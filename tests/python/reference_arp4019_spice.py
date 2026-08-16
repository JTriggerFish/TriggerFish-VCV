"""Compare the ARP 4019 input-pair reduction with a device-level proxy.

This manual reference sweeps a matched pair of manufacturer-model 2N3906 PNP
transistors at five tail currents, fits the effective thermal voltage of the
``tanh`` reduction, and reports normalized-current residuals. It writes only
generated netlists, logs, and data beneath ``build/arp4019-reference``.

The 2N3906 model is the Gummel-Poon model published in the onsemi data sheet:
https://www.onsemi.com/pdf/datasheet/2n3906-d.pdf
It is a modern device-level proxy for the 4019's matched TZ-581 pair, not a
reconstruction of the original transistor die or the complete 4019 module.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from reference_arp4019 import circuit_values

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "build" / "arp4019-reference"
TEMPERATURE_C = 27.0
PRODUCTION_THERMAL_VOLTAGE = 25.85e-3
SWEEP_LIMIT_VOLTS = 0.150
SWEEP_STEP_VOLTS = 0.00025

# Manufacturer-published 2N3906 model. Capacitances do not affect the DC sweep,
# but retaining the complete model makes the generated deck unambiguous.
MODEL = """.model Q2N3906 PNP(
+ Is=1.41f Xti=3 Eg=1.11 Vaf=18.7 Bf=180.7 Ne=1.5 Ise=0 Ikf=80m
+ Xtb=1.5 Br=4.977 Nc=2 Isc=0 Ikr=0 Rc=2.5 Cjc=9.728p Mjc=.5776
+ Vjc=.75 Fc=.5 Cje=8.063p Mje=.3677 Vje=.75 Tr=33.42n Tf=179.3p
+ Itf=.4 Vtf=4 Xtf=6 Rb=10)"""


def find_ngspice(explicit: Path | None) -> Path:
    if explicit is not None:
        candidate = explicit
    else:
        candidate_text = os.environ.get("NGSPICE") or shutil.which("ngspice")
        if candidate_text is None:
            raise RuntimeError(
                "ngspice was not found; set NGSPICE or pass --ngspice PATH"
            )
        candidate = Path(candidate_text)
    if not candidate.is_file():
        raise RuntimeError(f"ngspice executable does not exist: {candidate}")
    console_candidate = candidate.with_name(f"{candidate.stem}_con{candidate.suffix}")
    if console_candidate.is_file():
        candidate = console_candidate
    return candidate.resolve()


def netlist(tail_current: float, data_path: Path) -> str:
    return f"""ARP 4019 matched-PNP input-pair proxy
.temp {TEMPERATURE_C:g}
VCC supply 0 15
VLEFT collector_left 0 0
VRIGHT collector_right 0 0
VDIFF differential 0 0
BLEFT base_left 0 V=0.5*V(differential)
BRIGHT base_right 0 V=-0.5*V(differential)
QLEFT collector_left base_left emitters Q2N3906
QRIGHT collector_right base_right emitters Q2N3906
ITAIL supply emitters DC {tail_current:.17g}
{MODEL}
.control
set noaskquit
set wr_vecnames
set wr_singlescale
dc VDIFF {-SWEEP_LIMIT_VOLTS:.17g} {SWEEP_LIMIT_VOLTS:.17g} {SWEEP_STEP_VOLTS:.17g}
wrdata {data_path.name} i(VLEFT) i(VRIGHT)
quit
.endc
.end
"""


def run_sweep(ngspice: Path, tail_current: float, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    label = f"{tail_current * 1e6:.6f}uA".replace(".", "_")
    data_path = (output / f"pair-{label}.txt").resolve()
    netlist_path = output / f"pair-{label}.cir"
    log_path = output / f"pair-{label}.log"
    netlist_path.write_text(
        netlist(tail_current, data_path), encoding="ascii", newline="\n"
    )
    process = subprocess.run(
        [str(ngspice), "-b", str(netlist_path.resolve())],
        cwd=output,
        env={
            **os.environ,
            "PATH": f"{ngspice.parent}{os.pathsep}{os.environ.get('PATH', '')}",
        },
        text=True,
        capture_output=True,
        check=False,
    )
    log_path.write_text(process.stdout + process.stderr, encoding="utf-8", newline="\n")
    if process.returncode != 0 or not data_path.is_file():
        raise RuntimeError(f"ngspice failed; see {log_path}")

    data = np.loadtxt(data_path, skiprows=1)
    if data.ndim != 2 or data.shape[1] != 3:
        raise RuntimeError(
            f"unexpected ngspice output shape in {data_path}: {data.shape}"
        )
    differential = data[:, 0]
    left = data[:, 1]
    right = data[:, 2]
    collector_sum = np.abs(left) + np.abs(right)
    normalized = (left - right) / collector_sum
    if np.corrcoef(differential, normalized)[0, 1] < 0.0:
        normalized = -normalized
    return differential, normalized


def residuals(differential, normalized, thermal_voltage):
    reduced = np.tanh(differential / (2.0 * thermal_voltage))
    error = normalized - reduced
    return float(np.sqrt(np.mean(error**2))), float(np.max(np.abs(error)))


def fit_thermal_voltage(differential, normalized):
    result = minimize_scalar(
        lambda thermal_voltage: residuals(differential, normalized, thermal_voltage)[0],
        bounds=(0.015, 0.040),
        method="bounded",
        options={"xatol": 1.0e-14},
    )
    if not result.success:
        raise RuntimeError(f"thermal-voltage fit failed: {result.message}")
    return float(result.x)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ngspice",
        type=Path,
        help="ngspice executable; defaults to NGSPICE or PATH",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()

    ngspice = find_ngspice(arguments.ngspice)
    output = arguments.output.resolve()
    unity_current = circuit_values()["unity_control_current_amps"]
    tail_currents = unity_current * np.asarray((0.001, 0.01, 0.1, 1.0, 3.0))

    print(f"ngspice: {ngspice}")
    print(f"reference output: {output}")
    print("tail current | fitted VT | fitted RMS/max | fixed-25.85mV RMS/max")
    for tail_current in tail_currents:
        differential, normalized = run_sweep(ngspice, tail_current, output)
        fitted_voltage = fit_thermal_voltage(differential, normalized)
        fitted_rms, fitted_maximum = residuals(differential, normalized, fitted_voltage)
        fixed_rms, fixed_maximum = residuals(
            differential, normalized, PRODUCTION_THERMAL_VOLTAGE
        )
        print(
            f"{tail_current * 1e6:9.4f} uA | {fitted_voltage * 1e3:9.5f} mV | "
            f"{100 * fitted_rms:8.5f}/{100 * fitted_maximum:8.5f}% | "
            f"{100 * fixed_rms:8.5f}/{100 * fixed_maximum:8.5f}%"
        )


if __name__ == "__main__":
    main()
