"""Offline transistor-level reference for the Peterson Figure 11-9 module.

This deck is the authority for reducing the real-time solve; it is not a fitted
audio waveshaper.  Component connectivity and values come from Figure 11-9 of
the Fender Rhodes service manual.  The two 120725 parts are the same matched
PNP germanium device: substitution manuals identify it as DTG110B, whose short-
form data list hFE=65--300 at 1 A and fT around 500 kHz.

The manual does not specify T1.  Its magnetising inductance and winding
resistance are therefore explicit assumptions and must eventually be replaced
by impedance measurements from an original transformer.  Production uses a
linear core: the former unmeasured sinh saturation branch is retained only as
an audit variant.  Likewise, the Gummel-Poon model is constrained by the
DTG110B short-form limits but is not a manufacturer SPICE model.

Run with::

    uv run --with numpy python tests/python/reference_peterson_power_spice.py
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "build" / "peterson-power-reference"


def find_ngspice(explicit: Path | None) -> Path:
    candidates = (
        explicit,
        Path(os.environ["NGSPICE"]) if os.environ.get("NGSPICE") else None,
        Path(path) if (path := shutil.which("ngspice_con")) else None,
        Path(path) if (path := shutil.which("ngspice")) else None,
        Path(r"C:\msys64\ucrt64\bin\ngspice_con.exe"),
    )
    for candidate in candidates:
        if candidate is not None and candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError("ngspice was not found; pass --ngspice PATH")


def netlist(
    amplitude: float,
    frequency: float,
    output: Path,
    primary_inductance: float,
    variant: str = "full",
) -> str:
    duration = 60.0 / frequency
    start = 40.0 / frequency
    step = 1.0 / (frequency * 500.0)
    data = output.as_posix()
    core_branch = f"31m+V(flux)/{primary_inductance:.17g}"
    if variant == "nonlinear_core":
        core_branch += "+.012*sinh(V(flux)/.032)"
    winding_resistance = 2.7
    output_capacitors = (
        ""
        if variant == "no_output_caps"
        else """
CBCUP bu out .01u
CBCLOW bl cl .01u
"""
    )
    if variant == "resistive_load":
        load = "RLOAD out 0 16"
    else:
        load = """RCOIL out coil 12.8
LCOIL coil motion .55m
RMOTION motion 0 52.0047
CMOTION motion 0 163.265u
LMOTION motion 0 27.585m"""
    return f"""Peterson Figure 11-9 transistor reference
V25 v25 0 25
VP vp 0 35
VN vn 0 -35
VIN in 0 SIN(0 {amplitude:.17g} {frequency:.17g})
RIN in 0 10k
CIN in b1 6.4u
RBIAS b1 c1 68k
Q1 c1 b1 e1 SMALLNPN
RC1 v25 c1 6.8k
RE1 e1 0 100
RFB out e1 5.6k
CFB out e1 1200p
Q2 c2 c1 e2 SMALLNPN
RE2 e2 0 47
RCE2 c2 e2 22k
CCB2 c1 c2 100p
RCORE v25 c2 42k
* Ideal winding constraints plus explicit magnetising branch.  The production
* core is linear; nonlinear_core restores the former provisional sinh term.
* 31 mA term is the DC class-A primary current at the operating point.
EUP us vp VALUE={{V(v25,c2)/1.55}}
RSU us bu {winding_resistance:.17g}
ELOW out ls VALUE={{V(v25,c2)/1.55}}
RSL ls bl {winding_resistance:.17g}
FUP v25 c2 EUP -0.6451612903225806
FLOW v25 c2 ELOW -0.6451612903225806
BINT 0 flux I={{V(v25,c2)}}
CFLUX flux 0 1 IC=0
RFLUX flux 0 1G
BMAG v25 c2 I={{{core_branch}}}
QUP out bu eu DTG110B_PROXY
REUP vp eu .5
RBEUP eu bu 820
QLOW cl bl out DTG110B_PROXY
RCLOW cl vn .5
RBELOW out bl 820
{output_capacitors}

* Linear moving-coil equivalent of the provisional 16-ohm Suitcase load.
{load}
RBLEED out 0 270

.model SMALLNPN NPN(IS=4e-14 BF=180 BR=4 VAF=100
+ CJE=12p CJC=4p TF=.7n TR=30n)
.model DTG110B_PROXY PNP(IS=5e-4 BF=90 BR=3 VAF=35 IKF=3
+ ISE=1m NE=1.5 CJE=1n CJC=10n TF=.5u TR=6u)
.options reltol=1e-7 abstol=1e-10 vntol=1e-8 method=gear
.tran {step:.17g} {duration:.17g} {start:.17g} {step:.17g}
.control
set wr_singlescale
set wr_vecnames
run
wrdata {data} time v(in) v(out) i(VP) i(VN)
quit
.endc
.end
"""


def render(
    ngspice: Path,
    directory: Path,
    amplitude: float,
    frequency: float,
    primary_inductance: float,
    variant: str = "full",
) -> dict[str, float]:
    stem = f"{variant}-a{amplitude:g}-f{frequency:g}"
    circuit = directory / f"{stem}.cir"
    data = directory / f"{stem}.dat"
    circuit.write_text(
        netlist(amplitude, frequency, data, primary_inductance, variant),
        encoding="ascii",
        newline="\n",
    )
    result = subprocess.run(
        [str(ngspice), "-b", str(circuit)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=90,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    if not data.is_file() or data.stat().st_size == 0:
        raise RuntimeError(result.stdout + result.stderr)
    raw = np.loadtxt(data, skiprows=1)
    # wrdata may repeat the scale before each requested vector; take the final
    # value column for output and the first scale column for time.
    time = raw[:, 0]
    output = raw[:, -5] if raw.shape[1] >= 10 else raw[:, -3]
    query_count = 20 * 2048
    end = float(time[-1])
    query = np.linspace(end - 20.0 / frequency, end, query_count, endpoint=False)
    raw_signal = np.interp(query, time, output)
    dc_offset = float(np.mean(raw_signal))
    signal = raw_signal - dc_offset
    spectrum = np.fft.rfft(signal)
    fundamental_bin = 20
    fundamental = abs(spectrum[fundamental_bin])
    harmonics = {
        f"h{harmonic}_percent": 100.0
        * abs(spectrum[harmonic * fundamental_bin])
        / max(1.0e-30, fundamental)
        for harmonic in range(2, 8)
    }
    thd = math.sqrt(sum(value * value for value in harmonics.values()))
    return {
        "variant": variant,
        "amplitude_vpk": amplitude,
        "frequency_hz": frequency,
        "output_dc_v": dc_offset,
        "output_ac_rms_v": float(np.sqrt(np.mean(signal * signal))),
        "output_raw_min_v": float(np.min(raw_signal)),
        "output_raw_max_v": float(np.max(raw_signal)),
        "output_ac_min_v": float(np.min(signal)),
        "output_ac_max_v": float(np.max(signal)),
        "thd_percent": thd,
        **harmonics,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ngspice", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--primary-inductance", type=float, default=8.0)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--amplitudes", type=float, nargs="*")
    parser.add_argument(
        "--variants",
        nargs="*",
        choices=("full", "nonlinear_core", "resistive_load", "no_output_caps"),
    )
    parser.add_argument(
        "--ablations",
        action="store_true",
        help="compare topology reductions at 0.2 and 0.5 V peak",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    amplitudes = (
        tuple(args.amplitudes)
        if args.amplitudes
        else (
            (0.20, 0.50)
            if args.ablations
            else (
                (0.001, 0.01, 0.05)
                if args.quick
                else (
                    0.001,
                    0.01,
                    0.05,
                    0.10,
                    0.25,
                    0.50,
                )
            )
        )
    )
    variants = (
        tuple(args.variants)
        if args.variants
        else (
            (
                "full",
                "nonlinear_core",
                "resistive_load",
                "no_output_caps",
            )
            if args.ablations
            else ("full",)
        )
    )
    results = [
        render(
            find_ngspice(args.ngspice),
            args.output,
            amplitude,
            250.0,
            args.primary_inductance,
            variant,
        )
        for variant in variants
        for amplitude in amplitudes
    ]
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
