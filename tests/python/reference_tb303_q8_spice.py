"""Component-level reference for the TB-303 Q8 square-wave shaper.

The circuit and component values follow the Roland TB-303 service notes and
the open-source x0xb0x reconstruction.  ngspice supplies the Gummel-Poon BJT
solve; the 2SA733 parameters are constrained by the Renesas data sheet and the
documented P-rank current-gain range.

Run directly to compare the circuit with ``Tb303SquareShaper``::

    uv run --with numpy python tests/python/reference_tb303_q8_spice.py

Set ``NGSPICE`` when the executable is outside PATH.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
from scipy.optimize import least_squares

try:
    import _triggerfish_dsp as dsp
except ImportError as error:  # pragma: no cover - developer setup diagnostic
    raise ImportError(
        "Build the development binding with "
        "'uv sync --group dev --reinstall-package triggerfish-vcv-dsp'"
    ) from error

ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "build" / "tb303-q8-reference"


def find_ngspice() -> Path:
    candidates = [
        os.environ.get("NGSPICE"),
        shutil.which("ngspice_con"),
        shutil.which("ngspice"),
        r"C:\msys64\ucrt64\bin\ngspice_con.exe",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return Path(candidate)
    raise FileNotFoundError("ngspice was not found; set NGSPICE to its executable")


def _netlist(
    frequency: float,
    beta: float,
    temperature: float,
    output_path: Path,
) -> str:
    # Renesas specifies VBE=0.62 V typical at 1 mA.  The corresponding IS at
    # 27 C is about 4e-14 A.  fT and Cob constrain the high-frequency terms;
    # their effect is minute beside the two external RC networks in audio use.
    output = output_path.as_posix()
    settle_cycles = max(30, math.ceil(0.25 * frequency))
    return f"""TB-303 Q8 square shaper reference
.param frequency={frequency:.17g}
.param period={{1/frequency}}

V12 supply 0 12
Vbias bias 0 5.333
Vsaw saw 0 PULSE(11.25 5.75 0 {{period*0.9998}} 1n 1n {{period}})

R45 supply emitter 22k
C11 emitter 0 1u
Q8 square base emitter Q2SA733
R36 square bias 10k
R35 saw base 100k
C10 saw ac 10n
R34 ac base 10k

.model Q2SA733 PNP(
+ IS=4e-14 BF={beta:.17g} VAF=100 IKF=50m
+ ISE=4e-13 NE=1.5 BR=10 VAR=25 IKR=10m
+ RB=100 RE=1 RC=1
+ CJE=20p VJE=.75 MJE=.33 CJC=4.5p VJC=.55 MJC=.33
+ TF=.8n TR=50n)

.temp {temperature:.17g}
.options reltol=1e-8 abstol=1e-12 vntol=1e-9
.tran {{period/2048}} {{period*{settle_cycles + 1}}} {{period*{settle_cycles}}} {{period/2048}}

.control
set wr_singlescale
set wr_vecnames
run
wrdata {output} time v(saw) v(base) v(emitter) v(square)
quit
.endc

.end
"""


def render_spice(
    frequency: float,
    beta: float = 300.0,
    temperature: float = 27.0,
    samples_per_cycle: int = 4096,
) -> tuple[np.ndarray, np.ndarray]:
    BUILD.mkdir(parents=True, exist_ok=True)
    stem = f"q8-{frequency:g}hz-b{beta:g}-t{temperature:g}"
    circuit = BUILD / f"{stem}.cir"
    data = BUILD / f"{stem}.dat"
    circuit.write_text(_netlist(frequency, beta, temperature, data), encoding="utf-8")
    result = subprocess.run(
        [str(find_ngspice()), "-b", str(circuit)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)

    raw = np.loadtxt(data, skiprows=1)
    time = raw[:, 0]
    saw = raw[:, 2]
    square = raw[:, -1]
    period = 1.0 / frequency
    end = time[-1]
    phase = np.arange(samples_per_cycle, dtype=float) / samples_per_cycle
    query = end - period + phase * period
    return np.interp(query, time, saw), np.interp(query, time, square)


def render_spice_dc(
    beta: float = 300.0,
    temperature: float = 27.0,
) -> tuple[np.ndarray, np.ndarray]:
    BUILD.mkdir(parents=True, exist_ok=True)
    circuit = BUILD / f"q8-dc-b{beta:g}-t{temperature:g}.cir"
    data = BUILD / f"q8-dc-b{beta:g}-t{temperature:g}.dat"
    output = data.as_posix()
    circuit.write_text(
        f"""TB-303 Q8 DC transfer reference
V12 supply 0 12
Vbias bias 0 5.333
Vsaw saw 0 8.5
R45 supply emitter 22k
C11 emitter 0 1u
Q8 square base emitter Q2SA733
R36 square bias 10k
R35 saw base 100k
C10 saw ac 10n
R34 ac base 10k
.model Q2SA733 PNP(
+ IS=4e-14 BF={beta:.17g} VAF=100 IKF=50m
+ ISE=4e-13 NE=1.5 BR=10 VAR=25 IKR=10m
+ RB=100 RE=1 RC=1
+ CJE=20p VJE=.75 MJE=.33 CJC=4.5p VJC=.55 MJC=.33
+ TF=.8n TR=50n)
.temp {temperature:.17g}
.options reltol=1e-9 abstol=1e-13 vntol=1e-10
.dc Vsaw 5.75 11.25 .002
.control
set wr_singlescale
set wr_vecnames
run
wrdata {output} v(saw) v(square)
quit
.endc
.end
""",
        encoding="utf-8",
    )
    result = subprocess.run(
        [str(find_ngspice()), "-b", str(circuit)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    raw = np.loadtxt(data, skiprows=1)
    return raw[:, 1], raw[:, -1]


def fit_static_tanh(beta: float = 300.0, temperature: float = 27.0) -> dict[str, float]:
    saw, square = render_spice_dc(beta, temperature)
    normalized_saw = (saw - 8.5) / 2.75

    def residual(parameters: np.ndarray) -> np.ndarray:
        offset, gain, threshold, log_width = parameters
        width = math.exp(float(log_width))
        predicted = offset + gain * np.tanh((threshold - normalized_saw) / width)
        return predicted - square

    solution = least_squares(
        residual,
        np.asarray([7.0, 2.0, 0.0, math.log(0.055)]),
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
        max_nfev=2000,
    )
    error = residual(solution.x)
    span = float(np.max(square) - np.min(square))
    return {
        "threshold": float(solution.x[2]),
        "width": math.exp(float(solution.x[3])),
        "rms_percent_span": 100.0 * float(np.sqrt(np.mean(error**2))) / span,
        "max_percent_span": 100.0 * float(np.max(np.abs(error))) / span,
    }


def render_model(
    frequency: float,
    samples_per_cycle: int = 4096,
    settle_cycles: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    if settle_cycles is None:
        settle_cycles = max(30, math.ceil(0.25 * frequency))
    sample_rate = frequency * samples_per_cycle
    saw_cycle = 2.0 * np.arange(samples_per_cycle) / samples_per_cycle - 1.0
    saw = np.tile(saw_cycle, settle_cycles + 1)
    output = dsp.tb303_q8(saw, frequency, 0.0, sample_rate)
    return saw_cycle, output[-samples_per_cycle:]


def harmonic_amplitudes(signal: np.ndarray, count: int = 15) -> np.ndarray:
    centered = signal - np.mean(signal)
    phase = 2.0 * np.pi * np.arange(centered.size) / centered.size
    amplitudes = []
    for harmonic in range(1, count + 1):
        coefficient = np.sum(centered * np.exp(-1j * harmonic * phase))
        amplitudes.append(2.0 * abs(coefficient) / centered.size)
    return np.asarray(amplitudes)


def midpoint_duty(signal: np.ndarray) -> float:
    midpoint = 0.5 * (float(np.min(signal)) + float(np.max(signal)))
    return float(np.mean(signal > midpoint))


def normalized_harmonics(signal: np.ndarray, count: int = 15) -> np.ndarray:
    amplitudes = harmonic_amplitudes(signal, count)
    return amplitudes / amplitudes[0]


def compare_case(frequency: float, beta: float, temperature: float) -> dict[str, float]:
    _, spice = render_spice(frequency, beta, temperature)
    _, reduced = render_model(frequency)
    reference = normalized_harmonics(spice)
    candidate = normalized_harmonics(reduced)
    vector_error = float(
        np.linalg.norm(candidate[1:] - reference[1:])
        / max(np.linalg.norm(reference[1:]), np.finfo(float).eps)
    )
    harmonic_db = 20.0 * np.log10(
        np.maximum(candidate, 1e-12) / np.maximum(reference, 1e-12)
    )
    return {
        "spice_duty": midpoint_duty(spice),
        "reduced_duty": midpoint_duty(reduced),
        "harmonic_vector_error": vector_error,
        "h2_db": float(harmonic_db[1]),
        "h3_db": float(harmonic_db[2]),
        "h5_db": float(harmonic_db[4]),
        "h9_db": float(harmonic_db[8]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--frequencies", type=float, nargs="+", default=[10, 85, 135, 1000]
    )
    parser.add_argument("--betas", type=float, nargs="+", default=[200, 300, 400])
    parser.add_argument("--temperature", type=float, default=27.0)
    args = parser.parse_args()

    print("Static transfer fit to an affine tanh:")
    for beta in args.betas:
        fit = fit_static_tanh(beta, args.temperature)
        print(
            f"  beta {beta:5.0f}: threshold {fit['threshold']:+.4f}, "
            f"width {fit['width']:.4f}, RMS {fit['rms_percent_span']:.3f}% "
            f"of span, max {fit['max_percent_span']:.3f}%"
        )

    print("frequency  beta  duty SPICE/model  harmonic error  H2/H3/H5/H9 delta")
    for frequency in args.frequencies:
        for beta in args.betas:
            result = compare_case(frequency, beta, args.temperature)
            print(
                f"{frequency:8.1f} {beta:5.0f}  "
                f"{result['spice_duty']:.3f}/{result['reduced_duty']:.3f}       "
                f"{100*result['harmonic_vector_error']:6.2f}%       "
                f"{result['h2_db']:+6.2f}/{result['h3_db']:+6.2f}/"
                f"{result['h5_db']:+6.2f}/{result['h9_db']:+6.2f} dB"
            )


if __name__ == "__main__":
    main()
