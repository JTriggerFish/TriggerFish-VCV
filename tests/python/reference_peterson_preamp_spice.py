"""Offline transistor reference for the Figure 11-8 Peterson audio preamp.

The netlist includes the selected-2N3392 input pair, the actual bass/treble
bridge, the Darlington-connected selected-2N3392 feedback pair and its Q5
emitter follower.  It is the independent reference for
PetersonPreAmplifier.hpp plus PetersonTonePreAmplifier.hpp; in particular, no
post-EQ transfer curve or memoryless waveshaper stands in for tone-stage
overload.

Requires ngspice on PATH. Example:
    uv run --with numpy python tests/python/reference_peterson_preamp_spice.py
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "build" / "peterson-preamp-spice"


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
    data_path: Path,
    bass: float = 0.5,
    treble: float = 0.5,
) -> str:
    # The service sheet identifies selected 2N3392 devices but supplies no lot
    # curves. BF=220 lies inside the period datasheet's 150--300 range at 2 mA;
    # BR and Is match the real-time Ebers-Moll proxy and remain measurement gaps.
    pot_min = 1.0
    treble_top = pot_min + (1.0 - treble) * 100_000.0
    treble_bottom = pot_min + treble * 100_000.0
    bass_top = pot_min + (1.0 - bass) * 100_000.0
    bass_bottom = pot_min + bass * 100_000.0
    return f"""Peterson Figure 11-8 nonlinear audio preamp and tone feedback
VRAIL rail 0 25
VIN input 0 SIN(0 {amplitude:.12g} {frequency:.12g})
CIN input base 0.22u
RBIAS rail base 390k
RBG base 0 33k
RC1 rail collector 12k
RE1 emitter1 0 1.5k
CMILLER base collector 330p
Q1 collector base emitter1 Q2N3392
Q2 rail collector emitter2 Q2N3392
RE2 emitter2 0 4.7k
COUT emitter2 tonein 5u

* Figure 11-8 100k linear controls and frequency-dependent feedback bridge.
RTIN tonein tt 68k
RBIN tonein bt 6.8k
RTOP tt tw {treble_top:.12g}
RTBOT tw tb {treble_bottom:.12g}
R2B tb toneout 6.8k
CTREBLE tw tx 0.0047u
RBRIDGE tx bw 22k
CFBTREBLE tx fb 0.047u
RBASSTOP bt bw {bass_top:.12g}
RBASSBOT bw bb {bass_bottom:.12g}
CBASSTOP bw bt 0.1u
CBASSBOT bw bb 0.1u
R2A bb toneout 6.8k
RFBGROUND fb 0 390k
RFB toneout fb 1meg

* Q4/Q3 are the rotated Darlington common-emitter pair in Figure 11-8.
RC rail toneout 12k
Q4 toneout fb q3base Q2N3392
Q3 toneout q3base 0 Q2N3392

* Compensated Q5 follower and the dominant historical volume loading.
RZ toneout z 6.8k
RB5 z base5 68k
CZ z emitter5 0.0047u
CB5 base5 0 470p
Q5 rail base5 emitter5 Q2N3392
RE5 emitter5 0 6.8k
RVOLUME emitter5 0 100k
.model Q2N3392 NPN(Is=2e-14 Bf=220 Br=3 Nf=1 Nr=1)
.options abstol=1e-11 reltol=1e-6 vntol=1e-7
.tran 1.30208333333333u 0.24 0.08
.control
set wr_singlescale
set wr_vecnames
run
wrdata {data_path.as_posix()} time v(input) v(emitter5)
quit
.endc
.end
"""


def input_pair_netlist(amplitude: float, frequency: float, data_path: Path) -> str:
    """The Q1/Q2 reduction used by PetersonPreAmplifier.hpp in isolation."""
    return f"""Peterson Figure 11-8 input transistor pair
VRAIL rail 0 25
VIN input 0 SIN(0 {amplitude:.12g} {frequency:.12g})
CIN input base 0.22u
RBIAS rail base 390k
RBG base 0 33k
RC1 rail collector 12k
RE1 emitter1 0 1.5k
CMILLER base collector 330p
Q1 collector base emitter1 Q2N3392
Q2 rail collector emitter2 Q2N3392
RE2 emitter2 0 4.7k
COUT emitter2 output 5u
RLOAD output 0 64k
.model Q2N3392 NPN(Is=2e-14 Bf=220 Br=3 Nf=1 Nr=1)
.options abstol=1e-11 reltol=1e-6 vntol=1e-7
.tran 1.30208333333333u 0.24 0.08
.control
set wr_singlescale
set wr_vecnames
run
wrdata {data_path.as_posix()} time v(input) v(output)
quit
.endc
.end
"""


def harmonic_metrics(
    signal: np.ndarray, sample_rate: float, frequency: float
) -> dict[str, float]:
    signal = signal - np.mean(signal)
    phase = np.arange(signal.size) / sample_rate

    def magnitude(harmonic: int) -> float:
        omega = 2.0 * math.pi * frequency * harmonic * phase
        cosine = 2.0 * np.dot(signal, np.cos(omega)) / signal.size
        sine = 2.0 * np.dot(signal, np.sin(omega)) / signal.size
        return float(math.hypot(cosine, sine))

    fundamental = max(1.0e-20, magnitude(1))
    harmonics = [100.0 * magnitude(index) / fundamental for index in range(2, 8)]
    return {
        "fundamental_vpk": fundamental,
        "thd_percent": float(math.sqrt(sum(value * value for value in harmonics))),
        **{f"h{index}_percent": harmonics[index - 2] for index in range(2, 8)},
    }


def render(
    amplitude: float, frequency: float, ngspice: Path, section: str
) -> dict[str, float]:
    BUILD.mkdir(parents=True, exist_ok=True)
    stem = f"{section}-a{amplitude:g}-f{frequency:g}"
    circuit_path = BUILD / f"{stem}.cir"
    data_path = BUILD / f"{stem}.dat"
    circuit = (
        input_pair_netlist(amplitude, frequency, data_path)
        if section == "input"
        else netlist(amplitude, frequency, data_path)
    )
    circuit_path.write_text(circuit, encoding="utf-8")
    completed = subprocess.run(
        [str(ngspice), "-b", str(circuit_path)],
        cwd=BUILD,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0 or not data_path.exists():
        raise RuntimeError(
            f"ngspice failed ({completed.returncode})\n"
            f"{completed.stdout}\n{completed.stderr}"
        )
    raw = np.loadtxt(data_path, skiprows=1)
    # wrdata repeats the scale column for every saved vector. The final column
    # is v(output), independent of ngspice's repeated-scale convention.
    time = raw[:, 0]
    output = raw[:, -1]
    uniform_rate = 96000.0
    uniform_time = np.arange(time[0], time[-1], 1.0 / uniform_rate)
    uniform_output = np.interp(uniform_time, time, output)
    metrics = harmonic_metrics(uniform_output, uniform_rate, frequency)
    metrics.update(
        {
            "amplitude_vpk": amplitude,
            "frequency_hz": frequency,
            "output_ac_rms_v": float(np.std(uniform_output)),
            "output_peak_v": float(np.max(np.abs(uniform_output))),
        }
    )
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frequency", type=float, default=250.0)
    parser.add_argument("--ngspice", type=Path)
    parser.add_argument("--section", choices=("input", "full"), default="full")
    parser.add_argument(
        "--amplitudes", type=float, nargs="*", default=(0.01, 0.20, 0.50, 1.00, 1.50)
    )
    args = parser.parse_args()
    ngspice = find_ngspice(args.ngspice)
    print(
        json.dumps(
            [
                render(amplitude, args.frequency, ngspice, args.section)
                for amplitude in args.amplitudes
            ],
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
