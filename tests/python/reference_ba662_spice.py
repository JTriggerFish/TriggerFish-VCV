"""Run the published BA662-clone topology as an offline ngspice reference.

This is a research/benchmark tool, not part of the plugin build or automated
test suite. It downloads the current official Nexperia matched-transistor
models into ``build/ba662-reference`` and compares these variants:

* ideal: the clone topology with ideal, infinite-bandwidth BJTs;
* static: Nexperia's device models with junction/transit capacitances removed;
* full: the unmodified Nexperia device models;
* reduced: TriggerFish's level-matched, quasi-static differential-pair law.

The manufacturer model files are never copied into the repository. Pass an
explicit ngspice executable when it is not available on PATH, for example
``--ngspice C:/path/to/ngspice``.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from reference_ba662 import CURRENT_TRANSFER_EFFICIENCY, differential_pair_current

ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "tests" / "spice" / "ba662_clone.cir.in"
DEFAULT_OUTPUT = ROOT / "build" / "ba662-reference"
MODEL_URLS = {
    "PMP4201Y": "https://assets.nexperia.com/documents/spice-model/PMP4201Y.txt",
    "PMP5201Y": "https://assets.nexperia.com/documents/spice-model/PMP5201Y.txt",
}


@dataclass(frozen=True)
class Case:
    name: str
    frequency_hz: float
    differential_rms_volts: float
    control_current_amps: float
    output_load_ohms: float = 220_000.0


CASES = (
    Case("nominal", 1_000.0, 5.0e-3, 20.0e-6),
    Case("low_control", 1_000.0, 5.0e-3, 5.0e-6),
    Case("accented", 1_000.0, 5.0e-3, 40.0e-6),
    Case("driven", 1_000.0, 20.0e-3, 20.0e-6),
    Case("overdrive", 1_000.0, 100.0e-3, 20.0e-6),
    Case("datasheet", 1_000.0, 5.0e-3, 200.0e-6, 27_000.0),
    Case("high_frequency", 10_000.0, 5.0e-3, 20.0e-6),
)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--ngspice",
        type=Path,
        help="ngspice executable; defaults to NGSPICE or PATH",
    )
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_OUTPUT, help="ignored output directory"
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="run only the nominal and high-frequency cases",
    )
    return parser.parse_args()


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
    # The MSYS2 package also ships a GUI-subsystem executable that does not
    # provide useful batch status/output. Prefer its console sibling.
    console_candidate = candidate.with_name(f"{candidate.stem}_con{candidate.suffix}")
    if console_candidate.is_file():
        candidate = console_candidate
    return candidate.resolve()


def download_models(model_directory: Path) -> dict[str, Path]:
    model_directory.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, url in MODEL_URLS.items():
        path = model_directory / f"{name}.txt"
        if not path.exists():
            print(f"Downloading official {name} model from Nexperia")
            urllib.request.urlretrieve(url, path)
        text = path.read_text(encoding="ascii", errors="strict")
        if f".SUBCKT {name} 1 2 3" not in text:
            raise RuntimeError(f"unexpected {name} model format: {path}")
        paths[name] = path.resolve()
    return paths


def write_static_model(source: Path, destination: Path):
    """Remove only charge-storage terms while retaining all DC parameters."""

    text = source.read_text(encoding="ascii", errors="strict")
    dynamic_parameters = ("CJE", "CJC", "TF", "TR")
    for parameter in dynamic_parameters:
        text = re.sub(
            rf"(?im)^(\+\s*{parameter}\s*=\s*)[^\s]+",
            rf"\g<1>0",
            text,
        )
    destination.write_text(text, encoding="ascii", newline="\n")


def device_models(variant: str, models: dict[str, Path], output: Path) -> str:
    if variant == "ideal":
        return """.model NPN_IDEAL NPN(Is=1e-14 Bf=1e9 Br=1e9 Vaf=1e12 Var=1e12)
.model PNP_IDEAL PNP(Is=1e-14 Bf=1e9 Br=1e9 Vaf=1e12 Var=1e12)
.subckt NPN_DEVICE c b e
Q1 c b e NPN_IDEAL
.ends NPN_DEVICE
.subckt PNP_DEVICE c b e
Q1 c b e PNP_IDEAL
.ends PNP_DEVICE"""

    selected = models
    if variant == "static":
        static_directory = output / "models-static"
        static_directory.mkdir(parents=True, exist_ok=True)
        selected = {}
        for name, source in models.items():
            destination = static_directory / source.name
            write_static_model(source, destination)
            selected[name] = destination.resolve()

    npn = selected["PMP4201Y"].as_posix()
    pnp = selected["PMP5201Y"].as_posix()
    return f""".include \"{npn}\"
.include \"{pnp}\"
.subckt NPN_DEVICE c b e
X1 c b e PMP4201Y
.ends NPN_DEVICE
.subckt PNP_DEVICE c b e
X1 c b e PMP5201Y
.ends PNP_DEVICE"""


def analysis_deck(case: Case, output_path: Path) -> str:
    differential_peak = np.sqrt(2.0) * case.differential_rms_volts
    half_peak = 0.5 * differential_peak
    # Keep 200 source samples per period and simulate 100 settled periods. The
    # first 20 periods are omitted from metrics to remove startup/coupling DC.
    step = 1.0 / (200.0 * case.frequency_hz)
    duration = 100.0 / case.frequency_hz
    start = 20.0 / case.frequency_hz
    return f"""VCC vp 0 12
VBIAS bias 0 5.333333333333
VINP p3 bias SIN(0 {half_peak:.17g} {case.frequency_hz:.17g})
VINN p2 bias SIN(0 {half_peak:.17g} {case.frequency_hz:.17g} 0 0 180)
IABC bias p1 {case.control_current_amps:.17g}
XOTA p1 p2 p3 bias 0 ota ota buffer vp BA662
RO ota bias {case.output_load_ohms:.17g}
C38 buffer coupled 1u
RLOAD coupled 0 50k

.control
set noaskquit
set wr_vecnames
set wr_singlescale
tran {step:.17g} {duration:.17g} {start:.17g} {step:.17g}
wrdata {output_path.name} v(coupled) v(buffer) v(ota) v(p3,p2)
quit
.endc
.end"""


def run_variant(
    ngspice: Path,
    variant: str,
    case: Case,
    models: dict[str, Path],
    output: Path,
):
    case_directory = output / case.name
    case_directory.mkdir(parents=True, exist_ok=True)
    data_path = (case_directory / f"{variant}.txt").resolve()
    netlist_path = case_directory / f"{variant}.cir"
    template = TEMPLATE.read_text(encoding="utf-8")
    netlist = template.replace(
        "__DEVICE_MODELS__", device_models(variant, models, output)
    ).replace("__ANALYSIS_DECK__", analysis_deck(case, data_path))
    netlist_path.write_text(netlist, encoding="utf-8", newline="\n")
    process = subprocess.run(
        [str(ngspice), "-b", str(netlist_path.resolve())],
        cwd=case_directory,
        env={
            **os.environ,
            "PATH": f"{ngspice.parent}{os.pathsep}{os.environ.get('PATH', '')}",
        },
        text=True,
        capture_output=True,
        check=False,
    )
    (case_directory / f"{variant}.log").write_text(
        process.stdout + process.stderr, encoding="utf-8", newline="\n"
    )
    if process.returncode != 0 or not data_path.exists():
        raise RuntimeError(
            f"ngspice failed for {case.name}/{variant}; "
            f"see {case_directory / f'{variant}.log'}"
        )
    data = np.genfromtxt(data_path, names=True)
    names = data.dtype.names
    if names is None or len(names) != 5:
        raise RuntimeError(f"unexpected ngspice output columns in {data_path}: {names}")
    return tuple(np.asarray(data[name], dtype=float) for name in names)


def rms(signal):
    return float(np.sqrt(np.mean(np.square(signal))))


def level_matched_residual(reference, candidate):
    # Remove a fractional pure delay before gain matching. This prevents C38's
    # expected linear phase, or ngspice sample-grid placement, from being
    # reported as nonlinear model error.
    reference_spectrum = np.fft.rfft(reference)
    candidate_spectrum = np.fft.rfft(candidate)
    fundamental = int(np.argmax(np.abs(reference_spectrum[1:])) + 1)
    phase_per_bin = (
        np.angle(reference_spectrum[fundamental] / candidate_spectrum[fundamental])
        / fundamental
    )
    bins = np.arange(candidate_spectrum.size)
    candidate = np.fft.irfft(
        candidate_spectrum * np.exp(1j * phase_per_bin * bins), n=reference.size
    )
    candidate_energy = float(np.dot(candidate, candidate))
    gain = (
        float(np.dot(reference, candidate)) / candidate_energy
        if candidate_energy > 0.0
        else 0.0
    )
    return gain, rms(reference - gain * candidate) / rms(reference)


def harmonic_metrics(signal, time, frequency):
    # Resample adaptive SPICE output to a uniform grid containing an integer
    # number of periods before applying a rectangular-window Fourier series.
    period_count = 64
    sample_count = period_count * 256
    end = time[-1]
    start = end - period_count / frequency
    uniform_time = np.linspace(start, end, sample_count, endpoint=False)
    uniform = np.interp(uniform_time, time, signal)
    uniform -= np.mean(uniform)
    spectrum = np.fft.rfft(uniform)
    fundamental_bin = period_count
    harmonics = np.asarray(
        [abs(spectrum[harmonic * fundamental_bin]) for harmonic in range(1, 10)]
    )
    thd = np.linalg.norm(harmonics[1:]) / harmonics[0]
    even = np.linalg.norm(harmonics[1::2]) / harmonics[0]
    return uniform_time, uniform, float(thd), float(even)


def harmonic_projection(signal, period_count=64, harmonic_count=9):
    """Keep only the measured periodic harmonics used by the comparison.

    Adaptive SPICE output is interpolated onto a uniform grid. Comparing every
    FFT bin would count interpolation and solver noise as circuit-model error,
    even though a settled sine-driven circuit can contain only integer
    harmonics. Projecting both candidates onto the same harmonic basis makes
    the null measure the nonlinear transfer instead of the simulator floor.
    """

    spectrum = np.fft.rfft(signal)
    projected = np.zeros_like(spectrum)
    for harmonic in range(1, harmonic_count + 1):
        projected[harmonic * period_count] = spectrum[harmonic * period_count]
    return np.fft.irfft(projected, n=signal.size)


def analyze_case(ngspice, case, models, output):
    rendered = {}
    projected = {}
    metrics = {}
    full_stages = None
    for variant in ("ideal", "static", "full"):
        time, coupled, buffer, ota, differential = run_variant(
            ngspice, variant, case, models, output
        )
        uniform_time, uniform, thd, even = harmonic_metrics(
            coupled, time, case.frequency_hz
        )
        rendered[variant] = uniform
        projected[variant] = harmonic_projection(uniform)
        metrics[variant] = (rms(uniform), thd, even)
        if variant == "full":
            buffer_uniform = harmonic_metrics(buffer, time, case.frequency_hz)[1]
            ota_uniform = harmonic_metrics(ota, time, case.frequency_hz)[1]
            full_stages = (ota_uniform, buffer_uniform)

    reference_time = uniform_time
    differential = (
        np.sqrt(2.0)
        * case.differential_rms_volts
        * np.sin(2.0 * np.pi * case.frequency_hz * reference_time)
    )
    reduced = differential_pair_current(
        differential,
        np.full_like(differential, case.control_current_amps),
        efficiency=CURRENT_TRANSFER_EFFICIENCY,
    )
    reduced -= np.mean(reduced)
    reduced_thd, reduced_even = harmonic_metrics(
        reduced, reference_time, case.frequency_hz
    )[2:]
    ota_uniform, buffer_uniform = full_stages
    reduced_projected = harmonic_projection(reduced)
    gain, reduced_residual = level_matched_residual(
        harmonic_projection(ota_uniform), reduced_projected
    )
    _, static_residual = level_matched_residual(projected["full"], projected["static"])
    _, ideal_residual = level_matched_residual(projected["full"], projected["ideal"])
    _, buffer_residual = level_matched_residual(
        harmonic_projection(buffer_uniform), harmonic_projection(ota_uniform)
    )
    _, coupling_residual = level_matched_residual(
        projected["full"], harmonic_projection(buffer_uniform)
    )

    print(
        f"\n{case.name}: {case.frequency_hz:g} Hz, "
        f"{1e3 * case.differential_rms_volts:g} mV RMS differential, "
        f"{1e6 * case.control_current_amps:g} uA control, "
        f"{1e-3 * case.output_load_ohms:g} kohm load"
    )
    print("variant       RMS V       THD %      even %")
    for variant in ("ideal", "static", "full"):
        level, thd, even = metrics[variant]
        print(f"{variant:10s} {level:10.6f} {100 * thd:10.5f} {100 * even:10.5f}")
    print(
        f"reduced    gain={gain:10.3f} A->V  THD={100 * reduced_thd:8.5f}% "
        f"even={100 * reduced_even:8.5f}%"
    )
    print(
        "nine-harmonic residual vs full coupled output: "
        f"ideal={100 * ideal_residual:.5f}%  "
        f"static={100 * static_residual:.5f}%"
    )
    print(
        "reduced/full-OTA residual after harmonic/level/delay match: "
        f"{100 * reduced_residual:.5f}%"
    )
    print(
        "full-model stage residual after harmonic/level/delay match: "
        f"{100 * buffer_residual:.5f}%"
    )
    print(
        "C38 coupled/buffer residual after harmonic/level/delay match: "
        f"{100 * coupling_residual:.5f}%"
    )


def main():
    arguments = parse_arguments()
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    ngspice = find_ngspice(arguments.ngspice)
    models = download_models(output / "models")
    cases = (CASES[0], CASES[-1]) if arguments.quick else CASES
    print(f"ngspice: {ngspice}")
    print(f"reference output: {output}")
    for case in cases:
        analyze_case(ngspice, case, models, output)


if __name__ == "__main__":
    main()
