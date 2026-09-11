"""The optional inspector uses real packet preparation, with strict inputs."""

import json
from pathlib import Path
import subprocess

import pytest


@pytest.fixture(scope="module")
def executable():
    directory = Path(__file__).resolve().parents[2] / "build/workbench-api"
    for name in ("triggerfish_modal_inspect.exe", "triggerfish_modal_inspect"):
        if (directory / name).exists():
            return directory / name
    pytest.skip("Build the optional native workbench through dev.ps1 first")


def test_packet_state_budget_and_energy(executable):
    parameters = {f"resolved_level_{i}": -72 for i in range(32)}
    parameters.update(
        resolved_level_0=0,
        resolved_frequency_0=1000,
        field_satellite_density=1,
        field_turbulence=1,
        field_distribution=3,
        field_doublet_split=2,
    )
    request = "44100\n" + "".join(f"{k} {v}\n" for k, v in parameters.items())
    result = subprocess.run(
        [executable], input=request, text=True, capture_output=True, check=True
    )
    modes = json.loads(result.stdout)
    assert len(modes) == 512
    assert sum(m["input"] ** 2 for m in modes) == pytest.approx(1, abs=2e-6)
    assert all(m["packet"] == 0 and m["centre"] == 1000 for m in modes)
    assert all(1 <= m["frequency"] < 44100 * 0.49 for m in modes)


@pytest.mark.parametrize(
    "input_text", ["bad", "44100\nunknown 1", "44100\nfield_motion_depth 99"]
)
def test_invalid_inspection_is_not_silently_clamped(executable, input_text):
    result = subprocess.run(
        [executable], input=input_text, text=True, capture_output=True
    )
    assert result.returncode != 0
    assert result.stderr.strip()
