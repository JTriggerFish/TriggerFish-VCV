from types import SimpleNamespace
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from triggerfish_percussion.observation_energy_gate import ObservationEnergyGate


def fixture():
    rng = np.random.default_rng(714)
    columns = rng.normal(size=(3, 200))
    amplitudes = np.array([0.2, 0.3, 0.4])
    audio = amplitudes @ columns
    basis = SimpleNamespace(bases={None: (audio, columns)}, amplitudes=amplitudes)
    return basis, audio


def test_identity_and_literal_gain():
    basis, reference = fixture()
    gate = ObservationEnergyGate(basis, reference, 1000, bins=((0, 0.1), (0.1, 0.2)))
    margins, _ = gate.evaluate(basis.amplitudes)
    assert np.allclose(margins, 3.5)
    assert np.allclose(gate.actual_errors(reference * 2), 20 * np.log10(2))
    assert gate.evaluate(basis.amplitudes * 2)[0].min() < 0


def test_gate_jacobian_matches_finite_differences():
    basis, reference = fixture()
    gate = ObservationEnergyGate(basis, reference, 1000, bins=((0, 0.1), (0.1, 0.2)))
    x = basis.amplitudes * 1.2
    _, jacobian = gate.evaluate(x)
    for index in range(3):
        step = np.zeros(3)
        step[index] = 1e-6
        numeric = (gate.evaluate(x + step)[0] - gate.evaluate(x - step)[0]) / (2e-6)
        assert np.allclose(jacobian[:, index], numeric, rtol=1e-6, atol=1e-7)


def test_seeds_are_constrained_separately():
    basis, reference = fixture()
    audio, columns = basis.bases[None]
    basis.bases[2] = (-2 * audio, -2 * columns)
    gate = ObservationEnergyGate(basis, reference, 1000, bins=((0, 0.1),))
    margins, _ = gate.evaluate(basis.amplitudes)
    assert np.allclose(margins[:2], 3.5)
    assert margins[2] < 0


def test_invalid_bins_rejected():
    basis, reference = fixture()
    with pytest.raises(ValueError):
        ObservationEnergyGate(basis, reference, 1000, bins=((0.1, 0.3),))


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_nonfinite_inputs_fail_closed(value):
    basis, reference = fixture()
    gate = ObservationEnergyGate(basis, reference, 1000, bins=((0, 0.1),))
    amplitudes = basis.amplitudes.copy()
    amplitudes[0] = value
    with pytest.raises(ValueError):
        gate.evaluate(amplitudes)
    samples = reference.copy()
    samples[0] = value
    with pytest.raises(ValueError):
        gate.actual_errors(samples)
    with pytest.raises(ValueError):
        ObservationEnergyGate(basis, samples, 1000)


def test_empty_or_mismatched_inputs_rejected():
    basis, reference = fixture()
    for rate in (0, -1, np.nan):
        with pytest.raises(ValueError):
            ObservationEnergyGate(basis, reference, rate)
    with pytest.raises(ValueError):
        ObservationEnergyGate(basis, reference, 1000, bins=())
    with pytest.raises(ValueError):
        ObservationEnergyGate(basis, reference[:-1], 1000)
    basis.sample_rate = 2000
    with pytest.raises(ValueError):
        ObservationEnergyGate(basis, reference, 1000)
    basis.sample_rate = 1000
    basis.bases = {}
    with pytest.raises(ValueError):
        ObservationEnergyGate(basis, reference, 1000)


def test_invalid_basis_and_overflow_rejected():
    basis, reference = fixture()
    gate = ObservationEnergyGate(basis, reference, 1000, bins=((0, 0.1),))
    with pytest.raises(ValueError):
        gate.actual_errors(np.full(200, 1e300))
    with pytest.raises(ValueError):
        gate.actual_errors(reference[:-1])
    with pytest.raises(ValueError):
        gate.evaluate(np.ones(2))
    with pytest.raises(ValueError):
        ObservationEnergyGate(basis, np.full(200, 1e300), 1000)
    basis.bases[None][1][0, 0] = np.nan
    with pytest.raises(ValueError):
        ObservationEnergyGate(basis, reference, 1000)


def test_floor_jacobian_is_zero():
    basis, reference = fixture()
    gate = ObservationEnergyGate(basis, reference, 1000, bins=((0, 0.1),))
    margins, jacobian = gate.evaluate(np.zeros(3))
    assert np.isfinite(margins).all()
    assert np.all(jacobian == 0)


def fitting_tool():
    pytest.importorskip("torch")  # Only the optional optimizer needs Torch.
    path = Path(__file__).resolve().parents[2] / "tools/refine_workbench_attack_gate.py"
    spec = importlib.util.spec_from_file_location("attack_gate_tool", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_optimizer_rejects_out_of_bounds_start():
    tool = fitting_tool()
    basis, _ = fixture()
    basis.amplitudes[0] = 0
    with pytest.raises(ValueError, match="explicit"):
        tool.gated_weights(basis, None, None, 1)


def test_optimizer_retains_feasible_start():
    tool = fitting_tool()
    basis, reference = fixture()
    basis.keys = ["a", "b", "c"]
    gate = ObservationEnergyGate(basis, reference, 1000, bins=((0, 0.1),))
    weights, result = tool.gated_weights(basis, lambda x: (x * x).mean(), gate, 2)
    assert np.isfinite(weights).all()
    assert gate.evaluate(weights)[0].min() >= -1e-6
    assert np.isfinite(result["best_feasible_score"])


def test_optimizer_nonfinite_result_rejected(monkeypatch):
    tool = fitting_tool()
    basis, reference = fixture()
    basis.keys = ["a", "b", "c"]
    gate = ObservationEnergyGate(basis, reference, 1000, bins=((0, 0.1),))
    monkeypatch.setattr(
        tool, "minimize", lambda *a, **k: SimpleNamespace(x=np.array([np.nan]), fun=0)
    )
    with pytest.raises(ValueError, match="optimizer result"):
        tool.gated_weights(basis, None, gate, 1)


def test_tool_refuses_input_overwrite(tmp_path):
    tool = fitting_tool()
    args = SimpleNamespace(start=tmp_path, output=tmp_path)
    with pytest.raises(ValueError, match="separate"):
        tool.refine(args)


def contact_tool(monkeypatch):
    folder = Path(__file__).resolve().parents[2] / "tools"
    monkeypatch.syspath_prepend(str(folder))
    spec = importlib.util.spec_from_file_location(
        "contact_shape_tool", folder / "refine_workbench_contact_shape.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_contact_shape_finite_difference(monkeypatch):
    tool = contact_tool(monkeypatch)
    x = np.array([0.2, 0.3, 0.4])
    assert np.allclose(
        tool.derivative(x, lambda a: np.array([np.sum(a * a), np.sum(a)])),
        np.array([2 * x, np.ones(3)]),
    )


def test_contact_shape_invalid_audio(monkeypatch):
    tool = contact_tool(monkeypatch)
    for samples in (np.ones(10), np.full(200, np.nan), np.ones((2, 200))):
        with pytest.raises(ValueError):
            tool.early_power_db(samples, 1000)


def test_contact_shape_refuses_overwrite(monkeypatch, tmp_path):
    tool = contact_tool(monkeypatch)
    with pytest.raises(ValueError, match="separate"):
        tool.refine(SimpleNamespace(start=tmp_path, output=tmp_path))
