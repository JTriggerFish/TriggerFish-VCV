from types import SimpleNamespace

import numpy as np
import pytest

from triggerfish_percussion.observation_fit_basis import ObservationBasis


def renderer(nonlinear=False):
    def render(parameters, seconds, seed=None):
        gain = 10 ** (parameters["resolved_level_0"] / 20)
        return np.array([0.1, 0.2, 0.3]) + gain ** (2 if nonlinear else 1) * np.array(
            [0.4, -0.5, 0.6]
        )

    return SimpleNamespace(metadata={}, sample_rate=48000, render=render)


def test_polish_saves_actual_engine_audio(tmp_path):
    from triggerfish_percussion.observation_fit_polish import polish_observation
    from triggerfish_percussion.workbench_search import Search
    from triggerfish_percussion.audio_io import read_wav

    real = renderer()
    real.initial = dict(resolved_level_0=-6, resolved_level_1=0)
    real.metadata = {
        "descriptors": [dict(key=k, minimum=-72, maximum=6) for k in real.initial]
    }
    real.request = lambda **kwargs: {"fit": {"parameters": kwargs["parameters"]}}
    desired = real.render(dict(real.initial, resolved_level_0=-2), 1, None)
    loss = SimpleNamespace(
        residual=lambda samples, regions=range(5): 100 * (samples - desired),
        diagnostics=lambda samples: {},
    )
    search = Search(real, loss, tmp_path, seconds=1, seeds=(None,))
    polish_observation(search, 8)
    assert abs(search.parameters["resolved_level_0"] + 2) < 0.05
    actual = real.render(search.parameters, 1, None).astype(np.float32).astype(float)
    assert np.array_equal(read_wav(tmp_path / "candidate.wav").mono().samples, actual)


def test_basis_exact_and_only_observation_changes():
    initial = dict(resolved_level_0=-6, bloom_rate=1)
    real = renderer()
    basis = ObservationBasis(real, initial, ["resolved_level_0"], 1, [None])
    values = dict(initial, resolved_level_0=2)
    assert np.allclose(
        basis.render(values, 1), real.render(values, 1, None), atol=1e-12
    )
    with pytest.raises(ValueError, match="non-observation"):
        basis.render(dict(values, bloom_rate=2), 1)
    with pytest.raises(ValueError, match="active-mode"):
        basis.render(dict(values, resolved_level_0=-72), 1)
    with pytest.raises(ValueError, match="duration or seed"):
        basis.render(values, 2)


def test_rejects_nonlinear_observation():
    with pytest.raises(ValueError, match="not affine"):
        ObservationBasis(
            renderer(True), dict(resolved_level_0=-6), ["resolved_level_0"], 1, [None]
        )


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
@pytest.mark.parametrize("source", ["actual", "predicted"])
def test_validation_rejects_nonfinite_audio(bad, source):
    real = renderer()
    basis = ObservationBasis(
        real, dict(resolved_level_0=-6), ["resolved_level_0"], 1, [None]
    )
    invalid = lambda *args: np.full(3, bad)
    if source == "actual":
        real.render = invalid
    else:
        basis.render = invalid
    with pytest.raises(ValueError, match="finite audio"):
        basis.validate(None)
