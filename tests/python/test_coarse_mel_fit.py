from types import SimpleNamespace

import numpy as np
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("auraloss")

from triggerfish_percussion.coarse_mel_fit import polish_coarse_mel
from triggerfish_percussion.coarse_observation_fit import interpolation_weights
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.spectral_bloom_loss import SpectralBloomLoss


@pytest.mark.parametrize("exact_start", [False, True])
def test_six_coordinate_mel_fit_preserves_modes_and_improves_actual_audio(exact_start):
    torch.set_num_threads(1)
    frequencies = np.array([100, 200, 400, 800, 1400, 2000, 2600, 3000])
    knots = (100, 300, 600, 1200, 2200, 3000)
    t = np.arange(48000) / 8000
    signals = np.array(
        [
            np.sin(2 * np.pi * f * t) * np.exp(-t * (1 + i / 8))
            for i, f in enumerate(frequencies)
        ]
    )
    background = 0.002 * np.random.default_rng(81).normal(size=len(t)) * np.exp(-t / 4)

    class Renderer:
        def render(self, parameters, seconds, seed):
            weights = 10 ** (
                np.array([parameters[f"resolved_level_{i}"] for i in range(8)]) / 20
            )
            return background + weights @ signals

    renderer = Renderer()
    parameters = {f"resolved_level_{i}": -72 for i in range(32)}
    for i, f in enumerate(frequencies):
        parameters.update(
            {f"resolved_frequency_{i}": float(f), f"resolved_level_{i}": -20}
        )
    target = (
        background
        + (
            interpolation_weights(frequencies, knots)
            @ [0.07, 0.1, 0.14, 0.18, 0.16, 0.12]
        )
        @ signals
    )
    if exact_start:
        # A valid coarse curve below the optimizer's -45 dB lower bound.
        # Improving on the clipped start is not improving on the real start.
        for i in range(8):
            parameters[f"resolved_level_{i}"] = -60
        target = renderer.render(parameters, 6, 19)
    mel = AuralossMel(target, 8000)
    search = SimpleNamespace(
        renderer=renderer,
        parameters=parameters,
        seconds=6,
        seeds=(19,),
        loss=SpectralBloomLoss(target, 8000),
        history=[],
        save=lambda: None,
        audio=lambda p, seed: renderer.render(p, 6, seed),
    )
    before = mel.score(renderer.render(parameters, 6, 19))
    polish_coarse_mel(search, mel, knots)
    after = mel.score(renderer.render(search.parameters, 6, 19))
    if exact_start:
        assert after == before
        assert search.parameters == parameters
        assert not search.history[-1]["selected"]
    else:
        assert after < before
    assert search.history[-1]["stage"] == "6-coordinate Mel polish"
    assert {k: v for k, v in search.parameters.items() if "frequency" in k} == {
        k: v for k, v in parameters.items() if "frequency" in k
    }
    assert all(search.parameters[f"resolved_level_{i}"] == -72 for i in range(8, 32))
