"""The temporal objective must not hide delayed layers or wrong audition gain."""

import numpy as np
import pytest

from triggerfish_percussion.layered_band_loss import LayeredBandLoss


def example():
    rate = 32000
    t = np.arange(6 * rate) / rate
    low = 0.1 * np.exp(-t / 0.8) * np.sin(2 * np.pi * 360 * t)
    upper = 0.03 * np.maximum(t - 0.3, 0) * np.exp(-t / 0.3)
    upper *= np.sin(2 * np.pi * 10000 * t)
    return rate, low, upper


def test_reference_is_zero_and_gain_is_not_normalized():
    rate, low, upper = example()
    loss = LayeredBandLoss(low + upper, rate)
    assert loss.score(low + upper) < 1e-10
    assert loss.score(2 * (low + upper)) > 1


def test_layer_delay_and_shape_survive_level_independent_screen():
    rate, low, upper = example()
    loss = LayeredBandLoss(low + upper, rate)
    delay = rate // 3
    late = np.r_[np.zeros(delay), upper[:-delay]]
    assert loss.score_db(loss.envelopes(low + late), True) > 0.2


def test_invalid_audio_is_rejected():
    rate, low, upper = example()
    loss = LayeredBandLoss(low + upper, rate)
    with pytest.raises(ValueError, match="finite"):
        loss.envelopes(np.full_like(low, np.nan))
    with pytest.raises(ValueError, match="length"):
        loss.envelopes(low[:100])
    with pytest.raises(ValueError, match="six seconds"):
        LayeredBandLoss(low[:100], rate)
    with pytest.raises(ValueError, match="finite"):
        loss.score_db(np.full_like(loss.target, np.nan))


def test_reference_weighting_preserves_gain_and_delay_checks():
    rate, low, upper = example()
    loss = LayeredBandLoss(low + upper, rate, audibility=True)
    assert loss.score(low + upper) < 1e-10
    assert loss.score(2 * (low + upper)) > 1
    assert np.allclose(loss.weights.mean(axis=1), 1)
    frozen = loss.weights.copy()
    loss.score(0.01 * low + 100 * upper)
    assert np.array_equal(loss.weights, frozen)


@pytest.mark.parametrize("weighted", [False, True])
def test_attribution_exactly_decomposes_objective(weighted):
    rate, low, upper = example()
    loss = LayeredBandLoss(low + upper, rate, audibility=weighted)
    audio = low * 0.1 + upper * 3
    audit = loss.attribution(audio)
    cost = np.asarray(audit["cell_cost"])
    assert cost.sum() == pytest.approx(loss.score(audio) ** 2)
    assert cost.sum(axis=1) == pytest.approx(audit["band_cost"])
    assert cost.sum(axis=0) == pytest.approx(audit["region_cost"])
