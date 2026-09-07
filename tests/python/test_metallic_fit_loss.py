"""Known-answer checks before applying the objective to instrument fitting."""

import numpy as np
import pytest

from triggerfish_percussion.metallic_fit_loss import MetallicFitLoss
from triggerfish_percussion.instrument_fit_stages import metallic_stages, snare_stages


@pytest.fixture(scope="module")
def example():
    rate = 24000
    t = np.arange(rate * 6) / rate
    rng = np.random.default_rng(721)
    x = (
        0.15 * np.sin(2 * np.pi * 173 * t) * np.exp(-t / 1.2)
        + 0.03 * np.sin(2 * np.pi * 4300 * t) * np.exp(-t / 0.6)
        + 0.06 * rng.normal(size=len(t)) * np.exp(-t / 0.7)
    )
    return x, rate, MetallicFitLoss(x, rate)


def test_identity_and_fixed_level(example):
    x, rate, loss = example
    assert np.linalg.norm(loss.residual(x)) == 0
    quiet = np.linalg.norm(loss.residual(x * 0.5))
    assert quiet > 4
    assert np.linalg.norm(loss.residual(x * 0.25)) > quiet


def test_upper_ring_and_attack_are_not_hidden(example):
    x, rate, loss = example
    t = np.arange(len(x)) / rate
    ring = 0.08 * np.sin(2 * np.pi * 7200 * t) * np.exp(-t / 2)
    assert np.linalg.norm(loss.residual(x + ring)) > 1
    altered = x.copy()
    altered[: round(0.08 * rate)] *= 0.1
    assert np.linalg.norm(loss.residual(altered, (0,))) > 2


def test_decay_direction_known_answer(example):
    x, rate, loss = example
    t = np.arange(len(x)) / rate
    errors = [np.linalg.norm(loss.residual(x * np.exp(-k * t))) for k in (0, 0.2, 1)]
    assert errors[0] < errors[1] < errors[2]


def test_no_per_mode_decay_or_event_fitting():
    parameters = {f"resolved_level_{i}": -6 for i in range(32)}
    parameters.update({f"resolved_frequency_{i}": 100 * (i + 1) for i in range(32)})
    for stages in (metallic_stages(parameters), snare_stages(parameters)):
        keys = {key for _, bounds, _ in stages for key in bounds}
        assert not keys.intersection(
            {"strength", "location", "hardness", "implement", "constraint"}
        )
        assert not any(
            key.startswith("body_decay_")
            and key not in {"body_decay_seconds_0", "body_decay_seconds_7"}
            for key in keys
        )
