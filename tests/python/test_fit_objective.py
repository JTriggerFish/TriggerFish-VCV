"""Saved measurements must not silently change during continuation."""

import json
import numpy as np
import pytest
from types import SimpleNamespace

from triggerfish_percussion.fit_objective import check_objective, saved_metallic_loss
from triggerfish_percussion.metallic_balance_loss import MetallicBalanceLoss
from triggerfish_percussion.instrument_fit_configuration import configure_fit


def test_saved_erb_and_attack_configuration_is_preserved():
    rate = 24000
    reference = np.zeros(rate * 4)
    loss = MetallicBalanceLoss(
        reference, rate, contrast_weighting="erb", fast_attack=True
    )
    saved = json.loads(
        json.dumps(
            dict(
                objective=type(loss).__name__,
                objective_specification=loss.specification,
            )
        )
    )
    restored = saved_metallic_loss(saved, reference, rate)
    check_objective(saved, restored)
    assert restored.specification["contrast_weighting"] == "erb"
    assert restored.specification["fast_attack"]
    with pytest.raises(ValueError, match="objective changed"):
        check_objective(saved, MetallicBalanceLoss(reference, rate))
    saved["objective_specification"]["linear_scale"] = 42
    with pytest.raises(ValueError, match="objective changed"):
        saved_metallic_loss(saved, reference, rate)
    with pytest.raises(ValueError, match="Unsupported"):
        saved_metallic_loss({}, reference, rate)


def test_new_configuration_defaults_to_current_shared_method():
    renderer = SimpleNamespace(
        sample_rate=24000,
        reference=np.random.default_rng(42).normal(0, 0.01, 6 * 24000),
        metadata=dict(reference={}),
    )
    seconds, _, loss = configure_fit(renderer, "ride", None, {})
    assert seconds == 6
    assert loss.specification["contrast_weighting"] == "erb"
    assert loss.specification["fast_attack"]
    saved = dict(
        objective=type(loss).__name__,
        objective_specification=loss.specification,
        duration_seconds=seconds,
    )
    _, _, resumed = configure_fit(renderer, "ride", saved, {})
    check_objective(saved, resumed)
    for environment in (
        {"TF_FIT_FAST_ATTACK": "0"},
        {"TF_FIT_CONTRAST": "linear"},
        {"TF_FIT_OBJECTIVE": "trajectory-v1"},
        {"TF_FIT_SECONDS": "4"},
    ):
        with pytest.raises(ValueError, match="changed"):
            configure_fit(renderer, "ride", saved, environment)
    with pytest.raises(ValueError, match="0 or 1"):
        configure_fit(renderer, "ride", None, {"TF_FIT_FAST_ATTACK": "typo"})
