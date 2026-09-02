"""Compatibility facade for the staged crash calibration toolkit."""

from .crash_fit_acceptance import acceptance_diagnostics
from .crash_fit_causal import (
    CAUSAL_FIT_POLICIES,
    FIRST_100MS_TRADEOFF_POLICY,
    STRICT_CAUSAL_POLICY,
    CausalFitPolicy,
    causal_prefix_losses,
    fit_causal_model,
    prefix_gate_passes,
)
from .crash_fit_common import CrashFitCell
from .crash_fit_features import (
    FeatureSet,
    REGIONS,
    TRAJECTORY_TIMES_SECONDS,
    maximum_late_regrowth_db as _maximum_late_regrowth_db,
    perceptual_features,
)
from .crash_fit_influence import parameter_influence_diagnostics
from .crash_fit_prefix import (
    CAUSAL_PREFIX_SECONDS,
    CausalFeatures,
    causal_audio_quality,
    causal_feature_loss,
    causal_prefix_features,
    causal_prefix_quality,
)
from .crash_fit_modes import (
    fit_sparse_modes,
    fit_sparse_projection,
    modal_passes as _modal_passes,
    modal_tail as _modal_tail,
    reference_modal_residual,
    stable_modes as _stable_modes,
)

__all__ = [
    "CrashFitCell",
    "FeatureSet",
    "REGIONS",
    "TRAJECTORY_TIMES_SECONDS",
    "CAUSAL_PREFIX_SECONDS",
    "CausalFeatures",
    "CausalFitPolicy",
    "CAUSAL_FIT_POLICIES",
    "STRICT_CAUSAL_POLICY",
    "FIRST_100MS_TRADEOFF_POLICY",
    "acceptance_diagnostics",
    "causal_feature_loss",
    "causal_audio_quality",
    "causal_prefix_features",
    "causal_prefix_quality",
    "causal_prefix_losses",
    "fit_causal_model",
    "fit_sparse_modes",
    "fit_sparse_projection",
    "perceptual_features",
    "parameter_influence_diagnostics",
    "prefix_gate_passes",
    "reference_modal_residual",
]
