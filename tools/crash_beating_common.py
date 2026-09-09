"""Shared structured-crash objective and smooth observation parameterization."""

import numpy as np
from scipy.optimize import lsq_linear

from triggerfish_percussion.coarse_observation_fit import interpolation_weights
from triggerfish_percussion.perceptual_fit_losses import AuralossMel
from triggerfish_percussion.reference_floor_mel import ReferenceFloorMel
from fit_structured_metal_texture import Objective

KNOTS = (120, 400, 1000, 2500, 6500, 15000)


class CrashObjective(Objective):
    """Spectral/attack-led proposal score; texture is not allowed to dominate."""

    units = "Mel + 0.3 attack Mel + 0.05 bloom + 0.3 texture"

    def __init__(self, reference, rate, reference_floor_db=None):
        super().__init__(reference, rate)
        self.attack_frames = round(0.3 * rate)
        self.attack = AuralossMel(reference[: self.attack_frames], rate)
        self.specification.update(
            version="crash-low-blur-v1",
            weights=[1, 0.3, 0.05, 0.3],
            order=["mel", "attack_mel", "bloom", "texture"],
            attack=self.attack.specification,
            attack_seconds=0.3,
        )
        if reference_floor_db is not None:
            self.mel = ReferenceFloorMel(reference, rate, reference_floor_db)
            self.attack = ReferenceFloorMel(
                reference[: self.attack_frames], rate, reference_floor_db
            )
            self.specification.update(
                version="crash-reference-floor-v2",
                mel=self.mel.specification,
                attack=self.attack.specification,
            )

    def components(self, audio):
        return dict(
            super().components(audio),
            attack_mel=self.attack.score(audio[: self.attack_frames]),
        )

    def score(self, audio):
        return self.score_components(self.components(audio))

    @staticmethod
    def score_components(c):
        """Keep batch/seed audits on exactly the same declared weighting."""
        return c["mel"] + 0.3 * c["attack_mel"] + 0.05 * c["bloom"] + 0.3 * c["texture"]


def project_levels(parameters, source):
    """Project one smooth observation curve, not individual modal amplitudes."""
    keys = [i for i in range(32) if source[f"resolved_level_{i}"] > -71.99]
    weights = interpolation_weights(
        [source[f"resolved_frequency_{i}"] for i in keys], KNOTS
    )
    amplitudes = 10 ** (np.array([source[f"resolved_level_{i}"] for i in keys]) / 20)
    curve = lsq_linear(weights, amplitudes, bounds=(10 ** (-45 / 20), 10 ** (6 / 20))).x
    active = [i for i in range(32) if parameters[f"resolved_level_{i}"] > -71.99]
    projected = (
        interpolation_weights(
            [parameters[f"resolved_frequency_{i}"] for i in active], KNOTS
        )
        @ curve
    )
    return dict(
        parameters,
        **{
            f"resolved_level_{i}": float(20 * np.log10(a))
            for i, a in zip(active, projected)
        },
    )
