"""Fit observation prominence cheaply, then judge and save actual C++ renders."""

import numpy as np

from .observation_fit_basis import ObservationBasis
from .workbench_search import Search


def polish_observation(search, iterations=20):
    active = [
        key
        for key, value in search.parameters.items()
        if key.startswith("resolved_level_") and value > -71.99
    ]
    # The global observation gain is FIXED in this block, so a common offset
    # of the bars is identifiable. Anchoring an arbitrary bar would forbid a
    # necessary body/contact balance change after changing energy transport.
    keys = active
    if not keys:
        return
    before = float(np.linalg.norm(search.residual(search.parameters)))
    basis = ObservationBasis(
        search.renderer, search.parameters, keys, search.seconds, search.seeds
    )
    directory = search.output / "observation-polish"
    directory.mkdir(parents=True, exist_ok=True)
    trial = Search(
        basis, search.loss, directory, search.seconds, search.name, search.seeds
    )
    trial.stage(
        "validated affine observation amplitudes",
        {key: (-45, 6) for key in keys},
        iterations,
        parameter_scales={key: "amplitude_db" for key in keys},
        influence_threshold=0,
    )
    # The approximation is only float32 summation/subtraction error. Nevertheless,
    # only the exact engine may decide selection and supply published audio.
    after = float(np.linalg.norm(search.residual(trial.parameters)))
    selected = after < before
    if selected:
        search.parameters = dict(trial.parameters)
    search.history.append(
        dict(
            stage="exact-render observation polish",
            before=before,
            after=after,
            selected=selected,
            basis_validation=basis.validation,
            trial_history=trial.history,
            parameters=dict(search.parameters),
        )
    )
    search.save()
    trial.renderer = search.renderer
    trial.cache.clear()
    trial.save()
    return search.history[-1]
