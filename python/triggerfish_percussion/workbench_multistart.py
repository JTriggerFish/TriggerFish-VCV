"""Short refits before rejecting alternative modal layouts by their raw score."""

import numpy as np
from .workbench_search import Search


def refine_candidate_starts(search, candidates, bounds_for, *, count=3, iterations=10):
    before = best = float(np.linalg.norm(search.residual(search.parameters)))
    ranked = sorted(
        (float(np.linalg.norm(search.residual(values))), index, name, values)
        for index, (name, values) in enumerate(candidates)
    )
    records = []
    for score, index, name, values in ranked[:count]:
        directory = search.output / f"modal-restart-{index:02d}"
        directory.mkdir(parents=True, exist_ok=True)
        trial = Search(
            search.renderer,
            search.loss,
            directory,
            seconds=search.seconds,
            name=name,
            seeds=search.seeds,
        )
        trial.parameters = dict(values)
        trial.stage("short refit of modal alternative", bounds_for(values), iterations)
        after = float(np.linalg.norm(trial.residual(trial.parameters)))
        search.evaluations += trial.evaluations
        records.append(
            dict(
                name=name,
                raw_score=score,
                refined_score=after,
                parameters=trial.parameters,
                history=str(directory),
            )
        )
        if after < best:
            best, search.parameters = after, dict(trial.parameters)
    search.history.append(
        dict(
            stage="refitted modal alternatives",
            before=before,
            after=best,
            selected=best < before,
            trials=records,
            parameters=dict(search.parameters),
            evaluations=search.evaluations,
        )
    )
    search.save()
    return records
