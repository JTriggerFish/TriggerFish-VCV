"""Bounded population search; use fresh-seed validation before selecting a patch."""

import json

import numpy as np
from scipy.optimize import differential_evolution


def population_stage(search, bounds, generations=18):
    keys = list(bounds)
    low, high = np.asarray(list(bounds.values()), dtype=float).T
    descriptors = {
        item["key"]: item for item in search.renderer.metadata["descriptors"]
    }
    low = np.maximum(low, [descriptors[key]["minimum"] for key in keys])
    high = np.minimum(high, [descriptors[key]["maximum"] for key in keys])
    if np.any(high <= low):
        raise ValueError("Search bounds must overlap the UI range with positive width")
    initial = dict(search.parameters)
    start = np.clip(
        (np.array([initial[key] for key in keys]) - low) / (high - low), 0, 1
    )
    random = np.random.default_rng(812)
    population = random.uniform(0, 1, (max(5, 4 * len(keys)), len(keys)))
    population[: len(keys)] = np.clip(
        start + random.normal(0, 0.12, (len(keys), len(keys))), 0, 1
    )
    population[0] = start
    trials = {}
    generation = 0

    def unpack(x):
        return dict(initial, **dict(zip(keys, (low + x * (high - low)).tolist())))

    def objective(x):
        # One common realization makes exploration affordable. Selection below
        # re-ranks the leading candidates on BOTH training seeds.
        score = float(np.linalg.norm(search.loss.residual(search.audio(unpack(x)))))
        trials[tuple(x)] = score
        return score

    def progress(x, convergence):
        nonlocal generation
        generation += 1
        print(
            json.dumps(
                dict(
                    global_generation=generation,
                    best=min(trials.values()),
                    renderer_calls=search.evaluations,
                )
            ),
            flush=True,
        )

    before = float(np.linalg.norm(search.residual(initial)))
    differential_evolution(
        objective,
        [(0, 1)] * len(keys),
        init=population,
        maxiter=generations,
        rng=random,
        mutation=(0.4, 0.9),
        recombination=0.8,
        polish=False,
        tol=0.001,
        callback=progress,
    )
    best = before
    finalists = []
    for point, exploration_score in sorted(trials.items(), key=lambda item: item[1])[
        :12
    ]:
        parameters = unpack(np.array(point))
        score = float(np.linalg.norm(search.residual(parameters)))
        finalists.append(
            dict(
                exploration_error=exploration_score,
                training_error=score,
                parameters=parameters,
            )
        )
        if score < best:
            best, search.parameters = score, parameters
    search.history.append(
        dict(
            stage="population exploration with two-seed selection",
            before=before,
            after=best,
            evaluations=search.evaluations,
            objective_specification=search.loss.specification,
            finalists=finalists,
            parameters=dict(search.parameters),
        )
    )
    search.save()
    print(json.dumps(dict(global_before=before, global_after=best)), flush=True)
