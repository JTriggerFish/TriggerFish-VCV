"""Budgeted scalar optimization of the exact renderer, with logged differences."""

import time
import numpy as np
from scipy.optimize import minimize
from .fit_parameter_box import ParameterBox


class EvaluationBudget(Exception):
    pass


class ScalarTrial:
    def __init__(self, search, bounds, budget, step, method="L-BFGS-B"):
        if budget < 1 or not 0 < step < 0.1:
            raise ValueError(
                "Positive budget and small finite-difference step required"
            )
        self.search, self.bounds, self.budget, self.step = search, bounds, budget, step
        if method not in ("L-BFGS-B", "Powell"):
            raise ValueError("Unsupported scalar fitting method")
        self.method = method
        self.box = ParameterBox(
            search.parameters, bounds, search.renderer.metadata["descriptors"]
        )
        self.cache, self.trace, self.influence = {}, [], []
        self.best, self.best_x = float("inf"), self.box.initial.copy()
        self.began = time.perf_counter()

    def score(self, x):
        key = tuple(x)
        if key in self.cache:
            return self.cache[key]
        if len(self.cache) >= self.budget:
            raise EvaluationBudget()
        values = self.box.unpack(x)
        value = float(
            np.mean(
                [
                    self.search.loss.score(self.search.audio(values, seed))
                    for seed in self.search.seeds
                ]
            )
        )
        if not np.isfinite(value):
            raise ValueError("Non-finite fitting loss")
        self.cache[key] = value
        if value < self.best:
            self.best, self.best_x = value, x.copy()
            self.trace.append(
                dict(
                    evaluation=len(self.cache),
                    score=value,
                    elapsed_seconds=time.perf_counter() - self.began,
                )
            )
        if len(self.cache) % 100 == 0:
            print(
                f"{self.search.name}: {len(self.cache)}/{self.budget} evaluations; "
                f"best {self.best:.6g}",
                flush=True,
            )
        return value

    def gradient(self, x):
        slopes = []
        first = not self.influence
        for i, key in enumerate(self.box.keys):
            minus, plus = x.copy(), x.copy()
            minus[i], plus[i] = max(0, x[i] - self.step), min(1, x[i] + self.step)
            a, b = self.score(minus), self.score(plus)
            slopes.append((b - a) / (plus[i] - minus[i]))
            if first:
                self.influence.append(
                    dict(parameter=key, minus=a, plus=b, slope=slopes[-1])
                )
        return np.array(slopes)

    def run(self):
        before = self.score(self.box.initial)
        try:
            options = dict(maxiter=100, ftol=1e-9, gtol=1e-6, maxls=12)
            if self.method == "Powell":
                # Still record local sensitivity; do not use noisy finite
                # differences as derivatives in the bounded line searches.
                self.gradient(self.box.initial)
                options = dict(maxiter=100, ftol=0.0005, xtol=0.005)
            result = minimize(
                self.score,
                self.box.initial,
                method=self.method,
                jac=self.gradient if self.method == "L-BFGS-B" else None,
                bounds=[(0, 1)] * len(self.box.keys),
                options=options,
            )
            status = str(result.message)
        except EvaluationBudget:
            status = "explicit unique-parameter evaluation budget reached"
        return self.save(before, status)

    def save(self, before, status):
        self.search.parameters = self.box.unpack(self.best_x)
        record = dict(
            stage="matched scalar loss trial",
            before=before,
            after=self.best,
            budget=self.budget,
            parameter_evaluations=len(self.cache),
            renders=self.search.evaluations,
            elapsed_seconds=time.perf_counter() - self.began,
            solver=self.method,
            status=status,
            bounds=self.bounds,
            difference_step_fraction=self.step,
            logarithmic_parameters=[
                k for k, log in zip(self.box.keys, self.box.logarithmic) if log
            ],
            initial_parameters=self.box.start,
            influence=self.influence,
            trace=self.trace,
            parameters=self.search.parameters,
            selected=self.best < before,
        )
        self.search.history.append(record)
        self.search.save()
        return record


def refine_scalar(search, bounds, budget=1000, step=0.001, method="L-BFGS-B"):
    """Minimize the native scalar loss, not a rank-one least-squares surrogate."""
    return ScalarTrial(search, bounds, budget, step, method).run()
