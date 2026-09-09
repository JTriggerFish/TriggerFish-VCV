"""Bounded, logged searches over the exact workbench renderer."""

import json
from itertools import product

import numpy as np
from scipy.optimize import least_squares

from .audio_io import AudioBuffer, write_wav


class Search:
    def __init__(self, renderer, loss, output, seconds=6, name="Gong", seeds=(None,)):
        self.renderer, self.loss, self.output = renderer, loss, output
        self.parameters = dict(renderer.initial)
        self.evaluations = 0
        self.history = []
        self.cache = {}
        self.seconds, self.name = seconds, name
        self.seeds = seeds

    def audio(self, parameters, seed=None):
        key = (seed, tuple(sorted(parameters.items())))
        if key not in self.cache:
            self.cache[key] = (
                self.renderer.render(parameters, self.seconds)
                if seed is None
                else self.renderer.render(parameters, self.seconds, seed)
            )
            self.evaluations += 1
            if len(self.cache) > 80:
                del self.cache[next(iter(self.cache))]
        return self.cache[key]

    def residual(self, parameters, regions=range(5)):
        """Keep phase realizations separate, never average their audio."""
        return np.concatenate(
            [
                self.loss.residual(self.audio(parameters, seed), regions)
                for seed in self.seeds
            ]
        ) / np.sqrt(len(self.seeds))

    def stage(
        self,
        name,
        bounds,
        iterations,
        regions=range(5),
        difference_step=0.005,
        parameter_scales=None,
        influence_threshold=None,
    ):
        if not 0 < difference_step <= 0.1:
            raise ValueError("Finite-difference step must be in (0, .1]")
        keys = list(bounds)
        low, high = np.array(list(bounds.values()), dtype=float).T
        descriptors = {
            item["key"]: item for item in self.renderer.metadata["descriptors"]
        }
        low = np.maximum(low, [descriptors[key]["minimum"] for key in keys])
        high = np.minimum(high, [descriptors[key]["maximum"] for key in keys])
        if np.any(high <= low):
            raise ValueError(
                "Search bounds must overlap the UI range with positive width"
            )
        seed = dict(self.parameters)
        effective_bounds = dict(zip(keys, zip(low.tolist(), high.tolist())))
        physical_low, physical_high = low.copy(), high.copy()
        scales = [dict(parameter_scales or {}).get(key, "linear") for key in keys]
        if any(scale not in ("linear", "amplitude_db", "log") for scale in scales):
            raise ValueError("Unknown search parameter scale")
        # Exposed dB observation levels remain dB in every saved patch. Only
        # solver coordinates change: a quiet bar must receive a meaningful
        # amplitude probe rather than being declared dead at its current level.
        amplitude = np.array([scale == "amplitude_db" for scale in scales])
        logarithmic = np.array([scale == "log" for scale in scales])
        raw = np.array([seed[key] for key in keys], dtype=float)
        if np.any(low[logarithmic] <= 0):
            raise ValueError("Log search bounds must be strictly positive")
        raw[logarithmic] = np.log(np.maximum(raw[logarithmic], low[logarithmic]))
        low[logarithmic] = np.log(low[logarithmic])
        high[logarithmic] = np.log(high[logarithmic])
        low[amplitude] = 10 ** (low[amplitude] / 20)
        high[amplitude] = 10 ** (high[amplitude] / 20)
        raw[amplitude] = 10 ** (raw[amplitude] / 20)
        x = np.clip((raw - low) / (high - low), 0, 1)

        def unpack(values):
            physical = low + values * (high - low)
            physical[amplitude] = 20 * np.log10(physical[amplitude])
            physical[logarithmic] = np.exp(physical[logarithmic])
            # exp(log(maximum)) may round one ULP beyond a strict API bound.
            physical = np.clip(physical, physical_low, physical_high)
            return dict(seed, **dict(zip(keys, physical.tolist())))

        def residual(values):
            return self.residual(unpack(values), regions)

        influence = []
        # Compare with the ACTUAL previous patch, not a silently clamped start.
        # A narrower search box must never authorize a regression from that patch.
        first = self.residual(seed, regions)
        for index, key in enumerate(keys):
            minus, plus = x.copy(), x.copy()
            minus[index] = max(0, x[index] - 0.02)
            plus[index] = min(1, x[index] + 0.02)
            influence.append(
                dict(
                    parameter=key,
                    residual_change=float(
                        np.linalg.norm(residual(plus) - residual(minus))
                    ),
                    step=(unpack(plus)[key] - unpack(minus)[key]) / 2,
                )
            )
        print(
            json.dumps(
                dict(
                    stage=name,
                    baseline=float(np.linalg.norm(first)),
                    influence=influence,
                )
            ),
            flush=True,
        )
        active = [
            index
            for index, item in enumerate(influence)
            if item["residual_change"]
            >= (
                getattr(self.loss, "influence_threshold", 0.05)
                if influence_threshold is None
                else influence_threshold
            )
        ]
        if not active:
            return
        # Freeze effectively dead directions instead of allowing arbitrary
        # boundary values to appear in an otherwise equivalent fit.
        keys = [keys[index] for index in active]
        low, high, x = low[active], high[active], x[active]
        amplitude = amplitude[active]
        logarithmic = logarithmic[active]
        physical_low, physical_high = physical_low[active], physical_high[active]
        iteration = 0

        def jacobian(values):
            nonlocal iteration
            iteration += 1
            progress = dict(
                stage=name,
                iteration=iteration,
                score=float(np.linalg.norm(residual(values))),
                parameters=unpack(values),
                metadata=self.renderer.metadata,
                duration_seconds=self.seconds,
                training_seeds=self.seeds,
                status="iteration-checkpoint-not-reviewed",
            )
            self.output.mkdir(parents=True, exist_ok=True)
            pending = self.output / "stage-progress.pending.json"
            pending.write_text(json.dumps(progress, indent=2), encoding="utf8")
            pending.replace(self.output / "stage-progress.json")
            print(
                json.dumps({k: progress[k] for k in ("stage", "iteration", "score")}),
                flush=True,
            )
            columns = []
            for index in range(len(keys)):
                minus, plus = values.copy(), values.copy()
                minus[index] = max(0, values[index] - difference_step)
                plus[index] = min(1, values[index] + difference_step)
                columns.append(
                    (residual(plus) - residual(minus)) / (plus[index] - minus[index])
                )
            return np.array(columns).T

        # TRF's initial radius depends on ||x0||. Offset amplitude coordinates
        # so a bank starting near silence doesn't get a near-zero trust radius.
        offset = 1 if amplitude.any() else 0
        result = least_squares(
            lambda values: residual(values - offset),
            x + offset,
            jac=lambda values: jacobian(values - offset),
            bounds=(offset, 1 + offset),
            max_nfev=iterations,
            ftol=1e-7 if amplitude.any() else 0.002,
            xtol=1e-6 if amplitude.any() else 0.002,
            gtol=0.002,
        )
        candidate = unpack(result.x - offset)
        last = self.residual(candidate, regions)
        # Only select a search step here; this is never a listening approval.
        improved = np.linalg.norm(last) < np.linalg.norm(first)
        if improved:
            self.parameters = candidate
        record = dict(
            stage=name,
            renderer_sha256=self.renderer.metadata.get("rendererSha256"),
            objective_specification=getattr(self.loss, "specification", None),
            objective_units=getattr(self.loss, "units", "dB"),
            jacobian_step_fraction=difference_step,
            parameter_scales=dict(zip(list(bounds), scales)),
            bounds=effective_bounds,
            active_parameters=list(keys),
            fixed_parameters={
                key: value for key, value in seed.items() if key not in keys
            },
            max_nfev=iterations,
            regions=list(regions),
            solver_status=int(result.status),
            solver_message=result.message,
            before=float(np.linalg.norm(first)),
            after=float(np.linalg.norm(last)),
            selected=bool(improved),
            evaluations=self.evaluations,
            influence=influence,
            diagnostics=self.loss.diagnostics(self.audio(self.parameters)),
            parameters=dict(self.parameters),
        )
        self.history.append(record)
        self.save()
        print(
            json.dumps(
                {key: value for key, value in record.items() if key != "influence"}
            ),
            flush=True,
        )

    def screen_candidates(self, name, candidates):
        """Compare discrete starts with the same exact loss and seeds as refinement."""
        before = best = float(np.linalg.norm(self.residual(self.parameters)))
        scores = []
        for label, values in candidates:
            score = float(np.linalg.norm(self.residual(values)))
            scores.append(dict(name=label, score=score, parameters=values))
            if score < best:
                best, self.parameters = score, dict(values)
        self.history.append(
            dict(
                stage=name,
                before=before,
                after=best,
                selected=best < before,
                candidates=scores,
                parameters=dict(self.parameters),
                evaluations=self.evaluations,
            )
        )
        self.save()
        print(
            json.dumps(
                dict(
                    stage=name,
                    before=before,
                    after=best,
                    candidates=[dict(name=s["name"], score=s["score"]) for s in scores],
                )
            ),
            flush=True,
        )

    def explore_texture(self):
        # Cross the low-frequency zero-turbulence clamp deliberately. A local
        # Jacobian cannot discover a cleaner/diffuse regime on its other side.
        initial = dict(self.parameters)
        before = best = float(np.linalg.norm(self.loss.residual(self.audio(initial))))
        for index, (centre, slope, low_scale) in enumerate(
            product((400, 700, 1000, 1400), (0.15, 0.35, 0.6), (0.12, 0.5, 1))
        ):
            candidate = dict(
                initial,
                field_turbulence=initial["field_turbulence"] * (1000 / centre) ** slope,
                field_turbulence_slope=slope,
            )
            candidate.update({f"resolved_turbulence_{i}": low_scale for i in range(9)})
            score = float(np.linalg.norm(self.loss.residual(self.audio(candidate))))
            if score < best:
                best, self.parameters = score, candidate
            if index % 8 == 0:
                print(json.dumps(dict(exploration=index, best=best)), flush=True)
        self.history.append(
            dict(
                stage="low-band texture regimes",
                before=before,
                after=best,
                selected=best < before,
                evaluations=self.evaluations,
                diagnostics=self.loss.diagnostics(self.audio(self.parameters)),
                parameters=dict(self.parameters),
            )
        )
        self.save()
        print(json.dumps(self.history[-1]), flush=True)

    def save(self):
        audio = self.audio(self.parameters)
        write_wav(
            self.output / "candidate.wav", AudioBuffer(audio, self.renderer.sample_rate)
        )
        result = dict(
            parameters=self.parameters,
            history=self.history,
            metadata=self.renderer.metadata,
            listening_approved=False,
            training_seeds=self.seeds,
            objective=type(self.loss).__name__,
            objective_units=getattr(self.loss, "units", "dB"),
            objective_specification=getattr(self.loss, "specification", None),
            duration_seconds=self.seconds,
        )
        (self.output / "search.json").write_text(
            json.dumps(result, indent=2), encoding="utf8"
        )
        fit = self.renderer.request(
            command="snapshot",
            parameters=self.parameters,
            name=f"{self.name} fitting candidate — review required",
        )["fit"]
        (self.output / "candidate.fit.json").write_text(
            json.dumps(fit, indent=2), encoding="utf8"
        )
        if self.history:
            stage = f"stage-{len(self.history):02d}"
            (self.output / f"{stage}.fit.json").write_text(
                json.dumps(fit, indent=2), encoding="utf8"
            )
