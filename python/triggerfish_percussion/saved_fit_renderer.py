"""Render edits to a complete user snapshot, preserving its gesture and routing."""

from copy import deepcopy
from datetime import datetime, timezone
from uuid import uuid4


class SavedFitRenderer:
    def __init__(self, renderer, fit):
        self.renderer = renderer
        self.fit = deepcopy(fit)
        self.sample_rate = renderer.sample_rate
        self.initial = self.parameters(fit)
        reference = renderer.metadata["reference"]
        for key in ("sha256", "referenceGainDb", "sampleRate"):
            if fit["reference"].get(key) != reference.get(key):
                raise ValueError(f"Snapshot reference differs from renderer: {key}")
        onset = lambda item: (item.get("cell") or {}).get("onset_seconds", 0)
        if onset(fit["reference"]) != onset(reference):
            raise ValueError("Snapshot reference onset differs from renderer")
        if fit["instrument"]["recipe"] != renderer.metadata["recipeKey"]:
            raise ValueError("Snapshot recipe differs from renderer")
        if set(self.initial) != set(renderer.initial):
            raise ValueError("Snapshot controls differ from current renderer")

    @staticmethod
    def parameters(fit):
        result = {}
        for node in fit["instrument"]["nodes"]:
            for key, value in node["parameters"].items():
                if key in result:
                    raise ValueError(f"Duplicate parameter: {key}")
                result[key] = value
        return result

    def snapshot(self, parameters, name=None):
        if set(parameters) != set(self.initial):
            raise ValueError("Parameter surface changed")
        fit = deepcopy(self.fit)
        for node in fit["instrument"]["nodes"]:
            node["parameters"] = {key: parameters[key] for key in node["parameters"]}
        if name is not None:
            fit.update(
                id=str(uuid4()),
                parentId=self.fit["id"],
                name=name,
                createdAt=datetime.now(timezone.utc).isoformat(),
            )
        return fit

    def render(self, parameters, seconds, seed=None, event=None):
        fit = self.snapshot(parameters)
        fit["controls"]["event"].update(event or {})
        if seed is not None:
            fit["controls"]["event"]["seed"] = seed
        response = self.renderer.request(
            command="renderSnapshot", fit=fit, seconds=seconds
        )
        return self.renderer.decode(response["pcm"])
