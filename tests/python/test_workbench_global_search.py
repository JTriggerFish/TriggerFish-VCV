"""A primary-seed improvement must not worsen the multi-seed objective."""

from types import SimpleNamespace

import numpy as np

from triggerfish_percussion.workbench_global_search import population_stage
from triggerfish_percussion.workbench_search import Search


def test_population_selection_preserves_two_seed_optimum(tmp_path, monkeypatch):
    renderer = SimpleNamespace(
        initial={"gain": -0.4},
        metadata={"descriptors": [{"key": "gain", "minimum": -1, "maximum": 1}]},
        render=lambda parameters, seconds, seed=None: np.array(
            [parameters["gain"] + (0.8 if seed == 1 else 0)]
        ),
    )
    loss = SimpleNamespace(
        residual=lambda audio, regions=range(5): audio, specification={}
    )
    search = Search(renderer, loss, tmp_path, seeds=(None, 1))
    monkeypatch.setattr(search, "save", lambda: None)
    population_stage(search, {"gain": (-1, 1)}, generations=3)
    assert search.parameters["gain"] == -0.4
    assert search.history[-1]["after"] == search.history[-1]["before"]
