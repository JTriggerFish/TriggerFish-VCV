"""Only resume explicit-mode kick fits; no automatic legacy/EQ migration."""

import json
from pathlib import Path
from triggerfish_percussion.workbench_fit_baseline import check_reference


def load_start(search, path):
    saved = json.loads(Path(path).read_text())
    check_reference(saved["metadata"], search.renderer.metadata)
    values = saved["parameters"]
    if set(values) != set(search.parameters):
        raise ValueError(
            "This fit predates explicit kick modes; generate an editable modal start"
        )
    if values["equalizer_mode"] == 2:
        raise ValueError("Multiband-EQ fitting candidates are not kick fitting starts")
    search.parameters.update(values)
    search.history = saved["history"]
