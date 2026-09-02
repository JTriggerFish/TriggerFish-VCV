import json

import numpy as np

from triggerfish_percussion.capture_quality import _spearman, qualify_cells


def test_spearman_detects_order_and_reversal():
    values = np.asarray([1.0, 2.0, 3.0, 4.0])
    assert _spearman(values, values) == 1.0
    assert _spearman(values, values[::-1]) == -1.0


def test_empty_capture_grid_is_rejected_without_crashing(tmp_path):
    manifest = tmp_path / "cells.json"
    manifest.write_text(json.dumps({"cells": []}), encoding="utf-8")

    result = qualify_cells(manifest)

    assert not result["accepted"]
    assert not result["checks"]["cell_count"]
    assert result["cells"] == []


def test_spearman_rejects_an_unidentifiable_single_layer():
    assert _spearman(np.asarray([64]), np.asarray([0.5])) == 0.0
