"""Publishing a fit must not change its recorded comparison baseline."""

import importlib.util
import json
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location(
    "candidate_audit", Path(__file__).parents[2] / "tools/verify_metallic_candidate.py"
)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_baseline_survives_source_publication(tmp_path):
    source = tmp_path / "published.json"
    source.write_text(
        json.dumps({"instrument": {"nodes": [{"parameters": {"field_gain": 0.5}}]}})
    )
    parameters, digest = audit.audit_baseline("ride", tmp_path, source)
    report = dict(target="ride", baseline_parameters=parameters, baseline_sha256=digest)
    (tmp_path / "fresh-seed-audit.json").write_text(json.dumps(report))
    source.write_text("{}")
    assert audit.audit_baseline("ride", tmp_path, source) == (parameters, digest)
    with pytest.raises(ValueError, match="target"):
        audit.audit_baseline("gong", tmp_path, source)
    report["baseline_parameters"]["field_gain"] = 2
    (tmp_path / "fresh-seed-audit.json").write_text(json.dumps(report))
    with pytest.raises(ValueError, match="hash"):
        audit.audit_baseline("ride", tmp_path, source)
