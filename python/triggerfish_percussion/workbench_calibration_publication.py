"""Publish a verified fit as the actual workbench preset, not a second audition."""

import json

from .fit_provenance import verify_candidate


def publish_workbench_calibration(renderer, directory, source, served):
    """Keep full editable parameters and reference identity in one versioned JSON.

    The source copy survives rebuilds; the served copy updates the existing server.
    Already-open browser state is deliberately not reloaded or overwritten.
    """
    saved = verify_candidate(renderer, directory)
    if saved["metadata"]["recipeKey"] == "drum.kick.v1":
        from .kick_quality_checks import check_kick_candidate

        assessment = check_kick_candidate(directory, renderer=renderer)
        if not assessment["eligible"]:
            raise ValueError(
                "Kick shape/decay checks failed; workbench preset unchanged. "
                "See quality-checks.json; a lower aggregate score is not sufficient."
            )
    fit = json.loads((directory / "candidate.fit.json").read_text(encoding="utf8"))
    instrument = fit.get("instrument", {}).get("name", "Percussion")
    fit["name"] = f"{instrument} - current workbench calibration (review required)"
    payload = json.dumps(fit, indent=2, ensure_ascii=False) + "\n"
    for path in (source, served):
        pending = path.with_suffix(path.suffix + ".pending")
        pending.write_text(payload, encoding="utf8")
        pending.replace(path)
    print(f"Updated workbench calibration: {source}", flush=True)
