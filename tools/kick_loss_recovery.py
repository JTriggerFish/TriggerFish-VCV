"""Known C++ contact-noise recovery; these are not real-reference fit files."""

from contextlib import ExitStack
import json
from pathlib import Path

from triggerfish_percussion.scalar_fit_search import refine_scalar
from triggerfish_percussion.workbench_search import Search


def run_recovery(renderer, make_losses, output, remote):
    known = dict(
        renderer.initial,
        resonance_level=0,
        thump_level=0,
        contact_level=0.5,
        contact_noise_level=1.2,
        contact_noise_decay_seconds=0.18,
        contact_colour=0.7,
    )
    seed = renderer.metadata["event"]["seed"]
    reference = renderer.render(known, 1.2, seed)
    bounds = dict(contact_noise_level=(0.1, 3), contact_noise_decay_seconds=(0.03, 0.4))
    rows = []
    with ExitStack() as stack:
        losses = make_losses(reference, renderer.sample_rate, remote, stack)
        for name, loss in losses.items():
            search = Search(renderer, loss, output, 1.2, f"recovery-{name}", (seed,))
            # Do not emit a fit falsely identifying this synthetic reference as
            # the acoustic sample. This experiment saves its provenance below.
            search.save = lambda: None
            search.parameters = dict(
                known, contact_noise_level=0.5, contact_noise_decay_seconds=0.07
            )
            row = refine_scalar(search, bounds, 250)
            recovered = search.parameters
            passed = (
                abs(recovered["contact_noise_level"] - 1.2) < 0.03
                and abs(recovered["contact_noise_decay_seconds"] - 0.18) < 0.003
            )
            rows.append(
                dict(
                    objective=name,
                    fit=row,
                    passed=passed,
                    specification=loss.specification,
                )
            )
            print(
                json.dumps(
                    dict(
                        recovery=name,
                        passed=passed,
                        level=recovered["contact_noise_level"],
                        decay=recovered["contact_noise_decay_seconds"],
                    )
                ),
                flush=True,
            )
    report = dict(
        reference_kind="synthetic exact C++/Wasm contact",
        known=known,
        renderer_sha256=renderer.metadata["rendererSha256"],
        seed=seed,
        rows=rows,
        all_passed=all(r["passed"] for r in rows),
    )
    Path(output).mkdir(parents=True, exist_ok=True)
    (Path(output) / "synthetic-recovery.json").write_text(
        json.dumps(report, indent=2), encoding="utf8"
    )
