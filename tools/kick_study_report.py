"""One reference/candidate pair per ablation, plus separate playing checks."""

import json

from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from workbench_fit_report import write_report
from kick_playability import check_playability


def finish_study(renderer, reference, root):
    study = json.loads((root / "study.json").read_text())
    baseline = renderer.render(renderer.initial, 1.2)
    links = []
    for result in study["results"]:
        output = root / result["name"]
        for name, samples in (("reference", reference), ("baseline", baseline)):
            write_wav(
                output / f"{name}.wav", AudioBuffer(samples, renderer.sample_rate)
            )
        write_report(output, kind="kick", assets="../..", renderer=renderer)
        index = output / "index.html"
        html = index.read_text(encoding="utf8")
        label = result["name"].replace("-", " + ")
        html = html.replace(
            "Acoustic kick — reference vs candidate",
            f"{label} — reference vs candidate",
        )
        html += '<p><a href="../">Back to architecture study</a></p>'
        index.write_text(html, encoding="utf8")
        error = result["fitted_diagnostics"]["rms_error_db"]
        links.append(
            f'<li><a href="{result["name"]}/">{label}</a>: {error:.2f} dB engineering error</li>'
        )
    selected = study["results"][0]
    check_playability(renderer, selected["parameters"], root)
    html = '<!doctype html><meta charset="utf-8"><title>Kick architecture study</title>'
    html += "<style>body{background:#10141b;color:#ddd;font:17px system-ui;max-width:900px;margin:40px auto}a{color:#fcbf49}li{margin:16px}</style>"
    html += "<h1>Kick: contact + thump + adjustable resonance</h1><p>Same reference hit and fixed source gain. Each reduced structure was re-fitted after disabling its missing source. These are bounded local searches, not proof of a globally optimal fit.</p><ul>"
    html += (
        "".join(links)
        + "</ul><p>Each link shows only reference versus the named candidate. No factory presets were changed.</p>"
    )
    (root / "index.html").write_text(html, encoding="utf8")
