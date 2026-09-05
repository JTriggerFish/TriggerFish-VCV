"""Pin report links to immutable content, independent of the next search checkpoint."""

from hashlib import sha256
import re


def pin_report_assets(directory, html):
    """Call after current-renderer verification, with one writer per fit directory."""

    def pin(name):
        source = directory / name
        data = source.read_bytes()
        digest = sha256(data).hexdigest()[:20]
        target = f"published-{digest}-{name}"
        destination = directory / target
        if not destination.exists():
            destination.write_bytes(data)
        return target

    fit_name = pin("candidate.fit.json")
    html = html.replace('href="candidate.fit.json"', f'href="{fit_name}"')

    def audio(match):
        name = match.group(1)
        if not re.fullmatch(r"[a-zA-Z0-9_.-]+", name):
            raise ValueError("Invalid report audio name")
        return f'data-audio="{name}" data-file="{pin(name + ".wav")}"'

    return re.sub(r'data-audio="([^"]+)"', audio, html)


def publish_html(directory, html):
    temporary = directory / "index.pending.html"
    temporary.write_text(pin_report_assets(directory, html), encoding="utf8")
    temporary.replace(directory / "index.html")
