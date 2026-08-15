"""Render a Rack panel with its component-library controls, without running Rack.

The editable panel SVG remains the visual source. Widget positions are parsed
from the module's C++ constructor, and Rack's installed component SVGs are
embedded as data URLs in a standalone preview. A Chromium-family browser is
used for the optional PNG render when one is available.
"""

from __future__ import annotations

import argparse
import base64
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
SVG = "{http://www.w3.org/2000/svg}"

CONTROL_PATTERN = re.compile(
    r"add(?P<kind>Param|Input|Output)\s*\(\s*"
    r"create(?:Param|Input|Output)<(?P<type>[^>]+)>\s*\(\s*"
    r"Vec\(\s*(?P<x>[0-9.]+)\s*,\s*(?P<y>[0-9.]+)\s*\)\s*,\s*module\s*,\s*"
    r"TfDiodeLadderFilter::(?P<id>[A-Z0-9_]+)",
    re.MULTILINE,
)


def control_pattern(module_name: str) -> re.Pattern[str]:
    return re.compile(
        r"add(?P<kind>Param|Input|Output)\s*\(\s*"
        r"create(?:Param|Input|Output)<(?P<type>[^>]+)>\s*\(\s*"
        r"Vec\(\s*(?P<x>[0-9.]+)\s*,\s*(?P<y>[0-9.]+)\s*\)\s*,\s*module\s*,\s*"
        + re.escape(module_name)
        + r"::(?P<id>[A-Z0-9_]+)",
        re.MULTILINE,
    )


COMPONENTS = {
    "TfLargeAudioKnob": (
        "Davies1900hLargeBlack_bg.svg",
        "Davies1900hLargeBlack.svg",
    ),
    "TfAudioKob": ("Davies1900hBlack_bg.svg", "Davies1900hBlack.svg"),
    "TfCvKnob": ("RoundBlackKnob_bg.svg", "RoundBlackKnob.svg"),
    "TfSnapKnob": ("RoundBlackKnob_bg.svg", "RoundBlackKnob.svg"),
    "TfTrimpot": ("Trimpot_bg.svg", "Trimpot.svg"),
    "CKSS": ("CKSS_0.svg",),
    "PJ301MPort": ("PJ301M.svg",),
}

PANEL_GRAPHICS = ((ROOT / "res" / "TfDiodeConnectedTransistor.svg", 1.0, 44.0),)
MODULE_GRAPHICS = {
    "TfDiodeLadderFilter": PANEL_GRAPHICS,
    "Tf303Oscillator": (),
}

# Representative positions make the preview read naturally. They affect only
# the rendered knob indicators; parameter behavior remains defined in C++.
KNOB_ANGLES = {
    "CUTOFF": 0.0,
    "RESONANCE": -145.0,
    "DRIVE": 35.0,
    "ACCENT_SWEEP_MODE": 48.0,
    "BASS": -145.0,
    "ENV_AMOUNT": -50.0,
    "NORMAL_DECAY": -15.0,
    "ACCENT_AMOUNT": 0.0,
    "ACCENT_DECAY": -35.0,
    "VCA_DECAY": 0.0,
    "CV_AMOUNT": 0.0,
    "FM_AMOUNT": -145.0,
    "RES_AMOUNT": -145.0,
    "VCA_CV_AMOUNT": 145.0,
    "OCTAVE": 0.0,
    "TUNE": 0.0,
    "SLIDE_TIME": -25.0,
    "SHAPE": 0.0,
    "WAVE": -145.0,
    "TIME_AMOUNT": 0.0,
    "SHAPE_AMOUNT": 0.0,
    "WAVE_AMOUNT": 0.0,
    "FM_MODE": 0.0,
}

MODULE_KNOB_ANGLES = {
    "Tf303Oscillator": {
        "FM_AMOUNT": 0.0,
        "SLIDE_TIME": 43.0,
    }
}


def svg_dimensions(path: Path) -> tuple[float, float]:
    root = ET.parse(path).getroot()

    def pixels(value: str | None) -> float:
        if value is None:
            raise ValueError(f"SVG has no dimensions: {path}")
        return float(re.sub(r"[^0-9.+-]", "", value))

    width = root.get("width")
    height = root.get("height")
    if width is not None and height is not None:
        return pixels(width), pixels(height)

    view_box = root.get("viewBox")
    if view_box:
        _, _, view_width, view_height = (float(value) for value in view_box.split())
        return view_width, view_height
    raise ValueError(f"SVG has no dimensions: {path}")


def embedded_image(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/svg+xml;base64,{encoded}"


def image_element(
    asset: Path, x: float, y: float, *, angle: float | None = None
) -> str:
    width, height = svg_dimensions(asset)
    transform = ""
    if angle is not None:
        transform = (
            f' transform="rotate({angle:g} {x + width / 2:g} {y + height / 2:g})"'
        )
    return (
        f'  <image x="{x:g}" y="{y:g}" width="{width:g}" height="{height:g}"'
        f'{transform} href="{embedded_image(asset)}"/>\n'
    )


def find_browser() -> Path | None:
    configured = os.environ.get("PANEL_PREVIEW_BROWSER")
    candidates = [
        configured,
        shutil.which("msedge"),
        shutil.which("chrome"),
        shutil.which("chromium"),
        r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return Path(candidate)
    return None


def render_png(
    browser: Path, preview: Path, output: Path, width: float, height: float
) -> None:
    with tempfile.TemporaryDirectory(prefix="browser-", dir=output.parent) as profile:
        command = [
            str(browser),
            "--headless",
            "--disable-gpu",
            "--hide-scrollbars",
            "--no-first-run",
            "--force-device-scale-factor=1",
            f"--window-size={int(2 * width)},{int(2 * height)}",
            f"--user-data-dir={profile}",
            f"--screenshot={output.resolve()}",
            preview.resolve().as_uri(),
        ]
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode != 0 or not output.exists():
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"browser PNG render failed: {detail}")


def render_preview(
    rack_runtime: Path,
    output_directory: Path,
    png: bool,
    module_name: str = "TfDiodeLadderFilter",
) -> Path:
    if module_name not in MODULE_GRAPHICS:
        raise ValueError(f"Unsupported panel module: {module_name}")
    panel_source = ROOT / "res-src" / f"{module_name}.svg"
    widget_source = ROOT / "src" / f"{module_name}.cpp"
    component_directory = rack_runtime / "res" / "ComponentLibrary"
    if not component_directory.is_dir():
        raise FileNotFoundError(
            f"Rack component library not found at {component_directory}. "
            "Set RACK_RUNTIME_DIR to the Rack installation directory."
        )

    panel_width, panel_height = svg_dimensions(panel_source)
    panel = panel_source.read_text(encoding="utf-8")
    panel = re.sub(
        r'width="[0-9.]+" height="[0-9.]+"',
        f'width="{2 * panel_width:g}" height="{2 * panel_height:g}"',
        panel,
        count=1,
    )
    widgets = widget_source.read_text(encoding="utf-8")
    controls = list(control_pattern(module_name).finditer(widgets))
    if not controls:
        raise RuntimeError(f"No panel controls were found in {module_name}.cpp")

    overlays = ['\n<g id="rack-component-preview">\n']
    for asset, x, y in MODULE_GRAPHICS[module_name]:
        overlays.append(image_element(asset, x, y))

    screw = component_directory / "ScrewSilver.svg"
    for x, y in (
        (15.0, 0.0),
        (panel_width - 30.0, 0.0),
        (15.0, panel_height - 15.0),
        (panel_width - 30.0, panel_height - 15.0),
    ):
        overlays.append(image_element(screw, x, y))

    for control in controls:
        component_type = control.group("type")
        asset_names = COMPONENTS.get(component_type)
        if asset_names is None:
            raise KeyError(f"No preview assets configured for {component_type}")
        x = float(control.group("x"))
        y = float(control.group("y"))
        parameter_id = control.group("id")
        for index, asset_name in enumerate(asset_names):
            module_angles = MODULE_KNOB_ANGLES.get(module_name, {})
            angle = (
                module_angles.get(parameter_id, KNOB_ANGLES.get(parameter_id))
                if index == len(asset_names) - 1
                else None
            )
            overlays.append(
                image_element(component_directory / asset_name, x, y, angle=angle)
            )
    overlays.append("</g>\n")

    output_directory.mkdir(parents=True, exist_ok=True)
    preview = output_directory / f"{module_name}-preview.svg"
    preview.write_text(
        panel.replace("</svg>", "".join(overlays) + "</svg>"), encoding="utf-8"
    )

    if png:
        browser = find_browser()
        if browser is None:
            raise FileNotFoundError(
                "No Chromium-family browser found. Set PANEL_PREVIEW_BROWSER or use --no-png."
            )
        render_png(
            browser,
            preview,
            output_directory / f"{module_name}-preview.png",
            panel_width,
            panel_height,
        )
    return preview


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rack-runtime",
        type=Path,
        default=Path(os.environ.get("RACK_RUNTIME_DIR", ROOT.parent / "Rack2")),
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=ROOT / "build" / "panel-preview",
    )
    parser.add_argument("--no-png", action="store_true")
    parser.add_argument(
        "--module",
        choices=sorted(MODULE_GRAPHICS),
        default="TfDiodeLadderFilter",
    )
    arguments = parser.parse_args()
    preview = render_preview(
        arguments.rack_runtime,
        arguments.output_directory,
        not arguments.no_png,
        arguments.module,
    )
    print(preview)


if __name__ == "__main__":
    main()
