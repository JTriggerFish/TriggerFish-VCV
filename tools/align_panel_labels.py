"""Center SVG labels on their Rack controls using the real component widths.

Labels opt in with ``data-control="PARAM_OR_PORT_ID"``. Control positions are
read from the module widget's C++ source and widths from Rack's component SVGs,
so fractional asset dimensions cannot accumulate into visible alignment drift.
"""

from __future__ import annotations

import argparse
import html
from pathlib import Path
import re

from fontTools.pens.boundsPen import BoundsPen
from fontTools.ttLib import TTFont

from render_panel_preview import (
    COMPONENTS,
    component_asset,
    control_pattern,
    svg_dimensions,
)

ROOT = Path(__file__).resolve().parents[1]
TEXT_PATTERN = re.compile(r"<text\b(?P<attributes>[^>]*)>(?P<content>[^<]*)</text>")
X_PATTERN = re.compile(r'\bx="[0-9.+-]+"')
ATTRIBUTE_PATTERN = re.compile(r'(?P<name>[\w.-]+)="(?P<value>[^"]*)"')


def number(value: float) -> str:
    return f"{value:.4f}".rstrip("0").rstrip(".")


def control_centers(
    widget_source: Path, component_directory: Path, module_name: str
) -> dict[str, float]:
    source = widget_source.read_text(encoding="utf-8")
    centers = {}
    for control in control_pattern(module_name).finditer(source):
        component_type = control.group("type")
        asset_names = COMPONENTS.get(component_type)
        if asset_names is None:
            raise KeyError(f"No preview assets configured for {component_type}")
        widths = [
            svg_dimensions(component_asset(component_directory, asset_name))[0]
            for asset_name in asset_names
        ]
        centers[control.group("id")] = float(control.group("x")) + max(widths) / 2.0
    return centers


def attributes(source: str) -> dict[str, str]:
    return {
        match.group("name"): match.group("value")
        for match in ATTRIBUTE_PATTERN.finditer(source)
    }


class OpticalCentering:
    def __init__(self, font_path: Path):
        self.font = TTFont(font_path)
        self.glyph_set = self.font.getGlyphSet()
        self.cmap = self.font.getBestCmap()
        self.metrics = self.font["hmtx"].metrics
        self.units_per_em = self.font["head"].unitsPerEm

    def offset(self, text: str, font_size: float, letter_spacing: float) -> float:
        scale = font_size / self.units_per_em
        advance = 0.0
        ink_min = float("inf")
        ink_max = float("-inf")
        glyph_count = 0
        for character in text:
            glyph_name = self.cmap.get(ord(character), ".notdef")
            pen = BoundsPen(self.glyph_set)
            self.glyph_set[glyph_name].draw(pen)
            if pen.bounds is not None:
                ink_min = min(ink_min, advance + pen.bounds[0] * scale)
                ink_max = max(ink_max, advance + pen.bounds[2] * scale)
            advance += self.metrics[glyph_name][0] * scale + letter_spacing
            glyph_count += 1
        if glyph_count:
            advance -= letter_spacing
        if ink_min == float("inf"):
            return 0.0
        return advance / 2.0 - (ink_min + ink_max) / 2.0

    def close(self) -> None:
        self.font.close()


def aligned_source(
    source: str, centers: dict[str, float], centering: OpticalCentering
) -> tuple[str, int]:
    aligned = 0

    def replace(match: re.Match[str]) -> str:
        nonlocal aligned
        label_attributes = attributes(match.group("attributes"))
        control_id = label_attributes.get("data-control")
        fixed_center = label_attributes.get("data-center-x")
        if control_id is None and fixed_center is None:
            return match.group(0)
        if control_id is not None and control_id not in centers:
            raise KeyError(f"Panel label refers to unknown control {control_id}")
        opening_tag = f'<text{match.group("attributes")}>'
        if X_PATTERN.search(opening_tag) is None:
            raise ValueError("Optically centered panel label has no x coordinate")
        font_size = float(label_attributes.get("font-size", "12"))
        letter_spacing = float(label_attributes.get("letter-spacing", "0"))
        content = html.unescape(match.group("content"))
        optical_offset = centering.offset(content, font_size, letter_spacing)
        target_center = (
            centers[control_id] if control_id is not None else float(fixed_center)
        )
        aligned += 1
        opening_tag = X_PATTERN.sub(
            f'x="{number(target_center + optical_offset)}"',
            opening_tag,
            count=1,
        )
        return f'{opening_tag}{match.group("content")}</text>'

    return TEXT_PATTERN.sub(replace, source), aligned


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", default="Tf303Oscillator")
    parser.add_argument(
        "--rack-runtime",
        type=Path,
        default=ROOT.parent / "Rack2",
    )
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--font", type=Path)
    arguments = parser.parse_args()

    panel_source = ROOT / "res-src" / f"{arguments.module}.svg"
    widget_source = ROOT / "src" / f"{arguments.module}.cpp"
    component_directory = arguments.rack_runtime / "res" / "ComponentLibrary"
    font_path = (
        arguments.font or arguments.rack_runtime / "res" / "fonts" / "DejaVuSans.ttf"
    )
    centers = control_centers(widget_source, component_directory, arguments.module)
    original = panel_source.read_text(encoding="utf-8")
    centering = OpticalCentering(font_path)
    try:
        aligned, count = aligned_source(original, centers, centering)
    finally:
        centering.close()
    if count == 0:
        raise ValueError(f"No data-control labels found in {panel_source}")
    if arguments.check:
        if aligned != original:
            raise SystemExit(
                f"{panel_source} has labels that are not centered; run "
                "tools/align_panel_labels.py"
            )
    else:
        panel_source.write_text(aligned, encoding="utf-8", newline="\n")
        if aligned != original:
            print(f"Centered {count} labels in {panel_source}")


if __name__ == "__main__":
    main()
