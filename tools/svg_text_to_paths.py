"""Convert SVG text into paths for VCV Rack's NanoSVG renderer.

The editable source SVG may use text, but Rack does not render SVG text nodes.
This small build helper keeps the panel source readable while producing a
self-contained runtime asset. Font files are supplied explicitly and are not
copied into the repository.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

from fontTools.misc.transform import Transform
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

SVG_NAMESPACE = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NAMESPACE)


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("source", type=Path)
    result.add_argument("destination", type=Path)
    result.add_argument("--font", type=Path, required=True)
    result.add_argument("--bold-font", type=Path)
    return result


class Outliner:
    inherited_attributes = (
        "fill",
        "font-family",
        "font-size",
        "font-weight",
        "text-anchor",
    )

    def __init__(self, regular_path: Path, bold_path: Path | None):
        self.regular = TTFont(regular_path)
        self.bold = TTFont(bold_path) if bold_path else self.regular
        self.has_bold_font = bold_path is not None

    @staticmethod
    def number(value, default=0.0):
        if value is None:
            return default
        return float(value.removesuffix("px"))

    def convert_children(self, parent, inherited=None):
        inherited = dict(inherited or {})
        for key in self.inherited_attributes:
            if key in parent.attrib:
                inherited[key] = parent.attrib[key]

        for index, child in list(enumerate(list(parent))):
            child_inherited = dict(inherited)
            for key in self.inherited_attributes:
                if key in child.attrib:
                    child_inherited[key] = child.attrib[key]
            if child.tag == f"{{{SVG_NAMESPACE}}}text":
                parent.remove(child)
                parent.insert(index, self.outline(child, child_inherited))
            else:
                self.convert_children(child, child_inherited)

    def outline(self, text_node, attributes):
        content = "".join(text_node.itertext())
        weight = attributes.get("font-weight", "400")
        is_bold = weight == "bold" or (weight.isdigit() and int(weight) >= 600)
        font = self.bold if is_bold else self.regular
        glyph_set = font.getGlyphSet()
        cmap = font.getBestCmap()
        metrics = font["hmtx"].metrics
        units_per_em = font["head"].unitsPerEm
        font_size = self.number(attributes.get("font-size"), 12.0)
        scale = font_size / units_per_em
        letter_spacing = self.number(text_node.attrib.get("letter-spacing"), 0.0)

        glyphs = []
        advance = 0.0
        for character in content:
            glyph_name = cmap.get(ord(character), ".notdef")
            glyph_advance = metrics[glyph_name][0] * scale
            glyphs.append((glyph_name, advance))
            advance += glyph_advance + letter_spacing
        if glyphs:
            advance -= letter_spacing

        x = self.number(text_node.attrib.get("x"))
        y = self.number(text_node.attrib.get("y"))
        anchor = attributes.get("text-anchor", "start")
        if anchor == "middle":
            x -= advance / 2.0
        elif anchor == "end":
            x -= advance

        commands = []
        for glyph_name, offset in glyphs:
            pen = SVGPathPen(glyph_set)
            transform = Transform(scale, 0.0, 0.0, -scale, x + offset, y)
            glyph_set[glyph_name].draw(TransformPen(pen, transform))
            command = pen.getCommands()
            if command:
                commands.append(command)

        path = ET.Element(f"{{{SVG_NAMESPACE}}}path")
        path.set("d", " ".join(commands))
        path.set("fill", attributes.get("fill", "#000000"))
        if is_bold and not self.has_bold_font:
            path.set("stroke", attributes.get("fill", "#000000"))
            path.set("stroke-width", str(font_size * 0.025))
        if "opacity" in text_node.attrib:
            path.set("opacity", text_node.attrib["opacity"])
        return path


def main():
    arguments = parser().parse_args()
    tree = ET.parse(arguments.source)
    outliner = Outliner(arguments.font, arguments.bold_font)
    outliner.convert_children(tree.getroot())
    arguments.destination.parent.mkdir(parents=True, exist_ok=True)
    tree.write(arguments.destination, encoding="utf-8", xml_declaration=True)


if __name__ == "__main__":
    main()
