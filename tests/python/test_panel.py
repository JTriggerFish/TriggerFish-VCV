import importlib.util
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).parents[2]
SVG = "{http://www.w3.org/2000/svg}"

PREVIEW_SPEC = importlib.util.spec_from_file_location(
    "render_panel_preview", ROOT / "tools" / "render_panel_preview.py"
)
assert PREVIEW_SPEC is not None and PREVIEW_SPEC.loader is not None
PREVIEW_MODULE = importlib.util.module_from_spec(PREVIEW_SPEC)
PREVIEW_SPEC.loader.exec_module(PREVIEW_MODULE)
COMPONENTS = PREVIEW_MODULE.COMPONENTS
CONTROL_PATTERN = PREVIEW_MODULE.CONTROL_PATTERN
PANEL_GRAPHICS = PREVIEW_MODULE.PANEL_GRAPHICS
render_preview = PREVIEW_MODULE.render_preview
svg_dimensions = PREVIEW_MODULE.svg_dimensions


def test_diode_ladder_runtime_panel_outlines_all_editable_text():
    source = ET.parse(ROOT / "res-src" / "TfDiodeLadderFilter.svg").getroot()
    runtime = ET.parse(ROOT / "res" / "TfDiodeLadderFilter.svg").getroot()

    assert source.findall(f".//{SVG}text")
    assert not runtime.findall(f".//{SVG}text")
    assert len(runtime.findall(f".//{SVG}path")) > len(source.findall(f".//{SVG}path"))
    assert b"\r\n" not in (ROOT / "res" / "TfDiodeLadderFilter.svg").read_bytes()


def test_diode_ladder_documentation_includes_rendered_module_preview():
    preview = ROOT / "doc" / "TfDiodeLadderFilter.png"
    assert preview.is_file()
    assert 'src="doc/TfDiodeLadderFilter.png"' in (ROOT / "README.md").read_text(
        encoding="utf-8"
    )


def test_diode_ladder_uses_separate_diode_connected_transistor_artwork():
    source = ET.parse(ROOT / "res-src" / "TfDiodeLadderFilter.svg").getroot()
    assert source.find(f".//{SVG}g[@id='diode-connected-transistor']") is None

    asset = ROOT / "res" / "TfDiodeConnectedTransistor.svg"
    graphic = ET.parse(asset).getroot()
    assert graphic.attrib["viewBox"] == "0 0 77 77"
    assert graphic.attrib["width"] == graphic.attrib["height"] == "31"
    assert svg_dimensions(asset) == (31.0, 31.0)
    assert len(graphic.findall(f".//{SVG}path")) == 1
    paths = graphic.findall(f".//{SVG}path")
    assert [path.attrib.get("opacity") for path in paths] == ["0.48"]
    assert {element.attrib.get("stroke") for element in graphic.iter()} - {None} == {
        "#526870"
    }
    assert {element.attrib.get("fill") for element in graphic.iter()} - {None} == {
        "none"
    }
    assert len(PANEL_GRAPHICS) == 1
    assert PANEL_GRAPHICS[0] == (asset, 1.0, 44.0)

    widget_source = (ROOT / "src" / "TfDiodeLadderFilter.cpp").read_text(
        encoding="utf-8"
    )
    assert 'pluginInstance, "res/TfDiodeConnectedTransistor.svg"' in widget_source
    assert "transistorGraphic->box.pos = Vec(1, 44);" in widget_source


def test_diode_ladder_panel_uses_triggerfish_output_block_convention():
    source = ET.parse(ROOT / "res-src" / "TfDiodeLadderFilter.svg").getroot()
    assert source.attrib["width"] == "240"

    labels = {
        "".join(label.itertext()): label for label in source.findall(f".//{SVG}text")
    }
    assert labels["DRIVE"].attrib["x"] == "47.5"
    assert labels["BASS"].attrib["x"] == "189"

    output_blocks = source.find(f".//{SVG}g[@id='output-jack-blocks']")
    assert output_blocks is not None
    assert output_blocks.attrib["fill"] == "#545454"
    assert [
        (float(block.attrib["x"]), float(block.attrib["y"]))
        for block in output_blocks.findall(f"{SVG}rect")
    ] == [(147.0, 327.0), (195.0, 327.0)]


def test_panel_preview_embeds_rack_components_at_cpp_widget_positions(tmp_path):
    widget_source = (ROOT / "src" / "TfDiodeLadderFilter.cpp").read_text(
        encoding="utf-8"
    )
    controls = list(CONTROL_PATTERN.finditer(widget_source))
    assert len(controls) == 25

    by_id = {control.group("id"): control for control in controls}
    articulation_ids = (
        "ENV_AMOUNT",
        "NORMAL_DECAY",
        "ACCENT_AMOUNT",
        "ACCENT_DECAY",
        "VCA_DECAY",
    )
    modulation_ids = ("CV_AMOUNT", "FM_AMOUNT", "RES_AMOUNT", "VCA_CV_AMOUNT")
    assert [by_id[name].group("type") for name in articulation_ids] == ["TfCvKnob"] * 5
    assert [float(by_id[name].group("x")) for name in articulation_ids] == [
        10.0,
        58.0,
        106.0,
        154.0,
        202.0,
    ]
    assert [float(by_id[name].group("x")) for name in modulation_ids] == [
        34.0,
        82.0,
        130.0,
        178.0,
    ]
    assert [by_id[name].group("type") for name in modulation_ids] == ["TfCvKnob"] * 4

    ports = [control for control in controls if control.group("kind") != "Param"]
    assert [(float(port.group("x")), float(port.group("y"))) for port in ports] == [
        (12.0, 282.0),
        (60.0, 282.0),
        (108.0, 282.0),
        (156.0, 282.0),
        (204.0, 282.0),
        (12.0, 334.0),
        (60.0, 334.0),
        (108.0, 334.0),
        (156.0, 334.0),
        (204.0, 334.0),
    ]

    component_directory = tmp_path / "Rack2" / "res" / "ComponentLibrary"
    component_directory.mkdir(parents=True)
    asset_names = {"ScrewSilver.svg"}
    for names in COMPONENTS.values():
        asset_names.update(names)
    minimal_svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10" '
        'viewBox="0 0 10 10"><circle cx="5" cy="5" r="4"/></svg>'
    )
    for asset_name in asset_names:
        (component_directory / asset_name).write_text(minimal_svg, encoding="utf-8")

    preview = render_preview(tmp_path / "Rack2", tmp_path / "output", png=False)
    preview_root = ET.parse(preview).getroot()
    images = preview_root.findall(f".//{SVG}image")
    expected_images = (
        len(PANEL_GRAPHICS)
        + 4
        + sum(len(COMPONENTS[control.group("type")]) for control in controls)
    )
    assert len(images) == expected_images
    assert all(
        image.attrib["href"].startswith("data:image/svg+xml;base64,")
        for image in images
    )
    assert str(tmp_path) not in preview.read_text(encoding="utf-8")
