from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).parents[2]
SVG = "{http://www.w3.org/2000/svg}"


def test_diode_ladder_runtime_panel_outlines_all_editable_text():
    source = ET.parse(ROOT / "res-src" / "TfDiodeLadderFilter.svg").getroot()
    runtime = ET.parse(ROOT / "res" / "TfDiodeLadderFilter.svg").getroot()

    assert source.findall(f".//{SVG}text")
    assert not runtime.findall(f".//{SVG}text")
    assert len(runtime.findall(f".//{SVG}path")) > len(source.findall(f".//{SVG}path"))
