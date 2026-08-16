import importlib.util
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).parents[2]
SVG = "{http://www.w3.org/2000/svg}"

PREVIEW_SPEC = importlib.util.spec_from_file_location(
    "render_panel_preview", ROOT / "tools" / "render_panel_preview.py"
)
assert PREVIEW_SPEC is not None and PREVIEW_SPEC.loader is not None
PREVIEW_MODULE = importlib.util.module_from_spec(PREVIEW_SPEC)
PREVIEW_SPEC.loader.exec_module(PREVIEW_MODULE)
sys.modules["render_panel_preview"] = PREVIEW_MODULE
ALIGNMENT_SPEC = importlib.util.spec_from_file_location(
    "align_panel_labels", ROOT / "tools" / "align_panel_labels.py"
)
assert ALIGNMENT_SPEC is not None and ALIGNMENT_SPEC.loader is not None
ALIGNMENT_MODULE = importlib.util.module_from_spec(ALIGNMENT_SPEC)
ALIGNMENT_SPEC.loader.exec_module(ALIGNMENT_MODULE)
COMPONENTS = PREVIEW_MODULE.COMPONENTS
CONTROL_PATTERN = PREVIEW_MODULE.CONTROL_PATTERN
MODULE_NAMES = PREVIEW_MODULE.MODULE_NAMES
control_pattern = PREVIEW_MODULE.control_pattern
control_coordinates = PREVIEW_MODULE.control_coordinates
light_pattern = PREVIEW_MODULE.light_pattern
PANEL_GRAPHICS = PREVIEW_MODULE.PANEL_GRAPHICS
render_preview = PREVIEW_MODULE.render_preview
svg_dimensions = PREVIEW_MODULE.svg_dimensions
aligned_source = ALIGNMENT_MODULE.aligned_source


def test_303_voice_core_runtime_panel_outlines_all_editable_text():
    source = ET.parse(ROOT / "res-src" / "Tf303VoiceCore.svg").getroot()
    runtime = ET.parse(ROOT / "res" / "Tf303VoiceCore.svg").getroot()

    assert source.findall(f".//{SVG}text")
    assert not runtime.findall(f".//{SVG}text")
    assert len(runtime.findall(f".//{SVG}path")) > len(source.findall(f".//{SVG}path"))
    assert b"\r\n" not in (ROOT / "res" / "Tf303VoiceCore.svg").read_bytes()


def test_303_oscillator_runtime_panel_outlines_all_editable_text():
    source = ET.parse(ROOT / "res-src" / "Tf303Oscillator.svg").getroot()
    runtime = ET.parse(ROOT / "res" / "Tf303Oscillator.svg").getroot()

    assert source.attrib["width"] == "180"
    assert source.findall(f".//{SVG}text")
    assert not runtime.findall(f".//{SVG}text")
    assert len(runtime.findall(f".//{SVG}path")) > len(source.findall(f".//{SVG}path"))
    assert b"\r\n" not in (ROOT / "res-src" / "Tf303Oscillator.svg").read_bytes()
    assert b"\r\n" not in (ROOT / "res" / "Tf303Oscillator.svg").read_bytes()


def test_4072_voice_core_panel_matches_the_wide_dual_envelope_layout():
    source = ET.parse(ROOT / "res-src" / "Tf4072VoiceCore.svg").getroot()
    runtime = ET.parse(ROOT / "res" / "Tf4072VoiceCore.svg").getroot()
    assert source.attrib["width"] == "360"
    assert source.findall(f".//{SVG}text")
    assert not runtime.findall(f".//{SVG}text")
    assert b"\r\n" not in (ROOT / "res-src" / "Tf4072VoiceCore.svg").read_bytes()
    assert b"\r\n" not in (ROOT / "res" / "Tf4072VoiceCore.svg").read_bytes()

    widget_source = (ROOT / "src" / "Tf4072VoiceCore.cpp").read_text(encoding="utf-8")
    controls = list(control_pattern("Tf4072VoiceCore").finditer(widget_source))
    assert len(controls) == 34
    by_id = {control.group("id"): control for control in controls}
    envelope_params = (
        "FILTER_ATTACK",
        "FILTER_DECAY",
        "FILTER_SUSTAIN",
        "FILTER_RELEASE",
        "AMP_ATTACK",
        "AMP_DECAY",
        "AMP_SUSTAIN",
        "AMP_RELEASE",
    )
    assert [by_id[name].group("type") for name in envelope_params] == [
        "TfEnvelopeSlider"
    ] * 8
    assert [
        (float(by_id[name].group("x")), float(by_id[name].group("y")))
        for name in envelope_params
    ] == [
        (174.0, 52.0),
        (214.0, 52.0),
        (254.0, 52.0),
        (294.0, 52.0),
        (174.0, 177.0),
        (214.0, 177.0),
        (254.0, 177.0),
        (294.0, 177.0),
    ]
    assert svg_dimensions(ROOT / "res" / "TfSlider.svg") == (20.0, 100.0)
    assert svg_dimensions(ROOT / "res" / "TfSliderHandle.svg") == (7.0, 13.0)

    dividers = source.find(f".//{SVG}path[@id='section-dividers']")
    assert dividers is not None
    assert dividers.attrib["d"] == "M28 27h304"
    assert "v" not in dividers.attrib["d"]

    input_blocks = source.find(f".//{SVG}g[@id='input-jack-blocks']")
    output_blocks = source.find(f".//{SVG}g[@id='output-jack-blocks']")
    assert input_blocks is not None
    assert output_blocks is not None
    assert len(input_blocks.findall(f"{SVG}rect")) == 10
    assert len(output_blocks.findall(f"{SVG}rect")) == 4
    assert [
        block.attrib.get("fill", output_blocks.attrib["fill"])
        for block in output_blocks.findall(f"{SVG}rect")
    ] == ["#545454", "#545454", "#858585", "#858585"]
    assert [
        (float(by_id[name].group("x")), float(by_id[name].group("y")))
        for name in (
            "DRIVE",
            "FILTER_ENV_AMOUNT",
            "CUTOFF_CV_AMOUNT",
            "RES_CV_AMOUNT",
            "VCA_INITIAL_GAIN",
            "VCA_LINEAR_AMOUNT",
            "VCA_EXP_AMOUNT",
            "FILTER_RANGE",
        )
    ] == [
        (62.0, 112.0),
        (16.0, 173.0),
        (66.0, 173.0),
        (116.0, 173.0),
        (16.0, 231.0),
        (66.0, 231.0),
        (116.0, 231.0),
        (73.0, 65.0),
    ]
    assert "configSwitch(FILTER_ENV_MODE, 0.0f, 1.0f, 0.0f" in widget_source
    assert "configSwitch(AMP_ENV_MODE, 0.0f, 1.0f, 1.0f" in widget_source
    assert '{"2x (lower CPU)", "4x (default)"}' in widget_source
    assert "lightDivider.setDivision(512);" in widget_source
    assert "args.sampleTime * lightDivider.getDivision()" in widget_source
    assert widget_source.count("createLightParam<TfEnvelopeSlider>") == 8


def test_4072_voice_core_wrapper_preserves_normals_polyphony_and_runtime_state():
    source = (ROOT / "src" / "Tf4072VoiceCore.cpp").read_text(encoding="utf-8")

    # Every connected signal can establish the voice count, and monophonic
    # ports are explicitly broadcast by Rack's getPolyVoltage() contract.
    assert "for (int input = 0; input < NUM_INPUTS; ++input)" in source
    assert "channels = std::max(channels, inputs[input].getChannels());" in source
    assert source.count("getPolyVoltage(channel)") >= 10
    assert "for (int channel = channels; channel < activeChannels; ++channel)" in source
    assert "ResetChannel(channel);" in source
    assert "activeChannels = channels;" in source

    # The external filter and VCA controls replace their corresponding
    # normalled envelopes, while the audio VCA override selects its independent
    # oversampled path.
    assert "inputs[FILTER_ENV_CV_INPUT].isConnected() ?" in source
    assert "inputs[VCA_LINEAR_CV_INPUT].isConnected() ?" in source
    assert "const bool vcaOverride = inputs[VCA_AUDIO_INPUT].isConnected();" in source
    assert "if (vcaOverride)" in source

    # Runtime quality changes reset the newly selected DSP path rather than
    # resuming stale filter and resampler history.
    assert "if (oversampling != activeOversampling)" in source
    assert "activeOversampling = oversampling;" in source
    assert "ResetActivePath();" in source

    # Bypass follows the VCA input normalization instead of always routing the
    # filter input, and non-audio envelope outputs do not retain stale channels.
    assert "void processBypass(const ProcessArgs& args) override" in source
    assert "configBypass(" not in source
    assert (
        "route(inputs[VCA_AUDIO_INPUT].isConnected() ? VCA_AUDIO_INPUT :\n"
        "\t\t\tAUDIO_INPUT, VCA_OUTPUT);"
    ) in source
    assert "outputs[FILTER_ENV_OUTPUT].setChannels(0);" in source
    assert "outputs[AMP_ENV_OUTPUT].setChannels(0);" in source


def test_303_oscillator_panel_matches_widget_layout():
    widget_source = (ROOT / "src" / "Tf303Oscillator.cpp").read_text(encoding="utf-8")
    assert widget_source.count("oversampling = 1;") == 2
    assert widget_source.count("activeOversampling = 1;") == 2
    assert '{"2x (lower CPU)", "4x (default)"}' in widget_source
    controls = list(control_pattern("Tf303Oscillator").finditer(widget_source))
    assert len(controls) == 19
    by_id = {control.group("id"): control for control in controls}
    assert [
        by_id[name].group("type")
        for name in ("OCTAVE", "TUNE", "SLIDE_TIME", "SHAPE", "WAVE")
    ] == [
        "TfRotarySwitchKnob",
        "TfLargeAudioKnob",
        "TfAudioKob",
        "TfAudioKob",
        "TfAudioKob",
    ]
    assert [
        (float(by_id[name].group("x")), float(by_id[name].group("y")))
        for name in ("OCTAVE", "TUNE", "SLIDE_TIME", "SHAPE", "WAVE")
    ] == [
        (20.0, 49.0),
        (110.0, 45.0),
        (15.0, 120.0),
        (72.0, 120.0),
        (129.0, 120.0),
    ]
    amount_ids = ("FM_AMOUNT", "TIME_AMOUNT", "SHAPE_AMOUNT", "WAVE_AMOUNT")
    assert [by_id[name].group("type") for name in amount_ids] == ["TfCvKnob"] * 4
    assert [
        (float(by_id[name].group("x")), float(by_id[name].group("y")))
        for name in amount_ids
    ] == [(10.0, 178.0), (54.0, 178.0), (98.0, 178.0), (142.0, 178.0)]
    assert (
        float(by_id["FM_MODE"].group("x")),
        float(by_id["FM_MODE"].group("y")),
    ) == (83.0, 54.0)
    assert [
        (float(by_id[name].group("x")), float(by_id[name].group("y")))
        for name in ("VOCT_INPUT", "SLIDE_INPUT", "TIME_INPUT", "SYNC_INPUT")
    ] == [(12.0, 241.0), (56.0, 241.0), (100.0, 241.0), (144.0, 241.0)]
    assert [
        (float(by_id[name].group("x")), float(by_id[name].group("y")))
        for name in ("FM_INPUT", "SHAPE_INPUT", "WAVE_INPUT")
    ] == [(27.0, 284.0), (78.0, 284.0), (129.0, 284.0)]
    assert [
        (float(by_id[name].group("x")), float(by_id[name].group("y")))
        for name in ("CV_OUTPUT", "AUDIO_OUTPUT")
    ] == [(48.0, 344.0), (108.0, 344.0)]
    assert (
        float(by_id["SYNC_INPUT"].group("x")),
        float(by_id["SYNC_INPUT"].group("y")),
    ) == (144.0, 241.0)

    source = ET.parse(ROOT / "res-src" / "Tf303Oscillator.svg").getroot()
    detents = source.find(f".//{SVG}g[@id='octave-detents']")
    assert detents is not None
    assert len(detents.findall(f"{SVG}circle")) == 7
    dividers = source.find(f".//{SVG}path[@id='section-dividers']")
    assert dividers is not None
    assert dividers.attrib["d"] == "M28 27h124M28 222h124M28 324h124"
    labels = ["".join(label.itertext()) for label in source.findall(f".//{SVG}text")]
    assert labels.count("SL.TIME") == 1
    assert labels.count("SL.TIME CV") == 2
    labels_by_text = {
        "".join(label.itertext()): label for label in source.findall(f".//{SVG}text")
    }
    assert labels_by_text["TZ"].attrib["y"] == "50"
    assert labels_by_text["EXP"].attrib["y"] == "85"
    aligned_labels = [
        label
        for label in source.findall(f".//{SVG}text")
        if "data-control" in label.attrib
    ]
    centered_labels = [
        label
        for label in source.findall(f".//{SVG}text")
        if "data-control" in label.attrib or "data-center-x" in label.attrib
    ]
    assert len(aligned_labels) == 21
    assert len(centered_labels) == 23
    output_blocks = source.find(f".//{SVG}g[@id='output-jack-blocks']")
    assert output_blocks is not None
    assert output_blocks.attrib["fill"] == "#545454"
    blocks = output_blocks.findall(f"{SVG}rect")
    assert [
        (float(block.attrib["x"]), float(block.attrib["y"])) for block in blocks
    ] == [
        (39.0, 338.0),
        (99.0, 338.0),
    ]
    assert blocks[0].attrib["fill"] == "#858585"
    assert blocks[1].attrib.get("fill", output_blocks.attrib["fill"]) == "#545454"

    logo = PREVIEW_MODULE.MODULE_GRAPHICS["Tf303Oscillator"]
    assert logo == ((ROOT / "res" / "logo.svg", 16.0, 232.0, 148.0, 80.8, 0.12),)
    assert 'pluginInstance, "res/logo.svg"' in widget_source
    assert "logoGraphic->box.pos = Vec(16, 232);" in widget_source
    assert "logoGraphic->opacity = 0.12f;" in widget_source


def test_panel_label_alignment_uses_control_center_and_optical_offset():
    class FixedOpticalCentering:
        @staticmethod
        def offset(text, font_size, letter_spacing):
            assert text == "LABEL"
            assert font_size == 7.0
            assert letter_spacing == 0.0
            return 0.375

    source = '<text x="1" y="2" font-size="7" ' 'data-control="CONTROL">LABEL</text>'
    aligned, count = aligned_source(source, {"CONTROL": 12.5}, FixedOpticalCentering())
    assert count == 1
    assert 'x="12.875"' in aligned


def test_documentation_includes_every_rendered_module_preview_and_technical_report():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    for module_name in MODULE_NAMES:
        assert (ROOT / "doc" / f"{module_name}.png").is_file()
        assert f'src="doc/{module_name}.png"' in readme
    for report in (
        "Tf303Oscillator-technical-report.md",
        "Tf303VoiceCore-technical-report.md",
    ):
        assert (ROOT / "docs" / report).is_file()
        assert f"docs/{report}" in readme


def test_readme_and_technical_report_local_links_resolve():
    documents = [ROOT / "README.md", *sorted((ROOT / "docs").glob("*.md"))]
    link_pattern = re.compile(r"\[[^]]*\]\(([^)]+)\)|(?:src|href)=\"([^\"]+)\"")
    for document in documents:
        source = document.read_text(encoding="utf-8")
        for match in link_pattern.finditer(source):
            target = next(value for value in match.groups() if value)
            if target.startswith(("#", "http://", "https://")):
                continue
            local_path = target.split("#", 1)[0]
            assert (
                document.parent / local_path
            ).is_file(), f"Broken local link in {document.relative_to(ROOT)}: {target}"


def test_legacy_panel_coordinate_expressions_are_resolved():
    for module_name, expected in (
        ("TfSlop4", {"TRACK_SCALING4": (118.0, 223.0), "VOCT_OUTPUT4": (115.0, 319.0)}),
        ("TfVCA", {"LIN_CV_INPUT": (15.0, 313.0), "MAIN_OUTPUT": (141.0, 313.0)}),
    ):
        source = (ROOT / "src" / f"{module_name}.cpp").read_text(encoding="utf-8")
        controls = list(control_pattern(module_name).finditer(source))
        positions = {
            control.group("id"): control_coordinates(control, source)
            for control in controls
        }
        for control_id, position in expected.items():
            assert positions[control_id] == position


def test_vca_preview_includes_activity_light(tmp_path):
    component_directory = tmp_path / "Rack2" / "res" / "ComponentLibrary"
    component_directory.mkdir(parents=True)
    asset_names = {"ScrewSilver.svg", "MediumLight.svg"}
    for names in COMPONENTS.values():
        asset_names.update(names)
    minimal_svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10" '
        'viewBox="0 0 10 10"><circle cx="5" cy="5" r="4"/></svg>'
    )
    for asset_name in asset_names:
        (component_directory / asset_name).write_text(minimal_svg, encoding="utf-8")

    widget_source = (ROOT / "src" / "TfVCA.cpp").read_text(encoding="utf-8")
    lights = list(light_pattern("TfVCA").finditer(widget_source))
    assert len(lights) == 1
    assert lights[0].group("type") == "MediumLight<BlueLight>"
    assert control_coordinates(lights[0], widget_source) == (85.0, 250.0)

    preview = render_preview(
        tmp_path / "Rack2", tmp_path / "output", png=False, module_name="TfVCA"
    )
    preview_source = preview.read_text(encoding="utf-8")
    assert 'fill="#333"' in preview_source
    assert "MediumLight.svg" not in preview_source


def test_303_voice_core_uses_separate_diode_connected_transistor_artwork():
    source = ET.parse(ROOT / "res-src" / "Tf303VoiceCore.svg").getroot()
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
    assert PANEL_GRAPHICS[0] == (asset, 64.0, 266.0, 112.0, 112.0, 0.32)

    widget_source = (ROOT / "src" / "Tf303VoiceCore.cpp").read_text(encoding="utf-8")
    assert 'pluginInstance, "res/TfDiodeConnectedTransistor.svg"' in widget_source
    assert "transistorGraphic->box.pos = Vec(64, 266);" in widget_source
    assert "transistorGraphic->opacity = 0.32f;" in widget_source


def test_303_voice_core_panel_uses_triggerfish_output_block_convention():
    source = ET.parse(ROOT / "res-src" / "Tf303VoiceCore.svg").getroot()
    assert source.attrib["width"] == "240"

    dividers = source.find(f".//{SVG}path[@id='section-dividers']")
    assert dividers is not None
    assert dividers.attrib["d"] == "M28 27h184M28 264h184"

    labels = {
        "".join(label.itertext()): label for label in source.findall(f".//{SVG}text")
    }
    assert "303 VOICE CORE" in labels
    assert labels["DRIVE"].attrib["x"] == "47.5"
    assert labels["BASS"].attrib["x"] == "189"
    assert labels["RES RANGE"].attrib["data-control"] == "HIGH_RESONANCE"
    assert labels["STOCK"].attrib["data-control"] == "HIGH_RESONANCE"
    assert labels["HIGH"].attrib["data-control"] == "HIGH_RESONANCE"
    assert labels["RES RANGE"].attrib["y"] == "48"
    assert labels["HIGH"].attrib["y"] == "61"
    assert labels["STOCK"].attrib["y"] == "96"
    secondary_labels = source.find(f".//{SVG}g[@fill='#686868']")
    assert secondary_labels is not None
    assert ["".join(label.itertext()) for label in secondary_labels] == [
        "HIGH",
        "STOCK",
    ]

    output_blocks = source.find(f".//{SVG}g[@id='output-jack-blocks']")
    assert output_blocks is not None
    assert output_blocks.attrib["fill"] == "#545454"
    assert [
        (float(block.attrib["x"]), float(block.attrib["y"]))
        for block in output_blocks.findall(f"{SVG}rect")
    ] == [(147.0, 327.0), (195.0, 327.0)]


def test_panel_preview_embeds_rack_components_at_cpp_widget_positions(tmp_path):
    widget_source = (ROOT / "src" / "Tf303VoiceCore.cpp").read_text(encoding="utf-8")
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
