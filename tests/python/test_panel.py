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
MODULE_PANEL_FILES = PREVIEW_MODULE.MODULE_PANEL_FILES
control_pattern = PREVIEW_MODULE.control_pattern
control_coordinates = PREVIEW_MODULE.control_coordinates
centered_component_specs = PREVIEW_MODULE.centered_component_specs
light_pattern = PREVIEW_MODULE.light_pattern
PANEL_GRAPHICS = PREVIEW_MODULE.PANEL_GRAPHICS
render_preview = PREVIEW_MODULE.render_preview
screw_positions = PREVIEW_MODULE.screw_positions
switch_defaults = PREVIEW_MODULE.switch_defaults
svg_dimensions = PREVIEW_MODULE.svg_dimensions
aligned_source = ALIGNMENT_MODULE.aligned_source


def test_303_voice_core_runtime_panel_outlines_all_editable_text():
    source = ET.parse(ROOT / "res-src" / "Tf303VoiceCore.svg").getroot()
    runtime = ET.parse(ROOT / "res" / "Tf303VoiceCore.svg").getroot()

    assert source.findall(f".//{SVG}text")
    assert not runtime.findall(f".//{SVG}text")
    assert len(runtime.findall(f".//{SVG}path")) > len(source.findall(f".//{SVG}path"))
    assert b"\r\n" not in (ROOT / "res" / "Tf303VoiceCore.svg").read_bytes()


def test_prog_sequencer_has_three_valid_3u_widths_with_outlined_runtime_text():
    module_source = (ROOT / "src" / "TfProgSequencer.cpp").read_text(encoding="utf-8")
    parser_source = (ROOT / "src" / "tfseq_parser.cpp").read_text(encoding="utf-8")
    compiler_source = (ROOT / "src" / "tfseq_compiler.cpp").read_text(encoding="utf-8")
    runtime_source = (ROOT / "src" / "tfseq_runtime.cpp").read_text(encoding="utf-8")
    assert "ScrewSilver" not in module_source
    assert "MergeSelectionDocuments" in module_source
    assert "SplitPatternChildren" not in parser_source
    assert "ReadPatternNode" not in parser_source
    assert "SplitChordTones" not in compiler_source
    assert "NoteElement     <-" in parser_source
    assert "ExplicitVoicing" in parser_source
    assert "SlashSuffix" in parser_source
    assert "ParseArticulation(" not in compiler_source
    assert "CV1_OUTPUT" in module_source
    assert "CV2_OUTPUT" in module_source
    assert "CV3_OUTPUT" in module_source
    assert "static constexpr float MinimumHeight = 16.f" in module_source
    assert "wrappedTextHeight" in module_source
    assert "nvgTextBreakLines" in module_source
    assert "status->requiredHeight" in module_source
    assert "cvTraceLineSpan" in module_source
    assert "drawCvTrace(args, lane, now, span, sourcePositionsCurrent)" in module_source
    assert "visibleCvValues[index].store(cvOutputs[index].output" in module_source
    assert "stateTransferOrder" in runtime_source
    assert "activationCheckpointBeat" in module_source
    assert "activationNextStepBeat" in module_source
    assert "executionPulse.fetch_add(1, std::memory_order_release)" in module_source
    assert "cursorPulses[lane].fetch_add(1, std::memory_order_release)" in module_source
    assert "SchedulingLookaheadBeats" in module_source
    assert "guard++ < 64" not in module_source
    assert "delete pendingProgramPointer(pendingProgram.exchange(0" in module_source
    assert "PendingRestartBit" in module_source
    assert "publishSource(module->source, true)" in module_source
    assert "isKeyCommand(GLFW_KEY_SLASH, RACK_MOD_CTRL)" in module_source
    assert "isKeyCommand(GLFW_KEY_SPACE, RACK_MOD_CTRL)" in module_source
    assert "void onSelectText(const SelectTextEvent &event) override" in module_source
    assert "suppressShortcutSpace = true" in module_source
    assert "event.codepoint == static_cast<std::uint32_t>(' ')" in module_source
    assert (
        "isKeyCommand(GLFW_KEY_SPACE,\n                           RACK_MOD_CTRL | RACK_MOD_SHIFT)"
        in module_source
    )
    assert "isKeyCommand(GLFW_KEY_R, RACK_MOD_CTRL | RACK_MOD_SHIFT)" in module_source
    assert (
        "isKeyCommand(GLFW_KEY_BACKSPACE,\n                           RACK_MOD_CTRL | RACK_MOD_SHIFT)"
        in module_source
    )
    assert "getCompleteCablesOnPort(runPort)" in module_source
    assert "tftransport::RequestModuleCommand" in module_source
    assert "cable->outputId, command" in module_source
    assert '"QUEUED - activates on the next quarter beat"' in module_source
    assert "TransportStatus { Waiting, Playing, Paused, Stopped }" in module_source
    assert 'TransportStatus::Stopped ? "STOP"' in module_source
    assert (
        "transportStatus.load(std::memory_order_relaxed) !=\n"
        "          TransportStatus::Stopped" in module_source
    )
    assert "isKeyCommand(GLFW_KEY_D, RACK_MOD_CTRL)" in module_source
    assert "editorRunEnabled.load(std::memory_order_relaxed)" in module_source
    assert "change->oldCursor = previousCursor" in module_source
    assert "change->newCursor = cursor" in module_source
    assert "restoredEditorCursor.exchange(-1" in module_source
    assert screw_positions(module_source, 450.0, 380.0) == ()
    assert centered_component_specs("TfProgSequencer", module_source, 450.0) == (
        ("PJ301MPort", 409.0, 65.0),
        ("PJ301MPort", 437.0, 65.0),
        ("PJ301MPort", 423.0, 121.0),
        ("PJ301MPort", 409.0, 188.0),
        ("PJ301MPort", 437.0, 188.0),
        ("PJ301MPort", 409.0, 244.0),
        ("PJ301MPort", 437.0, 244.0),
        ("PJ301MPort", 409.0, 300.0),
        ("PJ301MPort", 437.0, 300.0),
        ("PJ301MPort", 409.0, 356.0),
        ("PJ301MPort", 437.0, 356.0),
    )
    for suffix, width in (("", "330"), ("-30", "450"), ("-38", "570")):
        source = ET.parse(ROOT / "res-src" / f"TfProgSequencer{suffix}.svg").getroot()
        runtime = ET.parse(ROOT / "res" / f"TfProgSequencer{suffix}.svg").getroot()
        assert source.attrib["width"] == width
        assert runtime.attrib["width"] == width
        assert source.attrib["height"] == "380"
        assert runtime.attrib["height"] == "380"
        assert source.findall(f".//{SVG}text")
        assert not runtime.findall(f".//{SVG}text")
        labels = {node.text: node for node in source.findall(f".//{SVG}text")}
        assert labels["TRIGGERFISH"].attrib["y"] == labels["PROG SEQUENCER"].attrib["y"]
        center = int(width) - 27
        left, right = center - 14, center + 14
        assert (labels["CLOCK"].attrib["x"], labels["RESET"].attrib["x"]) == (
            str(left),
            str(right),
        )
        assert labels["CLOCK"].attrib["y"] == labels["RESET"].attrib["y"] == "49"
        assert labels["RUN"].attrib["x"] == str(center)
        assert labels["RUN"].attrib["y"] == "103"
        assert labels["LOW = PAUSE"].attrib["x"] == str(center)
        assert labels["LOW = PAUSE"].attrib["y"] == "111"
        for first, second, y in (
            ("V/OCT", "GATE", "172"),
            ("TRIG", "VEL", "228"),
            ("ACC", "CV1", "284"),
            ("CV2", "CV3", "340"),
        ):
            assert (labels[first].attrib["x"], labels[second].attrib["x"]) == (
                str(left),
                str(right),
            )
            assert labels[first].attrib["y"] == labels[second].attrib["y"] == y


def test_preview_screws_follow_each_module_widget_specification():
    for module_name in MODULE_NAMES:
        widget_source = (ROOT / "src" / f"{module_name}.cpp").read_text(
            encoding="utf-8"
        )
        panel_name = MODULE_PANEL_FILES.get(module_name, f"{module_name}.svg")
        panel = ROOT / "res-src" / panel_name
        if not panel.is_file():
            panel = ROOT / "res" / panel_name
        width, height = svg_dimensions(panel)
        positions = screw_positions(widget_source, width, height)
        assert len(positions) == (0 if module_name == "TfProgSequencer" else 4)


def test_transport_panel_exposes_one_precise_fixed_rate_master_clock():
    source = (ROOT / "src" / "TfTransport.cpp").read_text(encoding="utf-8")
    panel = ET.parse(ROOT / "res-src" / "TfTransport.svg").getroot()
    labels = {node.text: node for node in panel.findall(f".//{SVG}text")}
    controls = {
        match.group("id") for match in control_pattern("TfTransport").finditer(source)
    }

    assert panel.attrib["width"] == "120"
    assert "FIXED 24 PPQN CLOCK" not in labels
    assert "120.00 BPM" not in labels
    assert {"RESTART", "PAUSE", "PLAY", "STOP"} <= labels.keys()
    assert controls == {
        "TEMPO",
        "PLAY_FROM_BEGINNING",
        "PAUSE",
        "PLAY",
        "STOP",
        "CLOCK_OUTPUT",
        "RUN_OUTPUT",
        "RESET_OUTPUT",
    }
    assert "rack::string::f(" in source
    assert '"%.2f BPM"' in source
    assert "nvgFontSize(args.vg, 9.f)" in source
    assert "NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE" in source
    assert "tempoParam->createContextMenu()" in source
    preview = (ROOT / "tools" / "render_panel_preview.py").read_text(encoding="utf-8")
    assert 'id="transport-tempo-preview"' in preview


def test_room_reverb_tooltips_use_physical_units():
    source = (ROOT / "src" / "TfReverb.cpp").read_text(encoding="utf-8")
    defaults = (ROOT / "src" / "tfdsp" / "reverb_defaults.hpp").read_text(
        encoding="utf-8"
    )

    assert "configParam<DecayQuantity>" in source
    assert '"Late decay", " s"' in source
    assert "configParam<LowCutQuantity>" in source
    assert re.search(r'"Wet low cut",\s*" Hz"', source)
    assert "configParam<HighCutQuantity>" in source
    assert re.search(r'"Wet high cut",\s*" Hz"', source)
    assert "struct HeightQuantity" not in source
    assert "RESERVED_HEIGHT" not in source
    assert "reverb_defaults::Height" not in source
    assert "configParam<AspectQuantity>" in source
    assert "reverb_defaults::HighCut" in source
    assert "inline constexpr float HighCut = MediumHall.highCut;" in defaults


def test_room_reverb_tooltip_descriptions_stay_concise():
    source = (ROOT / "src" / "TfReverb.cpp").read_text(encoding="utf-8")
    assignments = re.findall(
        r"description\s*=\s*((?:\s*\"[^\"]*\"\s*)+);", source, re.DOTALL
    )
    descriptions = [
        "".join(re.findall(r"\"([^\"]*)\"", block)) for block in assignments
    ]

    assert descriptions
    assert max(map(len, descriptions)) <= 72
    assert "Drag markers; double-click to reset." in source


def test_room_reverb_uses_a_two_dimensional_room_plan():
    source = (ROOT / "src" / "TfReverb.cpp").read_text(encoding="utf-8")
    defaults = (ROOT / "src" / "tfdsp" / "reverb_defaults.hpp").read_text(
        encoding="utf-8"
    )
    panel = ET.parse(ROOT / "res-src" / "TfReverb.svg").getroot()

    assert "struct TfRoomPlanWidget" in source
    assert "createWidget<TfRoomPlanWidget>" in source
    assert "move reverb source" in source
    assert "move reverb listener" in source
    assert "SourcePlanDefaults" in source
    assert "sourceXParamId(source)" in source
    assert "sourceYParamId(source)" in source
    assert "createWidget<TfRoomPlanWidget>(Vec(17, 27))" in source
    assert "roomPlan->box.size = Vec(206, 92)" in source
    assert panel.attrib["width"] == "240"
    assert "X_POSITION_INPUT" not in source
    assert "Y_POSITION_INPUT" not in source
    assert "Z_POSITION_INPUT" not in source
    draw = source.split("void draw(const DrawArgs &args) override", 1)[1].split(
        "void onButton", 1
    )[0]
    assert "module->inputs" not in draw
    assert "roomPlanSourceCount" in draw
    assert "sourceMarkerPosition" in draw
    assert "listenerMarkerPosition" in draw
    assert "roomPlanListenerPosition" in source
    assert "sourceAutoXParamId" in source
    assert "sourceXAutomatic(source)" in source
    assert source.index("SOURCE_8_Y") < source.index("SOURCE_1_AUTO_X")
    drag = source.split("void onDragMove", 1)[1].split("void onDoubleClick", 1)[0]
    reset = source.split("void onDoubleClick", 1)[1].split("void onDragEnd", 1)[0]
    assert "sourceAutoXParamId(activeSource)" in drag and "setValue(0.f)" in drag
    assert "sourceAutoXParamId(activeSource)" in reset and "setValue(1.f)" in reset
    assert "oldAutomaticX" in source and "history::ParamChange" in source
    assert '"sourceXPlacementVersion"' in source
    assert "Module::fromJson(root)" in source
    assert "std::abs(x - SourcePlanDefaults[source][0])" in source
    assert "including when that coordinate is exactly 0.5" in (
        ROOT / "docs" / "TfReverb-technical-report.md"
    ).read_text(encoding="utf-8")
    assert "reverb_defaults::Listener[0]" in source
    assert "reverb_defaults::Listener[1]" in source
    assert "{{0.50f, 0.682f, 0.45f}}" in defaults
    control_labels = {
        label.attrib["data-control"]
        for label in panel.findall(f".//{SVG}text")
        if "data-control" in label.attrib
    }
    assert "SOURCE_X" not in control_labels
    assert "SOURCE_Y" not in control_labels
    assert "LISTENER_X" not in control_labels
    assert "LISTENER_Y" not in control_labels
    width_label = panel.find(f".//{SVG}text[@data-control='WIDTH']")
    assert width_label is not None and width_label.text == "STEREO WIDTH"
    shimmer_label = panel.find(f".//{SVG}text[@data-control='SHIMMER']")
    assert shimmer_label is not None and shimmer_label.text == "SHIMMER"
    assert panel.find(f".//{SVG}text[@data-control='LUSH_INPUT_DIFFUSION']") is None
    assert "LUSH_INPUT_DIFFUSION" not in source
    assert panel.find(f".//{SVG}text[@data-control='HEIGHT']") is None
    assert '"Stereo width", "%"' in source
    assert '"Octave shimmer", "%"' in source

    preview_markup = PREVIEW_MODULE.module_preview_markup("TfReverb")
    assert 'id="reverb-room-plan-preview"' in preview_markup
    assert 'data-preview-source-count="4"' in preview_markup
    assert preview_markup.count('r="5.2" fill="#ffb032"') == 4
    assert 'fill="#36c8eb"' in preview_markup
    assert "SOURCES" in preview_markup
    assert "LISTENER" in preview_markup

    preview = ET.fromstring(
        f'<svg xmlns="http://www.w3.org/2000/svg">{preview_markup}</svg>'
    )
    source_markers = preview.findall(f".//{SVG}circle[@fill='#ffb032']")
    listener_marker = preview.find(f".//{SVG}circle[@fill='#36c8eb']")
    assert listener_marker is not None
    assert float(listener_marker.attrib["cy"]) > max(
        float(marker.attrib["cy"]) for marker in source_markers
    )
    source_x = [float(marker.attrib["cx"]) for marker in source_markers]
    assert source_x == sorted(source_x)
    assert source_x[0] < float(listener_marker.attrib["cx"]) < source_x[-1]

    title_lines = [
        "".join(label.itertext())
        for label in panel.findall(f".//{SVG}text")
        if "ROOM REVERB" in "".join(label.itertext())
    ]
    assert title_lines == ["TRIGGERFISH ROOM REVERB"]
    title = next(
        label
        for label in panel.findall(f".//{SVG}text")
        if "ROOM REVERB" in "".join(label.itertext())
    )
    title_parts = title.findall(f"{SVG}tspan")
    assert [part.text for part in title_parts] == ["TRIGGERFISH", " ROOM REVERB"]
    assert title_parts[0].attrib.get("font-weight") == "bold"
    assert title_parts[1].attrib.get("font-weight") is None

    controls = list(control_pattern("TfReverb").finditer(source))
    by_id = {control.group("id"): control for control in controls}
    assert [
        control_coordinates(by_id[name], source)
        for name in ("SPACE", "ASPECT", "PRE_DELAY", "DECAY", "DAMPING")
    ] == [(6.0, 130.0), (54.0, 130.0), (102.0, 130.0), (150.0, 130.0), (198.0, 130.0)]
    assert [
        control_coordinates(by_id[name], source)
        for name in ("DIFFUSION", "MODULATION", "SHIMMER", "WIDTH")
    ] == [(6.0, 184.0), (54.0, 184.0), (102.0, 184.0), (150.0, 184.0)]
    assert control_coordinates(by_id["BALANCE"], source) == (198.0, 184.0)
    assert [
        control_coordinates(by_id[name], source)
        for name in ("LOW_CUT", "HIGH_CUT", "MIX", "LEVEL")
    ] == [(33.826, 242.0), (81.826, 242.0), (129.826, 242.0), (177.826, 242.0)]
    knob_labels = [
        label
        for label in panel.findall(f".//{SVG}text")
        if label.attrib.get("data-control") in by_id
        and not label.attrib["data-control"].endswith(("INPUT", "OUTPUT"))
    ]
    assert sorted(
        sum(1 for label in knob_labels if label.attrib["y"] == row)
        for row in ("127", "181", "239")
    ) == [4, 5, 5]
    assert {
        control_coordinates(by_id[name], source)[1]
        for name in (
            "AUDIO_INPUT",
            "LISTENER_X_CV_INPUT",
            "DECAY_CV_INPUT",
            "DAMPING_CV_INPUT",
            "LEFT_OUTPUT",
        )
    } == {292.0}
    assert {
        control_coordinates(by_id[name], source)[1]
        for name in (
            "PRE_DELAY_CV_INPUT",
            "LISTENER_Y_CV_INPUT",
            "BALANCE_CV_INPUT",
            "MIX_CV_INPUT",
            "RIGHT_OUTPUT",
        )
    } == {339.0}
    assert "SPACE_CV_INPUT" not in source
    assert "ProgressiveSourceX(source, sourceCount)" in source
    assert "Margin + normalized.y *" in source
    assert "1.f - (screen.y - Margin)" not in source


def test_room_reverb_context_menu_exposes_complete_acoustic_presets():
    source = (ROOT / "src" / "TfReverb.cpp").read_text(encoding="utf-8")
    defaults = (ROOT / "src" / "tfdsp" / "reverb_defaults.hpp").read_text(
        encoding="utf-8"
    )

    for label in ("Medium Hall (default)", "Small Room", "Superlush"):
        assert f'"{label}"' in source
    assert "Centred dry" not in source
    assert 'createMenuLabel("Direct sound")' not in source
    for preset in ("MediumHall", "SmallRoom", "Superlush"):
        assert f"ReverbPreset {preset}" in defaults
    apply_body = source.split("void applyPreset", 1)[1].split("private:", 1)[0]
    assert "setParameter(LOW_CUT" in apply_body
    assert "setParameter(HIGH_CUT" in apply_body
    assert "setParameter(LEVEL" not in apply_body
    assert "source < SourcePlanDefaults.size()" in apply_body
    assert "setParameter(sourceXParamId(source), SourcePlanDefaults" in apply_body
    assert "setParameter(sourceYParamId(source), preset.source[1])" in apply_body
    assert "setParameter(sourceAutoXParamId(source), 1.f)" in apply_body


def test_303_oscillator_runtime_panel_outlines_all_editable_text():
    source = ET.parse(ROOT / "res-src" / "Tf303Oscillator.svg").getroot()
    runtime = ET.parse(ROOT / "res" / "Tf303Oscillator.svg").getroot()

    assert source.attrib["width"] == "180"
    assert source.findall(f".//{SVG}text")
    assert not runtime.findall(f".//{SVG}text")
    assert len(runtime.findall(f".//{SVG}path")) > len(source.findall(f".//{SVG}path"))
    assert b"\r\n" not in (ROOT / "res-src" / "Tf303Oscillator.svg").read_bytes()
    assert b"\r\n" not in (ROOT / "res" / "Tf303Oscillator.svg").read_bytes()


def test_wavefold_oscillator_panel_uses_compact_12hp_alive_layout():
    source = ET.parse(ROOT / "res-src" / "TfWavefoldOscillator.svg").getroot()
    runtime = ET.parse(ROOT / "res" / "TfWavefoldOscillator.svg").getroot()

    assert source.attrib["width"] == "180"
    assert runtime.attrib["width"] == "180"
    assert source.findall(f".//{SVG}text")
    assert not runtime.findall(f".//{SVG}text")

    dividers = source.find(f".//{SVG}path[@id='section-dividers']")
    assert dividers is not None
    assert dividers.attrib["d"] == "M18 27h144"

    widget_source = (ROOT / "src" / "TfWavefoldOscillator.cpp").read_text(
        encoding="utf-8"
    )
    controls = list(control_pattern("TfWavefoldOscillator").finditer(widget_source))
    by_id = {control.group("id"): control for control in controls}
    alive_params = ("ALIVE_SPEED", "MORPH_ALIVE", "FOLD_ALIVE", "SYMMETRY_ALIVE")
    assert [by_id[name].group("type") for name in alive_params] == ["TfTrimpot"] * 4
    assert [
        control_coordinates(by_id[name], widget_source) for name in alive_params
    ] == [
        (15.07, 171.0),
        (59.07, 171.0),
        (103.07, 171.0),
        (147.07, 171.0),
    ]
    assert [
        control_coordinates(by_id[name], widget_source)
        for name in ("FM_AMOUNT", "MORPH_AMOUNT", "FOLD_AMOUNT", "SYMMETRY_AMOUNT")
    ] == [
        (7.0, 257.0),
        (53.0, 257.0),
        (99.0, 257.0),
        (145.0, 257.0),
    ]
    assert all(control.group("type") != "TfSlider" for control in controls)
    assert by_id["OCTAVE"].group("type") == "TfSnapKnob"
    assert by_id["UNISON_VOICES"].group("type") == "TfSnapKnob"
    assert [
        control_coordinates(by_id[name], widget_source)
        for name in ("OCTAVE", "FM_MODE", "TUNE", "UNISON_VOICES")
    ] == [
        (9.83, 48.0),
        (59.0, 51.853),
        (99.07, 53.244),
        (135.83, 48.0),
    ]
    assert control_coordinates(by_id["CHARACTER"], widget_source) == (37.5, 205.0)
    assert control_coordinates(by_id["UNISON_SPREAD"], widget_source) == (
        121.0,
        205.0,
    )
    assert '"Alive drift time", " s"' in widget_source
    assert '"Wave CV amount", "%"' in widget_source
    assert '"Wave CV (sine / triangle morph)"' in widget_source
    assert (
        '"External audio to folder (replaces internal oscillator at folder input)"'
        in widget_source
    )
    assert '"Internal oscillator before folding"' in widget_source
    assert (
        '"Folder output (internal oscillator or external audio input)"' in widget_source
    )
    assert "TfSvgWatermark" not in widget_source


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
    assert len(controls) == 37
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
    assert input_blocks is None
    assert output_blocks is not None
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
            "FILTER_MOD_AMOUNT",
            "RES_CV_AMOUNT",
            "FILTER_MOD_MODE",
            "VCA_INITIAL_GAIN",
            "VCA_LINEAR_AMOUNT",
            "VCA_MOD_AMOUNT",
            "VCA_MOD_MODE",
        )
    ] == [
        (62.0, 108.0),
        (17.0, 161.0),
        (115.0, 161.0),
        (66.0, 161.0),
        (47.0, 249.0),
        (17.0, 205.0),
        (66.0, 205.0),
        (115.0, 205.0),
        (99.0, 249.0),
    ]
    assert "configSwitch(FILTER_ENV_MODE, 0.0f, 2.0f, 1.0f" in widget_source
    assert "configSwitch(AMP_ENV_MODE, 0.0f, 2.0f, 1.0f" in widget_source
    assert by_id["FILTER_ENV_MODE"].group("type") == "CKSSThree"
    assert by_id["AMP_ENV_MODE"].group("type") == "CKSSThree"
    assert "configSwitch(FILTER_MOD_ROUTING, 0.0f, 1.0f, 1.0f" in widget_source
    assert "configSwitch(VCA_MOD_ROUTING, 0.0f, 1.0f, 1.0f" in widget_source
    assert by_id["FILTER_MOD_ROUTING"].group("type") == "CKSS"
    assert by_id["VCA_MOD_ROUTING"].group("type") == "CKSS"
    assert by_id["FILTER_MOD_MODE"].group("type") == "CKSS"
    assert by_id["VCA_MOD_MODE"].group("type") == "CKSS"
    assert by_id["AMP_ENV_LAW"].group("type") == "CKSS"
    assert control_coordinates(by_id["FILTER_MOD_ROUTING"], widget_source) == (
        24.0,
        115.0,
    )
    assert control_coordinates(by_id["VCA_MOD_ROUTING"], widget_source) == (
        122.0,
        115.0,
    )
    assert by_id["ENVELOPE_CURVE"].group("type") == "TfCvKnob"
    assert (
        float(by_id["ENVELOPE_CURVE"].group("x")),
        float(by_id["ENVELOPE_CURVE"].group("y")),
    ) == (322.0, 143.0)
    assert control_coordinates(by_id["FILTER_ENV_MODE"], widget_source) == (
        329.0,
        78.0,
    )
    assert control_coordinates(by_id["AMP_ENV_MODE"], widget_source) == (
        329.0,
        194.0,
    )
    assert control_coordinates(by_id["AMP_ENV_LAW"], widget_source) == (
        329.0,
        241.0,
    )
    assert '{"2x (lower CPU)", "4x (default)"}' in widget_source
    assert "lightDivider.setDivision(512);" in widget_source
    assert "args.sampleTime * lightDivider.getDivision()" in widget_source
    assert widget_source.count("createLightParam<TfEnvelopeSlider>") == 8
    components = (ROOT / "src" / "components.hpp").read_text(encoding="utf-8")
    assert "TfEnvelopeSliderLight : RectangleLight<BlueLight>" in components
    assert "displayGain = 2.0f" in widget_source
    defaults = switch_defaults(widget_source)
    assert defaults["FILTER_ENV_MODE"] == 1.0
    assert defaults["AMP_ENV_MODE"] == 1.0
    assert defaults["FILTER_MOD_ROUTING"] == 1.0
    assert defaults["VCA_MOD_ROUTING"] == 1.0
    assert defaults["FILTER_MOD_MODE"] == 1.0
    assert defaults["VCA_MOD_MODE"] == 0.0
    assert defaults["AMP_ENV_LAW"] == 1.0
    assert [
        (float(by_id[name].group("x")), float(by_id[name].group("y")))
        for name in (
            "AUDIO_INPUT",
            "VOCT_INPUT",
            "FILTER_MOD_INPUT",
            "RES_CV_INPUT",
            "GATE_INPUT",
            "TRIGGER_INPUT",
        )
    ] == [(18.0 + 60.0 * index, 293.0) for index in range(6)]


def test_4072_voice_core_wrapper_preserves_normals_polyphony_and_runtime_state():
    source = (ROOT / "src" / "Tf4072VoiceCore.cpp").read_text(encoding="utf-8")

    # Every connected signal can establish the voice count, and monophonic
    # ports are explicitly broadcast by Rack's getPolyVoltage() contract.
    assert "for (int input = 0; input < NUM_INPUTS; ++input)" in source
    assert "channels = std::max(channels, inputs[input].getChannels());" in source
    assert source.count("getPolyVoltage(channel)") >= 8
    assert "for (int channel = channels; channel < activeChannels; ++channel)" in source
    assert "ResetChannel(channel);" in source
    assert "activeChannels = channels;" in source

    # An unpatched modulation input retains its internal envelope. A patched
    # input can add to it or replace it, and each destination selects its law.
    assert "const bool includeFilterEnvelope = !filterModConnected ||" in source
    assert "exponentialFilterModulation ? filterModulation : 0.0" in source
    assert "LinearFilterModulationHzPerVolt * filterModulation" in source
    assert "const bool exponentialAmpEnvelope =" in source
    assert "tfdsp::RouteArp4019ControlVoltages(" in source
    assert "addVcaModulationToEnvelope, exponentialAmpEnvelope," in source
    assert "exponentialVcaModulation);" in source

    # The audio VCA override selects its independent oversampled path.
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
    assert detents is None
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
    assert labels_by_text["−3"].attrib["x"] == "18"
    assert labels_by_text["+3"].attrib["x"] == "67"
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
    documented_modules = tuple(
        module_name
        for module_name in MODULE_NAMES
        if module_name != "TfUnisonOscillator"
    )
    for module_name in documented_modules:
        assert (ROOT / "doc" / f"{module_name}.png").is_file()
        assert f'src="doc/{module_name}.png"' in readme
    assert '<a href="#prog-sequencer"><img src="doc/TfProgSequencer.png"' in readme
    assert (
        "[Complete Prog Sequencer manual and language reference]"
        "(docs/TfProgSequencer-reference.md)" in readme
    )
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
        source = re.sub(r"```.*?```", "", source, flags=re.DOTALL)
        source = re.sub(r"`[^`\n]*`", "", source)
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


def test_legacy_pitch_and_vca_modules_expose_polyphonic_state_and_routing():
    slop = (ROOT / "src" / "TfSlop.cpp").read_text(encoding="utf-8")
    slop4 = (ROOT / "src" / "TfSlop4.cpp").read_text(encoding="utf-8")
    vdpo = (ROOT / "src" / "TfVDPO.cpp").read_text(encoding="utf-8")
    vca = (ROOT / "src" / "TfVCA.cpp").read_text(encoding="utf-8")

    # Slop mirrors its pitch cable. Hum is stepped once per frame, while drift
    # is indexed by channel so voices do not wobble as one rigid chord.
    assert "PORT_MAX_CHANNELS> _ou{};" in slop
    assert "inputs[VOCT_INPUT].getChannels()" in slop
    assert "outputs[VOCT_OUTPUT].setChannels(channels);" in slop
    assert "inputs[VOCT_INPUT].getPolyVoltage(channel)" in slop
    assert "_ou[channel].Step(_rng)" in slop

    # Slop4 retains one instrument-wide common process and adds a distinct
    # process for every Rack voice / oscillator-lane pair. Each lane keeps its
    # own input channel count instead of silently expanding sibling lanes.
    assert "PORT_MAX_CHANNELS> _ouIndividual{};" in slop4
    assert "_ouCommon.Step(_rng)" in slop4
    assert "inputs[oscillator].getChannels()" in slop4
    assert "outputs[oscillator].setChannels(channels);" in slop4
    assert "inputs[oscillator].getPolyVoltage(channel)" in slop4
    assert "_ouIndividual[channel][oscillator].Step(_rng)" in slop4

    # VDPO and VCA establish voice count from their widest input, broadcast
    # mono inputs, and own nonlinear/filter state for every possible channel.
    assert "std::array<std::unique_ptr<Oscillator>, PORT_MAX_CHANNELS>" in vdpo
    assert "inputs[DAMPING_INPUT].getChannels(), 1}" in vdpo
    assert "getPolyVoltage(channel)" in vdpo
    assert "outputs[OUTPUT].setChannels(channels);" in vdpo
    assert "_vdpHq[channel]->StepLogAngularFrequency" in vdpo

    assert "std::array<std::unique_ptr<Vca>, PORT_MAX_CHANNELS>" in vca
    assert vca.count("PORT_MAX_CHANNELS>") >= 3
    assert "inputs[EXP_CV_INPUT].getChannels(), 1}" in vca
    assert "getPolyVoltage(channel)" in vca
    assert "outputs[MAIN_OUTPUT].setChannels(channels);" in vca
    assert "_vcaTransi[channel]->StepControls" in vca
    assert "maximumControl = std::max(maximumControl, reconstructedCv);" in vca


def test_scene_pack_panel_preview_covers_every_control_and_runtime_text_is_outlined():
    widget_source = (ROOT / "src" / "TfScenePack4.cpp").read_text(encoding="utf-8")
    controls = list(control_pattern("TfScenePack4").finditer(widget_source))
    by_id = {control.group("id"): control for control in controls}

    assert len(controls) == 5
    assert len(by_id) == 5
    assert control_coordinates(by_id["SOURCE_1_INPUT"], widget_source) == (
        13.15,
        87.0,
    )
    assert control_coordinates(by_id["SOURCE_4_INPUT"], widget_source) == (
        53.15,
        149.0,
    )
    assert control_coordinates(by_id["AUDIO_OUTPUT"], widget_source) == (
        33.15,
        287.0,
    )
    assert "std::clamp(input.getChannels(), 0, MaximumSources)" in widget_source

    source_panel = ET.parse(ROOT / "res-src" / "TfScenePack4.svg").getroot()
    runtime_panel = ET.parse(ROOT / "res" / "TfScenePack4.svg").getroot()
    assert source_panel.attrib["width"] == runtime_panel.attrib["width"] == "90"
    assert source_panel.findall(f".//{SVG}text")
    assert not runtime_panel.findall(f".//{SVG}text")


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
        + sum(
            (
                1
                if control.group("type") in {"CKSS", "CKSSThree"}
                else len(COMPONENTS[control.group("type")])
            )
            for control in controls
        )
    )
    assert len(images) == expected_images
    assert all(
        image.attrib["href"].startswith("data:image/svg+xml;base64,")
        for image in images
    )
    assert str(tmp_path) not in preview.read_text(encoding="utf-8")
