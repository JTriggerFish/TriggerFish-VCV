import json
import math
from pathlib import Path
import re

import pytest

ROOT = Path(__file__).parents[2]
PLUGIN_VERSION = json.loads((ROOT / "plugin.json").read_text(encoding="utf-8"))[
    "version"
]
PATCH_PATHS = {
    "slop4": ROOT / "test-slop4.vcv",
    "vdpo": ROOT / "test-vdpo.vcv",
    "303_voice": ROOT / "test-303-voice.vcv",
}

# Canonical Rack 2.6 module tags. Rack accepts a few historical aliases, but
# manifests should publish the canonical names used by its module browser.
RACK_MODULE_TAGS = {
    "Arpeggiator",
    "Attenuator",
    "Blank",
    "Chorus",
    "Clock generator",
    "Clock modulator",
    "Compressor",
    "Controller",
    "Delay",
    "Digital",
    "Distortion",
    "Drum",
    "Dual",
    "Dynamics",
    "Effect",
    "Envelope follower",
    "Envelope generator",
    "Equalizer",
    "Expander",
    "External",
    "Filter",
    "Flanger",
    "Function generator",
    "Granular",
    "Hardware clone",
    "Limiter",
    "Logic",
    "Low-frequency oscillator",
    "Low-pass gate",
    "MIDI",
    "Mixer",
    "Multiple",
    "Noise",
    "Oscillator",
    "Panning",
    "Phaser",
    "Physical modeling",
    "Polyphonic",
    "Quad",
    "Quantizer",
    "Random",
    "Recording",
    "Reverb",
    "Ring modulator",
    "Sample and hold",
    "Sampler",
    "Sequencer",
    "Slew limiter",
    "Speech",
    "Switch",
    "Synth voice",
    "Tuner",
    "Utility",
    "Visual",
    "Vocoder",
    "Voltage-controlled amplifier",
    "Waveshaper",
}

EXPECTED_DEFAULTS = {
    "TfSlop": {0: 0.25, 1: 0.05, 2: 1.0, 3: -1.0},
    "TfSlop4": {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 0.1, 5: 0.05, 6: 0.05},
    "TfVDPO": {0: 0.5, 1: 0.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0},
    "TfVCA": {0: 0.5, 1: 1.0, 2: 1.0, 3: 0.5, 4: 50.0, 5: 1.0},
    "Tf303VoiceCore": {
        0: 0.9344246,
        1: 0.0,
        2: 0.0,
        3: 0.0,
        4: 0.0,
        5: 0.5321928,
        6: 0.0,
        7: 0.0,
        8: 1.0 / 3.0,
        9: math.log10(0.5),
        10: math.log10(0.2),
        11: 0.5,
        12: 0.5,
        13: 1.0,
        14: 2.0,
    },
    "Tf303Oscillator": {
        0: 0.0,
        1: 0.0,
        2: math.log10(0.060),
        3: 0.0,
        4: 0.0,
        5: 0.0,
        6: 0.0,
        7: 0.0,
        8: 0.0,
        9: 0.0,
    },
}


def test_manifest_uses_canonical_rack_module_tags():
    manifest = json.loads((ROOT / "plugin.json").read_text(encoding="utf-8"))
    published = {tag for module in manifest["modules"] for tag in module["tags"]}

    assert published <= RACK_MODULE_TAGS


def test_release_version_is_consistent_across_build_metadata():
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")

    make_version = re.search(r"^VERSION\s*=\s*([^\s]+)$", makefile, re.MULTILINE)
    cmake_version = re.search(
        r"^project\(TriggerFishDSP VERSION ([^\s]+) LANGUAGES CXX\)$",
        cmake,
        re.MULTILINE,
    )
    assert make_version is not None
    assert cmake_version is not None
    assert make_version.group(1) == PLUGIN_VERSION
    assert cmake_version.group(1) == PLUGIN_VERSION


def test_manifest_links_current_release_notes_and_module_guides():
    manifest = json.loads((ROOT / "plugin.json").read_text(encoding="utf-8"))

    assert manifest["changelogUrl"].endswith(f"/docs/releases/{PLUGIN_VERSION}.md")
    assert (ROOT / "docs" / "releases" / f"{PLUGIN_VERSION}.md").is_file()
    assert all(module.get("manualUrl") for module in manifest["modules"])


def test_303_modules_have_final_public_names_and_slugs():
    manifest = json.loads((ROOT / "plugin.json").read_text(encoding="utf-8"))
    modules_by_slug = {module["slug"]: module for module in manifest["modules"]}

    assert modules_by_slug["Tf303Oscillator"]["name"] == "303 Oscillator"
    assert modules_by_slug["Tf303VoiceCore"]["name"] == "303 Voice Core"
    assert "TfDiodeLadderFilter" not in modules_by_slug


def load_patch(name):
    return json.loads(PATCH_PATHS[name].read_text(encoding="utf-8"))


def modules(patch, model):
    return [module for module in patch["modules"] if module["model"] == model]


def param_values(module):
    return {param["id"]: param["value"] for param in module["params"]}


def has_cable(patch, source, output, target, input_):
    return any(
        cable["outputModuleId"] == source
        and cable["outputId"] == output
        and cable["inputModuleId"] == target
        and cable["inputId"] == input_
        for cable in patch["cables"]
    )


def test_smoke_patches_collectively_contain_every_triggerfish_module():
    models = {
        module["model"]
        for name in PATCH_PATHS
        for module in load_patch(name)["modules"]
        if module["plugin"] == "TriggerFish-Elements"
    }
    assert models == {
        "TfSlop",
        "TfSlop4",
        "TfVCA",
        "TfVDPO",
        "Tf303VoiceCore",
        "Tf303Oscillator",
    }


@pytest.mark.parametrize("name", PATCH_PATHS)
def test_smoke_patches_use_triggerfish_parameter_defaults(name):
    triggerfish_modules = [
        module
        for module in load_patch(name)["modules"]
        if module["plugin"] == "TriggerFish-Elements"
    ]
    for triggerfish_module in triggerfish_modules:
        assert (
            param_values(triggerfish_module)
            == EXPECTED_DEFAULTS[triggerfish_module["model"]]
        )


@pytest.mark.parametrize("name", PATCH_PATHS)
def test_smoke_patch_dependencies_and_cables_are_valid(name):
    patch = load_patch(name)
    module_ids = {module["id"] for module in patch["modules"]}

    assert len(module_ids) == len(patch["modules"])
    assert {module["plugin"] for module in patch["modules"]} <= {
        "Core",
        "Fundamental",
        "ImpromptuModular",
        "TriggerFish-Elements",
    }
    assert len({cable["id"] for cable in patch["cables"]}) == len(patch["cables"])
    assert all(cable["outputModuleId"] in module_ids for cable in patch["cables"])
    assert all(cable["inputModuleId"] in module_ids for cable in patch["cables"])


@pytest.mark.parametrize("name", PATCH_PATHS)
def test_smoke_patch_has_playable_control_and_stereo_audio_paths(name):
    patch = load_patch(name)
    audio = modules(patch, "AudioInterface")[0]
    adsrs = modules(patch, "ADSR")

    assert "Scope" not in {module["model"] for module in patch["modules"]}
    assert "Plateau" not in {module["model"] for module in patch["modules"]}
    if modules(patch, "MIDIToCVInterface"):
        midi = modules(patch, "MIDIToCVInterface")[0]
        assert all(has_cable(patch, midi["id"], 1, adsr["id"], 4) for adsr in adsrs)
    else:
        assert modules(patch, "Foundry")
        assert modules(patch, "Clocked")
    assert {
        cable["inputId"]
        for cable in patch["cables"]
        if cable["inputModuleId"] == audio["id"]
    } == {0, 1}


@pytest.mark.parametrize("name", PATCH_PATHS)
def test_smoke_patch_does_not_name_local_hardware(name):
    patch = load_patch(name)
    audio = modules(patch, "AudioInterface")[0]

    for midi_module in modules(patch, "MIDIToCVInterface"):
        assert midi_module["data"]["midi"]["driver"] == -1
        assert "deviceName" not in midi_module["data"]["midi"]
    assert audio["data"]["audio"]["driver"] == -1
    assert "deviceName" not in audio["data"]["audio"]


@pytest.mark.parametrize("name", PATCH_PATHS)
def test_smoke_patch_uses_current_plugin_version(name):
    triggerfish_modules = [
        module
        for module in load_patch(name)["modules"]
        if module["plugin"] == "TriggerFish-Elements"
    ]
    assert {module["version"] for module in triggerfish_modules} == {PLUGIN_VERSION}


def test_slop4_patch_exercises_all_outputs_and_tunes_octaves():
    patch = load_patch("slop4")
    slop4 = modules(patch, "TfSlop4")[0]
    vcos = modules(patch, "VCO")

    assert {
        cable["outputId"]
        for cable in patch["cables"]
        if cable["outputModuleId"] == slop4["id"]
    } == {0, 1, 2, 3}
    assert sorted(param_values(vco)[2] for vco in vcos) == [0.0, 0.0, 12.0, 12.0]


def test_slop4_patch_has_one_enveloped_vca_chain():
    patch = load_patch("slop4")
    adsr = modules(patch, "ADSR")[0]
    vca = modules(patch, "TfVCA")[0]
    voice_mixer = next(
        mixer for mixer in modules(patch, "VCMixer") if param_values(mixer)[0] == 0.7
    )

    assert has_cable(patch, voice_mixer["id"], 0, vca["id"], 0)
    assert has_cable(patch, adsr["id"], 0, vca["id"], 1)


def test_vdpo_patch_has_self_resonating_and_forced_enveloped_chains():
    patch = load_patch("vdpo")
    vdpos = modules(patch, "TfVDPO")
    vcas = modules(patch, "TfVCA")
    adsrs = modules(patch, "ADSR")
    forced_inputs = {
        cable["inputModuleId"]
        for cable in patch["cables"]
        if cable["inputId"] == 1
        and cable["inputModuleId"] in {vdpo["id"] for vdpo in vdpos}
    }

    assert len(vdpos) == len(vcas) == len(adsrs) == 2
    assert len(forced_inputs) == 1
    for vdpo in vdpos:
        audio_vca = next(
            vca for vca in vcas if has_cable(patch, vdpo["id"], 0, vca["id"], 0)
        )
        assert any(
            has_cable(patch, adsr["id"], 0, audio_vca["id"], 1) for adsr in adsrs
        )


def test_vdpo_patch_selects_sine_saw_or_square_forcing():
    patch = load_patch("vdpo")
    source_vco = modules(patch, "VCO")[0]
    selector = modules(patch, "SequentialSwitch2")[0]
    push = modules(patch, "Push")[0]
    forced_vdpo_id = next(
        cable["inputModuleId"]
        for cable in patch["cables"]
        if cable["outputModuleId"] == selector["id"] and cable["inputId"] == 1
    )

    assert param_values(selector) == {0: 1.0}
    assert selector["data"]["declick"] is True
    assert has_cable(patch, push["id"], 0, selector["id"], 0)
    assert {
        (cable["outputId"], cable["inputId"])
        for cable in patch["cables"]
        if cable["outputModuleId"] == source_vco["id"]
        and cable["inputModuleId"] == selector["id"]
    } == {(0, 2), (2, 3), (3, 4)}
    assert has_cable(patch, selector["id"], 0, forced_vdpo_id, 1)


def test_303_voice_patch_connects_sequencer_oscillator_filter_and_vca():
    patch = load_patch("303_voice")
    clock = modules(patch, "Clocked")[0]
    foundry = modules(patch, "Foundry")[0]
    oscillator = modules(patch, "Tf303Oscillator")[0]
    diode = modules(patch, "Tf303VoiceCore")[0]
    master = modules(patch, "VCMixer")[0]

    assert has_cable(patch, clock["id"], 1, foundry["id"], 6)
    assert has_cable(patch, foundry["id"], 0, oscillator["id"], 0)
    assert has_cable(patch, foundry["id"], 9, oscillator["id"], 1)
    assert has_cable(patch, oscillator["id"], 0, diode["id"], 1)
    assert has_cable(patch, oscillator["id"], 1, diode["id"], 0)
    assert has_cable(patch, foundry["id"], 8, diode["id"], 5)
    assert has_cable(patch, foundry["id"], 4, diode["id"], 6)
    assert has_cable(patch, diode["id"], 1, master["id"], 1)


def test_diode_ladder_default_cutoff_cv_range_reaches_fully_open():
    diode = modules(load_patch("303_voice"), "Tf303VoiceCore")[0]
    cutoff_pitch = diode["params"][0]["value"]
    cv_amount = diode["params"][5]["value"]
    cutoff_hz = 261.625565 * 2.0**cutoff_pitch
    envelope_peak_cutoff_hz = cutoff_hz * 2.0 ** (10.0 * cv_amount)

    assert cutoff_hz == pytest.approx(500.0, rel=1.0e-6)
    assert envelope_peak_cutoff_hz == pytest.approx(20_000.0, rel=1.0e-6)


def test_303_voice_foundry_pattern_has_accents_rests_legato_and_slide_gates():
    foundry = modules(load_patch("303_voice"), "Foundry")[0]
    data = foundry["data"]
    note_attributes = data["id0_attributes"][:16]
    slide_attributes = data["id1_attributes"][:16]

    assert data["id0_sequences"][0] == data["id1_sequences"][0] == 16
    assert sum(bool(value & (1 << 24)) for value in note_attributes) == 9
    assert {
        index for index, value in enumerate(note_attributes) if (value >> 28) == 5
    } == {4, 5, 11, 12}
    assert sum((value & 0xFF) == 200 for value in note_attributes) == 3
    assert {
        index for index, value in enumerate(slide_attributes) if value & (1 << 24)
    } == {5, 12}


@pytest.mark.parametrize("name", PATCH_PATHS)
def test_smoke_patch_has_quiet_master(name):
    patch = load_patch(name)
    audio = modules(patch, "AudioInterface")[0]
    master = next(
        mixer
        for mixer in modules(patch, "VCMixer")
        if math.isclose(20 * math.log10(param_values(mixer)[0]), -6.0, abs_tol=1e-6)
    )

    assert {
        cable["inputId"]
        for cable in patch["cables"]
        if cable["outputModuleId"] == master["id"]
        and cable["inputModuleId"] == audio["id"]
    } == {0, 1}
