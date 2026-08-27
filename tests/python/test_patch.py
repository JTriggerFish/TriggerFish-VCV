import json
import math
from pathlib import Path
import re
import tomllib

import pytest

ROOT = Path(__file__).parents[2]
PLUGIN_VERSION = json.loads((ROOT / "plugin.json").read_text(encoding="utf-8"))[
    "version"
]
PATCH_PATHS = {
    "slop4": ROOT / "test-slop4.vcv",
    "vdpo": ROOT / "test-vdpo.vcv",
    "303_voice": ROOT / "test-303-voice.vcv",
    "prog_303": ROOT / "test-prog-sequencer-303.vcv",
    "reverb_two_sources": ROOT / "test-room-reverb-two-sources.vcv",
    "4072_voice": ROOT / "test-4072-voice.vcv",
    "electric_piano": ROOT / "test-electric-piano.vcv",
    "wavefold": ROOT / "test-wavefold-oscillator.vcv",
    "unison": ROOT / "test-unison-oscillator.vcv",
    "scene_pack": ROOT / "test-scene-pack4.vcv",
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
    "TfProgSequencer": {},
    "TfTransport": {0: 120.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
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
    "Tf4072VoiceCore": {
        0: 0.9344246,
        1: 0.0,
        2: 0.0,
        3: 0.6,
        4: 0.0,
        5: 0.0,
        6: 0.0,
        7: 1.0,
        8: 1.0,
        9: math.log10(0.0014),
        10: math.log10(1.0),
        11: 0.5,
        12: math.log10(1.0),
        13: 1.0,
        14: 0.0,
        15: math.log10(0.0014),
        16: math.log10(1.0),
        17: 1.0,
        18: math.log10(1.0),
        19: 1.0,
        20: 1.0,
        21: 1.0,
        22: 1.0,
        23: 0.0,
        24: 1.0,
    },
    "TfElectricPiano": {
        0: 0.5,
        1: 1.0,
        2: 0.62,
        3: 0.52,
        4: 0.50,
        5: 0.52,
        6: 0.55,
        7: 0.48,
        8: 0.50,
        9: 0.24,
        10: 0.18,
        11: 0.32,
        12: 0.50,
        13: 0.50,
        14: 0.50,
        15: 0.32,
        16: 0.0,
    },
    "TfWavefoldOscillator": {
        0: 0.0,
        1: 0.0,
        2: 0.5,
        3: 0.4,
        4: 0.0,
        5: 0.0,
        6: 0.0,
        7: 0.0,
        8: 0.0,
        9: 0.0,
        10: 2.0,
        11: 0.5,
        12: 0.5,
        13: 0.5,
        14: 0.5,
        15: 1.0,
        16: 0.39841330778621553,
    },
    "TfUnisonOscillator": {
        0: 0.0,
        1: 0.0,
        2: 3.0,
        3: 0.0,
        4: 0.5,
        5: 0.0,
        6: 0.42,
        7: 0.39841330778621553,
        8: 0.65,
        9: 0.0,
        10: 0.0,
        11: 0.5,
        12: 0.1,
        13: 0.1,
        14: 0.15,
        15: 1.0,
        16: 1.0,
        17: 0.0,
        18: 0.0,
    },
    "TfScenePack4": {},
    "TfReverb": {
        0: 0.60,
        1: 0.5,
        2: 0.048,
        3: 0.5,
        4: 0.35,
        5: 0.5,
        6: 0.682,
        7: 0.64,
        8: 0.28,
        9: 0.82,
        10: 0.30,
        11: 0.6130368568946039,
        12: 0.5,
        13: 0.0,
        14: 0.829482217661603,
        15: 0.35,
        16: 0.0,
        17: 0.0,
        18: 0.50,
        19: 0.35,
        20: 0.50,
        21: 0.35,
        22: 0.50,
        23: 0.35,
        24: 0.50,
        25: 0.35,
        26: 0.50,
        27: 0.35,
        28: 0.50,
        29: 0.35,
        30: 0.50,
        31: 0.35,
    },
}


def test_manifest_uses_canonical_rack_module_tags():
    manifest = json.loads((ROOT / "plugin.json").read_text(encoding="utf-8"))
    published = {tag for module in manifest["modules"] for tag in module["tags"]}

    assert published <= RACK_MODULE_TAGS


def test_release_version_is_consistent_across_build_metadata():
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    uv_lock = tomllib.loads((ROOT / "uv.lock").read_text(encoding="utf-8"))

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
    assert pyproject["project"]["version"] == PLUGIN_VERSION
    triggerfish_lock = next(
        package
        for package in uv_lock["package"]
        if package["name"] == "triggerfish-vcv-dsp"
    )
    assert triggerfish_lock["version"] == PLUGIN_VERSION


def test_current_version_respects_the_declared_copyleft_license_boundary():
    manifest = json.loads((ROOT / "plugin.json").read_text(encoding="utf-8"))
    license_text = (ROOT / "LICENSE.txt").read_text(encoding="utf-8")
    notice = (ROOT / "NOTICE.md").read_text(encoding="utf-8")

    assert tuple(map(int, PLUGIN_VERSION.split("."))) >= (2, 4, 0)
    assert manifest["license"] == "GPL-3.0-or-later"
    assert "GNU GENERAL PUBLIC LICENSE" in license_text
    assert "Version 3, 29 June 2007" in license_text
    assert "version 2.4.0 and later" in notice
    assert "Releases through version 2.3.1" in notice
    assert (ROOT / "LICENSES" / "BSD-3-Clause.txt").is_file()


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
        "Tf4072VoiceCore",
        "TfElectricPiano",
        "TfWavefoldOscillator",
        "TfUnisonOscillator",
        "TfScenePack4",
        "TfReverb",
        "TfTransport",
        "TfProgSequencer",
    }


@pytest.mark.parametrize("name", PATCH_PATHS)
def test_smoke_patches_use_triggerfish_parameter_defaults(name):
    triggerfish_modules = [
        module
        for module in load_patch(name)["modules"]
        if module["plugin"] == "TriggerFish-Elements"
    ]
    for triggerfish_module in triggerfish_modules:
        actual = param_values(triggerfish_module)
        expected = EXPECTED_DEFAULTS[triggerfish_module["model"]]
        musical_overrides = set()
        if name == "prog_303":
            musical_overrides = {
                "Tf303Oscillator": {0, 4},
                "Tf303VoiceCore": {0, 1, 2, 3, 4, 9},
                "TfReverb": {16},
            }.get(triggerfish_module["model"], set())
        elif name == "scene_pack":
            if triggerfish_module["model"] == "TfReverb":
                musical_overrides = {16}
        elif name == "reverb_two_sources":
            musical_overrides = {
                "TfUnisonOscillator": {2, 7, 8, 9, 11, 12, 13, 14},
                "TfWavefoldOscillator": {2, 3, 4, 7, 10, 11, 12, 13, 14},
                "Tf4072VoiceCore": {
                    0,
                    1,
                    2,
                    3,
                    4,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                },
                "TfReverb": set(range(32)),
                "TfTransport": {0},
            }.get(triggerfish_module["model"], set())
        assert {
            param_id: value
            for param_id, value in actual.items()
            if param_id not in musical_overrides
        } == {
            param_id: value
            for param_id, value in expected.items()
            if param_id not in musical_overrides
        }


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

    if name in {"wavefold", "unison", "scene_pack", "reverb_two_sources"}:
        assert len(modules(patch, "Scope")) == 1
    else:
        assert "Scope" not in {module["model"] for module in patch["modules"]}
    assert "Plateau" not in {module["model"] for module in patch["modules"]}
    if name == "scene_pack":
        assert len(modules(patch, "VCO")) == 4
    elif name == "reverb_two_sources":
        assert len(modules(patch, "TfProgSequencer")) == 2
        assert len(modules(patch, "Tf4072VoiceCore")) == 2
        assert modules(patch, "TfTransport")
    elif modules(patch, "MIDIToCVInterface"):
        midi = modules(patch, "MIDIToCVInterface")[0]
        assert all(has_cable(patch, midi["id"], 1, adsr["id"], 4) for adsr in adsrs)
    elif name == "prog_303":
        assert modules(patch, "TfProgSequencer")
        assert modules(patch, "TfTransport")
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


def test_electric_piano_patch_is_16_voice_velocity_playable():
    patch = load_patch("electric_piano")
    midi_module = modules(patch, "MIDIToCVInterface")[0]
    piano = modules(patch, "TfElectricPiano")[0]

    assert midi_module["data"]["channels"] == 16
    assert has_cable(patch, midi_module["id"], 0, piano["id"], 0)
    assert has_cable(patch, midi_module["id"], 1, piano["id"], 1)
    assert has_cable(patch, midi_module["id"], 2, piano["id"], 2)
    assert {
        cable["outputId"]
        for cable in patch["cables"]
        if cable["outputModuleId"] == piano["id"]
    } == {1, 2}


def test_electric_piano_exposes_coupling_and_preserves_tails():
    module_source = (ROOT / "src" / "TfElectricPiano.cpp").read_text(encoding="utf-8")

    assert "COUPLING" in module_source
    assert "PreserveVoiceAsTail" in module_source

    preserve_tail = module_source.split("void PreserveVoiceAsTail", 1)[1].split(
        "tfdsp::ElectricPianoControls Controls", 1
    )[0]
    assert "std::swap" in preserve_tail
    assert "CreateVoice()" not in preserve_tail


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
    cutoff_lfo, resonance_lfo, fm_lfo = modules(patch, "LFO")
    oscillator = modules(patch, "Tf303Oscillator")[0]
    diode = modules(patch, "Tf303VoiceCore")[0]
    master = modules(patch, "VCMixer")[0]

    assert has_cable(patch, clock["id"], 1, foundry["id"], 6)
    assert has_cable(patch, clock["id"], 4, foundry["id"], 5)
    assert has_cable(patch, foundry["id"], 0, oscillator["id"], 0)
    assert has_cable(patch, foundry["id"], 9, oscillator["id"], 1)
    assert param_values(oscillator) == EXPECTED_DEFAULTS["Tf303Oscillator"]
    assert not has_cable(patch, oscillator["id"], 0, diode["id"], 1)
    assert has_cable(patch, oscillator["id"], 1, diode["id"], 0)
    assert has_cable(patch, cutoff_lfo["id"], 0, diode["id"], 2)
    assert has_cable(patch, resonance_lfo["id"], 0, diode["id"], 4)
    assert has_cable(patch, fm_lfo["id"], 0, diode["id"], 3)
    assert param_values(diode) == EXPECTED_DEFAULTS["Tf303VoiceCore"]
    assert param_values(cutoff_lfo)[0] == 0.0
    assert param_values(cutoff_lfo)[2] == pytest.approx(-4.00999546051025)
    assert param_values(resonance_lfo)[2] == pytest.approx(-4.1243371963501)
    assert param_values(fm_lfo)[2] == pytest.approx(3.80241107940674)
    assert has_cable(patch, foundry["id"], 8, diode["id"], 5)
    assert has_cable(patch, foundry["id"], 4, diode["id"], 6)
    assert has_cable(patch, diode["id"], 1, master["id"], 1)
    assert param_values(master)[1] == pytest.approx(0.838592648506165)


def test_303_voice_patch_preserves_saved_layout_and_view():
    patch = load_patch("303_voice")
    assert patch["zoom"] == pytest.approx(1.3599998950958252)
    assert patch["gridOffset"] == pytest.approx(
        [-12.460025787353516, -0.22162829339504242]
    )
    assert {module["id"]: module["pos"] for module in patch["modules"]} == {
        1: [-9, 0],
        2: [7, 0],
        3: [27, 0],
        4: [5, 1],
        5: [30, 1],
        6: [65, 0],
        7: [74, 0],
        8: [46, 1],
        9: [55, 1],
        10: [64, 1],
    }


def test_prog_303_replaces_foundry_with_the_concise_program():
    patch = load_patch("prog_303")
    sequencer = modules(patch, "TfProgSequencer")[0]
    clock = modules(patch, "TfTransport")[0]
    oscillator = modules(patch, "Tf303Oscillator")[0]
    voice = modules(patch, "Tf303VoiceCore")[0]
    reverb = modules(patch, "TfReverb")[0]
    audio = modules(patch, "AudioInterface")[0]

    assert modules(patch, "Foundry") == []
    assert has_cable(patch, clock["id"], 0, sequencer["id"], 0)
    assert has_cable(patch, clock["id"], 2, sequencer["id"], 1)
    assert has_cable(patch, clock["id"], 1, sequencer["id"], 2)
    assert has_cable(patch, sequencer["id"], 0, oscillator["id"], 0)
    assert has_cable(patch, sequencer["id"], 1, voice["id"], 5)
    assert has_cable(patch, sequencer["id"], 4, voice["id"], 6)
    assert not has_cable(patch, sequencer["id"], 0, oscillator["id"], 1)
    assert has_cable(patch, voice["id"], 1, reverb["id"], 0)
    assert has_cable(patch, reverb["id"], 0, audio["id"], 0)
    assert has_cable(patch, reverb["id"], 1, audio["id"], 1)
    assert param_values(reverb)[16] == -6.0
    assert modules(patch, "VCMixer") == []

    source = sequencer["data"]["source"]
    assert "subdiv 16n" in source
    assert "notes ^1 1!3 ^5 7 ^1 ~ 1 ~ ^8 >1, ~ 1 >^1, >1" in source
    assert "glide .5" in source
    assert sequencer["data"]["languageVersion"] == 1
    oscillator_values = param_values(oscillator)
    voice_values = param_values(voice)
    assert oscillator_values[0] == 1.0
    assert oscillator_values[4] == pytest.approx(0.17951805889606476)
    assert voice_values[0] == pytest.approx(0.6371597647666931)
    assert voice_values[1] == pytest.approx(0.48494037985801697)
    assert voice_values[2] == 1.0
    assert voice_values[3] == pytest.approx(8.949535369873047)
    assert voice_values[4] == pytest.approx(0.14879514276981354)
    assert voice_values[9] == pytest.approx(-0.06006646901369095)
    assert "articulation" not in source
    assert "\n  accent " not in source
    assert "gate .5" in source
    assert "velocity .5" in source
    assert "iv = acid |> modulate_degree 4 |> octave -1" in source
    assert "v  = acid |> modulate_degree 5 |> octave -1 |> scale major" in source
    assert "song = acid * 8 + iv * 4 + v * 4" in source
    assert len(source.splitlines()) == 14


def test_4072_voice_patch_connects_midi_oscillator_filter_envelopes_and_vca():
    patch = load_patch("4072_voice")
    midi = modules(patch, "MIDIToCVInterface")[0]
    oscillator = modules(patch, "VCO")[0]
    voice = modules(patch, "Tf4072VoiceCore")[0]
    master = modules(patch, "VCMixer")[0]

    assert has_cable(patch, midi["id"], 0, oscillator["id"], 0)
    assert has_cable(patch, midi["id"], 0, voice["id"], 2)
    assert has_cable(patch, midi["id"], 1, voice["id"], 6)
    assert has_cable(patch, oscillator["id"], 2, voice["id"], 0)
    assert has_cable(patch, voice["id"], 1, master["id"], 1)


def test_wavefold_patch_connects_midi_oscillator_envelope_and_vca():
    patch = load_patch("wavefold")
    midi = modules(patch, "MIDIToCVInterface")[0]
    adsr = modules(patch, "ADSR")[0]
    oscillator = modules(patch, "TfWavefoldOscillator")[0]
    vca = modules(patch, "TfVCA")[0]
    scope = modules(patch, "Scope")[0]

    assert has_cable(patch, midi["id"], 0, oscillator["id"], 0)
    assert has_cable(patch, midi["id"], 1, adsr["id"], 4)
    assert has_cable(patch, oscillator["id"], 1, vca["id"], 0)
    assert has_cable(patch, adsr["id"], 0, vca["id"], 1)
    assert has_cable(patch, oscillator["id"], 1, scope["id"], 0)
    assert has_cable(patch, oscillator["id"], 0, scope["id"], 1)
    assert has_cable(patch, oscillator["id"], 0, scope["id"], 2)
    assert param_values(scope)[4] == pytest.approx(-math.log2(0.010))
    assert param_values(scope)[7] == 0.0
    assert not any(
        cable["inputModuleId"] == oscillator["id"] and cable["inputId"] == 5
        for cable in patch["cables"]
    )


def test_unison_patch_connects_polyphonic_stereo_voice_and_scope():
    patch = load_patch("unison")
    midi = modules(patch, "MIDIToCVInterface")[0]
    adsr = modules(patch, "ADSR")[0]
    oscillator = modules(patch, "TfUnisonOscillator")[0]
    vcas = modules(patch, "TfVCA")
    scope = modules(patch, "Scope")[0]
    audio = modules(patch, "AudioInterface")[0]

    assert has_cable(patch, midi["id"], 0, oscillator["id"], 0)
    assert has_cable(patch, midi["id"], 1, adsr["id"], 4)
    for output, vca in zip((1, 2), vcas, strict=True):
        assert has_cable(patch, oscillator["id"], output, vca["id"], 0)
        assert has_cable(patch, adsr["id"], 0, vca["id"], 1)
    assert {
        cable["inputId"]
        for cable in patch["cables"]
        if cable["inputModuleId"] == audio["id"]
    } == {0, 1}
    assert has_cable(patch, oscillator["id"], 1, scope["id"], 0)
    assert has_cable(patch, oscillator["id"], 2, scope["id"], 1)
    assert has_cable(patch, oscillator["id"], 1, scope["id"], 2)


def test_scene_pack_patch_connects_four_sources_and_packed_outputs():
    patch = load_patch("scene_pack")
    oscillators = modules(patch, "VCO")
    scene_pack = modules(patch, "TfScenePack4")[0]
    reverb = modules(patch, "TfReverb")[0]
    scope = modules(patch, "Scope")[0]
    audio = modules(patch, "AudioInterface")[0]

    assert len(oscillators) == 4
    for lane, oscillator in enumerate(oscillators):
        assert has_cable(patch, oscillator["id"], 0, scene_pack["id"], lane)
    assert has_cable(patch, scene_pack["id"], 0, reverb["id"], 0)
    assert has_cable(patch, reverb["id"], 0, audio["id"], 0)
    assert has_cable(patch, reverb["id"], 1, audio["id"], 1)
    assert has_cable(patch, scene_pack["id"], 0, scope["id"], 0)
    assert has_cable(patch, reverb["id"], 0, scope["id"], 1)
    assert has_cable(patch, reverb["id"], 1, scope["id"], 2)
    assert modules(patch, "VCMixer") == []


def test_two_source_reverb_patch_builds_two_independently_sequenced_voices():
    patch = load_patch("reverb_two_sources")
    clock = modules(patch, "TfTransport")[0]
    sequencers = modules(patch, "TfProgSequencer")
    bass_sequencer = next(
        sequencer
        for sequencer in sequencers
        if "bass = sequence" in sequencer["data"]["source"]
    )
    arpeggio_sequencer = next(
        sequencer
        for sequencer in sequencers
        if "arpeggio = sequence" in sequencer["data"]["source"]
    )
    unison = modules(patch, "TfUnisonOscillator")[0]
    wavefolder = modules(patch, "TfWavefoldOscillator")[0]
    voices = modules(patch, "Tf4072VoiceCore")
    bass_voice = next(
        voice for voice in voices if has_cable(patch, unison["id"], 0, voice["id"], 0)
    )
    arpeggio_voice = next(
        voice
        for voice in voices
        if has_cable(patch, wavefolder["id"], 1, voice["id"], 0)
    )
    scene_pack = modules(patch, "TfScenePack4")[0]
    reverb = modules(patch, "TfReverb")[0]
    scope = modules(patch, "Scope")[0]
    audio = modules(patch, "AudioInterface")[0]

    assert len(sequencers) == len(voices) == 2
    assert param_values(clock)[0] == 96.0
    assert param_values(clock) == {
        0: 96.0,
        1: 0.0,
        2: 0.0,
        3: 0.0,
        4: 0.0,
    }
    assert clock["data"]["state"] == 2
    for sequencer in sequencers:
        assert has_cable(patch, clock["id"], 0, sequencer["id"], 0)
        assert has_cable(patch, clock["id"], 2, sequencer["id"], 1)
        assert has_cable(patch, clock["id"], 1, sequencer["id"], 2)
        assert sequencer["data"]["languageVersion"] == 1

    bass_source = bass_sequencer["data"]["source"]
    assert "tonic D@1" in bass_source
    assert "subdiv 2n" in bass_source
    assert "scale dorian" in bass_source
    assert "notes 1_2 7_2 4_2 5 7" in bass_source
    assert "~" not in bass_source
    assert "gate .96" in bass_source
    assert param_values(unison)[2] == 2.0
    assert param_values(unison)[3] == 0.0
    assert param_values(unison)[9] == pytest.approx(0.2965061068534851)
    assert param_values(unison)[13] == pytest.approx(0.29333335161209106)
    assert param_values(unison)[14] == pytest.approx(0.3920000493526459)
    assert param_values(bass_voice)[0] == pytest.approx(3.0938026905059814)
    assert param_values(bass_voice)[1] == pytest.approx(0.4920482337474823)
    assert param_values(bass_voice)[2] == pytest.approx(10.387954711914062)
    assert param_values(bass_voice)[10] == pytest.approx(math.log10(1.8))
    assert param_values(bass_voice)[17] == pytest.approx(0.85)
    assert has_cable(patch, bass_sequencer["id"], 0, unison["id"], 0)
    assert has_cable(patch, bass_sequencer["id"], 0, bass_voice["id"], 2)
    assert has_cable(patch, bass_sequencer["id"], 1, bass_voice["id"], 6)
    assert has_cable(patch, bass_sequencer["id"], 2, bass_voice["id"], 7)
    assert has_cable(patch, bass_voice["id"], 1, scene_pack["id"], 0)

    arpeggio_source = arpeggio_sequencer["data"]["source"]
    assert "tonic D@4" in arpeggio_source
    assert "notes 3' 1' 5 3 ; 1' 5 3 1 ; 2' 7 5 2 ; 7 5 2 1" in arpeggio_source
    assert "cv1 4 -4 |> interp linear |> rate 1/9" in arpeggio_source
    assert "cv2 env ad 5ms 16n depth 2 curve 0" in arpeggio_source
    assert arpeggio_sequencer["data"]["editorHeatmap"] == 5
    assert param_values(wavefolder)[3] == pytest.approx(0.5012047290802002)
    assert param_values(wavefolder)[7] == pytest.approx(0.75)
    assert param_values(wavefolder)[10] == 1.0
    assert param_values(wavefolder)[12] == pytest.approx(0.40800002217292786)
    assert param_values(wavefolder)[13] == pytest.approx(0.49066662788391113)
    assert param_values(wavefolder)[14] == pytest.approx(0.5093333721160889)
    assert param_values(arpeggio_voice)[0] == pytest.approx(0.6673828959465027)
    assert param_values(arpeggio_voice)[1] == 0.0
    assert param_values(arpeggio_voice)[2] == pytest.approx(3.5975842475891113)
    assert param_values(arpeggio_voice)[3] == pytest.approx(0.5493977069854736)
    assert param_values(arpeggio_voice)[4] == pytest.approx(0.8789158463478088)
    assert param_values(arpeggio_voice)[8] == pytest.approx(0.7283129692077637)
    assert param_values(arpeggio_voice)[10] == pytest.approx(0.04110236465930939)
    assert param_values(arpeggio_voice)[16] == pytest.approx(0.043807897716760635)
    assert param_values(arpeggio_voice)[7] == pytest.approx(1.0)
    assert has_cable(patch, arpeggio_sequencer["id"], 0, wavefolder["id"], 0)
    assert has_cable(patch, arpeggio_sequencer["id"], 0, arpeggio_voice["id"], 2)
    assert has_cable(patch, arpeggio_sequencer["id"], 1, arpeggio_voice["id"], 6)
    assert has_cable(patch, arpeggio_sequencer["id"], 2, arpeggio_voice["id"], 7)
    assert has_cable(patch, arpeggio_sequencer["id"], 5, arpeggio_voice["id"], 3)
    assert has_cable(patch, arpeggio_sequencer["id"], 6, wavefolder["id"], 3)
    assert has_cable(patch, arpeggio_voice["id"], 1, scene_pack["id"], 1)

    assert has_cable(patch, scene_pack["id"], 0, reverb["id"], 0)
    assert has_cable(patch, reverb["id"], 0, audio["id"], 0)
    assert has_cable(patch, reverb["id"], 1, audio["id"], 1)
    assert has_cable(patch, scene_pack["id"], 0, scope["id"], 0)
    assert has_cable(patch, reverb["id"], 0, scope["id"], 1)
    assert has_cable(patch, reverb["id"], 1, scope["id"], 2)
    reverb_values = param_values(reverb)
    superlush = {
        0: 0.78,
        1: 0.55,
        2: 0.112,
        3: 0.50,
        4: 0.32,
        5: 0.50,
        6: 0.68,
        7: 0.84,
        8: 0.26,
        9: 1.0,
        10: 1.0,
        11: 0.75,
        12: 0.75,
        13: 0.4580135308,
        14: math.log(3.0) / math.log(20.0),
        15: 0.40,
        16: -6.0,
        17: 0.45,
    }
    for param_id, value in superlush.items():
        assert reverb_values[param_id] == pytest.approx(value)
    assert 1000.0 * math.exp(reverb_values[14] * math.log(20.0)) == pytest.approx(
        3000.0
    )
    for source in range(8):
        assert reverb_values[3 if source == 0 else 18 + 2 * (source - 1)] == 0.5
        assert reverb_values[4 if source == 0 else 19 + 2 * (source - 1)] == 0.32


def test_diode_ladder_default_cutoff_cv_range_reaches_fully_open():
    cutoff_pitch = EXPECTED_DEFAULTS["Tf303VoiceCore"][0]
    cv_amount = EXPECTED_DEFAULTS["Tf303VoiceCore"][5]
    cutoff_hz = 261.625565 * 2.0**cutoff_pitch
    envelope_peak_cutoff_hz = cutoff_hz * 2.0 ** (10.0 * cv_amount)

    assert cutoff_hz == pytest.approx(500.0, rel=1.0e-6)
    assert envelope_peak_cutoff_hz == pytest.approx(20_000.0, rel=1.0e-6)


def test_303_voice_foundry_song_reproduces_gibber_pattern_and_transpositions():
    foundry = modules(load_patch("303_voice"), "Foundry")[0]
    data = foundry["data"]
    assert param_values(foundry)[78] == 0.0
    assert data["id0_songBeginIndex"] == data["id1_songBeginIndex"] == 0
    assert data["id0_songEndIndex"] == data["id1_songEndIndex"] == 2
    assert data["id0_phrases"][:3] == data["id1_phrases"][:3] == [1792, 769, 770]
    assert data["id0_sequences"][:3] == data["id1_sequences"][:3] == [16, 16, 16]
    assert data["id0_seqSaved"][:3] == data["id1_seqSaved"][:3] == [1, 1, 1]

    expected_semitones = [
        [-21, -21, -21, -21, -14, -11, -21, -21, -21, -21, -9, -33, -21, -21, -33, -21],
        [
            -28,
            -28,
            -28,
            -28,
            -21,
            -18,
            -28,
            -28,
            -28,
            -28,
            -16,
            -40,
            -28,
            -28,
            -40,
            -28,
        ],
        [
            -26,
            -26,
            -26,
            -26,
            -19,
            -15,
            -26,
            -26,
            -26,
            -26,
            -14,
            -38,
            -26,
            -26,
            -38,
            -26,
        ],
    ]
    actual_semitones = [
        [round(volts * 12) for volts in data["id0_cv"][start : start + 16]]
        for start in (0, 32, 64)
    ]
    assert actual_semitones == expected_semitones

    for start in (0, 32, 64):
        note_attributes = data["id0_attributes"][start : start + 16]
        slide_attributes = data["id1_attributes"][start : start + 16]
        assert sum(bool(value & (1 << 24)) for value in note_attributes) == 13
        assert all((value >> 28) == 0 for value in note_attributes)
        assert {
            index
            for index, value in enumerate(note_attributes)
            if (value & 0xFF) == 200
        } == {0, 4, 6, 10, 14}
        assert {
            index for index, value in enumerate(slide_attributes) if value & (1 << 24)
        } == set(range(8, 16))


@pytest.mark.parametrize("name", PATCH_PATHS)
def test_smoke_patch_has_safe_output_gain(name):
    patch = load_patch(name)
    audio = modules(patch, "AudioInterface")[0]
    if name == "electric_piano":
        piano = modules(patch, "TfElectricPiano")[0]
        masters = modules(patch, "VCMixer")
        assert len(masters) == 2
        assert all(
            [param_values(master)[0], param_values(master)[1]] == [1.0, 1.0]
            for master in masters
        )
        piano_routes = [
            cable
            for cable in patch["cables"]
            if cable["outputModuleId"] == piano["id"]
            and cable["inputModuleId"] in {master["id"] for master in masters}
        ]
        assert {cable["outputId"] for cable in piano_routes} == {1, 2}
        assert {cable["inputId"] for cable in piano_routes} == {1}
        assert {cable["inputModuleId"] for cable in piano_routes} == {
            master["id"] for master in masters
        }
        routed_inputs = {
            cable["inputId"]
            for cable in patch["cables"]
            if cable["outputModuleId"] in {master["id"] for master in masters}
            and cable["inputModuleId"] == audio["id"]
        }
        assert routed_inputs == {0, 1}
        return
    if name in {
        "reverb_two_sources",
        "prog_303",
        "scene_pack",
    }:
        reverb = modules(patch, "TfReverb")[0]
        assert param_values(reverb)[16] == -6.0
        routed_inputs = {
            cable["inputId"]
            for cable in patch["cables"]
            if cable["outputModuleId"] == reverb["id"]
            and cable["inputModuleId"] == audio["id"]
        }
        assert routed_inputs == {0, 1}
        return
    masters = [
        mixer
        for mixer in modules(patch, "VCMixer")
        if math.isclose(20 * math.log10(param_values(mixer)[0]), -6.0, abs_tol=1e-6)
    ]

    routed_inputs = {
        cable["inputId"]
        for cable in patch["cables"]
        if cable["outputModuleId"] in {master["id"] for master in masters}
        and cable["inputModuleId"] == audio["id"]
    }
    assert routed_inputs == {0, 1}
