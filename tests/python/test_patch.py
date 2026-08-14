import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[2]
PLUGIN_VERSION = json.loads((ROOT / "plugin.json").read_text(encoding="utf-8"))[
    "version"
]
PATCH_PATHS = {
    "slop4": ROOT / "test-slop4.vcv",
    "vdpo": ROOT / "test-vdpo.vcv",
}

EXPECTED_DEFAULTS = {
    "TfSlop": {0: 0.25, 1: 0.05, 2: 1.0, 3: -1.0},
    "TfSlop4": {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 0.1, 5: 0.05, 6: 0.05},
    "TfVDPO": {0: 0.5, 1: 0.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0},
    "TfVCA": {0: 0.5, 1: 1.0, 2: 1.0, 3: 0.5, 4: 50.0, 5: 1.0},
}


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
    assert models == {"TfSlop", "TfSlop4", "TfVCA", "TfVDPO"}


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
        "TriggerFish-Elements",
    }
    assert len({cable["id"] for cable in patch["cables"]}) == len(patch["cables"])
    assert all(cable["outputModuleId"] in module_ids for cable in patch["cables"])
    assert all(cable["inputModuleId"] in module_ids for cable in patch["cables"])


@pytest.mark.parametrize("name", PATCH_PATHS)
def test_smoke_patch_has_playable_midi_and_stereo_audio_paths(name):
    patch = load_patch(name)
    midi = modules(patch, "MIDIToCVInterface")[0]
    audio = modules(patch, "AudioInterface")[0]
    adsrs = modules(patch, "ADSR")

    assert "Scope" not in {module["model"] for module in patch["modules"]}
    assert "Plateau" not in {module["model"] for module in patch["modules"]}
    assert all(has_cable(patch, midi["id"], 1, adsr["id"], 4) for adsr in adsrs)
    assert {
        cable["inputId"]
        for cable in patch["cables"]
        if cable["inputModuleId"] == audio["id"]
    } == {0, 1}


@pytest.mark.parametrize("name", PATCH_PATHS)
def test_smoke_patch_does_not_name_local_hardware(name):
    patch = load_patch(name)
    midi = modules(patch, "MIDIToCVInterface")[0]
    audio = modules(patch, "AudioInterface")[0]

    assert midi["data"]["midi"]["driver"] == -1
    assert "deviceName" not in midi["data"]["midi"]
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
