import json
import math
from pathlib import Path

PATCH_PATH = Path(__file__).parents[2] / "test.vcv"


def load_patch():
    return json.loads(PATCH_PATH.read_text(encoding="utf-8"))


def test_smoke_patch_contains_every_triggerfish_module():
    patch = load_patch()
    models = {
        module["model"]
        for module in patch["modules"]
        if module["plugin"] == "TriggerFish-Elements"
    }
    assert models == {"TfSlop", "TfSlop4", "TfVCA", "TfVDPO"}


def test_smoke_patch_uses_triggerfish_parameter_defaults():
    expected_defaults = {
        "TfSlop": {0: 0.25, 1: 0.05, 2: 1.0, 3: -1.0},
        "TfSlop4": {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 0.1, 5: 0.05, 6: 0.05},
        "TfVDPO": {0: 0.5, 1: 0.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0},
        "TfVCA": {0: 0.5, 1: 1.0, 2: 1.0, 3: 0.5, 4: 50.0, 5: 1.0},
    }

    triggerfish_modules = [
        module
        for module in load_patch()["modules"]
        if module["plugin"] == "TriggerFish-Elements"
    ]
    for module in triggerfish_modules:
        values = {param["id"]: param["value"] for param in module["params"]}
        assert values == expected_defaults[module["model"]]


def test_smoke_patch_dependencies_and_cables_are_valid():
    patch = load_patch()
    modules = patch["modules"]
    cables = patch["cables"]
    module_ids = {module["id"] for module in modules}

    assert len(module_ids) == len(modules)
    assert {module["plugin"] for module in modules} <= {
        "Core",
        "Fundamental",
        "TriggerFish-Elements",
    }
    assert len({cable["id"] for cable in cables}) == len(cables)
    assert all(cable["outputModuleId"] in module_ids for cable in cables)
    assert all(cable["inputModuleId"] in module_ids for cable in cables)


def test_smoke_patch_has_playable_midi_and_stereo_audio_paths():
    patch = load_patch()
    modules = patch["modules"]
    cables = patch["cables"]
    ids_by_model = {module["model"]: module["id"] for module in modules}

    assert {
        "MIDIToCVInterface",
        "ADSR",
        "AudioInterface",
    } <= ids_by_model.keys()
    assert "Scope" not in ids_by_model
    assert "Plateau" not in ids_by_model

    midi_id = ids_by_model["MIDIToCVInterface"]
    adsr_id = ids_by_model["ADSR"]
    audio_id = ids_by_model["AudioInterface"]
    assert any(
        cable["outputModuleId"] == midi_id
        and cable["outputId"] == 1
        and cable["inputModuleId"] == adsr_id
        for cable in cables
    )
    assert {
        cable["inputId"] for cable in cables if cable["inputModuleId"] == audio_id
    } == {0, 1}


def test_smoke_patch_does_not_name_local_hardware():
    modules = load_patch()["modules"]
    midi = next(module for module in modules if module["model"] == "MIDIToCVInterface")
    audio = next(module for module in modules if module["model"] == "AudioInterface")

    assert midi["data"]["midi"]["driver"] == -1
    assert "deviceName" not in midi["data"]["midi"]
    assert audio["data"]["audio"]["driver"] == -1
    assert "deviceName" not in audio["data"]["audio"]


def test_smoke_patch_uses_current_plugin_version():
    triggerfish_modules = [
        module
        for module in load_patch()["modules"]
        if module["plugin"] == "TriggerFish-Elements"
    ]
    assert {module["version"] for module in triggerfish_modules} == {"2.1.0"}


def test_smoke_patch_exercises_all_slop4_outputs():
    patch = load_patch()
    slop4_id = next(
        module["id"] for module in patch["modules"] if module["model"] == "TfSlop4"
    )
    connected_outputs = {
        cable["outputId"]
        for cable in patch["cables"]
        if cable["outputModuleId"] == slop4_id
    }
    assert connected_outputs == {0, 1, 2, 3}

    vco_pitch_offsets = sorted(
        next(param["value"] for param in module["params"] if param["id"] == 2)
        for module in patch["modules"]
        if module["model"] == "VCO"
    )
    assert vco_pitch_offsets == [0, 0, 12, 12]


def test_smoke_patch_has_quiet_master_and_no_reverb():
    patch = load_patch()
    modules = patch["modules"]

    final_mixer = next(module for module in modules if module["id"] == 14)
    master_gain = next(
        param["value"] for param in final_mixer["params"] if param["id"] == 0
    )
    assert math.isclose(20 * math.log10(master_gain), -6.0, abs_tol=1e-6)

    audio = next(module for module in modules if module["model"] == "AudioInterface")
    master_to_audio = {
        cable["inputId"]
        for cable in patch["cables"]
        if cable["outputModuleId"] == final_mixer["id"]
        and cable["inputModuleId"] == audio["id"]
    }
    assert master_to_audio == {0, 1}
