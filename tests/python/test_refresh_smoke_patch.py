import json
import sys
from pathlib import Path

ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT))
from tools.refresh_smoke_patch import refresh_local_patch


def write_patch(path, modules):
    path.write_text(json.dumps({"modules": modules}), encoding="utf-8")


def test_refresh_creates_local_patch_from_portable_patch(tmp_path):
    portable_path = tmp_path / "test.vcv"
    local_path = tmp_path / "test.local.vcv"
    write_patch(portable_path, [{"id": 1, "plugin": "Core", "model": "Notes"}])

    assert refresh_local_patch(portable_path, local_path) is False
    assert json.loads(local_path.read_text(encoding="utf-8")) == json.loads(
        portable_path.read_text(encoding="utf-8")
    )


def test_refresh_updates_patch_but_retains_only_device_state(tmp_path):
    portable_path = tmp_path / "test.vcv"
    local_path = tmp_path / "test.local.vcv"
    portable_modules = [
        {
            "id": 10,
            "plugin": "Core",
            "model": "MIDIToCVInterface",
            "params": [{"id": 0, "value": 0.25}],
            "data": {
                "channels": 1,
                "midi": {"driver": -1, "channel": -1},
            },
        },
        {
            "id": 20,
            "plugin": "Core",
            "model": "AudioInterface",
            "data": {
                "audio": {"driver": -1, "sampleRate": 48_000},
                "dcFilter": True,
            },
        },
        {"id": 3, "plugin": "Core", "model": "Notes", "data": {"text": "new"}},
    ]
    local_modules = [
        {
            "id": 1,
            "plugin": "Core",
            "model": "MIDIToCVInterface",
            "params": [{"id": 0, "value": 0.75}],
            "data": {
                "channels": 16,
                "midi": {"driver": 7, "deviceName": "TEST_MIDI_DEVICE"},
            },
        },
        {
            "id": 2,
            "plugin": "Core",
            "model": "AudioInterface",
            "data": {
                "audio": {"driver": 8, "deviceName": "TEST_AUDIO_DEVICE"},
                "dcFilter": False,
            },
        },
        {"id": 3, "plugin": "Core", "model": "Notes", "data": {"text": "old"}},
    ]
    write_patch(portable_path, portable_modules)
    write_patch(local_path, local_modules)

    assert refresh_local_patch(portable_path, local_path) is True
    refreshed = json.loads(local_path.read_text(encoding="utf-8"))["modules"]

    assert refreshed[0]["params"] == [{"id": 0, "value": 0.25}]
    assert refreshed[0]["data"] == {
        "channels": 1,
        "midi": {"driver": 7, "deviceName": "TEST_MIDI_DEVICE"},
    }
    assert refreshed[1]["data"] == {
        "audio": {"driver": 8, "deviceName": "TEST_AUDIO_DEVICE"},
        "dcFilter": True,
    }
    assert refreshed[2]["data"] == {"text": "new"}
