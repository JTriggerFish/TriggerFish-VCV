import json
import sys
from pathlib import Path

ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT))
from tools.refresh_smoke_patch import prepare_local_patch


def test_prepare_creates_exact_local_copy(tmp_path):
    portable_path = tmp_path / "test.vcv"
    local_path = tmp_path / "test.local.vcv"
    portable_data = b'{"modules": []}\n'
    portable_path.write_bytes(portable_data)

    assert prepare_local_patch(portable_path, local_path) is False
    assert local_path.read_bytes() == portable_data


def test_prepare_preserves_existing_compressed_local_patch(tmp_path):
    portable_path = tmp_path / "test.vcv"
    local_path = tmp_path / "test.local.vcv"
    portable_path.write_bytes(b'{"modules": []}\n')
    compressed_local_data = b"\x28\xb5\x2f\xfd\x00user Rack archive"
    local_path.write_bytes(compressed_local_data)

    assert prepare_local_patch(portable_path, local_path) is True
    assert local_path.read_bytes() == compressed_local_data


def test_new_local_patch_can_inherit_private_device_selections(tmp_path):
    portable = tmp_path / "portable.vcv"
    local = tmp_path / "local.vcv"
    template = tmp_path / "template.local.vcv"
    portable.write_text(
        json.dumps(
            {
                "modules": [
                    {
                        "model": "MIDIToCVInterface",
                        "data": {"midi": {"driver": -1, "channel": -1}},
                    },
                    {
                        "model": "MIDICCToCVInterface",
                        "data": {"midi": {"driver": -1, "channel": -1}},
                    },
                    {
                        "model": "AudioInterface",
                        "data": {"audio": {"driver": -1, "blockSize": 256}},
                    },
                ]
            }
        ),
        encoding="utf-8",
    )
    template.write_text(
        json.dumps(
            {
                "modules": [
                    {
                        "model": "MIDIToCVInterface",
                        "data": {
                            "midi": {
                                "driver": 4,
                                "deviceName": "Private keyboard",
                                "channel": -1,
                            }
                        },
                    },
                    {
                        "model": "AudioInterface",
                        "data": {
                            "audio": {
                                "driver": 7,
                                "deviceName": "Private interface",
                                "blockSize": 64,
                            }
                        },
                    },
                ]
            }
        ),
        encoding="utf-8",
    )

    assert prepare_local_patch(portable, local, template) is False
    created = json.loads(local.read_text(encoding="utf-8"))
    assert created["modules"][0]["data"]["midi"]["deviceName"] == ("Private keyboard")
    assert created["modules"][1]["data"]["midi"]["deviceName"] == ("Private keyboard")
    assert created["modules"][2]["data"]["audio"]["deviceName"] == ("Private interface")


def test_refresh_replaces_topology_but_preserves_local_device_selection(tmp_path):
    portable = tmp_path / "portable.vcv"
    local = tmp_path / "local.vcv"
    portable.write_text(
        json.dumps(
            {
                "modules": [
                    {"model": "TfProgSequencer", "data": {"source": "new"}},
                    {
                        "model": "AudioInterface",
                        "data": {"audio": {"driver": -1, "blockSize": 256}},
                    },
                ]
            }
        ),
        encoding="utf-8",
    )
    local.write_text(
        json.dumps(
            {
                "modules": [
                    {"model": "Foundry", "data": {"old": True}},
                    {
                        "model": "AudioInterface",
                        "data": {
                            "audio": {
                                "driver": 7,
                                "deviceName": "Private interface",
                                "blockSize": 64,
                            }
                        },
                    },
                ]
            }
        ),
        encoding="utf-8",
    )

    assert prepare_local_patch(portable, local, refresh=True) is False
    refreshed = json.loads(local.read_text(encoding="utf-8"))
    assert [module["model"] for module in refreshed["modules"]] == [
        "TfProgSequencer",
        "AudioInterface",
    ]
    assert refreshed["modules"][1]["data"]["audio"]["deviceName"] == (
        "Private interface"
    )


def test_refresh_does_not_replace_selected_devices_with_disconnected_local_state(
    tmp_path,
):
    portable = tmp_path / "portable.vcv"
    local = tmp_path / "local.vcv"
    portable.write_text(
        json.dumps(
            {
                "modules": [
                    {
                        "model": "MIDIToCVInterface",
                        "data": {
                            "midi": {
                                "driver": 4,
                                "deviceName": "Checked-in keyboard",
                                "channel": -1,
                            }
                        },
                    },
                    {
                        "model": "AudioInterface",
                        "data": {
                            "audio": {
                                "driver": 7,
                                "deviceName": "Checked-in interface",
                                "blockSize": 64,
                            }
                        },
                    },
                ]
            }
        ),
        encoding="utf-8",
    )
    local.write_text(
        json.dumps(
            {
                "modules": [
                    {
                        "model": "MIDIToCVInterface",
                        "data": {"midi": {"driver": -1, "channel": -1}},
                    },
                    {
                        "model": "AudioInterface",
                        "data": {"audio": {"driver": -1, "blockSize": 256}},
                    },
                ]
            }
        ),
        encoding="utf-8",
    )

    assert prepare_local_patch(portable, local, refresh=True) is False
    refreshed = json.loads(local.read_text(encoding="utf-8"))
    assert refreshed["modules"][0]["data"]["midi"]["deviceName"] == (
        "Checked-in keyboard"
    )
    assert refreshed["modules"][1]["data"]["audio"]["deviceName"] == (
        "Checked-in interface"
    )
