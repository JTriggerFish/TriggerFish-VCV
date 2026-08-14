"""Generate the checked-in Rack smoke-test patches.

Run from the repository root with::

    uv run python tools/generate_test_patches.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORE_VERSION = "2.6.6"
FUNDAMENTAL_VERSION = "2.6.4"
TRIGGERFISH_VERSION = json.loads((ROOT / "plugin.json").read_text(encoding="utf-8"))[
    "version"
]

COLORS = ("#f3374b", "#ffcc00", "#0c8e15", "#0986ad", "#c91847", "#7a3ec8")

DEFAULT_PARAMS = {
    "ADSR": [0.2, 0.35, 0.75, 0.45, 0.0, 0.0, 0.0, 0.0, 0.0],
    "TfSlop": [0.25, 0.05, 1.0, -1.0],
    "TfSlop4": [1.0, 1.0, 1.0, 1.0, 0.1, 0.05, 0.05],
    "TfVCA": [0.5, 1.0, 1.0, 0.5, 50.0, 1.0],
    "TfVDPO": [0.5, 0.0, 1.0, 1.0, 1.0, 1.0],
    "VCO": [1.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0],
}


def params(values: list[float]) -> list[dict[str, float | int]]:
    return [{"id": index, "value": value} for index, value in enumerate(values)]


def module(
    module_id: int,
    plugin: str,
    model: str,
    pos: tuple[int, int],
    *,
    values: list[float] | None = None,
    data: dict | None = None,
) -> dict:
    versions = {
        "Core": CORE_VERSION,
        "Fundamental": FUNDAMENTAL_VERSION,
        "TriggerFish-Elements": TRIGGERFISH_VERSION,
    }
    result = {
        "id": module_id,
        "plugin": plugin,
        "model": model,
        "version": versions[plugin],
        "params": params(
            values if values is not None else DEFAULT_PARAMS.get(model, [])
        ),
    }
    if data is not None:
        result["data"] = data
    result["pos"] = list(pos)
    return result


def notes(module_id: int, text: str, pos: tuple[int, int] = (0, 0)) -> dict:
    return module(module_id, "Core", "Notes", pos, data={"text": text})


def midi(module_id: int, pos: tuple[int, int]) -> dict:
    return module(
        module_id,
        "Core",
        "MIDIToCVInterface",
        pos,
        data={"channels": 1, "pwRange": 2, "midi": {"driver": -1, "channel": -1}},
    )


def audio(module_id: int, pos: tuple[int, int]) -> dict:
    return module(
        module_id,
        "Core",
        "AudioInterface",
        pos,
        data={
            "audio": {
                "driver": -1,
                "sampleRate": 48000,
                "blockSize": 256,
                "inputOffset": 0,
                "outputOffset": 0,
            },
            "dcFilter": True,
        },
    )


def mixer(
    module_id: int,
    pos: tuple[int, int],
    levels: tuple[float, float, float, float, float],
) -> dict:
    return module(
        module_id,
        "Fundamental",
        "VCMixer",
        pos,
        values=list(levels),
        data={"chExp": False, "mixExp": False},
    )


class Patch:
    def __init__(
        self, zoom: float = 0.72, grid_offset: tuple[float, float] = (-1, -0.1)
    ):
        self.modules: list[dict] = []
        self.cables: list[dict] = []
        self.zoom = zoom
        self.grid_offset = grid_offset

    def add(self, item: dict) -> None:
        self.modules.append(item)

    def cable(self, source: int, output: int, target: int, input_: int) -> None:
        cable_id = 1000 + len(self.cables)
        self.cables.append(
            {
                "id": cable_id,
                "outputModuleId": source,
                "outputId": output,
                "inputModuleId": target,
                "inputId": input_,
                "color": COLORS[len(self.cables) % len(COLORS)],
            }
        )

    def write(self, filename: str) -> None:
        document = {
            "version": CORE_VERSION,
            "path": "",
            "zoom": self.zoom,
            "gridOffset": list(self.grid_offset),
            "modules": self.modules,
            "cables": self.cables,
        }
        with (ROOT / filename).open("w", encoding="utf-8", newline="\n") as output:
            output.write(json.dumps(document, indent=2) + "\n")


def generate_slop4_patch() -> None:
    patch = Patch()
    patch.add(
        notes(
            1,
            "TriggerFish Slop4 playable smoke test\n\n"
            "SETUP\nSelect MIDI and audio devices in MIDI-CV and Audio-8. "
            "Start with monitor volume low.\n\n"
            "MIDI pitch -> Slop4 -> four Fundamental VCOs -> VCA Mix -> "
            "TriggerFish VCA. MIDI gate drives the ADSR. Two oscillators are "
            "at unison and two are one octave higher.\n\n"
            "All TriggerFish parameters are at their declared defaults. The "
            "final VCA Mix is the master, set to -6 dB and routed to outputs 1/2.",
        )
    )
    patch.add(midi(2, (16, 0)))
    patch.add(module(3, "Fundamental", "ADSR", (32, 0)))
    patch.add(mixer(11, (68, 0), (0.5011872336, 0.7, 0.0, 0.0, 0.0)))
    patch.add(audio(12, (77, 0)))
    patch.add(module(4, "TriggerFish-Elements", "TfSlop4", (0, 1)))
    for module_id, x, octave, fine in (
        (5, 10, 0.0, 0.5),
        (6, 19, 0.0, 0.5),
        (7, 28, 12.0, 0.5),
        (8, 37, 12.0, 0.5),
    ):
        values = DEFAULT_PARAMS["VCO"].copy()
        values[2] = octave
        values[5] = fine
        patch.add(module(module_id, "Fundamental", "VCO", (x, 1), values=values))
    patch.add(mixer(9, (46, 1), (0.7, 0.35, 0.35, 0.3, 0.25)))
    patch.add(module(10, "TriggerFish-Elements", "TfVCA", (55, 1)))

    for channel, vco_id in enumerate((5, 6, 7, 8)):
        patch.cable(2, 0, 4, channel)
        patch.cable(4, channel, vco_id, 0)
        patch.cable(vco_id, channel, 9, channel + 1)
    patch.cable(2, 1, 3, 4)
    patch.cable(3, 0, 10, 1)
    patch.cable(9, 0, 10, 0)
    patch.cable(10, 0, 11, 1)
    patch.cable(11, 0, 12, 0)
    patch.cable(11, 0, 12, 1)
    patch.write("test-slop4.vcv")


def generate_vdpo_patch() -> None:
    patch = Patch(zoom=0.72, grid_offset=(-1, -0.1))
    patch.add(
        notes(
            1,
            "TriggerFish VDPO playable smoke test\n\n"
            "SETUP\nSelect MIDI and audio devices in MIDI-CV and Audio-8. "
            "Start with monitor volume low.\n\n"
            "VOICE 1 is a self-resonating Slop -> VDPO -> VCA chain. VOICE 2 "
            "is the same chain with its VDPO forced by a Fundamental VCO. Both "
            "VCAs have independent gate-driven ADSR envelopes.\n\n"
            "Click PUSH to cycle the forcing waveform through SINE, SAW, and "
            "SQUARE; the switch light shows the selected input. All TriggerFish "
            "parameters are at their declared defaults. The final VCA Mix is "
            "the master, set to -6 dB and routed to outputs 1/2.",
        )
    )
    patch.add(midi(2, (16, 0)))
    patch.add(
        module(
            13, "Fundamental", "Push", (32, 0), values=[0.0, 0.0], data={"hold": False}
        )
    )
    patch.add(
        module(
            12,
            "Fundamental",
            "SequentialSwitch2",
            (38, 0),
            values=[1.0],
            data={"declick": True},
        )
    )
    patch.add(module(11, "Fundamental", "VCO", (47, 0)))
    patch.add(mixer(14, (68, 0), (0.5011872336, 0.7, 0.7, 0.0, 0.0)))
    patch.add(audio(15, (77, 0)))

    patch.add(module(3, "Fundamental", "ADSR", (0, 1)))
    patch.add(module(4, "TriggerFish-Elements", "TfSlop", (10, 1)))
    patch.add(module(5, "TriggerFish-Elements", "TfVDPO", (16, 1)))
    patch.add(module(6, "TriggerFish-Elements", "TfVCA", (24, 1)))
    patch.add(module(7, "Fundamental", "ADSR", (0, 2)))
    patch.add(module(8, "TriggerFish-Elements", "TfSlop", (10, 2)))
    patch.add(module(9, "TriggerFish-Elements", "TfVDPO", (25, 2)))
    patch.add(module(10, "TriggerFish-Elements", "TfVCA", (33, 2)))

    patch.cable(2, 0, 4, 0)
    patch.cable(4, 0, 5, 0)
    patch.cable(2, 1, 3, 4)
    patch.cable(3, 0, 6, 1)
    patch.cable(5, 0, 6, 0)

    patch.cable(2, 0, 8, 0)
    patch.cable(8, 0, 9, 0)
    patch.cable(8, 0, 11, 0)
    patch.cable(2, 1, 7, 4)
    patch.cable(7, 0, 10, 1)
    patch.cable(9, 0, 10, 0)

    patch.cable(13, 0, 12, 0)
    for output, switch_input in ((0, 2), (2, 3), (3, 4)):
        patch.cable(11, output, 12, switch_input)
    patch.cable(12, 0, 9, 1)

    patch.cable(6, 0, 14, 1)
    patch.cable(10, 0, 14, 2)
    patch.cable(14, 0, 15, 0)
    patch.cable(14, 0, 15, 1)
    patch.write("test-vdpo.vcv")


if __name__ == "__main__":
    generate_slop4_patch()
    generate_vdpo_patch()
