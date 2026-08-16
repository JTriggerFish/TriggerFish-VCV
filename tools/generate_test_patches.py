"""Generate the checked-in Rack smoke-test patches.

Run from the repository root with::

    uv run python tools/generate_test_patches.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORE_VERSION = "2.6.6"
FUNDAMENTAL_VERSION = "2.6.4"
IMPROMPTU_VERSION = "2.5.0"
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
    "Tf303VoiceCore": [
        0.9344246,
        0.0,
        0.0,
        0.0,
        0.0,
        0.5321928,
        0.0,
        0.0,
        1.0 / 3.0,
        math.log10(0.5),
        math.log10(0.2),
        0.5,
        0.5,
        1.0,
        2.0,
    ],
    "Tf303Oscillator": [
        0.0,
        0.0,
        math.log10(0.060),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    "Tf4072VoiceCore": [
        0.9344246,
        0.0,
        0.0,
        0.6,
        0.0,
        0.0,
        0.0,
        1.0,
        1.0,
        math.log10(0.010),
        math.log10(0.300),
        0.5,
        math.log10(0.300),
        0.0,
        0.0,
        math.log10(0.005),
        math.log10(0.300),
        1.0,
        math.log10(0.300),
        1.0,
        0.0,
    ],
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
        "ImpromptuModular": IMPROMPTU_VERSION,
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


def foundry_step(*, gate: bool = False, velocity: int = 0, gate_type: int = 0) -> int:
    """Pack Foundry's public per-step fields into its patch representation."""
    value = velocity | (50 << 8) | (10 << 16) | (gate_type << 28)
    if gate:
        value |= 1 << 24
    return value


def foundry_track(
    track: int, pitches: list[float], attributes: list[int]
) -> dict[str, int | list[int] | list[float]]:
    """Build one 16-step Foundry track, padded to its 32-step storage size."""
    if len(pitches) != 16 or len(attributes) != 16:
        raise ValueError("Foundry smoke-test tracks must contain 16 steps")
    prefix = f"id{track}_"
    return {
        f"{prefix}pulsesPerStep": 1,
        f"{prefix}sequences": [16] + [32] * 63,
        f"{prefix}seqSaved": [1] + [0] * 63,
        f"{prefix}cv": pitches + [0.0] * 16,
        f"{prefix}attributes": attributes + [foundry_step() for _ in range(16)],
        f"{prefix}seqIndexEdit": 0,
    }


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


def generate_303_voice_patch() -> None:
    patch = Patch(zoom=0.68, grid_offset=(-1, -0.1))
    pitches = [
        -20 / 12,
        -20 / 12,
        -20 / 12,
        -20 / 12,
        -20 / 12,
        -8 / 12,
        -20 / 12,
        -20 / 12,
        -20 / 12,
        -20 / 12,
        -20 / 12,
        -12 / 12,
        -24 / 12,
        -20 / 12,
        -20 / 12,
        -20 / 12,
    ]
    slide_steps = {5, 12}
    legato_steps = slide_steps | {step - 1 for step in slide_steps}
    rest_steps = {1, 3, 6, 8, 10, 13, 15}
    accent_steps = {2, 9, 14}
    note_attributes = [
        foundry_step(
            gate=step not in rest_steps,
            velocity=200 if step in accent_steps else 0,
            gate_type=5 if step in legato_steps else 0,
        )
        for step in range(16)
    ]
    # Foundry track B is a full-step gate lane for our oscillator's Slide input.
    slide_attributes = [
        foundry_step(gate=step in slide_steps, gate_type=5) for step in range(16)
    ]
    foundry_data: dict = {
        "velocityMode": 0,
        "running": True,
        "stepIndexEdit": 0,
        "trackIndexEdit": 0,
    }
    foundry_data.update(foundry_track(0, pitches, note_attributes))
    foundry_data.update(foundry_track(1, [0.0] * 16, slide_attributes))

    patch.add(
        notes(
            1,
            "TriggerFish 303 voice smoke test\n\n"
            "SETUP\nSelect an output device in Audio-8 and start with monitor "
            "volume low. Impromptu Clocked and Foundry run the programmed "
            "16-step pattern automatically.\n\n"
            "Foundry A sends pitch, gate, and per-step accent. Track B marks "
            "the legato notes for Tf303Oscillator's Slide input. The oscillator's "
            "post-slide CV tracks 303 Voice Core cutoff, and VCA OUT feeds "
            "the -6 dB stereo master. Isolated repeated E notes alternate normal "
            "and accented articulation; two legato pairs slide by one octave.\n\n"
            "All TriggerFish parameters start at their declared defaults. Try "
            "WAVE and SHAPE on the oscillator, then CUTOFF, RES, ENV, ACCENT, "
            "DRIVE, and BASS on the filter.",
        )
    )
    patch.add(
        module(
            2,
            "ImpromptuModular",
            "Clocked",
            (16, 0),
            values=[125.0, 5.0, 0.0, 0.0] + [0.0] * 4 + [0.5] * 4 + [0.0] * 8,
            data={"running": True, "ppqn": 4},
        )
    )
    patch.add(
        module(
            3,
            "ImpromptuModular",
            "Foundry",
            (28, 0),
            data=foundry_data,
        )
    )
    patch.add(mixer(6, (65, 0), (0.5011872336, 0.7, 0.0, 0.0, 0.0)))
    patch.add(audio(7, (74, 0)))
    patch.add(module(4, "TriggerFish-Elements", "Tf303Oscillator", (28, 1)))
    patch.add(module(5, "TriggerFish-Elements", "Tf303VoiceCore", (40, 1)))

    patch.cable(2, 1, 3, 6)  # Clocked x4 -> Foundry track A clock
    patch.cable(2, 4, 3, 5)  # reset
    patch.cable(3, 0, 4, 0)  # track A pitch -> oscillator V/OCT
    patch.cable(3, 9, 4, 1)  # track B gate -> oscillator Slide
    patch.cable(4, 0, 5, 1)  # post-slide pitch -> filter V/OCT
    patch.cable(4, 1, 5, 0)  # oscillator audio -> filter input
    patch.cable(3, 8, 5, 5)  # track A gate -> filter Gate
    patch.cable(3, 4, 5, 6)  # track A CV2 -> filter Accent
    patch.cable(5, 1, 6, 1)
    patch.cable(6, 0, 7, 0)
    patch.cable(6, 0, 7, 1)
    patch.write("test-303-voice.vcv")


def generate_4072_voice_patch() -> None:
    patch = Patch(zoom=0.75, grid_offset=(-1, -0.1))
    patch.add(
        notes(
            1,
            "TriggerFish 4072 Voice Core smoke test\n\n"
            "Select MIDI and audio devices, then play from a keyboard. A "
            "Fundamental saw oscillator feeds the 4072 filter and normalled "
            "4019 VCA. MIDI pitch tracks both oscillator and filter; MIDI "
            "gate drives the two internal envelopes.\n\n"
            "All TriggerFish parameters are at their declared defaults. The "
            "final mixer is the master level, set to -6 dB and routed to "
            "outputs 1/2.",
        )
    )
    patch.add(midi(2, (16, 0)))
    patch.add(module(3, "Fundamental", "VCO", (28, 0)))
    patch.add(module(4, "TriggerFish-Elements", "Tf4072VoiceCore", (40, 0)))
    patch.add(mixer(5, (65, 0), (0.5011872336, 0.7, 0.0, 0.0, 0.0)))
    patch.add(audio(6, (74, 0)))

    patch.cable(2, 0, 3, 0)  # MIDI pitch -> oscillator
    patch.cable(2, 0, 4, 2)  # MIDI pitch -> filter 1V/oct
    patch.cable(2, 1, 4, 8)  # MIDI gate -> both internal envelopes
    patch.cable(3, 2, 4, 0)  # saw -> filter input
    patch.cable(4, 1, 5, 1)  # VCA output -> master
    patch.cable(5, 0, 6, 0)
    patch.cable(5, 0, 6, 1)
    patch.write("test-4072-voice.vcv")


if __name__ == "__main__":
    generate_slop4_patch()
    generate_vdpo_patch()
    generate_303_voice_patch()
    generate_4072_voice_patch()
