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
        math.log10(0.0014),
        math.log10(1.0),
        0.5,
        math.log10(1.0),
        1.0,
        0.0,
        math.log10(0.0014),
        math.log10(1.0),
        1.0,
        math.log10(1.0),
        1.0,
        1.0,
        1.0,
        1.0,
        0.0,
        1.0,
    ],
    "TfElectricPiano": [
        0.5,
        1.0,
        0.62,
        0.52,
        0.50,
        0.52,
        0.55,
        0.48,
        0.50,
        0.24,
        0.18,
        0.32,
    ],
    "TfWavefoldOscillator": [
        0.0,
        0.0,
        0.5,
        0.4,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        2.0,
        0.5,
        0.5,
        0.5,
        0.5,
        1.0,
        0.39841330778621553,
    ],
    "TfUnisonOscillator": [
        0.0,
        0.0,
        3.0,
        0.0,
        0.5,
        0.0,
        0.42,
        0.39841330778621553,
        0.65,
        0.0,
        0.0,
        0.5,
        0.1,
        0.1,
        0.15,
        1.0,
        1.0,
        0.0,
        0.0,
    ],
    "TfScenePack4": [],
    "TfTransport": [120.0, 0.0, 0.0, 0.0, 0.0],
    "TfReverb": [
        0.60,
        0.5,
        0.048,
        0.5,
        0.35,
        0.5,
        0.682,
        0.64,
        0.28,
        0.82,
        0.30,
        0.6130368568946039,
        0.5,
        0.0,
        0.829482217661603,
        0.35,
        0.0,
        0.0,
        0.50,
        0.35,
        0.50,
        0.35,
        0.50,
        0.35,
        0.50,
        0.35,
        0.50,
        0.35,
        0.50,
        0.35,
        0.50,
        0.35,
    ],
    "TfProgSequencer": [],
    "LFO": [0.0, 0.0, math.log2(0.07), 0.0, 0.0, 0.5, 0.0],
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


def midi(module_id: int, pos: tuple[int, int], *, channels: int = 1) -> dict:
    return module(
        module_id,
        "Core",
        "MIDIToCVInterface",
        pos,
        data={
            "channels": channels,
            "pwRange": 2,
            "midi": {"driver": -1, "channel": -1},
        },
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


def scope(module_id: int, pos: tuple[int, int]) -> dict:
    """Fundamental Scope configured for a triggered 10 ms display."""
    return module(
        module_id,
        "Fundamental",
        "Scope",
        pos,
        values=[0.0, 0.0, 0.0, 0.0, -math.log2(0.010), 0.0, 0.0, 0.0],
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
    track: int,
    pitch_sequences: list[list[float]],
    attribute_sequences: list[list[int]],
    phrase_repetitions: list[int],
) -> dict[str, int | list[int] | list[float]]:
    """Build a Foundry song from 16- or 32-step sequences."""
    if not pitch_sequences or len(pitch_sequences) != len(attribute_sequences):
        raise ValueError("Foundry pitch and attribute sequences must correspond")
    if len(phrase_repetitions) != len(pitch_sequences):
        raise ValueError("Each Foundry song phrase needs a repetition count")
    sequence_lengths = [len(sequence) for sequence in pitch_sequences]
    if any(length not in (16, 32) for length in sequence_lengths):
        raise ValueError("Foundry sequences must contain 16 or 32 steps")
    if any(
        len(attributes) != length
        for attributes, length in zip(attribute_sequences, sequence_lengths)
    ):
        raise ValueError("Foundry attributes must match their sequence length")

    saved_count = len(pitch_sequences)
    saved_cv: list[float] = []
    saved_attributes: list[int] = []
    for pitches, attributes in zip(pitch_sequences, attribute_sequences):
        padding = 32 - len(pitches)
        saved_cv.extend(pitches + [0.0] * padding)
        saved_attributes.extend(attributes + [foundry_step() for _ in range(padding)])

    phrases = [
        sequence + ((repetitions - 1) << 8)
        for sequence, repetitions in enumerate(phrase_repetitions)
    ]
    prefix = f"id{track}_"
    return {
        f"{prefix}pulsesPerStep": 1,
        f"{prefix}runModeSong": 0,
        f"{prefix}songBeginIndex": 0,
        f"{prefix}songEndIndex": len(phrases) - 1,
        f"{prefix}phrases": phrases + [0] * (99 - len(phrases)),
        f"{prefix}sequences": sequence_lengths + [32] * (64 - saved_count),
        f"{prefix}seqSaved": [1] * saved_count + [0] * (64 - saved_count),
        f"{prefix}cv": saved_cv,
        f"{prefix}attributes": saved_attributes,
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

    def cable(
        self,
        source: int,
        output: int,
        target: int,
        input_: int,
        *,
        color: str | None = None,
    ) -> None:
        cable_id = 1000 + len(self.cables)
        self.cables.append(
            {
                "id": cable_id,
                "outputModuleId": source,
                "outputId": output,
                "inputModuleId": target,
                "inputId": input_,
                "color": color or COLORS[len(self.cables) % len(COLORS)],
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
    patch = Patch(
        zoom=1.3599998950958252,
        grid_offset=(-12.460025787353516, -0.22162829339504242),
    )
    # The Gibber Acid playground pattern by fasttriggerfish and thecharlie.
    # Theory.degree changes both offset and mode: i is minor/Aeolian, -iv is
    # minor/Aeolian down a seventh, and -V is major/Ionian down a fifth.
    gibber_degrees = [0, 0, 0, 0, 4, 6, 0, None, 0, None, 7, -7, None, 0, -7, 0]
    root_d_sharp_2 = -21
    sections = [
        ({0: 0, 4: 7, 6: 10, 7: 12, -7: -12}, 0),
        ({0: 0, 4: 7, 6: 10, 7: 12, -7: -12}, -7),
        ({0: 0, 4: 7, 6: 11, 7: 12, -7: -12}, -5),
    ]
    pitch_sequences = [
        [
            (root_d_sharp_2 + degree_offset + mode.get(degree, 0)) / 12
            for degree in gibber_degrees
        ]
        for mode, degree_offset in sections
    ]
    # Gibber advances [1, 1, 100, 100] ms once per quarter note. The 303
    # oscillator's 2 ms floor makes the first half effectively instantaneous;
    # track B enables its 100 ms glide throughout the second half of each bar.
    slide_steps = set(range(8, 16))
    rest_steps = {step for step, degree in enumerate(gibber_degrees) if degree is None}
    accent_steps = {0, 4, 6, 10, 14}
    note_attributes = [
        foundry_step(
            gate=step not in rest_steps,
            velocity=200 if step in accent_steps else 0,
        )
        for step in range(len(gibber_degrees))
    ]
    # Foundry track B is a full-step gate lane for our oscillator's Slide input.
    slide_attributes = [
        foundry_step(gate=step in slide_steps, gate_type=5)
        for step in range(len(gibber_degrees))
    ]
    foundry_data: dict = {
        "velocityMode": 0,
        "running": True,
        "stepIndexEdit": 0,
        "phraseIndexEdit": 0,
        "trackIndexEdit": 0,
    }
    phrase_repetitions = [8, 4, 4]
    foundry_data.update(
        foundry_track(
            0,
            pitch_sequences,
            [note_attributes] * len(pitch_sequences),
            phrase_repetitions,
        )
    )

    # Retain the original pre-56ed864 smoke-patch parameter state while the
    # newer Acid arrangement exercises the current sequencer and reverb path.
    voice_params = DEFAULT_PARAMS["Tf303VoiceCore"].copy()
    oscillator_params = DEFAULT_PARAMS["Tf303Oscillator"].copy()
    foundry_data.update(
        foundry_track(
            1,
            [[0.0] * 16 for _ in pitch_sequences],
            [slide_attributes] * len(pitch_sequences),
            phrase_repetitions,
        )
    )

    patch.add(
        notes(
            1,
            "TriggerFish 303 voice smoke test\n\n"
            "SETUP\nSelect an output device in Audio-8 and start with monitor "
            "volume low. Impromptu Clocked and Foundry run the "
            "full 16-bar "
            "Gibber Acid playground line by fasttriggerfish and thecharlie at "
            "120 BPM: i for 8 bars, -iv for 4, then -V for 4.\n\n"
            "Foundry track A contains the notes, rests, and accents. Track B "
            "reproduces the repeating 1, 1, 100, 100 ms glide cycle; the first "
            "two quarters are effectively instantaneous. Independent slow "
            "bipolar sines modulate cutoff and resonance. A separate 14 Hz "
            "sine is connected "
            "to linear filter FM with its attenuverter at zero. VCA OUT feeds "
            "the -6 dB stereo master.\n\n"
            "The 303 oscillator and voice retain the original smoke patch's "
            "factory parameter state. Try "
            "WAVE and SHAPE on the oscillator, then CUTOFF, RES, ENV, ACCENT, "
            "DRIVE, and BASS on the filter.",
            pos=(-9, 0),
        ),
    )
    patch.add(
        module(
            2,
            "ImpromptuModular",
            "Clocked",
            (7, 0),
            values=[120.0, 5.0, 0.0, 0.0] + [0.0] * 4 + [0.5] * 4 + [0.0] * 8,
            data={"running": True, "ppqn": 4},
        )
    )
    foundry = module(
        3,
        "ImpromptuModular",
        "Foundry",
        (27, 0),
        data=foundry_data,
    )
    foundry["params"] = [{"id": 78, "value": 0.0}]  # Song mode
    patch.add(foundry)
    patch.add(
        mixer(
            6,
            (65, 0),
            (0.5011872336, 0.838592648506165, 0.0, 0.0, 0.0),
        )
    )
    patch.add(audio(7, (74, 0)))
    patch.add(
        module(
            4,
            "TriggerFish-Elements",
            "Tf303Oscillator",
            (5, 1),
            values=oscillator_params,
        )
    )
    patch.add(
        module(
            5,
            "TriggerFish-Elements",
            "Tf303VoiceCore",
            (30, 1),
            values=voice_params,
        )
    )
    patch.add(
        module(
            8,
            "Fundamental",
            "LFO",
            (46, 1),
            values=[0.0, 0.0, -4.00999546051025, 0.0, 0.0, 0.5, 0.0],
        )
    )
    patch.add(
        module(
            9,
            "Fundamental",
            "LFO",
            (55, 1),
            values=[0.0, 0.0, -4.1243371963501, 0.0, 0.0, 0.5, 0.0],
        )
    )
    patch.add(
        module(
            10,
            "Fundamental",
            "LFO",
            (64, 1),
            values=[0.0, 0.0, 3.80241107940674, 0.0, 0.0, 0.5, 0.0],
        )
    )

    patch.cable(2, 1, 3, 6)  # Clocked x4 -> Foundry track A clock
    patch.cable(2, 4, 3, 5)  # reset
    patch.cable(3, 0, 4, 0)  # track A pitch -> oscillator V/OCT
    patch.cable(3, 9, 4, 1)  # track B gate -> oscillator Slide
    patch.cable(4, 1, 5, 0)  # oscillator audio -> filter input
    patch.cable(8, 0, 5, 2, color="#f3374b")  # slow sine -> exponential cutoff
    patch.cable(9, 0, 5, 4, color="#ffb437")  # slow sine -> resonance
    patch.cable(10, 0, 5, 3)  # 14 Hz sine -> linear filter FM audition
    patch.cable(3, 8, 5, 5)  # track A gate -> filter Gate
    patch.cable(3, 4, 5, 6)  # track A CV2 -> filter Accent
    patch.cable(5, 1, 6, 1)
    patch.cable(6, 0, 7, 0)
    patch.cable(6, 0, 7, 1)
    patch.write("test-303-voice.vcv")


def generate_prog_sequencer_303_patch() -> None:
    """Replace the reverb-equipped smoke-303's Foundry with Prog Sequencer."""
    source_patch = _build_reverb_303_patch()
    patch = Patch(zoom=0.82, grid_offset=(-8.0, -0.15))
    patch.modules = [module for module in source_patch.modules if module["id"] != 3]
    patch.modules = [
        (
            module(
                2,
                "TriggerFish-Elements",
                "TfTransport",
                (7, 0),
                data={"state": 2},
            )
            if item["id"] == 2
            else item
        )
        for item in patch.modules
    ]
    for item in patch.modules:
        if item["id"] == 1:
            item["data"]["text"] = (
                "TriggerFish Prog Sequencer 303 smoke test\n\n"
                "This is the hand-tuned, reverb-equipped smoke-303 "
                "signal path with Foundry replaced by TfProgSequencer. The "
                "complete 16-bar "
                "Gibber Acid line is the short program visible in the "
                "sequencer: i for 8 bars, -iv for 4, then -V for 4.\n\n"
                "The test deliberately uses repetition, rests, slides, sparse "
                "accents, an octave mark, derived patterns, mode replacement, "
                "and concatenation as a syntax litmus test. Select an output "
                "device in Audio-8 and start with monitor volume low. The "
                "shared Transport can restart, pause, play, or stop the "
                "whole patch without changing the program."
            )

    program = """acid = sequence {
  subdiv 16n
  tonic D#@2
  scale minor
  notes ^1 1!3 ^5 7 ^1 ~ 1 ~ ^8 >1, ~ 1 >^1, >1
  velocity .5
  gate .5
  glide .5
}

iv = acid |> modulate_degree 4 |> octave -1
v  = acid |> modulate_degree 5 |> octave -1 |> scale major
song = acid * 8 + iv * 4 + v * 4
play song
"""
    patch.add(
        module(
            3,
            "TriggerFish-Elements",
            "TfProgSequencer",
            (23, 0),
            data={"source": program, "languageVersion": 1},
        )
    )

    # Musical settings captured from the Rack-saved Prog-303 audition patch.
    # Keep these overrides local to this patch: the standalone 303 smoke patch
    # deliberately retains its separate reference settings.
    prog_303_parameter_overrides = {
        4: {
            0: 1.0,  # oscillator octave
            4: 0.17951805889606476,  # saw / square morph
        },
        5: {
            0: 0.6371597647666931,  # cutoff pitch (about 407 Hz)
            1: 0.48494037985801697,  # resonance
            2: 1.0,  # high resonance range
            3: 8.949535369873047,  # drive in dB
            4: 0.14879514276981354,  # bass extension
            9: -0.06006646901369095,  # normal envelope decay (about 871 ms)
        },
    }
    for item in patch.modules:
        overrides = prog_303_parameter_overrides.get(item["id"], {})
        for parameter in item.get("params", []):
            if parameter["id"] in overrides:
                parameter["value"] = overrides[parameter["id"]]

    patch.cables = [
        cable
        for cable in source_patch.cables
        if cable["outputModuleId"] != 3 and cable["inputModuleId"] != 3
    ]
    for index, cable in enumerate(patch.cables):
        cable["id"] = 1000 + index
    patch.cable(2, 0, 3, 0)  # fixed 24 PPQN master clock
    patch.cable(2, 2, 3, 1)  # Transport reset -> Prog Sequencer reset
    patch.cable(2, 1, 3, 2)  # Transport RUN gate -> pause/resume
    patch.cable(3, 0, 4, 0)  # pitch, including programmed slides -> oscillator
    patch.cable(3, 1, 5, 5)  # gate -> 303 articulation
    patch.cable(3, 4, 5, 6)  # dedicated sparse accent -> 303 accent
    patch.write("test-prog-sequencer-303.vcv")


def _build_reverb_303_patch() -> Patch:
    """Build the reverb-equipped 303 path used by the Prog-303 smoke patch."""
    source = json.loads((ROOT / "test-303-voice.vcv").read_text(encoding="utf-8"))
    patch = Patch(zoom=0.82, grid_offset=(-8.0, -0.15))
    # Room Reverb already has equal-power Mix and an output Level control, so
    # replace smoke-303's output VCMixer instead of adding stereo master VCAs.
    # Also remove the three free-running filter/FM LFOs: this patch is a
    # reverb regression, so motion heard after transport stops must come from
    # the room engine rather than from an indefinitely sustained source.
    removed_module_ids = {6, 8, 9, 10}
    patch.modules = [
        module for module in source["modules"] if module["id"] not in removed_module_ids
    ]
    patch.cables = [
        cable
        for cable in source["cables"]
        if cable["outputModuleId"] not in removed_module_ids
        and cable["inputModuleId"] not in removed_module_ids
    ]

    note = next(module for module in patch.modules if module["model"] == "Notes")
    note["data"]["text"] = (
        "TriggerFish Room Reverb smoke test\n\n"
        "SETUP\nSelect an output device in Audio-8 and start with monitor "
        "volume low. The programmed 303 line feeds Room Reverb as a mono "
        "source; its independent LEFT and RIGHT outputs feed Audio-8 inputs "
        "1 and 2 directly. Use the reverb's LEVEL as the output trim.\n\n"
        "This reverb test deliberately removes the smoke-303 filter/FM LFOs "
        "so tail movement can be assessed without external modulation.\n\n"
        "Start by turning MIX up and sweep ER / TAIL from the geometric "
        "reflections through the inferred centre blend to the two-stage "
        "velvet late field. Drag the amber "
        "source and blue listener dots on the room plan, then try SIZE, "
        "ASPECT, PRE DELAY, DECAY, DAMPING, DIFFUSE, MOD, SHIMMER, "
        "WIDTH, LOW CUT, "
        "and HIGH CUT. The AUDIO input accepts Scene Pack 4's polyphonic "
        "source cable; each connected channel gets its own draggable marker."
    )
    note["pos"] = [-9, 0]

    audio_module = next(module for module in patch.modules if module["id"] == 7)
    audio_module["pos"] = [85, 0]
    reverb_values = [*DEFAULT_PARAMS["TfReverb"]]
    reverb_values[16] = -6.0
    patch.add(
        module(
            11,
            "TriggerFish-Elements",
            "TfReverb",
            (65, 0),
            values=reverb_values,
        )
    )
    # Reassign cable IDs after filtering so the checked-in patch is canonical.
    for cable_id, cable in enumerate(patch.cables, start=1000):
        cable["id"] = cable_id
    patch.cable(5, 1, 11, 0)  # 303 mono VCA output -> reverb
    patch.cable(11, 0, 7, 0)  # stereo left -> Audio-8 input 1
    patch.cable(11, 1, 7, 1)  # stereo right -> Audio-8 input 2
    return patch


def generate_two_source_reverb_patch() -> None:
    patch = Patch(zoom=0.66, grid_offset=(-7.0, -0.1))
    patch.add(
        notes(
            1,
            "TriggerFish musical two-source Room Reverb test\n\n"
            "At 96 BPM, two independent Prog Sequencers play an interlocking "
            "D-Dorian texture. Source 1 is a slow two-saw Unison "
            "Oscillator bass with a centred sub oscillator, analogue drift, "
            "and long 4072 filter and amplifier envelopes. Source 2 is a "
            "louder, driven descending folded arpeggio with a shorter, "
            "brighter 4072 voice. Its cutoff follows a slow 18-beat bipolar "
            "CV triangle across the four-beat note phrase, while a per-note "
            "AD envelope animates the wavefolder depth.\n\n"
            "The shared Transport drives Clock, Run, and Reset into both "
            "sequencers. RESTART plays from beat zero, PAUSE preserves the "
            "current position, PLAY continues it, and STOP returns to zero. "
            "Both voices therefore remain synchronized through every action.\n\n"
            "Scene Pack 4 keeps the voices as two independent room sources. "
            "Room Reverb uses the Superlush scene with its wet high cut "
            "lowered to 3 kHz. "
            "Drag their numbered amber markers to audition position, early "
            "reflections, and the late-field handoff; disconnect either Scene "
            "Pack input to confirm that only its marker disappears. The scope "
            "shows the packed sources and stereo result. Select an output "
            "device in Audio-8, then start with monitor volume low.",
            pos=(-11, 0),
        )
    )

    patch.add(
        module(
            2,
            "TriggerFish-Elements",
            "TfTransport",
            (5, 0),
            values=[96.0, 0.0, 0.0, 0.0, 0.0],
            data={"state": 2},
        )
    )

    bass_program = """bass = sequence {
  subdiv 2n
  tonic D@1
  scale dorian
  notes 1_2 7_2 4_2 5 7
  gate .96
}

play bass
"""
    arpeggio_program = """arpeggio = sequence {
  subdiv 16n
  tonic D@4
  scale dorian
  notes 3' 1' 5 3 ; 1' 5 3 1 ; 2' 7 5 2 ; 7 5 2 1
  gate .52
  cv1 4 -4 |> interp linear |> rate 1/9
  cv2 env ad 5ms 16n depth 2 curve 0
}

play arpeggio
"""
    patch.add(
        module(
            3,
            "TriggerFish-Elements",
            "TfProgSequencer",
            (25, 0),
            data={"source": bass_program, "languageVersion": 1},
        )
    )
    patch.add(
        module(
            4,
            "TriggerFish-Elements",
            "TfProgSequencer",
            (55, 0),
            data={
                "source": arpeggio_program,
                "languageVersion": 1,
                # tfui::HeatmapPalette::CrtGreen. The bass remains Magma so
                # both live programs are immediately distinguishable.
                "editorHeatmap": 5,
            },
        )
    )

    unison_values = [*DEFAULT_PARAMS["TfUnisonOscillator"]]
    unison_values[2] = 2.0  # two saw oscillators
    unison_values[7] = 0.55  # a wider, but still tonal, unison spread
    unison_values[8] = 0.45
    unison_values[9] = 0.2965061068534851  # centred sub oscillator
    unison_values[11] = 0.62
    unison_values[12] = 0.04
    unison_values[13] = 0.29333335161209106
    unison_values[14] = 0.3920000493526459
    patch.add(
        module(
            5,
            "TriggerFish-Elements",
            "TfUnisonOscillator",
            (29, 1),
            values=unison_values,
        )
    )

    bass_voice_values = [*DEFAULT_PARAMS["Tf4072VoiceCore"]]
    bass_voice_values[0] = 3.0938026905059814
    bass_voice_values[1] = 0.4920482337474823
    bass_voice_values[2] = 10.387954711914062
    bass_voice_values[3] = 0.55
    bass_voice_values[7] = 0.9
    bass_voice_values[9] = math.log10(0.025)
    bass_voice_values[10] = math.log10(1.8)
    bass_voice_values[11] = 0.12
    bass_voice_values[12] = math.log10(0.65)
    bass_voice_values[13] = 0.0
    bass_voice_values[14] = 0.25
    bass_voice_values[15] = math.log10(0.015)
    bass_voice_values[16] = math.log10(0.5)
    bass_voice_values[17] = 0.85
    bass_voice_values[18] = math.log10(0.6)
    bass_voice_values[19] = 0.0
    patch.add(
        module(
            6,
            "TriggerFish-Elements",
            "Tf4072VoiceCore",
            (41, 1),
            values=bass_voice_values,
        )
    )

    wavefold_values = [*DEFAULT_PARAMS["TfWavefoldOscillator"]]
    wavefold_values[2] = 0.28
    wavefold_values[3] = 0.5012047290802002
    wavefold_values[4] = 0.08
    wavefold_values[7] = 0.75  # CV2 AD adds a 0.3 fold-depth sweep at 2 V
    wavefold_values[10] = 1.0  # Hinge character
    wavefold_values[11] = 0.58
    wavefold_values[12] = 0.40800002217292786
    wavefold_values[13] = 0.49066662788391113
    wavefold_values[14] = 0.5093333721160889
    patch.add(
        module(
            7,
            "TriggerFish-Elements",
            "TfWavefoldOscillator",
            (65, 1),
            values=wavefold_values,
        )
    )

    arpeggio_voice_values = [*DEFAULT_PARAMS["Tf4072VoiceCore"]]
    arpeggio_voice_values[0] = 0.6673828959465027
    arpeggio_voice_values[1] = 0.0
    arpeggio_voice_values[2] = 3.5975842475891113
    arpeggio_voice_values[3] = 0.5493977069854736
    arpeggio_voice_values[4] = 0.8789158463478088
    arpeggio_voice_values[7] = 1.0
    arpeggio_voice_values[8] = 0.7283129692077637
    arpeggio_voice_values[9] = math.log10(0.004)
    arpeggio_voice_values[10] = 0.04110236465930939
    arpeggio_voice_values[11] = 0.08
    arpeggio_voice_values[12] = math.log10(0.16)
    arpeggio_voice_values[13] = 0.0
    arpeggio_voice_values[14] = 0.35
    arpeggio_voice_values[15] = math.log10(0.002)
    arpeggio_voice_values[16] = 0.043807897716760635
    arpeggio_voice_values[17] = 0.12
    arpeggio_voice_values[18] = math.log10(0.15)
    arpeggio_voice_values[19] = 0.0
    patch.add(
        module(
            8,
            "TriggerFish-Elements",
            "Tf4072VoiceCore",
            (77, 1),
            values=arpeggio_voice_values,
        )
    )

    patch.add(module(9, "TriggerFish-Elements", "TfScenePack4", (85, 0)))

    reverb_values = [*DEFAULT_PARAMS["TfReverb"]]
    # Superlush, with a deliberately darker 3 kHz wet high cut. LEVEL remains
    # a separate output trim, matching the preset menu's behaviour.
    reverb_values[0] = 0.78
    reverb_values[1] = 0.55
    reverb_values[2] = 0.112
    reverb_values[3] = 0.50
    reverb_values[4] = 0.32
    reverb_values[5] = 0.50
    reverb_values[6] = 0.68
    reverb_values[7] = 0.84
    reverb_values[8] = 0.26
    reverb_values[9] = 1.0
    reverb_values[10] = 1.0
    reverb_values[11] = 0.75
    reverb_values[12] = 0.75
    reverb_values[13] = 0.4580135308
    reverb_values[14] = math.log(3.0) / math.log(20.0)
    reverb_values[15] = 0.40
    reverb_values[16] = -6.0
    reverb_values[17] = 0.45
    for source in range(1, 8):
        reverb_values[18 + 2 * (source - 1)] = 0.50
        reverb_values[19 + 2 * (source - 1)] = 0.32
    patch.add(
        module(
            10,
            "TriggerFish-Elements",
            "TfReverb",
            (91, 0),
            values=reverb_values,
        )
    )
    patch.add(audio(11, (107, 0)))
    patch.add(scope(12, (117, 0)))

    for sequencer_id in (3, 4):
        patch.cable(2, 0, sequencer_id, 0)  # fixed 24 PPQN master clock
        patch.cable(2, 2, sequencer_id, 1)  # shared reset
        patch.cable(2, 1, sequencer_id, 2)  # shared run/pause gate

    patch.cable(3, 0, 5, 0)  # bass pitch -> Unison Oscillator
    patch.cable(3, 0, 6, 2)  # bass pitch -> 4072 cutoff tracking
    patch.cable(3, 1, 6, 6)  # bass gate -> envelopes
    patch.cable(3, 2, 6, 7)  # bass trigger -> envelope retrigger
    patch.cable(5, 0, 6, 0)  # unison saw/sub mix -> 4072 filter
    patch.cable(6, 1, 9, 0)  # enveloped bass -> source 1

    patch.cable(4, 0, 7, 0)  # arpeggio pitch -> folding oscillator
    patch.cable(4, 0, 8, 2)  # arpeggio pitch -> 4072 cutoff tracking
    patch.cable(4, 1, 8, 6)  # arpeggio gate -> envelopes
    patch.cable(4, 2, 8, 7)  # arpeggio trigger -> envelope retrigger
    patch.cable(4, 5, 8, 3)  # slow CV1 triangle -> filter modulation
    patch.cable(4, 6, 7, 3)  # per-note CV2 AD -> wavefolder Fold modulation
    patch.cable(7, 1, 8, 0)  # folded oscillator -> 4072 filter
    patch.cable(8, 1, 9, 1)  # enveloped arpeggio -> source 2

    patch.cable(9, 0, 10, 0)
    patch.cable(10, 0, 11, 0)
    patch.cable(10, 1, 11, 1)
    patch.cable(9, 0, 12, 0)
    patch.cable(10, 0, 12, 1)
    patch.cable(10, 1, 12, 2)
    patch.write("test-room-reverb-two-sources.vcv")


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
    patch.cable(2, 1, 4, 6)  # MIDI gate -> both internal envelopes
    patch.cable(3, 2, 4, 0)  # saw -> filter input
    patch.cable(4, 1, 5, 1)  # VCA output -> master
    patch.cable(5, 0, 6, 0)
    patch.cable(5, 0, 6, 1)
    patch.write("test-4072-voice.vcv")


def generate_electric_piano_patch() -> None:
    patch = Patch(zoom=0.82, grid_offset=(-1, -0.1))
    patch.add(
        notes(
            1,
            "TriggerFish Electric Piano playable prototype\n\n"
            "Select MIDI and audio devices, then play from a velocity-sensitive "
            "keyboard. MIDI-CV is set to 16-channel polyphonic mode; pitch, gate, "
            "and velocity feed the instrument directly.\n\n"
            "Start with VEL CURVE and DYNAMICS to match your controller. BODY and "
            "BELL balance the resonances, HAMMER changes physical tip stiffness, "
            "and COUPLING moves from a short isolated tine to a sustained "
            "common-base tuning fork. "
            "TONE changes vertical tine/pickup alignment; "
            "GAP changes pickup proximity and dynamic response. DECAY scales the "
            "natural tine/tone-bar ring time; RELEASE sets damper speed. MECHANICS "
            "adds synchronized key noise, and DRIVE overloads "
            "the shared amplifier. DIRECT POLY bypasses that final amp.\n\n"
            "The -6 dB master protects the first listen. PEDAL accepts a gate; "
            "patch your preferred MIDI CC-to-CV module there for sustain.",
        )
    )
    patch.add(midi(2, (17, 0), channels=16))
    patch.add(module(3, "TriggerFish-Elements", "TfElectricPiano", (29, 0)))
    patch.add(mixer(4, (45, 0), (0.5011872336, 0.7, 0.0, 0.0, 0.0)))
    patch.add(mixer(5, (54, 0), (0.5011872336, 0.7, 0.0, 0.0, 0.0)))
    patch.add(audio(6, (63, 0)))

    patch.cable(2, 0, 3, 0)  # MIDI pitch -> 1V/oct
    patch.cable(2, 1, 3, 1)  # MIDI gate -> key gate
    patch.cable(2, 2, 3, 2)  # MIDI velocity -> strike velocity
    patch.cable(3, 1, 4, 1)  # shared amplifier left -> left master
    patch.cable(3, 2, 5, 1)  # shared amplifier right -> right master
    patch.cable(4, 0, 6, 0)
    patch.cable(5, 0, 6, 1)
    patch.write("test-electric-piano.vcv")


def generate_wavefold_patch() -> None:
    patch = Patch(zoom=0.78, grid_offset=(-1, -0.1))
    patch.add(
        notes(
            1,
            "TriggerFish Wavefold Oscillator smoke test\n\n"
            "Select MIDI and audio devices, then play from a keyboard. MIDI "
            "pitch drives the internal sine/triangle oscillator. Master "
            "channel 1 receives FOLD OUT directly for uncoloured evaluation. "
            "Channel 2 contains an ADSR-controlled TriggerFish VCA and starts "
            "muted; raise it and mute channel 1 for gated playing.\n\n"
            "The module starts with the Lockhart folder selected. Compare "
            "LOCKHART, HINGE, and SERGE, then try SINE-TRI, FOLD, and SYMMETRY. "
            "SPEED sets the timescale of three independent drift processes; "
            "the WAVE, FOLD, and SYM sliders set their individual depths. "
            "The scope compares FOLD OUT with OSC OUT and triggers from OSC "
            "OUT, so the traces remain stable while the waveform changes. "
            "The oscillator uses 4x oversampling by default; a 2x lower-CPU "
            "mode is available from its context menu. "
            "For an intentionally severe alias test, patch OSC OUT back to "
            "IN, set FOLD fully clockwise, and play around C8; the external "
            "path leaves fold depth untapered. "
            "The OSC OUT jack provides the waveform before folding. Patching "
            "another nominal +/-5 V source into IN replaces the internal "
            "oscillator at the folder input.",
        )
    )
    patch.add(midi(2, (16, 0)))
    patch.add(module(3, "Fundamental", "ADSR", (28, 0)))
    patch.add(
        module(
            4,
            "TriggerFish-Elements",
            "TfWavefoldOscillator",
            (38, 0),
        )
    )
    patch.add(module(5, "TriggerFish-Elements", "TfVCA", (50, 0)))
    patch.add(mixer(6, (59, 0), (0.5011872336, 0.7, 0.0, 0.0, 0.0)))
    patch.add(audio(7, (68, 0)))
    patch.add(scope(8, (76, 0)))

    patch.cable(2, 0, 4, 0)  # MIDI pitch -> oscillator 1V/oct
    patch.cable(2, 1, 3, 4)  # MIDI gate -> ADSR gate
    patch.cable(4, 1, 5, 0)  # folded audio -> VCA
    patch.cable(3, 0, 5, 1)  # ADSR -> VCA linear CV
    patch.cable(4, 1, 6, 1)  # direct folded reference -> master channel 1
    patch.cable(5, 0, 6, 2)  # enveloped voice -> muted master channel 2
    patch.cable(6, 0, 7, 0)
    patch.cable(6, 0, 7, 1)
    patch.cable(4, 1, 8, 0)  # folded output -> scope channel 1
    patch.cable(4, 0, 8, 1)  # oscillator output -> scope channel 2
    patch.cable(4, 0, 8, 2)  # oscillator output -> scope trigger
    patch.write("test-wavefold-oscillator.vcv")


def generate_unison_patch() -> None:
    patch = Patch(zoom=0.76, grid_offset=(-1, -0.1))
    patch.add(
        notes(
            1,
            "TriggerFish Unison Oscillator smoke test\n\n"
            "Select MIDI and audio devices, then play from a keyboard. MIDI "
            "pitch drives a seven-voice stereo saw stack; the ADSR controls "
            "independent left and right TriggerFish VCAs. Two -6 dB masters "
            "feed Audio-8 outputs 1 and 2.\n\n"
            "Try VOICES, SPREAD, and WIDTH first. Switch WAVE to PULSE; PWM at "
            "centre leaves PULSE WIDTH fixed, while moving it applies the "
            "internal LFO at the selected RATE. A PW CV cable replaces that LFO. "
            "Raise SUB LEVEL and "
            "compare CENTRE with STACK; SUB is also available separately. "
            "The SLOP controls add shared hum, common drift, independent drift, "
            "and fixed per-oscillator tracking differences. The scope shows "
            "the stereo stack with stable triggering.",
        )
    )
    patch.add(midi(2, (16, 0)))
    patch.add(module(3, "Fundamental", "ADSR", (28, 0)))
    patch.add(module(4, "TriggerFish-Elements", "TfUnisonOscillator", (38, 0)))
    patch.add(module(5, "TriggerFish-Elements", "TfVCA", (60, 0)))
    patch.add(module(6, "TriggerFish-Elements", "TfVCA", (69, 0)))
    patch.add(mixer(7, (78, 0), (0.5011872336, 0.7, 0.0, 0.0, 0.0)))
    patch.add(mixer(8, (87, 0), (0.5011872336, 0.7, 0.0, 0.0, 0.0)))
    patch.add(audio(9, (96, 0)))
    patch.add(scope(10, (104, 0)))

    patch.cable(2, 0, 4, 0)  # MIDI pitch -> oscillator 1V/oct
    patch.cable(2, 1, 3, 4)  # MIDI gate -> ADSR gate
    patch.cable(4, 1, 5, 0)  # left -> left VCA
    patch.cable(4, 2, 6, 0)  # right -> right VCA
    patch.cable(3, 0, 5, 1)
    patch.cable(3, 0, 6, 1)
    patch.cable(5, 0, 7, 1)
    patch.cable(6, 0, 8, 1)
    patch.cable(7, 0, 9, 0)
    patch.cable(8, 0, 9, 1)
    patch.cable(4, 1, 10, 0)  # left -> scope channel 1
    patch.cable(4, 2, 10, 1)  # right -> scope channel 2
    patch.cable(4, 1, 10, 2)  # left -> scope trigger
    patch.write("test-unison-oscillator.vcv")


def generate_scene_pack4_patch() -> None:
    patch = Patch(zoom=0.76, grid_offset=(-1, -0.1))
    patch.add(
        notes(
            1,
            "TriggerFish Scene Pack 4 smoke test\n\n"
            "Four Fundamental oscillators feed the four source inputs. Scene "
            "Pack compacts every connected mono or polyphonic channel into "
            "one polyphonic output, which feeds Room Reverb.\n\n"
            "The reverb room plan shows four numbered amber sources. Drag "
            "them to place each oscillator independently; disconnect a Scene "
            "Pack input and its marker disappears. Feed another Scene Pack's "
            "poly output into one input to exercise the eight-source path. "
            "Start with monitor volume low.",
        )
    )
    pitches = (-2.0, -1.5, -1.0, -0.5)
    for lane, pitch in enumerate(pitches):
        values = DEFAULT_PARAMS["VCO"].copy()
        values[2] = pitch
        patch.add(
            module(2 + lane, "Fundamental", "VCO", (16 + 11 * lane, 0), values=values)
        )
    patch.add(module(6, "TriggerFish-Elements", "TfScenePack4", (60, 0)))
    reverb_values = [*DEFAULT_PARAMS["TfReverb"]]
    reverb_values[16] = -6.0
    patch.add(
        module(
            7,
            "TriggerFish-Elements",
            "TfReverb",
            (68, 0),
            values=reverb_values,
        )
    )
    patch.add(audio(8, (89, 0)))
    patch.add(scope(9, (97, 0)))

    for lane in range(4):
        patch.cable(2 + lane, 0, 6, lane)
    patch.cable(6, 0, 7, 0)
    patch.cable(7, 0, 8, 0)
    patch.cable(7, 1, 8, 1)
    patch.cable(6, 0, 9, 0)
    patch.cable(7, 0, 9, 1)
    patch.cable(7, 1, 9, 2)
    patch.write("test-scene-pack4.vcv")


if __name__ == "__main__":
    generate_slop4_patch()
    generate_vdpo_patch()
    generate_303_voice_patch()
    generate_two_source_reverb_patch()
    generate_prog_sequencer_303_patch()
    generate_4072_voice_patch()
    generate_electric_piano_patch()
    generate_wavefold_patch()
    generate_unison_patch()
    generate_scene_pack4_patch()
