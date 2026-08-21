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
    "TfScenePack4": [
        3.5,
        3.5,
        4.5,
        4.5,
        3.5,
        4.5,
        5.5,
        3.5,
        4.5,
        6.5,
        3.5,
        4.5,
    ],
    "TfReverb": [
        0.5,
        0.5,
        0.28,
        0.0,
        0.5,
        0.35,
        0.5,
        0.682,
        0.55,
        0.18,
        0.75,
        0.4,
        0.6130368568946039,
        0.0,
        0.0,
        0.0,
        0.9039693650225663,
        0.35,
        0.0,
        0.0,
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

    # Auditioned settings saved with the musical smoke patch.
    voice_params = DEFAULT_PARAMS["Tf303VoiceCore"].copy()
    voice_params[0] = 0.128505349159241
    voice_params[1] = 0.408434182405472
    voice_params[2] = 1.0
    voice_params[4] = 0.479517936706543
    voice_params[5] = 0.406198114156723
    voice_params[6] = 0.2313252389431
    voice_params[7] = 0.60891592502594
    voice_params[8] = 1.0
    voice_params[9] = -1.18070936203003
    voice_params[11] = 0.589156568050385
    voice_params[12] = 0.497590363025665
    oscillator_params = DEFAULT_PARAMS["Tf303Oscillator"].copy()
    oscillator_params[2] = math.log10(0.100)
    oscillator_params[4] = 0.296385675668716
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
            "volume low. Impromptu Clocked and Foundry run the full 16-bar "
            "Gibber Acid playground line by fasttriggerfish and thecharlie at "
            "120 BPM: i for 8 bars, -iv for 4, then -V for 4.\n\n"
            "Foundry track A contains the notes, rests, and accents. Track B "
            "reproduces the repeating 1, 1, 100, 100 ms glide cycle; the first "
            "two quarters are effectively instantaneous. Independent slow "
            "bipolar sines modulate cutoff and resonance. A separate 14 Hz "
            "sine is connected "
            "to linear filter FM with its attenuverter at zero. VCA OUT feeds "
            "the -6 dB stereo master.\n\n"
            "The cutoff, resonance, CV depths, and slide time are set for this "
            "line; other TriggerFish parameters use their declared defaults. Try "
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
    source_patch = json.loads(
        (ROOT / "test-room-reverb.vcv").read_text(encoding="utf-8")
    )
    patch = Patch(zoom=0.82, grid_offset=(-8.0, -0.15))
    patch.modules = [module for module in source_patch["modules"] if module["id"] != 3]
    for item in patch.modules:
        if item["id"] == 1:
            item["data"]["text"] = (
                "TriggerFish Prog Sequencer 303 smoke test\n\n"
                "This is the tuned, reverb-equipped smoke-303 signal path with "
                "Foundry replaced by TfProgSequencer. The complete 16-bar "
                "Gibber Acid line is the short program visible in the "
                "sequencer: i for 8 bars, -iv for 4, then -V for 4.\n\n"
                "The test deliberately uses repetition, rests, slides, sparse "
                "accents, an octave mark, derived patterns, mode replacement, "
                "and concatenation as a syntax litmus test. Select an output "
                "device in Audio-8 and start with monitor volume low."
            )

    program = """acid = sequence {
  cycle 16
  tonic D#@2
  scale minor
  notes 1!4 5 7 1!2 8 1, 1 1, 1
  articulation x!7 ~ > ~ >!2 ~ >!3
  velocity .5
  accent + .!3 + . + . + .!2 + .
  gate .5
  glide .8
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
            data={"source": program, "languageVersion": 4},
        )
    )

    patch.cables = [
        cable
        for cable in source_patch["cables"]
        if cable["outputModuleId"] != 3 and cable["inputModuleId"] != 3
    ]
    for index, cable in enumerate(patch.cables):
        cable["id"] = 1000 + index
    patch.cable(2, 1, 3, 0)  # Clocked x4 -> Prog Sequencer clock
    patch.cable(2, 4, 3, 1)  # Clocked reset -> Prog Sequencer reset
    patch.cable(3, 0, 4, 0)  # pitch, including programmed slides -> oscillator
    patch.cable(3, 1, 5, 5)  # gate -> 303 articulation
    patch.cable(3, 4, 5, 6)  # dedicated sparse accent -> 303 accent
    patch.write("test-prog-sequencer-303.vcv")


def generate_room_reverb_patch() -> None:
    """Insert Room Reverb into the full smoke-303 patch."""
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
        "Start by turning MIX up and use the EARLY and TAIL trims to isolate "
        "the geometric reflections and two-stage velvet late field. Drag the amber "
        "source and blue listener dots on the room plan, then try SIZE, "
        "ASPECT, PRE DELAY, DECAY, DAMPING, DIFFUSE, MOD, SHIMMER, "
        "WIDTH, LOW CUT, "
        "and HIGH CUT. The AUDIO/X/Y/Z inputs also "
        "accept the aligned polyphonic scene cables from Scene Pack 4."
    )
    note["pos"] = [-9, 0]

    audio_module = next(module for module in patch.modules if module["id"] == 7)
    audio_module["pos"] = [85, 0]
    reverb_values = [*DEFAULT_PARAMS["TfReverb"]]
    reverb_values[18] = -6.0
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
    patch.write("test-room-reverb.vcv")


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
            "Four Fundamental oscillators feed the four local mono inputs. "
            "Scene Pack compacts them into aligned AUDIO/X/Y/Z polyphonic "
            "outputs. AUDIO feeds both -6 dB masters; the first packed source "
            "is therefore audible on outputs 1 and 2.\n\n"
            "Move each source's X, Y, and Z controls and inspect the packed "
            "AUDIO and X cables on the scope. Chain a second Scene Pack into "
            "the four BUS inputs to exercise the eight-source path. Start "
            "with monitor volume low.",
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
    patch.add(scope(7, (75, 0)))
    patch.add(mixer(8, (87, 0), (0.5011872336, 0.7, 0.0, 0.0, 0.0)))
    patch.add(mixer(9, (96, 0), (0.5011872336, 0.7, 0.0, 0.0, 0.0)))
    patch.add(audio(10, (105, 0)))

    for lane in range(4):
        patch.cable(2 + lane, 0, 6, 4 + lane)
    patch.cable(6, 0, 8, 1)
    patch.cable(6, 0, 9, 1)
    patch.cable(8, 0, 10, 0)
    patch.cable(9, 0, 10, 1)
    patch.cable(6, 0, 7, 0)
    patch.cable(6, 1, 7, 1)
    patch.cable(2, 0, 7, 2)
    patch.write("test-scene-pack4.vcv")


if __name__ == "__main__":
    generate_slop4_patch()
    generate_vdpo_patch()
    generate_303_voice_patch()
    generate_room_reverb_patch()
    generate_prog_sequencer_303_patch()
    generate_4072_voice_patch()
    generate_wavefold_patch()
    generate_unison_patch()
    generate_scene_pack4_patch()
