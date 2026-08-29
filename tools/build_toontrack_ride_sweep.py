"""Build deterministic MIDI and onset metadata for Toontrack ride captures.

The generated MIDI file can be imported into Superior Drummer 3 or EZdrummer 3.
Render it once for each reference cymbal, with humanization and processing disabled.
The adjacent JSON file records the exact onset, velocity, articulation, and MIDI note
of every hit so that reference-analysis tools do not have to rediscover onsets.
"""

from __future__ import annotations

import argparse
import json
import struct
from dataclasses import asdict, dataclass
from pathlib import Path

TICKS_PER_BEAT = 480
TEMPO_US_PER_BEAT = 500_000  # 120 BPM; 960 ticks per second.
TICKS_PER_SECOND = TICKS_PER_BEAT * 1_000_000 / TEMPO_US_PER_BEAT

ARTICULATIONS = {
    "sd3": (
        ("bow-tip", 51),
        ("bow-shank", 29),
        ("bell-tip", 30),
        ("bell-shank", 53),
        ("edge", 59),
    ),
    "ezd3": (
        ("bow", 51),
        ("bell", 53),
        ("edge", 59),
    ),
}


@dataclass(frozen=True)
class Hit:
    index: int
    onset_seconds: float
    articulation: str
    midi_note: int
    velocity: int
    repeat: int


def variable_length(value: int) -> bytes:
    if value < 0:
        raise ValueError("MIDI delta time cannot be negative")
    encoded = bytearray((value & 0x7F,))
    value >>= 7
    while value:
        encoded.insert(0, 0x80 | (value & 0x7F))
        value >>= 7
    return bytes(encoded)


def meta_event(kind: int, payload: bytes) -> bytes:
    return bytes((0xFF, kind)) + variable_length(len(payload)) + payload


def parse_velocities(value: str) -> tuple[int, ...]:
    try:
        velocities = tuple(int(part.strip()) for part in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "velocities must be comma-separated integers"
        ) from error
    if not velocities or any(velocity < 1 or velocity > 127 for velocity in velocities):
        raise argparse.ArgumentTypeError("velocities must be in the range 1..127")
    return velocities


def build_sequence(
    target: str,
    velocities: tuple[int, ...],
    repeats: int,
    gap_seconds: float,
    group_gap_seconds: float,
    lead_seconds: float,
) -> list[Hit]:
    hits: list[Hit] = []
    onset = lead_seconds
    for articulation_index, (articulation, note) in enumerate(ARTICULATIONS[target]):
        if articulation_index:
            onset += group_gap_seconds
        for velocity in velocities:
            for repeat in range(1, repeats + 1):
                hits.append(
                    Hit(
                        index=len(hits),
                        onset_seconds=onset,
                        articulation=articulation,
                        midi_note=note,
                        velocity=velocity,
                        repeat=repeat,
                    )
                )
                onset += gap_seconds
    return hits


def build_midi(hits: list[Hit], note_seconds: float) -> bytes:
    # Each tuple is (absolute tick, ordering at equal tick, event bytes).
    events: list[tuple[int, int, bytes]] = [
        (0, 0, meta_event(0x03, b"Toontrack ride calibration sweep")),
        (0, 1, meta_event(0x51, TEMPO_US_PER_BEAT.to_bytes(3, "big"))),
        (0, 2, meta_event(0x58, bytes((4, 2, 24, 8)))),
    ]
    note_ticks = max(1, round(note_seconds * TICKS_PER_SECOND))
    for hit in hits:
        tick = round(hit.onset_seconds * TICKS_PER_SECOND)
        label = (f"{hit.articulation}/v{hit.velocity:03d}/r{hit.repeat:02d}").encode(
            "ascii"
        )
        events.append((tick, 3, meta_event(0x06, label)))
        events.append((tick, 4, bytes((0x99, hit.midi_note, hit.velocity))))
        events.append((tick + note_ticks, 5, bytes((0x89, hit.midi_note, 0))))

    events.sort(key=lambda event: (event[0], event[1]))
    track = bytearray()
    previous_tick = 0
    for tick, _, event in events:
        track.extend(variable_length(tick - previous_tick))
        track.extend(event)
        previous_tick = tick
    track.extend(variable_length(round(2.0 * TICKS_PER_SECOND)))
    track.extend(meta_event(0x2F, b""))

    header = b"MThd" + struct.pack(">IHHH", 6, 0, 1, TICKS_PER_BEAT)
    return header + b"MTrk" + struct.pack(">I", len(track)) + track


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", choices=tuple(ARTICULATIONS), default="sd3")
    parser.add_argument(
        "--velocities",
        type=parse_velocities,
        default=parse_velocities("16,32,48,64,80,96,112,127"),
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--gap-seconds", type=float, default=4.0)
    parser.add_argument("--group-gap-seconds", type=float, default=4.0)
    parser.add_argument("--lead-seconds", type=float, default=1.0)
    parser.add_argument("--note-seconds", type=float, default=0.08)
    arguments = parser.parse_args()
    if arguments.repeats < 1:
        parser.error("--repeats must be at least one")
    if (
        min(
            arguments.gap_seconds,
            arguments.group_gap_seconds,
            arguments.lead_seconds,
            arguments.note_seconds,
        )
        <= 0
    ):
        parser.error("timing arguments must be positive")

    hits = build_sequence(
        arguments.target,
        arguments.velocities,
        arguments.repeats,
        arguments.gap_seconds,
        arguments.group_gap_seconds,
        arguments.lead_seconds,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_bytes(build_midi(hits, arguments.note_seconds))
    manifest_path = arguments.output.with_suffix(".json")
    manifest = {
        "schema": 1,
        "target": arguments.target,
        "midi_file": arguments.output.name,
        "midi_channel": 10,
        "ticks_per_beat": TICKS_PER_BEAT,
        "tempo_bpm": 120,
        "capture_guidance": {
            "usage": "private listening/validation unless vendor grants model-development permission",
            "disable_humanization": True,
            "disable_internal_processing": True,
            "recommended_outputs": ["direct", "overhead", "room"],
        },
        "hits": [asdict(hit) for hit in hits],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    duration = hits[-1].onset_seconds + 2.0 if hits else 0.0
    print(f"wrote {arguments.output} ({len(hits)} hits, {duration:.1f} seconds)")
    print(f"wrote {manifest_path}")


if __name__ == "__main__":
    main()
