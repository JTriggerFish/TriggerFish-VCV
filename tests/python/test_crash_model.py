from pathlib import Path
import sys

import numpy as np

import _triggerfish_dsp as native
from triggerfish_percussion.audio_io import AudioBuffer
from triggerfish_percussion.crash_corpus import (
    isolate_capture,
    load_fit_cells,
    write_cells_manifest,
)
from triggerfish_percussion.crash_model import (
    CrashEvent,
    CrashFit,
    render_crash,
    render_crash_sequence,
)
from triggerfish_percussion.report import (
    ReportPair,
    _common_magnitude_db,
    write_comparison_report,
)

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))
import build_toontrack_cymbal_sweep as sweep  # noqa: E402


def test_native_crash_parameters_round_trip_and_render():
    fit = CrashFit(strength_gamma=2.0, output_gain=2.0)
    first = render_crash(fit, 0.1, strength=0.7, seed=12)
    repeated = render_crash(fit, 0.1, strength=0.7, seed=12)
    assert first.shape == (4800,)
    assert np.array_equal(first, repeated)
    assert np.isfinite(first).all()
    assert np.max(np.abs(first)) > 0


def test_native_crash_sequence_retriggers_one_stateful_body():
    events = (
        CrashEvent(0.0, 0.4, seed=1),
        CrashEvent(0.05, 0.8, location=0.0, seed=2),
    )
    sequence = render_crash_sequence(CrashFit(), 0.1, events)
    assert sequence.shape == (4800,)
    assert np.isfinite(sequence).all()
    assert np.max(np.abs(sequence[2400:])) > 0


def test_cymbal_sweep_uses_crash2_articulations():
    hits = sweep.build_sequence("sd3-crash2", (32, 96), 2, 1.0, 2.0, 0.5)
    assert len(hits) == 20
    assert {(hit.articulation, hit.midi_note) for hit in hits} == {
        ("edge", 49),
        ("bow-tip", 27),
        ("bow-shank", 92),
        ("bell-tip", 93),
        ("bell-shank", 28),
    }
    midi = sweep.build_midi(hits, 0.08)
    assert midi.startswith(b"MThd")
    assert b"edge/v032/r01" in midi


def test_report_contains_only_reference_and_current_model(tmp_path):
    time = np.arange(4096) / 48000.0
    reference = AudioBuffer(np.sin(2 * np.pi * 1000 * time) * np.exp(-8 * time), 48000)
    synthesis = AudioBuffer(np.sin(2 * np.pi * 1100 * time) * np.exp(-9 * time), 48000)
    assets = tmp_path / "report-assets"
    assets.mkdir()
    (assets / "stale-reference.wav").write_bytes(b"stale")
    path = write_comparison_report(
        (ReportPair("one", reference, synthesis),), tmp_path / "report.html", "test"
    )
    html = path.read_text(encoding="utf-8")
    assert ">Reference<" in html
    assert ">TriggerFish<" in html
    assert "old model" not in html.lower()
    assert not (assets / "stale-reference.wav").exists()
    assert (tmp_path / "report-assets/1-one-reference.wav").is_file()
    assert (tmp_path / "report-assets/1-one-triggerfish.wav").is_file()


def test_report_spectrograms_share_one_level_reference():
    louder, quieter = _common_magnitude_db(
        np.asarray([[1.0, 0.1]]), np.asarray([[0.1, 0.01]]), -100.0
    )
    assert louder[0, 0] == 0.0
    assert quieter[0, 0] == -20.0


def test_capture_isolation_preserves_level_and_labels_cells(tmp_path):
    sample_rate = 48000
    samples = np.zeros(2 * sample_rate)
    samples[round(0.5 * sample_rate)] = 0.25
    samples[round(1.5 * sample_rate)] = 0.75
    sweep_manifest = {
        "hits": [
            {
                "index": 0,
                "onset_seconds": 0.5,
                "articulation": "edge",
                "velocity": 32,
                "repeat": 1,
            },
            {
                "index": 1,
                "onset_seconds": 1.5,
                "articulation": "bell-shank",
                "velocity": 96,
                "repeat": 2,
            },
        ]
    }
    output = tmp_path / "cells"
    cells = isolate_capture(
        AudioBuffer(samples, sample_rate), sweep_manifest, output, 0.2
    )
    manifest = output / "cells.json"
    write_cells_manifest(manifest, cells, tmp_path / "capture.wav")
    loaded = load_fit_cells(manifest)
    assert len(loaded) == 2
    assert loaded[0].location == 1.0
    assert loaded[1].location == 0.0
    assert loaded[1].hardness > loaded[0].hardness
    assert np.max(np.abs(loaded[0].reference.samples)) == np.float32(0.25)
    peak = np.argmax(np.abs(loaded[0].reference.samples))
    assert abs(peak - round(0.05 * sample_rate)) < round(0.001 * sample_rate)
    assert cells[0].cell_onset_seconds == 0.05
