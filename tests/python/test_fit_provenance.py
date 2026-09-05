"""An old WAV or fit must never be advertised as current-renderer output."""

import json
from copy import deepcopy
from types import SimpleNamespace
import numpy as np
import pytest
from triggerfish_percussion.audio_io import AudioBuffer, write_wav
from triggerfish_percussion.fit_provenance import verify_candidate
from triggerfish_percussion.fit_publication import publish_html


def fixture(tmp_path):
    audio = np.array([0.25, -0.5])
    metadata = dict(
        recipeKey="fixture",
        rendererSha256="current",
        descriptors=[dict(key="gain")],
        event=dict(seed=1),
        reference=dict(sha256="ref", sampleRate=2, referenceGainDb=0),
    )
    renderer = SimpleNamespace(
        metadata=metadata,
        initial={"gain": 1},
        sample_rate=2,
        reference=audio,
        render=lambda p, seconds: audio,
        decode=lambda value: value,
        request=lambda **kw: dict(pcm=audio * kw["fit"]["gain"]),
    )
    saved = dict(
        metadata=deepcopy(metadata), parameters={"gain": 1}, duration_seconds=1
    )
    (tmp_path / "search.json").write_text(json.dumps(saved))
    (tmp_path / "candidate.fit.json").write_text('{"gain":1}')
    for name in ("reference", "candidate"):
        write_wav(tmp_path / f"{name}.wav", AudioBuffer(audio, 2))
    return renderer, saved


def test_verified_report_fixture(tmp_path):
    renderer, saved = fixture(tmp_path)
    assert verify_candidate(renderer, tmp_path) == saved


@pytest.mark.parametrize(
    "key,value",
    [
        ("rendererSha256", "old"),
        ("recipeKey", "old"),
        ("descriptors", [dict(key="old")]),
    ],
)
def test_old_model_is_rejected_before_render(tmp_path, key, value):
    renderer, saved = fixture(tmp_path)
    saved["metadata"][key] = value
    (tmp_path / "search.json").write_text(json.dumps(saved))
    with pytest.raises(ValueError, match="Stale"):
        verify_candidate(renderer, tmp_path)


@pytest.mark.parametrize("name", ["reference", "candidate"])
def test_changed_audio_is_rejected(tmp_path, name):
    renderer, _ = fixture(tmp_path)
    write_wav(tmp_path / f"{name}.wav", AudioBuffer(np.zeros(2), 2))
    with pytest.raises(ValueError):
        verify_candidate(renderer, tmp_path)


def test_ui_fit_must_reproduce_candidate(tmp_path):
    renderer, _ = fixture(tmp_path)
    (tmp_path / "candidate.fit.json").write_text('{"gain":0}')
    with pytest.raises(ValueError, match="UI fit"):
        verify_candidate(renderer, tmp_path)


def test_published_audio_does_not_follow_later_search_checkpoints(tmp_path):
    import re

    renderer, _ = fixture(tmp_path)
    verify_candidate(renderer, tmp_path)
    publish_html(
        tmp_path,
        '<a href="candidate.fit.json" download>Fit</a><button data-audio="candidate">Play</button>',
    )
    html = (tmp_path / "index.html").read_text(encoding="utf8")
    audio_name = re.search('data-file="([^"]+)"', html).group(1)
    original = (tmp_path / audio_name).read_bytes()
    write_wav(tmp_path / "candidate.wav", AudioBuffer(np.zeros(2), 2))
    assert (tmp_path / audio_name).read_bytes() == original
    assert (tmp_path / "index.html").read_text(encoding="utf8") == html
