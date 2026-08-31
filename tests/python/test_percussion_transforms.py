from pathlib import Path

import numpy as np
import pytest

from triggerfish_percussion.audio_io import AudioBuffer
from triggerfish_percussion.comparison import ComparisonConfig, compare_audio
from triggerfish_percussion.erb import ErbFilterbank, erb_rate, frequency_from_erb_rate
from triggerfish_percussion.transform_cache import TransformCache
from triggerfish_percussion.transforms import StftConfig, multiresolution_stft, stft


def test_erb_rate_round_trip():
    frequencies = np.array([0.0, 40.0, 1_000.0, 8_000.0, 20_000.0])
    assert frequency_from_erb_rate(erb_rate(frequencies)) == pytest.approx(frequencies)


def test_stft_finds_a_bin_centred_sine_at_calibrated_level():
    sample_rate = 48_000
    fft_samples = 1_024
    frequency = 48 * sample_rate / fft_samples
    time = np.arange(8_192) / sample_rate
    signal = 0.7 * np.sin(2 * np.pi * frequency * time)
    result = stft(signal, sample_rate, StftConfig(fft_samples, 256, fft_samples))
    middle = result.magnitude[:, result.magnitude.shape[1] // 2]
    peak_bin = int(np.argmax(middle))
    assert result.frequencies_hz[peak_bin] == frequency
    assert middle[peak_bin] == pytest.approx(0.7, rel=2.0e-3)
    assert np.max(result.magnitude_db()) == pytest.approx(0.0)


def test_erb_projection_conserves_selected_spectral_energy():
    frequencies = np.fft.rfftfreq(2_048, 1 / 48_000)
    bank = ErbFilterbank.create(frequencies, 80.0, 18_000.0, 32)
    generator = np.random.default_rng(4)
    power = np.square(generator.normal(size=(frequencies.size, 23)))
    bands = bank.apply_power(power)
    assert np.sum(bands, axis=0) == pytest.approx(
        np.sum(power[bank.selected_bins], axis=0), rel=1.0e-12
    )
    assert np.sum(bank.weights[:, bank.selected_bins], axis=0) == pytest.approx(1.0)


def test_multiresolution_and_versioned_cache_share_canonical_transform(tmp_path: Path):
    audio = AudioBuffer(np.sin(2 * np.pi * 440 * np.arange(4_096) / 48_000), 48_000)
    configs = (StftConfig(256, 64, 512), StftConfig(1_024, 256, 2_048))
    results = multiresolution_stft(audio.samples, audio.sample_rate, configs)
    assert len(results) == 2
    cache = TransformCache(tmp_path)
    first = cache.stft(audio, configs[0])
    second = cache.stft(audio, configs[0])
    assert second.spectrum == pytest.approx(first.spectrum)
    assert len(list(tmp_path.glob("*.npz"))) == 1


def test_comparison_contract_aligns_once_and_reuses_axes_for_losses():
    sample_rate = 48_000
    reference = np.zeros(8_192)
    event = np.hanning(512) * np.sin(2 * np.pi * 2_000 * np.arange(512) / sample_rate)
    reference[1_000:1_512] = event
    candidate = np.zeros_like(reference)
    candidate[1_027:1_539] = event
    config = ComparisonConfig((StftConfig(512, 128, 1_024),), level_mode="preserve")
    comparison = compare_audio(
        AudioBuffer(reference, sample_rate), AudioBuffer(candidate, sample_rate), config
    )
    assert comparison.alignment.candidate_lag_samples == pytest.approx(-27, abs=2)
    assert comparison.resolutions[0].reference.times_seconds == pytest.approx(
        comparison.resolutions[0].candidate.times_seconds
    )
    assert max(loss.value for loss in comparison.losses) < 1.0e-10
    assert {loss.name.split("/")[1] for loss in comparison.losses} == {
        "contact",
        "bloom",
        "early_body",
    }


def test_comparison_regions_require_stft_window_support():
    sample_rate = 48_000
    reference = np.zeros(48_000)
    onset = 2_000
    contact_samples = round(0.015 * sample_rate)
    reference[onset : onset + contact_samples] = np.hanning(contact_samples)
    candidate = reference.copy()
    bloom_start = onset + contact_samples
    bloom_samples = round(0.08 * sample_rate)
    time = np.arange(bloom_samples) / sample_rate
    candidate[bloom_start : bloom_start + bloom_samples] += np.sin(
        2 * np.pi * 4_000 * time
    )
    config = ComparisonConfig(
        (
            StftConfig(256, 64, 512),
            StftConfig(1_024, 256, 2_048),
        )
    )
    comparison = compare_audio(
        AudioBuffer(reference, sample_rate), AudioBuffer(candidate, sample_rate), config
    )
    names = {loss.name for loss in comparison.losses}
    assert "log_spectral/contact/w256" in names
    assert "log_spectral/contact/w1024" not in names
    contact = next(
        loss.value
        for loss in comparison.losses
        if loss.name == "log_spectral/contact/w256"
    )
    bloom = next(
        loss.value
        for loss in comparison.losses
        if loss.name == "log_spectral/bloom/w256"
    )
    assert contact == pytest.approx(0.0)
    assert bloom > 1.0
