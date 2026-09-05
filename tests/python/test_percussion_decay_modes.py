import numpy as np
import pytest

from triggerfish_percussion.audio_io import AudioBuffer
from triggerfish_percussion.crash_fit_t60 import t60_diagnostics
from triggerfish_percussion.decay import (
    band_decay_fits,
    energy_decay_db,
    fit_decay,
    pre_onset_noise_power,
)
from triggerfish_percussion.erb import erb_rate, frequency_from_erb_rate
from triggerfish_percussion.descriptors import (
    contact_descriptors,
    spectral_trajectories,
)
from triggerfish_percussion.distances import (
    decay_time_distance,
    energy_matching_gain,
    erb_trajectory_distance,
    log_spectral_distance,
    log_ratio_relationship_distance,
    modal_distance,
)
from triggerfish_percussion.modes import (
    DampedMode,
    estimate_damped_modes,
    estimate_mode_count,
    match_modes,
    refit_mode_amplitudes,
    resynthesize_modes,
)
from triggerfish_percussion.modal_evidence import (
    EspritPass,
    estimate_modal_evidence,
    estimate_subband_modes,
)
from triggerfish_percussion.t60_envelope import (
    interpolate_t60,
    measure_band_t60,
    recover_two_point_t60,
)


def _coloured_noise_with_t60(
    low_seconds: float, high_seconds: float, tilt_db_per_octave: float
) -> tuple[np.ndarray, int]:
    """Build independent noise bands with the production decay curve."""
    sample_rate = 48_000
    sample_count = 6 * sample_rate
    time = np.arange(sample_count) / sample_rate
    generator = np.random.default_rng(731)
    edge_rates = np.linspace(float(erb_rate(40.0)), float(erb_rate(22_000.0)), 41)
    edges = frequency_from_erb_rate(edge_rates)
    frequencies = np.fft.rfftfreq(sample_count, 1.0 / sample_rate)
    result = np.zeros(sample_count)
    point_frequencies = np.linspace(0.0, sample_rate / 2, 8)
    point_seconds = np.linspace(low_seconds, high_seconds, 8)
    active = (True, False, False, False, False, False, False, True)
    for lower, upper in zip(edges, edges[1:]):
        selected = (frequencies >= lower) & (frequencies < upper)
        spectrum = np.zeros(frequencies.size, dtype=np.complex128)
        spectrum[selected] = generator.normal(
            size=np.count_nonzero(selected)
        ) + 1j * generator.normal(size=np.count_nonzero(selected))
        noise = np.fft.irfft(spectrum, sample_count)
        noise /= max(float(np.sqrt(np.mean(np.square(noise)))), 1.0e-12)
        center = float(
            frequency_from_erb_rate(0.5 * (erb_rate(lower) + erb_rate(upper)))
        )
        t60 = float(
            interpolate_t60(
                center,
                point_frequencies,
                point_seconds,
                active,
                sample_rate,
            )
        )
        colour = 10.0 ** (
            tilt_db_per_octave * np.log2(max(center, 40.0) / 1000.0) / 20.0
        )
        result += colour * noise * 10.0 ** (-3.0 * time / t60)
    return result / np.sqrt(len(edges) - 1), sample_rate


def test_energy_decay_recovers_a_known_t60():
    sample_rate = 1_000
    time = np.arange(2_000) / sample_rate
    energy = 10.0 ** (-6.0 * time / 0.5)
    decay = energy_decay_db(energy)
    fit = fit_decay(time, decay, -5.0, -25.0, robust=False)
    assert fit.t60_seconds == pytest.approx(0.5, rel=0.01)
    assert fit.r_squared > 0.999
    fits = band_decay_fits(time, np.vstack((energy, 10.0 ** (-6.0 * time))))
    assert [item.t60_seconds for item in fits] == pytest.approx([0.5, 1.0], rel=0.03)


@pytest.mark.parametrize("tilt_db_per_octave", (-4.0, 4.0))
def test_two_point_t60_fit_recovers_coloured_noise_envelope(tilt_db_per_octave):
    samples, sample_rate = _coloured_noise_with_t60(3.2, 0.75, tilt_db_per_octave)
    onset_sample = round(0.1 * sample_rate)
    noise = 1.0e-5 * np.random.default_rng(93).normal(size=onset_sample + samples.size)
    samples = noise + np.pad(samples, (onset_sample, 0))
    if tilt_db_per_octave > 0:
        samples *= 0.02
    fitted = recover_two_point_t60(
        samples,
        sample_rate,
        band_count=24,
        onset_sample=onset_sample,
        minimum_r_squared=0.8,
    )
    assert fitted.low_seconds == pytest.approx(3.2, rel=0.12)
    assert fitted.high_seconds == pytest.approx(0.75, rel=0.12)
    assert fitted.log_rmse < 0.12
    assert fitted.band_frequencies_hz.size >= 20


def test_t60_acceptance_gate_rejects_the_wrong_decay():
    samples, sample_rate = _coloured_noise_with_t60(3.2, 0.75, -2.0)
    onset_sample = round(0.1 * sample_rate)
    reference = np.pad(samples, (onset_sample, 0))
    time = np.arange(samples.size) / sample_rate
    too_short = samples * 10.0 ** (-3.0 * time / 0.45)
    matching = t60_diagnostics(
        AudioBuffer(reference, sample_rate), AudioBuffer(samples, sample_rate)
    )
    rejected = t60_diagnostics(
        AudioBuffer(reference, sample_rate), AudioBuffer(too_short, sample_rate)
    )
    assert matching["passed"]
    assert not rejected["passed"]


def test_t60_gate_does_not_invent_evidence_for_an_unmeasurable_reference():
    sample_rate = 48_000
    silence = AudioBuffer(np.zeros(round(0.1 * sample_rate)), sample_rate)
    diagnostics = t60_diagnostics(silence, silence)
    assert not diagnostics["passed"]
    assert not diagnostics["evaluated"]
    assert diagnostics["status"] == "insufficient_reference_evidence"


def test_band_decay_rejects_a_recording_too_short_for_the_fit_interval():
    sample_rate = 1_000
    time = np.arange(sample_rate) / sample_rate
    energy = 10.0 ** (-6.0 * time / 4.0)
    fit = band_decay_fits(time, energy[None, :], np.array([0.0]))[0]
    assert fit.status == "recording_too_short"
    assert np.isnan(fit.t60_seconds)


def test_band_decay_uses_pre_onset_noise_without_treating_tail_as_noise():
    sample_rate = 1_000
    time = np.arange(2_000) / sample_rate
    noise_power = 1.0e-7
    pre_onset = np.full((1, 100), noise_power)
    measured_noise = pre_onset_noise_power(pre_onset)
    energy = 10.0 ** (-6.0 * time / 1.0) + noise_power
    fit = band_decay_fits(time, energy[None, :], measured_noise)[0]
    assert fit.status == "measured"
    assert fit.t60_seconds == pytest.approx(1.0, rel=0.03)
    assert fit.noise_limit_seconds > 0.5


def test_silent_decay_is_finite_but_not_a_measurement():
    time = np.arange(100, dtype=float) / 1_000
    assert np.isfinite(energy_decay_db(np.zeros(time.size))).all()
    fit = band_decay_fits(time, np.zeros((1, time.size)))[0]
    assert fit.status == "insufficient_energy"
    assert np.isnan(fit.t60_seconds)


def test_t60_measurement_rejects_window_leakage_as_unobservable():
    sample_rate = 48_000
    time = np.arange(4 * sample_rate) / sample_rate
    signal = np.sin(2 * np.pi * 220.0 * time) * 10.0 ** (-3.0 * time / 2.0)
    frequencies, fits = measure_band_t60(
        signal, sample_rate, band_count=24, start_seconds=0.1
    )
    measured = np.asarray([item.status == "measured" for item in fits])
    rejected = np.asarray(
        [item.status == "below_relative_energy_floor" for item in fits]
    )
    assert np.any(measured & (frequencies < 500.0))
    assert np.all(rejected[frequencies > 2_000.0])


def test_esprit_recovers_close_damped_sinusoids():
    sample_rate = 8_000
    time = np.arange(2_048) / sample_rate
    signal = np.sin(2 * np.pi * 347.0 * time) * np.exp(-time / 0.5) + 0.4 * np.sin(
        2 * np.pi * 371.0 * time + 0.3
    ) * np.exp(-time / 0.2)
    modes = estimate_damped_modes(signal, sample_rate, 2, pencil_samples=256)
    assert [mode.frequency_hz for mode in modes] == pytest.approx(
        [347.0, 371.0], abs=0.02
    )
    assert [mode.decay_seconds for mode in modes] == pytest.approx(
        [0.5, 0.2], rel=1.0e-3
    )
    assert [mode.amplitude for mode in modes] == pytest.approx([1.0, 0.4], rel=1.0e-3)
    residual = signal - resynthesize_modes(modes, sample_rate, signal.size)
    assert np.sqrt(np.mean(residual**2)) < 1.0e-8


def test_multiresolution_modal_evidence_requires_repeatable_candidates():
    sample_rate = 8_000
    time = np.arange(2_048) / sample_rate
    clean = np.sin(2 * np.pi * 347.0 * time) * np.exp(-time / 0.5) + 0.4 * np.sin(
        2 * np.pi * 371.0 * time + 0.3
    ) * np.exp(-time / 0.2)
    generator = np.random.default_rng(9)
    hits = tuple(clean + 1.0e-4 * generator.normal(size=clean.size) for _ in range(2))
    passes = (
        EspritPass(300.0, 430.0, 2, 128),
        EspritPass(300.0, 430.0, 2, 256),
    )
    evidence = estimate_modal_evidence(hits, sample_rate, passes, merge_cents=10.0)
    accepted = [
        item for item in evidence if item.hit_count == 2 and item.observation_count == 4
    ]
    assert [item.mode.frequency_hz for item in accepted] == pytest.approx(
        [347.0, 371.0], abs=0.15
    )
    structural = tuple(item.mode for item in accepted)
    fitted = refit_mode_amplitudes(clean, sample_rate, structural)
    residual = clean - resynthesize_modes(fitted, sample_rate, clean.size)
    assert np.sqrt(np.mean(residual**2)) < 1.0e-3


def test_mdl_model_order_is_evidence_for_known_modes_in_white_noise():
    sample_rate = 8_000
    time = np.arange(2_048) / sample_rate
    generator = np.random.default_rng(12)
    signal = (
        np.sin(2 * np.pi * 347.0 * time) * np.exp(-time / 0.5)
        + 0.4 * np.sin(2 * np.pi * 711.0 * time + 0.3) * np.exp(-time / 0.2)
        + 0.01 * generator.normal(size=time.size)
    )
    evidence = estimate_mode_count(signal, 6, pencil_samples=256)
    assert evidence.selected_mode_count == 2
    assert evidence.mdl_values.shape == evidence.candidate_mode_counts.shape
    lowpass_modes = estimate_subband_modes(
        signal, sample_rate, EspritPass(0.0, 500.0, 1, 128)
    )
    assert lowpass_modes[0].frequency_hz == pytest.approx(347.0, abs=0.1)


def test_hungarian_mode_matching_reports_errors_and_unmatched_modes():
    reference = (
        DampedMode(300.0, 1.0, 1.0, 0.0),
        DampedMode(600.0, 0.5, 0.5, 0.0),
    )
    candidate = (
        DampedMode(606.0, 0.45, 0.45, 0.0),
        DampedMode(303.0, 1.1, 0.9, 0.0),
        DampedMode(1_200.0, 0.2, 0.1, 0.0),
    )
    matching = match_modes(reference, candidate)
    assert {
        (match.reference_index, match.candidate_index) for match in matching.matches
    } == {
        (0, 1),
        (1, 0),
    }
    assert matching.extra_candidate == (2,)
    assert modal_distance(matching)[-1].value == 1.0

    rejected = match_modes(
        (DampedMode(300.0, 1.0, 1.0, 0.0),),
        (DampedMode(1_200.0, 1.0, 1.0, 0.0),),
    )
    assert not rejected.matches
    assert rejected.missing_reference == (0,)
    assert rejected.extra_candidate == (0,)


def test_named_distances_are_zero_only_for_matching_descriptors():
    reference = np.array([[1.0, 0.5], [0.25, 0.125]])
    assert energy_matching_gain(reference, 0.5 * reference) == pytest.approx(2.0)
    assert log_spectral_distance(reference, reference).value == pytest.approx(0.0)
    level, change = erb_trajectory_distance(reference, reference)
    assert level.value == pytest.approx(0.0)
    assert change.value == pytest.approx(0.0)
    assert decay_time_distance(
        np.array([1.0, 2.0]), np.array([1.0, 2.0])
    ).value == pytest.approx(0.0)
    assert log_ratio_relationship_distance(
        np.array([1.0, 2.0]),
        np.array([2.0, 3.0]),
        np.array([0.5, 1.0]),
        np.array([1.0, 1.5]),
        "velocity_relation",
    ).value == pytest.approx(0.0)


def test_contact_and_dense_descriptors_have_declared_units_and_shapes():
    sample_rate = 8_000
    time = np.arange(256) / sample_rate
    contact = np.hanning(256) * np.sin(2 * np.pi * 1_000 * time)
    descriptor = contact_descriptors(contact, sample_rate)
    assert descriptor.energy > 0
    assert descriptor.spectral_centroid_hz == pytest.approx(1_000, abs=5)
    frequencies = np.array([100.0, 1_000.0, 4_000.0])
    power = np.array([[1.0, 0.5], [2.0, 2.0], [0.5, 1.0]])
    trajectories = spectral_trajectories(frequencies, power)
    assert trajectories.centroid_hz.shape == (2,)
    assert np.isfinite(trajectories.flatness).all()
