import numpy as np
import pytest

import _triggerfish_dsp as dsp
from reference_ba662 import (
    BASE_CONTROL_CURRENT,
    c38_transfer,
    differential_pair_current,
    small_signal_gm,
    static_tb303_vca,
)
from reference_ba662_spice import find_ngspice, harmonic_projection


def test_spice_reference_honours_ngspice_environment_variable(tmp_path, monkeypatch):
    executable = tmp_path / "custom-ngspice"
    executable.write_bytes(b"")
    monkeypatch.setenv("NGSPICE", str(executable))

    assert find_ngspice(None) == executable.resolve()


def test_spice_residual_projection_removes_nonperiodic_solver_floor():
    period_count = 64
    size = period_count * 256
    sample = np.arange(size)
    fundamental = np.sin(2.0 * np.pi * period_count * sample / size)
    third = 0.1 * np.sin(2.0 * np.pi * 3 * period_count * sample / size)
    off_harmonic_noise = 0.5 * np.sin(2.0 * np.pi * 17 * sample / size)

    projected = harmonic_projection(fundamental + third + off_harmonic_noise)

    assert projected == pytest.approx(fundamental + third, abs=1.0e-13)


def test_ota_core_matches_exact_differential_pair_reference():
    differential = np.linspace(-0.25, 0.25, 20_001)
    control = np.linspace(0.0, 500.0e-6, differential.size)

    actual = dsp.ota_vca_current(differential, control)
    expected = differential_pair_current(differential, control)

    assert actual == pytest.approx(expected, rel=2.0e-14, abs=2.0e-18)


def test_ota_core_small_signal_gain_and_current_limit():
    differential = np.asarray((-1.0e-7, 1.0e-7))
    control = np.full(2, BASE_CONTROL_CURRENT)
    current = dsp.ota_vca_current(differential, control)
    measured_gm = (current[1] - current[0]) / (differential[1] - differential[0])

    assert measured_gm == pytest.approx(
        small_signal_gm(BASE_CONTROL_CURRENT), rel=1e-10
    )
    saturated = dsp.ota_vca_current(
        np.asarray((-10.0, 10.0)), np.full(2, BASE_CONTROL_CURRENT)
    )
    assert np.max(np.abs(saturated)) <= 0.85 * BASE_CONTROL_CURRENT


def test_tb303_wrapper_matches_static_reference_above_coupling_corner():
    sample_rate = 48_000.0
    time = np.arange(int(sample_rate)) / sample_rate
    audio = 0.1 * np.sin(2.0 * np.pi * 1000.0 * time)
    base = np.ones_like(audio)
    accent = np.zeros_like(audio)

    actual = dsp.tb303_vca(audio, base, accent, sample_rate)
    expected = static_tb303_vca(audio, base, accent)
    settled = slice(int(0.25 * sample_rate), None)

    # The production wrapper includes the 1 uF C38 output coupling capacitor;
    # this static reference deliberately does not. At 1 kHz the remaining
    # phase error is small but measurable.
    assert np.sqrt(np.mean((actual[settled] - expected[settled]) ** 2)) < 3.0e-4


def test_c38_and_50k_volume_control_set_the_output_coupling_corner():
    sample_rate = 48_000.0
    frequency = 1.0 / (2.0 * np.pi * 50_000.0 * 1.0e-6)
    time = np.arange(int(4.0 * sample_rate)) / sample_rate
    audio = 0.1 * np.sin(2.0 * np.pi * frequency * time)
    base = np.ones_like(audio)
    accent = np.zeros_like(audio)

    actual = dsp.tb303_vca(audio, base, accent, sample_rate)
    expected = static_tb303_vca(audio, base, accent)
    settled = slice(int(2.0 * sample_rate), None)
    ratio = np.sqrt(np.mean(actual[settled] ** 2)) / np.sqrt(
        np.mean(expected[settled] ** 2)
    )

    assert ratio == pytest.approx(1.0 / np.sqrt(2.0), rel=0.015)


def test_c38_nonidealities_do_not_require_another_audio_band_state():
    frequencies = np.asarray((20.0, 1000.0))
    ideal = np.abs(c38_transfer(frequencies))
    conservative = np.abs(
        c38_transfer(
            frequencies,
            source_ohms=1000.0,
            esr_ohms=10.0,
            leakage_ohms=1.0e6,
        )
    )
    # Remove the fixed level change at 1 kHz; Rack output calibration already
    # absorbs it. The remaining 20 Hz-to-1 kHz shape error is below 0.01 dB.
    relative_db = 20.0 * np.log10(
        (conservative[0] / ideal[0]) / (conservative[1] / ideal[1])
    )
    assert abs(relative_db) < 0.01


def test_tb303_wrapper_is_closed_at_zero_control_and_accent_is_additive():
    sample_rate = 48_000.0
    time = np.arange(int(0.25 * sample_rate)) / sample_rate
    audio = 5.0 * np.sin(2.0 * np.pi * 1000.0 * time)
    zero = np.zeros_like(audio)
    one = np.ones_like(audio)

    closed = dsp.tb303_vca(audio, zero, zero, sample_rate)
    base = dsp.tb303_vca(audio, one, zero, sample_rate)
    accented = dsp.tb303_vca(audio, one, one, sample_rate)

    assert np.max(np.abs(closed)) == 0.0
    assert np.sqrt(np.mean(accented**2)) > 1.8 * np.sqrt(np.mean(base**2))


def test_tb303_wrapper_nominal_distortion_matches_transistor_reference():
    sample_rate = 48_000.0
    time = np.arange(int(sample_rate)) / sample_rate
    # The wrapper maps a 5 V-peak Rack oscillator to 5 mV RMS differential
    # drive. The 20 uA control current and 220k load are the TB-303-oriented
    # operating point; the separate 200 uA/27k modern-datasheet case lives in
    # the offline transistor-level reference harness.
    audio = 5.0 * np.sin(2.0 * np.pi * 1000.0 * time)
    output = dsp.tb303_vca(
        audio, np.ones_like(audio), np.zeros_like(audio), sample_rate
    )[int(0.25 * sample_rate) :]
    windowed = output * np.hanning(output.size)
    spectrum = np.abs(np.fft.rfft(windowed))
    frequencies = np.fft.rfftfreq(output.size, 1.0 / sample_rate)
    bins = [
        np.argmin(np.abs(frequencies - harmonic * 1000.0)) for harmonic in range(1, 10)
    ]
    thd = np.sqrt(np.sum(spectrum[bins[1:]] ** 2)) / spectrum[bins[0]]

    assert 0.001 <= thd <= 0.003


def test_stock_articulation_timings_and_tied_gate():
    sample_rate = 48_000.0
    gate = np.ones(int(0.1 * sample_rate)) * 10.0
    accent = np.zeros_like(gate)
    resonance = np.zeros_like(gate)
    output = dsp.tb303_articulation(gate, accent, resonance, sample_rate=sample_rate)

    main, _, volume, _ = output.T
    assert volume[int(0.003 * sample_rate)] == 0.0
    assert volume[int(0.0075 * sample_rate)] > 0.99
    assert np.all(np.diff(main) <= 0.0)


def test_stock_release_holds_then_closes_linearly():
    sample_rate = 48_000.0
    gate = np.concatenate(
        (np.ones(int(0.02 * sample_rate)) * 10.0, np.zeros(int(0.03 * sample_rate)))
    )
    accent = np.zeros_like(gate)
    resonance = np.zeros_like(gate)
    volume = dsp.tb303_articulation(gate, accent, resonance, sample_rate=sample_rate)[
        :, 2
    ]
    release = int(0.02 * sample_rate)

    assert volume[release + int(0.007 * sample_rate)] == pytest.approx(
        volume[release], rel=2.0e-4
    )
    assert volume[release + int(0.017 * sample_rate)] == 0.0


@pytest.mark.parametrize("sample_rate", (44_100.0, 48_000.0, 96_000.0, 192_000.0))
def test_accent_branch_uses_47k_33nf_and_sweep_retains_memory(sample_rate):
    note = int(0.025 * sample_rate)
    gap = int(0.015 * sample_rate)
    gate = np.concatenate(
        (
            np.ones(note) * 10.0,
            np.zeros(gap),
            np.ones(note) * 10.0,
        )
    )
    accent = np.ones_like(gate) * 10.0
    resonance = np.ones_like(gate)
    output = dsp.tb303_articulation(gate, accent, resonance, sample_rate=sample_rate)
    filter_accent = output[:, 1]
    vca_accent = output[:, 3]

    tau_sample = round(47_000.0 * 33.0e-9 * sample_rate)
    assert vca_accent[tau_sample] == pytest.approx(1.0 - np.exp(-1.0), abs=0.02)
    assert filter_accent[note + gap] > 0.0


def test_vca_decay_control_spans_short_decay_stock_decay_and_sustain():
    sample_rate = 48_000.0
    gate = np.ones(int(sample_rate)) * 10.0
    accent = np.zeros_like(gate)
    resonance = np.zeros_like(gate)

    short = dsp.tb303_articulation(
        gate, accent, resonance, vca_decay=0.0, sample_rate=sample_rate
    )[-1, 2]
    stock = dsp.tb303_articulation(
        gate, accent, resonance, vca_decay=0.5, sample_rate=sample_rate
    )[-1, 2]
    sustain = dsp.tb303_articulation(
        gate, accent, resonance, vca_decay=1.0, sample_rate=sample_rate
    )[-1, 2]

    assert short < 1.0e-6
    assert 0.7 < stock < 0.8
    assert sustain > 0.999


def test_devil_fish_release_is_immediate_instead_of_stock_hold():
    sample_rate = 48_000.0
    note = int(0.02 * sample_rate)
    release = int(0.002 * sample_rate)
    gate = np.concatenate((np.ones(note) * 10.0, np.zeros(release)))
    accent = np.zeros_like(gate)
    resonance = np.zeros_like(gate)

    stock = dsp.tb303_articulation(gate, accent, resonance, sample_rate=sample_rate)[
        :, 2
    ]
    devil_fish = dsp.tb303_articulation(
        gate, accent, resonance, sample_rate=sample_rate, devil_fish=True
    )[:, 2]

    assert stock[-1] == pytest.approx(stock[note], rel=2.0e-4)
    assert devil_fish[-1] < 0.25 * devil_fish[note]


def test_accent_without_gate_does_not_trigger_either_envelope():
    samples = 4096
    output = dsp.tb303_articulation(
        np.zeros(samples), np.ones(samples) * 10.0, np.ones(samples)
    )

    assert np.max(np.abs(output)) == 0.0


def test_accent_sweep_panel_modes_follow_documented_behaviour():
    sample_rate = 48_000.0
    note = int(0.01 * sample_rate)
    gap = int(0.01 * sample_rate)
    gate = np.concatenate((np.ones(note) * 10.0, np.zeros(gap), np.ones(note) * 10.0))
    accent = np.ones_like(gate) * 10.0
    resonance = np.ones_like(gate)

    off = dsp.tb303_articulation(
        gate, accent, resonance, sample_rate=sample_rate, accent_sweep=0
    )[:, 1]
    fast = dsp.tb303_articulation(
        gate, accent, resonance, sample_rate=sample_rate, accent_sweep=1
    )[:, 1]
    assert np.max(np.abs(off)) == 0.0
    assert fast[0] > fast[note + gap]

    long_gate = np.ones(int(0.5 * sample_rate)) * 10.0
    long_accent = np.ones_like(long_gate) * 10.0
    long_resonance = np.ones_like(long_gate)
    normal = dsp.tb303_articulation(
        long_gate,
        long_accent,
        long_resonance,
        normal_decay=3.0,
        accent_decay=3.0,
        sample_rate=sample_rate,
        accent_sweep=2,
    )[:, 1]
    slow = dsp.tb303_articulation(
        long_gate,
        long_accent,
        long_resonance,
        normal_decay=3.0,
        accent_decay=3.0,
        sample_rate=sample_rate,
        accent_sweep=3,
    )[:, 1]
    assert slow[-1] > normal[-1]
