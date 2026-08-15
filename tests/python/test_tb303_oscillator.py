import numpy as np

import _triggerfish_dsp as dsp

SAMPLE_RATE = 48_000.0


def render(
    function,
    pitch,
    *,
    seconds=2.0,
    slide=0.0,
    fm=0.0,
    shape=0.0,
    wave=0.0,
    slide_time=0.060,
    linear_fm=False,
):
    samples = int(SAMPLE_RATE * seconds)

    def signal(value):
        if np.isscalar(value):
            return np.full(samples, value, dtype=np.float64)
        return np.asarray(value, dtype=np.float64)

    return function(
        signal(pitch),
        signal(slide),
        signal(fm),
        signal(shape),
        signal(wave),
        sample_rate=SAMPLE_RATE,
        slide_time=slide_time,
        linear_fm=linear_fm,
    )


def dominant_frequency(signal, low, high):
    windowed = signal * np.hanning(signal.size)
    spectrum = np.abs(np.fft.rfft(windowed))
    frequencies = np.fft.rfftfreq(signal.size, 1.0 / SAMPLE_RATE)
    selection = (frequencies >= low) & (frequencies <= high)
    return frequencies[selection][np.argmax(spectrum[selection])]


def test_x4_oscillator_tracks_one_volt_per_octave():
    rendered = render(dsp.tb303_oscillator_x4, 0.0)
    frequency = dominant_frequency(rendered[-48_000:, 0], 200.0, 320.0)
    assert abs(frequency - 261.625565) < 0.5


def test_stock_slide_uses_the_22_ms_pitch_cv_time_constant():
    samples = 48_000
    pitch = np.ones(samples)
    pitch[:4_800] = 0.0
    slide = np.zeros(samples)
    slide[4_800:] = 1.0
    rendered = render(
        dsp.tb303_oscillator_x4,
        pitch,
        seconds=1.0,
        slide=slide,
        slide_time=0.060,
    )
    after_sixty_ms = rendered[4_800 + 2_880 - 1, 3]
    expected = 1.0 - np.exp(-0.060 / 0.022)
    assert abs(after_sixty_ms - expected) < 2e-3


def test_wave_morph_reaches_both_waveform_endpoints():
    saw = render(dsp.tb303_oscillator_x4, 0.0, wave=0.0)
    square = render(dsp.tb303_oscillator_x4, 0.0, wave=1.0)
    settled = slice(24_000, None)
    np.testing.assert_allclose(saw[settled, 2], saw[settled, 0], atol=1e-6)
    np.testing.assert_allclose(square[settled, 2], square[settled, 1], atol=1e-6)


def test_wave_morph_does_not_cancel_between_endpoints():
    fundamental_levels = []
    rms_levels = []
    for wave in np.linspace(0.0, 1.0, 11):
        mixed = render(dsp.tb303_oscillator_x4, np.log2(85.0 / 261.625565), wave=wave)[
            -48_000:, 2
        ]
        window = np.hanning(mixed.size)
        spectrum = np.abs(np.fft.rfft(mixed * window))
        frequencies = np.fft.rfftfreq(mixed.size, 1.0 / SAMPLE_RATE)
        fundamental_levels.append(spectrum[np.argmin(np.abs(frequencies - 85.0))])
        rms_levels.append(np.sqrt(np.mean(mixed**2)))

    assert np.all(np.diff(fundamental_levels) > 0.0)
    assert min(rms_levels) > 0.95 * min(rms_levels[0], rms_levels[-1])


def test_linear_fm_runs_backwards_through_zero():
    # C4 with -2 V at the 200 Hz/V linear input produces -138.37 Hz. A
    # clamp-to-zero implementation would collapse to an almost-static signal.
    rendered = render(
        dsp.tb303_oscillator_x4,
        0.0,
        fm=-2.0,
        linear_fm=True,
    )
    saw = rendered[-48_000:, 0]
    frequency = dominant_frequency(saw, 100.0, 180.0)

    assert abs(frequency - (400.0 - 261.625565)) < 0.5
    assert np.std(saw) > 2.0


def test_stock_square_duty_decreases_with_pitch():
    collector_high_duties = []
    for frequency in (10.0, 85.0, 1_000.0):
        pitch = np.log2(frequency / 261.625565)
        square = render(dsp.tb303_oscillator_x4, pitch, wave=1.0)[-48_000:, 1]
        # The Rack-facing square is polarity-aligned with the saw for a
        # cancellation-free morph. Its negative portion corresponds to Q8's
        # collector-high portion used by the hardware measurements.
        collector_high_duties.append(float(np.mean(square < 0.0)))

    assert (
        collector_high_duties[0] > collector_high_duties[1] > collector_high_duties[2]
    )
    assert 0.62 < collector_high_duties[0] < 0.75
    assert 0.43 < collector_high_duties[1] < 0.53
    assert 0.36 < collector_high_duties[2] < 0.47


def test_shape_control_extends_the_square_bias_in_both_directions():
    low = render(dsp.tb303_oscillator_x4, 0.0, shape=-0.7, wave=1.0)
    high = render(dsp.tb303_oscillator_x4, 0.0, shape=0.7, wave=1.0)
    low_duty = np.mean(low[-48_000:, 1] < 0.0)
    high_duty = np.mean(high[-48_000:, 1] < 0.0)
    assert high_duty > low_duty + 0.3


def test_oversampling_variants_remain_finite_at_extended_pitch():
    pitch = np.log2(12_000.0 / 261.625565)
    for function in (
        dsp.tb303_oscillator_x1,
        dsp.tb303_oscillator_x2,
        dsp.tb303_oscillator_x4,
    ):
        rendered = render(function, pitch, seconds=0.1, shape=1.0, wave=1.0)
        assert np.isfinite(rendered).all()
        assert np.max(np.abs(rendered[:, :3])) < 12.0
