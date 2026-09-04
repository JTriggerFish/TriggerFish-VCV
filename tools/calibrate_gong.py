"""Audit and bootstrap the constructive gong model against one reference hit.

This intentionally fits the temporal spectral trajectory before absolute level:
the defining target is a low, tonal onset whose stored energy progressively
reappears as a broader high-frequency response.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

from triggerfish_percussion.audio_io import (
    AudioBuffer,
    read_wav,
    resample_audio,
    write_wav,
)
from triggerfish_percussion.crash_model import (
    CrashEvent,
    CrashFit,
    render_crash_components,
    render_crash_sequence,
)
from triggerfish_percussion.transforms import StftConfig, stft

BANDS_HZ = ((40, 300), (300, 1000), (1000, 3000), (3000, 8000), (8000, 16000))
TIMES_SECONDS = (0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 2.5, 4.0)
WORKBENCH_OUTPUT_CALIBRATION = 174.0
GONG_MODEL_LEVEL_DB = -35.0


def default_reference() -> Path:
    return Path(os.environ.get("LOCALAPPDATA", "")) / (
        "Bitwig Studio/installed-packages/5.0/Bitwig/"
        "Acoustic Drums and Percussion/Percussion/Gongs/Gong Dresden 03.wav"
    )


def power_transform(audio: AudioBuffer):
    return stft(
        audio.samples,
        audio.sample_rate,
        StftConfig(window_samples=4096, hop_samples=512),
    )


def trajectory(audio: AudioBuffer) -> dict[str, object]:
    transform = power_transform(audio)
    power = transform.power
    frequencies = transform.frequencies_hz
    total_peak = max(float(np.max(np.sum(power, axis=0))), np.finfo(float).tiny)
    rows = []
    for target_time in TIMES_SECONDS:
        if target_time > transform.times_seconds[-1]:
            continue
        frame = int(np.argmin(np.abs(transform.times_seconds - target_time)))
        frame_power = power[:, frame]
        total = max(float(np.sum(frame_power)), np.finfo(float).tiny)
        band_db = []
        for low, high in BANDS_HZ:
            selected = (frequencies >= low) & (frequencies < high)
            value = float(np.sum(frame_power[selected])) / total_peak
            band_db.append(10.0 * np.log10(max(value, 1.0e-15)))
        rows.append(
            {
                "time_seconds": target_time,
                "band_db_relative_to_peak": band_db,
                "centroid_hz": float(np.sum(frequencies * frame_power) / total),
            }
        )
    band_peak_times = []
    for low, high in BANDS_HZ:
        selected = (frequencies >= low) & (frequencies < high)
        envelope = np.sum(power[selected], axis=0)
        band_peak_times.append(float(transform.times_seconds[np.argmax(envelope)]))
    return {"band_peak_times_seconds": band_peak_times, "rows": rows}


def low_ridges(audio: AudioBuffer) -> list[dict[str, float]]:
    transform = power_transform(audio)
    power = transform.power
    selected_time = (transform.times_seconds >= 0.03) & (
        transform.times_seconds <= 0.25
    )
    spectrum = np.mean(power[:, selected_time], axis=1)
    frequencies = transform.frequencies_hz
    selected_frequency = (frequencies >= 40.0) & (frequencies <= 2500.0)
    indices = np.flatnonzero(selected_frequency)
    db = 10.0 * np.log10(np.maximum(spectrum[indices], np.finfo(float).tiny))
    peaks, properties = find_peaks(db, prominence=3.0, distance=3)
    ordered = sorted(
        zip(peaks, properties["prominences"]),
        key=lambda item: db[item[0]],
        reverse=True,
    )[:12]
    return [
        {
            "frequency_hz": float(frequencies[indices[index]]),
            "level_db_relative": float(db[index] - np.max(db)),
            "prominence_db": float(prominence),
        }
        for index, prominence in ordered
    ]


def bootstrap_fit() -> CrashFit:
    frequencies = (
        128.9,
        245.0,
        304.7,
        375.0,
        421.9,
        550.8,
        621.1,
        668.0,
        726.6,
        878.9,
        1100.0,
        1380.0,
        1699.2,
        1980.5,
        2400.0,
        2900.0,
        3500.0,
        4200.0,
        5000.0,
        6000.0,
        7200.0,
        8800.0,
        11000.0,
        15000.0,
    )
    levels_db = (
        -23.0,
        -17.0,
        -10.0,
        0.0,
        -14.0,
        4.0,
        -12.8,
        -17.3,
        -13.7,
        -14.3,
        -10.0,
        -8.0,
        -6.0,
        -6.0,
        -5.0,
        -4.0,
        -3.0,
        -2.0,
        -1.0,
        0.0,
        0.0,
        -2.0,
        -6.0,
        -12.0,
    )
    return replace(
        CrashFit(),
        sparse_frequency_hz=frequencies,
        sparse_amplitude=tuple(10.0 ** (level / 20.0) for level in levels_db),
        field_turbulence_scale=tuple(
            0.12 if frequency < 500.0 else 0.3 if frequency < 3500.0 else 1.2
            for frequency in frequencies
        ),
        body_decay_seconds=(12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 1.7),
        body_tilt_db_per_octave=-10.0,
        body_tilt_centre_hz=550.0,
        bloom_rate_octaves_per_second=5.0,
        bloom_energy_dependence=1.0,
        bloom_phase_diffusion=0.8,
        body_excitation_gain=1.0,
        field_gain=0.08,
        field_turbulence=0.72,
        field_turbulence_slope_per_octave=0.22,
        field_turbulence_centre_hz=1200.0,
        field_packet_spread_erb=5.5,
        field_phase_bandwidth_erb=0.8,
        field_exchange=0.5,
        contact_duration_scale=2.2,
        contact_pulse_gain=0.79,
        contact_chirp_gain=1.426585,
        contact_chirp_frequency_scale=0.15,
        contact_noise_gain=0.309017,
        contact_noise_duration_scale=2.2,
        contact_micro_gain=0.309017,
        contact_micro_duration_scale=2.2,
        direct_gain=0.002,
        output_gain=WORKBENCH_OUTPUT_CALIBRATION * 10.0 ** (GONG_MODEL_LEVEL_DB / 20.0),
        direct_high_cut_hz=3000.0,
        body_low_cut_hz=25.0,
        body_colour_frequency_hz=8500.0,
        body_colour_gain_db=2.0,
        body_high_cut_hz=12000.0,
        velocity_brightness_db_per_octave=5.0,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, default=default_reference())
    parser.add_argument("--output", type=Path, default=Path("build/gong-calibration"))
    parser.add_argument("--seconds", type=float, default=6.0)
    arguments = parser.parse_args()

    reference = resample_audio(read_wav(arguments.reference).mono(), 48_000)
    count = min(
        reference.samples.size, round(arguments.seconds * reference.sample_rate)
    )
    reference = AudioBuffer(reference.samples[:count], reference.sample_rate)
    fit = bootstrap_fit()
    components = render_crash_components(
        fit,
        arguments.seconds,
        sample_rate=reference.sample_rate,
        strength=0.76,
        location=0.55,
        hardness=0.35,
        implement=0.5,
        contact_spread=0.3,
        seed=3,
    )
    synthesis = AudioBuffer(components[:, 3], reference.sample_rate)
    repeated = AudioBuffer(
        render_crash_sequence(
            fit,
            arguments.seconds,
            tuple(
                CrashEvent(
                    onset,
                    strength=0.76,
                    location=0.55,
                    hardness=0.35,
                    implement=0.5,
                    contact_spread=0.3,
                    seed=3 + index,
                )
                for index, onset in enumerate((0.0, 1.25, 2.5, 3.75))
                if onset < arguments.seconds
            ),
            reference.sample_rate,
        ),
        reference.sample_rate,
    )

    arguments.output.mkdir(parents=True, exist_ok=True)
    write_wav(arguments.output / "reference.wav", reference)
    write_wav(arguments.output / "candidate.wav", synthesis)
    write_wav(arguments.output / "candidate-repeated.wav", repeated)
    write_wav(
        arguments.output / "candidate-contact.wav",
        AudioBuffer(components[:, 0], reference.sample_rate),
    )
    write_wav(
        arguments.output / "candidate-body.wav",
        AudioBuffer(components[:, 2], reference.sample_rate),
    )
    result = {
        "reference": str(arguments.reference),
        "reference_low_ridges": low_ridges(reference),
        "candidate_low_ridges": low_ridges(synthesis),
        "reference_trajectory": trajectory(reference),
        "candidate_trajectory": trajectory(synthesis),
        "candidate_contact_trajectory": trajectory(
            AudioBuffer(components[:, 0], reference.sample_rate)
        ),
        "candidate_body_trajectory": trajectory(
            AudioBuffer(components[:, 2], reference.sample_rate)
        ),
        "fit": asdict(fit),
    }
    (arguments.output / "audit.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
