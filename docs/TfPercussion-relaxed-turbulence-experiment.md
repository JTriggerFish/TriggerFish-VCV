# Relaxed turbulence: crash and gong experiment

This is an opt-in control-mapping experiment, not a replacement calibration.
All four metallic factory fits explicitly keep `field_relaxed_turbulence: 0`.
The saved user Ride remains unchanged in sound and gesture.

## Mapping

Classic remains available and retains its existing mapping, including its
global ceiling before the local multiplier. In the relaxed mapping, intensity
at a painted handle frequency is

$$I_i=L\left(\frac{f_i}{f_c}\right)^s m_i.$$

$L$ is the global level, $s$ the slope, $f_c$ the centre and $m_i$ the existing
local multiplier. The global slider covers 0–1 in classic and 0–4 in relaxed
mode. In relaxed mode +1 slope doubles intensity per octave; classic slope
remains additive. Global zero or local zero means no turbulence. There is no
clamp at intensity one, before or after the local multiplier.

This is a different curve family: identical knob values are deliberately
included as a diagnostic, not asserted to be equivalent settings.

| Quantity | Relaxed mapping |
|---|---|
| Fraction of drive energy assigned to satellites | $q_i=1-\exp[-\ln(10)I_i^2]$ |
| Centre/satellite drive amplitudes | $\sqrt{1-q_i}$ and $\sqrt{q_i/(2N_i)}$ |
| Packet spread | $I_i$ × existing ERB spread control |
| Phase linewidth | $I_i^2$ × existing ERB bandwidth and existing centre/satellite scale |
| Exchange and optional smooth-drift weight | $q_i$, bounded between zero and one |

At intensity one, satellites carry 90% of their packet's drive energy; the
fraction approaches 100% smoothly. With no satellites allocated, all drive
energy stays at the centre. Excitation remains unit-normalized regardless of
packet allocation, and the total state budget is still 512. This is a drive
partition, not a promise that interference or repeated coherent excitation
cannot change observed loudness.

Spread/linewidth can continue increasing while the energy fraction remains
bounded. Existing positive-frequency/Nyquist limits still apply. At extreme
settings, packets can crowd those limits; this is not an unlimited useful
range or a guarantee of realistic sound. Wider satellites also sample different
positions on the unchanged T60 curve, so the observed decay need not be identical.

The small mapping helper runs during preparation, not per audio sample.
The modal editor now uses the same slope, centre and mapping for its drawn
widths and width dragging. Previously its drawing ignored slope and centre.
Switching from relaxed to classic explicitly caps a global value above one
back to one; the slider has no inactive upper range in classic mode.

## Comparisons

Existing workbench URLs:

- `/?audition=turbulence-check/crash/manifest.json`
- `/?audition=turbulence-check/gong/manifest.json`

Each has the real reference and four synthetic variants, with single and
four repeated hits at 500 ms intervals. Synthesized clips are eight seconds
at 48 kHz. Reference files are SHA-verified, mono-folded, onset-aligned and
given only their saved reference gain. No loudness/peak normalization is used.
Playback uses the usual master and limiter; no new server is started.

1. Saved settings, classic mapping.
2. Saved settings, relaxed mapping.
3. Stronger contrast settings, classic mapping.
4. The same stronger contrast settings, relaxed mapping.

| Contrast setting | Level | Slope | Centre |
|---|---:|---:|---:|
| Crash | 0.7 | 0.65 | 2500 Hz |
| Gong | 0.4 | 1 | 2000 Hz |

Only these three controls and the mapping switch vary. Painted handles/local
multipliers, T60 curve, contact, cascade, density, exchange knob, observation
and output gains remain fixed. Smooth drift remains off. The paired contrast
trials distinguish a benefit of the new mapping from ordinary knob adjustment.
These are manually chosen probes, not optimized fits.

## Inspection and acceptance

The baseline crash renders match the previous classic renderer bit-for-bit.
Within each classic/relaxed pair, eight-second RMS changes are under 0.06 dB
for single hits and under 0.4 dB for repeated hits; no correction is applied.

The reference-scaled STFT plots do **not** establish a successful fit:

- Gong: relaxed mapping at the saved knobs makes the low region more diffuse.
  Stronger contrast restores clearer low ringing with upper wash, but the
  reference's bloom development and decay are still different. Simply enabling
  relaxed mapping is not an improvement we should accept automatically.
- Crash: the relaxed contrast trial redistributes upper texture, but the
  lower/middle spectral structure and decay remain mismatched. It does not
  remove the need to tune the instrument, nor diagnose the remaining harshness.

The judgement pending is whether the mapping provides a more useful audible
range and easier low/high separation. No trial is promoted to a factory fit.

Tests cover continuity/monotonicity, bounded energy allocation, clean modes,
normalized modal drive across density/level settings, finite passive tails,
prepared-Wasm replay, real compiled preset descriptors and editor curve values.
All 20 percussion tests, both native workbench API tests and all 10 Wasm/site
CTest checks pass, including the native/Wasm signature comparison.

Private scripts, exact snapshots/overrides, renderer hashes, WAVs, measurements
and inspected Plotly STFTs are in `build/turbulence-v12/`.
