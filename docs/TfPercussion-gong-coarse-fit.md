# Gong: protected core and coarse fitting

Current main-workbench audition candidate: `gong_calibration.fit.json`, named
**Gong - protected harmonic core, coarse fit**. This replaces the previous
individually adjusted candidate; it is not a claim of improved realism.
Exact checkpoint: `build/gong-protected-final`. Previous parameters/audio are
archived in `build/gong-protected-coarse/baseline`.

## Constraints

- Generated layout: base **121.8 Hz**, **24 modes**, **harmonic core 4**,
  **stretch .7**, using the shared workbench generator. First four centres:
  121.8, 243.6, 365.4 and 487.2 Hz.
- **No individual frequency adjustments.** All 24 centres remain bit-for-bit
  equal to the generated double-precision values in the saved parameter file.
- Levels follow **four broad amplitude controls**, at 120, 600, 3000 and
  15000 Hz. Linear positive amplitudes interpolate in log frequency. There
  are no independent bar-level fits or corrections for individual ridges.
- All local noisiness multipliers stay at one. Only shared noisiness controls
  move. Damping retains the two existing endpoints; no per-mode decay fits.
- Same reference, gesture, excitation gain, output gains and radiation as the
  previous trial. No per-file normalization or velocity compression.

The four level coordinates are offline fitting constraints, not hidden DSP
parameters. The resulting ordinary observation bars are fully visible/editable.
Their smooth shape can be reconstructed to floating-point precision from the
four coordinates; there are no inherited irregular upper-bar offsets.

## Actual fitting procedure

`tools/fit_coarse_gong.py` compares entire layouts with counts 16/24/32 and
stretch .4/.7/1, all using the same 121.8-Hz root and four-mode core. Layouts
outside 15 kHz are reported and skipped. All bars initially start at −20 dB;
none of the old detailed bar pattern is copied.

Each valid layout receives a four-coordinate level fit against the absolute
STFT-band onset/bloom/decay objective. The best two layouts, ranked by full Mel
plus .05 times the bloom residual norm, receive two alternating passes:

1. Bounded finite-difference fitting of shared excitation tilt/centre,
   diffusion strength/energy dependence, noisiness/slope/spread/phase bandwidth
   and the two T60 endpoints, using actual Wasm renders.
2. Refit only the four broad observation amplitudes.

The selected layout then receives auraloss Mel refinement of those same four
coordinates over seeds 1675 and 1776. A 1-dB per-cell guard limits deterioration
of the already fitted onset/bloom/decay. Neither stage moves individual modes.

`coarse_observation_fit.py` builds an exact four-column affine basis from actual
renders. Independent amplitude probes validate it before fitting; cached STFT
cross-power retains interference. Analytic derivatives and interpolation are
unit-tested. The perceptual pass also checks its predicted score against direct
Wasm renders. No differentiable replacement synthesizer is used.

## Result and limitations

The reference remains Gong Dresden 03, first six seconds, mono downmix at
44.1 kHz, existing −6 dB family gain, strength .76 and the standard mallet event.

| Four-seed mean diagnostic; lower is better | Previous detailed fit | Coarse fit |
|---|---:|---:|
| Full Mel | 1.2398 | 1.2985 |
| First-400-ms Mel | 1.2425 | 1.3367 |
| First-400-ms linear-bin MR-STFT | 1.9825 | 2.1081 |

The standard-seed bloom residual is 6.5576 versus 6.3454 previously. This is a
deliberate simplicity trade-off, **not a lower-error calibration**. Fixed-scale
spectrogram, difference and band-envelope plots were inspected. The low decay
and delayed high bloom remain broadly aligned; the irregular 300–700-Hz
reference ridges do not match the protected harmonic series, and the upper tail
is somewhat too persistent. No isolated mode edits were added to conceal this.

Final broad levels are approximately −20.01, −14.00, −16.91 and +5.70 dB at the
four curve coordinates. T60 endpoints are 7.289 and 2.960 seconds. Diffusion
strength is 3.526, energy dependence .142 and phase-noise bandwidth .158.

Eight quarter-note hits peak at −7.03 dBFS; eight rapid full-strength hits peak
at −0.73 dBFS, before master/limiter. Output stays finite. Strength .25/.5/.75/1
checks produce increasing levels and retain delayed upper energy. These are
robustness checks, not fits to additional recorded velocities. Listening
approval remains outstanding.

Reproduce with `fit_coarse_gong.py --output ...`, then
`polish_coarse_gong.py SOURCE/candidate --output ...`. Review with
`review_metal_refit.py`, `audit_metal_strength.py` and the existing fixed-scale
plot tools. Build/serve the main workbench through `dev.ps1`; analysis remains
optional for normal Rack builds.
