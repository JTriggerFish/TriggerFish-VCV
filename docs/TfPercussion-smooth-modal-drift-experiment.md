# Smooth modal drift: crash texture experiment

## Question and controlled comparison

Does slowly correlated frequency motion sound less synthetic than the current
independent per-sample phase kicks? This is not a new calibration or a claim
that modulation fixes the crash. The preceding 1x/4x audition did not identify
oversampling as a useful direction, so all these renders use 48 kHz, 1x.

The existing workbench serves the comparison at
`/?audition=drift-check/manifest.json`. The saved graded-turbulence crash is the
starting point. Each variant has a seven-second single strike and four strikes
500 ms apart, with persistent state. The actual reference is also available.

| Variant | Phase diffusion | Smooth drift depth | Drift speed |
|---|---:|---:|---:|
| Current | 0.35 ERB | 0% | irrelevant |
| Phase diffusion off (control) | 0 ERB | 0% | irrelevant |
| Smooth, moderate | 0 ERB | 2% | 8 knots/s |
| Smooth, stronger | 0 ERB | 6% | 8 knots/s |

All modal centre frequencies, bars, T60 points, contact, cascade, neighbour
exchange, packet allocation, turbulence profile, observation filters and gains
remain fixed. No loudness normalization is used. Playback goes through the
usual master and limiter. The audition buttons play fixed renders, not the
current edited patch. Single-hit RMS differs by less than 0.44 dB between
synthetic variants; repeated-hit RMS differs by less than 0.22 dB. Radiated
levels need not be identical when phase relationships change.

## Implementation and explicit controls

`SmoothModalDrift` supplies independent random trajectories per modal state,
not a shared vibrato. Uniform targets in [-1, 1] are joined using

$$
h(t)=6t^5-15t^4+10t^3,\qquad 0\leq t\leq 1.
$$

Values and their first two derivatives are continuous at target boundaries.
Independent initial positions stagger those boundaries. The random stream is
separate from existing phase/exchange randomness and resets with the field,
not with each strike. The process has zero mean in the ensemble; an individual
short hit can have a temporary positive or negative frequency offset.

For centre frequency $f_i$, global depth $d$ in percent, packet turbulence
$T_i$, and random curve $u_i(t)$, the frequency deviation is

$$
\Delta f_i(t)=f_i\frac{d}{100}T_i^2u_i(t).
$$

This reuses the existing squared-turbulence weighting already used for modal
diffusion/exchange. Thus a clean packet remains fixed, whereas turbulent
upper packets wander more. Excursion is bounded near Nyquist. Depth describes
frequency deviation, **not spectral linewidth**, and speed is random targets
per second, **not a periodic LFO frequency**.

The added phase increment rotates the modal recurrence coefficients using a
unit-length Cayley rotation. A bounded small-angle approximation avoids
per-sample sine/cosine calls. Pole radii are unchanged: drift adds no independent
energy source or decay envelope. It can change interference, excitation energy
accumulation and the observed spectrum without changing ordinary modal damping.

Two controls are exposed in the modal field's advanced section and serialize
through the ordinary parameter/JSON path:

- `field_drift_depth`: 0–10%, default 0 (disabled).
- `field_drift_rate`: 0.1–40 knots/s, default 8.

`field_phase_bandwidth` remains independent. To reproduce the smooth replacement
trials, set it to zero; otherwise drift is added to the existing diffusion.
No saved instrument preset is changed by this experiment. These are experimental
controls, not a decision to permanently expand the instrument's control surface.

## Checks and observations

- Default-off single and repeated renders are bit-identical to the previous
  Wasm output (all 336,000 samples per clip).
- Primitive tests cover bounded/continuous trajectories, clean-mode protection,
  deterministic reset, passive free decay and declared T60.
- Prepared native/Wasm replay and block-size tests cover the enabled path.
- All 19 percussion tests and both native workbench API tests pass. Wasm smoke,
  calibration and diagnostic-player tests pass. A disposable browser tab
  reached Ready, displayed both controls and all nine audition buttons; no
  test triggered audio or modified the user's open patch.
- Inspection of 0.2–1.5 s and 1.5–3 s spectra shows that removing phase kicks
  exposes narrow spikes; smooth drift broadens them while reducing the broad
  upper-frequency skirts of the existing phase diffusion. This does not prove
  that the resulting timbre is preferable. Lower/middle reference mismatches
  remain and are not addressed by this experiment.

Private reproduction files, the exact snapshot/overrides, renderer hash, audio,
measurements and inspected Plotly figure are in `build/crash-drift-v11/`.
The saved reference WAV there is already onset-aligned; do not remove its
original 50 ms pre-roll a second time. No new server is started.

### Preset-loader regression correction

The initial browser check above only tested page initialization; it missed
factory preset selection. The new controls increased the compiled surface
from 149 to 151 parameters, while the four metallic factory files still had
149. The strict calibration loader correctly rejected those incomplete files.
All four now explicitly contain drift depth 0 and speed 8; existing sound
parameters are unchanged. The new `calibration_surface_tests.mjs` validates
all six factory targets against real Wasm descriptors, rather than deriving
mock descriptors from the same JSON being tested. It is registered in CTest.
The disposable browser check now selects all four metallic targets as well.

## Inspiration, not a copied cymbal model

[Lokki and Hiipakka, DAFx 2001](https://users.aalto.fi/~ktlokki/Publs/lokki_dafx01.pdf)
describe continuous modulation to move resonant peaks in a reverberation
enhancement system, including different modulation trajectories across
channels. That supports testing time variation against persistent coloration.
The smooth random modal-frequency construction above is our experimental
translation of that principle, not their algorithm or evidence of a physically
correct cymbal model.
