# Gong: pitched body followed by metallic sizzle

## Correction: preserve the approved bass

The initially published fit below was rejected by listening: its aggregate
reference scores hid a destructive change to the approved low body. The 120 Hz
bar had fallen from −14.86 to −54.73 dB, and the 240 Hz bar from −17.34 to
−22.71 dB. Raising 360–480 Hz did not preserve the same pitched body. Calling
this a successful body-preserving refinement was wrong.

The current workbench preset is **Gong — restored low body and metallic sizzle**.
Only the two low observation levels change from that rejected fit: 120 Hz is
−16.886 dB and 240 Hz is −19.683 dB. These account for the new excitation and
transport: measured 80–180/180–300 Hz power over .05–.5 seconds matches the
previous approved preset within .09 dB on both tested seeds. They are deliberately
matched to the approved sound, not suppressed to improve the reference score.
All upper controls, frequencies, other bars and damping remain unchanged.

Upper band-power changes above 3.2 kHz over .5–4 s stay below .013 dB on the
standard render. At 48 kHz the restored single hit peaks at −3.53 dBFS and four
half-second-spaced hits at −2.99 dBFS, without limiting at Master 0 dB.
Artifacts: `build/gong-low-body-restored`; correction tool:
`tools/restore_gong_low_body.py`. The main workbench save/load probe passes.

The temporary below-300 Hz bar lock was an inadequate workaround and has been
removed. [The objective/search-space review](TfPercussion-fitting-objective-failure.md)
reproduces why the fitter removed the fundamental: broad-band pitch substitution,
hard-coupled smooth correction coordinates, and a changed fitting target. The
restored preset above remains unchanged during that investigation. The historical
scores below describe the rejected balance, not approval of the current sound.

## Original refinement record

This pass tunes the existing engine. It does **not** add a delayed noise layer,
change modal frequencies, change the reference gesture, or turn up the master.
The target is the standard Dresden Gong03 reference at its saved −6 dB gain.
The saved gesture is strength 0.76, location 0.55, hardness 0.35, implement 0.5,
contact spread 0.3 and seed 1675. The user's stretched harmonic centres remain
fixed; painted observation amplitudes may change smoothly across frequency.

## What is fitted

| Stage | Variables | Solver |
|---|---|---|
| Temporal layering | Diffusion strength, concentration dependence, energy sensitivity | Coarse and local actual-WASM screens |
| Broad spectral balance | Seven or eight smooth frequency-curve values, converted to explicit painted levels | Bounded Powell search on a validated observation basis |
| Upper texture | Packet spread/density, beat controls, bounded movement and phase blur | Two-round coordinate comparison with timing/body/ridge guards |
| Decay check | Only the two active T60 endpoints, jointly with transport | Actual-WASM screen, then observation refit |
| Upper energy headroom | Existing Body excitation and transport controls | Actual-WASM screen, then observation refit |

No per-mode decay multiplier, extra damping knot, frequency placement, EQ or
output gain is fitted. The curve coordinates are optimizer coordinates only:
the JSON and UI contain the resulting ordinary painted levels, with no hidden
curve applied by C++. The absolute-curve alternative is also tested against
smooth adjustments to the incoming bars; a lower score alone is not acceptance.

## Temporal and spectral objective

`LayeredBandLoss` measures an actual 4096-point STFT with 10 ms hops. Frequency
edges are 80, 250, 500, 800, 1250, 2000, 3200, 5000, 7000, 9000, 11500 and
15000 Hz. Time boundaries are 0, .05, .1, .2, .3, .4, .5, .65, .8, 1, 1.25,
1.5, 2, 3, 4 and 6 seconds. Each cell is mean band power, converted to dB.
An analysis floor is fixed at 45 dB below **the reference's** peak in that band.
Neither audio signal is normalized.

The loss combines absolute cell errors and half-weight within-band envelope
shape errors. Shape removes a band's mean error; it does not rescale audition
audio. A shape-only first screen is explicitly an approximation: later fitting
must realize the levels with legal observation amplitudes. Upper bars reaching
their +6 dB limit are a reason to examine energy delivery, not to invent gain.

Both uniform and reference-weighted cells are tested. For the weighted variant,
with reference cell level $R_{bt}$ and band peak $P_b$:

$$
w_{bt}=\max\left(0.05,10^{(R_{bt}-P_b)/20}\right).
$$

Normalize weights to mean one within each band. Fit weighted absolute errors
and weighted, mean-removed shape errors. This keeps very quiet onset/noise-floor
cells from dominating an audible peak, while retaining an onset penalty. It is
a stated fitting heuristic, **not a validated perceptual-equivalence loss**.
Scores from different weighting variants must not be compared directly.

## Efficient fitting without a surrogate synthesizer

For fixed dynamics, observation is affine in the positive modal amplitudes.
`ObservationBasis` obtains its columns from the actual WASM renderer and validates
mixed-level predictions. `SpectralBloomBasis` then stores each cell's STFT
cross-power matrix. A candidate amplitude vector gives exact cached cell power
through a quadratic form, so thousands of optimizer evaluations require no DSP
render. Frequencies, excitation, damping, allocation and transport cannot change
inside that solve.

Every selected result is rendered again. Predicted versus actual cell levels
must agree within 0.02 dB. Reused caches require matching non-observation
parameters and validated renderer/reference provenance. Saved JSON must reproduce
the actual candidate WAV through the current renderer before publication.

## Texture and acceptance checks

The texture search uses the reference's upper spectral centre and 8–32/32–128 Hz
envelope fluctuation power near 12 kHz, after dividing out a 100 ms-smoothed
trend **for analysis only**. This is a useful local diagnostic, not sufficient
coverage by itself. Check additional auditory bands and fresh random seeds
afterwards. Fine spectral contrast in 3–7 and 7–14 kHz guards against replacing
moving ridges with smooth noise. Timing and low-body guards reject texture
changes that undo the fitted layering.

Inspect absolute band-envelope plots, same-scale reference/model STFTs and their
difference. Half-peak-power crossings are rise landmarks, not physical onsets.
Independent Mel, multi-resolution attack and joint time-frequency scattering
audits check for metric-specific overfitting. The current scattering audit runs
at 16 kHz and **cannot assess sizzle above 8 kHz**. Also check repeated strikes,
velocity changes, finite output and the actual browser limiter.

Initial screens use seeds 1675/1982. Final observation fitting uses **1675 only**:
averaging random realizations into the one recorded target compromised its
actual standard-strike balance. Seed 1982 remains a validation render. Texture
checks also use 6841/7913; a separate
broader audit uses 9047/9851/11029. These are random realizations of one gesture,
not additional independent reference recordings. Only listening can establish
whether the remaining mismatch is acceptable.

## Tools and artifacts

- `tools/refine_gong_layered_timing.py`: transport screens, smooth observation
  fitting, optional reference weighting and verified cache reuse.
- `tools/refine_gong_upper_movement.py`: guarded texture coordinate comparisons.
- `tools/refine_gong_layered_energy.py`: energy delivery versus observation
  headroom check, without changing output gain.
- `tools/plot_layered_fit.py`: fixed-level band envelopes and rise landmarks.
- Existing `audit_gong_sizzle.py`, `review_metal_refit.py`, STFT difference and
  limiter tools supply independent checks.

Artifacts are under `build/gong-layered-*`; these are developer diagnostics.
The audition destination remains the **main workbench's Gong reference target**.

## Published result — 10 September 2026

`workbench/web/gong_calibration.fit.json` now contains **Gong — pitched body and
metallic sizzle**, ID `5117b722-dc9e-413a-93a6-f316a6f8dd07`. The selected checkpoint
is `build/gong-layered-excitation/trial-0`. All 184 saved parameters reproduce the
validated render. The main workbench's target-to-editable-patch-to-saved-fit test
passes in an isolated browser tab. Existing user tabs/snapshots are not modified.

Material changes from the incoming preset:

- Diffusion strength 6 → 4; concentration dependence .2 → .05. A lower exponent
  allows weaker upper spectral energy to continue spreading instead of requiring
  excessive observation gain. Increasing speed alone was not the answer.
- Excitation centre 2118 → 1500 Hz, with the same dark excitation tilt. Body
  excitation rounds from 2.990 to 3; it is not an added output gain.
- Packet spread 2.717 → 1.8 ERB. More concentrated oscillator groups produce
  stronger local upper fluctuations without increasing random phase blur.
- Beat-rate tilt .25 → .5, retaining beat depth .15 and the low-frequency rate.
- Smooth broad modal-level rebalancing, with every centre frequency unchanged.
  The two-point damping curve stays approximately 10/2.14 seconds: shorter/longer
  endpoint screens did not justify a larger decay redesign or extra knots.

Phase blur (.035 ERB, tilt −.5), bounded movement (1 rad, 100 changes/s, sharing
.25), density, output EQ, Body observation, Model level and reference gain stay
unchanged. The engine and UI control set are unchanged in this tuning pass.

For the saved standard strike:

| Diagnostic | Incoming | Published | Reference |
|---|---:|---:|---:|
| 3–15 kHz centroid, .5–1.5 s | 4829 Hz | 5933 Hz | 5995 Hz |
| 9–14 kHz half-peak rise | .863 s | .678 s | .688 s |
| 5–9 kHz peak power | −33.48 dB | −33.35 dB | −33.43 dB |
| 12 kHz detrended envelope power, 8–32 Hz | .00253 | .00763 | .00839 |
| Same, 32–128 Hz | .01260 | .03030 | .03869 |
| Fine ridge contrast, 7–14 kHz | 4.62 dB | 4.63 dB | 5.10 dB |

The same reference-weighted band objective drops from 4.581 to 1.538; this is not
an auditory acceptance threshold. Same-scale STFT/difference and absolute band
plots were inspected. **Remaining mismatch:** the 2.5–5 kHz half-peak rise is
still .234 s versus .484 s, the highest peak remains about 1.9 dB quiet, and the
lowest band's decay/pitch detail is not an exact match. The upper spectral centre
also varies across random realizations (about 5.1–5.9 kHz in the texture audit).

The final independent audit improves all five reported measures on all four
tested seeds (1675/9047/9851/11029): Mel, attack Mel, attack ridge MR-STFT,
regional modulation texture, and JTFS below 8 kHz. On the standard seed, Mel
falls 1.736 → 1.185, texture .443 → .360 and JTFS .166 → .074. These checks were
not the final observation-fit objective and still do not certify a listening
match. Full results are in the selected checkpoint's `refit-audit.json`.

Velocity tests at .3/.5/.76/1 produce increasing energy and increasing upper-to-
body energy ratio, with no new velocity remapping. At 48 kHz, the tested standard
single hit and four half-second-spaced strikes peak at −6.57 dBFS and require no
limiting even at Master 0 dB. Eight rapid .125-second-spaced hits exceed 0 dBFS
before browser gain/limiting at 44.1 kHz; the limiter remains necessary protection,
not part of the fitted timbre. There is no claim that arbitrary sequences cannot
reach it.
