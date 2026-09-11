# Gong refinement — 9 September 2026

## Scope and starting point

This pass starts from `Gong — irregular gentle ringing` at revision `50d53e6`.
The target remains the standard Dresden gong reference and its saved mallet
gesture. Render six seconds through the **actual workbench WASM**, with fixed
source, alignment and output gains. No sample normalization or engine changes
are part of this fitting pass.

The immediate problem is excessively regular low/mid ringing, not merely the
amount of noise. A lower aggregate error is not sufficient evidence of a fix.

## Diagnosis

Fixed-scale STFTs and one-second spectral slices show that the old 240, 360 and
480 Hz handles dominate the lower spectrum. The reference has a fundamental
near 120 Hz, a cluster around 285–300 Hz, strong components around 345–375 Hz,
and a cluster around 535–545 Hz. Their relative prominence changes with time.
These are **clusters**, not proof of stationary, independently identifiable
physical modes. Zero padding the spectrum does not improve its true resolution.

The previous modulation summary covered only 0.5–12 Hz. It missed prominent
20–60 Hz envelope modulation in several model bands. The analysis now reports
12–80 Hz flutter separately, and the motion plot shows the wider range. A
synthetic 30 Hz amplitude-modulation regression test exercises this blind spot.
The reference also contains fast modulation: minimizing all motion is wrong.

## What is fitted, and how

1. **Ablations:** screen packet layout, pair depth, spread, pitch wander and
   frequency-dependent phase blur. Preserve every other parameter. More wander
   alone produces little improvement. Aggressively cleaning the first packets
   loses too much low/mid energy.
2. **Low-core trial:** compare two explicit four-handle arrangements against the
   untouched baseline, with 0/4/8 dB attenuation of handles 2–4. The remaining
   28 centres and levels stay unchanged. This is a small reference-guided
   exception to series placement, not an unconstrained ridge-fitting search.
3. **Shared dynamics:** bounded Powell search over bloom rate, energy
   nonlinearity, excitation tilt, packet spread, pitch-wander depth, and the
   two active T60 endpoints. No extra decay knots or per-mode decay factors.
   Local finite differences record influence; the budget is 200 parameter
   evaluations, each rendered at two fixed seeds (1675 and 1982).
4. **Validation:** inspect fixed-scale STFT/difference and low-band motion plots;
   compare additional seeds, attack, bloom, decay and repeated strikes. Verify
   that the saved JSON reproduces the exact candidate WAV before publication.

The shared-dynamics scalar ranking is:

$$
L=L_{\mathrm{Mel},60\,\mathrm{dB}}+0.04L_{\mathrm{bloom}}+0.025L_{\mathrm{decay}}.
$$

The first term is Auraloss-based Mel comparison with a **reference-fixed**
floor. The second measures absolute band envelopes and onset-to-bloom contrast;
the third measures relative band decay. The explicit weights are pragmatic
search trade-offs, not a validated psychoacoustic scale. Keep the constituent
scores and visual checks separate. No scalar score confers listening approval.

## Reproduction and artifacts

Use the existing development Python environment and `EMSDK_NODE`; build the
workbench through `dev.ps1` when required. The server is not involved in fitting.

- `tools/refine_gong_slow_beating.py`: named screens and four-seed validation.
- `tools/refine_gong_shared_envelope.py`: shared refinement (originally seven
  coordinates; now eight after separating concentration and energy sensitivity).
- `tools/plot_spectral_difference.py`: reference-fixed STFT and differences.
- `tools/plot_modal_motion.py`: band motion and modulation spectra.
- `tools/review_metal_refit.py`: additional-seed and repeated-hit audit.

Local experiment directories start with `build/gong-` and end with `50d53e6`.
Each screen archives its starting parameters; checkpoints include renderer
hashes, reference metadata, complete parameters and exact audio. These private
audio artifacts are not committed. Audition remains in the main workbench.

## Published candidate and limitations

The main workbench's standard Gong now loads **Gong — retuned low core**.
The selected checkpoint is
`build/gong-post-envelope-ring-50d53e6/candidate`. The exact original fit is
archived alongside it as `original-workbench.fit.json`; the original remains in
Git history as well. Browser target selection and save/reload reproduced all
180 parameters and the reference metadata. The published patch also reproduces
the verified WASM candidate sample-for-sample.

Changes relative to the original:

- Lowest centres: 120/240/360/480 → 120/285/350/540 Hz. Reduce the latter three
  observation bars by 4 dB. Keep the upper 28 handles unchanged.
- Paired central ringing, depth 0.15, nominal split 2 Hz; low four packet
  noisiness multipliers 0.65. This replaces doublets within the satellite cloud,
  so the old and new split numbers do **not** describe the dominant total beat
  rate. Independent pitch wander is 1.5 Hz at one target per second.
- Phase blur 0.035 ERB with a −0.5 frequency tilt; spread 2.72. No new controls.
- Shared bloom rate 2.84, excitation tilt −47.90 dB/octave; T60 endpoints
  9.89/2.14 seconds. No extra knots, per-mode damping, gain or EQ adjustments.

Four-seed audit (two fitting seeds and two additional seeds): attack-ridge
MR-STFT error improves about 9%, and band-decay shape error about 28%.
Reference-fixed 60 dB-floor Mel error **worsens about 4%**. The candidate is
selected as a pitch/decay and ringing trade-off, not as a universal loss winner.
The intermediate shared-envelope candidate has better Mel error but less
improvement in slow-motion depth and decay.

Fixed-scale STFT and motion plots were inspected. The low/mid stripes are less
dominated by the old harmonic centres, but the 90–180 Hz modulation is still
more concentrated near 17–20 Hz than in the reference, and irregular slow
motion remains deficient in other bands. The upper bloom also remains an
approximation. **This is not a completed perceptual match or listening approval.**
Audio was rendered and measured; final auditory acceptance remains with the
user. Quarter-note and rapid-hard repeat renders were checked for finite output;
the browser's existing limiter remains the listening safety stage.

## Subsequent control-surface split

Concentration dependence and total-energy sensitivity are now separate visible
parameters. The preset gains the explicit energy exponent
`bloom_energy_sensitivity = 2 * bloom_energy_acceleration`, preserving the old
law and bringing the saved surface to 181 parameters. This is a parameter
conversion, not another gong refit; the six-second migration render was checked
sample-for-sample. The current refinement tool can search the two exponents
independently. Its extra degree of freedom must not be attributed retrospectively
to the seven-coordinate experiment and scores above.

## Pitched attack refinement — 10 September

The standard Gong now loads **Gong — clearer pitched attack**. This is a
partial pitch improvement, not a completed match of the first half-second.
The starting revision is `e2a7eb6`; the selected local checkpoint is
`build/gong-pitched-balanced-e2a7eb6`. No DSP or control-surface changes.

### Diagnosis and scope

Fixed-level spectra over 40–200 ms and 200–500 ms show that simply boosting
the fundamental is the wrong correction: the roughly 121 Hz component already
has comparable strength. The reference has a prominent 344/374 Hz pair,
whereas the previous model overemphasizes 285 Hz and has no corresponding
374 Hz centre. These are measured low-core hypotheses, not an attempt to
identify and individually fit the dense upper spectrum.

Retune the first five handles to 121/288/344/374/537 Hz. Their observation
level changes are 0/−7/+2/+2/−4 dB; the fifth local noisiness becomes 0.65,
matching the other four. Leave the upper 27 handles unchanged. Keep both T60
endpoints, gains, EQ, event velocity and implement unchanged. There is no
per-mode damping adjustment, extra decay knot or audio normalization.

### Search and checks

1. `refine_gong_pitched_attack.py` screens explicit texture, four-core,
   five-core and concentration/rate hypotheses using the actual WASM renderer.
   It archives the original full parameter surface before searching.
2. `polish_gong_pitched_attack.py` refines only rate, concentration exponent,
   energy exponent and excitation tilt with bounded Powell search and recorded
   finite-difference influence checks. The selected 70-evaluation run starts
   from five-core trial `pair-local-0.65-balance-2-rate-1`.
3. The primary objective is existing Auraloss attack MR-STFT on the first
   0.5 s after identical causal sixth-order 1.5 kHz low-pass filtering of both
   signals. Its FFT sizes are 4096/8192/16384. Whole-six-second reference-fixed
   Mel, band-envelope and bloom-rise errors remain separate diagnostics.
4. Ranking divides attack error by the original two-seed mean and adds hinge
   penalties of weight 10 for Mel/envelope/rise exceeding 1.05/1.10/1.10 times
   their original means. These are soft search penalties, not hard guarantees
   for every seed, and not perceptual acceptance thresholds.
5. Fit seeds are 1675/1982; additional checks use 2586/3276. Saved parameters
   and snapshot audio are checked against the exact candidate render.
   `plot_gong_attack.py` displays absolute-level low spectra and band envelopes;
   `plot_spectral_difference.py` supplies reference-fixed STFT differences.
   Both plots were inspected, including the remaining mismatches.

Final shared parameters are rate 3.0380, concentration exponent 0.04693,
energy exponent 0.11856 and excitation tilt −47.9417 dB/octave. Exact values
are in the preset and checkpoint. On the four-seed mean, low-attack error
improves about 8.5% and Mel about 3%; **band-envelope error worsens about 9%**.
More aggressive low-packet cleaning scored better on the nominal attack but
regressed other seeds and the whole envelope; it was rejected.

### Remaining limitation

The reference's 300–450 Hz body holds up through roughly 0.4 s; ours falls
away during that interval. The reference's 3–12 kHz wash has a deeper early
dip and blooms later; ours fills that region too early. Changing rate and
concentration alone in the tested ranges did not recover both behaviours.
The next fitting target is this **sustained pitched body versus developing
wash**, not another blanket fundamental boost or a claim that the attack
loss alone establishes success. This pass rendered and inspected measurements;
it does not claim auditory approval.

Publication checks: the served target loads/saves all 181 parameters and the
correct reference; its instrument, gesture and reference objects exactly match
the verified checkpoint. Quarter-note and rapid-hard repeat renders are finite,
with peaks −5.47 and −2.14 dBFS respectively, before the browser safety limiter.
The five attack-study tests and the `dev.ps1 build-workbench` build pass.

## User harmonic body plus metallic bloom — 10 September

The main target now loads **Gong — tuned body and bloom**, starting from user
snapshot `477215bd-6dc9-40ef-897b-15981756f73c` (the limiter-audit gong).
This is deliberately a **sound-design compromise**, not an exact reference fit:
keep the user's stronger 120/240 Hz pitched body and add the reference-like
upper swell. Every one of the user's 32 centre frequencies remains unchanged.
No new controls, DSP changes, direct noise layer or artificial onset delay.

### Method and selected values

- Measure five STFT bands (80–300, 300–900, 900–3000, 3000–7000 and
  7000–14000 Hz) in twelve time regions through six seconds. Use the same
  4096-point window/10 ms hop for reference and model, with a fixed −100 dB
  analysis floor and no audio normalization.
- Explicit hybrid target: the lowest band's envelope is the user's minus
  16 dB, leaving it approximately 7 dB stronger than the reference early on.
  Other bands retain the reference's saved level. The scalar objective is
  twice the lowest-band RMS dB error, plus middle-band, high-band and high
  onset-to-bloom contrast errors. This low-band target is a design choice,
  not something inferred from the reference or automatically applied by the UI.
- Screen 48 rate/concentration/energy-sensitivity combinations. A subsequent
  broad-tilt joint search was stopped after about 50 evaluations because its
  observation shape hollowed the middle. Retain the screen, not that abandoned
  optimizer. The scripts now separate dynamics screening and observation fitting.
- Fit only four observation-curve coordinates, at 240/600/3000/12000 Hz, with
  smooth log-frequency interpolation. Use the existing validated exact-render
  observation basis and STFT cross-power cache. No independent ridge placement
  or 32-dimensional level optimization. Verify the predicted band envelopes
  against complete WASM renders before accepting the result.
- Compare two dynamics starts and screen the existing high T60 endpoint.
  Longer damping settings did not improve the composite result. Keep the curve
  at approximately 9.89/2.14 seconds, still only two active endpoints; no
  per-mode damping changes.
- Selected diffusion strength **6**, concentration dependence **0.2**, energy
  sensitivity **0.3**. Excitation, packet texture, phase blur and EQ stay as
  in the user's snapshot. Raising Body observation to **4** gives upper bars
  more output range; compensate by lowering low/mid observation bars. This is
  the visible existing gain, not an extra gain or energy injection.

The final four bar adjustments relative to the user's levels are
−20.857/−6.907/+9.793/+47.9999 dB, interpolated and clipped to the existing
bar range. The exact resulting bars are stored in JSON; no hidden interpolation
curve exists in the instrument. With the observation-gain increase included,
the lowest pair is about 13.5 dB quieter than the user's very loud snapshot,
but remains clearly stronger than the published gong before this pass.
Model level stays 0 dB and Master is not changed.

### Verification and limitations

Fit seeds: 1675/1982; additional seeds: 2586/3276. Four-seed mean high-band
envelope error drops from 18.17 to 5.49 dB, and high-band rise error from 29.43
to 6.78 dB. The selected hybrid low-envelope error is 1.99 dB. These are
diagnostics against the stated hybrid target, not a perceptual approval score.
STFT/difference and band-envelope plots were inspected. The 3–12 kHz band
now rises about 25 dB; the user's original upper band barely rose at all.
The candidate still blooms somewhat early, lacks some very-high-frequency
energy, and intentionally retains excess harmonic bass versus the reference.

The affine cache's largest tested band-envelope discrepancy was 0.0008 dB.
Saved snapshot audio reproduces the exact WASM render. At 48 kHz the saved
single hit peaks at −4.97 dBFS and four half-second-spaced strikes at −4.12 dBFS;
neither invokes the actual browser limiter even at Master 0 dB. This is not
a guarantee for arbitrary stronger gestures or longer sequences.

Artifacts are under `build/gong-layered-*`; the selected checkpoint is
`build/gong-layered-final`. Reusable steps are
`refine_gong_layered_bloom.py`, `polish_gong_layered_observation.py` and
`refine_gong_layered_decay.py`; limiter checks use `audit_fit_limiter.mjs`.
The prior user snapshot remains untouched. Audition is in the main workbench,
not a separate report page; auditory acceptance still requires listening.

The later [layered gong refinement](TfPercussion-gong-layered-refinement.md)
supersedes this hybrid-target tuning as the published starting point. It fits
the real reference throughout, separates rise/level/texture checks, and records
the remaining timing mismatch rather than treating a lower scalar score as a
complete match.
