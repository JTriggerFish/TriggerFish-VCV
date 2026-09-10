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
