# Crash: perceptual losses and recoverability

## Question and scope

Is the poor crash match caused by the search, its distance measure, or the
instrument's capacity? These are different questions. A failed local fit cannot
establish that a different DSP architecture is necessary.

The 2026-09-08 experiment keeps the C++ engine, reference, level and performance
event fixed: Private crash A, edge, velocity 72, repeat 1, onset 50 ms, six-second
fit window. No playback normalization, new DSP coefficients, per-mode damping or
extra T60 knots. The successful kick preset and other instruments are untouched.
Private artifacts are under `build/crash-perceptual-v5/`.

## Which perceptual representations were tested?

- **auraloss mel MR-STFT:** the library loss at full 48 kHz; windows
  512/2048/8192, 64 mel bands, spectral convergence plus log-magnitude error.
  Scale invariance is disabled. An A-weighted variant is audited separately.
- **JTFS:** the existing WaveSpin transform, evaluated on MLBox. Both signals
  are resampled to 16 kHz, retaining the time axis with 16 ms averaging. This
  represents spectrotemporal modulation; it is not a learned listener score.
  Its 8 kHz bandwidth is a significant limitation for cymbals.
- **Control:** the existing ERB-weighted metallic spectral/envelope objective.

See the [auraloss implementation](https://github.com/csteinmetz1/auraloss) and
[JTFS paper](https://arxiv.org/abs/2204.08269). The precise library versions,
transform settings and reference floor are recorded in `audit.json`.

Identity returns zero. Increasing pitch-and-speed changes, faster decay and
added envelope-following wash receive increasing penalties in all three audited
perceptual variants. Pitch resampling also changes duration; it is not an
isolated pitch test. These checks establish sensitivity, not human agreement.

The previous crash trials were cross-scored before fitting. Plain mel preferred
the published baseline to every trial. JTFS gave one trial only a small advantage.
This confirms that the earlier custom-score gains were not a general improvement.

## Known-answer test on the full C++ crash

Render the published crash with a fixed seed. Change four controls in the start:
tune +6%, low T60 -25%, high T60 +30%, phase bandwidth -50%. All remaining
controls and the target's noise realization are fixed. The target is generated
by the actual full model, so an exact solution demonstrably exists.

Fit with mel MR-STFT, bounded L-BFGS-B and central differences of 0.001 of each
normalized control range. Frequency/time coordinates use their declared UI
scales. After 200 parameter evaluations:

| Control | Remaining relative error |
| --- | ---: |
| Tune | -0.0006% |
| Low T60 | -9.76% |
| High T60 | +33.89% |
| Phase bandwidth | +1.84% |

Holding the recovered tune and phase bandwidth, a separate 120-evaluation
two-endpoint damping stage reduces the T60 errors to **+0.11% and -0.61%**.
The loss falls 0.933 → 0.287 → 0.132. The remaining phase-bandwidth error means
this is not waveform-exact recovery. No true target parameter was substituted
into the recovered patch between stages.

This demonstrates a search-conditioning/staging problem. It does **not** prove
that the real crash recording lies within the model's capabilities.

## Matched real-reference searches

Each objective starts at the published patch and fits the same six controls:
tune [0.8, 1.2], low T60 [5, 28] s, high T60 [0.2, 5] s, phase bandwidth
[0, 1.5] ERB, cascade [0, 1.5] oct/s and excitation tilt [-12, 8] dB/oct.
Other controls, including modal observation levels, stay fixed for this test.
Each gets 180 unique parameter evaluations and averages two separate seed scores.
All exhaust that budget; none is declared converged.

| Objective | Own objective before → after |
| --- | ---: |
| Custom, squared residual norm | 87.155 → 86.281 |
| Mel MR-STFT | 1.3692 → 1.3591 |
| JTFS | 0.12875 → 0.11445 |

Numbers in different rows have different units and cannot be compared directly.
JTFS's improvement still leaves excessive late wash in the inspected difference
plot, and worsens independent envelope and ridge measurements. Substituting a
perceptual representation therefore does not by itself solve the crash.

## Changes to the fitting procedure

The actual-render affine observation fitter can now differentiate the library
mel loss directly. It validates the library score and its amplitude derivative,
and still verifies mixed-gain C++ renders. No second synthesizer is trained.

An optional `PerceptualEnvelopeGuard` constrains **each seed independently**:

- absolute band-level RMS error and relative decay-shape RMS error in each
  usable reference band, from 0.2 s onward;
- raw-power errors in 0–1, 1–3, 3–10, 10–30 and 30–100 ms.

Each error must not exceed its explicit comparator's error plus 0.05 dB. This
is a non-regression allowance, not an accuracy target or playback correction.
The guard uses the canonical 4096/512 STFT and the reference-only contamination
mask from `BandDecayShapeLoss`. Gradients are checked numerically. SLSQP retains
its best feasible evaluated candidate if its last step is infeasible. A prepared
body with different frequencies or damping may have no feasible amplitude-only
solution against the old body; that is not proof that the whole model lacks one.

Additional trials retain broad packets while using the six spare handles for
cleaner measured resonances. Those use the already exposed local turbulence
scales, not extra decay parameters. A single-realization fit is also compared
with the two-seed fitting procedure; other seeds remain validation data. These
are explicit ablations, not reasons to waive reference/decay/texture checks.

The early guarded-run history records square roots of native mel scores, as the
old scalar-residual adapter did. The observation fitter now records native mel
units consistently; rerendered native scores, not differently scaled history
numbers, must be used for cross-run comparisons.

## Final audit and publication decision

The guarded amplitude-only searches did not find a feasible improvement within
35 SLSQP iterations, including a run starting from the feasible published patch.
The 0.05 dB allowance is deliberately tight. This is a failed bounded experiment,
not evidence that these constraints are always appropriate or that a solution
cannot exist. The guard remains optional, not the default search policy.

The mixed-packet experiment retains the 26 existing active handles and adds six
cleaner centres at 436.52, 1520.51, 1573.24, 2088.87, 2885.74 and 3723.63 Hz,
each with local turbulence scale 0.25. The two existing clean low handles move
to 407.23 and 125.98 Hz. No per-mode damping is fitted. After the two-seed trial,
a single-seed observation/damping/observation sequence produces the following
standard-strike measurements:

| Measurement | Published patch | Mixed single-strike trial |
| --- | ---: | ---: |
| Composite engineering score, mixed units | 9.472 | 8.959 |
| Envelope component | 3.542 | 3.405 |
| Linear-spectrum component | 15.169 | 14.224 |
| Contrast component | 7.371 | 7.077 |
| Attack component | 6.022 | 6.051 |

This is a useful local improvement, not an accepted calibration. In particular,
the first 1 ms RMS excess increases from 1.31 to 3.44 dB. Three fresh seed
composite scores change 9.49 → 9.38, 9.75 → 9.89 and 10.18 → 10.21. Contrast
improves in all three, but envelope error worsens in two. Thus the improvement
does not consistently generalize to other realizations.

The inspected difference plots still show incorrect resonance structure and
decay. The independent ten-second audit finds the 1–4 kHz band **3.93 dB too
quiet at 0.12–0.5 s and 4.95 dB too quiet at 0.5–1.5 s**; the late 4–16 kHz
band is 2.15 dB too loud at 3–6 s. A single static gain cannot repair that
time evolution. Snapshot reload reproduces the render exactly. Repeated-hit
renders were also generated, but have no listening approval; the full-strength
sequence reaches a raw peak of 1.76, before browser output protection.

Artifacts: `mixed-single/search.json`, `fresh-seed-audit.json`,
`full-review.json` and `difference.png` under the experiment directory. The
published workbench preset is **unchanged**. No other instrument is altered.

The next targeted experiment should isolate the 1–4 kHz build-up and decay,
varying existing excitation/energy-transfer controls with observation levels
fixed before refitting those levels. Keep two damping endpoints and independent
attack, late-tail and fresh-seed checks. If that restricted experiment cannot
produce the required time evolution, inspect the relevant DSP mechanism before
adding controls. The successful kick is a useful reference for the fitting
workflow, not proof that the same objective/search is sufficient for a crash.

The [following dynamics and modal-placement pass](TfPercussion-crash-dynamics-review.md)
tests these hypotheses and adds a reusable exact quadratic observation-energy
fitter. Its results and publication status are recorded separately.
