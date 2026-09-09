# Bloom and damping: what can be made independent?

September 2026 investigation of the user's **Candidate crash**, snapshot
`e7f00f01-d8ab-4602-a181-8be49ce0f5d4`. This document separates measured
results from the implemented **Hold decay** design-time UI assistance. No
audio-rate DSP compensation, hidden envelope or gain matching has been added.

## Why the controls interact

The current spectral diffusion step conserves stored modal energy. Each mode
also has ordinary damping from the shared T60 curve. In energy units:

$$
\dot E_i=F_i(E)-\lambda_iE_i,\qquad
\sum_iF_i=0,\qquad \lambda_i=\frac{6\ln 10}{T_{60,i}}.
$$

Consequently $\dot E_{\mathrm{total}}=-\sum_i\lambda_iE_i$: moving energy
into faster-damped modes accelerates total loss even though the transport
step itself removes nothing. A band's decay also includes its inflow/outflow.
Painted prominence and observation EQ further change what is heard as energy
moves. Stored energy, output loudness and per-band decay are not interchangeable.

With uniform damping, conservative redistribution cannot change total stored
energy decay. With frequency-dependent damping, complete independence is not
generally possible. Cancelling arbitrary band outflow would require negative
damping once outflow exceeds the requested decay. We should not do that.

This interaction is consistent with the phenomenological plate literature:
[Humbert et al., equations 3–4 and section 4](https://arxiv.org/html/1709.09884)
explicitly couple spectral transfer and frequency-dependent dissipation.
The paper is not evidence for a particular compensation algorithm.

## Why this snapshot is especially difficult

Diffusion nonlinearity is 0.001. Our conductivity scales approximately as
$\rho^{2a}$. If energy density falls by 40 dB without changing spectral shape,
conductivity at $a=0.001$ retains $10^{-8a}\simeq0.982$ of its value. It is
almost linear diffusion throughout the audible tail, not just near the hit.
Local gradients still evolve; this calculation isolates the energy scaling.

The candidate's low T60 endpoint was already 27.87 seconds. Raising that alone
cannot prevent continuing transport out of the low body. The plots show low
and mid tails falling too quickly, while the upper spectrum is sustained too
long early in the decay. One overall gain or decay multiplier cannot fix both.

The user snapshot has strength 0.960 / location 0.512, whereas the selected
reference cell is strength 0.567 / edge. Both auditions are archived, but
reference fitting uses the latter gesture. These are not interchangeable tests.

## Existing-control fit

`tools/refine_user_crash_decay.py` first fits only two T60 endpoints, then
jointly fits those endpoints, diffusion rate and nonlinearity. Frequencies,
painted prominence, packet texture, excitation, observation, gain and gesture
stay fixed. No per-mode damping or extra curve knots are added.

The objective separately measures reference-relative decay shape in 0.5–1,
1–2.5 and 2.5–5.5 seconds. Equal-band regional weighting prevents the longest
region from dominating. A soft penalty above 1.5 dB per time-frequency cell
protects the user's first 450 ms at the reference gesture. This is not a hard
guarantee of perceptual equality. The analysis subtracts one constant per band
for decay shape only; no playback level is normalized.

Log-parameter finite differences use ±0.02, with bounded trust-region least
squares (28 function evaluations maximum; Jacobian probes are additional).
Bounds are T60 low 0.1–30 s, T60 high 0.1–10 s, rate 0.2–16, nonlinearity
0.001–0.5. Actual Wasm renders are used throughout.

| Variant | Decay shape error | Front change from user's patch |
|---|---:|---:|
| User parameters, reference gesture | 5.73 dB | 0 dB |
| Two T60 endpoints only | 5.46 dB | 0.52 dB |
| Joint existing-control fit | 3.17 dB | 0.88 dB |

Selected values: rate 6.5966, nonlinearity 0.2260, T60 endpoints 30.000 and
0.9287 seconds. A higher nominal rate offsets the stronger energy dependence
during bloom; late transfer then subsides more substantially. Across four
seeds, decay-shape errors fall from 5.73–5.82 to 3.17–3.21 dB. Reference-floor
Mel and the independent texture metric also improve on these seeds.

At the user's exact stronger gesture, the first 450 ms changes by 1.60 dB RMS.
That is an important limitation of fitting at one velocity. The original
snapshot remains archived. Plots still show residual ridge mismatch, excess
upper-band energy early in decay, and an imperfect low tail. This is an
audition candidate, not listening approval.

Artifacts: `build/crash-user-decay/before`, `damping-only`, `coupled`.
The `coupled` directory contains standard/strong-strike renders, seed and
repeated-hit audits, and fixed-scale spectrogram/decay diagnostics. The main
workbench—not a separate report—is the audition surface.

## Hold decay — implemented workbench assistance

Use **design-time compensation**, not an audio-rate AGC, output envelope,
parallel tail or another energy reservoir:

1. Capture the patch's late band envelopes at the start of a bloom edit.
2. Keep the user's requested bloom-rate change fixed. Estimate sensitivities
   of the existing damping curve and, when needed, energy dependence using
   actual renderer probes.
3. Find a small bounded correction preserving the late envelopes while
   retaining the intended early bloom change. Penalize unnecessarily large
   parameter changes. Validate with an actual render, not just a Jacobian.
4. Update the visible controls/curve and save those exact ordinary parameters.
   Anchor the whole drag to its starting patch so compensation does not drift.
   Report limits or reject a correction if it cannot preserve decay adequately.

A compact **Hold decay** lock beside Bloom is preferable to more synthesis
knobs. The default-on checkbox sits above Bloom timing. After drag release a
dedicated worker runs actual C++/Wasm renders, with elapsed time, render count
and Cancel. A new edit, reference/strike change, reset or controls rebuild
cancels the job. No user tab is reloaded automatically. Ordinary editing and
double-click reset remain available. Only visible parameter values are saved;
the checkbox and gesture anchor are UI state, not sound parameters.

An initial test in `tools/audit_bloom_decay_coupling.py` raised diffusion rate
25%. The incoming patch's late envelope changed by 2.27 dB RMS. A bounded
two-endpoint correction reduced this to 1.30 dB (predicted 1.34 dB), but the
low T60 hit 30 seconds. This demonstrates useful **partial** compensation,
not a complete solution. This was a preliminary T60-only measurement on the
incoming user patch, not the same experiment as the shipped four-control fit.

### Solver and acceptance

`decay_hold_measurement.mjs` measures 24 log-frequency bands (80 Hz to 16 kHz,
limited by sample rate), with fixed 4096-point Hann analysis. Display settings
do not affect this measurement. Five front regions cover 0–450 ms; six late
regions cover 1–6 s. Late cells below the starting peak minus 55 dB are excluded.
The tail target is the **starting patch**, not the reference sample. The front
target is the **requested edited patch**, so compensation does not simply undo
the desired bloom change. The 450–1000 ms transition and tails after six seconds
are not constrained; this is deliberately partial assistance, not an exact
output-decay law or a reference fitter.

`decay_hold_solver.mjs` / `decay_hold_fit.mjs` use three bounded Gauss–Newton
iterations with symmetric finite differences, regularized least squares and
actual-render backtracking. Coordinates are log T60 and linear nonlinearity.
Only existing active damping knots move (at most a factor of two), plus
nonlinearity within ±0.12 **unless the user directly edited it**. No new knots,
gain, mode positions, prominence, routing or strike changes are allowed.

Acceptance requires at least 10% late-error reduction, at most 1.5 dB late RMS,
1.5 dB front RMS and 4.5 dB worst front-cell deviation. A second strike seed
must also improve (or remain below 0.2 dB), with late/front RMS below 1.75 dB.
Rejected/limited corrections leave the user's raw edit intact and explain why.
Accepted residuals above 0.5 dB explicitly say **Decay partly held**. Limits are
reported; negative damping and output-level restoration are never used.

### Tests

On the published refined crash, increasing diffusion strength 25% produces
0.967 dB late-envelope change. Compensation reduces it to 0.723 dB; front
deviation from the requested edit is 0.284 dB. The second seed improves from
1.005 to 0.689 dB. The low T60 is already at its 30 s limit. This run took
25 renders / about 14 s on this machine; it is not an audio-thread operation.

`tools/audit_decay_hold.mjs` reproduces this with the actual Wasm engine.
`workbench/tests/decay_hold_tests.mjs` covers the small solver with a known test
function, unchanged gain/requested edit, bounds, rejection and analysis.
The separate silent browser integration test verifies real-worker acceptance,
exact saved parameter values, reset, cancellation and stale-result rejection.
It also checks that the refined crash is named in the preset menu. Normal Rack
builds do not depend on any of this optional tooling.

Do not silently redefine the modal T60 curve as an exact output-envelope
target: energy transport, beating and observation prevent that interpretation.
