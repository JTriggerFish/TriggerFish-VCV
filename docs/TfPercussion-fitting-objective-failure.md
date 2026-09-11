# Why the gong fitter removed the low body

This is a diagnosis of a rejected fit, not another calibration or a new DSP
model. The current workbench sound is unchanged. The temporary below-300 Hz
observation lock has been removed: it concealed the cause rather than fixing it.

## The actual failure

Three problems compounded one another.

1. **The fitting coordinates coupled the wrong changes.** A smooth *correction*
   curve had anchors at 120 and 360 Hz, but no independent direction at 240 Hz.
   Raising 360 therefore raised 240. Suppressing 120 was how the optimizer
   counteracted that unwanted increase while retaining the stronger 360 mode.
2. **The objective allowed one pitch to substitute for another.** It sums power
   across 80–250 Hz before comparing envelopes. It cannot directly distinguish
   the 120 and 240 Hz components inside that band. Increasing FFT size alone
   cannot restore information discarded by this sum.
3. **The task's target changed.** The preceding `LayeredBloomLoss` explicitly
   combined the user's stronger low body with the real reference's upper body.
   `LayeredBandLoss` instead targets the recording everywhere, while retaining
   the user's harmonic centres. That is a different task. Matching the recording
   would legitimately reduce some bass, but did not justify almost removing it.

The renderer/cache checks proved that the saved sound was the one being scored.
They did **not** prove that the score or permitted parameter directions described
the intended sound. Aggregate Mel, MR-STFT and JTFS improvements also failed to
establish that each important feature survived. Publication was the additional
review failure: those improvements were accepted despite the low-body regression.

## The parameter coupling is exact, not speculative

For log-frequency smoothstep interpolation between 120 and 360 Hz, the weight
of the 360 Hz anchor at 240 Hz is 0.691906. Therefore:

$$
\Delta L_{240}=0.308094\,\Delta L_{120}+0.691906\,\Delta L_{360}.
$$

The actual accumulated changes were:

| Bar | Level change |
|---|---:|
| 120 Hz | −39.8704 dB |
| 240 Hz | −5.3705 dB |
| 360 Hz | +9.9917 dB |

Substitution reproduces the 240 Hz change to floating-point precision. This
wasn't a request from the engine to eliminate bass. It was a restriction imposed
by our optimizer's coordinates. Smooth *adjustments* are not necessarily a
sensible smooth *result*, especially when repeatedly applied to uneven bars.

## Controlled measurements

All measurements use the actual workbench WASM, the same saved gesture and
reference gain, and six seconds of mono audio. No limiter, normalization or
frequency lock is introduced. Reference and saved renders are checked against
the current renderer before the experiment.

| Test | Weighted layered score (lower is better) |
|---|---:|
| Incoming sound | 4.5806 |
| Only suppress its 120 Hz bar as the rejected fit did | 4.0932 |
| Rejected complete fit | 1.5379 |
| Re-optimize only its first two amplitudes independently, same old loss | 1.4394 |

For the single-bar suppression, 80–180 Hz power over 40–200 ms falls from
−17.56 to −44.79 dB; the reference is −27.24 dB. So the original was too strong
for an exact reference match, but the suppressed version is **17.55 dB too weak**.
The reference's early spectrum contains a peak near 122 Hz. Its strongest early
components are near 346 and 375 Hz; these are spectral observations, not proof
of a single harmonic fundamental or exact physical mode identification.

Allowing the first two amplitudes to vary independently gives approximately
−26.44 and −35.88 dB, from both tested starting points. The first was −54.73 dB
in the rejected fit. **Even the old loss restores that mode when it is given a
less restrictive search space.** This is why a loss-only fix is insufficient.

A second experiment uses equal-log-frequency regional spectral errors over the
whole audible range, rather than the broad-band sum. With the same two free
amplitudes it gives approximately −29.42 and −29.86 dB, again from both starts.
This is an experimental comparison, not a perceptual certificate or a published
fit. No new rule tells either optimizer to retain a particular note.

## Decay is a separate issue, not an excuse for the missing note

The broad lowest-band envelope-shape error worsened from 1.89 to 2.81 dB even
as that band's total objective contribution improved. Level and shape must be
reported separately. Lowering an observation bar does not repair its decay.

One-factor tests shortening the low T60 endpoint from 9.89 to 7, 5 or 3 seconds
worsen the full objective at fixed other parameters. That is consistent with
damping also reducing energy available for later bloom; it is **not** evidence
that the original damping is correct, or that T60 should be replaced by gain.
Joint dynamics/observation tests are needed to distinguish these trade-offs.

## What changes in the methodology

- **State the target before fitting:** matching the recording, or refining a
  user's designed sound against selected reference attributes. Record that
  distinction explicitly; do not silently replace one with the other.
- **Keep frequency-family parameters distinct from prominence fitting.** Start
  from the regular/stretched series. Test root/stretch when pitch differs;
  don't force a fixed incorrect series to match through amplitude deletion.
- **Fit with directions the actual controls can express.** Use actual modal
  amplitudes for observation optimization, with a *soft* complexity penalty
  where needed. Do not hard-interpolate corrections through a few anchors that
  prohibit a required peak or trough. This need not introduce new UI controls,
  individual frequency placement or per-mode damping parameters.
- **Retain frequency information and separate time regions.** The broad-band
  objective is useful for energy delivery, not sufficient for pitch/ridge
  preservation. Inspect regional spectral deficits/excesses alongside band
  bias, envelope shape and texture. Strong high-frequency ridge identification
  is not assumed; distribution/spacing remains preferable there.
- **Validate optimizer behaviour on counterexamples.** Tests now demonstrate
  the same missing-tone failure at three frequency ranges. The diagnostic has
  a reference-fixed floor/mask and detects the missing component in all three.
  The exact-render two-amplitude experiment tests search-space bias separately
  from loss bias. These are stronger evidence than a reference-versus-itself
  zero-loss test, but are not a complete perceptual validation suite.

The new regional diagnostic is intentionally **not** silently added with an
arbitrary weight to every fitter. Its use as an optimization residual remains
experimental. Next is a regularized actual-amplitude fit, jointly checked with
the shared dynamics, plus known-parameter recovery and other-instrument checks
before promoting a replacement objective. The restored workbench preset has
not been overwritten by these experiments.

## Tools and reproduction

`tools/audit_gong_fit_objective.py` produces additive loss attribution, regional
spectra and one-factor actual-render tests. `tools/profile_observation_objectives.py`
compares the two objectives in the same validated observation subspace.

```powershell
$env:EMSDK_NODE='D:/dev/emsdk/node/24.19.0_64bit/node.exe'
$env:OPENBLAS_NUM_THREADS='1'
.venv/Scripts/python.exe tools/audit_gong_fit_objective.py --incoming build/gong-layered-timing/baseline --rejected build/gong-layered-excitation/trial-0 --output build/gong-objective-audit
.venv/Scripts/python.exe tools/profile_observation_objectives.py --target gong-standard --source build/gong-layered-excitation/trial-0 --output build/gong-objective-audit --keys resolved_level_0 resolved_level_1
```

Artifacts are local analysis files, not a new audition site or shipped preset.
The JSON includes the source paths, renderer/reference metadata and exact-render
verification. No browser audio is played.

## Literature: relevant, but not a substitute for these tests

[Torres, Peeters and Richard, ICASSP 2024](https://arxiv.org/html/2312.14507v2)
discuss why pointwise spectral losses can poorly guide oscillator-frequency
estimation and investigate spectral optimal transport for harmonic synthesis.
That supports investigating a frequency-displacement comparison for root/stretch;
it does not validate a gong loss, nor explain our amplitude-coordinate defect.
Transport compares normalized spectral distributions, so a separate absolute
level comparison remains necessary. We have not added it to the production fit.
