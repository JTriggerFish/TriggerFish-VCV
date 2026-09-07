# Kick: controlled perceptual-loss comparison

This is an offline experiment on the existing `drum.kick.v1` voice, **not a DSP
redesign**. Contact noise is present throughout. No result is published by this
command. The main workbench remains the destination for selected listening fits.

## Question and fixed experiment

Does changing the objective improve the audible kick match, when the renderer,
reference, starting patches, editable controls and search budget stay fixed?

Use the medium-centre acoustic-kick-oak reference, velocity 64 / take 1. Trim
the recorded onset (1.315 ms), then pad to 1.2 s. Keep reference monitoring at
+2 dB and model level at -12 dB. Performance inputs are unchanged (strength,
location, hardness and implement 0.5). No waveform level matching or limiter
participates. Source hash, Wasm hash and all parameters are saved in
`build/kick-loss-comparison/experiment.json`.

The topology stays contact-direct + contact-driven modal body + independent
pitched thump. The source noise can feed the body even if direct contact is
muted. Previous fitting nearly muted that direct output; it did not remove the
noise generator. Neither of the two experimental starts mutes it.

## Objectives

| Label | Actual computation | Important limitation |
| --- | --- | --- |
| `region` | Previous band/time and regional-spectrum residual squared norm | Engineering comparison, not a perceptual model |
| `mel` | auraloss mel MR-STFT, spectral convergence + log-magnitude L1 | Frequency aggregation can conceal narrow modes |
| `mel_a` | Same library loss with its A-weighting option | Can undervalue bass relative to attack noise |
| `jtfs` | Squared distance between log-compressed joint time-frequency scattering coefficients | Not a trained listening score; finite frequency/time resolution |

For auraloss: Hann windows 512 / 2048 / 8192, hops 128 / 512 / 2048, 64 mel
bands at the original 44.1 kHz, no scale invariance, spectral convergence and
log-magnitude terms both weight 1, power clamp `1e-10`. This invokes the library
loss directly; it is not a custom STFT approximation bearing its name.

For JTFS: resample both signals to 16 kHz with the same polyphase filter. Use
the DAFx2022 authors' WaveSpin implementation, commit
`5ff1b72785703cd09d4b1bf4b5f52ffcc8a926ae`: J=10, Q=8, T=256, J_fr=3,
Q_fr=1, F=4, zero boundary padding, temporal and frequency averaging enabled.
Retain the time axis (16 ms averaging scale), not a whole-hit global average.
Compress coefficients as `log1p(abs(c)/epsilon)`, with epsilon fixed to 0.001
times the reference's largest coefficient (minimum `1e-10`). Compare by mean
squared feature distance. This tests JTFS features directly; it does **not**
reproduce the PNP paper's inverse-network training or parameter-space metric.
Above 8 kHz is excluded from this objective, but not from saved full-rate audio.

These objectives retain temporal information; attack (0–30 ms), early decay
(30–100 ms), body (100–250 ms) and tail are also inspected separately using
the fixed diagnostic tools. Native loss numbers are incomparable across columns.

## Controlled faults before fitting

Each objective is evaluated on identity and three severities of:

- Reduced low frequencies using a zero-phase frequency-domain shelf; its
  gain is exactly `1 - amount/(1+(f/180)^8)`. This is an offline test, not a
  causal synthesizer stage. Subtracting a causal low-pass would confound the
  test with phase-dependent boosts, so it is deliberately not used here.
- Added 700 Hz damped ringing, at increasing amplitude.
- Added high-passed/band-limited noise, with fixed initial amplitude and increasing decay.
- Resampling upward by 1, 2 and 3 semitones. **This changes duration too:** it
  is explicitly a pitch-and-speed test, not an isolated pitch test.

These deliberately changed signals are not estimates of the reference's latent
sources. They are useful sanity checks, not human-validated perceptual rankings.
Identity should score zero; increasing each fault should increase its penalty.
Separate synthetic tests change only a damped sinusoid's pitch while retaining
its duration, amplitude and envelope, avoiding the resampling test's confound.
The first run passes these checks for all four objectives, so none is ruled out
at this stage. This is not evidence that any objective produces a good fit.

## Matched fitting procedure

Both starts activate the same 16 existing modal slots: retain the six previous
centres and fill spare slots with the documented 200 Hz–8 kHz coverage layout.
They differ only in direct contact level and source duration/noise balance:

- `direct-long`: direct gain 0.5; existing force width and source noise settings.
- `direct-short`: direct gain 1, force width 2 ms, noise gain 0.3, base noise T60 35 ms.

Use the same bounds for all eight trials. Fit the contact controls, clean-thump
pitch/envelope/gain, modal frequencies and relative prominence, shared damping,
body tension and the existing simple output low-pass. The low-pass starts at
2.5 kHz; HP stays at 10 Hz and colour gain at zero. **No multiband EQ, per-mode
decay, routing changes, extra mode capacity or hidden DSP variables.** One modal
prominence anchor is fixed to remove a redundant common gain direction.

Optimize the native scalar objective with bounded L-BFGS-B. Frequency and time
parameters use logarithmic coordinates; others use linear coordinates. Compute
central finite differences at 0.001 of each transformed range, one-sided at a
boundary. Log the initial directional sensitivity for every control. Average
the two separate seed scores (1449 and 1460), never their audio. Each trial has
the same cap of 1000 unique parameter evaluations, including finite differences,
or earlier solver convergence. Save the best evaluated bounded candidate.

This is a bounded local experiment, not an exhaustive search or proof of model
capacity. It changes the optimizer from the earlier least-squares run; the
`region` control trial uses this same optimizer, so the **within-experiment**
loss comparison is not confounded by that change.

Finally, cross-score every result with all four objectives and run the existing
absolute band/shape checks on the primary strike and held-out seeds 1450–1452.
Inspect shared-scale plots against the reference. A lower own-objective score
does not authorize publication or establish listening quality.

## Running and optional dependencies

`perceptual-fit` is an explicit development dependency group in `pyproject.toml`.
It is not needed by normal VCV, official release, or workbench builds. The
experiment uses auraloss 0.4.0, librosa 0.11.0 and the pinned WaveSpin commit.
Local tests use torch 2.8.0. Record actual environment versions with results.

```powershell
# Add packages to the existing development venv; no native plugin build:
uv pip install --python .venv/Scripts/python.exe --group perceptual-fit
.\dev.ps1 test-perceptual-losses
.\dev.ps1 compare-kick-losses  # controlled-fault audit only, local JTFS
$env:TF_KICK_LOSS_FIT = '1'
.\dev.ps1 compare-kick-losses  # audit plus eight matched trials
```

Completed trial directories are not overwritten by another fit run. Choose a
new `TF_KICK_LOSS_OUTPUT` directory for a different experiment. Audit-only reruns
verify the saved experiment identity before refreshing diagnostics.

The current machine can set `TF_KICK_LOSS_REMOTE=1` for the explicitly prepared
MLBox worker. It scores transient PCM through a persistent SSH connection on
the GPU; **the actual C++/Wasm synthesis still runs locally**. The worker uses
an isolated dependency directory under `/tmp/tf-kick-perceptual-20260905`, not a
replacement of the existing reverb environment. It writes no audio. This remote
path is optional; local JTFS runs the same transform without SSH.

## Results: 2026-09-05

All eight trials exhausted their 1000-parameter-evaluation budgets (16,000
optimization renders in total, plus independent checks). **None is converged
by a solver stopping criterion, and none passes the publication checks.**
The long-noise start wins within each objective; that is not evidence that long
excitation is physically correct. Both starts were allowed to change duration.

For each objective's better start, the independent primary-strike results are:

| Objective | Worst evaluated band/time error | Regional spectrum P90 error | Worst band error across held-out seeds |
| --- | ---: | ---: | ---: |
| Region control | 5.4 dB | 6.3 dB | 9.8 dB |
| Mel MR-STFT | 11.8 dB | 10.9 dB | 12.8 dB |
| A-weighted mel MR-STFT | 17.4 dB | 8.8 dB | 18.4 dB |
| JTFS | 15.6 dB | 8.6 dB | 16.8 dB |

These diagnostics are not a perceptual league table. In particular, the region
control is directly optimized for the measurements in that table. The complete
cross-objective scores are retained in `trials.json`, with published-preset
scores in `baselines.json` so changed starting patches cannot manufacture an
improvement claim.

JTFS's primary score improves from the published preset's 0.01905 to 0.01544
(about 19%). Yet it still has a roughly 6 dB deficit at 1–2 kHz in the first
30 ms, and 10–16 dB excess at 2–8 kHz during 30–100 ms. Its 120–250 Hz body is
about 12 dB deficient during 100–250 ms. Those are not acceptable tradeoffs.
Plain mel improves its own score from 1.450 to 1.196 while worsening the JTFS
score from 0.01905 to 0.02374. A-weighting is not an overall improvement either.

The offline shared-scale spectrograms were inspected. The reference has a
different low-body decay structure; the trials either retain a pronounced low
ridge or shorten the body too much, and their upper attack/decay distributions
still differ. These observations do not constitute listening approval.
The workbench preset is unchanged; no experimental sound is silently selected.

### Exact-model recovery and implementation checks

Separately render the **actual C++ contact generator** with known noise level
1.2 and base decay 180 ms, then start at 0.5 and 70 ms. Fit only those two
parameters with each loss and the same 250-evaluation cap, keeping the same
noise seed. All recover the known settings:

| Objective | Recovered noise level | Recovered base decay |
| --- | ---: | ---: |
| Region | 1.20104 | 179.915 ms |
| Mel | 1.19968 | 180.004 ms |
| A-weighted mel | 1.20063 | 179.907 ms |
| JTFS | 1.20000 | 180.010 ms |

This validates a small identifiable subproblem, **not the full 47-parameter
joint fit**. Its provenance is saved as `synthetic-recovery.json`, deliberately
not as a fit file pretending to use the acoustic reference. Run it with
`TF_KICK_LOSS_RECOVERY=1` and `TF_KICK_LOSS_FIT=0` through the same launcher.

Eight optional loss/search tests and 49 existing fitting-tool tests pass.
The CPU/GPU JTFS score discrepancy on a non-identical signal is approximately
4.4e-7 relative. The local check uses torch 2.8.0; MLBox uses 2.13.0+cu130.
Existing locked reverb dependencies were preserved rather than downgraded to
the local test version.

### What this establishes—and does not

The established loss implementations work and can recover existing contact
noise parameters. Merely substituting them has not solved this real kick under
the common search budget. It does not establish that the DSP lacks capacity,
that JTFS is ineffective generally, or that contact noise should be removed.
The next investigation should address joint-search conditioning/recoverability
and initialization with this frozen model, rather than another exciter redesign
or untested mixture of objective weights.

## References

- [auraloss implementation and usage](https://github.com/csteinmetz1/auraloss).
- [Differentiable Time-Frequency Scattering on GPU, DAFx2022](https://arxiv.org/abs/2204.08269)
  and [authors' implementation](https://github.com/OverLordGoldDragon/wavespin/tree/dafx2022-jtfs).
- [Learning to Solve Inverse Problems for Perceptual Sound Matching](https://arxiv.org/abs/2311.14213).
- [Evaluating Sound Similarity Metrics for Differentiable, Iterative Sound-Matching](https://arxiv.org/abs/2506.22628):
  loss effectiveness depends on the synthesizer, motivating this controlled comparison.
