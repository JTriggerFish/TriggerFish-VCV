# Crash: separating timing, linewidth and modal placement

This continues the [perceptual-loss review](TfPercussion-crash-perceptual-review.md).
The later [shared T60 curve review](TfPercussion-crash-decay-review.md) records
the current decay-focused audition and supersedes this page's final preset.
The reference remains Private crash A, edge, velocity 72, repeat 1. Fit duration
is six seconds at 48 kHz; independent review includes the ten-second recording
and repeated strikes. Kick and other presets are out of scope. No DSP changes
or per-mode damping multipliers are used. One extra shared T60 knot was tested
late in the pass and rejected; the audition patch retains two endpoints.

Private reproducibility artifacts are in `build/crash-perceptual-v6/`. Every
candidate includes the actual parameter vector, reference/event metadata,
renderer hash, snapshot and WAV. An optimizer result is not listening approval.

## Timing subproblem

`RegionalEnergyLoss` integrates squared output from causal second-order
Butterworth band-pass filters over explicit, non-overlapping time regions.
Reference and candidate use the same filters, region boundaries and fixed
reference floor (70 dB below the largest reference regional power). Nothing
normalizes the playback waveform or removes a fitted gain error.

The broad timing experiment uses bands 40–125, 125–300, 300–700, 700–1500,
1500–3000, 3000–6000 and 6000–16000 Hz, with time boundaries 0, 0.01, 0.03,
0.12, 0.3, 0.7, 1.5, 3 and 6 seconds. Each band/region has equal squared-error
weight. These coarse measurements cannot certify linewidth, tuning or texture.

At fixed bars, wide cascade probes (0–8 oct/s), excitation tilt/centre, energy
acceleration, transfer diffusion, neighbour exchange and both damping endpoints
only modestly improve the timing error, 3.262 → 3.075 dB. Increasing cascade
speed alone drains the lower bands too rapidly. Lower-excitation-centre/faster-
transfer starts also fail, even after observation amplitudes are refitted.

Increasing the exposed body excitation from 1 to 4 while reducing the exposed
body observation gain from 1.0724 to 0.2681 tests stronger internal energy
without a trivial fourfold output increase. This is an explicit pair of stored
UI parameters, not automatic level matching. With acceleration 1 and cascade
0.3 oct/s, observation fitting reaches 1.896 dB regional error. However, its
independent spectral error is worse than the published patch. The difference
plot still shows excess diffuse low-mid energy and missing narrow ridges.

## Efficient, exact observation-energy fitting

The existing `ObservationBasis` validates an affine combination of actual C++
renders before optimization. If the observation amplitudes are $a_i$, write

$$
y(t)=b(t)+\sum_i a_i c_i(t).
$$

The fixed intercept $b$ includes contact observation and the correction for the
basis's original amplitudes. For each linear band filter and time region,
`RegionalEnergyBasis` caches the matrix of filtered signal cross-products:

$$
P_{b,r}(a)=v^\mathsf{T}Q_{b,r}v,\qquad v=(1,a_1,\ldots,a_n)^\mathsf{T}.
$$

This retains interference terms; summing individual modal powers would not be
equivalent. The cache yields analytic amplitude derivatives without rerendering
or filtering each optimizer step. It is not a surrogate of the nonlinear DSP.
Quadratic predictions are checked against waveform measurements; mixed-gain
basis predictions and final snapshots are checked against the actual renderer.

`RegionalEnergyGuard` bounds each absolute band/time error relative to an
explicit comparator plus a declared dB allowance. It implements both signed
inequalities, avoiding an absolute-value corner in the derivative. The 0.4 dB
trial protects the improved timing while full-spectrum amplitude fitting runs;
it improves spectral error but does not resolve the incorrect ridge layout.

Reusable development-only command:

```powershell
.venv/Scripts/python.exe tools/fit_regional_observation.py crash `
  build/my-candidate build/envelope-trial measurements.json
```

`measurements.json` contains `bands` (Hz pairs) and `regions` (seconds pairs),
as used above. `EMSDK_NODE` must point to the configured Emscripten Node runtime.
The output directory must be new. This tool fits only positive observation
amplitudes at fixed body settings and the standard recorded strike. It never
publishes a preset; validate spectra and other seeds separately. All of this
remains optional analysis tooling, outside normal VCV builds.

## Modal-placement and linewidth experiments

Globally narrowing phase diffusion to 0.02 ERB and disabling neighbour exchange
exposes incorrectly placed resonances. Increasing density does not repair them.
These candidates improve coarse envelopes while worsening spectral contrast.

The next start therefore uses 23 measured peak centres and nine broader packets
at 180, 450, 800, 1300, 2200, 4000, 7000, 11000 and 15000 Hz. Peak proposals use
a 32768/4096 STFT, averaged over 0.15–1.5 s, prominence at least 4 dB, with
log-power parabolic interpolation. They are placement proposals, not physical
mode identification. Core local turbulence is 0.08; broad packets use 1. Global
turbulence is 1, spread 4 ERB, phase bandwidth 0.8 ERB, exchange 0.02 and density
1. All values are explicit existing controls. No hidden secondary resonator
bank is introduced.

Full-spectrum amplitude fitting followed by two-endpoint damping and contact
refinement gives a standard composite score of 8.653 versus 9.472 published.
Three fresh seed scores also improve (9.49 → 9.27, 9.93 → 9.10, 10.05 → 9.89).
Full-rate mel MR-STFT improves 1.3696 → 1.3457. Nevertheless, the midrange remains
too quiet during the bloom, and some attack/envelope components regress.

## Measuring the ridges' time evolution

The following timing fit adds bands 400–450, 1495–1540, 1550–1600, 2070–2110,
2850–2910 and 3690–3750 Hz to the broad bands. The first region is 0–30 ms;
the rest follow the same boundaries after 30 ms. Narrow-band filters have their
own causal settling time, so these are like-for-like filtered observations,
not instantaneous modal-energy estimates or physical T60 measurements.

Excitation centre, tilt, cascade and two T60 endpoints are fitted with bars
fixed, then the full spectral observation objective is restored. This gives
8.137 composite error, with envelope 3.512, linear-spectrum 12.444 and contrast
6.831 (published: 3.542, 15.169 and 7.371). Fresh seed composite scores improve
9.61 → 8.97, 10.00 → 9.09 and 9.91 → 8.81. The initial 1–3 ms remain too strong;
contact refinement and final audition/publication review follow this result.

The subsequent contact fit reaches 8.117 composite error. It improves the
0–1 ms level difference to +1.10 dB, but the 1–3 ms difference remains +4.68 dB.
Mel MR-STFT is 1.3916, slightly worse than the published patch despite the
better engineering score. A constrained mel polish is therefore tested rather
than declaring agreement between the objectives.

One diagnostic explicitly tests the UI blend's fixed half-sine pulse floor.
Native C++ and workbench baseline waveforms agree to relative RMS 0.0000699
before the ablation. Removing only that pulse changes the 1–3 ms excess from
4.68 to 4.53 dB and worsens overall envelope error from 3.54 to 6.50. This does
not support blaming the pulse floor for that spike or changing the model on
that assumption. These raw-C++ variants are marked non-UI-compatible and are
not published as fitted patches.

## Search bounds were confounding the contact comparison

Testing contact tone/noise blends of 0.15, 0.35 and 0.6 initially suggested that
the quieter-noise start could not reproduce the ridges. Several important bars
were actually pinned at the optimizer's +6 dB limit while the exposed body
observation gain was only 0.268. This was a search-coordinate limitation, not
evidence of insufficient synthesis capacity.

The corrected experiment raises the explicit observation gain to 4 and lowers
all bars by the corresponding 23.47 dB before searching. Starting waveforms
agree within 0.0003 relative RMS. The optional amplitude search range is widened
downward to −70 dB, still inside the UI's active range; its normal default stays
−45 to +6 dB. Bounds are recorded in the search history. With the freed headroom,
the 0.15 trial's composite score improves from 8.788 to 8.226. There is no
automatic waveform normalization and no new DSP parameter.

For presentation, the selected 0.35 candidate is then rescaled once so its
largest painted bar is 0 dB and the exposed observation gain is 0.6632. Their
products remain unchanged (waveform relative RMS difference 0.000000277).
Body excitation remains 4; master level remains −0.76224 dB and reference gain
42 dB. These are stored settings, not a level-matching process on trigger.

## Perceptual checks and rejected refinements

The native full-rate auraloss mel MR-STFT objective was used both for evaluation
and constrained amplitude optimization. The differentiable adapter is checked
against the library scorer and finite differences. The constrained result can
reduce mel error while worsening attack, or regress on fresh synthesis seeds.
It is therefore not promoted solely because its training loss falls.

JTFS was independently evaluated through the existing isolated MLBox worker.
This uses a 16 kHz analysis copy; full-rate mel and spectrogram checks are still
needed for the upper cymbal spectrum. JTFS is a structured time/frequency
representation, not a human listening verdict or a learned listener model.

The following comparison uses the original standard strike:

| Measurement (lower is better) | Previous workbench | Audition candidate |
| --- | ---: | ---: |
| Engineering composite (mixed units, not dB) | 9.4722 | 8.1296 |
| Full-rate mel MR-STFT | 1.3696 | 1.3846 |
| JTFS | 0.1282 | 0.0977 |

The same three development-audit seeds give composite scores of
9.8933 / 9.6959 / 9.4496 before and 9.5898 / 8.7082 / 9.4004 after. Mel errors
are worse on those seeds. These seeds have been inspected during selection:
they are **not** an untouched holdout set. All tests use one reference recording,
not three different reference strikes. An ensemble-trained variant is steadier
across seeds but does not improve all standard-strike measurements either.

Full-tail inspection identifies a substantial unresolved timing error:
250–1000 Hz is about 8–9 dB too quiet in the first 120 ms, but slightly too loud
after 1.5 s. As a final sparse damping experiment, a single 450 Hz T60 knot is
tested at 5, 7 and 9 s. The low endpoint is adjusted to preserve damping at
125 Hz, and observation amplitudes are refitted for each trial. Composite
scores are 9.522, 8.577 and 8.182; envelope errors are 7.045, 4.727 and 3.659.
None justifies retaining the extra point. Simply adding damping flexibility
does not resolve the discrepancy in this experiment.

## Workbench handoff and remaining question

`workbench/web/crash_calibration.fit.json` now contains the **unfinished audition
fit**, sourced from `build/crash-perceptual-v6/audition/`. Its parent ID identifies
the previous committed preset. No other instrument preset or DSP implementation
changes. Snapshot reloading is exact; the six-second fitting window is also
reviewed over the full ten-second reference. Eight repeated strikes at 250 ms
spacing are finite: raw pre-limiter peaks are 0.835 at the standard velocity and
2.244 at full strength. The browser safety limiter must remain enabled; these
are not claims of limiter-free headroom or listening approval.

This is evidence of **methodology problems**, particularly coarse-loss tradeoffs,
modal placement and search conditioning. It does not yet establish a necessary
model replacement. The signed plot still shows an early low-mid deficit and
excess diffuse energy between ridges. The next useful model-capacity test is a
joint excitation/timing/linewidth fit with those discrepancies explicitly
guarded, not another acceptance decision based on an aggregate score alone.

## Audition-guided bloom refinement

The user's next audition judged tuning substantially improved, but requested
less noise, more bloom and faster high-frequency decay. This is a **fixed-tuning**
pass: all 32 modal centres, their membership, contact parameters, strike event,
reference gain and master level stay unchanged. Private artifacts are under
`build/crash-bloom-v7/`.

One-control probes separate cascade rate, the upper T60 endpoint, phase
diffusion, transfer diffusion, turbulence and packet spread. Faster transfer
alone also prolongs the upper output tail. Shorter upper damping alone suppresses
too much bloom. The reference comparison further distinguishes deficient
4–8 kHz body energy from excess 8–15 kHz persistence; a global treble boost
cannot correct both.

The joint trial screens three upper damping endpoints, three cascade rates and
two phase bandwidths. During this screen the lower damping endpoint compensates
to preserve the existing 1 kHz T60, with both endpoint values stored explicitly.
It then fits the two endpoints, cascade rate and four upper packet observation
levels using finite differences through the actual workbench renderer. The
search records control influence before optimizing. Modal frequencies, lower
bar levels and per-mode damping are excluded. No intermediate T60 knot is added.

The trial objective combines regional absolute band/time energy (75%) and the
existing ERB-weighted ridge-contrast error (25%). Bands are 40–250, 250–1000,
1000–2000, 2000–4000, 4000–8000 and 8000–15000 Hz; region boundaries are 0,
0.03, 0.12, 0.3, 0.7, 1.5, 3 and 6 s. These weights express this refinement's
focus, not a validated perceptual distance. Full-spectrum, mel, signed plots,
fresh-seed comparisons and the ten-second/restrike audit remain separate checks.

The new audition uses cascade 0.1603 oct/s (previously 0.06385), phase bandwidth
0.35 ERB (0.8) and transfer diffusion 0.35 (0.6915). The upper T60 endpoint is
1.437 s (2.920); the lower endpoint rises to 29.416 s (16.832), compensating
for increased transport out of the lower packets. The four upper packet levels
also change; all other parameters remain exactly fixed. These endpoint values
are damping settings, not a claim that observed band tails follow those T60s
when energy is arriving from lower frequencies.

At the standard strike, 4–8 kHz energy during 0.3–0.7 s moves from −2.99 to
−1.68 dB relative to the reference. The 8–15 kHz excess during 1.5–3 s falls
from +4.48 to +1.91 dB. Fresh seed offsets 70001, 71011 and 72019 all reproduce
the shorter upper persistence (excess falls from +4.94…+6.36 to +1.59…+3.20 dB).
They are not used to tune the candidate. Mel error improves on these three
seeds, but changes from 1.3846 to 1.3897 on the standard strike. Regional error
improves on the standard strike and regresses slightly on those seeds.

The signed plot was inspected before publication. Spectral magnitude error
still regresses; less diffuse texture and the improved upper-band timing are
the reason to audition this version, not a claim that all metrics improved.
The ten-second and repeated-strike audit passes with exact snapshot reload.
Raw repeated-hit peaks are 0.806 at standard strength and 2.054 at full strength,
before the browser limiter. This replaces only the Crash workbench preset;
the prior audition is recoverable from `build/crash-perceptual-v6/audition/`.
The new artifact set is `build/crash-bloom-v7/audition/`.
