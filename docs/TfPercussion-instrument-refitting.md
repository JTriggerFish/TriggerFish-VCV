# Reference-specific instrument refitting

This pass keeps the user-accepted kick unchanged and fits one standard reference
each for snare, crash, ride, half-open hi-hat and gong. It uses the **actual C++
WebAssembly engine**, not a surrogate synthesizer. A better numerical fit is a
candidate for audition, not proof that the instrument sounds right.

The subsequent parallel crash/ride/gong pass is documented in
[Shared metallic fitting v2](TfPercussion-metallic-fitting-v2.md). It replaces the
objective for those experiments, adds validated analysis autograd, and corrects
quiet-mode search coordinates. The original pass below remains recorded for
provenance rather than being relabelled as the new method.

## What is fitted

The reference-target catalog fixes the sample, onset, corpus gain, strike
strength, implement, hardness, location and pedal constraint. Those performance
inputs are **not optimization variables**. Audio is neither peak-normalized nor
time-warped. The visible, attenuation-only model level can be fitted explicitly;
there is no automatic output correction at playback time.

For snare, bounded blocks adjust body pitch/character and damping, wire response,
contact/body balance, ring frequency and the existing radiation bandwidth. Wire
noise balance is adjusted separately with its modal mix anchored; master gain,
wire level and both mix gains are not simultaneously free.

For metallic instruments, the variables are modal frequencies and observation
levels, initial energy distribution, passive upward energy travel, turbulence,
contact presentation, and the **two endpoint T60 values only**. Interior decay
knots and per-mode decay multipliers are not fitted. Global observation gain is
held fixed while painted modal levels are fitted, avoiding a gain nullspace
without preventing a necessary change in overall body/contact balance. No DSP topology or velocity
curve is changed by this procedure.

## Starting modal layouts

The factory crash's first handle is 421 Hz. A sample with defined ringing below
that needs lower handles, not merely more turbulence or longer damping. Candidate
layouts include both:

- retaining the current layout and adding up to four measured lower peaks;
- a smaller measured-peak layout with explicitly identified broadband coverage
  handles at higher frequencies.

Peak proposals use the existing long-window, prominence-ranked spectral tool.
They are **not a physical modal decomposition**, and measured spectral power is
not treated as a modal coefficient. Alternative layouts receive a short bounded
refit of observation gain, excitation tilt and endpoint damping before comparison.
Rejecting an untuned layout by its raw score can reject the better architecture.
Selected frequencies subsequently receive local refinement.
Local frequency refinement is limited to eight salient handles below 3 kHz;
small moves of every broad high-frequency packet create weak, seed-sensitive
search directions. Broadband colour is handled by observation levels and texture.

## Objective and search

Snare uses `DrumBalanceLoss`: multiresolution spectra plus identically filtered,
causal band-power envelopes, with separate contact, attack/body and tail regions.
Its optimization horizon is 1.2 seconds.

Metallic instruments use `MetallicFitLoss`, over six seconds:

| Share | Measurement | Purpose |
| --- | --- | --- |
| 60% | ERB-band trajectories and resolved low/mid spectra | Spectral energy travel, body ringing and decay |
| 20% | 512-sample STFT, first 120 ms | Contact colour without long-window attack smearing |
| 20% | 4096-sample regional spectra, 1–16 kHz | Upper ringing that the older low-frequency-only ridge term missed |

All shares refer to squared residuals; each uses a reference-derived floor and
salience weighting. Metallic regions are 0–0.12, 0.12–0.5, 0.5–1.5, 1.5–3 and
3–6 seconds. Region errors and peak times remain available individually. Scores
from different objectives or instruments are **not directly comparable**.

`Search.stage` uses bounded SciPy least squares, parameter-range-scaled central
finite differences, explicit sensitivity probes and cached exact renders. Dead
directions remain fixed. A step must improve on the actual previous parameter
vector, not a silently clamped replacement. Short staged passes are warm-started
when the updated contact or texture changes the useful body solution.

Painted observation amplitudes have a faster, separately validated path. With
active membership, excitation, frequency, damping and texture fixed, they do not
affect the stored state, so

$$
y(t;g)=c(t)+\sum_i g_i s_i(t),\qquad g_i=10^{d_i/20}.
$$

The contributions $s_i$ are obtained from differences of **actual C++ renders**,
not synthesized in Python. Two mixed-gain probes per training seed check the
affine prediction against new full renders before optimization. Turning a mode
off or changing another parameter invalidates the basis. This permits rapid
observation-level iterations without repeating the nonlinear state simulation.
Selection and saved WAVs always use a fresh full render; prediction roundoff is
never passed off as an exact saved-engine result. The kick code is untouched.

Two fixed random seeds train each candidate. Three different seeds are evaluated
after fitting; stochastic audio is never averaged before measurement. Known-answer
tests cover exact identity, literal gain, altered attack, added upper ringing and
decay direction, alongside the existing coloured-noise T60 recovery tests.

## Review and publication

Inspect same-scale reference/candidate spectrograms, attack and body spectra, and
four broad-band power envelopes. Inspect remaining differences explicitly: an
aggregate improvement cannot waive a missing bass mode, delayed attack or noisy
ring. Also verify finite audio, repeated strikes and exact snapshot reload.

Full fitted JSON includes the reference identity, fixed event and every editable
parameter. `reference_calibration_library.mjs` connects reviewed candidates to the
existing **Reference targets** toolbar. The main workbench is the audition UI;
there are no competing report servers or alternate “new” sounds. Existing browser
edits are not overwritten; reselect the target after reloading to audition a new
published candidate. User listening acceptance is tracked separately from checks.

## Reproduction

Use the optional developer launcher:

```powershell
$env:TF_FIT_TARGETS = 'snare' # or crash, ride, hihat, gong; comma-separated
$env:TF_FIT_ITERATIONS = '12'
./dev.ps1 fit-instruments
```

Output defaults to `build/instrument-refits-v1/<instrument>`. It contains the
immutable baseline, reference WAV, candidate WAV/fit, full parameter bounds and
stage history, sensitivity results, seed identities and renderer hash. Set
`TF_FIT_RESUME=1` to resume a provenance-checked search; `TF_FIT_LAYOUTS=1` requests
another comparison of modal layouts. Set `TF_FIT_OUTPUT` for a separate experiment.
`TF_FIT_PHASE=surface` resumes after layout/damping selection to refine observation,
texture, local frequencies and contact. `./dev.ps1 polish-instruments` performs
only the validated observation-level refinement on completed metallic candidates.
The fitter does **not** automatically publish. Private offline Plotly inspection
images can be generated using `tools/instrument_fit_plots.py` and
`tools/capture_fit_plot.mjs`; these do not create another audition page.

Python, browser tooling and WebAssembly remain optional development dependencies,
not requirements of normal Rack builds.

## Earlier pass results and limits — 2026-09-07

This section is a historical record, **not the current metallic calibration**.
The [subsequent shared crash/ride/gong pass](TfPercussion-metallic-fitting-v2.md)
supersedes those three rows. In particular, its reference-only onset audit found
the ride marker early, so these earlier ride timing scores must not be reused.
The two-endpoint statement below describes this earlier pass; the later ride
experiment documents a single explicit, measured damping-knot exception.

At the end of this earlier pass, all five candidates were saved in the main workbench. These were improvements on
their previous presets, **not five accepted matches at the kick's quality**.
The plots were inspected; listening acceptance remains pending.

| Target | Previous objective | Candidate objective | Remaining mismatch |
| --- | ---: | ---: | --- |
| Snare | 20.71 | 5.03 | Body/ring detail; strongest improvement alongside hi-hat |
| Crash | 21.58 | 6.92 | Midrange too diffuse; upper ringing persists too long |
| Ride | 19.25 | 8.71 | Attack too weak, stable ridges insufficiently defined, tail too short |
| Half-open hi-hat | 10.49 | 5.67 | Low-frequency detail and late decay |
| Gong | 21.55 | 6.95 | Initial low mode clearer, but upper bloom trajectory still wrong |

Compare scores only within a row. They are objective units, not perceptual grades.
Three held-out random seeds also improve on each corresponding baseline. That
checks stochastic robustness, not generalization to other velocities or samples.

Full-reference-duration audits extend beyond the optimization horizons, including
the ride's 22.4-second recording. This matters: a six-second objective cannot
certify its later decay. The ride's reference/model peaks are approximately
0.77/0.18 at the fixed source gain, so its overall score also fails to express an
important attack-level discrepancy sufficiently strongly. No playback gain
correction conceals that difference.

Further experiments tested jointly coherent texture settings, measured late-ridge
layouts, and low-fed stronger bloom, including short refits before selection.
None supplanted the retained candidates. One explicit ride experiment added a
single 1 kHz damping knot as a final refinement; it was rejected too. **Every
published metallic candidate still uses just two T60 endpoints**, with no fitted
per-mode decay multipliers. These failed experiments do not establish that the
DSP cannot fit the samples; they show that these particular starts/search blocks
did not find a better solution.

The developer-only polishing switches `TF_FIT_TEXTURE_SCREEN`,
`TF_FIT_ENERGY_SCREEN`, `TF_FIT_LATE_RIDGES`, `TF_FIT_LOW_FED_ALTERNATIVE` and
`TF_FIT_RIDGE_DAMPING` reproduce those experiments when set to `1`. Use them
individually in a separate `TF_FIT_OUTPUT`; named alternatives resume an existing
provenance-checked child search rather than discarding it. The damping experiment
requires the preceding coherent-tail-layout result. These switches are not
instrument controls or part of ordinary fitting defaults.

Verification includes full-render snapshot reproduction, standard/max-strength
repeated hits, 438 Python tests, native/Wasm parity, and browser checks that each
Reference target saves back the exact candidate parameters, event and reference
hash. Repeated maximum-strength strikes can exceed unity before the browser
limiter; no new limiter or normalization was added to the core. The accepted kick
preset and C++ DSP are unchanged.
