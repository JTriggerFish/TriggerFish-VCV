# Current percussion fitting method

See [reusable fitting safeguards](TfPercussion-fitting-lessons.md) for the shared
placement, normalization, search logging and current-renderer publication rules.

The current experiments fit one reference at a time through the exact workbench
renderer. The gong experiment below uses one Dresden mallet hit. The newer
[acoustic-kick review](TfPercussion-kick-fitting.md) uses a short-drum objective,
source-regime exploration and two-seed refinement; do not apply the gong's
six-second windows to it. Neither experiment fits a library, infers velocity
laws, or declares perceptual success from an optimizer score.

## What was wrong with the previous procedure

- `tools/calibrate_gong.py` constructs a hard-coded `CrashFit`, then renders and
  measures it. It is a bootstrap audit, **not an optimizer**.
- That script uses another event seed and resamples to 48 kHz rather than using
  the actual UI sample rate. Its manually copied contact macro expansion can
  drift from the workbench.
- Its band plots divide each signal by its own peak, hiding overall gain error.
- Its low-mode estimates use an 85 ms FFT window and broad early-time averages.
  Close frequency-bin agreement is not evidence of matching modal trajectories.
- The focused spectral refiner stops at 100 ms. The general staged crash fitter
  does include a 600 ms bloom stage and a 4 s tail stage, but only damping is
  varied in that final tail stage. Its transport search therefore stops before
  this gong's roughly one-second upper bloom peak. It is not accurate to say
  the whole earlier toolkit ignored the tail.
- Historical numerical claims predate gain/topology revisions and are withdrawn.

## Renderer and data contract

`dev.ps1 fit-gong-start` builds the optional workbench, then starts a persistent
Node/Wasm renderer and a Python optimizer. No compiler or fitting dependency is
added to normal Rack builds.

`workbench_fit_bridge.mjs` imports the actual workbench engine, parameter adapter,
reference decoder and snapshot serializer. Every parameter is checked against
the published C++ descriptor. No surrogate renderer or independent macro mapping
is used. Each evaluation resets state and uses the selected event unchanged.
The reference source hash, gain, onset and complete parameter vector are saved.
Reference gain is fixed throughout. The optimizer does not change velocity,
implement, model level, body excitation, or the velocity-response law.

## Objective: what is compared

Both signals retain absolute amplitude. No per-render or per-band normalization,
time warping, waveform-phase matching, or silent output-level correction occurs.

Five regions receive equal weight: 0–120 ms, 120–500 ms, 0.5–1.5 s, 1.5–3 s,
and 3–6 s. Each region contains two equally weighted residual blocks:

1. **Time-varying band power:** 36 ERB-spaced bands from 40 Hz to 16 kHz,
   a 2048-sample Hann STFT with 512-sample hop and 32 ms power smoothing.
   This describes where spectral energy starts, grows and disappears.
2. **Low-mode spectrum:** an 8192-sample Hann STFT with 1024-sample hop,
   region-averaged power from 40 Hz to 3 kHz with 0.6-bin smoothing.
   This prevents broad-band agreement from excusing the wrong low ringing.

Power differences are measured in dB. The floor and salience weights are derived
from the **reference only**, then held fixed. Every block is divided by its
total weight before aggregation, so more FFT bins or longer tails do not silently
dominate. Regional errors and five broad-band peak times are reported separately.
The objective is an engineering discrepancy, not a validated perceptual metric.

## What is fitted, and how

The initial experiment preserves the existing 17 centre frequencies and active
handles. It changes no topology and adds no control. The sequence is:

1. Fit the 17 visible painted levels across the full response.
2. Fit cascade rate, initial excitation tilt/centre and global turbulence
   distribution/diffusion/exchange controls.
3. Fit the **two existing T60 endpoints**, jointly with cascade rate. These are
   intrinsic damping, not copied observed decay times; band energy can also
   leave through transport. No interior knot or per-mode decay is fitted.
4. Refine the painted levels against the changed dynamics.
5. Test the existing contact/radiation controls against the remaining onset error.

Search uses bounded trust-region least squares. Coordinates are scaled to their
declared search ranges and clipped to UI bounds. Before each stage, symmetric
2%-of-range perturbations measure influence. Central finite differences at
0.5% of the range provide the Jacobian; default machine-epsilon differences are
inappropriate for float DSP with stochastic phase sensitivity. Each search is
bounded in iterations and logs the actual number of renderer evaluations.
Controls changing the residual by less than 0.05 dB for the influence probe
are frozen, avoiding arbitrary changes along nearly dead directions.
Numerical selection of an experiment does not promote it to a factory preset.

An additional 36-cell texture scan crosses the low-frequency zero-turbulence
clamp: four turbulence centres (400/700/1000/1400 Hz), three slopes
(0.15/0.35/0.6 per octave), and three common low-anchor scalers (0.12/0.5/1).
These are expanded into the existing visible controls, not extra DSP parameters.
This checks a regime boundary that a local Jacobian may not cross.

## Verification and remaining limitations

Fixtures verify identity, retained gain errors, wrong-pitch rejection, late-tail
errors despite a correct attack, and recovery of a known coloured-noise/tone
decay. The final candidate is also rendered with three held-out random seeds.
The local report presents reference versus one candidate with shared colour
limits, full audio, band trajectories, and downloadable workbench JSON.

This first search does not yet model microphone response separately, estimate
individual ridge lifetimes, or fit movement of modal centres. It may encounter
allocator discontinuities or local minima in the turbulence controls. A measured
improvement is not enough to claim these issues are solved. Inspect the plots,
regional errors, seed variation and audio, then decide what needs another pass.
If the model cannot span the observed bloom with the current controls, report
that limitation rather than adding delays, hidden gain, or extra decay knots.

Generated reference audio, PCM derivatives, reports and fits remain in ignored
`build/workbench-wasm/site/gong-review/`. They are not committed or redistributed.

## Current candidate: measured improvement, not acceptance

The review pass reduced the aggregate engineering error from about 19.5 to
7.0 dB. This is **not** a perceptual score. In the 0.5–1.5 s bloom region,
band error fell from 25.6 to about 7.1 dB; early 0–120 ms band error remains
about 8 dB. The final 1–3 kHz and 3–8 kHz peaks are around 0.57 and 0.84 s,
versus 0.49 and 0.77 s in the reference.

The directly inspected plots still show excessive initial >8 kHz energy and
low ridges that are too stationary and persistent. Several upper painted levels
reach the current +6 dB bound, so the search is constrained there. Do not expand
those bounds or add decay knots silently to make the score improve.

The two modal-loss endpoints changed from 12/1.1 s to approximately 7.24/2.01 s;
cascade rate rose from 0.5 to 1.17 octaves/s. Initial excitation became more
low-weighted, while low-packet texture was broadened using existing scalers.
Source probing found direct gain and direct high-cut effectively inactive at
this setting; their inconsequential optimizer changes were explicitly reverted.
Velocity, implement, body excitation, model level, frequencies and active
handle count remain unchanged. No synthesis code was modified by this fit.

The exported candidate was reloaded through the UI's JSON validation/adapter
and rendered sample-identically. Three fresh stochastic seeds are measured in
`heldout-seeds.json`. The factory preset is unchanged; audition the candidate
from the report and load its JSON explicitly if it is a useful starting point.
