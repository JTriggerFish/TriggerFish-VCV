# Crash refinement: low phase blur and structured spectra

This pass follows the user's positive audition of the beating-packet gong.
The gong preset, DSP topology and UI control count are unchanged. It tests
existing degrees of freedom before adding another synthesis mechanism.

The subsequent [audition-directed refinement](TfPercussion-crash-audition-refinement.md)
records the currently published broader-shimmer trial; results below describe
the earlier low-blur starting point.

## Research and working hypotheses

[Skare and Abel (DAFx 2019)](https://www.dafx.de/paper-archive/2019/DAFx2019_paper_48.pdf)
describe high-Q modal banks whose adjacent resonances produce evolving beats.
They distinguish this from reproducing the nonlinear attack/bloom; their
examples use thousands of modes. This supports testing low stochastic blur,
but does not demonstrate that our 512-mode budget is sufficient.

[Nguyen and Touzé (JASA 2019)](https://perso.ensta-paris.fr/~touze/PDF/2019_BaoCT_TaperedCymbals_JASA.pdf)
show how thickness variation helps edge strikes develop high frequencies very
rapidly, whereas gong-like development can be much slower. They also emphasize
frequency-dependent damping. Their [simulation examples](https://perso.ensta.fr/~touze/tapercymbals.html)
are useful qualitative references, not measurements of our target cymbal.
We test faster transport together with initial excitation and damping; we do
not equate a larger diffusion exponent with faster transport, nor introduce
a hidden delay or independently generated high-frequency noise tail.

The stretched harmonic family is a low-dimensional design prior, **not** a
claim that crash-cymbal eigenfrequencies form one harmonic series. Root, upper
stretch and protected-core size select a whole family. Fitting individual
unidentifiable upper ridges remains excluded.

## What is fitted, and how

`tools/refine_crash_beating.py` uses the actual workbench Wasm. The target is the
existing medium edge strike with its stored onset, velocity, implement and
reference gain. Model level, body drive, observation gain, contact parameters,
output filters and local allocation/noisiness weights initially stay fixed.

1. Archive the currently published preset and render as `baseline`.
2. Screen three layouts, densities 0.45/0.85/1, and phase-blur coefficients
   0/0.001/0.004/0.01 with unchanged geometry and observation levels.
3. Screen 32-centre families with roots 110/120/130/145 Hz, protected harmonic
   cores 2/4/6, and upper stretch computed by the **same UI generator** to
   place the last centre at 13.5 kHz. No individual frequency fitting.
4. Fit six smooth observation-amplitude coordinates (120, 400, 1000, 2500,
   6500, 15000 Hz) using a verified affine audio basis and exact STFT
   cross-power derivatives. This writes the visible modal bars, not hidden EQ.
5. Refine shared excitation, transport, noisiness/spread/blur and two T60
   endpoints using bounded Powell search, preceded by finite-difference
   sensitivity checks. Save the actual bounds and evaluation trace.
6. Apply differentiable Mel observation polishing, constrained not to give
   back the fitted band envelopes/rise beyond the recorded tolerance. Reject
   the polish if the full proposal score worsens.

The proposal score is Mel MRSTFT + 0.3 attack Mel (first 300 ms) + 0.05 spectral
bloom norm + 0.3 envelope-modulation texture distance. This increases the
importance of attack/spectral balance compared with the previous texture-led
trial. It remains an engineering ranking, **not a validated perceptual loss**.
Final checks must report its separate components, additional random seeds,
band decay, paired/difference plots and repeated strikes. A low score does not
constitute listening approval.

The phase-blur slider is a coefficient, not the final linewidth in Hz. The
existing implementation multiplies it by squared local noisiness and the
frequency's ERB bandwidth (the coherent centre also has its existing factor).
A small coefficient can therefore still blur noisy upper packets appreciably.
Tests must consider the resulting sound and linewidth, not compare the raw
gong and crash slider numbers as if they meant the same bandwidth.

## Reproduction

Use the MinGW/optional-Wasm build through `dev.ps1`, then the development Python
environment with `EMSDK_NODE` pointing to the configured Emscripten Node.

```powershell
.\dev.ps1 build-workbench
.\.venv\Scripts\python.exe tools/refine_crash_beating.py screen --output build/crash-low-blur
.\.venv\Scripts\python.exe tools/refine_crash_beating.py refine --source build/crash-low-blur/series-130-4 --output build/crash-low-blur-refine-130 --budget 480
```

All stages produce unpublished checkpoints. Workbench publication is a separate
reviewed action and must reproduce the saved candidate audio exactly. None of
these Python tools or analysis dependencies are needed by normal Rack builds.

## Identified loss bias and corrective test

The ordinary six-second Mel MRSTFT loss overweights the quiet upper tail in
this reference. Two controlled edits of the **reference itself** demonstrate
the issue (fourth-order zero-phase 6-kHz high-pass component, smooth 100-ms
crossfades; these are diagnostics, not edits to the fitting target):

| Reference perturbation | Ordinary Mel | Reference-floor Mel, 60 dB |
|---|---:|---:|
| Reduce upper tail after 3 s by 20 dB; difference RMS −74.8 dBFS | 0.33563 | 0.000325 |
| Reduce upper content during the first second by 6 dB; difference RMS −49.7 dBFS | 0.08076 | 0.05992 |

The old ranking favours preserving quiet late hiss over a much larger early
spectral difference. `ReferenceFloorMel` retains auraloss's Mel projection,
spectral convergence and log-magnitude terms, but derives **each resolution's
linear FFT power floor from the reference peak**, 60 dB below it. The floor
is fixed before comparing any candidates and before the Mel projection. FFT
window scales differ, so each resolution needs its own threshold. Neither the
reference/candidate audio nor display normalization changes. No gain matching.

This is a declared analysis dynamic range, not a measured absolute hearing
threshold or a full masking model. The 50/60/70-dB audit exposes sensitivity
to this choice. Unit tests cover the late-noise/early-colour ordering, retained
gain sensitivity, fixed thresholds, invalid references and agreement of the
double-precision autograd adapter and numerical derivatives.

Pass `--reference-floor-db 60` to use this comparison in the crash fitter.
Its objective specification changes to `crash-reference-floor-v2`; scores
from the ordinary and reference-floor versions are **not comparable**.
`review_metal_refit.py --floor-audit` reports both, alongside the existing
decay/modulation diagnostics. The new loss still cannot establish sound quality.

Additional experiments:

- `contact`: five shared exciter/presentation controls, same playing gesture.
- `dynamics`: shared excitation/diffusion/two-endpoint damping, with the
  chosen packet geometry/noisiness/phase blur held fixed.
- `stable`: zero phase blur, full oscillator budget, doublet layout; no
  compensating increase in stochastic linewidth is possible.
- `decay`: test one shared 600-Hz knot after the main fit, retain only if
  the proposal score improves by more than 0.02. No per-mode decay fitting.
- `screen_stable_crash_texture.py`: existing smooth frequency drift and a
  whole frequency-dependent allocation curve, not individually fitted weights.
- `refine_crash_front.py`: faster, explicitly energy-dependent initial
  transport as a research-motivated starting hypothesis, not a new DSP model.

The nominal blur coefficient alone was an inadequate fitting constraint: one
trial reduced it while raising packet noisiness, leaving the calculated 10-kHz
satellite linewidth at 332 Hz versus the starting 483 Hz. The final low-blur
experiments hold packet texture fixed while fitting other controls, preventing
that compensation. This is why parameter-value changes alone are not evidence
that the intended acoustic change occurred.

## Published audition trial: 9 September 2026

The main workbench's **Crash — low-blur structured trial** is the verified
`build/crash-linewidth9/candidate` snapshot, not the lowest score from every
experiment. Its parent is the archived `build/crash-low-blur/baseline`.
Gong is unchanged. This is a partial result awaiting listening approval.

Retained configuration:

- Original 24-centre stretched series: 120-Hz root, four protected harmonics,
  stretch 0.7. No individual frequency adjustments or local allocation fits.
- Original scattered layout, density 0.85, packet spread and noisiness curve.
  Smooth frequency drift remains zero. Doublets are useful for the gong but
  did not establish a better crash fit in this experiment.
- Phase-blur coefficient 0.0842177 → 0.0015: calculated satellite linewidth at
  10 kHz falls from 483.1 to 8.6 Hz. This is a phase-process parameter, not an
  independently measured linewidth of the complete rendered spectrum.
- Shared transport rate 1.661 → 2.361, energy exponent 0.0639 → 0.1008;
  initial excitation tilt −4.88 → −7.20 dB/oct and centre 929 → 2225 Hz.
- Contact mix, width, chirp pitch, noise tilt and direct observation gain
  retuned together; six broad modal-observation coordinates retuned afterward.
- Two active T60 endpoints only: 25.83 s and 0.655 s. No per-mode decay fits.
  Model level, body gain, reference gain and playing gesture remain fixed.

Mean independent audit results over the standard seed and three additional
seeds are below. Lower is better for every row, but units differ by metric.
The additional seeds were not used for the Powell/observation fits.

| Diagnostic | Previous preset | Published trial |
|---|---:|---:|
| Reference-floor Mel, 50 dB | 0.909 | 0.874 |
| Reference-floor Mel, 60 dB | 1.007 | 0.965 |
| Reference-floor Mel, 70 dB | 1.091 | 1.030 |
| Band-decay shape, dB | 3.543 | 3.232 |
| Ordinary attack Mel | 1.313 | 1.274 |
| Ordinary full Mel | 1.599 | 2.111 |
| Envelope-modulation texture | 0.280 | 0.352 |
| JTFS, restricted below 8 kHz | 0.1124 | 0.1141 |

These are **mixed**, not universal improvement: decay improves on three of
four seeds; reference-floor Mel improves on two, is almost unchanged on one,
and worsens on one. Texture worsens. The 17.2-Hz linewidth variant has nearly
identical reference-floor scores and slightly better texture scores; the
8.6-Hz trial is exposed for the user's low-blur listening preference, not
because every numerical metric prefers it.

Visual inspection of absolute band envelopes and paired/difference STFTs
shows closer upper-band decay, but excess around 200–300 Hz, deficient early
400-Hz energy and missing detailed reference ridges remain. Late reference
high-band power plateaus while the model dies away. The fixed analysis floor
limits the leverage of this quiet plateau; it does not erase it from plots
or audio. Neither matching the plateau nor reducing this loss proves a
realistic crash. No claim of having listened to the render is made.

Repeated quarter-note strikes peak at −3.36 dBFS before master volume;
rapid hard strikes reach +8.89 dBFS in floating-point synthesis (−3.11 dBFS
after the normal −12-dB browser master). The browser safety limiter remains
in place. There is no hidden normalization or synth clipping introduced here.

The saved snapshot, candidate WAV and fixed reference reproduce exactly with
the rebuilt workbench. All 551 Python tests, 14 optional Wasm tests and two
native workbench API tests pass. A silent disposable browser test loads crash,
gong, ride and hi-hat without changing the user's tab or playing audio.
An initial float32-versus-float64 loss
discrepancy correctly stopped polishing; scoring now uses float64 like the
autograd adapter. The agreement tolerance was not relaxed.
