# Constrained attack-energy fitting

This optional developer tool refines painted observation amplitudes in an
existing, verified workbench fit. It is **not** a new synthesis model, a playback
normalizer, or a declaration that a fit sounds correct. Normal Rack builds do
not use Python, Torch, this optimizer, or the browser renderer.

## Why use an explicit constraint?

An average spectral/envelope loss can improve the sustained body while allowing
the first milliseconds to become much too weak. Five separate RMS measurements
cover 0–1, 1–3, 3–10, 10–30 and 30–100 ms after the authoritative reference onset.
Every training seed must satisfy each interval, rather than passing after
averaging seeds. The default limit is ±3.5 dB. This is an engineering rejection
limit, not a measured perceptual threshold.

Total RMS is insufficient on its own: an incorrect low-frequency impulse can
supply the missing attack power. Always inspect short-region **band energies and
spectra** alongside these constraints. A passing RMS gate cannot excuse false
bass, missing high-frequency content, or incorrect modal decays.

## Exact procedure

1. Verify the starting JSON, source audio identity, source gain/onset, renderer
   hash and exact actual-Wasm render. Never align the candidate independently.
2. Render the affine observation basis with the actual C++ engine for the
   standard event seed and that seed plus 101. Validate the basis using fresh
   actual renders and validate Torch/NumPy measurement parity.
3. Optimize active painted observation amplitudes only, in linear-amplitude
   solver coordinates. Exposed JSON/UI values remain dB. Bounds are −45 to
   +6 dB; an active starting value outside these bounds is rejected. Inactive
   handles stay inactive. All excitation, transfer, damping, global output,
   strike inputs and reference gains are fixed.
4. SLSQP minimizes the shared `MetallicBalanceLoss` with ERB-weighted contrast
   and fast-attack measurements. Gradients pass through the analysis and the
   exact rendered affine basis, **not** through a substitute synthesizer.
   The five RMS inequalities have analytic amplitude Jacobians.
5. Retain the best feasible evaluation, including the starting point when
   feasible. Reaching the iteration limit does not mean convergence. No fit is
   saved if no feasible point was observed, or if results are nonfinite.
6. Rerender both training seeds through C++, check the constraints again
   (0.02 dB numerical tolerance), then save parameters, bounds, fixed controls,
   objective specification, validation and actual errors. Refuse to overwrite
   an existing search directory or the starting fit.
7. Independently inspect fresh seeds, full-duration band decays, first-30-ms
   spectra, shared-scale plots, and repeated strikes. Report failed held-out
   intervals rather than tuning those seeds until they pass. Audition the
   result in the main workbench before calling it a successful calibration.

## Invocation

Build the optional workbench through `dev.ps1 build-workbench` and initialize
the existing Emscripten environment (`EMSDK_NODE`). From the repository root:

```powershell
.venv\Scripts\python.exe tools/refine_workbench_attack_gate.py ride START_DIRECTORY OUTPUT_DIRECTORY --seconds 12 --iterations 25
```

The process explicitly limits OpenBLAS and Torch to one CPU thread to avoid
oversubscribing concurrent instrument experiments. It does not start a server,
change source presets, or publish additional report pages.

`tests/python/test_observation_energy_gate.py` checks literal gain, independent
seed constraints, finite-difference Jacobian agreement, the silence floor,
invalid/overflowing data, input bounds, best-feasible selection, and overwrite
protection. Source implementation lives in
`python/triggerfish_percussion/observation_energy_gate.py` and the CLI above.

## Ride: final contact-observation check

The September 2026 ride trial exposed precisely the RMS loophole described
above. A constrained body fit recovered the initial total energy but its
first-30-ms 40–250 Hz band was 21.32 dB too strong. Increasing the **existing
contact high-pass** from 40 Hz to 1 kHz removed this false bass without changing
the stored body; simply accepting that change would have weakened the short
RMS bins too much.

A separate bounded contact-observation fit therefore used the same two
training seeds and fixed every contact-generator and body parameter:

| Exposed control | Bounds | Selected value |
| --- | --- | --- |
| Contact observation gain | 0.5–2 | 0.8426893 |
| Contact high-pass | 400–1000 Hz | 702.878 Hz |
| Contact shelf frequency | 1500–12000 Hz | 10902.382 Hz |
| Contact shelf gain | −8–12 dB | 12 dB |

For each seed, the objective was the mean squared error in dB across the five
short RMS bins **and** the four first-30-ms band powers: 40–250, 250–1000,
1000–4000 and 4000–16000 Hz. Bands use second-order Butterworth band-pass
sections, causal filtering from zero state, then mean squared output. These
are the same measurements as `tools/review_instrument_fit.py`, not frequency
bins silently cropped from a long STFT.

SLSQP constrained all 18 signed errors independently to ±3.5 dB. Controls were
normalized to their stated bounds; central finite differences used ±0.005 of
each normalized span, clipped at bounds. The initial point was gain 1.4,
high-pass 1000 Hz, shelf 4000 Hz/+3 dB. Maximum iterations were 30 and `ftol`
was 1e−6. Each evaluation rendered **0.2 seconds through the actual Wasm
engine**, caching identical control vectors; the final selected parameters
were rerendered for 12 seconds and scored with the shared full objective.
The best feasible evaluation was retained, not an infeasible final iterate.

The selected standard event has first-30-ms band errors
`[+0.91, −1.14, +2.02, −0.30] dB` and short RMS errors
`[−0.28, +1.30, −3.42, −1.84, +2.42] dB`. Its shared full score improves
from 13.041 to 9.333. These are improvements, not equivalence: held-out seeds
still show some larger errors, and the 6–12-second low-mid modal tail remains
too weak. No per-mode damping or additional damping knots were introduced
by either final constrained pass.

The reusable source CLI `tools/refine_workbench_contact_shape.py` implements
this contact-only procedure separately from the observation-amplitude tool.
Its bounds and diagnostic starting point are explicit in the source and saved
history; it does not silently clamp the starting patch. The following command
was executed successfully to reproduce the selected parameters exactly:

```powershell
.venv\Scripts\python.exe tools/refine_workbench_contact_shape.py ride build/instrument-refits-v3/ride/corrected-onset/constrained-final build/instrument-refits-v3/ride/contact-tool-verification --seconds 12 --iterations 30
```

These local build paths are optional experiment outputs, not required source
assets. Supply any verified compatible input directory and a fresh output
directory. The script saves only after a fresh full-duration render satisfies
both spectral and energy constraints; it does not update the workbench preset.
