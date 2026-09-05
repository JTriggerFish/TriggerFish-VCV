# Crash cymbal calibration baseline

## Status

The current `CrashCymbal` graph is an experimental, deterministic instrument
built from tested percussion primitives. It is a useful starting point for
interactive listening and ablation, but it is **not a calibrated crash model**.
Earlier numerical checkpoints did not survive direct listening and are not
accepted baselines. In particular, a passing aggregate or prefix metric must
never be described as a perceptual fit.

Calibration is now fitting-by-ear first. Numerical analysis remains a set of
diagnostic views and regression checks; it does not approve sound quality.
Automated search code is retained as experimental developer tooling. It may
produce a neutral starting state for the ear-fitting workbench, but it never
approves a fit. The workbench base now contains only documented neutral DSP
defaults and a fixed output-unit calibration; it contains no corpus-specific
modal phases, velocity exponents, contact fit, or instrument preset.

## Failure review and corrected fit contract

An earlier automated pass accepted several stages on relative improvement
alone. The unified initial-body, bloom, and tail stages had their absolute
quality and acceptance gates disabled, so a candidate could be less wrong than
an already-wrong seed while remaining audibly unrelated to the reference.
Automatic level matching then rewrote the model-level parameter, sometimes
above 0 dB, and concealed the underlying output-calibration error. Those were
methodological bugs, not evidence that the synthesis graph was calibrated.

The corrected contract is:

1. align the physical onset before measuring, while retaining pre-onset audio;
2. estimate level offset separately and never include it in shape fitting;
3. divide the response into contact (`0..15 ms`), initial body (`15..120 ms`),
   bloom (`120..600 ms`), and tail (`600 ms..end`) regions;
4. optimize only controls whose finite differences materially affect the active
   region, while protecting every earlier region from regression;
5. require absolute regional quality gates at every promotable unified stage;
6. preserve the complete vector of failures—envelope, ERB trajectory,
   fine-spectrum, centroid/rolloff, flatness, crest, and ridge/noise balance;
7. audition the candidate and relevant branch solos before it may be named a
   calibration.

Modal T60 is handled by a narrower contract. Factory and calibration starts use
only the fixed 40 Hz and 15 kHz curve positions; the six interior slots are
inactive. Contact, 100 ms colour, and bloom stages cannot change T60. The tail
schedule is the only generic schedule allowed to expose the two endpoint
values, but its result is not accepted from the generic spectrogram objective:
it must also pass the dedicated ERB-band decay measurement. Interior knots are
a sparse, last-step exception after other topology and colour errors have been
resolved, not extra degrees of freedom for reducing a spectrogram loss. The
controlled coloured-noise recovery fixture and numerical procedure are
documented in `TfPercussion-analysis-toolkit.md`.

The 4 ms contact-only stage remains diagnostic because individual stochastic
waveforms are not phase-comparable. Every later stage that can alter a saved
render now requires both relative improvement and absolute-quality acceptance.
Failed searches may be stored as experiments, but factory or reference starts
must never import them silently.
The offline fitter therefore writes `promotion_status` and an explicit
`listening_approval: false`; even a candidate that passes numerical gates is
not labelled a calibration or promoted by the tool.

Multi-resolution STFT spectral convergence and log-magnitude terms are useful
search objectives because different windows expose contact and sustained
structure at different time scales. Perceptually spaced ERB-band trajectories
then expose broad colour and decay without pretending that unrelated stochastic
phase should match. Neither family is a perceptual approval metric: direct,
level-controlled listening is the final gate. PEAQ is not used as the primary
loss because it was standardized for audible impairment of a reference signal,
not for judging whether two independently generated stochastic percussion
events describe the same instrument.

Relevant methodological references are the multi-resolution STFT loss in
[Parallel WaveGAN](https://arxiv.org/abs/1910.11480), its perceptually weighted
extension in [Yamamoto et al.](https://arxiv.org/abs/2101.07412), and the scope
of [ITU-R BS.1387 (PEAQ)](https://www.itu.int/rec/R-REC-BS.1387/).

A fresh search after limiting editable anchor centres to 15 kHz was also kept
as a rejected experiment rather than promoted automatically. It improved the
search's final causal-prefix objective, but failed the independent acceptance
gate: contact/bloom/early ERB shape errors remained roughly 5.6/6.5/8.7 dB,
fine-spectrum error was 6.7 dB, and the persistent-mode presentation-level
error was 14.1 dB. This is exactly the boundary between an optimizer seed and
an accepted listening fit.

## Reference object

The first object is the consistent 18-inch crash grid described by
`data/crash-calibration/private-corpus-a-crash-v1.json`. The public repository
uses the neutral name **private corpus A**. Exact vendor, product, instrument,
and filesystem details belong in an ignored local companion manifest. Source
audio and recoverable derivatives remain local and are never redistributed.

The grid contains edge, bow-tip, bow-shank, bell-tip, and bell-shank strikes,
multiple velocity layers, and repeated hits. `tools/build_private_cymbal_sweep.py`
generates deterministic MIDI and onset metadata. Capture qualification checks
that requested articulations are genuinely distinct, repeats have not
collapsed, onsets are valid, and velocity ordering is credible before samples
are offered to the workbench.

Other licensed and open cymbal collections are validation sources, not cells
silently pooled into this one-object fit. They are introduced after one object
works across its own strike regions and velocities.

### Gong bootstrap fit

`tools/calibrate_gong.py` holds one deliberately sparse, listening-oriented
gong bootstrap against the local `Gong Dresden 03` reference. It uses 17 active
painted handles, one normalized excitation shelf with a movable centre, passive
upward transport, and only the fixed 40 Hz/15 kHz T60 boundary values; all six
interior decay knots remain disabled. The saved workbench preset is `gong-v1`.

After separating excitation from observation, the first nine measured low
ridges are 375, 551, 129, 305, 422, 246, 621, 727, and 879 Hz in both renders;
their relative levels agree within 0.25 dB. Full-render RMS is 0.08148 versus
0.08131 in the reference, with a negative model-level setting. The model's five
ERB-band peak times are approximately 0.01, 0.20, 0.73, 0.82, and 0.90 seconds,
versus 0.01, 0.07, 0.48, 0.87, and 0.96 seconds in the reference.

The reference's directly measured band-decay curve is about 6.25 seconds at
40 Hz and 3.43 seconds at 15 kHz. Those are observed envelopes, not the modal
loss values: passive upward transport is another reason energy leaves a source
band. Joint fitting therefore uses intrinsic endpoints of 12 and 1.1 seconds;
the cascade plus these two losses matches the 4-second five-band levels within
about 3.3 dB, without adding an interior T60 knot.

This is an auditionable bootstrap, not a listening-approved final fit. Its
known structural mismatch is the 300 Hz--3 kHz envelope: its two band peaks
still arrive roughly 0.13 and 0.25 seconds too late. Modal painting now affects
observation only and the two-point T60 is already reserved for true loss, so
neither may be abused to conceal this transport-kinetics limitation.

## Current synthesis graph

```text
contact direct -------------------------------------> contact observation --+
       |
       `-> body excitation -> projected force -> one stochastic modal field -+
                                      |                               |
                                      `-> intrinsic upward cascade    |
                                                         body observation
                                                                  +--> output
```

The experimental body is one stored state with a 512-state ceiling: zero to 32
paintable centre handles reserve one state each, and deterministic stochastic
sideband pairs are allocated from the remaining shared pool. Global turbulence
control transfers normalized excitation energy from centres to satellites and
increases ERB spread, phase diffusion, and passive local exchange. A paintable
per-anchor `0..2x` scaler lets selected ridges stay clean. All anchors share one
frequency-dependent T60 curve and one mono radiation/output path. The older
sparse bank, statistical cloud, and separate turbulent residual primitives are
retained as reusable DSP modules, but their crash graph, state, and controls are
disconnected. They cannot alter the active renderer.

Contact writes directly into that one stored modal state. Strike location and
velocity colour redistribute a normalized new-force vector without recolouring
energy already circulating. Painted bars form an independently normalized
observation-prominence vector; they do not change stored strike energy.
The visible body-excitation gain sets the nonlinear body's drive; the separate
body-observation gain cannot alter stored energy or cascade behaviour.
Bloom is a passive state-level cascade through a fixed half-octave transport
stencil interpolated onto the available higher packets;
it has no separate signal, delay, latch, output gain, or feedback T60. Its rate,
energy acceleration, and destination phase diffusion are explicit. Higher strike
energy accelerates the upward transfer, while the separately visible velocity
brightness control also couples new force more strongly into high modes.
The renderer is mono; stereo presentation belongs after the instrument state.
Passive mute can only remove stored or future energy, and a zero-strength event
is a strict no-op.

This graph is deliberately decomposable: each branch can be soloed, bypassed,
and compared. Its architecture may change when listening shows that a component
cannot span the reference behaviour. Keeping it is not a claim that all current
parameterizations are perceptually useful.

## Manual fitting workflow

The browser workbench will compile the same C++ renderer to WebAssembly and
will be a separate, default-off developer build target. A session proceeds as:

1. Select one qualified reference hit and preserve its source level and
   provenance.
2. Audition the reference and synthesis at one shared, safe monitor gain.
3. Adjust a small set of perceptual macro controls while watching synchronized
   waveform, spectrogram, difference, and branch-ablation views.
4. Exercise repeated strikes and a two-dimensional strike pad so the fit is
   not valid only for one isolated render.
5. Snapshot promising states with all controls, renderer/schema versions,
   reference identity and hash, source transform, sample rate, and notes.
6. Promote a snapshot to a named fit only after direct listening at normal and
   level-matched gain. Keep its parent so a later regression can be reversed.
7. Warm-start the neighbouring velocity or articulation, then inspect which
   values are object constants and which form meaningful control curves.

The first pass fits one sample at a time. It does not average unlike hits or
infer size/location curves prematurely. Once several independently saved fits
exist, common object parameters can be frozen and input-dependent mappings can
be designed from the observed trends.

## Minimal control policy

The workbench exposes perceptual macros by default, not every implementation
parameter. A macro is accepted only when its range is useful, its action is
reasonably monotonic, and its label predicts what is heard. Raw parameters,
branch solos, and curve editors remain in an advanced panel for diagnosis.
Advanced means disclosed on demand, not inaccessible: no active DSP parameter
is silently fixed by the focused UI. Redundant controls are removed from both
UI and renderer only after controlled ablation.

The initial macro groups are:

- impact: strength, hardness/contact width, tonal-to-noise contact balance;
- object: 0-32 centre frequencies/levels up to 15 kHz, global turbulence, per-anchor
  turbulence response, packet spread, phase bandwidth, passive neighbour
  exchange, and broad body colour;
- evolution: intrinsic bloom rate/energy response and the shared decay shape;
- strike projection: radial location and implement blend; and
- constraint/presentation: mute, master level, and safe level matching.

This list is provisional. Finite-difference and sweep renders can reveal dead,
redundant, discontinuous, or dangerously sensitive controls, but listening
decides whether the parameterization is intelligible.

## Visual and numerical diagnostics

The workbench must show the unaltered reference and current synthesis with
explicit analysis settings. Required views include synchronized waveform and
high-resolution STFT, reference/synthesis split or mirror layouts, signed dB
difference, ERB-band energy trajectories, and branch solos. Zooming changes the
displayed region without silently changing the underlying audio or analysis
definition. Window, hop, FFT size, scale, floor, and colour range remain visible
and reproducible.

Diagnostics may flag clipping, onset error, level mismatch, missing bandwidth,
wrong decay, excessive tonality, or unstable repeated-hit accumulation. They
must not collapse those independent failures into a green “fit passed” label.

## Acceptance boundary

A first-object fit is credible only when:

- reference and synthesis remain distinguishable but plausibly describe the
  same cymbal across the fitted velocities and strike regions;
- attack, bloom, body, and tail each survive isolated listening;
- repeated hits accumulate and decay without instability or limiter pumping;
- the result works with the safety limiter bypassed in an offline diagnostic;
- saved parameters reproduce exactly with the recorded renderer version; and
- held-out references expose understandable model limitations rather than a
  memorized recording or microphone curve.

The always-on listening limiter is protection, never part of the fitted
instrument. Analysis uses the pre-limiter signal unless a view is explicitly
labelled otherwise.

### Modal decay policy

Ordinary fitting places modal frequencies and levels only. The active
40 Hz/15 kHz modal-T60 envelope is the sole decay model; no resolved-mode decay
parameter exists. Ridge-specific decay values may be measured as diagnostics,
but adding such a control requires a future explicit model extension after the
global T60 acceptance gate has passed. They must never be inherited from
another instrument calibration.
