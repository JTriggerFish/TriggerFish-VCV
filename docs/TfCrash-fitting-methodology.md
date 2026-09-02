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
Automated search code is retained as experimental developer tooling so useful
feature extraction and parameter-influence work is not lost, but it is not the
active fitting workflow.

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

## Current synthesis graph

```text
contact direct --------------------------------------> contact observation
       |
       `-> body drive -----> sparse arbitrary modal bank -> modal observation
                 |
                 `-> nonlinear dispersion -> dense modal cloud ---------+
                                          `-> turbulent energy residual -+-> residual observation
```

The sparse modal bank owns resolved persistent structure. The deterministic
modal cloud and optional turbulent residual own dense metallic wash. The raw
dispersion signal is an excitation/analysis tap, not an automatically audible
layer. The renderer is mono; any stereo presentation belongs after the physical
instrument state. Passive mute can only remove stored or future energy.

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

The initial macro groups are:

- impact: strength, hardness/contact width, tonal-to-noise contact balance;
- object: sparse pitch scale, sparse/dense balance, spectral colour;
- evolution: bloom amount/time and low/mid/high decay shape;
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
