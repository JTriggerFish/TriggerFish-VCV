# Ride fitting methodology

## Status

This document specifies the next calibration system. The previous ride model,
renderer, objectives, fit scripts, and fitted parameters have been removed.
Nothing from that system is an accepted baseline or target implementation.

The first instrument is one ride cymbal because its reference corpus and
articulation metadata are already available. Hi-hat, snare, and kick work starts
only after the reusable ride components and calibration gates are credible.

## Preconditions

Calibration cannot begin until every component in the candidate graph:

- has a small, explicit API and declared units;
- passes the analytic, numerical-reference, automation, and stress tests in
  `TfPercussion-quality-and-calibration.md`;
- exposes deterministic intermediate renders;
- reports latency and valid sample-rate/control ranges; and
- can be bypassed for an ablation render.

A fitter must never compensate for an implementation defect.

## Frozen first target

The first target is the University of Iowa 21-inch ride selected by
`data/ride-calibration/iowa-21ride-mf-ff-v1.json`:

- bow MF and FF;
- bell MF and FF;
- one instrument, recording chain, and stick family;
- source levels preserved across the four cells.

Intermediate velocities, edge/shoulder, mallet, brush, tune/size, and other
cymbals are held out. They are introduced only after the four-cell graph passes.
Bitwig, Salamander, VCSL, manufacturer clips, and private Toontrack renders
remain independent validation or later generalization sources.

## Analysis representation

Each hit is represented by separate, overlapping perceptual regions. Nominal
boundaries are contact 0--15 ms, bloom 15--120 ms, early body 120--600 ms, and
tail after 600 ms. The stored manifest records data-derived boundaries for each
cell; the nominal values are not hard-coded fitting truths.

The target contains four families of measurements:

1. Contact: onset envelope, crest, direct spectral distribution, and energy.
2. Evolution: ERB-band energy trajectories, time-to-peak, flux, and upward
   spreading.
3. Persistent structure: stable low/mid damped components with uncertainty,
   salience, decay, and cross-velocity evidence.
4. Dense response: band decay, spectral envelope, occupancy, crest/flatness,
   modulation, and autocorrelation after masking accepted persistent structure.

No FFT peak list is called a modal target. A mode estimator must first pass
synthetic close-pair, damping, drift, coloured-residual, and noise-floor tests,
and its resynthesis must be auditioned.

## Objectives and acceptance

The calibration result is a vector of named losses and hard gates, not one
opaque score. At minimum it reports:

- contact/envelope loss;
- multiresolution ERB or log-spectral loss for each perceptual region;
- band-energy trajectory and band-decay losses;
- salience-weighted persistent-component frequency, level, and decay errors;
- dense-response tonality, occupancy, and modulation errors;
- bow/bell contrast and MF-to-FF level/brightness relationships; and
- pre-limiter peak, energy accounting, and stability diagnostics.

The modal, spectrotemporal, decay, relationship, and listening gates are joined
with AND. Improvement in wash cannot excuse wrong dominant resonances; matching
resonances cannot excuse a stiff or incorrectly evolving tail. Raw spectrogram
MSE and per-hit normalization are insufficient acceptance criteria.

## Fit order

Fit only parameters whose mechanism is active and observable:

1. Validate the contact, delay, allpass, nonlinear, shifter, resonator, loss,
   coupling, and radiation components independently.
2. Assemble the smallest ride graph and fit bow MF contact and broad response.
3. Add persistent resonator structure and frequency-dependent damping.
4. Add dispersion/bloom and dense response, with isolated branch renders.
5. Fit bow FF while preserving the MF object; velocity changes excitation and
   nonlinear operating point, not arbitrary object tuning.
6. Fit bell MF/FF through contact/radiation projections and only physically or
   perceptually justified regional structure.
7. Test held-out velocities and the other corpora before adding edge, implement,
   mute, tune/size, or additional cymbals.

Every stage produces reference-versus-current one-shots, fixed quarter-note
sequences, component solos, and ablations using one recorded gain convention.

## Optimization

Begin with bounded derivative-free searches while topology and discrete branch
choices are unsettled. Use parameter sweeps and sensitivity analysis to reject
unidentifiable controls before optimization. A differentiable PyTorch renderer
or surrogate is appropriate only when it reproduces the tested C++ graph and
the loss vector; gradients do not validate a bad objective.

Each run records the commit, dataset hashes, graph version, parameter schema,
seed, sample rate, renderer build, descriptor version, objectives, constraints,
and generated audio hashes. A result without this manifest is exploratory.

## Tool layout

The replacement tooling will be built as small modules with one responsibility:

```text
tools/percussion_fit/
  audio_io.py       loading, channel policy, resampling, hashes
  dataset.py        manifest validation and cell metadata
  segmentation.py   onset and perceptual-region boundaries
  descriptors.py    named observable measurements
  objectives.py     explicit loss terms and hard gates
  optimizer.py      search adapter only
  report.py         tables, plots, and manifests
  audition.py       deterministic A/B and sweep assembly
```

Instrument rendering remains separate from analysis. Individual files and
functions should stay focused; orchestration code may compose components but
must not contain their DSP or measurement algorithms.
