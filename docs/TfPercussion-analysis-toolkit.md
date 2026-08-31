# Percussion analysis toolkit

## Scope

`python/triggerfish_percussion/` is the presentation-neutral numerical layer
used by fitting, regression tests, and the future interactive report. Plotly,
HTTP serving, and WebAssembly synthesis are deliberately not implemented yet.
They will consume these same arrays and named losses rather than recomputing a
second interpretation of the audio.

The current implementation provides:

- level-preserving integer and floating-point WAV I/O, explicit channel
  policy, polyphase resampling, and content hashes;
- high-frequency impact-onset detection, sample-only pair alignment, and four
  explicit perceptual regions;
- canonical SciPy `ShortTimeFFT` transforms at several resolutions;
- an ERB-rate triangular energy partition whose band powers sum exactly to the
  selected FFT-bin power;
- noise-corrected Schroeder energy-decay curves, robust or least-squares decay
  fits, and per-ERB-band decay estimation;
- contact energy, envelope crest, centroid, bandwidth, spectral flatness,
  crest, and flux trajectories;
- a damped-sinusoid ESPRIT baseline, complex-amplitude reconstruction, and
  Hungarian one-to-one mode matching;
- settled linear-phase FIR subbands, multiresolution/repeated-hit evidence
  aggregation, uncertainty summaries, full-signal amplitude/phase refitting,
  and deterministic modal resynthesis/residuals;
- singular-value and MDL model-order evidence, retained as a diagnostic rather
  than an automatic acceptance decision;
- named log-spectral, ERB-level, ERB-change, decay, modal, and cross-cell
  relationship distances; and
- a versioned transform cache plus a report-neutral reference/synthesis
  comparison object.

Losses and future plots share the `StftResult`, `ErbFilterbank`, and
`Comparison` data structures. A visual report must never apply different
normalization, alignment, masking, or transforms from the numerical objective.

## ERB representation

The current ERB transform is an energy-conserving aggregation of STFT power,
not a reconstructing gammatone filterbank. Frequency bins are linearly
interpolated between uniformly ERB-rate-spaced centres, and the weights for
each selected bin sum to one. It is appropriate for efficient fitting losses
and auditable band-energy trajectories.

A Hohmann reconstructing gammatone bank, for example through `pyfar`, remains
a useful independent auditory cross-check. It should not silently replace the
canonical loss transform: the two representations have different delay,
overlap, and energy conventions and must be named separately.

## What the tests prove

Synthetic tests currently establish:

- WAV full-scale conversion and floating-point level preservation;
- resampling duration and pass-band level;
- impact alignment without time warping;
- calibrated STFT sine magnitude and common reference/synthesis axes;
- exact ERB partition energy conservation;
- transform-cache identity and version separation;
- known T60 recovery in two bands;
- separation of two close damped modes with known frequencies, amplitudes, and
  decay constants;
- correct optimal assignment and unmatched-mode reporting; and
- zero distance for identical spectral, ERB, decay, relational, and modal
  descriptors.

These are synthetic numerical contracts, not proof that real cymbal modes have
been estimated correctly. Before fitting the first ride, the ESPRIT estimator
still needs adverse coloured-residual fixtures and validation on the DAFx
modal-estimation fixtures described in
`TfCymbal-analysis-methods.md`. Real modes must meet minimum pass and repeated-
hit evidence thresholds; an MDL minimum alone cannot make one a target.

## Future interactive report

The future served Plotly report will include reference and synthesized audio,
gain-preserved and explicitly gain-matched A/B modes, common-scale STFT and ERB
heatmaps, signed differences, band decays, modal matches, descriptors, and
component ablations. Plotly is the primary plotting library so time and
frequency regions can be inspected interactively.

When the first instrument renderer exists, its actual C++ implementation will
also be compiled to WebAssembly. Sliders and a clickable strike surface will
therefore run the same DSP as the Rack module. Static deterministic renders
remain the fallback and regression format.
