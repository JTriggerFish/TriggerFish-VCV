# Percussion analysis toolkit

## Scope

`python/triggerfish_percussion/` is the numerical layer used by fitting,
regression tests, and reporting. A first static Plotly A/B report and a small
local HTTP server are implemented. WebAssembly synthesis is not yet
implemented. All views consume the same aligned, level-preserving arrays rather
than recomputing a second interpretation of the audio.

This entire layer is advanced development tooling. Normal VCV Rack builds use
the Rack SDK Makefile and require none of Python, SciPy, Plotly, HTTP tooling,
or WebAssembly. The standalone CMake build also leaves Python bindings disabled
by default. Developers install the Python analysis dependencies explicitly;
future browser/WebAssembly targets must remain opt-in as well.

The current implementation provides:

- level-preserving integer and floating-point WAV I/O, explicit channel
  policy, polyphase resampling, and content hashes;
- high-frequency impact-onset detection, onset-only reference/synthesis
  alignment, optional waveform-copy correlation, and four explicit perceptual
  regions;
- canonical SciPy `ShortTimeFFT` transforms at several resolutions;
- an ERB-rate triangular energy partition whose band powers sum exactly to the
  selected FFT-bin power;
- pre-onset-noise-corrected Schroeder energy-decay curves with tail correction,
  explicit noise-limited/short-recording status, robust or least-squares fits,
  and per-ERB-band decay estimation;
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
- rejection of raw-waveform correlation for unrelated impact phases;
- calibrated STFT sine magnitude and common reference/synthesis axes;
- exact ERB partition energy conservation;
- transform-cache identity and version separation;
- known T60 recovery in two bands;
- explicit rejection of noise-limited and truncated T60 measurements;
- STFT region membership based on complete window support rather than frame
  centres alone;
- separation of two close damped modes with known frequencies, amplitudes, and
  decay constants;
- correct optimal assignment and unmatched-mode reporting; and
- zero distance for identical spectral, ERB, decay, relational, and modal
  descriptors.

These are synthetic numerical contracts, not proof that real cymbal modes have
been estimated correctly. Before fitting the first ride, the ESPRIT estimator
still needs adverse coloured-residual fixtures and validation on the DAFx
modal-estimation fixtures described below. Real modes must meet minimum pass
and repeated-hit evidence thresholds; an MDL minimum alone cannot make one a
target.

## Canonical analysis and fitting method

Use a hybrid analysis rather than one universal ridge detector:

1. Model the initial contact separately.
2. Estimate resolvable low- and mid-frequency components as exponentially
   damped sinusoids, including close pairs that produce beating.
3. Treat the dense upper spectrum statistically. Compare its spectral
   envelope, density, modulation, and decay on a perceptual filterbank rather
   than interpreting every FFT peak as a mode.
4. Validate both representations with resynthesis and held-out instruments. A
   fit is not accepted merely because it lowers one spectrogram error.

This split reflects the high modal overlap and nonlinear wide-band energy
cascade measured in cymbals. The modal/statistical crossover must be inferred
from overlap and reconstruction residual, not fixed globally; 3--8 kHz is only
a starting search range.

### Modal-estimation references

The implemented baseline follows sub-band ESPRIT. Sirdey et al., [Modal
analysis of impact sounds with ESPRIT in Gabor
transforms](https://www.dafx.de/paper-archive/2011/Papers/61_e.pdf), estimate
exponentially damped sinusoids in Gabor sub-bands. Ege, Boutillon, and David's
[high-resolution modal analysis](https://arxiv.org/abs/0909.0885) combines
sub-band conditioning, ESTER model-order selection, and ESPRIT for high modal
overlap.

[SAMPLE and BeatsDROP](https://github.com/LIMUNIMI/SAMPLE) remain useful
MIT-licensed independent baselines, particularly for interrupted tracks and
close pairs. They are corroboration, not ground truth. Reassigned STFTs,
[ssqueezepy](https://github.com/OverLordGoldDragon/ssqueezepy), and
[loristrck](https://github.com/gesellkammer/loristrck) are diagnostic ridge or
deterministic/noise views; they do not replace high-resolution parameter
identification. The GPL status of Loris/loristrck also keeps it out of the
shipped plugin.

For resolvable modes, run overlapping sub-bands at several time/frequency
resolutions. Merge duplicates in log-frequency, retain reconstructed energy,
and require stability across adjacent passes and repeated hits or microphones.
Store frequency, decay, amplitude, phase where meaningful, uncertainty, and
evidence count. After selecting structural modes, refit complex amplitudes and
phases against the full signal before resynthesis and residual analysis.

### Pre-processing and dense residual

- Preserve source level unless a named experiment explicitly gain-matches.
- Align from the high-frequency energy derivative and retain a pre-onset
  noise-floor region; do not align on the waveform peak. Reference/synthesis
  pairs align by their independently detected onsets. Raw-waveform correlation
  is reserved for explicitly identified copies or repeated recordings.
- Analyze close, overhead, and room microphones independently. Concatenating
  channels can invent modes.
- Keep contact, bloom, early-body, and tail measurements separate. Select exact
  boundaries from the recording rather than assuming nominal boundaries are
  physically exact. An STFT frame belongs to a region only when its declared
  window support satisfies that region's overlap threshold.
- Estimate band noise from the retained pre-onset frames, integrate only to the
  declared noise limit, and report noise-limited or recording-too-short status
  instead of extrapolating an unsupported T60 target.
- Subtract only the accepted deterministic resynthesis. On the residual,
  measure ERB-band energy and robust decay, spectral centroid/bandwidth,
  flatness and crest, temporal modulation, autocorrelation, late diffuseness,
  and velocity-dependent energy transfer.

### Comparison and acceptance gates

Match modal sets one-to-one in log-frequency with Hungarian assignment, then
report cents, log-decay, amplitude, and unmatched-mode errors. Compare the
dense residual with multiresolution log-STFT and ERB losses, compare band decay
curves directly, and compare onset/contact descriptors separately. Report
distributions across repeated hits rather than only one aggregate.

The acceptance rule is an AND: salience-weighted modal structure,
broadband spectrotemporal evolution, and matched listening must all pass.
Improved flatness or broad-band T20 cannot compensate for wrong dominant
ridges, and matching sparse modes cannot compensate for the wrong wash. Keep
absolute velocity-to-level curves for dynamics experiments; use one small
gain-only alignment for explicitly timbral comparisons, never independent
time-frequency normalization.

### Dataset and challenge validation

The fitting corpus includes Bitwig, Iowa, Salamander, VCSL, and controlled
Toontrack renders that are in scope. Fit trends on multiple sources and hold
out entire cymbals, rotating the held-out instruments, so microphone or library
processing cannot masquerade as a physical trend. Separate cymbal-level
structure, strike-location amplitudes/contact spectrum, velocity-dependent
contact and nonlinear transfer, and implement contact duration/spectrum.

Before real estimates become targets, validate on synthetic impacts spanning
close-pair spacing, decay, SNR, coloured wash, drift, contact duration, and room
tail. Score precision/recall, cents, decay, amplitude, and reconstruction. The
[DAFx Parameter Estimation Challenge](https://github.com/LOGUNIVPM/1st-DAFx-Challenge)
provides an external metal-plate generator and evaluation procedure.

### Nonlinear wash research boundary

Prior static 64-line VFM and 4,096-mode upper-field experiments could improve
individual density metrics but did not jointly reproduce flatness, attack,
phase predictability, band envelopes, and upward bloom. Static FDN coupling
also made persistent 2--8 kHz poles more exposed. These rejected experiments
are evidence against equating more linear modes with realistic wash.

Nonlinear thin-shell work by Ducceschi and Touze, Cirio et al., Skare and Abel,
and Conan et al. instead motivates energy-dependent upward transfer and phase
decorrelation. Any future implementation must expose transfer and intermediate
audio, return toward the stable body at low energy, and be evaluated as a
hypothesis rather than accepted from a single aggregate score.

## Report status and future interaction

The served Plotly report currently includes reference and synthesized audio
with one shared non-clipping audition gain, aligned waveforms, and common-scale
STFT heatmaps. It deliberately contains no old-model comparisons. ERB heatmaps,
signed differences, band decays, modal matches, descriptor/loss tables, and
component ablations are the next report views. Plotly remains the primary
library so time and frequency regions can be inspected interactively.

When the first instrument renderer exists, its actual C++ implementation will
also be compiled to WebAssembly. Sliders and a clickable strike surface will
therefore run the same DSP as the Rack module. Static deterministic renders
remain the fallback and regression format.
