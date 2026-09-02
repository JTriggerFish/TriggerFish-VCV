# Percussion analysis toolkit

## Scope

`python/triggerfish_percussion/` is the numerical layer used by fitting,
regression tests, and reporting. A static Plotly A/B report remains available
for regression work. The primary listening tool is now a local browser
workbench backed by the same C++ graph through WebAssembly. All views consume
the same aligned, level-preserving arrays rather than silently peak-normalizing
or recomputing a second interpretation of the audio.

This entire layer is advanced development tooling. Normal VCV Rack builds use
the Rack SDK Makefile and require none of Python, SciPy, Plotly, HTTP tooling,
or WebAssembly. The standalone CMake build also leaves Python bindings disabled
by default. Developers install the Python analysis dependencies explicitly;
the browser/WebAssembly target is also explicit and opt-in.

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
- Preserve the same decomposition in synthesis: direct body drive excites the
  sparse modal branch, while dispersed drive excites the parallel dense
  residual. Do not force the residual through the modal bank or expose the raw
  dispersion tap to hide an error in either decomposition.

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

The fitting loss uses four complementary views rather than treating one
spectrogram distance as perception:

1. Multiresolution log-magnitude STFT/ERB trajectories retain attack timing
   and spectral evolution. This follows the multiresolution STFT loss used by
   [Parallel WaveGAN](https://arxiv.org/abs/1910.11480) and its
   [perceptually weighted extension](https://arxiv.org/abs/2101.07412), while
   using ERB pooling to avoid fitting random fine-bin fluctuations.
2. A normalized 96-band ERB spectrum constrains fine spectral shape, while
   short-time centroid and roll-off constrain brightness independently of
   total band energy.
3. Within-band flatness, crest, and a median-filter ridge/residual ratio reject
   tonal substitutes for noise-like wash. The ridge split is based on the
   time/frequency median filtering used in Fitzgerald's
   [HPSS method](http://dafx10.iem.at/papers/DerryFitzGerald_DAFx10_P15.pdf).
4. Explicit persistent-mode matching constrains stable resolved ridges. This
   is the deterministic half of the established sinusoidal-plus-residual
   decomposition represented by UPF's
   [SMS Tools](https://www.upf.edu/web/mtg/sms-tools).

For 30 ms and longer residual texture, cochlear-envelope modulation and
cross-band correlation are the next validation layer. McDermott and
Simoncelli's
[auditory texture model](https://mcdermottlab.mit.edu/bib2php/papers/McDermott_Simoncelli_2011_sound_texture_synthesis.pdf)
shows why band power and sparsity alone are insufficient. Joint time-frequency
scattering is a particularly relevant held-out metric: the
[perceptual-neural-physical sound-matching study](https://arxiv.org/abs/2301.02886)
evaluates it directly on nonstationary inharmonic percussion. Neither a
speech-trained learned metric nor exact waveform/phase error is a primary
gate for this stochastic synthesizer.

Schema 5 adds the brightness and texture terms to both optimizer penalties and
hard acceptance. The initial numerical limits are conservative rejection
thresholds, not claims of perceptual equivalence. Final tolerances must be
estimated from independent repeated recordings of the same articulation and
velocity, then validated by blinded A/B listening. A report cannot say that a
fit passes merely because one aggregate number fell.

### Dataset and challenge validation

The fitting corpus includes independent licensed libraries, Iowa, Salamander,
VCSL, and controlled private-corpus-A renders that are in scope. Fit trends on
multiple sources and hold
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

## Report status and interactive direction

The served Plotly report includes like-for-like reference and synthesized audio
with one shared non-clipping audition gain. That gain targets the mean dB peak
of the distinct source velocity layers, weighting each velocity once rather
than each round robin, and is capped against the loudest exported clip. It is
used only for report WAVs; analysis and stored calibration audio retain their
source levels. The report also provides aligned waveforms, separate
128-sample attack heatmaps, shared-scale full STFT heatmaps, signed spectral
differences, band decays, named acceptance gates, and the cumulative causal-fit
table. It deliberately contains no old-model comparisons or derived-component
pairs presented as references. ERB heatmaps, modal-match overlays, and component
ablations remain later report views. Plotly remains the primary library so time
and frequency regions can be inspected interactively.

The current static report is a regression fallback, not the main fitting UI.
The next tool compiles the actual C++ renderer to WebAssembly and adds live
perceptual controls, a clickable strike surface, synchronized spectrogram
comparison, and versioned snapshots. It is a separate, default-off developer
target and is never a dependency of Rack, package, or official release builds.
