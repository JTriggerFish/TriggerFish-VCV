# Percussion DSP quality and calibration

This document defines the acceptance contract for reusable DSP components,
instrument integration, rendering, and calibration.

## DSP quality contract

Every reusable primitive and component must pass a quality gate before an
instrument is calibrated around it. *Artifact-free* means free of avoidable
numerical, interpolation, control, and implementation artifacts over the
declared operating range. It does not mean removing intentional FM sidebands,
nonlinear bloom, comb coloration, collisions, or other parts of the model.

### Common invariants

Every component must:

- remain finite for every documented control value and supported sample rate;
- have deterministic reset and fixed-seed behaviour;
- produce exact silence from silent state and silent input;
- avoid denormal slowdowns and audio-thread allocation;
- state its latency, safe frequency range, and minimum delay explicitly;
- preserve gain, energy, or passivity where its contract claims to do so;
- separate deliberate loss or drive from accidental numerical attenuation;
- expose an unprocessed or reference mode where that is meaningful;
- survive rapid control movement without clicks, zipper noise, stale state, or
  unintended pitch jumps;
- crossfade or transform stored state safely when a live topology change cannot
  be interpreted continuously.

An output limiter is a last safety boundary, not proof of internal stability.
Tests must inspect pre-limiter states and fail if a supposedly passive or
contractive component relies on clipping to remain bounded.

### Sample-rate and reference testing

Static and automated renders are compared at 44.1, 48, 88.2, 96, and 192 kHz.
After resampling and compensating declared latency, equivalent physical
parameters should retain their frequency, decay, envelope, and level within
component-specific tolerances. A higher-rate or longer-kernel offline renderer
serves as the numerical reference for fractional delay, modulation, FM/PM, and
nonlinear processing.

Thresholds must be derived from audibility and the component's intended use,
then stored beside the tests. They must not be loosened merely because a later
instrument-level loss hides the error.

### Delay and feedback quality

Fractional-delay tests cover integer and non-integer positions throughout the
declared range, including continuous sweeps across every integer boundary.
Measurements include impulse position, magnitude response, phase/group-delay
error, interpolation noise, modulation sidebands, and feedback-loop pole and
RT60 error. Short high-feedback tests are distinct from room-scale delay tests.

Every feedback component is exercised with impulses, sustained bounded input,
silence after excitation, maximum decay, maximum modulation, and adversarial
automation. Orthogonal scattering is tested for energy preservation before
loss filters; complete loops are tested for the declared contraction and
frequency-dependent decay. No hidden interpolation filter may silently replace
the requested damping law.

### Nonlinear and modulation quality

FM/PM, self-phase delay, waveshaping, collision, and saturation paths are
compared against an oversampled offline reference. Tests measure aliased energy,
DC generation, intermodulation, level dependence, and recovery to the linear
or silent state. Oversampling and modulation-band limiting are selected per
component rather than applied blindly to the whole instrument.

Zero drive must reduce exactly to the documented linear topology. Increasing
drive must change only declared behaviours; for example, phase-distortion drive
may increase spectral complexity but must not introduce an unexplained output
gain jump or change the centre delay.

### Controls, events, and stochastic processes

Knob and CV paths share one documented mapping and smoothing law. Tests sweep
every control slowly and rapidly, exercise simultaneous control changes, and
verify trigger, gate, retrigger, and fractional-event timing. Mute, pedal gap,
and other nominally passive controls are checked for injected energy.

Random or granular contact processes use deterministic seeds and bounded event
energy. Density transitions must not expose individual clicks unless discrete
collisions are the requested regime. Repeated hits must vary in permitted
contact detail without randomizing stable object tuning or creating unbounded
loudness variation.

### Required validation levels

A component is not production-ready until it passes all applicable levels:

1. **Analytic unit tests:** known frequency, delay, gain, energy, and decay
   relationships.
2. **Numerical reference tests:** comparison with a deliberately expensive
   offline implementation.
3. **Automation and stress tests:** parameter extremes, retriggers, accumulated
   state, sample-rate changes, and long renders.
4. **Integration ablations:** component bypassed, isolated, and inserted into
   each relevant feedback position.
5. **Perceptual regression renders:** fixed listening files which make sweeps,
   transitions, tails, and failures easy to hear.

Calibration begins only after the active topology passes this contract. A fit
must never compensate for a known implementation artifact.

## Component and integration tests

Unit tests establish DSP contracts rather than perceived realism:

- contact duration, energy, spectrum, and monotonic velocity response;
- allpass magnitude and group delay;
- fractional-delay accuracy, continuity, and feedback stability;
- self-phase distortion reducing exactly to a constant delay at zero drive;
- bounded nonlinear excursion and controlled aliasing;
- predictable feedback-to-RT60 mapping;
- exact wet-only resonator routing;
- energy preservation of the coupling matrix;
- membrane mode ratios, split-pair behaviour, and state-dependent tuning;
- passive bidirectional transfer through cavities and other energy couplers;
- one-sided distributed contact with no force outside compression;
- negligible energy injection from slow contact-gap changes;
- observation-path delay, polarity, EQ, and source isolation;
- zero direct leakage from the dispersion tap;
- passive mute and constraint changes.

Deterministic render tests then expose every intermediate component and the
complete instrument. Fixed-seed WAV suites should contain even quarter-note
hits and slow, one-parameter sweeps rather than musical phrases that obscure
the comparison.

## Calibration system

Calibration operates on a complete, hashed, onset-aligned reference manifest.
It must preserve absolute level relationships across velocities and locations;
per-hit normalization can be used for a diagnostic but not for accepting a
model.

A percussion hit is divided into four perceptual regions whose losses remain
separate:

| Region | Primary constraints |
| --- | --- |
| Contact, approximately 0--15 ms | waveform/envelope, direct spectrum, crest |
| Build or bloom, approximately 15--120 ms | band time-to-peak, flux, spreading |
| Early body, approximately 120--600 ms | resonator balance and decay slopes |
| Tail, approximately 600 ms onward | band T20, density, tonality, modulation |

The exact boundaries are estimated from the reference. Useful objectives
include multiresolution ERB/log-spectrogram distance, band-energy trajectories,
spectral centroid and flux, bandwise energy-decay curves, peak prominence,
spectral crest, autocorrelation/cepstral periodicity, and tail occupancy. Modal
ridge estimates remain diagnostics; they cannot certify a perceptually wrong
sound by themselves.

Fitting is staged but ends with a joint validation because parameters cross
window boundaries:

1. one neutral instrument, one location, one medium velocity;
2. matched soft, medium, and hard velocities of that instrument;
3. additional locations of the same object;
4. held-out intermediate velocities and locations;
5. implement and constraint changes;
6. other sizes and constructions;
7. regression or interpolation between fitted instrument parameter sets.

Subsystem renders and ablations must accompany every fit. A candidate is
rejected if any explicit tuning, decay, density, transient, or level gate gets
worse beyond tolerance, even if an aggregate perceptual loss improves. Human
level-matched comparison against the same real sample is the final acceptance
step.

Derivative-free bounded optimization is the first choice while topology,
branch assignment, and feedback structure remain discrete or rough. A
differentiable PyTorch surrogate becomes useful only after the topology and
losses have survived listening tests. Multiple objectives should be retained
rather than concealed in one spectrogram-MSE score.
