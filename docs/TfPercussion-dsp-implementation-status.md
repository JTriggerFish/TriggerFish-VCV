# Percussion DSP implementation status

This note records what the clean-slate component library implements and what
its tests currently prove. It is an engineering status report, not a claim
that the ride model is calibrated or ready for listening.

## Implemented foundations

Reusable code lives in `src/tfdsp/percussion/`; no component is hidden in an
instrument or Rack module.

| Family | Components | Current analytic coverage |
| --- | --- | --- |
| Control | fixed-capacity linear/geometric breakpoint trajectory, asymmetric passive-loss controller | exact endpoints, continuous retrigger, sanitization, fast loss/slow release, ordered attenuation-only gains |
| Contact and compact body | half-sine force pulse, tonal chirp, enveloped noise, finite micro-contact burst, gated/finite renewal micro-contact process, 2x correlated FM burst with 4x reference, explicit direct/body router | duration, endpoints, strength-to-energy mapping, deterministic seeds, gate release, FM spectral convergence, bounds, routing |
| Delay and diffusion | 12-tap/2048-phase moving sinc delay, static Thiran delay, shared cubic Lagrange delay, fractional Schroeder all-pass | tone gain, low-frequency delay, polynomial exactness, integer-boundary continuity, impulse energy, five sample rates |
| Spectral motion | 255-tap antisymmetric FIR Hilbert transformer, phase-continuous signed SSB frequency shifter, and fourth-order translation-band guards | wanted level, image rejection, exact zero shift, through-zero automation, DC/Nyquist translated-content rejection, five sample rates |
| Resonance and loss | complementary three-band T60 loss, passive constraint gains, orthogonal Givens mixer, projected wet-only fractional-comb network, explicit output submix | exact identity at zero coupling, scattering energy, T60 recurrence, passive attenuation, excitation isolation, group routing |
| Cymbal bloom | slow stochastic delay, 2x oversampled bounded self-phase delay with a 4x reference implementation, serial two-all-pass dispersion loop with explicit outer feedback | 1x/2x comparison against 4x, causal onset and declared recurrence delay, zero-drive linear reduction, no hidden feedback sample, long contractive stress |
| Radiation | guarded TDF2 biquad designs and a static high-pass/colour/low-pass chain | centre gain, pass/reject bands, sample-rate sweep, non-finite recovery |
| Observation | zero-capable static fractional delay and per-source gain, polarity, delay, and radiation paths | exact zero and integer delay, source isolation, polarity, non-finite recovery |

All sample loops are allocation-free. Tables and delay storage are created in
`Prepare`; filter and recurrence coefficients are built during preparation,
static configuration, hit creation, or quantized/smoothed frequency-shift
boundary updates. The active chirp evaluates one sine and
the self-phase delay evaluates its intentional nonlinear map at 2x per sample.
The reverb's velvet mixer advances its smoothed rotation algebraically, without
per-sample trigonometric calls. The
release micro-benchmark command, `dev.ps1 benchmark-percussion`, reports
component cost through the MinGW launcher without imposing machine-dependent
CI limits.

Rack calls `rack::system::resetFpuFlags()` on the main engine thread and every
engine worker. On x64 this enables FTZ and DAZ; on ARM64 it enables FTZ, with
DAZ supplied by the architecture. These flags are thread-local, so setting them
once during plugin initialization would not protect worker threads. Reusable
DSP still flushes non-finite and subnormal values at recursive filter states,
feedback delay writes, and public audio boundaries. This keeps the standalone
tests and non-Rack consumers safe without changing a caller's floating-point
environment from inside a component.

## Deliberate distinctions

`EnvelopedNoiseBurst`, `MicroContactBurst`, and `MicroContactProcess` are
separate. The first is a
compact oscillator-plus-noise source for tutorial-style cymbal/hi-hat graphs;
zero tilt is exactly white before its envelope. The second represents many
small implement contacts in one finite hit. The process adds finite clusters
and gated streams driven by a Poisson scheduler and smooth overlapping contact
windows. Density and incident amplitude remain separate controls. None is
mixed into an instrument output implicitly.

`CorrelatedFmBurst` is a general compact body/contact source rather than a kick
voice. Amplitude, carrier frequency, and frequency-deviation trajectories are
independent. The modulator continuously blends seeded band-limited irregular
motion with a periodic oscillator, and optional perturbation is explicit.

The dispersion return is stored in its base propagation delay. The signal is
read, passed serially through slow delay, two all-passes, self-phase delay and
loss, then written back with the new body drive. This avoids an accidental
extra sample at the outer feedback sum. The dispersion tap remains an analysis
and resonator-drive output, not audible dry leakage.

The coupled resonator is one superset: coupling zero gives independent combs;
nonzero coupling applies an orthogonal matrix before line loss. It is not the
room-reverb velvet FDN or its coefficient set.

`RadiationFilter` is an object/observation voicing stage. Its final biquad is a
two-pole low-pass, so the asymptotic electrical slope is 12 dB/octave. It is
not, by itself, a physical distance-dependent air-absorption model: atmospheric
loss belongs in an optional propagation/observation stage and must vary with
distance, humidity, and frequency when that distinction matters.

## Not yet accepted for instrument calibration

These tests are the first analytic quality level only. The following gates are
still open:

- add smoothed, state-preserving live retuning for size and tune; use dynamic
  excitation projections for location and integrate the implemented passive
  constraint controller for mute without clearing stored state;
- extend the existing five-rate, non-finite, denormal, retrigger and automation
  tests to rapid simultaneous control changes across the assembled graph;
- produce isolated ablation WAVs and then assemble the first neutral ride graph;
- fit and validate that graph against the frozen real ride dataset.

Calibration must not begin by optimizing around any of these open numerical or
control-path questions.

The Python numerical-analysis implementation and its remaining real-data gates
are recorded in `TfPercussion-analysis-toolkit.md`. No Plotly report or browser
renderer is part of the present component pass.

## Deferred structured-model extensions

`EnergyCoupler`, `AcousticCavity`, and `DistributedContactCoupler` remain
documented future components rather than prerequisites for the compact
video-derived instruments. Existing orthogonal mixing covers passive exchange
inside one resonator bank. Cross-body energy coupling waits for a selected
membrane/body model with energy-normalized ports; distributed collision waits
until compact hi-hat or snare fits demonstrate a state-dependent interaction
that passive loss and driven stochastic contact cannot reproduce.
