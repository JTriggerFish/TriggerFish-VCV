# Percussion DSP implementation status

This note records what the clean-slate component library implements and what
its tests currently prove. It is an engineering status report, not a claim
that the ride model is calibrated or ready for listening.

## Implemented foundations

Reusable code lives in `src/tfdsp/percussion/`; no component is hidden in an
instrument or Rack module.

| Family | Components | Current analytic coverage |
| --- | --- | --- |
| Contact | half-sine force pulse, tonal chirp, enveloped noise with neutral-white and shelf-tilt settings, stochastic micro-contacts, explicit direct/body router | duration, endpoints, strength-to-energy mapping, deterministic seeds, bounds, routing |
| Delay and diffusion | 12-tap/2048-phase moving sinc delay, static Thiran delay, shared cubic Lagrange delay, fractional Schroeder all-pass | tone gain, low-frequency delay, polynomial exactness, integer-boundary continuity, impulse energy, five sample rates |
| Spectral motion | 255-tap antisymmetric FIR Hilbert transformer, phase-continuous signed SSB frequency shifter, and fourth-order translation-band guards | wanted level, image rejection, exact zero shift, through-zero automation, DC/Nyquist translated-content rejection, five sample rates |
| Resonance and loss | complementary three-band T60 loss, orthogonal Givens mixer, wet-only coupled fractional-comb network | exact identity at zero coupling, scattering energy, T60 recurrence, dry isolation |
| Cymbal bloom | slow stochastic delay, 2x oversampled bounded self-phase delay with a 4x reference implementation, serial two-all-pass dispersion loop with explicit outer feedback | 1x/2x comparison against 4x, causal onset and declared recurrence delay, zero-drive linear reduction, no hidden feedback sample, long contractive stress |
| Radiation | guarded TDF2 biquad designs and a static high-pass/colour/low-pass chain | centre gain, pass/reject bands, sample-rate sweep, non-finite recovery |

All sample loops are allocation-free. Tables and delay storage are created in
`Prepare`; filter and recurrence coefficients are built during preparation,
static configuration, or hit creation. The active chirp evaluates one sine and
the self-phase delay evaluates its intentional nonlinear map at 2x per sample.
The reverb's velvet mixer advances its smoothed rotation algebraically, without
per-sample trigonometric calls. The
release micro-benchmark command, `dev.ps1 benchmark-percussion`, reports
component cost through the MinGW launcher without imposing machine-dependent
CI limits.

## Deliberate distinctions

`EnvelopedNoiseBurst` and `MicroContactBurst` are separate. The first is a
compact oscillator-plus-noise source for tutorial-style cymbal/hi-hat graphs;
zero tilt is exactly white before its envelope. The second represents many
small implement contacts and is appropriate for mallets and brushes. Neither
is mixed into an instrument output implicitly.

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

- add smoothed, state-preserving live retuning for size, tune, location, mute
  and other controls; current delay-network retuning is explicitly static and
  clears state;
- extend the existing five-rate, non-finite, denormal, retrigger and automation
  tests to rapid simultaneous control changes across the assembled graph;
- produce isolated ablation WAVs and then assemble the first neutral ride graph;
- fit and validate that graph against the frozen real ride dataset.

Calibration must not begin by optimizing around any of these open numerical or
control-path questions.
