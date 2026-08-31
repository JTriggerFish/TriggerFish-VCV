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
| Delay and diffusion | 12-tap/2048-phase moving sinc delay, static Thiran delay, fractional Schroeder all-pass | tone gain, low-frequency delay, integer-boundary continuity, impulse energy |
| Spectral motion | 255-tap antisymmetric FIR Hilbert transformer and phase-continuous signed SSB frequency shifter | wanted level, image rejection, signed/zero shift, through-zero automation, five sample rates |
| Resonance and loss | complementary three-band T60 loss, orthogonal Givens mixer, wet-only coupled fractional-comb network | exact identity at zero coupling, scattering energy, T60 recurrence, dry isolation |
| Cymbal bloom | slow stochastic delay, bounded signed-audio self-phase delay, serial two-all-pass dispersion loop with explicit outer feedback | zero-drive linear reduction, declared recurrence delay, no hidden feedback sample, long finite stress |
| Radiation | guarded TDF2 biquad designs and a static high-pass/colour/low-pass chain | centre gain, pass/reject bands, sample-rate sweep, non-finite recovery |

All sample loops are allocation-free. Tables and delay storage are created in
`Prepare`; filter and recurrence coefficients are built during preparation,
static configuration, or hit creation. The active chirp evaluates one sine and
self-phase delay evaluates its intentional nonlinear map per sample. The
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

## Not yet accepted for instrument calibration

These tests are the first analytic quality level only. The following gates are
still open:

- compare fractional delay, SSB shift and self-phase modulation with expensive
  offline/high-rate numerical references;
- measure frequency-shifter DC/Nyquist transition bands and translated-content
  rejection explicitly;
- choose and validate 2x production oversampling for self-phase distortion
  against the documented 4x offline reference;
- add smoothed, state-preserving live retuning for size, tune, location, mute
  and other controls; current delay-network retuning is explicitly static and
  clears state;
- stress every component across rapid simultaneous automation, retriggers,
  denormals and sample-rate changes;
- produce isolated ablation WAVs and then assemble the first neutral ride graph;
- fit and validate that graph against the frozen real ride dataset.

Calibration must not begin by optimizing around any of these open numerical or
control-path questions.
