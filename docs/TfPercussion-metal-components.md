# Metallic percussion DSP components

This note defines the reusable metallic-percussion components and distinguishes
the active cymbal recipe from retained experiments. The complete current body
is documented in
[the nonlinear modal resonator architecture](TfPercussion-nonlinear-resonator-architecture.md).

## Active cymbal components

The metallic-plate recipe contains four nodes:

```text
contact -> unified modal field -> observation -> mono output
    `----------------------------> observation
```

The three pre-output connections are Boolean route switches. All gains are
named node parameters visible through the same C++, JSON, and workbench
descriptor surface. The recipe contains no hidden edge coefficients.

### Contact exciter

The contact model combines a half-sine force pulse, damped chirp, enveloped
tilted noise, and dense micro-contact gesture. Direct and body projections are
explicit. Implement and strength derive hit-local coefficients but do not add
saved parameters or recolour old body state.

### Unified stochastic modal field

Twenty-four paintable anchors each prepare one coherent centre and sixteen
satellites. One field therefore spans defined ridges, overlapping metallic
texture, and dense wash. Increasing turbulence redistributes normalized input
energy into satellites, expands ERB spread, adds magnitude-preserving phase
diffusion, and enables local orthogonal exchange.

All states share one paintable frequency/T60 curve and one mono body
observation. The former separate sparse bank, modal cloud, and filtered-noise
turbulent residual are reusable library primitives only; they are disconnected
from the active cymbal renderer.

### Modal energy cascade

The current bloom mechanism is state-level upward transport. It moves a bounded
fraction of stored energy from each lower packet to its adjacent higher packet,
preserving total energy before declared modal damping. The pairs are processed
high-to-low, so energy crosses at most one boundary per sample.

Three controls are sufficient:

- rate in octaves per second;
- dependence of rate on lower-packet stored energy; and
- magnitude-preserving destination phase diffusion.

This creates progressive spectral travel rather than a delayed high-frequency
burst. Higher-energy strikes move farther upward. The renderer separately
provides visible velocity-dependent excitation brightness, so the two effects
can be tested and fitted independently.

### Passive constraints and observation

Mute adds smoothed multiplicative loss to stored modal state. It never drives a
silent body. Observation provides independent direct/body levels and static
high-pass, broad colour, and low-pass filtering. It does not alter the body T60
or stored-energy accounting.

## Retained dispersion-loop experiment

`DispersionLoop` remains a tested reusable nonlinear transformation, but it is
not part of `CrashCymbal` or the active metallic patch. It contains one explicit
outer feedback loop:

```text
body drive -> base delay -> moving delay -> AP1 -> AP2 -> AP3 -> AP4
                 ^                                      |
                 `------ loss <- self-phase delay <-----+
```

The base and moving delays provide causal propagation. Each Schroeder allpass
has direct feedthrough when its coefficient is nonzero; its configured delay is
a group-delay scale, not a latency claim. The loop consequently tracks both:

- minimum causal propagation, used for recurrence reasoning; and
- nominal/group propagation, used to calibrate loop decay.

The self-phase stage uses bounded signed audio to modulate a fractional-delay
read. It runs at 2x in production and has a 4x test reference. Slow modulation
changes recurrence relationships independently of signal level. These remain
useful for other synthetic bodies, but the earlier idea of feeding a delayed
dispersion tap into the cymbal modal bank was rejected: it produced an
arbitrary late bloom instead of intrinsic upward energy flow.

Feedback components should not receive an automatic extra unit delay. A delay
network already has explicit state. A genuinely algebraic nonlinear loop must
instead be solved or redesigned explicitly.

## Other reusable resonator and texture components

### Modal bank and statistical cloud

`ModalBank` supplies directly placed two-pole modes with frequency, T60,
excitation projection, and output projection. `StatisticalModalCloud` prepares
a deterministic dense set from smooth distributions. They remain useful for
other recipes and offline decomposition, even though the cymbal combines their
roles in one field.

### Turbulent residual

`TurbulentResidual` is a three-band passive energy reservoir with a correlated
noise readout. It is a TriggerFish perceptual heuristic, not the wave-turbulence
solver of Cirio et al. It is disconnected from the cymbal because the unified
packet field plus state cascade currently offers a simpler, more controllable
representation of the same ridge-to-wash continuum.

### Coupled feedback resonator

The wet-only fractional-comb network is a coupled resonator/FDN superset:

```text
A = I                    independent feedback combs
A = orthogonal A(theta)  passive coupled feedback network
```

Its poles depend jointly on delay lengths, loss, and the feedback matrix. It is
not an orthogonal modal parameterization, and it is not the room-reverb velvet
FDN or its coefficient set. Use it only when an instrument needs a delay-body
signature that directly placed modal states cannot supply.

### Frequency shifting

Directly moving modal frequencies provides arbitrary inharmonicity for a modal
body. The single-sideband frequency shifter remains useful for translated
residual submixes and effects. It uses a latency-matched Hilbert pair and
phase-continuous signed oscillator. It should stay outside a feedback loop
unless repeated frequency translation is specifically intended.

Pitch shift, frequency shift, and delay modulation are distinct:

```text
pitch ratio r:      f -> r f
frequency offset d: f -> f + d
delay modulation:   time-varying phase and sidebands
```

None should be hidden behind a generic detune control.

## Fractional-delay requirements

Reusable delays declare their valid range, interpolation method, and latency.
Tests cover tone gain, group delay, integer-boundary continuity, modulation
sidebands, impulse response, recurrence T60, and supported sample rates.
Interpolation choice matters inside feedback: a convenient moving Lagrange
read may add frequency-dependent magnitude error, while a static Thiran allpass
preserves magnitude but approximates group delay.

## Component acceptance rules

Every stored-energy component must state:

- what enters and leaves the state;
- what can create, preserve, transfer, and remove energy;
- whether its controls act at preparation, hit time, or every sample;
- its causal latency and recurrence delay where applicable;
- behavior under retrigger and live constraints; and
- finite, denormal, sample-rate, and maximum-control behavior.

No component is accepted into an instrument merely because it sounds plausible
in isolation. It must improve reference listening and diagnostics without
creating redundant or hidden control dimensions.

## Deferred components

`DistributedContactCoupler`, cross-body `EnergyCoupler`, and a dedicated
`AcousticCavity` remain deferred. A cavity may initially be synthesized from a
modal/delay body and observation filters; a new primitive is justified only by
a behavior those components cannot express. Distributed plate contact is most
likely to become relevant when the hi-hat's two-body interaction is revisited.
