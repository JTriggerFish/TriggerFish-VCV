# Metallic percussion DSP components

This document specifies the reusable cymbal and metallic-resonator components.
The overall design and code-organization rules are in
[the percussion architecture](TfPercussion-synthesis-architecture.md).

## Cymbal dispersion and bloom

The proposed cymbal core follows a perceptual excitation--dispersion--resonance
partition. A short contact enters a serial processor chain contained inside one
outer feedback loop:

```text
                          +-----------------------------------------------+
                          |                                               |
body drive -> sum -> base delay -> slow delay -> AP1 -> AP2 -> AP3 -> AP4 |
              ^                                           |                |
              |                                           v                |
              +-- loop loss/filter <- self-phase delay <------------------+
                                                        |
                                                        +-> dense residual
```

The base delay sets the principal circulation time. The serial allpasses spread
each circulation into several arrivals. The slow modulated delay supplies
signal-independent irregularity, preventing a static loop from settling into
obvious periodic combs. The self-phase delay supplies signal-dependent
nonlinearity and therefore velocity-dependent bloom. The dispersion signal is
available as an analysis tap but is not mixed directly into the instrument
output.

The outer return shown above is a real, essential feedback connection: every
circulation traverses the delays, allpasses, self-phase stage, and loop loss.
It is not shorthand for a feed-forward cascade or repeated offline processing.
Its explicit propagation delays make it causal without a separate hidden unit
delay at the summing junction.

All of these stages are serial from the perspective of the outer loop, although
each allpass contains its own internal feedback. Parallel dispersion branches
are a possible later extension, not part of the first identifiable model.

### Slow modulation versus self-phase modulation

The two moving delays have different control signals:

```text
d_slow[n] = D_slow + depth_slow * smooth_random[n]
d_self[n] = D_centre + depth_self * f(audio[n])
```

Slow modulation continually changes the recurrence pattern even in the quiet
tail. Self-phase modulation is strongest while the circulating signal is
strong and complex and becomes progressively more linear as it decays. Slow
modulation before self-phase distortion continually presents different phase
relationships to the nonlinear stage; placing it after would instead smear the
components already generated there.

Neither the slow delay nor this ordering is accepted on description alone.
Component and listening ablations must compare allpasses alone, each modulation
mechanism separately, both together, and both possible orders.

### Self-phase-delay implementation

The phrase *amplitude-controlled delay* means that signed instantaneous audio
controls phase; it is not ordinary amplitude modulation. A causal initial
implementation is:

```text
s[n]       = signal after the allpasses
e[n]       = amplitude follower(s[n])
control[n] = tone_lowpass(s[n])
control[n] = control[n] / (epsilon + e[n])       # optional normalization
offset[n]  = maximum_excursion * tanh(drive * control[n])
delay[n]   = centre_delay + offset[n]
y[n]       = fractional_delay.read(delay[n])
fractional_delay.push(s[n])
```

The control is signed audio, not an absolute-value envelope. The envelope is
used only for optional normalization. Partial rather than complete
normalization should remain available because complete normalization removes
useful velocity dependence. `Tone` limits the bandwidth of the modulation,
`tanh` bounds its excursion, and the centre delay must leave the interpolator a
safe causal margin in both directions.

This construction approximates self-phase modulation: a component
`sin(omega*n)` is read at a delay which depends on the signal itself, generating
correlated sidebands. It is an operational model, not a claim to reproduce the
private implementation of any commercial effect. The audio-rate nonlinear
stage should initially run at 2x, with a 4x offline reference used to assess
aliasing.

### Feedback causality and zero-delay loops

Do not insert a one-sample delay into every feedback path merely to make the
code causal. At 48 kHz one sample contributes 37.5 degrees at 5 kHz and 75
degrees at 10 kHz, enough to retune a short high-feedback metallic structure.
Choose the solution from the intended topology:

- Delay networks, combs, FDN lines, and the dispersion loop already contain an
  explicit propagation delay. Read the prior buffer state before writing the
  current input; no additional unit delay or Newton solve is required.
- Allpass filters use their defined internal state. If an allpass with direct
  feedthrough is enclosed by a loop that contains no other delay, solve its
  linear feedthrough algebraically or revise the topology explicitly.
- Analog-derived filters or instantaneous nonlinear feedback stages may create
  a genuine algebraic loop. Use a topology-preserving/zero-delay formulation;
  solve a linear loop in closed form and a nonlinear loop with a bounded Newton
  or equivalent robust solver, with convergence tests and a defined fallback.
- A state-controlled fractional delay is explicit when its modulation is
  derived from already available state and its read stays behind the write
  head. Making its delay depend implicitly on its own current output is a
  different topology and requires a solve.

The proposed dispersion chain begins with a base propagation delay and its
moving reads retain a declared causal margin. Consequently its outer loop is
not algebraic and should not acquire another one-sample delay. Zero-delay
solving is reserved for any later instantaneous nonlinear/filter substage whose
equations actually demand it. Component tests must measure the complete loop's
pole/recurrence shift so an accidental sample of latency cannot pass unnoticed.

## Sparse modes and dense residual

### Experimental unified modal field

The complete current instrument graph, including contact, bloom, modal-state
equations, energy accounting, control mapping, and observation, is documented
in the self-contained
[nonlinear resonator architecture](TfPercussion-nonlinear-resonator-architecture.md).
This section keeps
the reusable metallic-component rationale and legacy comparison.

The crash workbench can replace the three body renderers below with one
stochastic modal field. Twenty-four editable anchors, limited to a useful
40 Hz--15 kHz design range, each generate one exact centre mode and sixteen
symmetric ERB-domain satellites. Anchor energy is normalized
before it is divided between the centre and satellites, so increasing
turbulence does not implicitly increase body level. The same contact and bloom
forces excite one stored field and one radiation path.

Each mode can receive independent signed phase kicks. For a requested
decorrelation bandwidth `B`, the fixed rotation uses
`cos(theta) = exp(-pi * B / sampleRate)` and a randomly signed `sin(theta)`.
The rotation has unit magnitude: it broadens the expected spectral line without
changing modal energy or the separately declared T60. At zero bandwidth it is
exactly the ordinary deterministic modal recurrence.

Neighbouring states also undergo small randomly signed Givens rotations.
Alternating pairings exchange energy within a packet and across adjacent packet
boundaries. These rotations are orthogonal, so they cannot inject energy. The
main Turbulence macro transfers energy from each coherent centre to its
satellites and increases packet spread, phase bandwidth, and exchange. Every
anchor also has a paintable `0..2x` response scaler. Zero keeps that packet
coherent, one follows the global trajectory, and values above one reach its
diffuse endpoint earlier. Neighbour exchange uses the geometric mean of the two
local amounts, so a clean packet cannot be randomized through its neighbour.
The last three global controls remain visible as advanced workbench ablations
while the turbulence trajectory is calibrated.

The displayed bell around an anchor is a design metaphor for this stochastic
neighbourhood, not a claim that the cymbal contains a literal Gaussian family
of physical eigenmodes. This is a constructive perceptual approximation, not the wave-turbulence solver
of Cirio et al. Its acceptance criteria are smooth resolved-to-diffuse
interpolation, invariant level and T60, stable restrikes, useful location
projections, and a better match to reference ridge width and occupancy than the
legacy branches. The workbench retains the legacy path for direct A/B testing;
it is not part of the proposed final instrument topology.

The bloom loop separately exposes all-pass diffusion. Scaling every Schroeder
coefficient to zero leaves the matched pure delays in place, preserving causal
loop length while removing frequency-dependent group-delay dispersion. This is
kept independent from bloom level, development, and nonlinear self-phase.

### Legacy three-branch comparator

The following topology remains implemented only for workbench A/B diagnosis.
It is not the current proposed crash body.

A cymbal output is divided into resolvable persistent structure and unresolved
high-density response. These are parallel audible representations of the same
dispersed excitation, not a serial effects chain:

```text
body drive ----+-----------> arbitrary modal bank -> modal radiation --------+
               |                    ^                                        |
               `-> dispersion -----+-> statistical cloud -> wash radiation -+-> output
                                `----> three-band turbulent residual
                                                   -> turbulence radiation --+
```

The raw dispersion tap is not mixed into the output. An independent bloom-route
gain sends it to both modal representations and the turbulent residual; at zero
the loop may still evolve for analysis but cannot excite the body. Resolved modes additionally retain
a strong direct-body feed for the immediate coherent response; the dense cloud
uses only a weak direct bypass so its audible development follows the bloom.
Keeping separate radiation paths allows their spectral balance and spatial
presentation to be fitted without misusing modal gains as post-EQ. Turbulence
is not scaled by the resolved-to-wash balance; its own amount controls that
parallel source.

A local contact and the dispersed bloom enter each modal bank through separate
excitation projections but write the same modal state. The contact projection
carries the current strike's location and velocity colour. Bloom uses a fixed
body-wide projection. Consequently a retrigger adds force to stored energy
without changing how an earlier strike's circulating bloom excites the bank.
The self-phase stage likewise derives its excursion from the instantaneous loop
signal rather than a latched latest-strike velocity.

The sparse bank uses independently parameterized two-pole modes. Each mode has
an explicit frequency, T60/decay, excitation gain, and observation gain. This
gives direct inharmonic mode placement and a substantially better-conditioned
fit than using coupled comb delays as the primary pole model. Location and
implement select smooth excitation/observation projections over a shared mode
set; they do not create unrelated objects. Size, profile, and tune are
regularized macros over modal frequency, decay, and gain curves and need not
move every mode by the same ratio.

The dense residual is fitted only after accepted modal resynthesis is removed
from the reference. Its target is ERB-band evolution, occupancy, crest and
flatness, autocorrelation, temporal modulation, and band decay rather than
individual ridge frequency. The first renderer is a deterministic statistical
modal cloud. Its individual frequencies, phases, gains, and decay variations
come from a fixed object seed; the fitted controls are smooth frequency-density,
gain-envelope, and decay distributions. This preserves strike correlation and
stateful repeated-hit accumulation without presenting independent noise grains.
The current cloud uses a 33-knot ERB-domain object gain profile and a six-knot
ERB-domain decay envelope. Older gain grids are interpolated into the current
profile; the decay endpoints remain the named low/high decay controls and only
four internal decay offsets are free. Both curves are interpolated at
instrument preparation, so no per-sample table lookup or independently exposed
modal gain is added. The normalized gain profile is applied to modal excitation,
so it describes stored spectral energy. Per-path level and radiation remain a
separate observation stage. This factorization is transfer-equivalent for the
present independent linear modes but remains meaningful if coupling or
state-dependent loss is added later.

For strong, noise-like metallic responses, a second dense representation may
approximate unresolved nonlinear energy transfer. The turbulent residual splits
the dispersion signal into complementary low/mid/high bands, integrates squared
band energy into independently damped passive reservoirs, and reads those
reservoirs with bands from one deterministic correlated noise stream. Its T60s
describe macro energy loss, not contact-grain duration. Repeated strikes add
energy to the existing reservoirs, while mute applies the same per-sample modal
damping gains to their state. Zero band gains make it exactly inaudible. This
path is retained only because a private-corpus-A ablation improved envelope and flatness;
the current anchor still fails fine-spectrum acceptance, so it is not yet a
calibrated default.

This three-band reservoir is a TriggerFish perceptual heuristic, not an
implementation of a published wave-turbulence algorithm. The closest physical
reference, [Cirio et al. (2018)](https://www.cs.columbia.edu/cg/waveturb/),
detects chaos in a nonlinear low-frequency shell simulation and evolves a
time-varying spectral-energy distribution with a frequency-domain diffusion
model before adding the synthesized turbulent texture to the shell sound.
[Skare and Abel (2019)](https://www.dafx.de/paper-archive/2019/DAFx2019_paper_48.pdf)
instead retain a very large modal bank and approximate nonlinear energy transfer
by activating or coupling measured high-energy modes. Both keep structured modal
sound alongside the nonlinear contribution. A future residual replacement should
therefore test a multi-band energy cascade driving dense modes; the present
filtered-noise readout must not be described as either paper's method.

A coupled wet-only fractional-comb network remains an ablation if the modal
cloud lacks measured temporal correlation. Its branches are parallel when
their feedback matrix is the identity; replacing that identity with an
energy-preserving scattering matrix turns the same bank into a coupled feedback
delay network:

```text
A = I                   independent parallel combs
A = orthogonal A(theta) coupled resonator FDN
```

The implementation remains a useful FDN superset with continuously
controllable coupling. Products of Givens rotations introduce passive energy
exchange. However, its delay lengths, coupling, feedback, and output projection
jointly determine its poles; they are not called an orthogonal modal
parameterization. Coupling is retained only when ablation shows improved dense
evolution over an uncoupled residual with no exposed comb lattice or unstable
late-energy regrowth.

This is not the room-reverb FDN transplanted into a cymbal. The room network's
sixteen lines and delayed velvet scattering were optimized for a diffuse,
spectrally flat late field. A cymbal network must retain controllable pole
structure, frequency-dependent decay, and correlated metallic evolution. Some
reverb delay and scattering primitives are reusable, but its coefficient set
and complete topology are not a starting assumption.

Each residual line needs:

- a fractional delay or explicitly parameterized resonant frequency;
- feedback expressed through decay time rather than an arbitrary gain;
- frequency-dependent loss;
- wet-only output and a deterministic input/output projection;
- optional branch delay or allpass;
- optional frequency shifting of selected groups to break exact harmonic comb
  relationships, if an ablation justifies it.

Frequency shifting adds a fixed hertz offset and is distinct from pitch
shifting, which preserves harmonic ratios. Direct placement supplies all
inharmonicity required by the sparse modal bank. A shifter therefore remains
outside the feedback loop and is only an optional coloration of residual
submixes; it must beat the unshifted residual in matched listening and residual
metrics.

### Pitch shift, frequency shift, and delay modulation

These are three different operations and must not share one vaguely named
"detune" control:

```text
pitch shift by ratio r:       f_k -> r f_k
frequency shift by offset d:  f_k -> f_k + d
delay modulation:             time-varying phase and sidebands
```

A pitch shifter preserves the ratios among partials. It can transpose an
external source or provide an intentional effect, but it cannot turn the exact
harmonic series of a feedback comb into an inharmonic cymbal series. The
repository's existing `WindowedPitchShifter` is an eight-grain, fixed
one-octave-up shifter designed for a diffuse reverb-feedback branch. It is not
part of the neutral cymbal graph and is not a general implementation of the
tutorial's frequency-shifter stage.

A frequency shifter translates every positive-frequency component by a fixed
signed number of hertz, so harmonic ratios are broken. This is the candidate
operation after selected resonator groups. The first implementation should be
a single-sideband shifter:

```text
q[n] = Hilbert{x[n]}
phi[n+1] = phi[n] + 2 pi d[n] / sampleRate
up[n]   = delayed(x[n]) cos(phi[n]) - q[n] sin(phi[n])
down[n] = delayed(x[n]) cos(phi[n]) + q[n] sin(phi[n])
```

The real path must be delayed by exactly the Hilbert transformer's latency.
For the quality reference, use a linear-phase odd-symmetric FIR Hilbert
transformer with a declared valid band. A lower-latency polyphase-allpass
quadrature implementation may replace it only after nulling against that
reference and meeting image-rejection and phase-error limits. The oscillator
must be phase-continuous under CV, accept signed offsets for through-zero
shifting, and never reset on a parameter change. Pre/post band limiting must
discard translated content beyond DC or Nyquist rather than folding it back.

Apply this shifter to a wet resonator-group output before group summation and
stereo radiation. Do not shift the direct contact transient by default. Keep it
outside feedback during identification: placing it in a loop creates a
frequency-shifting feedback network in which every circulation moves again,
which is a separate, harder-to-calibrate effect.

Delay modulation in the dispersion loop is neither of the above. It produces
vibrato and correlated sidebands by continuously changing propagation delay;
its rate and depth describe motion rather than a fixed spectral remapping.

The instrument `tune` or `size` control is also neither of the above. It is a
macro over modal frequencies/delay lengths, damping, dispersion, coupling, and
radiation. Its fitted response may be nonuniform across mode groups. A global
audio-domain pitch shifter would preserve the wrong relationships and would
transpose contact noise and radiation artifacts along with the object.

The general component library may ultimately contain both shifter types, but
their initial roles are deliberately unequal:

- `SingleSidebandFrequencyShifter`: optional dense-residual coloration for
  breaking exact comb relationships;
- variable-ratio pitch shifter: optional source transformation or creative
  effect, not required by the neutral cymbal model;
- `WindowedPitchShifter`: retain specifically as the existing octave shimmer
  implementation until it receives a separate generalization and quality
  audit.

Frequency-shifter tests must include single tones over the declared band,
unwanted-sideband/image rejection, signed offsets, zero-offset equivalence to
the latency-matched input, through-zero sweeps, DC/Nyquist boundaries, impulse
latency and ringing, stereo coherence, and sample-rate equivalence. Pitch-
shifter tests must additionally verify ratio/cents accuracy on multitone input,
harmonic-ratio preservation, transient smearing, grain-boundary continuity,
live ratio sweeps, and upward-shift alias rejection.

## Fractional-delay requirements

The repository currently has three different designs, none of which is yet a
fully validated general short metallic delay:

- `LateReverb::FractionalDelay` accepts a floating-point delay and uses
  four-point Lagrange interpolation with a two-sample minimum. It was designed
  and tested for room-scale delay lines.
- `VelvetFeedbackMatrix::TransitionFractionalDelay` uses integer nominal taps;
  fractional reads occur during modulation.
- `MetalFeedbackField` uses a first-order Thiran allpass for short fractional
  delays and currently enforces a seven-sample configured minimum.

Lagrange interpolation is convenient for moving reads but its magnitude is not
exactly unity, so repeated high-feedback circulation changes the intended
frequency-dependent decay. A Thiran allpass preserves magnitude for a static
fractional delay, but its group-delay approximation is imperfect at high
frequency and naive time variation can produce artifacts at integer-delay
boundaries.

A reusable `MetallicFractionalDelay` should consequently provide a static
allpass mode and a smoothly modulated mode, state its safe minimum delay, and
be validated against a high-rate reference. Tests must cover fractional delay
and group-delay error, loop pole frequency, RT60, integer-boundary sweeps,
modulation sidebands, long high-feedback stability, and sample-rate equivalence.
Existing tests which merely show retained treble energy are insufficient to
claim artifact-free short-delay operation.

## Radiation and live losses

The audible instrument output contains direct contact radiation, sparse-modal
radiation, and dense-residual radiation, followed by object-specific
observation filtering and conservative output conditioning. The residual path
is a first-class model output, not leakage. The raw dispersion or stochastic
body-drive tap must not leak into the output merely to improve spectral
flatness.

Body synthesis remains mono. A mono observation is always available. Optional
stereo presentation applies two observation projections to the same direct,
modal, and residual taps; it may use small propagation delays and different
radiation views, but it does not duplicate or decorrelate the underlying body
state. Stereo width is therefore an output-presentation control, not an
instrument-size or modal-coupling parameter.

Location affects both excitation and radiation. For a ride, a bell strike has
stronger tonal contact and cup coupling, bow balances contact and plate bloom,
and edge contact is broader and favours high-mode wash. Velocity primarily
changes incident energy, excitation bandwidth, and nonlinear operating point;
it must not substantially retune the structural object. Size and profile are
macros over delay times, damping, dispersion, radiation, and regional coupling,
not a global linear pitch shift.

Mute, hi-hat pedal contact, and other constraints modify stored energy and
future loss continuously. They require explicit energy accounting: a slow
constraint change injects negligible energy, while a sufficiently fast new
contact may have its own distributed exciter.
