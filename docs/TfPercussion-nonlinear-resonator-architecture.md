# TriggerFish nonlinear modal resonator architecture

This document describes the current constructive cymbal body. It is a mono,
stateful instrument made from contact excitation, one unified stochastic modal
field, passive loss, and observation filtering. It is not a finite-element
model and it does not claim that each control is a measurable physical property.
Its controls are intended to be audible, reasonably orthogonal, and fit-able to
recordings.

For the current oscillator allocation, stable/doublet layouts, phase-blur
controls and regional texture fitting, see [Beating packets](TfPercussion-beating-packets.md).
The [control-surface review](TfPercussion-visible-controls-and-slow-beating.md)
records visible-by-default grouping, a proposed noisiness-parameter reduction,
and the slower gong-beating audit.
The current [modal-wander and EQ refinement](TfPercussion-modal-wander-and-eq.md)
covers fixed-pivot noisiness, irregular narrowband ringing, graphical EQ and
the live full-output spectrum.
That refinement supersedes the earlier width-limited allocator; the mono energy
path and spectral-diffusion law described here are unchanged.

Signal units and all gain stages are defined in
[the gain-staging contract](TfPercussion-gain-staging.md).

The defining design choice is that attack ridges, dense wash, and nonlinear
bloom are different behaviours of one modal state. There is no delayed bloom
sample, hidden latch, separately audible noise tail, or second resonator bank in
the active cymbal recipe.

## Active signal flow

![Current cymbal signal flow](TfNonlinear-resonator-signal-flow.svg)

```text
hit event
   |
   v
contact exciter -----------------------> direct radiation --+
   |                                                        |
   `---- body force ---> unified stochastic modal field ----+--> mix --> final EQ --> mono output
                              |       |       |
                              |       |       `-- smooth pitch wander
                              |       `---------- phase blur
                              `------------------ spectral-energy diffusion

mute/constraint ----------------------> additional modal loss only
```

One hit computes a new excitation projection and starts a contact gesture. Each
sample then executes, in this order:

```text
contact = ContactExciter::Process()
body = StochasticModalField::ProcessExcitedPair(contact.bodyDrive, 0, muteLoss)
mix = contactLevel * contact.directRadiation + bodyLevel * body
output = modelLevel * (eqEnabled ? RadiationFilter::Process(mix) : mix)
```

Inside the modal field, recurrence happens before spectral-energy diffusion.
The entire audio loop is fixed C++; it performs no graph traversal,
allocation, or JSON processing.

## Contact excitation

`ContactExciter` combines four finite gestures:

| Primitive | Perceptual role |
| --- | --- |
| Half-sine force | coherent compression and release |
| Damped chirp | hard-tip ping or bell articulation |
| Enveloped tilted noise | broadband collision |
| Dense micro-contact burst | brushes and rough contact |

It exposes a direct-radiation port and a body-force port. A noisy force can
therefore excite resonances without necessarily appearing as an unrelated dry
noise layer. Brush, mallet, and stick are families on one performance control;
their character control changes bristle stiffness, mallet firmness, or tip
hardness respectively.

Strike strength is not compressed by a hidden velocity curve. It scales the
incident force linearly and also changes physically plausible contact
properties: duration, chirp frequency, micro-contact contribution, and noise
bandwidth. A visible `Velocity brightness` control changes the high-frequency
modal coupling of new force. It never recolours energy already ringing.

## One body, represented by modal packets

![Modal-packet preparation and audio processing](TfNonlinear-resonator-modal-packet.svg)

The editor contains between zero and 32 active centre handles. Each handle has:

- a centre frequency;
- a relative observation prominence; and
- a local multiplier for the global turbulence field.

At preparation time every active handle reserves one coherent centre state.
Deterministic sideband pairs are then allocated from one shared pool of 512
complex modal states. Local turbulence and ERB spread request useful coverage;
the global satellite-density control scales that request, and overlap between
neighbouring packet widths reduces redundant satellites. Painted centres are
never removed by this heuristic. Allocation is performed outside the audio
loop and the resulting states are packed into one flat SoA bank.

For mode $m$, before cross-mode processing, the recurrence is

$$
z_m[n+1] = r_m e^{j\omega_m[n]}z_m[n]
           + q_m x[n],
$$

where $r_m$ comes from the shared T60 curve, $q_m$ contains the excitation
tilt, packet weight, and current strike projection, and $\omega_m[n]$ may
contain a small stochastic phase perturbation. The complete drive vector is
renormalized after location and velocity colour are applied:

$$
\sum_m q_m^2 = 1.
$$

This normalizes spatial coupling, not the work done by an entire contact.
Delivered energy also depends on contact waveform and existing modal phase.
Observation is

$$
y[n]=G_\mathrm{body}\sum_m g_m\operatorname{Re}(z_m[n]),
$$

where the painted anchor levels define literal $g_m$ amplitudes
and $G_\mathrm{body}$ is the explicit body-observation level. Painted levels do
not affect $q_m$, stored energy, or cascade activation. A handle at the editor's
silence boundary is absent from both vectors. With no drive or constraint, the
pole radius is the only loss in this recurrence.

Global noisiness and its per-anchor multiplier distribute excitation energy
between the coherent centre (or pair) and surrounding modes. Packet spread
sets their frequency range, density/local allocation sets the number of
oscillators, and phase blur sets stochastic coherence loss. No random neighbour
exchange is enabled in the current workbench recipe.

The noisiness field is defined at a fixed 1 kHz pivot:

$$
t(f)=t_{1\mathrm{kHz}}\left(\frac{f}{1000\,\mathrm{Hz}}\right)^s
t_{\mathrm{local}}.
$$

The relaxed response maps this smoothly to a bounded satellite-energy fraction;
it does not add energy. See [the packet design](TfPercussion-beating-packets.md)
and [current movement controls](TfPercussion-modal-wander-and-eq.md) for
allocation, paired-ring coefficients, phase blur and independent Hz wander.

## Intrinsic bloom: spectral-energy diffusion

The current workbench uses **bidirectional diffusion down the spectral-energy
density gradient**, not a prescribed upward velocity. A low-heavy strike
usually drives energy upward initially; direction follows the state and can
reverse. The old one-way cascade remains a generic library primitive, not the
active recipe's transfer law.

With normalized coordinate $x=f/(1000\,\mathrm{Hz})$, packet cell width $w_i$,
stored energy $e_i$ and reference input energy $E_{\mathrm{ref}}$:

$$
\rho_i=\frac{e_i}{E_{\mathrm{ref}}w_i}.
$$

A semi-implicit tridiagonal solve redistributes these densities using
nonnegative, energy-dependent conductances. Closed boundary fluxes preserve
total stored energy and keep it nonnegative. Packet amplitudes are rescaled;
previously silent packets receive exactly their allocated energy through the
same normalized centre/satellite weights used for excitation. There is no
separate noise signal or timed high-frequency injection.

The complete conductivity, boundary and discretization equations are maintained
in [spectral-energy diffusion](TfPercussion-spectral-diffusion.md).
Diffusion strength is a coefficient, **not octaves per second**. Nonlinearity
sets the energy dependence: zero is linear diffusion; positive values reduce
conductivity as energy fades. It is not a guarantee of a fixed bloom time.

Phase blur rotates modal states without changing their energy. Smooth pitch
wander changes instantaneous frequency without changing the pole radius.
Neither replaces diffusion or T60 damping. Random neighbour exchange and
arrival-phase randomization are disconnected in this recipe.

## Decay and mute

All modes use one shared frequency-dependent local-T60 curve. The ordinary
curve has only two active boundary values at the modal design limits, 40 Hz
and 15 kHz. Up to six interior knots
are available for sparse last-stage correction, but fitting may not use them to
hide errors in excitation, modal placement, turbulence, or bloom.

T60 is interpolated in ERB rate and log seconds and held at its boundary value
outside the design range. Fixed frequency boundaries keep a saved patch's loss
law invariant when sample rate changes. The curve supplies each mode's pole
radius

$$
r_m = 10^{-3/(F_s T_{60}(f_m))}.
$$

Bloom and exchange do not define another hidden low/mid/high decay system.
Moving energy upward can change the observed envelope because the destination
modes use the T60 appropriate to their frequencies; that is the intended
interaction between spectral transport and loss.

Mute is a smoothed multiplicative loss applied to the modal recurrence. It can
only remove energy. A change of mute on a silent body remains silent, and a slow
constraint movement injects no energy.

## Constructive colour controls

`Initial excitation tilt` and `Excitation centre` form a smooth shelf that
shapes where a strike deposits energy across the modal field. For mode
frequency $f$, centre $f_c$, and high-side slope $s$ in dB/octave, its
unnormalized gain is

$$
q(f)=\left(1+\left(\frac{f}{f_c}\right)^2\right)^{s/(40\log_{10}2)}.
$$

The response is flat below $f_c$ and approaches slope $s$ above it. A centre
therefore remains meaningful after normalization, unlike the pivot of a pure
power law. Painted modal levels control what is audible; this shelf controls
what is initially driven. The complete shelf/location/velocity
projection is energy-normalized, while the painted observation vector is
normalized independently. Consequently modal painting changes spectral
prominence, `Body observation level` changes overall audible level, and neither
changes stored strike energy. The tilt range is deliberately wide enough to
move from low-dominated gong starts to bright cymbal starts.

`Body excitation` is the explicit gain between the contact body port and this
modal field. It changes stored energy and therefore the operating point of the
energy-dependent cascade. `Body observation level` is downstream and changes
only the audible readout. There is no independent hidden graph-edge gain.

The turbulence level/slope/centre controls describe coherence and density, not
spectral amplitude. Body tune scales modal centre frequencies, while contact
chirp pitch changes only the impact. Neither is a global audio-domain pitch
shifter.

The final observation mixes independent direct and body levels, then applies
one shared high-pass/colour/low-pass EQ. Observation filtering changes recorded
colour without changing the body's stored energy or T60. The reusable observation
delay primitive remains available for future presentation work; this voice does
not instantiate it or any per-path EQ. Stereo remains
an output-presentation extension; the synthesized object is mono.

## Explicit non-features

The active metallic recipe does not contain:

- a dispersion-loop node or delayed secondary body input;
- a separately audible turbulent-noise reservoir;
- independent resolved and dense modal banks;
- trigger-relative bloom latches, timers, or envelopes;
- hidden connection gains; or
- per-mode fitted decay multipliers.

The tested allpass/self-phase `DispersionLoop` remains a reusable library
component for other recipes and experiments. It is not part of the current
cymbal renderer and must not be presented as its bloom mechanism.

## Verification contract

Automated tests establish structural properties, not perceptual calibration.
They cover:

- deterministic rendering and finite state at supported sample rates;
- monotonic output energy across a velocity sweep;
- increased high-frequency response for stronger strikes;
- passive diffusion with nonnegative energy and closed boundaries;
- packet-coordinate and coincident-anchor consistency;
- preservation of normalized centre/satellite excitation weights;
- bounded phase blur and smooth frequency movement;
- passive mute and zero-strength no-op;
- additive restrikes when nonlinear transport is disabled; and
- bounded long rendering at maximum bloom settings.

Reference spectrograms and listening remain necessary to choose rates,
turbulence, modal placement, decay, and observation settings.

## References and provenance

- Zion Jaymes,
  [cymbal-synthesis walkthrough](https://www.youtube.com/watch?v=netcpYINyBQ).
  This motivated separating contact, developing metallic texture, resonators,
  and observation. The current state cascade is a revision, not a literal copy
  of the tutorial's feedback/allpass graph.
- Travis Skare and Jonathan Abel,
  [*Real-Time Modal Synthesis of Crash Cymbals with Nonlinear Approximations,
  using a GPU*](https://www.dafx.de/paper-archive/2019/DAFx2019_paper_48.pdf),
  DAFx-19. This supports dense complex modal recurrence, persistent state, and
  energy-dependent extensions to a static bank.
- Gabriel Cirio, Ante Qu, George Drettakis, Eitan Grinspun, and Changxi Zheng,
  [*Multi-Scale Simulation of Nonlinear Thin-Shell Sound with Wave
  Turbulence*](https://www.cs.columbia.edu/cg/waveturb/), SIGGRAPH 2018. Its
  frequency-domain energy cascade motivates directed spectral-energy transport.
  TriggerFish does not implement its shell or diffusion solvers.
- Michele Ducceschi and Cyril Touzé,
  [*Modal approach for nonlinear vibrations of damped impacted plates:
  Application to sound synthesis of gongs and
  cymbals*](https://doi.org/10.1016/j.jsv.2015.01.029), 2015. This establishes
  the relevance of modal state, impact excitation, frequency-dependent loss,
  and nonlinear coupling for struck plates.
- Quoc Bao Nguyen and Cyril Touzé,
  [*Nonlinear vibrations of thin plates with variable thickness: Application
  to sound synthesis of
  cymbals*](https://doi.org/10.1121/1.5091013), 2019. This supports treating
  profile and strike region as modal-coupling changes rather than one pitch
  control.
- Dan Stowell,
  [*Cymbal synthesis tutorial*](https://mcld.co.uk/cymbalsynthesis/), an
  independent real-time, spectrogram-guided constructive approach.

The 17-state packet, normalized centre/satellite trajectory, state-level upward
cascade, stochastic phase rotation, and local Givens exchange are TriggerFish
designs assembled from these requirements. They are evaluated by their declared
invariants and listening results, not attributed to any single source.

## Source map

| Responsibility | Source |
| --- | --- |
| Instrument/event composition | `src/tfdsp/percussion/crash_cymbal.cpp` |
| Fit expansion and strike projections | `src/tfdsp/percussion/crash_cymbal_parameters.cpp` |
| Contact gesture | `src/tfdsp/percussion/contact_exciter.hpp` |
| Unified recurrence and local exchange | `src/tfdsp/percussion/stochastic_modal_field.hpp` |
| Directed state-energy transport | `src/tfdsp/percussion/modal_energy_cascade.hpp` |
| Passive live damping | `src/tfdsp/percussion/passive_constraint.hpp` |
| Final output EQ | `src/tfdsp/percussion/radiation_filter.hpp` |

Related documents:

- [Metallic percussion DSP components](TfPercussion-metal-components.md)
- [Crash fitting methodology](TfCrash-fitting-methodology.md)
- [Percussion analysis toolkit](TfPercussion-analysis-toolkit.md)
- [Percussion ear-fitting workbench](TfPercussion-ear-fitting-workbench.md)
