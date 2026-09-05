# TriggerFish nonlinear modal resonator architecture

This document describes the current constructive cymbal body. It is a mono,
stateful instrument made from contact excitation, one unified stochastic modal
field, passive loss, and observation filtering. It is not a finite-element
model and it does not claim that each control is a measurable physical property.
Its controls are intended to be audible, reasonably orthogonal, and fit-able to
recordings.

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
   `---- body force ---> unified stochastic modal field ----+--> mono output
                              |       |       |
                              |       |       `-- local passive exchange
                              |       `---------- phase diffusion
                              `------------------ upward energy cascade

mute/constraint ----------------------> additional modal loss only
```

One hit computes a new excitation projection and starts a contact gesture. Each
sample then executes, in this order:

```text
contact = ContactExciter::Process()
body = StochasticModalField::ProcessExcitedPair(contact.bodyDrive, 0, muteLoss)
output = ObservationModel::Process({contact.directRadiation, body})
```

Inside the modal field, recurrence happens before the upward cascade and local
exchange. The entire audio loop is fixed C++; it performs no graph traversal,
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

Global turbulence and its per-anchor multiplier control a normalized
centre-to-satellite trajectory. Increasing turbulence:

1. transfers excitation weight from the centre to its satellites;
2. widens their spread on the ERB-rate axis;
3. increases magnitude-preserving phase diffusion; and
4. increases passive exchange with compatible neighbours.

The squared excitation weights are normalized, so turbulence is not an
unlabelled gain control. A frequency-dependent turbulence field is

$$
t(f)=\operatorname{clamp}\!\left(
  t_c+s\log_2\!\frac{f}{f_c},\;0,\;1
\right)t_{\mathrm{local}},
$$

where $t_c$ is `Turbulence`, $s$ is `Turbulence slope`, and $f_c$ is
`Turbulence centre`. Positive slope keeps low packets more defined while making
high packets progressively wash-like. This is deliberately independent of the
painted modal prominence.

## Intrinsic bloom: stored-energy transport

Bloom is passive, one-way transport between frequency-ordered packets. It acts
directly on their complex states. A fixed half-octave transport stencil is
interpolated onto the available painted anchors. Adding an anchor between two
existing frequencies therefore refines the audible modal field; it does not
insert another serial bloom stage or change the meaning of the rate control.
For source energy $E_l$, one sample computes a bounded transfer fraction

$$
q = 1-\exp\!\left(-\frac{R\,a(E_l)}{F_s\,\Delta o}\right),
$$

where $R$ is `Bloom rate` in octaves per second, $\Delta o$ is the half-octave
stencil distance (shortened only at the upper field boundary), and $F_s$ is
sample rate. `Energy acceleration` controls

$$
a(E)=1+7d\frac{E}{E+E_{\mathrm{ref}}}.
$$

Here $E$ and $E_{\mathrm{ref}}$ are totals for the whole field, rather than
per-packet values. Consequently the activation does not change when an artist
inserts an intermediate anchor. At $d=0$, transport rate is level-independent.
At $d>0$, stored energy accelerates transport from the declared baseline $R$
toward $8R$. This deliberately broad range allows a forceful onset to bloom
quickly while a quiet late tail travels slowly. The cascade may therefore drain
lower packets upward; the shared
frequency-dependent T60 curve remains the final decay shaper. This is in
addition to, and testable
separately from, the visible velocity-brightness excitation tilt.
With the factory mapping, every increase in strike strength therefore both
injects more total energy and increases high-frequency excitation; it cannot
become darker merely because the input became stronger.

The event energy is deposited into the one or two anchors bracketing the target
log-frequency. Linear interpolation in octaves preserves the requested mean
travel distance. If a field is sparser than the stencil, part of the event
remains at the source and the rest enters the next available anchor. The state
update is therefore, schematically,

$$
E_l'=(1-q)E_l+w_lqE_l, \qquad
E_u'=E_u+w_uqE_l, \qquad w_l+w_u=1.
$$

It is implemented by magnitude scaling of the actual complex packet states, so
apart from floating-point error,

$$
E_l'+E_u'=E_l+E_u.
$$

No signal is added to the output. If the upper packet was silent it is seeded
with exactly the transferred energy, distributed between its coherent centre
and stochastic satellites using the same normalized weights as a direct hit;
otherwise its existing state is scaled.
`Bloom phase diffusion` then rotates destination phasors by deterministic
signed angles. Rotation changes correlation and texture but preserves
magnitude.

All transfers use packet energies captured at the beginning of the sample.
Newly received energy therefore cannot cascade again in that sample. The bloom
is continuous from the beginning of a strike and progresses upward at a
declared rate; it cannot suddenly appear after an arbitrary timer. Repeated
hits add force to the state already present, so a crash can build and swell.

This mechanism is a constructive approximation to nonlinear energy transfer,
not a solver for thin-shell equations. It intentionally models the perceptual
property the control names: low-frequency stored energy progressively populates
higher, denser modal regions.

## Local diffusion and exchange

Two separate energy-neutral processes control density:

- phase diffusion rotates individual complex modal states;
- local exchange applies signed Givens rotations to neighbouring state pairs.

For two real state components,

$$
\begin{bmatrix}x_a'\\x_b'\end{bmatrix}=
\begin{bmatrix}\cos\theta&-\sin\theta\\
                \sin\theta& \cos\theta\end{bmatrix}
\begin{bmatrix}x_a\\x_b\end{bmatrix}.
$$

The same rotation is applied to the imaginary components, so paired energy is
unchanged. Unlike the upward cascade, local exchange is symmetric and does not
create a preferred spectral direction.

These operations have distinct jobs:

| Control | Changes |
| --- | --- |
| Packet spread | prepared satellite frequencies |
| Phase bandwidth | time-varying coherence/linewidth |
| Local exchange | passive short-range state mixing |
| Bloom rate | directed upward state-energy travel |
| Energy acceleration | how strongly travel accelerates with stored strike energy |
| Bloom phase diffusion | correlation of newly populated upper states |

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

The final observation has independent direct and body gains and static
radiation filters. Observation filtering changes recorded colour without
changing the body's stored energy or T60. An optional relative propagation
delay exists for future presentation work and is off by default. Stereo remains
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
- stronger energy-dependent upward transport with velocity brightness disabled;
- exact passivity of an isolated cascade within tolerance;
- invariant high-band arrival after inserting intermediate painted anchors;
- preservation of a packet's centre-to-sideband energy ratio on cascade entry;
- no same-sample retransmission of newly arrived cascade energy;
- energy-preserving phase diffusion and Givens exchange;
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
| Observation | `src/tfdsp/percussion/observation_model.hpp` |

Related documents:

- [Metallic percussion DSP components](TfPercussion-metal-components.md)
- [Crash fitting methodology](TfCrash-fitting-methodology.md)
- [Percussion analysis toolkit](TfPercussion-analysis-toolkit.md)
- [Percussion ear-fitting workbench](TfPercussion-ear-fitting-workbench.md)
