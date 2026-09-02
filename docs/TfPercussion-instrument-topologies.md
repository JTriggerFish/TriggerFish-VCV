# Percussion instrument topologies

This document maps the reusable component library onto snare, kick, cymbal,
and hi-hat topologies. The structures are candidates to validate, not retained
implementations.

## Lessons from the cymbal, snare, and kick studies

The following tutorial studies are useful perceptual prototypes rather than
authoritative physical derivations:

- Zion Jaymes,
  [cymbal-synthesis walkthrough](https://www.youtube.com/watch?v=netcpYINyBQ);
- Zion Jaymes, [Synthesize SNARES that sound REAL using the power of
  FEEDBACK](https://www.youtube.com/watch?v=1Db9rGbth_o);
- [the snare follow-up](https://www.youtube.com/watch?v=LxV5J6UpR8Y);
- [Synthesize Kick drums that sound REAL using the power of
  FM](https://www.youtube.com/watch?v=ndG-6-vONNc).

The snare prototype divides the sound into direct stick contact, a filtered
feedback-delay body, an independent bright snare component, and room response.
Delay sets the main feedback pitch, loss-filter magnitude sets decay, and the
loss filter's phase also changes the loop tuning. The follow-up demonstrates a
comb implementation, amplitude-followed filter movement as a tension effect,
and several delay paths in the common feedback loop. It also finds that adding
approximately 2 kHz contact energy and clipping makes the hit sound firmer, and
that randomizing the noisy part of the contact is more natural than strongly
randomizing its tonal component.

The kick prototype separates a close-microphone sound from a far-field sound.
Its close path is a very short falling sine whose phase or frequency is
modulated by a noise burst: the low impulse and beater click are consequently
correlated rather than two unrelated layers. Its far path is a longer,
approximately 100--200 Hz decaying carrier with a longer noise-FM envelope,
different EQ, and a small propagation delay. The quiet far path materially
changes the apparent size even when the close path dominates.

These demonstrations support several reusable components, but their convenient
DAW routing must not be mistaken for a complete physical topology. A physical
snare contains two membranes, an enclosed acoustic field, a rigid shell, and
distributed one-sided collisions between the lower membrane and the wires.
Bilbao's full time-domain model treats all of those explicitly
([JASA paper](https://www.pure.ed.ac.uk/ws/portalfiles/portal/11222120/1.3651240.pdf)).
Likewise, a bass drum's batter and resonant membranes exchange energy through
the enclosed air and shell; measured double-headed drums exhibit coupled
membrane/cavity behaviour rather than two unrelated sine generators
([bass-drum synthesis overview](https://www.soundonsound.com/techniques/synthesizing-drums-bass-drum)).

### Additional reusable components

#### Correlated FM burst

A reusable `CorrelatedFmBurst` contains:

- a carrier waveform and amplitude envelope;
- a carrier pitch trajectory;
- a noise or irregular oscillator modulator;
- an independent FM-index envelope;
- optional deterministic per-hit perturbation;
- oversampled generation and a radiation filter.

It can create a hard kick contact, a far-field low correlated burst, a
synthetic contact supplement, or the complete compact body response of a
perceptual drum model. Its role is determined by routing and envelope length,
not by a restriction that FM may occur only in the transient. Velocity may
change FM index and duration as well as amplitude. Tonal carrier phase and
tuning remain mostly stable while the noisy modulator is allowed greater seeded
variation.

The cymbal tutorial motivates self-PM by observing that sinusoidal FM/PM
sideband amplitudes follow Bessel functions and that circular membrane
solutions also involve Bessel functions. This is a useful creative hint, but
not a physical equivalence. In sinusoidal FM, `J_n(index)` sets the amplitudes
of equally spaced spectral sidebands. In a circular membrane, Bessel functions
describe spatial mode shapes and their boundary-condition zeros set modal
frequency ratios. An FM oscillator does not automatically reproduce those
ratios, strike-position projections, coupled heads, or modal decays merely
because both calculations contain Bessel functions.

FM/PM remains a strong candidate because it efficiently produces correlated,
energy-dependent spectral complexity. It should therefore be tested as a
selectable `BodyResponse` implementation alongside feedback-delay and arbitrary
modal implementations, without claiming that its special-function expansion
is itself a membrane solution.

#### Membrane resonator

A `MembraneResonator` must support arbitrary membrane mode ratios, split mode
pairs, location-dependent mode-shape projections, frequency-dependent loss,
and state-dependent tension. A circular membrane is not well represented by
one harmonic comb. The initial efficient implementation can be an arbitrary
modal bank; a circular digital waveguide mesh remains a later alternative.

State-dependent tension should be expressed explicitly as energy-dependent
modal frequency or propagation delay. An amplitude follower moving the cutoff
of a feedback-loop filter is an effective tutorial shortcut because the
filter's phase retunes the loop, but it entangles damping and pitch. Our model
should expose those two effects separately and permit the shortcut only as a
candidate reduced implementation.

#### Passive resonator coupling and acoustic cavity

**Status: deferred structured-model extension.** None of the compact
video-derived ride, hi-hat, snare, or kick graphs requires a generic
bidirectional body coupler. The existing orthogonal resonator mixer already
covers passive exchange among lines within one resonator bank. Do not implement
`EnergyCoupler` or `AcousticCavity` until a coupled membrane/body candidate is
selected and exposes energy-normalized contact ports.

An `EnergyCoupler` exchanges state bidirectionally without creating energy. It
connects two membranes through an `AcousticCavity`, whose reduced model may
contain low acoustic modes, propagation delay, frequency-dependent loss, and a
port or vent. The same abstraction can connect:

- snare batter head, cavity, and snare-side head;
- kick batter head, cavity, and resonant head;
- tom heads and shell;
- cymbal bell and plate;
- sympathetic drums in a kit, at a much lower coupling level.

The cavity is not room reverberation. It is part of the instrument and affects
its poles. The external room remains a separate renderer.

#### Distributed one-sided contact

**Status: deferred higher-fidelity extension.** It is not required by the
compact video-derived graphs or by the first ride and hi-hat implementations.
First test pedal-controlled passive loss plus driven stochastic contact for the
hi-hat, and the compact bright residual for the snare. Implement distributed
contact only if matched references expose state-dependent plate or wire
interaction that those smaller graphs cannot reproduce.

A `DistributedContactCoupler` converts relative displacement and velocity
between resonating objects into many bounded, one-sided collisions. It should
support contact gap, compliance, damping, spatial distribution, and a family of
small resonators or a correlated residual. It is reusable for snare wires,
partly open hi-hat plates, rattles, chains, and objects resting on a membrane.

Snare-wire sound should therefore be driven primarily by motion of the lower
head, not by an independent noise envelope launched directly from HIT. The
relative prominence of wires in a ghost note can still increase because the
body, contact nonlinearity, radiation, and bottom microphone scale differently.
A static touching state must not generate energy, and a slow change of contact
gap must approach passive behaviour.

#### Observation and microphone model

The synthesizer needs explicit source-state outputs before any room effect:

- direct contact;
- each resonating body or resonator group;
- cavity or port radiation;
- distributed-contact radiation;
- a mono instrument sum and optional stereo radiation view.

An `ObservationModel` combines those states using delay, polarity, radiation
EQ, directivity, and distance. This produces kick-in, kick-out, snare-top,
snare-bottom, overhead-like, and far-field views without rerunning the physical
core. Microphone coloration may be a fitted low-order filter; external room
response is handled by the existing reverb or measured convolution during
calibration.

The tutorial's two independent kick synthesizers are a useful approximation,
and form a valid compact `FmBodyResponse` candidate. The more physical candidate
derives close and far views from shared membrane and cavity state. A hybrid can
add a far-field FM residual when it is driven by the same contact event, shares
deterministic timing and variation, and improves held-out measurements without
concealing an incorrect body model.

#### Bounded nonlinear output processing

**Status: deferred pending feature design.** A single saturator is too narrow
an abstraction for the intended optional production stage. The later design
should consider oversampled saturation, compression, trigger-aware envelope
shaping, filtering, and deliberately high-quality sample-rate or bit-depth
reduction as one composable character processor. It must be exactly bypassable
during body calibration, expose every gain change, and never let an optimizer
use processing to disguise incorrect contact levels or resonances.

### Candidate snare topology

```text
stick/brush contact -------------------------------> direct/top radiation
        |
        v
batter membrane <-> cavity/shell <-> snare-side membrane
                                          |
                                          v
                              distributed wire contacts
        |                 |                |
        +-----------------+----------------+-> observation model
                                                   |-> top
                                                   |-> bottom
                                                   `-> room send
```

The feedback-delay body demonstrated in the videos is a useful reduced
`MembraneResonator` candidate. Several fractionally tuned delay lines with
frequency-dependent feedback and controlled coupling are preferable to one
raw loop. They must be tested against an arbitrary modal membrane bank, using
the same contact, cavity, wire-contact, and radiation components.

### Candidate kick topology

```text
beater contact / correlated FM supplement
        |
        v
batter membrane <-> cavity/shell <-> resonant membrane/port
        |                 |                  |
        +-----------------+------------------+-> observation model
                                                    |-> kick in
                                                    |-> kick out
                                                    `-> far/room send
```

The falling pitch trajectory belongs partly to the contact impulse and partly
to state-dependent membrane tension. The two contributions must remain
separately observable. Resonant-head tuning is a real second subsystem, not
only a narrow EQ boost, although the latter remains a useful reduced-model
ablation.

### Calibration consequences for membrane drums

Multi-microphone sample libraries are especially valuable. Matched exports of
one MIDI performance from a multi-output drum library can provide kick-in,
kick-out, snare-top, snare-bottom, overhead, and room channels for the same
physical event. This permits a staged identification:

1. fit contact and batter-head state primarily from close/top channels;
2. fit head-to-head and cavity transfer from onset lag, phase, and low-mode
   splitting across close channels;
3. fit snare-wire activation from the bottom channel and snares-on/off data;
4. fit radiation and microphone transfer jointly across close and overhead
   views;
5. fit external room response only after the dry instrument passes.

Additional objectives include time-varying fundamental/pole trajectories,
cross-channel onset and coherence, delayed lower-head or resonant-head energy,
wire activation density and decay, port/cavity band energy, and close/far level
relationships. Reference manifests must record microphone channel, polarity,
bleed configuration, processing, head tuning, damping, beater or stick, strike
location, and velocity.

The first snare dataset should include matched top and bottom channels,
multiple velocities, centre and off-centre strikes, and preferably snares-on
and snares-off takes. The first kick dataset should include matched in, out,
and room channels across velocity, with one fixed beater and damping state.
Other drum sizes, head tensions, implements, and rooms remain held out until
one instrument works across those initial cells.

### Candidate hi-hat topologies

The tutorial explicitly removes the ride/crash dispersion-bloom loop for
hi-hats. Its compact topology is:

```text
short pitched oscillator + shaped trailing noise
                    |                 |
                    +-> direct output |
                    `-> parallel comb resonators -> filtered resonator output
```

The oscillator supplies the short fundamental/contact articulation and uses a
pitch envelope. The noise supplies the sizzly tail. Both drive the same
parallel comb bank which supplies the characteristic modes. For a foot chick,
the source receives a small attack and the noise decays faster than the tonal
component. The tutorial treats shaped noise plus oscillator as a complete
replacement for the cymbal bloom/dispersion stage, not as an additional layer
on top of it.

This compact topology should be the first hi-hat ablation because it is simple
and identifiable. The advanced topology remains available:

```text
contact -> top and bottom resonator networks
                         <-> distributed plate contact
                                      |
                                      `-> shaped contact drive into resonators
```

The advanced graph explains pedal-continuous damping, reopening, repeated
inter-plate collisions, and energy already stored in the two plates. It still
does not require the ride/crash bloom loop: the distributed-contact component
can reduce to the compact shaped-noise/oscillator response when detailed
two-body motion is unnecessary. Compact, hybrid, and coupled versions must use
the same event/control interface and be compared against the same pedal,
velocity, and location dataset.

### Model fidelity is a topology choice

The component library deliberately permits several models behind the same
semantic role:

| Role | Compact perceptual option | Structured option | Coupled option |
| --- | --- | --- | --- |
| Body response | FM/PM oscillator and noise | feedback-delay or modal bank | several coupled bodies |
| Cymbal propagation | shaped stochastic response | nonlinear dispersion loop | regionally coupled loops |
| Hi-hat interaction | closure-shaped noise/oscillator | one resonator plus contact state | two bodies with distributed collisions |
| Drum interior | fitted EQ/resonance | reduced cavity modes | reciprocal heads/cavity/port |

Greater structural detail is not accepted automatically. Start with the
smallest topology that exposes the required controls, then add a component only
when held-out measurements or listening identify a failure that the new state
can explain.

## Tutorial technique coverage audit

The generic graph can reproduce every synthesis technique demonstrated in the
three Zion Jaymes tutorials, while keeping the tutorial recipes separate from
claims about physical mechanism. A technique being representable does not make
it mandatory in every calibrated instrument.

### Snare and follow-up

| Tutorial technique | Generic representation | Architectural status |
| --- | --- | --- |
| 4--15 ms oscillator-plus-noise stick transient | `ContactExciter` pulse, chirp, and noisy contact | Core |
| Transient heard directly and sent into the head | separate direct-radiation and body-drive ports | Core |
| Sample delay sets head pitch | `FeedbackDelayBody` fractional main delay | Compact body option |
| Feedback gain sets ring time | decay-time parameter converted to loop loss | Core resonator feature |
| EQ/filter inside feedback supplies damping | per-circulation frequency-dependent loss | Core resonator feature |
| Filter phase also retunes the loop | measured complete loop phase; optional coupled damping/tune shortcut | Supported ablation |
| Limiter protects arbitrary DAW feedback | bounded loop gain, finite-state guard, optional safety limiter | Safety requirement |
| Comb-filter plugin as Karplus--Strong body | one-line or multiline `FeedbackDelayBody` | Compact body option |
| Several delays in the common feedback bus | parallel/coupled delay network or explicit serial delay graph | Supported topology |
| Envelope follower moves a high-pass filter to imitate tension | energy-followed loop-filter mode | Compact nonlinear option |
| Explicit membrane tension trajectory | energy-dependent resonator frequency/delay | Structured alternative |
| Separate bright snare-noise track | `ShapedStochasticResponse` routed directly | Tutorial-compatible compact option |
| Lower membrane drives snare wires | `DistributedContactCoupler` plus wire resonators/residual | Structured alternative |
| Approximately 2 kHz contact emphasis | contact-radiation filter | Voicing parameter |
| Heavy clipping for harder impact | bypassable oversampled nonlinear output stage | Production/observation option |
| Randomize transient noise more than tonal phase | seeded micro-contact variation with stable structural state | Core variation rule |
| Ghost notes reduce body/contact more than apparent wires | velocity-dependent contact/coupling/radiation maps | Calibration target |
| Lower drums use longer snare decay | fitted size/tension-to-contact-decay macro | Calibration hypothesis |
| Room supplies size and power | instrument taps sent to external reverb/room observation | Core separation |
| Pre-fader sends and low DAW buffer size | explicit graph routing and sample-accurate internal feedback | Implementation detail |

The independent snare-noise track remains available because it is part of the
tutorial and may win a compact-model comparison. It is not silently identified
with physical wire collision; the structured candidate must earn its extra
complexity against matched top/bottom microphone data.

### Kick

| Tutorial technique | Generic representation | Architectural status |
| --- | --- | --- |
| Very short falling sine for the close thump | `CorrelatedFmBurst` carrier pitch and amplitude envelopes | Compact body/contact option |
| Noise FM creates a correlated beater click | band-limited irregular modulator with short FM-index envelope | Compact body/contact option |
| Roughly 1 ms hold and 2--10 ms FM decay | independent FM-index envelope | Fit parameter |
| Noise/irregular waveform changes attack character | selectable seeded modulator spectrum | Implement parameter |
| Velocity changes FM amount and duration | velocity-to-index and velocity-to-time maps | Core mapping |
| Close-microphone bass/mid/high EQ | `ObservationModel` close-path radiation/microphone filter | Observation parameter |
| Narrow 60--80 Hz resonant-head boost | low resonant observation filter | Reduced-model ablation |
| Explicit resonant head | second `MembraneResonator` coupled through `AcousticCavity` | Structured alternative |
| Longer falling carrier for electronic kicks | longer FM-body amplitude/pitch envelopes | Creative range |
| Far sound uses a 100--200 Hz decaying carrier | second or shared-state `FmBodyResponse` | Compact far-body option |
| Longer noise-FM burst makes boom and rattling | far-path FM-index envelope | Compact far-body option |
| FM index has a boom-to-burst sweet spot | bounded fitted FM-index control | Calibration target |
| Far path has less sub and different EQ | far-field radiation/observation filter | Observation parameter |
| Small far-path delay | physical observation propagation delay | Core observation feature |
| Haas widening | optional stereo observation delays | Production option |
| Quiet far path materially increases apparent size | shared event feeding calibrated close/far observation mix | Calibration target |
| Pure noise shaped with bass boost, distortion, high-pass, and dynamics | arbitrary source, filters, waveshaper, and envelope/dynamics nodes | Creative compact topology |
| Supersaw or arbitrary recordings used instead of noise | generic external/oscillator source node | Creative topology |

The compact two-FM-layer kick is a first-class candidate, not merely an effect
added to a mandatory physical kick. The coupled-head/cavity graph is a second
candidate using the same controls and observation taps. A hybrid is accepted
only when ablation shows that each branch explains a distinct measurement.

### Cymbal and hi-hat

| Tutorial technique | Generic representation | Architectural status |
| --- | --- | --- |
| Very short sine with extreme FM or pitch envelope | `ContactExciter` tonal chirp | Core contact option |
| Contact is heard directly | direct contact-radiation port | Core routing |
| Contact also drives dispersion | body-drive port into `DispersionLoop` | Ride/crash topology |
| Base/sample delay inside feedback | dispersion-loop base delay | Core dispersion stage |
| Serial allpasses smear each circulation | allpass stack inside dispersion feedback | Core dispersion stage |
| Slow vibrato creates irregular frequency spreading | slow modulated fractional delay | Candidate dispersion stage |
| Phase distortion maps signed audio to delay | self-phase fractional-delay stage | Candidate nonlinear stage |
| Tone filters the phase-modulation signal | modulation-path filter | Core self-PM control |
| Normalize reduces level dependence | continuously variable modulation normalization | Optional self-PM control |
| Drive controls velocity-sensitive bloom | incident energy plus self-PM drive | Calibration target |
| Feedback loop contains every dispersion processor | one explicit outer-loop graph | Core routing |
| Dispersion is not heard directly | analysis-only dispersion tap | Core routing |
| Dispersion drives the body renderer | secondary excitation port of the unified stochastic modal field | Core routing |
| Plugin dry signal is polarity-cancelled | native wet-only implementation; no dry-null hack required | Equivalent implementation |
| Approximately twelve parallel combs | optional parallel/coupled dense-residual lines | Residual candidate |
| A few dark lower-order groups plus messy high groups | clustered resonator parameters and group filters | Calibration structure |
| Per-resonator filters | per-line/group radiation and loss filters | Core resonator feature |
| Frequency shifter breaks residual comb relationships | optional fixed-hertz residual-group shift; sparse modes are placed directly | Candidate coloration |
| Additional comb over the summed output | optional `SpectralColorationComb` node | Candidate coloration |
| Different delays/allpasses stop branch flanging | optional per-branch delay/allpass | Candidate correction |
| Resonator buses are freely mixed | explicit group output gains | Core observation routing |
| Serial comb mode gives a special ringy hat | serial `ResonatorNetwork` graph mode | Creative/special articulation |
| Whole-output EQ | radiation/observation filter | Core voicing |
| Hi-hat replaces bloom with oscillator plus shaped noise | `ShapedStochasticResponse`; no `DispersionLoop` | Compact hi-hat topology |
| Hi-hat source is both direct and sent to combs | direct/source and resonator-drive ports | Core compact-hat routing |
| Foot chick adds attack and changes relative source decays | closure-velocity contact event and source envelopes | Compact pedal articulation |
| Ride uses bright direct stick and low dispersion drive | ride graph macro over contact and dispersion | Instrument calibration |
| Ride size/hammering change delay, vibrato, allpasses, and filtering | graph-level size/profile macro, not global pitch | Instrument calibration |
| Small direct stick feed may also reach resonators | optional contact-to-resonator bypass | Candidate routing |
| Bell may use a second dispersion and resonator chain | optional region-specific subgraph with passive/shared coupling candidate | Structured ride option |
| Crash uses stronger nonlinear dispersion and longer feedback | crash graph macro over shared components | Instrument calibration |
| Limiting protects the nonlinear loop | bounded loop, finite guard, and optional oversampled safety limiter | Safety requirement |

The tutorial's special-effect routings are supported without forcing them into
the neutral acoustic graph. In particular, serial combs, Haas delay, clipping,
and arbitrary-source kick processing are selectable graph nodes and remain
bypassed during first-pass body identification.
