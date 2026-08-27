# TriggerFish Elements

Circuit-modelled sound generators and processors, plus pitch utilities for VCV Rack 2.

[![CI](https://github.com/JTriggerFish/TriggerFish-VCV/actions/workflows/ci.yml/badge.svg)](https://github.com/JTriggerFish/TriggerFish-VCV/actions/workflows/ci.yml)

[Release notes for 2.4.0](docs/releases/2.4.0.md)

<table>
  <tr>
    <td align="center"><a href="#slop-and-slop-4"><img src="doc/TfSlop.png" height="260" alt="Slop module"><br><strong>Slop</strong></a></td>
    <td align="center"><a href="#slop-and-slop-4"><img src="doc/TfSlop4.png" height="260" alt="Slop 4 module"><br><strong>Slop 4</strong></a></td>
    <td align="center"><a href="#vdpo"><img src="doc/TfVDPO.png" height="260" alt="VDPO module"><br><strong>VDPO</strong></a></td>
  </tr>
  <tr>
    <td align="center"><a href="#vca"><img src="doc/TfVCA.png" height="260" alt="VCA module"><br><strong>VCA</strong></a></td>
    <td align="center"><a href="#303-oscillator"><img src="doc/Tf303Oscillator.png" height="260" alt="303 Oscillator module"><br><strong>303 Oscillator</strong></a></td>
    <td align="center"><a href="#303-voice-core"><img src="doc/Tf303VoiceCore.png" height="260" alt="303 Voice Core module"><br><strong>303 Voice Core</strong></a></td>
  </tr>
  <tr>
    <td align="center"><a href="#4072-voice-core"><img src="doc/Tf4072VoiceCore.png" height="260" alt="4072 Voice Core module"><br><strong>4072 Voice Core</strong></a></td>
    <td align="center"><a href="#electric-piano"><img src="doc/TfElectricPiano.png" height="260" alt="Electric Piano module"><br><strong>Electric Piano</strong></a></td>
    <td align="center"><a href="#wavefold-oscillator"><img src="doc/TfWavefoldOscillator.png" height="260" alt="Wavefold Oscillator module"><br><strong>Wavefold Oscillator</strong></a></td>
  </tr>
  <tr>
    <td align="center"><a href="#unison-oscillator"><img src="doc/TfUnisonOscillator.png" height="260" alt="Unison Oscillator module"><br><strong>Unison Oscillator</strong></a></td>
    <td align="center"><a href="#scene-pack-4"><img src="doc/TfScenePack4.png" height="260" alt="Scene Pack 4 module"><br><strong>Scene Pack 4</strong></a></td>
    <td align="center"><a href="#room-reverb"><img src="doc/TfReverb.png" height="260" alt="Room Reverb module"><br><strong>Room Reverb</strong></a></td>
  </tr>
  <tr>
    <td align="center"><a href="#transport"><img src="doc/TfTransport.png" height="260" alt="Transport module"><br><strong>Transport</strong></a></td>
    <td align="center"><a href="#prog-sequencer"><img src="doc/TfProgSequencer.png" height="260" alt="Prog Sequencer module"><br><strong>Prog Sequencer</strong></a></td>
    <td></td>
  </tr>
</table>

## Polyphony

| Module | Handling | Output voice count |
| --- | --- | --- |
| Slop | Polyphonic, with shared hum and independent drift per channel | Channel count of `1V/OCT`, up to 16 |
| Slop 4 | Four polyphonic oscillator paths, with instrument-wide hum/common drift and independent drift per channel and oscillator | Each output follows its corresponding input, up to 16 |
| VDPO | Fully polyphonic, with independent oscillator and resampling state per voice; mono inputs are broadcast | Widest connected input, up to 16 |
| VCA | Fully polyphonic, with independent nonlinear, resampling, and DC-rejection state per voice; mono inputs are broadcast | Widest connected input, up to 16 |
| 303 Oscillator | Fully polyphonic, with independent DSP state per voice; mono inputs are broadcast | Widest connected input, up to 16 |
| 303 Voice Core | Fully polyphonic, with independent filter, envelope, accent, and VCA state; mono controls are broadcast | Channel count of `IN`, up to 16 |
| 4072 Voice Core | Fully polyphonic, with independent filter, dual-envelope, and VCA state; mono controls are broadcast | Widest connected input, up to 16 |
| Electric Piano | Sixteen independent modal tine/pickup note channels, sixteen internal release-tail slots, and one shared nonlinear amplifier | Widest pitch, gate, or velocity input, up to 16; `DIRECT POLY` preserves the originating channels |
| Wavefold Oscillator | Fully polyphonic, with independent oscillator, folder, and resampling state; mono controls are broadcast | Widest connected input, up to 16 |
| Unison Oscillator | Fully polyphonic, with an independent oscillator stack and drift state per channel; mono controls are broadcast | Widest connected input, up to 16 |
| Scene Pack 4 | Concatenates four mono/poly source bundles in jack order | One polyphonic audio output, up to 8 channels |
| Room Reverb | Mono by default; accepts and places up to eight source channels | Independent stereo left/right outputs |
| Transport | Monophonic shared clock, run state, and reset | Clock pulse, Run gate, and Reset pulse |
| Prog Sequencer | Polyphonic chord-capable program with independently cycling control lanes and optional Notes-pass edge alignment | Up to 16 pitch/gate/trigger/velocity/accent channels plus CV1-CV3 |

## Modules

### Electric Piano

Electric Piano is a sample-free, velocity-sensitive instrument voice loosely
inspired by a Rhodes Mk1. A finite-mass hammer is coupled bidirectionally to a
modal tine/tone-bar system through a Hunt-Crossley contact law. Compression and
relative hammer/tine velocity determine force and hysteretic loss at every
oversample, so contact duration, rebound, and modal energy distribution emerge
from the collision rather than from a prescribed brightness envelope. Two
transverse tine coordinates then traverse a finite-pole magnetic-flux field.
The sixteen pickup voices meet in one shared, stateful amplifier, so chords
compress and overload differently from isolated notes.

The resonator state represents displacement. An initial cantilever-scaling
profile derives modal mass from each key's exact pitch; modal mass and frequency
together determine displacement per impulse. A separate per-key pickup-
sensitivity curve balances small-signal level without removing the larger
displacement of low tines. Repeated strikes add impulses to existing motion.
When a Rack MIDI channel is reassigned, its released note moves into one of
sixteen internal tail slots instead of being silenced or retuned.

- **Vel Curve** curves incoming 0–10 V velocity to suit the controller, while
  **Dynamics** sets the range of physical hammer speed. That speed changes both
  loudness and timbre through the coupled collision and pickup excursion. Its
  default is the full incoming range; reduce it only for a less dynamic
  controller response. An inaudible sub-MIDI hammer-speed floor treats tiny
  positive CV residue as zero before contact begins; active collisions are never
  ended by an elapsed-time cutoff.
- **Body** and **Bell** balance the coupled fundamental body and the signed,
  short-lived inharmonic modes. **Hammer** changes neoprene stiffness,
  hysteretic loss, and contact duration rather than directly boosting Bell or
  treble; as a strike property, changes apply to newly played notes. **Coupling**
  changes the common-base stiffness of a genuine two-coordinate tine/tone-bar
  system. Its normal-mode frequencies, hammer projections, modal masses, decay,
  and tine pickup motion are derived together. Weak coupling lets the tine drive
  energy into the lossy mounting block; strong anti-phase tine/bar motion cancels
  that reaction, raising Q and adding a short tone-bar sub-fundamental bloom.
  This changes body and sustain rather than crossfading levels.
- **Tone** changes the tine's vertical alignment to the pickup pole and hence
  the balance of fundamental and curvature-generated harmonics. **Gap** changes
  front-to-back pickup proximity and therefore both level and dynamic response.
  `TONE CV` accepts polyphonic modulation at 10% of the knob range per volt.
- **Decay** scales natural modal decay, **Release** controls damper speed, and
  **Mechanics** ranges from a clean pickup sound to deliberately emphasized,
  differentiated keybed, key-release, and damper events. The calibrated default
  remains subordinate to the pitched voice. `PEDAL` holds released notes while
  its gate is high.
- **Drive** controls level into the post-sum electronics. The amplifier combines
  harp/load voicing, asymmetric feedback saturation, level-dependent supply
  recovery, and restrained cabinet filtering rather than a static waveshaper.
  `DIRECT POLY` exposes each pickup channel; `LEFT` and `RIGHT` expose the shared
  amplified output (currently dual mono).

The default curves vary spectrum and decay by velocity and keyboard position.
There is no dedicated bark detector or bass enhancement: low-note bark must
arise from the resonator motion, nonlinear magnetic pickup, and shared
amplifier. Its presence and onset are therefore calibration evidence for the
physical model rather than a directly imposed feature.

The magnetic pickup and shared amplifier nonlinearities run at 4x sample rate.
High resonator modes taper smoothly before Nyquist. Live timbre controls use
short smoothing; Decay uses a slower transition and staggered control-rate
coefficient updates so polyphonic adjustment remains continuous and inexpensive.

### Slop and Slop 4

Slop adds slow pitch drift and 60 Hz power-supply hum to a 1 V/octave signal.
The drift can be proportional, measured in cents, or linear, measured in hertz.
On a polyphonic cable, hum remains common to the instrument while every channel
has its own drift process. The tracking control also allows small
oscillator-tracking errors to be modelled.

Slop 4 treats each Rack polyphonic channel as one synthesizer voice and its four
paths as four oscillators within that voice. Proportional thermal drift and hum
are shared by the complete instrument; linear drift is independent for every
voice/oscillator pair. This produces coherent ensemble motion together with
evolving voice-to-voice and oscillator beating. Each path retains its own
tracking trim, and each output follows the channel count of its corresponding
input.

### VDPO

VDPO is a driven [Van der Pol oscillator](https://en.wikipedia.org/wiki/Van_der_Pol_oscillator).
It can self-oscillate, or lock imperfectly to a signal patched into `IN`, moving
from harmonic locking through inharmonic and chaotic responses.

- **Self freq** sets the oscillator's natural frequency.
- **Damping** sets the nonlinear damping and therefore affects both tuning and
  waveform complexity.
- **In** controls the forcing-signal level.
- **Damp** accepts bipolar damping modulation.

The stiff differential equation is evaluated by a structure-aware, prewarped
split integrator at 4x sample rate. Adaptive work is confined to the demanding
high-frequency, high-damping region. Audio, pitch, and damping controls are
reconstructed at 4x; pitch is converted to frequency inside that path.

### VCA

VCA is a transistor-inspired amplifier with separate nonlinear audio and control
paths, followed by a saturating output stage. The nonlinear amplifier core runs
at 2x sample rate. Linear and exponential CV are reconstructed separately before
the exponential control law and limiting are applied. CV bleed follows the
reconstructed control; final DC rejection runs at Rack's sample rate.

- **Drive** increases input saturation while compensating most of the associated
  level change.
- **Lin** and **Exp** set the depth of the two 0–10 V control inputs. Both default
  to 100%.
- **Curve** changes the exponential CV response, while **Bleed** adds a small,
  high-passed control transient to the output.
- **Output** sets the final level and the amount driven into the output stage.

For conventional AR and ADSR envelopes, the linear input is usually the most
predictable starting point.

### 303 Oscillator

303 Oscillator combines a band-limited saw core with a circuit-modelled TB-303
Q8 saw-to-square shaper. The module supports polyphony,
voltage-controlled waveform morphing, hard sync, extended tuning, and two FM
responses.

- **Octave** selects an integer offset from −3 to +3 octaves; **Tune** covers
  ±7 semitones.
- `1V/OCT` is the pitch input. `SLIDE` enables the internal pitch glide while
  high, and **SL.Time** sets its nominal duration from 2 to 360 ms.
- **Shape** shifts the square-wave switching point. **Wave** morphs continuously
  from saw to square with consistent level through the midpoint.
- `FM` is exponential in **EXP** mode and signed linear through-zero FM in **TZ**
  mode. The attenuverter sets depth in either mode.
- `SYNC` performs fractionally timed, minBLEP-corrected hard sync.
- `CV OUT` exposes the post-slide pitch signal for driving another oscillator or
  filter; `OUT` is the audio output.

The widest connected input sets the output channel count. Each voice has its own
oscillator, shaper, slide, sync, and resampling state; mono modulation inputs are
shared across the active voices.

The oscillator uses 4x oversampling by default. A 2x mode is available from the
context menu for dense polyphonic patches where CPU use matters more than the
extra alias suppression at high pitch or under complex FM and sync.
The [303 Oscillator technical report](docs/Tf303Oscillator-technical-report.md)
describes the model and its antialiasing in detail.

### Wavefold Oscillator

Wavefold Oscillator combines a band-limited sine/triangle oscillator with three
selectable folding characters. **Hinge** is a smooth four-stage design;
**Lockhart** and **Serge** follow the circuit models analysed by
[Esqueda et al.](https://doi.org/10.3390/app7121328). Lockhart is the brighter,
more harmonically active response, while Serge has a smoother, darker contour.

- **Octave** selects −3 to +3 octaves and **Tune** covers ±7 semitones.
- **Sine–Tri** morphs the source waveform. **Fold** controls input drive through
  the selected cascade, and **Symmetry** offsets its operating point.
- `WAVE`, `FOLD`, and `SYM` each have a bipolar CV attenuverter. FM can use
  exponential or signed through-zero response.
- **Alive Speed** sets a common 0.5–120 s mean-reversion time for three
  independent drift processes. The **Wave**, **Fold**, and **Sym** trimmers set
  their individual excursions; all three depths default to 50%.
- **Voices** selects one to four true-unison oscillator voices. **Spread** sets
  outer-voice offsets up to ±50 cents while preserving the centre pitch.
- `OSC OUT` provides the source before folding. Patching a nominal ±5 V signal
  into `EXT IN` replaces the internal oscillator at the folder input;
  `FOLD OUT` provides the processed signal.

The nonlinear path runs at 4x sample rate by default. The internal oscillator
gradually reduces extreme fold depth at high notes to keep its harmonic density
musically useful. External signals retain the requested fold depth. A 2x mode
is available from the context menu. Each polyphonic voice has independent
oscillator, folder, resampling, and Alive drift state.

### Unison Oscillator

Unison Oscillator builds a normalized stack of one to sixteen band-limited saw
or pulse voices for each polyphonic input channel. Voice-count changes fade and
reposition the stack smoothly, avoiding abrupt level and stereo jumps.

- **Octave** selects −3 to +3 octaves, **Tune** covers ±7 semitones, and
  **Voices** selects the stack size from one to sixteen.
- **Wave** selects saw or pulse. **P Width** sets pulse width from 5% to 95%; the
  bipolar **PWM** amount uses `PW CV` when patched, or the internal
  0.03–10 Hz **Rate** oscillator when unpatched.
- **Spread** distributes voices with outer offsets up to ±50 cents. **Width**
  distributes the stack across the stereo outputs. Both have bipolar CV
  attenuverters.
- **Sub Level** mixes a pulse sub-oscillator one octave down into all outputs.
  **Centre** uses one oscillator at the central pitch; **Stack** gives every
  unison voice its own detuned and panned sub.
- The built-in Slop section adds common 60 Hz hum, common proportional drift,
  individual drift in cents or hertz, and per-voice tracking error. **Speed**
  sets the drift time from 0.5 to 120 seconds.
- `MONO` is the normalized stack mix. `LEFT` and `RIGHT` use equal-power
  panning and match the mono level when **Width** is zero.

The widest connected input sets the output channel count, up to 16. Each Rack
channel owns an independent one-to-sixteen-voice oscillator stack and
individual drift processes; common drift and hum are shared across channels.

### 303 Voice Core

303 Voice Core is a circuit-modelled TB-303 voice back end. It includes the
four-stage diode ladder, coupling and resonance networks, filter and volume
envelopes, accent path, and BA662-style OTA VCA. `LP OUT` exposes the filter directly;
`VCA OUT` provides the complete articulated signal.

- **Cutoff** spans 10 Hz to 20 kHz. `1V/OCT` tracks directly, while `EXP CV`
  provides a second exponential input with an attenuverter.
- **FM** applies bipolar, AC-coupled linear cutoff modulation. **Res** modulates
  resonance, and **Res Range** selects the stock or extended feedback range.
- **Drive** sets the level entering the nonlinear ladder. **Bass** continuously
  restores low-frequency response lost in the original coupling network.
- `GATE` drives the internal filter and volume envelopes. `ACC` adds the accent
  response. The five envelope controls set filter depth, normal and accented
  decay, accent amount, and VCA decay/sustain.
- **Accent Sweep** provides Off, Fast, Normal, and Slow responses. Normal is the
  stock setting; Fast and Slow extend the timing in the spirit of Devil Fish.
- Patching `VCA CV` replaces the internal volume envelope. Its attenuator controls
  the external 0–10 V signal, while accent remains additive.

> **Level warning:** Bass is an extended boost control. With high resonance,
> large Bass settings can exceed Rack's nominal audio range and clip downstream
> modules. Reduce Bass, Drive, or mixer gain if this occurs.

The `IN` cable sets the filter's polyphonic voice count. Each audio channel has
independent filter, envelope, accent, and VCA state; mono control inputs are
shared across those voices.

The filter uses 4x oversampling by default, with a 2x context-menu option. The
context menu also selects stock or Devil Fish volume-envelope timing. Circuit
analysis, equations, calibration, and numerical validation are collected in the
[303 Voice Core technical report](docs/Tf303VoiceCore-technical-report.md).

### 4072 Voice Core

4072 Voice Core combines a circuit-modelled ARP 4072 four-pole low-pass filter,
two switchable AR/AD/ADSR envelope generators, and an ARP 4019-style discrete
VCA.
The filter's outer differential pair preserves the original unequal audio and
resonance-return levels, while the four locally fed-back filter stages operate
around their circuit bias points.

- **Cutoff** and `F 1V/OCT` set the filter frequency. The internal filter
  envelope uses the original exponential cutoff law through `ENV>CUT`.
- **Input** covers attenuation through +24 dB of filter drive. The cutoff
  control extends to 20 kHz using ARP's published correction for the original
  4072 high-frequency limitation.
- The two A/D/S/R slider rows generate 0–10 V envelopes. Each row can be switched
  among AR, AD, and ADSR. The playable default uses AD for both envelopes, with
  minimum Attack and one-second Decay: a bright, gate-length-independent
  starting point derived from the [ARP patch book's](https://www.korg.com/us/support/download/manual/0/877/4471/)
  plucked-instrument settings. Sustain is inactive in AD, but defaults to 50%
  for the filter and 100% for the amplifier so selecting ADSR gives an immediate
  sustained sound. Release defaults to one second for AR and ADSR operation.
  **Curve** varies both envelope shapes around the original circuit response.
- The amplifier envelope's `LIN / EXP` switch selects the 4019 VCA response for
  all three envelope modes. It defaults to EXP for the sharper 10 dB/V contour
  used by classic plucked patches.
- `CUT MOD` and `VCA MOD` each have one input, one amount control, and a
  `LIN / EXP` switch. Filter LIN is signed 200 Hz/V frequency modulation;
  filter EXP is 1 V/oct at 100%. VCA LIN follows the 4019 linear-gain input,
  while VCA EXP follows its 10 dB/V control law.
- The adjacent `+ / EXT` switches add each external modulation signal to its
  internal envelope or replace that envelope. An unpatched modulation input
  keeps the internal envelope active in either position. Both default to `+`.
- `VCA IN` replaces the normalled filter-to-VCA audio connection. `LP OUT`,
  `VCA OUT`, and both envelope outputs expose each stage independently.

#### Envelope modes

- **AR** attacks when Gate rises, holds its peak while Gate remains high, and
  releases when Gate falls. It is useful for straightforward note-length
  articulation, sustained sounds, and gate-following filter movement. AR uses
  the Board-4 response and ignores Trigger.
- **AD** is a retriggerable Attack–Decay one-shot independent of Gate length.
  It is useful for plucks, percussion, sequencer steps, and repeatable filter
  sweeps.
- **ADSR** attacks, decays to Sustain, holds while Gate remains high, and then
  releases. It is useful for keyboard articulation and sweeps with distinct
  peak and held levels. At maximum Sustain it can reproduce an AR level
  sequence; AR retains its different attack curve and Gate-only triggering.

The amplifier envelope law is independent of these modes. LIN gives direct
0–10 V gain control; EXP applies the original 10 dB/V VCA response for sharper
attack, decay, and sustain contours.

The original 2600 normally routes ADSR to filter cutoff. Its VCA receives AR
through the linear control input and ADSR through the exponential input. Its
plucked patches obtain an AD-like contour from ADSR with Sustain at zero. The
module's dedicated AD mode gives this contour timing independent of gate
length; select AR and LIN for the original gate-following AR path.

The complete filter/VCA path runs at 4x sample rate by default. A 2x mode is
available from the context menu for larger polyphonic patches.
Circuit analysis, equations, calibration, and numerical validation are in the
[4072 Voice Core technical report](docs/Tf4072VoiceCore-technical-report.md).

### Scene Pack 4

Scene Pack 4 combines four mono or polyphonic source bundles into one
polyphonic audio cable. Connected channels are concatenated in jack order,
preserving channel order within each input, and the result is capped at Room
Reverb's eight-source maximum. The polyphonic output of another Scene Pack can
feed any input to chain bundles without separate position cables.

Source placement belongs to Room Reverb: every channel present on its AUDIO
input appears as a numbered draggable marker in the room plan. Disconnected or
zero-channel Scene Pack inputs consume no output channel and therefore create
no marker.

The topology, signal model, equations, design rationale, latency semantics, and
implementation invariants are described in the
[TfReverb technical report](docs/TfReverb-technical-report.md).

### Room Reverb

Room Reverb combines a rectangular three-dimensional image-source model with a
16-line, two-stage paraunitary velvet feedback network. SIZE derives all three
room dimensions. The main delays and both velvet delay banks scale together
from the resulting mean free time, while DECAY remains an independent RT60.
The geometric response provides directional, frequency-dependent early reflections;
the calibrated late field develops underneath it and decays to independent stereo
outputs. FIR construction and room movement run on a rate-limited background worker,
while the audio thread uses allocation-free prepared convolution banks.

SIZE defines the room's floor scale and derives a plausible ceiling height;
ASPECT reshapes the floor at constant area. The top-view placement pad moves up
to eight numbered amber sources and the blue listener in two dimensions. Its asymmetric listener
default avoids the coincident reflection paths produced by the exact room centre;
it remains draggable because no one receiver position is optimal for every source
and room shape. Only channels present on the polyphonic AUDIO input are shown;
one unedited source is centred, while additional unedited channels spread
progressively across the room width. Manually placed positions remain stored
when disconnected. Source and listener heights use
the calibrated factory values. PRE DELAY is shared by the complete wet response. The engine measures
each source's four-band ER handoff and automatically chooses its late send. It
renders every source's immediate direct sound with constant-power azimuth
panning and a bounded room-normalized distance gain. ER / TAIL balances from
the geometric response through the automatically inferred handoff to the late
field; its centre preserves the inferred blend. Direct sound is always positioned. LOW CUT and
HIGH CUT filter the wet output, and DECAY,
DAMPING, DIFFUSE, MOD, SHIMMER, STEREO WIDTH, MIX, and LEVEL control the
remaining response. The late tail sums handoff-aligned sources into one fixed
normalized 16-line FDN input vector and uses a fixed orthonormal stereo decoder;
source/listener geometry may set handoff gain and onset but never steers the
late stereo field. MOD applies slow, sample-rate-invariant random motion to the
main and velvet delays; its subtle default prevents a completely static tail,
while the upper range reaches audible, lush movement. SHIMMER sends four
orthogonal late-field buses through independently seeded eight-grain octave
shifters and returns them at a bounded gain inside the production velvet loop.
This creates successive octave bloom without adding a direct shifted layer or
touching the early reflections and dry signal.
The CV inputs cover PRE DELAY, DECAY, DAMPING, ER / TAIL, MIX, and LISTENER
X/Y. Listener CV is intended for gestures or slow LFOs: direct placement follows
the smoothed control while the background worker transitions the corresponding
early-reflection scene.
The factory scene is a medium hall with a 2.30 second decay and 12 ms wet
pre-delay. Its wet filters default to 20 Hz and 12 kHz. STEREO WIDTH is an output-image control—0%
collapses the complete wet field to mono, 100% preserves its native stereo image,
and values up to 150% exaggerate the side component—whereas ASPECT changes the
actual room width/depth geometry and therefore its reflection times.
The right-click **Presets** menu provides the default Medium Hall, a restrained
Small Room ambience, and a filtered Superlush synth texture. Applying a preset
preserves the module's output LEVEL.

The [two-source Room Reverb test patch](test-room-reverb-two-sources.vcv) uses
a slow two-saw/sub bass and a descending folded arpeggio, each driven by its
own Prog Sequencer and shaped by a 4072 voice. Scene Pack 4 preserves them as
independent left/right room-plan sources for auditioning the spatial direct,
early, and late fields with a musical input. The patch uses Superlush with a
3 kHz wet high cut; arpeggio CV1 drives an 18-beat bipolar cutoff triangle
across the four-beat note phrase, and CV2 sends a per-note AD envelope to the
wavefolder Fold input. One Transport supplies Clock, Run, and Reset
to both sequencers so all transport operations act on them together.

### Transport

Transport is the shared master clock for TriggerFish sequencers. TEMPO covers
30 to 300 BPM. The unsnapped knob supports Rack's fine adjustment and the
numeric display shows hundredths of a BPM; clicking the display opens Rack's
parameter menu for direct value entry. CLOCK has one fixed 24 PPQN rate:
each quarter-note beat contains 24 one-millisecond pulses. The fixed rate gives
a newly started sequencer a
usable period measurement after 1/24 of a beat, including at slow tempos, and
supports straight and triplet subdivisions with the same clock setting in every
patch. The phase accumulator runs at the audio sample rate.

The four buttons have distinct transport meanings:

- RESTART pulses Reset, waits 1.5 ms so downstream modules see it, then raises
  Run and emits the first Clock pulse at beat zero;
- PAUSE lowers Run and freezes the master-clock phase;
- PLAY raises Run and continues from the frozen phase; and
- STOP lowers Run, pulses Reset, and returns the phase to beat zero.

CLOCK is a 10 V pulse stream, RUN is 10 V only while playing, and RESET is a
one-millisecond 10 V pulse. Fan all three outputs from one Transport to every
Prog Sequencer that belongs to the same song. This gives every sequencer the
same pulse edges and makes restart, pause, play, and stop sample-aligned.

### Prog Sequencer

[Complete Prog Sequencer manual and language reference](docs/TfProgSequencer-reference.md)

Prog Sequencer is a 24 PPQN externally clocked live-coding sequencer with
up to 16 polyphonic V/OCT, Gate, Trigger, Velocity, and Accent channels, plus
three monophonic CV lanes. Edit the program
directly on the module. Ctrl+`.` compiles the complete document; Ctrl+Enter
evaluates the top-level statement containing the selection or current line.
Ctrl+Shift+`.` compiles the document and restarts it on the next quarter beat;
Ctrl+Space toggles local transport. With RUN connected directly to a
TriggerFish Transport, Ctrl+Shift+Space pauses or plays the shared Transport,
Ctrl+Shift+R restarts it, and Ctrl+Shift+Backspace stops it. Ctrl+`/` toggles
line comments and Ctrl+D duplicates the current line or selection without
evaluating the draft. Ctrl+Z undoes editor changes; Ctrl+Shift+Z and Ctrl+Y
redo them.
Other edited statements remain inactive until separately evaluated; unrelated
malformed draft text does not block a complete valid selected statement. A
successful edit swaps in on the next quarter beat while preserving the
arrangement phase and the lane cursors of sequences whose names still exist.
Successfully executed code flashes in the editor. The final executed
`play NAME` statement selects the looping sequence or arrangement. Diagnostics
wrap in the status strip, and the last valid program keeps playing after an
error.

While it plays, a soft highlight and progress underline remain on the active
note, rest, or tie. The RUN light flashes on each quarter beat, and each
active CV lane can show a six-second scrolling trace with a fixed -5 V to +5 V
range in the free space to the right of its source line. A CV lane can generate
an AD, AR, or ADSR envelope on its own, or add one to a sequenced CV contour;
the trace follows the final combined voltage.

Its language is inspired by TidalCycles, Gibber, and other live-coding and
pattern-sequencing systems. It is intended to make musical ideas quick to type,
combine, and reshape during live coding or conventional sequencing.

Prog Sequencer is currently in beta. Its syntax, semantics, and feature set are
subject to change as musical workflows are refined.

Double-click selects a word and triple-click selects a complete row. Rack
requires every module to be exactly one 3U row high, but the module context
menu's Editor section offers 22, 30, and 38 HP widths and the complete heatmap
selection. New modules default to 30 HP, and the chosen width is saved with the
patch. The Examples section can load the Acid bassline, Slow bassline, and
Descending arpeggio programs into the editor. Loading an example participates
in Rack undo; evaluate it with Ctrl+`.` when ready. The editor uses a thin outer
margin and a compact right-side I/O strip at every width.

```text
riff = sequence {
  subdiv 16n
  tonic D@3
  scale dorian
  glide 32n
  notes ^1 1!2 x1 [3 4] >5{stacc} ~ 6{quiet} 7{ten}
}

fill = sequence {
  subdiv 16n
  tonic D@3
  scale harmonic_minor
  notes [5 6 7] ^8*2 ~ V7
}

song = riff * 3 + fill
seed 42
play song
```

`riff` and `fill` are independently defined sequences with their own musical
settings and lengths. The `song` assignment arranges complete passes of both;
an arrangement can freely mix independent sections with derived variations.

`subdiv` uses standard musical note values: `4n` is a quarter note, `8n` an
eighth, and `16n` a sixteenth. Append `d` for dotted values and `t` for
triplets, as in `8nd` and `8nt`. Other tempo-relative controls accept the same
notation or bare quarter-note beats.

Articulation normally lives on the note itself. `^3` and `^^3` assert Accent
and raise Velocity to their two built-in emphasis levels. Velocity is a
continuous performance value; Accent remains separate so voices such as the
303 can change timbre and envelope behaviour rather than merely loudness.
`x3` is a short quiet ghost, while `{quiet}` only halves Velocity. `{stacc}`
and `{ten}` select quarter-span and 95-percent Gate without changing the event
span. `_` ties the preceding pitch, `~` rests, and `>3` slides from its
predecessor. `[8 9]` subdivides one parent span, `7*3` ratchets within a span,
`7!3` repeats across three cells, and
`7(3,8,1)` distributes that event over a rotated eight-cell Euclidean rhythm.
Inactive Euclidean cells are rests. `7_`, `7_3`, and `7.` make a note doubled,
three times as long, and dotted respectively. Sparse numerical lanes use `.`
as a typed no-op. A velocity such as `.5` means one half; CV lane numbers are
volts. There is no separate numerical Accent lane. A semicolon between pattern
terms is a visual separator with no effect on timing.

`3?0.35` plays degree 3 with 35 percent probability and otherwise produces a
rest for the same span. `3??0.35` instead omits the event and its span on a
miss, shortening that Notes pass. Free numerical lanes advance only for
surviving events and continue cycling when Notes wraps. Lanes using `...` align
to the beginning, end, or both edges of each surviving Notes pass.

`$` produces a seeded random note from the active scale. `${1,8}` selects an
inclusive scale-degree range, `$n{4,1.5}` uses a normal distribution in degree
space, and `$c{0,11}` selects unquantized chromatic semitone offsets. Numerical
lanes use forms such as `velocity $u{.4,.9}` and `cv1 $n{0,2}`.

Scale degrees refer to the selected scale, so degree `3` in `scale minor` is
already the minor third; `b3` lowers that result by one more semitone. Degree
numbers continue according to scale cardinality: `8` is the next tonic in a
seven-note scale, `6` in a pentatonic scale, and `9` in an octatonic scale.
`@` gives a pitch or chord an unambiguous absolute register. `D@4` is the
individual note D4, while `D7` is a dominant-seventh chord and `D7@3` roots
that chord in octave 3. A trailing apostrophe raises a pitch or chord by one
relative octave and a comma lowers it: `1'`, `1''`, `1,`, and `Cm7'`. The
marks compose, while signed forms such as `C@-1` still name absolute octaves.
An optional `octave` lane selects absolute octaves and loops independently like
every other control lane, so `octave 3 3 4` can displace against an eight-note
melody. Named sections
concatenate with `+` and repeat with `*`. `1!4` repeats a note in its lane, and
`.!3` repeats a sparse default.
`gate .5` holds Gate for half the event; only `_` is a semantic tie.
`glide .8` slews a `>` target for .8 quarter-note beats, capped by that
target event's span.

`//` begins a line comment and can quickly shorten a lane while retaining the
unused material: `notes 1 2 3 // 4 5 6 7 8` compiles as the three-note pattern
`1 2 3`. Whole comment lines and trailing comments on commands are also valid.
Comment text is shown in a dimmer shade of the display's existing colour
scheme. `#` remains available for sharp degrees such as `#4`.

Existing sequences can be varied without copying lanes:

```text
iv = acid |> modulate_degree 4 |> octave -1
v  = acid |> modulate_degree 5 |> octave -1 |> scale major
song = acid * 8 + iv * 4 + v * 4
```

`modulate_degree 4` moves the tonal centre by the selected scale's fourth and
`modulate_degree 5` by its fifth, preserving the riff's chromatic shape;
negative degree numbers move downward. `shift_degree 4` instead remaps every
scale degree four indices upward (tonic to fifth), staying in the current key;
this diatonic operation may change the semitone shape of the riff. In harmonic
minor, for example, that in-key shift produces the major V without a separate
`scale major`. `transpose_semitone` remains the separate operation for exact
chromatic movement.
Supported scale names include `harmonic_minor`, `major_pentatonic`,
`minor_pentatonic`, `octatonic_whole_half`, and `octatonic_half_whole`.

Parentheses are simultaneous scale-degree voicings. Conventional jazz symbols
and Roman-degree chords are accepted directly; uppercase Roman degrees imply
major quality, lowercase implies minor, and degrees stop at VII:

```text
notes I i iim7 bVII (1 b3 5) Cm7 D7 Bbm7b5@3 / D@2 Cmaj9 C7#9
```

Jazz chord syntax supports major, minor, diminished, augmented, suspended,
power, sixth, seventh, ninth, eleventh, and thirteenth chords, plus `add` and
`b`/`#` alterations. Close-position voices are emitted as Rack polyphonic
channels on Pitch, Gate, Trigger, Velocity, and Accent. A slash bass is channel
1 and is placed below the chord unless its octave is explicit. Rack's native
16-channel cable capacity is the chord-size limit.
Inside an explicit voicing, bare named pitches are unambiguous, so
`(C E G)@4` shares octave 4 across the three tones while
`(C@3 E@4 G@4)` remains available for a spread voicing.

Whole-sequence timing transformations are explicitly quantified:

```text
groove = riff
  |> fast 2
  |> swing .58 32n
  |> late 128n
  |> early random 6ms
```

`slow 3/2` is also valid. Swing `.5` is straight; larger values give the first
member of each equal subdivision pair more of the pair. With no second argument
the event subdivision supplies the grid. Explicit grids and other
tempo-relative times accept note values such as `8n`, `16nd`, and `8nt`, while
bare numbers remain quarter-note beats. `early` and `late` also accept `ms`,
and `random AMOUNT` chooses a deterministic
amount from zero to the stated maximum. Subdivision density can already vary
with note groups such as `[1 2] [3 4 5]`. When a separate rhythmic cycle is
deliberate, independent `ratchet` and `duration` lanes are also available. A
patternable `offset` lane accepts signed beats, note values, or `ms`; negative
values are early and positive values late. Its optional numeric `rate`, as in
`offset -10ms!2 8ms |> rate 1/2`, changes only that lane's score-time phase.
CV1-CV3 use the same rate and sparse-lane rules and support `step`, `linear`,
`smooth`, and `power P` interpolation. Ellipsis-aligned CV is stepped in this
version; continuous modes use free score-time lanes.

During playback, heatmap cursors show the active arrangement term and every
independently advancing lane. See the
[complete Prog Sequencer reference](docs/TfProgSequencer-reference.md) for all
syntax, additional musical examples, and current limits.

## Contributing

Issues, pull requests, and suggestions are welcome.

## License

TriggerFish Elements 2.4.0 and later are licensed under the
[GNU General Public License v3.0 or later](LICENSE.txt). Releases through
2.3.1 remain available under the BSD 3-Clause terms under which they were
published. See [NOTICE.md](NOTICE.md) for the license transition and retained
notices; vendored dependencies retain their own licenses.

## Development

See [DEVELOPMENT.md](DEVELOPMENT.md) for Windows and Linux setup and the build,
test, package, install, panel-rendering, and Rack workflows.
