# TriggerFish Elements

Circuit-inspired sound generators, processors, and pitch utilities for VCV Rack 2.

[![CI](https://github.com/JTriggerFish/TriggerFish-VCV/actions/workflows/ci.yml/badge.svg)](https://github.com/JTriggerFish/TriggerFish-VCV/actions/workflows/ci.yml)

[Release notes for 2.3.0](docs/releases/2.3.0.md)

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
    <td colspan="3" align="center"><a href="#4072-voice-core"><img src="doc/Tf4072VoiceCore.png" height="260" alt="4072 Voice Core module"><br><strong>4072 Voice Core</strong></a></td>
  </tr>
</table>

## Polyphony

| Module | Handling | Output voice count |
| --- | --- | --- |
| Slop | Monophonic; reads channel 1 of `1V/OCT` | One |
| Slop 4 | Four independent monophonic signal paths | One per output jack |
| VDPO | Monophonic; reads channel 1 of each input | One |
| VCA | Monophonic; reads channel 1 of each input | One |
| 303 Oscillator | Fully polyphonic, with independent DSP state per voice; mono inputs are broadcast | Widest connected input, up to 16 |
| 303 Voice Core | Fully polyphonic, with independent filter, envelope, accent, and VCA state; mono controls are broadcast | Channel count of `IN`, up to 16 |
| 4072 Voice Core | Fully polyphonic, with independent filter, dual-envelope, and VCA state; mono controls are broadcast | Widest connected input, up to 16 |

## Modules

### Slop and Slop 4

Slop adds slow pitch drift and 60 Hz power-supply hum to a 1 V/octave signal.
The drift can be proportional, measured in cents, or linear, measured in hertz.
The tracking control also allows small oscillator-tracking errors to be modelled.

Slop 4 applies shared proportional drift and hum to four independent monophonic
paths, then adds a separate linear drift to each voice. This produces coherent
ensemble motion with stable beating across the keyboard. Each path has its own
tracking trim.

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

303 Oscillator combines a band-limited saw core with a circuit-informed model of
the TB-303 Q8 saw-to-square shaper. The module supports polyphony,
voltage-controlled waveform morphing, hard sync, extended tuning, and two FM
responses.

- **Octave** selects an integer offset from −3 to +3 octaves; **Tune** covers
  ±7 semitones.
- `1V/OCT` is the pitch input. `SLIDE` enables the internal pitch glide while
  high, and **SL.Time** sets its nominal duration from 2 to 360 ms.
- **Shape** shifts the square-wave switching point. **Wave** morphs continuously
  from saw to square without a level notch around the middle.
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

### 303 Voice Core

303 Voice Core is the back end of a TB-303-style voice. It models the four-stage
diode ladder, coupling and resonance networks, filter and volume envelopes,
accent path, and BA662-style OTA VCA. `LP OUT` exposes the filter directly;
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

The `IN` cable sets the filter's polyphonic voice count. Each audio channel has
independent filter, envelope, accent, and VCA state; mono control inputs are
shared across those voices.

The filter uses 4x oversampling by default, with a 2x context-menu option. The
context menu also selects stock or Devil Fish volume-envelope timing. Circuit
analysis, equations, calibration, and numerical validation are collected in the
[303 Voice Core technical report](docs/Tf303VoiceCore-technical-report.md).

### 4072 Voice Core

4072 Voice Core combines a circuit-scaled ARP 4072 four-pole low-pass filter,
two switchable ADSR/AR envelope generators, and an ARP 4019-style discrete VCA.
The filter's outer differential pair preserves the original unequal audio and
resonance-return levels, while the four locally fed-back filter stages operate
around their circuit bias points.

- **Cutoff** and `1V/OCT` set the filter frequency. `CUTOFF CV`, `RES CV`, and
  the internal filter envelope have independent amount controls.
- **Input** covers attenuation through +24 dB of filter drive. **Range** selects
  the stock 4072 high-frequency limit or an extended 20 kHz range.
- The two A/D/S/R slider rows generate 0–10 V envelopes. Each row can be switched
  to AR; filter ADSR and amplifier AR are the defaults. **Curve** varies both
  envelope shapes around the original circuit response.
- `ENV CV` replaces the internal filter-envelope normalization. `LIN CV`
  replaces the amplifier-envelope normalization, while `EXP CV` accesses the
  4019's 10 dB/V exponential control path.
- `VCA IN` replaces the normalled filter-to-VCA audio connection. `LP OUT`,
  `VCA OUT`, and both envelope outputs expose each stage independently.

The complete filter/VCA path runs at 4x sample rate by default. A 2x mode is
available from the context menu for larger polyphonic patches.
Circuit analysis, equations, calibration, and numerical validation are in the
[4072 Voice Core technical report](docs/Tf4072VoiceCore-technical-report.md).

## Contributing

Issues, pull requests, and suggestions are welcome.

## Development

See [DEVELOPMENT.md](DEVELOPMENT.md) for Windows and Linux setup and the build,
test, package, install, panel-rendering, and Rack workflows.
