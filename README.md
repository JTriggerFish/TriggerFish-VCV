
# TriggerFish Elements plugins

<img src="doc/modules.png" width="600">

[![CI](https://github.com/JTriggerFish/TriggerFish-VCV/actions/workflows/ci.yml/badge.svg)](https://github.com/JTriggerFish/TriggerFish-VCV/actions/workflows/ci.yml)

## Modules
- [Slop and Slop 4](#slop-and-slop-4)
- [VDPO](#vdpo)
- [VCA](#vca)
- [Diode Ladder Filter](#diode-ladder-filter)


### Slop and Slop 4
Slop and Slop4 are utilities to add drift and hum to V/oct signals in order to add pleasant detuning to VCOs.

Slop can add linear detuning ( Hz mode ) or proportional detuning ( cents mode ).

In Slop4, the common detuning is in cents ( i.e. proportional ) and the individual detunings are linear to give a more stable and pleasant beating accross octaves.

Hum adds 60hz frequency modulation to the V/Oct signals to replicate the power supply bleed and cross modulation in analog oscillators.


### VDPO
VDPO or Van Der Pol oscillator is based on the classic differential equation [wikipedia](http://en.wikipedia.org/wiki/Van_der_Pol_oscillator)

While it can self oscillate, best results are obtained by feeding it another oscillator at the input and playing with the self-freq, damping and input level to go from harmonic to inharmonic and chaotic.

**self-freq** controls the natural resonance frequency of the oscillator

**damping** controls the non linear damping of the oscillator. Note that it will affect the tuning as well as the harmonic and inharmonic content of the output waveform.

**in** controls the level of the (forcing) input which will force the oscillator to run at the same frequency as the input signal.

The higher the level the more the VDPO will follow the input, however if self-freq is too low the oscillator is too slow to follow it and all sorts of fun stuff happens.

The stiff differential equation is solved with a structure-aware, prewarped
split integrator and 4x oversampling. This improves pitch tracking across the
audio range while using adaptive work only for high damping at high frequency.

## VCA
TriggerFish VCA is an analog modelled VCA that is loosely based on the minimoog's circuit.
Just like the original it includes 3 non linearities, one on the audio, one on the CV and one on the output.

The input non linearities are tanh-like but with a 1 pole feedback loop, resulting in some amount of slew limiting. The output nonlinearity is also a tanh but with no feedback loop, it also serves to limit the output to +-12v.

**Drive** controls the level-compensated input drive. Its extended range mostly increases saturation and nonlinear slew rather than acting as another output-volume control.
At very low levels more pink noise will be heard on the output, and with the knob fully counterclockwise the input will be cut out.

**lin** and **exp** are controls for the gain of the linear and exponential cv inputs. 
Exponential is more snappy but linear is typically better for normal enveloppe signals. Higher gains will result in more saturation and slew limiting of the CV.

For standard useage with enveloppes ADSR or AR enveloppes I recommend using the linear input.

**cv curve** controls the curve of the exponential input - higher is more exponential and snappier, lower is more linear.

**output** controls the output level, higher values will cause saturation as the level approach +-12v

**bleed** will send part of the input CV to the output to make it more clicky.

The model use antialiased integration for the nonlinearities and the whole module runs at 2x oversampling for low aliasing. 
Because of this the CPU useage is relatively high.

## Diode Ladder Filter

The complete circuit analysis, equations, numerical method, calibration, and
reference measurements are documented in the
[TfDiodeLadderFilter technical report](docs/TfDiodeLadderFilter-technical-report.md).

`TfDiodeLadderFilter` is a circuit-structured four-stage diode ladder inspired
by the Roland TB-303 filter and its extended high-resonance, drive, filter-FM,
and bass modifications. The unbuffered stages are solved as one coupled
nonlinear system; the AC-coupled resonance path is included in the implicit
solve rather than approximated by a separate output high-pass filter.

**Cutoff** spans 10 Hz to 20 kHz, defaults to the service-note 500 Hz centre
calibration, and has dedicated 1 V/octave and attenuverted exponential CV inputs.
The exponential CV attenuverter defaults to 53.2%, so a standard 0--10 V
envelope sweeps the default cutoff to 20 kHz; at 100% it follows the usual
1 V/octave Rack mapping.
**FM** is an AC-coupled linear-Hz input intended for audio-rate modulation. At
full amount, a nominal +/-5 V Rack signal sweeps the cutoff by +/-1 kHz.
**Resonance** has an attenuverted CV input and a Stock/High switch; High mode can
self-oscillate in the upper part of the knob travel at suitable mid/high cutoff
settings. Resonance-dependent output makeup prevents the source from vanishing,
while retaining the level thinning and drive-dependent resonance quenching of
the circuit.

**Drive** ranges from silence through a marked 0 dB stock point to 66.6 times
stock level. **Bass** continuously morphs the calibrated stock coupling response
to the extended response while retaining the circuit's internal
frequency-dependent feedback behaviour.

Input drive is calibrated from the stock oscillator and resistor divider in the
Roland schematic. Output gain is calibrated independently so that a stock-driven,
open filter keeps a normal +/-5 V oscillator in the same practical Rack voltage
range; it is not an algebraic reciprocal of the nonlinear input scaling.

The module supports Rack polyphony and uses 4x oversampling by default. A 2x
mode is available in the context menu for dense patches where lower CPU use is
more important than extreme high-frequency nonlinear accuracy.

`GATE` drives separate TB-303-style filter and volume envelopes. `ACCENT` adds
the short main envelope to both the filter sweep and BA662 control-current
paths; the filter sweep retains capacitor memory across repeated accents.
Tied/high Gates do not retrigger, so pitch slides can preserve the current
articulation. The panel's snapped **Accent Sweep** control selects
Off/Fast/Normal/Slow behavior; Normal is the circuit-derived stock default.
The context menu selects stock or Devil Fish gate behavior. Fast and Slow
follow the published Devil Fish behavior descriptions, since their changed
component values have not been published.

`LP OUT` remains the raw filter output for ordinary modular use. `VCA OUT`
passes the same signal through a reduced BA662-style OTA model and its output
coupling. `VCA CV` is normalled to the internal volume envelope, but a patched
0--10 V signal replaces that base control; its attenuator defaults to 100%.
The accent current remains additive, so turning `ACCENT` down gives a fully
independent externally controlled VCA.



## Contributing

Issues, pull requests and suggestions for improvement are all very welcome.

## Development

See [DEVELOPMENT.md](DEVELOPMENT.md) for the Windows toolchain and the
one-command build, test, package, install, and Rack workflows. The old
machine-specific Visual Studio/Python wrapper has been replaced by CMake,
pybind11, scikit-build-core, and uv.
