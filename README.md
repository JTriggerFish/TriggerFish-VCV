
# TriggerFish Elements plugins

This packaga currently includes 4 modules.

Slop and Slop4 are utilities to add drift and hum to V/oct signals in order to add pleasant detuning to VCOs.
Slop can add linear detuning ( Hz mode ) or proportional detuning ( cents mode ).
In Slop4, the common detuning is in cents and the individual detuning is linear which gives a more stable and pleasant beating accross octaves.

VDPO or Van Der Pol oscillator is based on the classic differential equation.
While it can self oscillate, best results are obtained by feeding it another oscillator at the input and playing with the self-freq, damping and input level to go from harmonic to inharmonic and chaotic.
Note that this modules solves a stiff differential equation in real time using backward differentiation and 4x oversampling and is thus quite CPU heavy despite using vectorized arithmetic.

VCA is an analog modelled VCA that is loosely based on the minimoog VCA.
Just like the original it includes 3 non linearities, once on the audio, one on the CV and one on the output.
These use antialiased integration and the whole module 2x oversampling for low aliasing. Again CPU useage is relatively high.

See https://vcvrack.com/manual/PluginDevelopmentTutorial.html for a development tutorial.

## Contributing

I welcome Issues and Pull Requests to this repository if you have suggestions for improvement.

