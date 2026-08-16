# 303 Oscillator technical report

<p align="center"><img src="../doc/Tf303Oscillator.png" height="520" alt="Tf303Oscillator module"></p>

The [module guide](../README.md#303-oscillator) describes the panel controls and
patching interface.

## Overview

`Tf303Oscillator` follows the signal path of the TB-303 oscillator: a sawtooth
core feeds the Q8 transistor stage that produces the square waveform. The model
captures the frequency-dependent pulse width and asymmetric edge shaping of that
stage, together with the original pitch-slide time constant. VCV extensions add
polyphony, continuous wave morphing, shape CV, through-zero linear FM, hard sync,
and a wider pitch range.

The audio path runs at 4x sample rate by default. A 2x option is available from
the module context menu for lower CPU use.

## Pitch and slide

The octave switch, fine tuning, and post-slide pitch form the pitch value $p$ in
octaves relative to C4. Exponential FM adds to this value before conversion:

$$
f = 261.625565\,2^{p + 0.2a_\mathrm{FM}v_\mathrm{FM}} .
$$

Here $a_\mathrm{FM}$ is the bipolar FM amount and $v_\mathrm{FM}$ is the input
voltage. At full depth, 5 V changes the exponential pitch by one octave.

While `SLIDE` is high, the pitch state approaches the input pitch through a
matched one-pole update:

$$
p_{n+1}=p_n+\alpha(p_\mathrm{target}-p_n),
$$

$$
\alpha=-\operatorname{expm1}\!\left(-\frac{1}{f_s\tau}\right),
\qquad
\tau=T_\mathrm{slide}\frac{22}{60}.
$$

The stock 60 ms setting therefore uses the 22 ms pitch-CV time constant associated
with the TB-303 slide circuit. Releasing `SLIDE` transfers the pitch immediately.
`CV OUT` exposes this post-slide state before octave and fine-tune offsets.

In through-zero mode, FM is applied in hertz:

$$
f = 261.625565\,2^p + 200a_\mathrm{FM}v_\mathrm{FM}.
$$

Negative frequency reverses phase accumulation. A smooth knee begins at 40% of
the host sample rate and approaches a 45% limit, keeping extreme modulation
finite without introducing an abrupt frequency clamp.

## Saw oscillator

The saw phase advances at the internal oversampled rate. Ordinary phase wraps
use a short polynomial BLEP correction, which is efficient for the predictable
once-per-cycle discontinuity. The Rack-facing saw is scaled to approximately
10 V peak-to-peak.

Hard-sync edges can arrive anywhere inside a host sample. A Schmitt trigger
estimates the crossing fraction by linear interpolation, then maps that event to
the corresponding oversampled segment. The oscillator advances to the crossing,
resets its phase, advances through the remainder of the segment, and inserts a
minimum-phase band-limited step correction with the measured discontinuity
magnitude.

The shared minBLEP kernel is built from an eight-zero-crossing Blackman–Harris
windowed sinc. Cepstral minimum-phase reconstruction moves most of its energy near
the event, and a table oversampling factor of 32 supports fractional lookup. The
resulting correction occupies 16 internal samples.

## Q8 square-wave shaper

The square shaper is a reduced dynamic model of the Q8 common-emitter transistor
stage. Its switching threshold varies with frequency to reproduce the measured
trend of the hardware waveform: a broad positive pulse at very low frequency,
approximately equal duty near 85 Hz, and a narrower positive pulse at higher
audio frequencies. This behaviour follows Richie Burnett's
[measurements of an original TB-303](https://synth-diy.org/pipermail/synth-diy/2015-October/150013.html).

With frequency $f$ and $o=\log_2(f/85)$, the stock threshold is

$$
q(f)=
\begin{cases}
-0.1359o, & o < 0,\\
-0.0281o, & o \ge 0,
\end{cases}
$$

limited to $[-0.15,0.45]$. The Shape control adds $0.55s$ and extends the final
threshold over $[-0.82,0.82]$.

The transistor transition uses

$$
c=\frac{1}{2}\left[1+\tanh\!\left(
\frac{q_\mathrm{effective}-x_\mathrm{saw}}{0.055}\right)\right].
$$

An emitter-bias state follows $c$ with the R104/C11 time constant of 2.2 ms and
feeds a small history-dependent shift back into the threshold. The collector
state follows $2c-1$ with a 100 µs rising time and 22 µs falling time. These two
states produce the characteristic top droop and unequal edge shapes. At extended
VCV pitches the collector time is smoothly bounded to 8% of a cycle so the
square response remains useful.

The output removes the threshold-dependent DC component. Its polarity is aligned
with the saw before morphing, and its gain reflects the hardware’s lower square
level relative to the saw. Linear interpolation then provides a continuous Wave
control without a cancellation notch near the centre.

## Oversampling and control interpolation

V/oct pitch, logarithmic Slide Time, FM, Shape, and Wave CV are interpolated
into the oversampled domain. Pitch is converted to frequency after
reconstruction, while Shape and Wave limits are applied at the internal rate.
The saw core, sync correction, Q8 shaper, and waveform mix are evaluated there,
followed by separate seventh-order decimators for the saw, square, and mixed
signals.
Keeping the nonlinear shaper inside the oversampled path suppresses aliases from
its steep transitions. The 4x default gives additional margin for high-pitch FM
and sync sounds. The 2x mode approximately halves the oscillator's DSP cost for
dense polyphonic patches.

## Controls and ports

| Control or port | Function |
| --- | --- |
| Octave | Integer pitch offset from −3 to +3 octaves |
| Tune | Continuous offset over ±7 semitones |
| SL.Time | Nominal slide duration, 2–360 ms; 60 ms default |
| Shape | Square switching bias, extended around the stock response |
| Wave | Continuous saw-to-square morph |
| FM mode | Exponential or linear through-zero response |
| FM / SL.Time CV / Shape CV / Wave CV | Bipolar modulation amounts |
| 1V/OCT | Primary pitch input |
| Slide | Gate-controlled pitch glide |
| Sync | Fractionally timed rising-edge hard sync |
| CV OUT | Post-slide pitch CV |
| OUT | Mixed audio output |

All signal and modulation inputs support Rack polyphony. The widest connected
input determines the channel count, and mono inputs are shared across channels.

## Validation

The native and Python tests cover 1 V/octave tracking, the 22 ms stock slide
constant, waveform endpoints and morph level, frequency-dependent square duty,
Shape range, through-zero reversal, fractional sync convergence against a
16-times-rate reference, and bounded output under combined sync and through-zero
FM. Oversampling comparisons can be reproduced with:

```console
uv run pytest tests/python/test_tb303_oscillator.py
uv run python tests/python/benchmark_tb303_oscillator.py --seconds 2 --cpu-samples 960000 --cpu-repeats 9
```

The published run uses a two-second analysis window. It compares each production
mode at 48 kHz with a 768 kHz render followed by high-rejection FIR decimation.
The spectral metric includes aliasing and passband-magnitude differences;
more-negative values are closer to the reference.

| Scenario | 2x | 4x | 4x improvement |
| --- | ---: | ---: | ---: |
| 1 kHz square | −42.66 dB | −43.30 dB | 0.64 dB |
| 6 kHz saw | −17.21 dB | −40.48 dB | 23.27 dB |
| 6 kHz square | −21.65 dB | −25.68 dB | 4.03 dB |
| 12 kHz square | −14.58 dB | −19.59 dB | 5.01 dB |
| 6.1 kHz saw, 997 Hz sync | −25.49 dB | −28.07 dB | 2.58 dB |
| 1.5 kHz saw, exponential FM | −35.48 dB | −36.55 dB | 1.07 dB |
| 1.2 kHz saw, through-zero FM | −34.59 dB | −36.60 dB | 2.01 dB |
| Through-zero FM, sync, and 50% morph | −19.69 dB | −26.34 dB | 6.65 dB |

On an AMD Ryzen 9 PRO 8945HS with a GCC 16.2 release build, the median cost was
175–183 ns per sample at 2x and 342–349 ns at 4x. The 4x path therefore used
1.91–1.95 times the CPU in this harness. Absolute timings include the Python
binding and output-buffer write and should be treated as machine-specific; the
ratio is the useful comparison.

The implementation is in
[`src/models/Tb303Oscillator.hpp`](../src/models/Tb303Oscillator.hpp), with the
generic phase, minBLEP, and fractional-trigger components in
[`src/tfdsp/oscillator.hpp`](../src/tfdsp/oscillator.hpp),
[`src/tfdsp/minblep.hpp`](../src/tfdsp/minblep.hpp), and
[`src/tfdsp/control.hpp`](../src/tfdsp/control.hpp).

## References

1. Roland Corporation, [TB-303 Service Notes](https://www.synfo.nl/servicemanuals/Roland/ROLAND_TB-303_SERVICE_NOTES.pdf): oscillator, waveform switch, pitch CV, and slide circuitry.
2. Robin Whittle, [TB-303 unique accent and envelope characteristics](https://www.firstpr.com.au/rwi/dfish/303-unique.html): note, slide, and articulation behaviour.
3. Robin Whittle, [Devil Fish user manual](https://www.firstpr.com.au/rwi/dfish/Devil-Fish-Manual.pdf): extended slide controls and operating ranges.
4. VCV Rack, [Voltage Standards](https://vcvrack.com/manual/VoltageStandards): Rack pitch, audio, gate, and modulation conventions.
5. Richie Burnett, [Roland TB-303 “square” wave variation](https://synth-diy.org/pipermail/synth-diy/2015-October/150013.html): measured pitch-dependent Q8 waveform shape and the approximately 85 Hz equal-duty point.
