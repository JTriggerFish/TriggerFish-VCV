# TfDiodeLadderFilter technical report

## Overview

`TfDiodeLadderFilter` models the Roland TB-303 diode ladder, its surrounding
AC-coupling and resonance network, and selected Devil Fish extensions. The
module adds a TB-303-style filter envelope, volume envelope, accent paths, and
a reduced BA662-style OTA VCA. Separate outputs expose the filter and the
post-VCA voice.

The circuit values and transfer functions come from the Roland service notes
and Tim Stinchcombe's analyses. Devil Fish ranges and articulation behaviour
come from Robin Whittle's published documentation. The OTA reduction was
checked against the Open Music Labs BA662 clone populated with manufacturer
models for modern matched transistors. Open303 and Soundpipe were consulted for
software resonance-level comparisons.

## Diode ladder

The ladder has four mutually loading stages. Its first capacitor is 18 nF and
the remaining three are 33 nF. Define

$$
r = \frac{33}{18}, \qquad
c = r^{-1/4} = 0.8593887047640296, \qquad
\Omega = 2\pi f_c c .
$$

Let $x_0,\ldots,x_3$ be the ladder states, $u$ the forward-path signal,
$F(x_3)$ the resonance-path signal, and $k$ the feedback amount. The normalized
nonlinear junction currents are

$$
\begin{aligned}
j_0 &= \tanh\!\left(u-kF(x_3)\right),\\
j_1 &= \tanh(x_0-x_1),\\
j_2 &= \tanh(x_1-x_2),\\
j_3 &= \tanh(x_2-x_3),\\
j_4 &= \tanh(x_3).
\end{aligned}
$$

The state equations are

$$
\begin{aligned}
\dot{x}_0 &= \Omega r(j_0-j_1),\\
\dot{x}_1 &= \Omega (j_1-j_2),\\
\dot{x}_2 &= \Omega (j_2-j_3),\\
\dot{x}_3 &= \Omega (j_3-j_4).
\end{aligned}
$$

Linearization around zero gives the normalized denominator

$$
D(p)=p^4+(r+6)cp^3+(5r+10)c^2p^2
 +(6r+4)c^3p+1,
\qquad p=\frac{s}{2\pi f_c}.
$$

For the idealized ratio $r=2$, this becomes Stinchcombe's polynomial

$$
D(p)=p^4+2^{11/4}p^3+10\sqrt{2}\,p^2+2^{13/4}p+1.
$$

The transition region is shallower than four buffered identical poles, while
the asymptotic stop-band slope remains 24 dB/octave.

## Coupling and resonance network

The input, output, and resonance paths contain several coupling sections. The
complete fitted forward path is

$$
\begin{aligned}
H_f(s)={}&1.06
\frac{s}{s+578.1}
\frac{s}{s+97.5}
\frac{s}{s+38.5}\\
&\times\frac{s+109.9}{s+20.0}
\frac{s+34.0}{s+4.45}.
\end{aligned}
$$

The resonance path is

$$
\begin{aligned}
H_r(s)={}&18.7
\frac{s}{s+578.1}
\frac{s}{s+97.5}
\frac{s}{s+38.5}
\frac{s}{s+20.0}\\
&\times\frac{s+46.5}{s+7.41}
\frac{s+4.40}{s+4.45}.
\end{aligned}
$$

The constants are angular frequencies in radians per second. A common
$(s+7.41)$ factor in the published forward numerator and coupling denominator
has been cancelled. With $L(s)=1/D(s)$ denoting the linearized ladder, the
closed-loop response is

$$
H(s)=\frac{L(s)H_f(s)}{1+kL(s)H_r(s)}.
$$

The stock and extended feedback mappings are

$$
k_{\mathrm{stock}}=0.78q,
\qquad
k_{\mathrm{high}}=2(0.78q),
$$

where $0\le q\le1$ is the panel resonance. High mode follows the published
Devil Fish doubled-feedback range and reaches self-oscillation.

The output applies a resonance-dependent calibration

$$
G_{\mathrm{out}}=2.75\left(1+qm\right),
\qquad
m=\begin{cases}1,&\text{stock},\\2,&\text{high}.
\end{cases}
$$

The fixed factor 2.75 preserves a practical Rack level. The second factor is a
software compensation informed by Open303 and `tbvcf`; it retains part of the
source-level reduction produced by the circuit.

## Bass extension

The Devil Fish modification increases two output-coupling capacitors by a
factor of ten. The software uses a calibrated two-shelf reduction

$$
H_b(s,b)=\left(\frac{s+\omega_b}
 {s+\omega_b10^{-b}}\right)^2,
\qquad
\omega_b=2\pi(24.66),
\qquad 0\le b\le1.
$$

At $b=0$ the shelves cancel. At $b=1$ their poles move down one decade. The
result adds approximately 4 dB at 32 Hz relative to the stock setting. The
24.66 Hz corner and squared-shelf form are fitted software parameters. A 10 ms
smoother prevents abrupt coefficient changes.

## Discrete-time realization

### Analog ratio sections

Each coupling factor

$$
H_a(s)=\frac{s+z}{s+p}
$$

uses a topology-preserving one-pole state. At sample rate $f_s$,

$$
g=\tan\!\left(\frac{p}{2f_s}\right),
\qquad
a=\frac{g}{1+g},
$$

followed by

$$
\begin{aligned}
\ell &= ax+(1-a)s_1,\\
y &= x+\left(\frac{z}{p}-1\right)\ell,\\
s_1' &= 2\ell-s_1.
\end{aligned}
$$

Before updating the state, the section output is affine in its input:

$$
y=\left[1+\left(\frac{z}{p}-1\right)a\right]x
 +\left(\frac{z}{p}-1\right)(1-a)s_1.
$$

Composing these gain and offset terms through the resonance cascade gives

$$
F(x_3)=Ax_3+B.
$$

This expression enters the current nonlinear solve and preserves the
instantaneous resonance-path response.

### Ladder update and Newton solve

The ladder junction differences use an implicit midpoint update. The resonance
signal uses the current-sample affine TPT preview. For the previous and next
state vectors,

$$
x_m=\frac{x_n+x_{n+1}}{2},
$$

Define the dimensionless junction-flow vector

$$
\Phi(x)=
\begin{bmatrix}
r(j_0-j_1) & j_1-j_2 & j_2-j_3 & j_3-j_4
\end{bmatrix}^{\mathsf T}.
$$

The residual is

$$
R(x_{n+1})=x_{n+1}-x_n
-\Gamma\Phi(x_m,x_{3,n+1}),
$$

with the prewarped coefficient

$$
\Gamma=2c\tan\!\left(\frac{\pi f_c}{f_{s,\mathrm{internal}}}\right).
$$

Newton iterations solve

$$
J(x_i)\Delta_i=-R(x_i),
\qquad
x_{i+1}=x_i+d_i\Delta_i.
$$

The analytic $4\times4$ Jacobian uses

$$
\frac{d}{dz}\tanh z=1-\tanh^2z.
$$

Its ladder terms form a tridiagonal matrix, with an additional row-0,
column-3 term from the resonance preview. Gaussian elimination uses partial
pivoting. Updates larger than one normalized state unit are damped. The solve
allows eight iterations and terminates when

$$
\max_i |R_i| < 10^{-11}.
$$

A finite bounded final iterate is retained when the iteration limit is reached,
and a diagnostic counter records the event. Non-finite states or magnitudes
above 100 reset the affected channel.

### Oversampling

The module uses the repository's polyphase IIR half-band resampler, derived
from the Valenzuela--Constantinides elliptic half-band construction. The
seventh-order design runs at 2x or 4x, with 4x as the default.

One interpolator feeds the input network, ladder, Bass correction, and OTA VCA.
Independent decimators give both outputs the same resampling phase; the VCA's
C38 coupling section contributes its intended low-frequency phase response.
VCA CV is interpolated before the OTA stage. The envelope and accent state
machines run at the host rate because their bandwidth is low.

A coherent high-drive sine test measured the following non-harmonic residuals.
The 16x result uses the ladder without an internal resampler followed by a
high-rejection FIR decimator.

| Input | 2x, order 7 | 4x, order 7 | 16x + FIR |
|---:|---:|---:|---:|
| 5 kHz | -18.26 dBc | -23.84 dBc | -37.74 dBc |
| 7 kHz | -14.89 dBc | -19.33 dBc | -36.32 dBc |
| 9 kHz | -10.74 dBc | -16.33 dBc | -28.01 dBc |
| 11 kHz | -8.21 dBc | -13.28 dBc | -32.15 dBc |

The ninth-order repository design was also evaluated. At 4x its residual moved
between 0.3 dB worse and 1.2 dB better across the same cases. Its flatter
passband reduced the 20 kHz round-trip loss from 0.097 dB to below 0.001 dB,
while group delay increased from 3.06 to 3.44 host samples at low frequency and
from 6.34 to 7.80 samples at 20 kHz. The seventh-order design gives the better
balance of alias rejection and phase dispersion for this module.

Changing mode resets the selected filter, resampler, and VCA coupling states.
Gate and accent articulation continue from their current state.

## Input and modulation calibration

The service schematic shows an approximately 6.5 V peak-to-peak saw attenuated
by 220 kohm and 10 kohm before the input pair. For a 10 V peak-to-peak Rack
oscillator,

$$
S_{\mathrm{in}}=
\frac{6.5\,[10/(220+10)]}{2(0.02585)(10)}
\approx0.547
$$

normalized units per Rack volt. Drive applies

$$
u_{\mathrm{normalized}}
=0.547V_{\mathrm{in}}10^{D/20},
$$

where $D$ spans silence through the marked 0 dB stock point to 36.47 dB, or
66.6 times stock.

Cutoff follows the Rack pitch convention

$$
f_c=f_{\mathrm{C4}}2^P.
$$

The knob spans 10 Hz to 20 kHz and defaults to 500 Hz. The V/OCT input adds
directly to $P$. The bipolar exponential CV attenuverter is 1 V/octave at 100%
and defaults to 53.22%, allowing a 0--10 V envelope to move 500 Hz through
5.3219 octaves to 20 kHz.

The linear FM input is AC-coupled at approximately 5 Hz and applies

$$
\Delta f_c=200a_{\mathrm{FM}}V_{\mathrm{FM}}\ \mathrm{Hz}.
$$

At full depth, a $\pm5$ V signal produces a $\pm1$ kHz excursion.

Key calibration values are summarized below.

| Value | Function | Basis |
|---|---|---|
| 18 nF / 33 nF | Ladder capacitor ratio | Roland schematic |
| 0.547 | Rack input scale | Schematic signal and divider estimate |
| 0.78 | Stock feedback scale | Software fit |
| 2x feedback | High resonance | Published Devil Fish range |
| 2.75 | Rack output scale | Level calibration |
| $1+qm$ | Resonance makeup | Software calibration |
| 24.66 Hz, two shelves | Bass response | Calibrated reduction |
| 200 Hz/V | Linear FM depth | Modular range |
| 53.22% | Default exponential CV depth | Modular range |

## Filter and volume envelopes

On a Gate rising edge, the main filter envelope is set to one and decays as

$$
e_{n+1}=e_n\exp\!\left(-\frac{1}{f_s\tau_e}\right).
$$

Accent interpolates geometrically between the normal and accented decay times:

$$
\tau_e=\exp\!\left[(1-a)\ln\tau_{\mathrm{normal}}
 +a\ln\tau_{\mathrm{accent}}\right].
$$

Both controls span 30 ms to 3 s. Filter-envelope pitch modulation is

$$
P_{\mathrm{env}}=
6a_{\mathrm{env}}(e-0.3137)
+2a_{\mathrm{accent}}e_{\mathrm{accent}}.
$$

The 0.3137 pivot represents the original bias network: increasing envelope
depth opens the attack and lowers the tail.

The stock volume envelope uses a 4 ms onset delay, a 3 ms attack, exponential
decay, an 8 ms release hold, and an 8 ms linear close. Devil Fish mode uses a
0.5 ms onset delay and an immediate exponential release with time constant

$$
\tau_{\mathrm{release}}=1.1581186\ \mathrm{ms},
$$

which reaches approximately -60 dB in 8 ms. A continuously high Gate preserves
the current envelope states for tied-note slides.

The first half of the VCA Decay control maps logarithmically from 16 ms to
3.5 s. The second half keeps the 3.5 s decay and raises sustain from zero to
one.

## Accent sweep

The normal accent sweep uses the published 47 kohm, 100 kohm, and 1 uF network.
Its direct gain, capacitor target, attack time, and release time are

$$
G_d=\frac{100}{147},
\qquad
G_c=\frac{100}{247},
$$

$$
\tau_a=(147\ \mathrm{k\Omega}\parallel100\ \mathrm{k\Omega})(1\ \mu\mathrm{F})
=59.5\ \mathrm{ms},
\qquad
\tau_r=100\ \mathrm{ms}.
$$

Resonance interpolates between the direct and capacitor paths, preserving the
capacitor state across repeated accents. Off disables this contribution. Fast
and Slow implement the published Devil Fish performance descriptions; their
switched component values remain unpublished.

The VCA accent branch follows the 47 kohm and 33 nF time constant

$$
\tau_{\mathrm{VCA\ accent}}=(47\ \mathrm{k\Omega})(33\ \mathrm{nF})
=1.551\ \mathrm{ms}.
$$

## BA662-style OTA VCA

The reusable OTA core uses the matched-pair large-signal law

$$
i_{\mathrm{out}}=\eta I_{\mathrm{ABC}}
\tanh\!\left(\frac{v_d}{2V_T}\right),
$$

with

$$
\eta=0.85,
\qquad
V_T=25.85\ \mathrm{mV}.
$$

Its small-signal transconductance is

$$
g_m=\frac{\eta I_{\mathrm{ABC}}}{2V_T}.
$$

The wrapper receives each Rack-scaled oversampled filter value after resonance
makeup and applies

$$
v_d=\sqrt{2}(10^{-3})V_{\mathrm{Rack}},
$$

$$
I_{\mathrm{ABC}}=
(20\ \mu\mathrm{A})e_{\mathrm{volume}}
+(20\ \mu\mathrm{A})e_{\mathrm{accent}}.
$$

The physical transimpedance is 220 kohm. A Rack output calibration of
9.8181818 restores a practical modular level. The efficiency value is
calibrated against modern 662-family data and the modern-device clone
reference. The VCA output has its own half-band decimator, so both module
outputs traverse one resampling round trip.

### C38 output coupling

C38 is 1 uF and drives the 50 kohm volume load, giving

$$
f_{\mathrm{HP}}=
\frac{1}{2\pi(50\ \mathrm{k\Omega})(1\ \mu\mathrm{F})}
=3.183\ \mathrm{Hz}.
$$

The implementation subtracts an exact-exponential low-pass state:

$$
\alpha=1-\exp\!\left(-\frac{2\pi f_{\mathrm{HP}}}{f_s}\right),
$$

$$
\ell_n=\ell_{n-1}+\alpha(x_n-\ell_{n-1}),
\qquad
y_n=x_n-\ell_n.
$$

## Numerical experiments

### Filter reference

The independent Python model represents every coupling section as a
continuous one-pole state and integrates the complete nonlinear system with
SciPy DOP853 at tight tolerances. Regression tests compare the production
small-signal response from 30 Hz to 6 kHz, host rates from 44.1 to 192 kHz, DC
rejection, resonance thresholds, extreme drive, and 2x/4x agreement.

In a representative heavy-drive run, the first five odd harmonics agreed with
the continuous-time reference by approximately 0.02% to 2.4%. The automated
limit is 3.5% per harmonic.

Randomized abrupt parameter stress at 44.1, 48, 96, and 192 kHz produced finite
outputs and zero Newton failures in both oversampling modes. The largest
observed iteration count was eight.

### OTA transistor reference

The ngspice reference expands the Open Music Labs BA662 clone into individual
transistors and uses Nexperia PMP4201Y/PMP5201Y models. It compares ideal BJTs,
models with dynamic capacitances removed, full manufacturer models, and the
reduced OTA law. The deck includes the 12 V and 5.333 V supplies, 220 kohm load,
emitter-follower buffer, C38, and 50 kohm load.

Residual calculations retain the first nine steady-state harmonics before
level and delay matching. The reduced law is compared at the OTA node, while
the buffer and C38 are measured separately.

| Case | Reduced/full OTA | Full/static coupled |
|---|---:|---:|
| 1 kHz, 5 mV RMS, 20 uA | 0.00757% | 0.00529% |
| 1 kHz, 5 mV RMS, 5 uA | 0.00481% | 0.00450% |
| 1 kHz, 5 mV RMS, 40 uA | 0.01582% | 0.01010% |
| 1 kHz, 20 mV RMS, 20 uA | 0.03913% | 0.02387% |
| 1 kHz, 100 mV RMS, 20 uA | 0.20432% | 0.13442% |
| 10 kHz, 5 mV RMS, 20 uA | 0.05329% | 0.05142% |

At the nominal point, the full model produces 0.15434% THD and the reduced law
produces 0.15514%. The nominal buffer residual is 0.01616%; C38 contributes
0.00139%. A C38 sensitivity case with 1 kohm source resistance, 10 ohm ESR, and
1 Mohm leakage changes the gain-matched 20 Hz-to-1 kHz shape by less than
0.01 dB.

The reduced wrapper measured about 40.7 million samples/s in a representative
Windows release build. The older TriggerFish OTA-inspired and shipped
transistor VCA cores measured 2.16 and 1.42 million samples/s. Their
level-matched residuals are approximately 40% and 55%, reflecting their
additional nonlinear poles and saturation. A 9 kHz host-rate render of the
reduced VCA differed from an 8x reference by approximately 0.152% RMS; the
production voice therefore evaluates it at the selected 2x or 4x rate.

## Model boundaries

The production filter omits component tolerance, temperature drift, individual
diode mismatch, power-supply coupling, and transistor parasitics. The Bass
response is a calibrated low-order reduction. Fast and Slow accent modes are
behavioural implementations of published descriptions. The resonance makeup,
FM range, and Rack output gains are software calibrations.

The OTA reference uses a published clone topology with modern matched devices.
Production omits device mismatch, Early effect, and parasitic capacitances at
the measured operating levels.

## Reproducing the results

```console
uv run pytest
uv run python tests/python/benchmark_diode_ladder.py
uv run python tests/python/benchmark_ba662.py
uv run python tests/python/reference_ba662_spice.py --ngspice /path/to/ngspice
```

The SPICE harness stores downloaded manufacturer models and generated data in
the ignored `build/ba662-reference` directory.

## Source map

- `src/models/DiodeLadderFilter.hpp`: coupling network, ladder equations,
  Newton solve, bass response, and oversampling.
- `src/models/OtaVca.hpp`: reusable matched-pair OTA law.
- `src/models/Tb303Voice.hpp`: envelopes, accent, VCA wrapper, and C38.
- `src/TfDiodeLadderFilter.cpp`: Rack controls, ports, modulation, and
  polyphony.
- `tests/python/reference_diode_ladder.py`: analytic and continuous-time filter
  references.
- `tests/python/reference_ba662_spice.py`: transistor-level OTA reference.

## References

1. Roland Corporation, [TB-303 Service Notes](https://www.synfo.nl/servicemanuals/Roland/ROLAND_TB-303_SERVICE_NOTES.pdf).
2. Tim Stinchcombe, [Analysis of the Moog Transistor Ladder and Derivative Filters](https://www.timstinchcombe.co.uk/synth/Moog_ladder_tf.pdf).
3. Tim Stinchcombe, [A Comprehensive TB-303 Diode Ladder Filter Model](https://www.timstinchcombe.co.uk/index.php?pge=diode2).
4. Robin Whittle, [Devil Fish modifications and documentation](https://www.firstpr.com.au/rwi/dfish/) and [Devil Fish user manual](https://www.firstpr.com.au/rwi/dfish/Devil-Fish-Manual.pdf).
5. Robin Whittle, [TB-303 unique accent and envelope characteristics](https://www.firstpr.com.au/rwi/dfish/303-unique.html).
6. Open Music Labs, [BA662 clone](http://www.openmusiclabs.com/projects/ba662-clone/index.html) and [BA662 reverse engineering](http://wiki.openmusiclabs.com/wiki/BA662).
7. ROHM, [BA6110 voltage-controlled operational amplifier datasheet](https://experimentalistsanonymous.com/diy/Datasheets/BA6110.pdf).
8. Nexperia, [PMP4201Y matched NPN pair](https://www.nexperia.com/product/PMP4201Y) and [PMP5201Y matched PNP pair](https://www.nexperia.com/product/PMP5201Y).
9. Robin Schmidt, [Open303 TeeBeeFilter](https://github.com/RobinSchmidt/Open303/blob/313bf0d9ade7c1dcb6b3a74f5ea1780a29d70074/Source/DSPCode/rosic_TeeBeeFilter.h) and [Open303 envelope/VCA processing](https://github.com/RobinSchmidt/Open303/blob/313bf0d9ade7c1dcb6b3a74f5ea1780a29d70074/Source/DSPCode/rosic_Open303.h).
10. Paul Batchelor, [Soundpipe `tbvcf`](https://github.com/PaulBatchelor/Soundpipe/blob/3efb43bdabd0ed23b17c694292b5a79f1692a3ea/modules/tbvcf.c).
11. Vadim Zavalishin, [The Art of VA Filter Design](https://www.discodsp.net/VAFilterDesign.pdf).
12. VCV Rack, [Voltage Standards](https://vcvrack.com/manual/VoltageStandards) and [Plugin Development Guide](https://vcvrack.com/manual/PluginGuide).
