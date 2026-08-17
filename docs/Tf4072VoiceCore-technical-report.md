# 4072 Voice Core technical report

<p align="center"><img src="../doc/Tf4072VoiceCore.png" height="520" alt="4072 Voice Core module"></p>

The [module guide](../README.md#4072-voice-core) describes the panel controls and
patching interface. This report describes the circuits represented by the
module and the equations evaluated by its DSP model.

## 1. Module architecture

`Tf4072VoiceCore` combines four functional blocks:

1. an ARP 4072 four-pole resonant low-pass filter;
2. an ARP 4019 voltage-controlled amplifier;
3. two envelope generators based on the ARP 4020 ADSR and Board-4 AR circuits;
4. Rack-facing control routing, polyphony, and sample-rate conversion.

The audio path is:

```text
FILTER IN -> 4072 filter -> LP OUT
                         -> normalled 4019 VCA -> VCA OUT

VCA IN, when patched, replaces the normalled filter-to-VCA connection.
```

`GATE` and `TRIG` drive both internal envelopes. The filter envelope is routed
to cutoff through `ENV→CUT`; the amplifier envelope is routed to the VCA through
`ENV→VCA`. `CUT MOD` and `VCA MOD` can add to or replace those internal
envelopes. Both envelope outputs remain available as 0–10 V signals.

Each envelope can operate as AR, AD, or ADSR. ADSR follows the 4020 timing
structure, AR follows the Board-4 gated envelope, and AD is a retriggerable
one-shot extension.

The model has the following dynamic state per polyphonic channel:

| Block | State | Runtime model |
|---|---|---|
| 4072 filter | Four LM3900 integrator voltages | Nonlinear implicit state equations |
| 4019 VCA | LM301 output-pole voltage | Nonlinear differential pair followed by one pole |
| Filter envelope | Stage, level, segment phase | Event-driven AR/AD/ADSR state machine |
| Amplifier envelope | Stage, level, segment phase | Event-driven AR/AD/ADSR state machine |

The filter and VCA run at 4x sample rate by default, with a 2x option. The
envelopes run at the Rack sample rate.

### 1.1 Nominal voltage mapping

The model uses the same numerical voltage on each side of the Rack-to-circuit
boundary. This one-to-one mapping follows the original 2600 calibration:

- the oscillator audio outputs are approximately 10 V peak-to-peak;
- the open 4072 is adjusted for unity gain from oscillator input to filter
  output;
- the fully open 4019 is adjusted to reproduce a 10 V peak-to-peak audio
  signal;
- the envelope generators produce 0–10 V control signals;
- filter pitch tracking uses 1 V/octave.

These levels coincide with Rack's
[voltage conventions](https://vcvrack.com/manual/VoltageStandards): ±5 V audio,
0–10 V unipolar CV, and 1 V/octave pitch. Drive and modulation attenuators
operate around that mapping. Supply headroom is handled separately because the
2600 uses ±15 V internal rails while Rack cable voltages are based on ±12 V
Eurorack rails.

## 2. Circuit analysis

### 2.1 ARP 4072 filter

The 4072 contains an input differential pair, four locally fed-back nonlinear
integrators, and an output level shifter. The Resonance potentiometer and its
return network are on the ARP 2600 main board.

<p align="center">
<a href="https://till.com/arptech/pdf/4072.pdf"><img src="Tf4072-original-schematic.png" width="100%" alt="Published ARP 4072 circuit schematic"></a>
</p>

<p align="center"><em>Figure 1. Published ARP 4072 schematic. Q9/Q10 form the
input differential pair. Q1–Q8 and the four LM3900 sections form the filter
stages. Z2B drives output pin 5. Click the image for the source PDF.</em></p>

<p align="center">
<a href="https://synthfool.com/docs/Arp/Arp2600/Arp2600ServiceManual.pdf"><img src="Tf4072-board3-wiring.png" width="100%" alt="ARP 2600 Board 3 wiring for the 4072 filter"></a>
</p>

<p align="center"><em>Figure 2. Board-3 wiring around the 4072. R160 is the
100 kΩ Resonance potentiometer. Its wiper returns to pin 2 through R161
(150 kΩ), with C62 (50 pF) in parallel. Click the image for the source service
manual.</em></p>

#### Input differential pair and resonance loop

Audio reaches the Q9 base through 100 kΩ. Resonance feedback reaches Q10 from
the Z2B output, the R160 wiper, and R161. Q9/Q10 share an approximately fixed
tail current. Their collector-current difference follows the matched-pair law

$$
i_\Delta=I_L\tanh\!\left(\frac{v_Q}{2V_T}\right),
$$

where $v_Q$ is the Q9-to-Q10 base-voltage difference, $I_L$ is the pair tail
current, and $V_T$ is the thermal voltage.

At 27 °C,

$$
V_T=25.85\ \mathrm{mV}.
$$

The 56 kΩ emitter resistor gives

$$
I_L=\frac{15-0.65}{56\ \mathrm{k\Omega}}
=256.25\ \mu\mathrm{A}.
$$

At low frequency, the input and feedback resistor networks apply the fractions

$$
a=\frac{220}{100000+220}=0.00219517,
$$

$$
b=\frac{2200}{150000+2200}=0.01445466.
$$

For external input $v_{in}$, output $v_o$, and Resonance-pot fraction $\rho$,

$$
v_Q=a v_{in}-\rho b v_o.
$$

A 5 V input therefore applies 10.98 mV to Q9. At maximum Resonance, a 5 V
output applies 72.27 mV to Q10. The ratio

$$
\frac{b}{a}=6.5848
$$

shows that the resonance return drives the pair much harder than the audio
input. The `tanh` current law then limits the error signal smoothly during high
resonance and self-oscillation.

C62 increases the feedback near the top of the audio band. Its small value and
the changing source resistance at the R160 wiper are omitted from the reduced
DSP network, which uses the low-frequency ratio $b$.

#### Conversion into the first filter stage

The Q9/Q10 collectors drive unequal base impedances. One side sees 220 Ω. The
other sees 220 Ω loaded by two 12.1 kΩ paths:

$$
R_b=\left(\frac{1}{220}+\frac{1}{12100}+\frac{1}{12100}\right)^{-1}
=212.2807\ \Omega.
$$

The mean differential load is

$$
R_L=\frac{220+R_b}{2}=216.1404\ \Omega.
$$

Expressed in the voltage domain used by the first stage, the nominal maximum
drive is

$$
L_{nom}=I_L(12.1\ \mathrm{k\Omega})\frac{R_L}{R_b}
=3.1570\ \mathrm{V}.
$$

Connector pin 10 of the 4072 is labelled `GAIN TRIM`; it is an internal service
connection. Board trimmer R163 is adjusted for unity gain with the filter open.
At low signal level the input pair contributes the slope $L/(2V_T)$, the audio
divider contributes $a$, and the output amplifier contributes $G$. Unity gain
therefore requires

$$
G\frac{La}{2V_T}=1,
$$

which gives the calibrated value

$$
L=\frac{2V_T}{Ga}=3.06172\ \mathrm{V}.
$$

#### Four current-integrator stages

Q1/Q2, Q3/Q4, Q5/Q6, and Q7/Q8 drive the four sections of the LM3900. An LM3900
is a Norton, or current-differencing, amplifier. Each section integrates the
pair's collector-current difference into its capacitor voltage and buffers that
state from the following stage.

Each driven transistor base receives the preceding-stage voltage and local
negative feedback through equal 12.1 kΩ resistors. The loaded 220 Ω base shunt
converts their signed voltage sum $v_s$ to the pair's base-voltage difference:

$$
v_{BE,\Delta}=\frac{R_b}{12.1\ \mathrm{k\Omega}}v_s.
$$

The normalized differential current is therefore

$$
\tanh(k_s v_s),
$$

with

$$
k_s=\frac{R_b}{(12.1\ \mathrm{k\Omega})(2V_T)}
=0.3393396\ \mathrm{V}^{-1}.
$$

The constant $k_s$ is the stage-pair voltage sensitivity. A larger value makes
the pair enter its curved region at a smaller voltage. Its inverse,
$1/k_s=2.9469$ V, is the signed stage-voltage sum that produces a `tanh`
argument of one. Local feedback keeps normal signals in the gentler part of
this curve; large stage voltages progressively limit the integrating current.

#### Output amplifier

The fourth stage feeds the inverting Z2B level shifter. R41 and R42 set its gain
magnitude:

$$
G=\frac{100}{13}=7.69231.
$$

Z2B drives both output pin 5 and the external resonance network.

#### Cutoff range

All four stages receive the same exponential cutoff-control current. The late
4072 has an unintended upper-frequency limitation. ARP documented changing the
four stage-frequency resistors from 4.7 kΩ to 2.2 kΩ to extend the range. The
module incorporates this correction and provides a 10 Hz–20 kHz cutoff control.

### 2.2 ARP 4019 VCA

The 4019 is a DC-coupled VCA built from a matched PNP audio pair, current
mirrors, linear and exponential control-current converters, and an LM301
transimpedance output stage.

<p align="center">
<a href="https://till.com/arptech/pdf/4019.pdf"><img src="Tf4019-schematic.png" width="100%" alt="ARP 4019 voltage-controlled amplifier schematic"></a>
</p>

<p align="center"><em>Figure 3. ARP 4019 schematic. Q9/Q10 are the audio
differential pair. The linear and exponential control circuits are at left;
the LM301 output stage is at lower right. Click the image for the source
PDF.</em></p>

<p align="center">
<a href="https://synthfool.com/docs/Arp/Arp2600/Arp2600ServiceManual.pdf"><img src="Tf4019-board4-wiring.png" width="100%" alt="ARP 2600 Board 4 wiring around the 4019 VCA"></a>
</p>

<p align="center"><em>Figure 4. Board-4 wiring around the 4019. The external
audio resistors, Initial Gain, linear AR control, and exponential ADSR control
are shown. Click the image for the source service manual.</em></p>

#### Audio path

Each audio input passes through a 100 kΩ Board-4 resistor, a 1 kΩ resistor
inside the 4019, and a 220 Ω base shunt. The input divider is

$$
S=\frac{220}{100000+1000+220}=0.0021734835.
$$

For external differential audio voltage $v_d$ and tail current $I_t$, Q9/Q10
produce

$$
i_o=I_t\tanh\!\left(\frac{Sv_d}{2V_T}\right).
$$

The current mirrors feed the LM301 transimpedance stage. Its low-frequency
target voltage is

$$
v_{target}=R_f i_o,
$$

where $R_f=56$ kΩ. This large-signal equation gives gradual compression as the
audio-pair base difference increases.

#### Gain controls

The ARP specification gives linear gain $V_L/10$ and exponential sensitivity
of 10 dB/V. With Initial Gain $g_0$, linear voltage $V_L$, and exponential
voltage $V_E$,

$$
g_L=\frac{V_L}{10},
$$

$$
g_E=10^{(V_E-10)/2},
$$

$$
g=g_0+g_L+g_E.
$$

The linear path reaches unity at 10 V. The exponential path reaches unity at
10 V and $10^{-5}$ gain at 0 V, covering the specified 100 dB range.

#### Output bandwidth

R26 and C5 form the LM301 feedback pole:

$$
f_o=\frac{1}{2\pi(56\ \mathrm{k\Omega})(100\ \mathrm{pF})}
=28.4205\ \mathrm{kHz}.
$$

Although its corner lies above 20 kHz, this pole affects the upper audio band:
its magnitude is -0.51 dB at 10 kHz, -1.07 dB at 15 kHz, and -1.75 dB at
20 kHz. It also attenuates ultrasonic harmonics generated by the audio
differential pair before decimation. The pole follows that nonlinearity, so it
shapes the generated spectrum without changing the signal driving the pair.

### 2.3 Envelope circuits

The ARP 2600 contains a 4020 ADSR generator and a simpler AR generator on
Board 4.

<p align="center">
<a href="https://synthfool.com/docs/Arp/Arp2600/Arp2600ServiceManual.pdf"><img src="Tf2600-board4-envelope-schematic.png" width="100%" alt="ARP 2600 board 4 AR, ADSR, and VCA schematic"></a>
</p>

<p align="center"><em>Figure 5. Corrected Board-4 schematic. The AR generator
is at left, the 4020 ADSR is at upper right, and the 4019 VCA is below them.
Click the image for the source service manual.</em></p>

#### 4020 ADSR

<p align="center">
<a href="https://till.com/arptech/pdf/4020.pdf"><img src="Tf4020-schematic.png" width="100%" alt="ARP 4020 ADSR envelope generator schematic"></a>
</p>

<p align="center"><em>Figure 6. ARP 4020 ADSR schematic. Gate enters pin 2,
Trigger enters pin 1, and Q3/Q4 with C1/C2 form the timing network. Click the
image for the source PDF.</em></p>

Gate establishes the held or releasing state. A Trigger pulse starts Attack.
During Attack, the timing capacitor charges toward +15 V; A1 switches to Decay
when the envelope reaches +10 V. Decay approaches the Sustain setting, Sustain
holds while Gate remains high, and Gate removal selects Release.

Normalizing the +10 V peak to one gives the zero-level Attack trajectory

$$
v_A(t)=1.5\left(1-e^{-t/\tau_A}\right).
$$

It reaches the peak at

$$
T_A=\tau_A\log 3.
$$

#### Board-4 AR

<p align="center">
<a href="https://synthfool.com/docs/Arp/Arp2600/Arp2600ServiceManual.pdf"><img src="Tf2600-board4-ar-schematic.png" width="65%" alt="ARP 2600 Board 4 AR envelope circuit"></a>
</p>

<p align="center"><em>Figure 7. Board-4 AR circuit. Attack and Release set the
charge and discharge paths around C69, Q5 responds to Gate, and A20 buffers the
output. Click the image for the source service manual.</em></p>

A high Gate charges C69 through the Attack resistance and then holds the peak.
A low Gate discharges C69 through the Release resistance. The resulting sequence
is Attack, Hold, Release.

## 3. DSP implementation

### 3.1 Sample-rate architecture

At each Rack sample, the module computes the two envelopes and the Rack-facing
control values for every active polyphonic channel. Matching seventh-order
polyphase IIR interpolators reconstruct the following signals at the internal
2x or 4x rate:

- filter audio input;
- logarithmic cutoff pitch;
- linear cutoff FM;
- resonance;
- VCA linear CV;
- VCA exponential CV;
- patched VCA audio input, when present.

Logarithmic cutoff is converted to hertz after interpolation. This preserves
the 1 V/octave trajectory between Rack samples. The normalled filter-to-VCA
path remains at the internal rate, avoiding an intermediate decimation and
interpolation.

One 2x stage uses the seventh-order Chebyshev-derived half-band design in
`tfdsp/sampleRate`. The 4x path cascades two of those stages. Matching
decimators suppress frequencies above the Rack Nyquist limit before returning
to the host rate. Filtering the modulation paths as well as audio prevents
host-rate control steps from being repeated unchanged across the internal
samples.

### 3.2 Nonlinear 4072 filter model

#### Variables and constants

For one internal sample:

| Symbol | Meaning | Value or range |
|---|---|---:|
| $u$ | Reconstructed filter input after Drive | Rack volts |
| $q$ | Reconstructed normalized resonance | 0–1 |
| $x_0,\ldots,x_3$ | Four LM3900 integrator states | volts |
| $y$ | Filter output before supply compliance | volts |
| $V_T$ | Thermal voltage | 25.85 mV |
| $a$ | Audio input divider | 0.00219517 |
| $b$ | Feedback divider | 0.01445466 |
| $L$ | Calibrated conversion from input-pair current to first-stage drive | 3.06172 V |
| $k_s$ | Local stage-pair voltage sensitivity | 0.3393396 V$^{-1}$ |
| $G$ | Output level-shifter gain | 7.69231 |
| $f_c$ | Cutoff after modulation and smooth bounds | hertz |

The four states are AC voltages; the large negative quiescent LM3900 voltages
have been removed analytically.

The signal flow represented by the equations is:

```text
audio u ----> input differential pair ----> equivalent drive L*ell
                  ^                                  |
                  |                                  v
output y -- resonance q ----------------------> stage 0
                                                   |
                                                   v
                                                stage 1
                                                   |
                                                   v
                                                stage 2
                                                   |
                                                   v
                                                stage 3 ----> gain G ----> y
```

Each stage also feeds its own state voltage back into its transistor pair.
The sign convention below places that feedback inside $z_i$ and gives the
state derivative a leading minus sign.

#### Continuous-time large-signal equations

The output is

$$
y=Gx_3.
$$

The input differential pair compares the audio and resonance-return voltages:

$$
v_Q=a\,u-q\,b\,y,
$$

$$
\ell=\tanh\!\left(\frac{v_Q}{2V_T}\right).
$$

The four local stage drives are

$$
z_0=L\ell+x_0,
$$

$$
z_1=x_0+x_1,
$$

$$
z_2=x_1+x_2,
$$

$$
z_3=x_2+x_3.
$$

With $\omega_c=2\pi f_c$, the four state equations are

$$
\dot{x}_i=-\frac{\omega_c}{k_s}\tanh(k_s z_i),
\qquad i=0,1,2,3.
$$

The LM3900 timing capacitor integrates the differential-pair collector
current. The common cutoff-control current scales that current equally in all
four stages. Writing the scale as $\omega_c/k_s$ makes the low-level slope

$$
\dot{x}_i\simeq-\omega_c z_i,
$$

because $\tanh(k_s z_i)\simeq k_s z_i$ near zero. Each stage therefore has the
requested small-signal pole frequency, while the complete hyperbolic-tangent
law limits the available charging current at larger signal levels.

These are the model's signal equations. The input-pair and
local-stage `tanh` functions are evaluated for every internal sample. Their
slopes change continuously with signal level, producing the low-level filter
response, driven-stage compression, resonance limiting, and self-oscillation
from the same set of equations.

#### Cutoff and resonance controls

Let $p_{knob}$ be the Cutoff-knob position in octaves relative to C4,
$V_{F\,1V/oct}$ the filter tracking input in volts, $E_f$ the 0–10 V filter
envelope, $a_e$ the bipolar `ENV→CUT` amount, and $c_f$ the internal-envelope
routing flag. The flag is one when the internal envelope is active. The base
pitch is

$$
p_i=p_{knob}+V_{F\,1V/oct}+0.5\,c_f\,a_e\,E_f.
$$

At 100% envelope amount, a 10 V envelope sweeps five octaves. For external
`CUT MOD` voltage $V_m$ and amount $a_m$,

$$
m_f=a_m\,V_m.
$$

An unpatched `CUT MOD` input leaves $c_f=1$. With a patched input, the `+`
route also gives $c_f=1$, while `EXT` gives $c_f=0$.

EXP mode gives

$$
f_{requested}=f_{C4}2^{p_i+m_f}.
$$

Here $f_{C4}=261.626$ Hz. The filter envelope always enters this exponential
pitch path, matching the normal ARP filter-envelope routing. The `CUT MOD`
switch selects only the law of the external modulation input.

LIN mode gives

$$
f_{requested}=f_{C4}2^{p_i}+(200\ \mathrm{Hz/V})m_f.
$$

The smooth lower knee is centred near 0 Hz. The upper limit is the lower of
20 kHz and $0.45$ times the Rack sample rate. For Resonance-knob value
$q_{knob}$, bipolar amount $a_q$, and external voltage $V_{RES\,MOD}$,
resonance is

$$
q=\operatorname{clamp}\!\left(
q_{knob}+a_q\frac{V_{RES\,MOD}}{10},0,1
\right).
$$

Drive converts its dB control to a linear input gain, with a maximum of +24 dB.

#### Implicit midpoint update

Resonance closes a feedback loop from $x_3$ to the input differential pair.
An explicit update would evaluate that loop from an older state, adding a
sample of numerical delay. The implicit midpoint method instead finds all four
next-state voltages together, with the nonlinear current laws evaluated halfway
between the accepted and candidate states.

The model is integrated at internal sample rate $f_{s,i}$. The prewarped
coefficient is

$$
\gamma=2\tan\!\left(\frac{\pi f_c}{f_{s,i}}\right).
$$

Implicit midpoint integration maps an analog one-pole frequency through a
bilinear transform. The tangent term compensates that mapping so the
linearized stage reaches the requested cutoff at the internal sample rate.

Let $x_n$ be the accepted four-state vector and $r$ a candidate for
$x_{n+1}$. Evaluate the nonlinear equations at the midpoint

$$
m=\frac{x_n+r}{2}.
$$

At that midpoint,

$$
y_m=Gm_3,
$$

$$
\ell_m=\tanh\!\left(\frac{au-qby_m}{2V_T}\right),
$$

$$
\phi_0=\tanh\!\left[k_s(L\ell_m+m_0)\right],
$$

$$
\phi_i=\tanh\!\left[k_s(m_{i-1}+m_i)\right],
\qquad i=1,2,3.
$$

The four residual components are

$$
R_i(r)=r_i-x_{n,i}+\frac{\gamma}{k_s}\phi_i.
$$

The accepted state satisfies

$$
R(r)=0.
$$

#### Newton solution

For any candidate $r$, the residual $R(r)$ measures how far that candidate is
from satisfying the midpoint update. Newton's method locally approximates the
residual by the Jacobian matrix and solves for a voltage correction that drives
the residual toward zero.

Newton iteration starts from the previous state:

$$
r^{(0)}=x_n.
$$

For iteration $k$, the correction $\Delta r_k$ solves

$$
J_k\Delta r_k=-R\!\left(r^{(k)}\right),
$$

where

$$
J_k=\left.\frac{\partial R}{\partial r}\right|_{r^{(k)}}.
$$

Define

$$
s_L=1-\ell_m^2,
$$

$$
s_i=1-\phi_i^2,
$$

$$
d=-\frac{Lq bG}{4V_T}s_L.
$$

The analytic Jacobian is

$$
J_k=
\begin{bmatrix}
1+\frac{\gamma}{2}s_0 & 0 & 0 & \gamma s_0d \\
\frac{\gamma}{2}s_1 & 1+\frac{\gamma}{2}s_1 & 0 & 0 \\
0 & \frac{\gamma}{2}s_2 & 1+\frac{\gamma}{2}s_2 & 0 \\
0 & 0 & \frac{\gamma}{2}s_3 & 1+\frac{\gamma}{2}s_3
\end{bmatrix}.
$$

Partial pivoting solves this $4\times4$ system. The largest correction component
sets the Newton step scale:

$$
D_k=\max_i|\Delta r_{k,i}|.
$$

For $D_k\le2$ V,

$$
\lambda_k=1.
$$

For $D_k>2$ V,

$$
\lambda_k=\frac{2\ \mathrm{V}}{D_k}.
$$

The estimate is updated explicitly:

$$
r^{(k+1)}=r^{(k)}+\lambda_k\Delta r_k.
$$

The 2 V step limit keeps a single local Jacobian from crossing most of the
$1/k_s=2.9469$ V stage-pair transition range. Ordinary iterations use
$\lambda_k=1$.

Convergence requires

$$
\max_i|R_i|<10^{-11}.
$$

The solver performs at most eight iterations. A converged estimate becomes the
new filter state. On failure, the previous accepted state is held for one
internal sample.

#### Filter pseudocode

```text
for each Rack sample and polyphonic channel:
    compute envelope and Rack-facing filter controls
    reconstruct audio, pitch, linear FM, and resonance at 2x or 4x

    for each internal sample:
        convert pitch to hertz
        add linear FM and apply the smooth cutoff bounds
        gamma = 2 * tan(pi * cutoff / internal_sample_rate)

        accepted = filter_state
        candidate = accepted
        converged = false

        repeat up to 8 times:
            midpoint = (accepted + candidate) / 2
            output_midpoint = G * midpoint[3]
            input_pair = tanh(
                (a * audio - q * b * output_midpoint) / (2 * VT)
            )

            z[0] = L * input_pair + midpoint[0]
            z[1] = midpoint[0] + midpoint[1]
            z[2] = midpoint[1] + midpoint[2]
            z[3] = midpoint[2] + midpoint[3]
            phi[i] = tanh(stage_sensitivity * z[i])
            residual[i] = (
                candidate[i] - accepted[i]
                + (gamma / stage_sensitivity) * phi[i]
            )

            if max_absolute_component(residual) < 1e-11:
                converged = true
                stop iterating

            build the analytic Jacobian above
            correction = solve_with_partial_pivoting(Jacobian, -residual)
            maximum_delta = max_absolute_component(correction)
            step_scale = 1
            if maximum_delta > 2 V:
                step_scale = 2 V / maximum_delta
            candidate = candidate + step_scale * correction

        if converged:
            filter_state = candidate
        else:
            filter_state = accepted

        circuit_filter_output = supply_compliance(G * filter_state[3])
        rack_filter_sample = rack_output_adapter(circuit_filter_output)

    decimate rack_filter_sample to the Rack sample rate
    apply the post-decimation Rack safety stage
```

#### Output compliance

The Z2B output is linear through ±10 V and approaches ±13.5 V smoothly. For
input $v$, define $h=3.5$ V. The compliance function is

$$
C(v)=v,\qquad |v|\le10\ \mathrm{V},
$$

$$
C(v)=\operatorname{sgn}(v)\left[
10\ \mathrm{V}+h\tanh\!\left(\frac{|v|-10\ \mathrm{V}}{h}\right)
\right],\qquad |v|>10\ \mathrm{V}.
$$

This is the internal 2600 circuit domain. The normalled filter-to-VCA
connection carries $C(v)$ directly, so Rack output protection cannot change
the level driving the VCA.

Rack-facing outputs use the shared soft adapter

$$
S_R(v;k,l)=v,\qquad |v|\le k,
$$

$$
S_R(v;k,l)=\operatorname{sgn}(v)\left[
k+(l-k)\frac{e}{\sqrt{1+e^2}}
\right],\qquad |v|>k,
$$

where

$$
e=\frac{|v|-k}{l-k}.
$$

Before decimation, $k=10.5$ V and $l=11.5$ V. The 10.5 V knee leaves Rack's
±10 V full-scale range unchanged and gives the overload transition 1 V of
headroom. After decimation, the same function uses $k=11.5$ V and $l=11.7$ V.
This final stage is normally an identity operation; it catches IIR overshoot
and keeps every finite Rack cable voltage below ±11.7 V without a hard clamp.

#### Small-signal calibration

This subsection analyzes the slope of the nonlinear equations around their
zero-input equilibrium. The audio path continues to evaluate the
large-signal equations above.

Using $\tanh(z)\simeq z$, one stage has transfer

$$
P(s)=-\frac{\omega_c}{s+\omega_c}.
$$

The complete linearized response is

$$
H(s)=\frac{GAP(s)^4}{1+qK_0P(s)^4},
$$

with

$$
A=\frac{La}{2V_T},
$$

$$
K_0=\frac{LbG}{2V_T}=6.5848.
$$

The unity-gain calibration makes $GA=1$. Four identical linearized poles reach
the self-oscillation condition near

$$
q=\frac{4}{K_0}=0.6075.
$$

Above this threshold, the full Q9/Q10 `tanh` law determines the oscillation
amplitude. The linearized response is used for gain, cutoff, and resonance
validation.

### 3.3 Nonlinear 4019 VCA model

#### Control routing

Let $E_a$ be the 0–10 V amplifier envelope, $a_a$ its amount, and $c_a$ the
internal-envelope routing flag. Let $V_m$ be the external `VCA MOD` voltage and
$a_m$ its amount:

$$
e_a=c_a\,a_a\,E_a,
$$

$$
m_a=a_m\,V_m.
$$

Let $q_a$ select the internal envelope law and $q_m$ the external modulation
law, with zero for LIN and one for EXP. The two 4019 control voltages are

$$
V_L=(1-q_a)e_a+(1-q_m)m_a,
$$

$$
V_E=q_a e_a+q_m m_a.
$$

An unpatched `VCA MOD` input leaves $c_a=1$. With a patched input, the `+`
position also gives $c_a=1$, while `EXT` gives $c_a=0$.

The amplifier-envelope law switch selects $q_a$. The `VCA MOD` law switch
selects $q_m$ and therefore affects only the external modulation input.

#### Runtime equations

The audio divider, feedback resistance, and output pole are

$$
S=0.0021734835,
$$

$$
R_f=56\ \mathrm{k\Omega},
$$

$$
f_o=28.4205\ \mathrm{kHz}.
$$

The unity tail current follows from the zero-input slope of the transistor
pair:

$$
I_{t,u}=\frac{2V_T}{R_fS}=424.7625\ \mu\mathrm{A}.
$$

Let $g_0$ be Initial Gain in the range 0–1. The combined Initial Gain and
linear-CV contribution uses a smooth positive knee with $\epsilon=10^{-6}$:

$$
g_{lin}=\epsilon\log\!\left[
1+\exp\!\left(\frac{g_0+V_L/10}{\epsilon}\right)
\right].
$$

The implementation evaluates this expression with stable asymptotic branches.
The exponential contribution is

$$
g_E=10^{(V_E-10)/2}.
$$

Their sum is limited smoothly to 16 times unity:

$$
g_{sum}=g_{lin}+g_E,
$$

$$
g=\frac{g_{sum}}{
\left[1+(g_{sum}/16)^4\right]^{1/4}}.
$$

The tail current and differential output current are

$$
I_t=I_{t,u}g,
$$

$$
i_o=I_t\tanh\!\left(\frac{S(v_+-v_-)}{2V_T}\right).
$$

The LM301 target is

$$
v_{target}=R_fi_o.
$$

For internal interval $\Delta t=1/f_{s,i}$, the exact one-pole coefficient is

$$
\beta=1-e^{-2\pi f_o\Delta t}.
$$

The single VCA state advances as

$$
v_o[n+1]=v_o[n]+\beta\left(v_{target}-v_o[n]\right).
$$

The same internal supply-compliance function $C(v)$ used by the filter is
applied to the VCA output. Rack output adaptation follows at the boundary
described above.

#### VCA pseudocode

```text
for each Rack sample and polyphonic channel:
    form the internal-envelope and external-modulation voltages
    route each voltage to the linear or exponential control input
    use the oversampled filter output, or reconstruct patched VCA audio
    reconstruct linear CV and exponential CV at 2x or 4x

    for each internal sample:
        linear_gain = smooth_positive(initial_gain + linear_cv / 10)
        exponential_gain = 10 ** ((exponential_cv - 10) / 2)
        gain_sum = linear_gain + exponential_gain
        gain = smooth_limit_to_16(gain_sum)
        tail_current = unity_tail_current * gain

        base_difference = S * differential_audio
        output_current = tail_current * tanh(base_difference / (2 * VT))
        target = 56k * output_current
        output_state += beta * (target - output_state)
        circuit_vca_output = supply_compliance(output_state)
        rack_vca_sample = rack_output_adapter(circuit_vca_output)

    decimate rack_vca_sample to the Rack sample rate
    apply the post-decimation Rack safety stage
```

The unity-current equation and the R26/C5 pole set low-level gain and bandwidth.
The complete `tanh` equation determines compression and harmonic distortion at
larger audio levels.

### 3.4 Envelope model

#### Modes and events

Each generator stores its mode, current stage, normalized level $y$, normalized
segment phase $u$, and Gate/Trigger Schmitt states. Gate and Trigger switch high
at 1 V and low at 0.1 V.

| Mode | Start or retrigger event | Stage sequence | Gate falling |
|---|---|---|---|
| ADSR | Trigger while Gate is high; Gate when Trigger is unpatched | Attack → Decay → Sustain | Starts Release |
| AR | Gate rising | Attack → Hold | Starts Release |
| AD | Gate or Trigger rising | Attack → Decay → Idle | Active cycle continues |

AD is a software extension derived from the AR circuit's ordinary RC Attack
shape. It uses the Decay slider for the automatic return to zero.

#### Curve control and segment shape

The shared **Curve** knob changes the shape of both envelopes while preserving
the Attack, Decay, and Release times selected by their sliders. Negative
settings make each segment approach a straight line. Positive settings make
the segment change more rapidly at first and approach its target more
gradually. The same rule applies in both directions: a curved rise has a fast
initial rise, and a curved fall has a fast initial fall.

This shape control is independent of the amplifier envelope's **LIN / EXP**
switch, which selects how the completed envelope signal controls the 4019 VCA.

At the centre position, the model uses two curve shapes derived from the ARP
circuits:

- The 4020 ADSR Attack charges its timing capacitor towards $+15\,\mathrm{V}$
  and ends when the capacitor reaches $+10\,\mathrm{V}$. This gives the curve
  coefficient $c_A=\ln 3$.
- Decay, Release, and the AR circuit follow an ordinary RC response. Their
  displayed duration represents a 95% RC transition, giving
  $c_R=\ln 20$.

AD uses the ordinary RC curve for both of its segments. The AR and AD Attack
curves consequently differ from the thresholded 4020 ADSR Attack even when the
Curve knob is centred.

A physical RC response approaches its target asymptotically. Each slider time
$T$ is defined as the exact end of the corresponding model segment. Normalizing
the RC response retains its shape while keeping the endpoints fixed as the
Curve knob moves.

For normalized segment time $u=t/T$, define

$$
R(u,c)=\frac{\operatorname{expm1}(-cu)}
{\operatorname{expm1}(-c)},
\qquad 0\le u\le1.
$$

At $c=0$,

$$
R(u,0)=u.
$$

For segment start $y_0$ and target $y_1$,

$$
y(u)=y_0+(y_1-y_0)R(u,c).
$$

Since $R(0,c)=0$ and $R(1,c)=1$, every curve starts at $y_0$ and reaches $y_1$
at exactly $t=T$. Curve adjustment therefore changes the motion within the
segment without changing its duration.

Let $q_c\in[-1,1]$ be the Curve-knob position, with zero at the centre, and let
$h$ be the applicable centre coefficient. The coefficient is interpolated
geometrically so that equal knob travel produces an even proportional change
in curvature. For $q_c<0$,

$$
c=h\left(\frac{h}{c_{min}}\right)^{q_c}.
$$

For $q_c\ge0$,

$$
c=h\left(\frac{c_{max}}{h}\right)^{q_c}.
$$

The resulting ranges are:

| Segment family | Minimum, $q_c=-1$ | Centre, $q_c=0$ | Maximum, $q_c=1$ |
|---|---:|---:|---:|
| 4020 ADSR Attack | $0.1$ | $\ln 3$ | $\ln 1000$ |
| Ordinary RC segments | $0.25$ | $\ln 20$ | $8$ |

The control spans positive coefficients. Values near zero approach linear
motion; larger values produce progressively stronger fast-then-slow curvature.

#### Discrete segment update

For envelope sample rate $f_s$ and selected duration $T$, advance phase by

$$
u_{new}=\min\!\left(1,u_{old}+\frac{1}{f_sT}\right).
$$

Define

$$
r_{old}=R(u_{old},c),
$$

$$
r_{new}=R(u_{new},c).
$$

The fraction of the remaining distance travelled during this sample is

$$
\eta=\frac{r_{new}-r_{old}}{1-r_{old}}.
$$

The state update is

$$
y\leftarrow y+\eta(y_{target}-y).
$$

When Attack is retriggered from a nonzero level, its starting phase is recovered
from the inverse curve

$$
R^{-1}(r,c)=-\frac{1}{c}
\log\!\left[1+r\operatorname{expm1}(-c)\right].
$$

This makes retriggering continuous.

#### Envelope pseudocode

```text
for each Rack sample and polyphonic channel:
    update Gate and Trigger Schmitt states
    detect Gate rising, Gate falling, and Trigger rising

    if mode is ADSR:
        Trigger rising while Gate is high starts or retriggers Attack
        Gate rising starts Attack when Trigger is unpatched
        Gate falling starts Release
    else if mode is AR:
        Gate rising starts Attack
        Gate falling starts Release
    else if mode is AD:
        Gate rising or Trigger rising starts or retriggers Attack

    select target, duration, and curve from the active stage:
        ADSR Attack -> 1, Attack time, 4020 Attack curve
        AR/AD Attack -> 1, Attack time, ordinary RC curve
        ADSR Decay -> Sustain, Decay time, ordinary RC curve
        AD Decay -> 0, Decay time, ordinary RC curve
        Release -> 0, Release time, ordinary RC curve
        Sustain/Hold -> keep the held level

    for a moving segment:
        advance normalized phase
        evaluate the old and new normalized curve values
        move level by the corresponding remaining-distance fraction

    at the segment endpoint:
        enter Decay, Sustain, Hold, or Idle according to the mode

    output 10 * level volts
```

#### Ranges and defaults

| Segment | Module range | Published ADSR range | Published AR range |
|---|---:|---:|---:|
| Attack | 1.4 ms–5 s | 1.4 ms–1.5 s | 20 ms–5 s |
| Decay | 6.4 ms–6 s | 6.4 ms–6 s | — |
| Release | 0.52 ms–6 s | 0.52 ms–6 s | 2.5 ms–5 s |

Both envelopes default to AD with 1.4 ms Attack, 1 s Decay, and 1 s Release.
The stored Sustain values are 50% for the filter envelope and 100% for the
amplifier envelope, ready for ADSR mode. The amplifier envelope defaults to
the VCA's exponential control path.

## 4. Combined oversampled signal path

When `VCA IN` is unpatched, every circuit-domain filter output sample is passed
directly to the 4019 VCA model at the internal rate. A separate copy passes
through the Rack output adapter for `LP OUT`. The VCA output enters its own
Rack adapter, and the two Rack-facing paths then use independent matched
decimators. This preserves the 2600 filter-to-VCA headroom while keeping both
cable outputs within Rack's voltage range.

When `VCA IN` is patched, the filter and VCA use separate input interpolators.
The filter still produces `LP OUT`, while the patched signal drives the VCA.

The default internal rate is 4x the Rack sample rate. A 2x context-menu option
reduces CPU use. Audio and modulation inputs entering the oversampled path are
interpolated to the internal rate; this avoids holding one host-rate value
across all internal samples during fast modulation.

## 5. Numerical validation

### 5.1 Filter

An independent continuous-time reference integrates the four nonlinear state
equations with SciPy's DOP853 solver. The comparison uses relative tolerance
$2\times10^{-11}$, absolute tolerance $2\times10^{-12}$, and a maximum solver
step of $1/(8f_s)$. At a 48 kHz Rack sample rate, the measured small-signal
magnitude errors were:

| Cutoff | Test frequency | 2x error | 4x error |
|---:|---:|---:|---:|
| 8 kHz | 8 kHz | +0.302 dB | +0.075 dB |
| 12 kHz | 12 kHz | +0.688 dB | +0.169 dB |
| 8 kHz | 18 kHz | -1.518 dB | -0.339 dB |
| 12 kHz | 18 kHz | -0.799 dB | -0.178 dB |

Large-signal tests compare harmonic amplitudes. At 4x, a 1 kHz, 50 V-peak
stress input with a 2 kHz cutoff and 0.55 resonance matched the continuous
reference within 2.9% through the ninth harmonic; the largest relative error
was on a seventh harmonic below 12 mV. An 80 V-peak case with an 8 kHz cutoff
matched within 4.6% through the ninth harmonic.

With an 800 Hz requested cutoff and full resonance, the continuous model
self-oscillated at 735.619 Hz with a 6.63671 V peak. The 4x implementation
produced 735.628 Hz and 6.63721 V peak.

These tests serve different purposes. The small-signal sweep measures the
discretization of the four-pole response. The harmonic tests exercise the
differential-pair nonlinearities. The self-oscillation test exercises the
closed feedback loop without an input signal.

### 5.2 VCA

The VCA differential-pair reduction was checked against ngspice using the
manufacturer's 2N3906 Gummel-Poon model as a modern matched-PNP proxy. Five
tail currents spanning 0.001 to 3 times the unity-gain current were tested.
After fitting the effective thermal voltage, the largest normalized-current
difference from the hyperbolic-tangent law was 0.019% at unity current. Using
the model's fixed value $V_T=25.85$ mV gave a maximum difference of 0.207%.

The 4x implementation's magnitude error relative to the analog 56 k$\Omega$ /
100 pF output pole was 0.0003 dB at 1 kHz, 0.0371 dB at 10 kHz, and 0.0889 dB
at 20 kHz. A 1 kHz large-signal comparison gave:

| Input peak | 4x implementation THD | 16x continuous-time proxy THD |
|---:|---:|---:|
| 0.1 V | 0.000147% | 0.000147% |
| 1 V | 0.014655% | 0.014650% |
| 5 V | 0.362525% | 0.362409% |
| 10 V | 1.404421% | 1.403971% |

The comparison supports the matched-pair `tanh` reduction and the modeled
output pole. A transistor-by-transistor solve would add substantial work per
sample while changing the nominal transfer by much less than component and
trim variation in an original module.

### 5.3 Envelopes

The ADSR Attack equation follows the nominal 4020 rail and threshold values.
It reaches 36.0% of its final value at one quarter of the selected Attack time
and 63.4% at one half.

Decay and Release use an endpoint-normalized $\log 20$ curve. Compared over one
displayed segment with an unbounded RC decay, the RMS difference is 0.0384 of
full scale. At the end of that interval the physical RC retains 5% of its
initial error, whereas the endpoint-normalized DSP segment reaches its target
exactly. This finite endpoint makes envelope stage timing deterministic while
retaining the characteristic RC curvature.

Retrigger tests cover every envelope phase and verify continuity at the event
sample. Mode tests cover Gate and Trigger combinations for AR, AD, and ADSR.

## 6. Model boundaries

The 4072 model includes the input differential pair, the four local
differential-pair integrator stages, the component-derived audio and feedback
ratios, the calibrated output gain, and smooth output-supply compliance. It
does not include transistor mismatch, finite current gain, Early effect,
LM3900 current-mirror error, amplifier slew and offset, temperature drift,
component tolerance, or parasitic capacitance. Capacitor C62 and the varying
R160 wiper impedance are reduced to the low-frequency feedback ratio $b$.

The 4019 model includes the audio differential pair, linear and exponential
control laws, summed control current, the R26/C5 output pole, and output-supply
compliance. Offset, control feedthrough, noise, current-mirror asymmetry, and
LM301 dynamics beyond the feedback pole are omitted. Those effects vary with
unit adjustment and may contribute VCA thump in original hardware.

The envelope model includes the nominal 4020 Attack threshold, ordinary RC
segments, Gate and Trigger behavior, and continuous retriggering. Component
tolerance, the short high-frequency feature measured on one modern 2600
clone, and the original keyboard's delayed Trigger pulse are outside the
model.

The following software provisions extend the component-level circuit models:
the 16x VCA control-gain limit, smooth $\pm13.5$ V internal output compliance,
the shared Rack output adapter, AD envelope mode, expanded cutoff range,
continuously variable envelope curvature, and the extended shared envelope
time ranges. The adapter's 10.5 V knee is a transition-quality calibration;
its 11.7 V limit follows Rack's protected-rail convention.

## 7. Reproducing the comparisons

Build the Python extension and run the regression suite:

```powershell
.\dev.ps1 python-test
```

Print the filter and VCA comparison tables:

```powershell
uv run python tests/python/benchmark_arp4072.py --full
uv run python tests/python/benchmark_arp4019.py
```

Run the device-level 2N3906 pair comparison with ngspice:

```powershell
uv run python tests/python/reference_arp4019_spice.py --ngspice C:/path/to/ngspice
```

Generated decks, logs, and sweep data are written beneath
`build/arp4019-reference`. The independent reference equations are in
`tests/python/reference_arp4072.py`, `tests/python/reference_arp4019.py`, and
`tests/python/reference_arp_envelope.py`.

The implementations are in
[`Arp4072Filter.hpp`](../src/models/Arp4072Filter.hpp),
[`Arp4019Vca.hpp`](../src/models/Arp4019Vca.hpp), and
[`ArpEnvelope.hpp`](../src/models/ArpEnvelope.hpp). Rack-facing rail handling
is centralized in
[`rail.hpp`](../src/tfdsp/rail.hpp).

## 8. References

1. [ARP Models 2600-2606 Service Manual](https://synthfool.com/docs/Arp/Arp2600/Arp2600ServiceManual.pdf) — module wiring, calibration, envelope operation, VCA control laws, and corrected schematics.
2. [ARP 4072 schematic](https://till.com/arptech/pdf/4072.pdf) — filter component values and ARP's documented high-frequency modification.
3. [ARP 4019 schematic](https://till.com/arptech/pdf/4019.pdf) — VCA input divider, matched pair, current mirrors, and output network.
4. [ARP 4020 schematic](https://till.com/arptech/pdf/4020.pdf) — ADSR switching, timing capacitors, and panel-control connections.
5. [Alan R. Pearlman and Dennis P. Gillette, *Dynamic Filter*](https://patents.google.com/patent/US4011466A/en) — filter-stage operation, local feedback, biasing, and differential-pair design.
6. [Rasmus Booberg, *Analysing and emulating the 2600 envelope and filter modules*](https://www.kth.se/social/files/650471ca9fc6eda417c3181b/dt2217-booberg-2023-analysing-and-emulating-the-2600-envelope-and-filter-modules.pdf) — modern-clone envelope measurements and filter analysis.
7. [Antonus 2600 technical features](https://www.antonus-synths.com/antonus2600/2600%20tech%20features.pdf) — measured envelope ranges.
8. [TI LM3900 datasheet](https://www.ti.com/lit/ds/symlink/lm3900.pdf) — Norton-amplifier operation and electrical limits.
9. [TI LM301A datasheet](https://www.ti.com/lit/ds/symlink/lm301a-n-mil.pdf) — output-stage electrical limits.
10. [onsemi 2N3906 datasheet and SPICE model](https://www.onsemi.com/pdf/datasheet/2n3906-d.pdf) — modern transistor proxy used in the VCA pair comparison.
11. [Vadim Zavalishin, *The Art of VA Filter Design*](https://www.discodsp.net/VAFilterDesign.pdf) — implicit virtual-analog integration and numerical background.
12. [VCV Rack Voltage Standards](https://vcvrack.com/manual/VoltageStandards) — Rack audio, CV, gate, and output-saturation conventions.
