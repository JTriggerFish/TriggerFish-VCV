# TfReverb technical report

## 1. Purpose and rendering model

TfReverb accepts a mono or polyphonic input and produces a stereo room response.
A rectangular-room image-source model renders early reflections, and a 16-line
feedback delay network (FDN) renders the late field. Up to eight sources may
have independent X/Y positions in the two-dimensional floor plan. X coordinates
support automatic horizontal spreading; Y coordinates remain individually
stored. Source and listener heights use fixed values in the Rack module.

A room response has three useful regions. The **direct sound** is the first
arrival from source to listener. **Early reflections** are the following,
individually distinguishable arrivals from nearby walls, floor, and ceiling.
After enough bounces, the arrivals become too dense to associate with separate
surfaces; this is the **late field** or reverb tail. Source direction is most
useful in the direct and early regions. The late field is heard mainly through
its decay rate, spectral colour, echo density, and stereo spread.

TfReverb is a hybrid reverb. It calculates the sparse part from room geometry,
then uses a recursive delay network for the dense tail. Calculating every
physical reflection for several seconds would require a rapidly growing number
of image sources. The FDN produces that density with a fixed amount of work per
audio sample.

The renderer targets conventional two-speaker playback from a single listener
point. Listener coordinates specify that point in the room. Constant-power
amplitude panning represents horizontal direction; the model omits binaural
head filtering and interaural time differences.

In this report, **direct** names the positioned source branch, **wet** names the
early-plus-late room response, and **dry/wet mix** names the final crossfade
between those two stereo signals. A polyphonic input channel is treated as one
source with its own X/Y position.

Timing uses the input sample as the direct-arrival reference:

- the input sample is the logical arrival time of the direct sound;
- the direct path is rendered on that same sample and has no delay;
- geometric reflections use their excess travel time relative to the direct
  path;
- Pre-delay shifts the complete wet path by an artistic delay.

The immediate direct branch supports live playing while each reflection retains
its physically meaningful offset from direct arrival. Audio-interface buffering
still contributes to the computer's overall input/output latency; changing that
buffer does not change the sample offsets between TfReverb's direct sound and
its reflections.

## 2. End-to-end topology

The diagram follows one input sample from left to right. Each source first
splits into an immediate direct branch and a delayed room branch. The room
branch splits again: the image-source finite impulse response (FIR) produces
early reflections, while a source-specific delay and gain feed the shared FDN.
The early and late results meet at the ER / TAIL balance before wet filtering
and the final dry/wet mix.

```text
up to eight polyphonic audio channels + stored normalized XY positions
        (source Z uses a fixed normalized height)
        |
        +-- per-source distance gain and constant-power pan
        |       -> stereo direct -----------------------------------------+
        |
        +-- per-source wet Pre-delay                                      |
                |                                                         |
                +-- image-source stereo early-reflection FIR --+          |
                |                                             |          |
                +-- source-specific handoff delay and gain     |          |
                        -> mono source sum                      |          |
                        -> fixed unit-norm 16-line injection    |          |
                        -> 16 main delays                       |          |
                        -> three-band decay filters             |          |
                        -> optional four-bus octave return      |          |
                        -> two-stage velvet feedback matrix ----+          |
                                      ^                        |          |
                                      +-- slow delay modulation |          |
                                                               |          |
                       fixed orthonormal stereo decoder <-------+          |
                                                               |          |
        early reflections + late tail -> balance -> size calibration       |
            -> wet width -> wet low/high cuts -> constant-power dry/wet mix +
            -> output level -> stereo output
```

Direct sound and sparse early reflections carry source direction. Once the
response becomes a dense, statistically mixed late field, source azimuth is
discarded and all sources feed one position-independent room network. FDN line
indices represent statistical modes, while localization cues remain in the
direct and early fields.

The FDN receives a mono sum because its job is to synthesize the common room
tail. A fixed stereo decoder converts its internal modes into two decorrelated
loudspeaker signals. Section 8 derives this routing after introducing the FDN.

## 3. Notation and common control curve

Normalized controls are clamped to $[0,1]$. Several perceptual mappings use
the smoothstep function

$$
s(x)=x^2(3-2x).
$$

Its zero slope at both endpoints gives fine adjustment near minimum and maximum
settings. Bold lower-case symbols denote vectors and uppercase letters denote
matrices. The sample rate is $f_s$, the
speed of sound is $c=343\ \mathrm{m/s}$, and $V$ and $S$ are room volume
and total surface area.

The notation $z^{-m}$ means a delay of $m$ samples in the $z$ domain. A
diagonal matrix of such terms represents several independent delay lines. The
function $\operatorname{clamp}(x,a,b)$ limits $x$ to $[a,b]$. Amplitude ratios
use $20\log_{10}$ when expressed in decibels; power ratios use $10\log_{10}$.

## 4. Room geometry

### 4.1 Size and Aspect

Write the room dimensions as $\mathbf D=(D_x,D_y,D_z)$, with width, depth, and
height measured in metres. For a rectangular room,

$$
V=D_xD_yD_z,\qquad
S=2(D_xD_y+D_xD_z+D_yD_z).
$$

Volume $V$ controls the overall acoustic scale. Surface area $S$ also matters
because a long, narrow room meets boundaries more frequently than a compact room
of the same volume.

Across the Size control range, the minimum and maximum dimensions in metres are

$$
\mathbf D_{\min}=(2.8,3.5,2.4),\qquad
\mathbf D_{\max}=(18,25,8).
$$

Size $x$ log-interpolates each dimension:

$$
D_j=\exp\!\left(\log D_{\min,j}+
s(x)(\log D_{\max,j}-\log D_{\min,j})\right).
$$

At $x=0$ the dimensions equal $\mathbf D_{\min}$; at $x=1$ they equal
$\mathbf D_{\max}$. Log interpolation makes equal control movements correspond
more closely to equal perceived scale changes. A single Size control keeps the
geometric response and recursive time scale tied to the same room dimensions.

Aspect $a_c$ produces

$$
a=\exp((2a_c-1)\log 1.8),\qquad
D_x\leftarrow D_x\sqrt a,\quad D_y\leftarrow D_y/\sqrt a.
$$

The square-root factors multiply to one, so floor area, volume, and height stay
fixed as the room becomes wider and shallower or narrower and deeper. The Aspect
range is $1/1.8$ to $1.8$.

### 4.2 Coordinates

Source and listener coordinates are normalized fractions of the current room
dimensions. Listener fractions are clamped to $[0.02,0.98]$, and source
fractions to $[0.001,0.999]$, placing every point inside the room boundaries.
XY means lateral and front/back position. The module stores an XY pair for each
of the eight possible audio channels and displays the subset present at the
polyphonic AUDIO IN jack. Source and listener Z remain at normalized fractions
0.42 and 0.45 respectively. Their similar heights suit the floor-plan model;
the small offsets from the vertical centre avoid exact midplane symmetry.

For example, normalized source coordinates $(x_s,y_s,z_s)$ correspond to the
physical point $(x_sD_x,y_sD_y,z_sD_z)$ metres. The normalized coordinates stay
attached to the same relative part of the room when Size or Aspect changes.

The room plan uses top-view coordinates: increasing X moves from the left edge
to the right edge, and increasing Y moves from the top toward the bottom. A
source drawn left of the listener has negative lateral displacement and favours
the left speaker under the panning law below. The factory listener is drawn
below the sources.

An unedited single source starts at X=0.5. With multiple connected channels,
unedited source X positions are distributed evenly around the centre. The
half-span grows from 0.12 for two sources by 0.04 per added source, capped at
0.35, so a larger ensemble progressively occupies more of the room width.
Dragging a marker switches that channel to an explicitly stored X coordinate,
including when that coordinate is exactly 0.5. Resetting the marker switches it
back to automatic spreading. Patch data stores the automatic/manual state
independently of the coordinate value.

The factory Medium Hall has dimensions approximately
$9.35\times12.51\times5.24\ \mathrm m$, source position
$(0.50,0.35,0.42)$, and listener position
$(0.50,0.682,0.45)$.

## 5. Direct path

The direct branch performs two operations for each source: it chooses a stereo
pan from horizontal direction and applies a bounded distance gain. It introduces
zero samples of delay. This gives the dry sound the same apparent direction as
the first arrivals in the room response.

For a source-to-listener displacement in metres
$\mathbf d=(d_x,d_y,d_z)$, define

$$
q=\begin{cases}
d_x/\sqrt{d_x^2+d_y^2},&d_x^2+d_y^2>0,\\
0,&\text{otherwise},
\end{cases}
\qquad p=\frac{q+1}{2}.
$$

The ratio $q$ is the lateral component of the horizontal unit direction. It is
-1 for a source directly to the left, 0 for a source directly in front of or
behind the listener, and +1 for a source directly to the right. The remapping
$p=(q+1)/2$ puts that direction in the interval $[0,1]$ required by the panner.

The left and right speaker gains are

$$
g_L=\cos\left(\frac{\pi p}{2}\right),\qquad
g_R=\sin\left(\frac{\pi p}{2}\right).
$$

Since $g_L^2+g_R^2=1$, the panning law has constant power through the centre.
The endpoints produce $(g_L,g_R)=(1,0)$ and $(0,1)$, and the centre produces
$(1/\sqrt2,1/\sqrt2)$.
In two-speaker stereo, the front/back coordinate affects horizontal arrival
angle, distance, and the reflection pattern.

Distance gain is normalized to the room's characteristic length. Let

$$
\ell=\sqrt[3]{V},\qquad \rho=\frac{\|\mathbf d\|}{\ell},
$$

and let $\rho_0$ be the same quantity for the factory room and placement.
Then

$$
g_d=\operatorname{clamp}\left(
\frac{\rho_0}{\max(\rho,0.05)},\ 0.25,\ 2.0\right).
$$

The dimensionless distance $\rho$ measures source separation relative to the
room size. The reference value $\rho_0$ makes the factory placement unity gain.
The bounds are -12 dB and +6 dB.
Normalizing distance by $\sqrt[3]{V}$ decouples the Size control from dry level.
Every source is panned and scaled before summation, so polyphonic positions
remain independent.

The patched sample establishes direct arrival and is emitted in the same
processing call. Source placement changes direct gain and pan while retaining
the interchannel timing supplied by the patch.

## 6. Geometric early reflections

An impulse response is the output produced by a unit impulse. Convolving a
source signal with that response applies every stored reflection delay and gain
to the source. TfReverb builds one left and one right early-reflection impulse
response for each connected source.

### 6.1 Image-source construction

The image-source method turns a reflected path into an equivalent straight
path. For a single wall, mirror the source across the wall and draw a straight
line from the mirrored source to the listener. That line has the same length and
arrival angle as the physical one-bounce reflection. Repeating the mirror
operation tiles space with virtual rooms and accounts for higher reflection
orders. Allen and Berkley's image method gives the standard rectangular-room
construction [1].

In a rectangular room, repeated mirror images replace reflected paths with
straight paths from virtual sources. For image index $k$, source coordinate
$x_s$, and room dimension $D$, define

$$
p_k=|k|\bmod 2,\qquad n_k=\frac{k+p_k}{2},
$$

$$
x_{\mathrm{image}}=2n_kD+(1-2p_k)x_s.
$$

Applying this independently on all three axes enumerates image index triples
$(k_x,k_y,k_z)$. The triple $(0,0,0)$ belongs to the direct branch; the
early-reflection response uses the remaining triples.

In the one-dimensional formula, $k$ identifies a virtual room in the tiled
space. Its parity $p_k$ determines whether the source is mirrored, and $n_k$
determines which translated room contains it. Applying the same rule to X, Y,
and Z gives every virtual source position.

Let $d_0$ be direct distance and $d_i$ the distance from image $i$ to the
listener. Its stored delay is

$$
t_i=\max\left(0,\frac{d_i-d_0}{c}\right).
$$

The value $t_i$ is the excess path delay. Applying the wet Pre-delay $t_p$
makes the audible reflection time $t_p+t_i$ after the undelayed direct sample.

The image-to-listener vector supplies the horizontal arrival direction, which
uses the same constant-power stereo law as the direct sound. The
single-listener amplitude-panning model produces coherent stereo early
reflections.

### 6.2 Frequency-dependent path gain

The response uses four complementary bands separated at 250 Hz, 1 kHz, and
4 kHz. For band $b$, path $i$ has pressure gain

$$
G_{i,b}=\frac{1\ \mathrm m}{d_i}
10^{-d_i\alpha_b/20}
\prod_r R_{r,b}^{N_{i,r}}.
$$

Here $\alpha_b=(0,0.001,0.006,0.020)\ \mathrm{dB/m}$ models air absorption,
$R_{r,b}$ is a surface reflection-amplitude coefficient, and $N_{i,r}$
counts encounters with surface group $r$. The four surface groups are floor,
ceiling, side walls, and front/rear walls. The inverse-distance factor models
spherical pressure spreading.

“Complementary” means that the four band filters sum back to the original
full-band signal when their gains are equal. Each successive low-pass is
subtracted from the remaining signal, leaving four adjacent bands. A path can
therefore receive a different wall-loss product in each band.

Each factor in $G_{i,b}$ has a direct physical meaning. The first term reduces
pressure with distance. The second accumulates frequency-dependent air loss.
The product accumulates one reflection coefficient for every encounter with a
surface group. A long, high-order path is consequently quieter and darker than
a short, low-order path.

The Damping control sets the frequency dependence of wall absorption. It
smooth-log-interpolates each $R_{r,b}$ between bright and absorptive material
tables. Cascaded amplitude losses then vary smoothly in decibels over repeated
wall encounters. Minimum Damping selects the bright material table, which still
contains modest high-frequency absorption; higher settings make high-order
paths progressively darker.

Each fractional sample delay is deposited with four-point Lagrange
interpolation. The calculated arrival generally falls between two sample
indices; cubic Lagrange weights distribute its gain among four neighbouring FIR
coefficients. A causal form handles paths shorter than one sample and retains
both sub-sample geometry and causal FIR timing.

### 6.3 Choosing the early/late boundary

The response should hand over to a statistical late model when individual
echoes have become dense. A room-based first estimate is

$$
t_{\mathrm{avg}}=10^{-3}\left(20\frac VS+12\right),
$$

and a conservative estimate is

$$
t_{\mathrm{cons}}=10^{-3}(0.0117V+50.1).
$$

The average is limited to 20--150 ms. The conservative value is limited to the
same response horizon and kept at least 20 ms after the average. Image paths
are generated through the conservative value plus a 20 ms analysis margin,
with an absolute 150 ms ceiling.

These two empirical predictors provide a search interval. The $V/S$ term tracks
the average distance between boundary encounters, while the volume term gives a
later upper estimate for large rooms. The generated response then supplies the
actual handoff point for the current source and listener placement.

The estimates are refined from the generated response. A 20 ms window moves
in 2.5 ms hops and is accepted as diffuse when normalized echo density lies in
$[0.9,1.1]$ and absolute excess kurtosis is below 0.5. The condition must
remain true for 15 ms. Bounds around the predicted interval constrain the
handoff selected from these statistics.

Both statistics compare the window with noise-like late reverberation. For the
left and right samples in one window, let $\mu$, $\sigma^2$, and $\mu_4$ be the
mean, variance, and fourth central moment. The normalized echo density used here
is

$$
\eta=\frac{\Pr\{|x-\mu|>\sigma\}}{0.3173105},
$$

Here $\Pr$ is the fraction of pooled left and right samples satisfying the
condition. The denominator 0.3173105 is the corresponding probability for a
Gaussian signal, so a Gaussian-like dense window has $\eta\approx1$. Excess
kurtosis is

$$
\kappa=\frac{\mu_4}{\sigma^4}-3.
$$

Sparse impulses produce larger $|\kappa|$; a Gaussian signal has $\kappa=0$.
Requiring both tests to remain inside their limits for several windows stops a
single dense cluster from ending the early response prematurely.

For final handoff start $t_a$ and end $t_b$, both stereo FIRs receive the
same taper

$$
w(t)=\begin{cases}
1,&t\le t_a,\\
\cos\!\left(\frac\pi2 s\!\left(\frac{t-t_a}{t_b-t_a}\right)\right),
&t_a<t<t_b,\\
0,&t\ge t_b.
\end{cases}
$$

Using one envelope for both channels preserves arrival direction during the
handoff.

The cosine of smoothstep gives value and slope continuity at both ends of the
taper. Early energy is unchanged before $t_a$, fades over $[t_a,t_b]$, and has
reached zero by $t_b$. The FDN input described in Section 7 is aligned with the
same interval.

### 6.4 Zero-latency FIR rendering

For an FIR $h[k]$, direct convolution computes

$$
y[n]=\sum_{k=0}^{K-1}h[k]x[n-k].
$$

Its cost grows with the number of stored taps. Fast Fourier transform (FFT)
convolution is more efficient for a long response, but processing fixed-size
FFT blocks normally delays the result by one partition. TfReverb combines the
two methods so that early taps remain sample-synchronous and the long tail
remains affordable.

The first 128 coefficients of each source/channel FIR are evaluated by a
direct-form head. Coefficients from sample 128 onward use uniform 128-sample
FFT partitions, shifted so the partition latency restores their intended
absolute tap times. A reflection inside the first 128 samples appears at its
correct direct-relative sample. The convolution scheme has zero block and
partition latency.

For example, a coefficient whose geometric arrival is sample 160 belongs to an
FFT partition that itself takes 128 samples to produce an output. Storing that
coefficient 128 samples earlier within the partition makes its output appear at
sample 160. The direct-form head covers the interval in which that compensation
is unavailable. This split is an implementation detail: the rendered impulse
response retains the delays calculated from geometry.

Geometry and FIR construction run on a background worker capped at 20 scene
requests per second. Building a scene allocates memory and performs enough work
to disrupt the audio callback, so the worker prepares it away from the audio
thread. The worker publishes a complete, immutable bank containing the FIRs and
matching handoff data. The audio thread crossfades old and new banks over 100 ms,
which avoids clicks and keeps every item of scene metadata synchronized.

## 7. Early-to-late handoff

The early-reflection renderer and the FDN describe adjacent parts of one room
response. Their join needs a compatible start time and level. Starting the FDN
too early produces a dense cloud while the geometric echoes are still sparse;
starting it too late leaves a gap. A large level mismatch makes the join sound
like a second effect switching on or off.

For each source, the worker measures the untapered early-response power around
the handoff selected in Section 6.3. It measures the same four frequency bands
used for wall absorption and combines them with weights

$$
(0.16,0.29,0.34,0.21).
$$

The weights sum to one and give the two middle bands slightly greater influence,
where a level discontinuity is especially apparent. If $P_h$ is the weighted
power, the source's late-send amplitude is

$$
g_h=\operatorname{clamp}\left(
\sqrt{\frac{P_h}{2.5\times10^{-5}}},\ 0.35,\ 3.0\right).
$$

The denominator is a calibration power measured around the factory handoff. The
square root converts the power ratio into an amplitude ratio. The limits keep
unusual placements from creating a vanishing or excessive tail. This gain
changes the amount injected into the FDN. Separate filters inside the feedback
loop determine how long that energy remains in the tail.

The FDN begins with sixteen unequal main delays. Let $t_m$ denote their mean
duration and let $r_i$ denote the dimensionless ratio between line $i$ and that
mean. Section 8.2 derives their values. The earliest possible FDN output is
approximately

$$
t_{\mathrm{FDN},0}=t_m\min_i r_i,
$$

A source-specific alignment delay

$$
t_{h}=\max(0,t_a-t_{\mathrm{FDN},0})
$$

precedes the late input. The aligned, gain-matched sources are then summed to
mono. Thus a source whose geometric response becomes dense later also enters
the shared late field later. A scene bank stores its FIR, handoff gain, and
handoff delay under one sequence number so all three change together. Existing
energy remains in the FDN while a new scene fades in.

Explicit source direction ends at this summation. A mature room field has
undergone many reflections and is represented by diffuse stereo statistics.

## 8. Late field: 16-line FDN

### 8.1 FDN background

An FDN is a bank of delay lines connected in a loop [2]. On every sample, the
renderer reads the end of each line, attenuates and mixes the sixteen readings,
adds any new room input, and writes the results back to the line beginnings.
An impulse therefore circulates many times. Each circulation produces more
arrivals at the outputs, giving a long response with a fixed amount of work per
sample.

One delay line with feedback produces a regularly spaced comb of resonances.
Several unequal delays place those resonances at different frequencies. Mixing
the line outputs on every trip distributes energy among them, producing a much
denser set of resonant modes. Here a **mode** means one of the frequencies at
which energy can persist in the feedback system. A sparse or overly regular set
of modes can sound pitched or metallic.

The mixing transforms are orthogonal: for a vector of sixteen samples, their
output has the same sum of squared amplitudes as their input. They redistribute
energy among the lines while separate filters provide the intended loss.

For the time-invariant core, with fixed delay lengths and zero octave-up return,
the TfReverb loop can be written

$$
\mathbf d(z)=\mathbf Z_m(z)
\left[\mathbf F(z)\mathbf d(z)+\gamma\mathbf b\,u(z)\right],
\qquad
\mathbf y(z)=\mathbf C\mathbf d(z),
$$

where $u$ is the mono late input, $\mathbf d\in\mathbb R^{16}$ is the line
state, $\mathbf Z_m$ contains the 16 main delays, $\gamma=0.42$ is input
calibration, $\mathbf b$ is a unit-norm injection vector, and
$\mathbf C\in\mathbb R^{2\times16}$ is the stereo decoder.

The bracketed term is the signal written into the main delays: it combines the
returning line state $\mathbf F\mathbf d$ with the new input distributed along
$\mathbf b$. The diagonal matrix $\mathbf Z_m$ applies a different main delay
to every component. The decoder $\mathbf C$ forms left and right output samples
from the resulting sixteen line values.

The feedback operator is

$$
\mathbf F(z)=
\mathbf U_2\mathbf L_2(z)\mathbf D_2(z)
\mathbf U_1\mathbf L_1(z)\mathbf D_1(z)
\mathbf U_0\mathbf L_m(z).
$$

$\mathbf D_1,\mathbf D_2$ are velvet delay banks, the $\mathbf U_j$ are
orthogonal transforms, and each $\mathbf L$ applies the decay appropriate to
the immediately preceding physical delay segment.

In time order, a returning sample passes through the main-segment loss
$\mathbf L_m$, transform $\mathbf U_0$, the first velvet delays and their loss,
then $\mathbf U_1$, the second velvet delays and their loss, and finally
$\mathbf U_2$. The main delays $\mathbf Z_m$ complete the loop. Section 9
explains why the two extra delay banks increase diffusion.

### 8.2 Room-scaled main delays

For a rectangular room,

$$
t_m=\operatorname{clamp}\left(\frac{4V}{cS},\ 3\ \mathrm{ms},\ 78\ \mathrm{ms}\right).
$$

The expression $4V/(cS)$ is the mean free time between boundary encounters
for an isotropic ray field: a larger volume lengthens the average trip, while a
larger boundary area increases the frequency of encounters. The clamp keeps the
recursive time scale useful for very small and very large control settings.

Main delay $i$ is $r_i t_m$. The ratios $r_i$ are all different, so the lines
have different recurrence periods. Their arithmetic mean is one, making $t_m$
the actual mean delay. This normalization separates two design choices: room
geometry sets the common time scale, and the ratio set arranges the modes around
that scale.

### 8.3 Base and Optimized delay ratios

The context menu offers **Base FDN** and **Optimized FDN**. These selections use
the same sixteen-line architecture. They also share the velvet delays,
permutations, signs, orthogonal transforms, decay filters, moving-delay
processing, octave-up return, input vector, and stereo decoder. The selection
changes the sixteen main-delay ratios shown below.

| Line | Base | Optimized | Line | Base | Optimized |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.566877 | 0.566204 | 9 | 1.025879 | 1.026567 |
| 2 | 0.629627 | 0.629611 | 10 | 1.077396 | 1.078809 |
| 3 | 0.679982 | 0.679534 | 11 | 1.136272 | 1.135451 |
| 4 | 0.729562 | 0.730948 | 12 | 1.200184 | 1.201339 |
| 5 | 0.789987 | 0.788747 | 13 | 1.262934 | 1.262491 |
| 6 | 0.850413 | 0.851401 | 14 | 1.328007 | 1.327143 |
| 7 | 0.920135 | 0.920101 | 15 | 1.386496 | 1.386828 |
| 8 | 0.968165 | 0.967278 | 16 | 1.448084 | 1.447550 |

**Base** is the hand-curated reference set. Its delay distribution was derived
from the prime-number FDN family used in the velvet-noise FDN literature [3],
then normalized and tuned for this implementation. Distinct, irregular ratios
reduce repeated coincidences between line recurrences.

**Optimized** was produced by the repository's PyTorch model and optimization
script. PyTorch differentiates a loss computed from the static FDN frequency
response, and normalized steepest descent with adaptive backtracking updates the
main-delay ratios. “Static” here means that the analysis holds every delay
length fixed and excludes the octave-up feedback path. The response is
evaluated from 20 Hz to 20 kHz over 135 combinations of room geometry, decay
time, frequency-dependent loss, and short-delay mixing strength and span. The
model scores both representative input/output paths and the full internal
sixteen-line response, so an internal resonance still contributes when a
particular input/output combination hides it. The loss penalizes upward spectral
peaks above a local spectral average measured over 2, 8, 32, and 128 Hz
neighbourhoods. It combines average performance with a smooth approximation of
the worst case. The total objective adds 0.12 times a regularizer that separates
nearby delays and discourages repeated pairwise delay differences. Candidate
selection checks the training grid, a held-out grid, and a second independent
validation grid. The exported ratios use a 50% blend between the
gradient-descent result and Base.

Both sets have unit mean. Their largest relative ratio change is about 0.19%,
and the root-mean-square relative change is about 0.095%. The expected audible
change is therefore subtle: individual modal peaks and beating patterns move
slightly, with the Optimized set chosen for a smoother numerical resonance
measure across the tested controls. Listening preference can still depend on
source material and settings. Optimized is the factory selection.

Changing the selection crossfades between two fixed main-delay reads over
50 ms. This retains the circulating tail and avoids the pitch sweep produced by
continuously moving a delay tap.

### 8.4 Fixed input and stereo output vectors

Let the normalized Walsh vector for mask $m$ be

$$
W_m[n]=\frac{1}{4}(-1)^{\operatorname{popcount}(n\mathbin{\&}m)},
\qquad n=0,\ldots,15.
$$

Here $n$ and $m$ are four-bit integers, $n\mathbin{\&}m$ is their bitwise AND,
and $\operatorname{popcount}$ counts the resulting one-bits. This construction
produces sixteen Walsh vectors whose pairwise dot products are zero. A dot
product is the sum of component-by-component products; zero means that the two
directions are orthogonal in the sixteen-dimensional line space.

The injection vector is $\mathbf b=\mathbf W_0$. Its sixteen entries are all
$1/4$, so its squared norm is $16(1/4)^2=1$. One mono input sample is therefore
distributed equally across the lines with unit total squared amplitude. Define
$\mathbf m=\mathbf W_1$ and $\mathbf s=\mathbf W_2$; the decoder rows are

$$
\mathbf c_L=\frac{\mathbf m+\mathbf s}{\sqrt2},\qquad
\mathbf c_R=\frac{\mathbf m-\mathbf s}{\sqrt2}.
$$

The labels $\mathbf m$ and $\mathbf s$ suggest mid-like and side-like sign
patterns inside the FDN. Since they are orthogonal and unit norm, both decoder
rows also have unit norm and their dot product is zero. The left and right
outputs consequently have equal nominal energy and use different combinations
of line modes. A fixed factor of 1.30 calibrates their output amplitude. These
vectors give the tail a stable stereo field; source and listener direction has
already been represented by the direct and early branches.

## 9. Velvet feedback matrix and Diffusion

**Diffusion** controls how quickly a small number of arrivals becomes a dense,
noise-like response. In this design it changes the mixing strength and the span
of short delays inside the FDN feedback path. The following construction shows
how those two quantities are related.

An ordinary FDN mixes the line values at one instant and then sends them into
the main delays. TfReverb's **velvet feedback matrix** (VFM) adds two banks of
short, irregular delays between three mixing transforms. One sample returning
through the loop is consequently redistributed into many nearby sub-arrivals
before it reaches the main lines again. Repetition over many loop trips raises
echo density and softens regular comb patterns. This follows the
velvet-noise FDN construction described by Fagerström and colleagues [3].

Each transform $\mathbf U_j$ mixes every line with the other lines in four
layers. The layers pair indices separated by 1, 2, 4, and 8, so information can
reach all sixteen lines after four pairwise operations. Each pair uses the
butterfly matrix

$$
\begin{bmatrix}c&s\\s&-c\end{bmatrix},
\qquad c=\cos\theta,\quad s=\sin\theta.
$$

Multiplying this matrix by its transpose gives the identity for every
$\theta$. It therefore changes the distribution of a pair's energy while
preserving $x_1^2+x_2^2$. A fixed permutation and sign pattern follows each
four-layer transform. Permuting values and changing signs also preserves their
sum of squares, so every complete $\mathbf U_j$ is orthogonal.

The two diagonal delay banks place a different short delay on every line. Their
unscaled ranges are 0.125--9.75 ms and 0.375--10.125 ms. The name “velvet”
refers to the sparse, irregular timing pattern of these taps. With Diffusion
control $x$, let $u=s(x)$. The pairwise mixing angle and delay scale are

$$
\theta=\frac\pi{16}+u\left(\frac\pi4-\frac\pi{16}\right),
$$

$$
q_D=(0.20+1.30u)q_R,
\qquad
q_R=\operatorname{clamp}\left(\frac{t_m}{t_{m,\mathrm{factory}}},0.35,2.25\right).
$$

Here $q_R$ scales the short delays with the room's mean delay, limited to a
practical range, and $q_D$ adds the user-controlled amount. Every velvet delay
is multiplied by $q_D$. Diffusion therefore changes two related properties:
the amount transferred between paired lines and the time span over which that
transfer is scattered. Its minimum retains mild coupling and a short temporal
spread. At maximum, each butterfly gives equal-magnitude contributions from its
two inputs, and the velvet span is 1.5 times the room-scaled nominal span.

With fixed velvet-delay lengths, each $\mathbf U_j$ is orthogonal and each
$\mathbf D_j(z)$ is a bank of pure delays. The resulting velvet operator is
**paraunitary**: at every frequency, its matrix frequency response preserves the
total squared magnitude of a sixteen-channel signal [5]. Diffusion redistributes
that energy across lines and time. The filters in Section 10 determine its
loss.

## 10. Decay and Damping

Reverberation time, written RT60 or $T_{60}$, is the time required for the
reverberant amplitude to fall by 60 dB after excitation stops. A 60 dB
amplitude reduction is a factor of $10^{-3}$. TfReverb assigns a separate
target time to low, middle, and high frequencies, then derives each feedback
gain from the actual delay travelled [4].

The Damping control introduced for early-reflection materials in Section 6.2
also sets these late-field frequency-dependent decay times. This gives the
sparse reflections and dense tail a related spectral character.

The Decay control $d$ defines the middle-band target:

$$
T_{60,\mathrm{mid}}=0.25\,32^d\ \mathrm s,
$$

giving 0.25--8 seconds. For $h=s(\text{Damping})$, the other targets are

$$
T_{60,\mathrm{low}}=T_{60,\mathrm{mid}}\,0.72^h,
\qquad
T_{60,\mathrm{high}}=T_{60,\mathrm{mid}}\,0.20^h.
$$

For the late field, minimum Damping makes all three targets equal to the
middle-band value. As Damping increases, the low-frequency target reaches 72%
of it and the high-frequency target reaches 20%. The tail therefore loses
treble fastest while retaining a longer low and middle body.

Every physical delay segment—each main delay and both velvet banks—is followed
by a complementary three-band attenuation filter. Its approximate crossover
frequencies are 220 Hz and $\min(3500\ \mathrm{Hz},0.2f_s)$. “Complementary”
means that adding the three band outputs with equal gains reconstructs the
input. For a segment of duration $\tau$, the target amplitude gain in band $b$
is

$$
g_b=10^{-3\tau/T_{60,b}}.
$$

Suppose a signal travels through segments with durations
$\tau_1,\ldots,\tau_k$. Their gains multiply, and the exponents add:

$$
\prod_{j=1}^{k}10^{-3\tau_j/T_{60,b}}
=10^{-3(\tau_1+\cdots+\tau_k)/T_{60,b}}.
$$

After a total travel time $T_{60,b}$, the amplitude is therefore $10^{-3}$,
or -60 dB, however that time was divided among delay segments. This is why the
loss filters follow every segment: changing Size or Diffusion changes segment
lengths while preserving the requested decay time.

## 11. Modulation

Fixed delay lengths give the FDN fixed modal frequencies. Slowly varying those
lengths makes the modes drift by small amounts, reducing persistent ringing and
adding the gentle motion associated with a lush tail. Excessive depth or speed
is heard as pitch wobble or chorus, so the control spans subtle movement at its
lower settings and an intentional effect at its maximum.

For normalized Modulation control $m$, the depth factor is

$$
A=m^2.
$$

The square law gives the full control travel a useful purpose: small settings
have fine resolution, 50% gives 25% of maximum depth, and 100% reaches the full
depth. Any value above zero produces some motion.

The sixteen main lines use independent deterministic smooth-random modulators
with nominal rates from 0.061 to 0.229 Hz. A smooth-random modulator chooses a
sequence of seeded random targets and interpolates smoothly between them. It
lacks the obvious periodic cycle of a sine LFO while producing the same result
on every run. The maximum main-line excursion is
$\pm0.50\ \mathrm{ms}\,A$.

Both velvet banks have separately seeded smooth-random motion. At each sample,
the mean modulation value of a complete bank is subtracted from its sixteen
line values. Some taps move longer while others move shorter, preserving the
bank's average delay and reducing a common pitch shift of the whole late field.

Velvet depth is bounded per tap by both a fraction of its base delay and an
absolute limit:

- stage 1: at most 8% of the base tap and 0.35 ms;
- stage 2: at most 12% of the base tap and 0.65 ms.

A moving delay almost always falls between two stored samples. Four-point
Lagrange interpolation estimates the waveform at that fractional position from
four neighbouring samples, giving smoother motion than rounding to the nearest
sample. A causal margin limits very short taps so every required neighbour is
already in the history buffer. At zero Modulation the exact integer-delay read
is used. During modulation the operator varies with time; its delay centres and
nominal RT60 law remain fixed while the interpolation introduces small
instantaneous energy variation.

## 12. Shimmer

Shimmer creates a rising, pitched component by returning an octave-up copy of
part of the late field to the feedback loop. Because the copy recirculates, a
sound may be shifted again on a later trip, producing progressively higher
octaves. Filtering and a bounded return gain keep this process controlled.

After main-delay attenuation, the sixteen-line state is projected onto four
orthonormal Walsh coordinates. A projection takes the dot product with a unit
Walsh vector, reducing the sixteen line values to one bus signal. Four such
buses provide a balanced sample of the tank while limiting pitch-shifter cost.
Each bus is high-pass filtered at 250 Hz, shifted up one octave, and darkened by
two one-pole low-pass stages at
$\min(6.5\ \mathrm{kHz},0.2f_s)$.

The pitch shifter reads short overlapping pieces, called **grains**, from a
history buffer at twice the normal read rate. Each grain uses a 120 ms Hann
window, whose gain rises smoothly from zero and returns to zero. Eight grains
start at uniformly staggered phases, so their windows overlap to an almost
constant total gain and conceal individual grain boundaries. An eighth-order
low-pass before the history buffer limits aliasing; seeded pitch wander of at
most 3 cents and small read-position jitter reduce repetitive texture.

For Shimmer control $z$, the return gain is

$$
g_s=0.85s(z).
$$

The four shifted bus signals are projected back into sixteen lines with the
transpose of the same Walsh basis and added to the unshifted state. This mapping
preserves the bus energy and confines the shifted return to four known
directions in the tank. The resulting sixteen values pass through the complete
velvet matrix and main delays before reaching the stereo decoder.

The in-loop route makes successive octaves bloom and diffuse with the room tail.
At zero Shimmer the return gain is zero; the history and grain phases continue
to advance so automation can raise the control without an uninitialized start.

## 13. Wet output and final mix

The geometric FIR and FDN have already been timed and level-matched at their
handoff. The centre of the ER / TAIL control preserves that inferred
relationship. Turning the control toward either end attenuates the opposite
branch until the output contains a pure early response or a pure tail. For
normalized balance $b$, the branch multipliers are

$$
g_E=\min(1,2(1-b)),\qquad g_T=\min(1,2b).
$$

At $b=0$ the output contains the geometric early response. At $b=0.5$, both
multipliers are unity, preserving the inferred early response and the
geometry-matched late send described in Section 7. At $b=1$ the output contains
the FDN tail. Each branch gain is capped at its calibrated centre level.

Before width processing, both branches receive a room-size calibration

$$
g_V=\operatorname{clamp}\left(
\sqrt[3]{\frac{V}{V_0}},\ 0.30,\ 3.0\right),
$$

where $V_0$ is the factory-room volume. The cube root of volume has dimensions
of length and serves as the characteristic room scale. Image-source pressure
falls approximately with inverse distance, so multiplying by this relative
length prevents a larger room from becoming quiet merely because its paths are
longer.

For width control $w_c$, define $w=1.5s(w_c)$ and use the orthonormal
mid/side transform

$$
M=\frac{L+R}{\sqrt2},\quad S=\frac{L-R}{\sqrt2},\quad S\leftarrow wS.
$$

The transform separates a stereo signal into a common component $M$ and a
difference component $S$. Scaling $S$ changes stereo separation while retaining
the centre content. Displayed 100% width corresponds to $w=1$; the range is
0--150%. Width acts on the wet field. Direct-source placement has already set
the direct stereo image and proceeds to the final mix.

Wet filtering uses a one-pole low-pass followed by a one-pole high-pass formed
by subtracting a low-passed copy. The normalized controls map logarithmically:

$$
f_{\mathrm{low}}=20\,50^{x_L}\ \mathrm{Hz},
\qquad
f_{\mathrm{high}}=\min(1000\,20^{x_H},0.45f_s)\ \mathrm{Hz}.
$$

One-pole filters have a gradual 6 dB-per-octave slope and a small state, making
them suitable for broad tone shaping. Logarithmic mapping gives roughly equal
control travel to each frequency ratio. Both filters operate on the wet branch;
the direct branch retains the source spectrum.

Finally, for Mix control $x$, $q=s(x)$, and output level $G$,

$$
g_{\mathrm{dry}}=G\cos\frac{\pi q}{2},\qquad
g_{\mathrm{wet}}=G\sin\frac{\pi q}{2}.
$$

At the ends, one gain is one and the other is zero. At the centre, both are
$1/\sqrt2$, and their squared gains sum to one throughout the travel. This
constant-power law reduces the level dip that a linear crossfade creates for
uncorrelated signals. The positioned stereo direct signal is then mixed
channel-for-channel with the wet stereo signal.

## 14. Latency and time-varying controls

**Implementation latency** is a delay imposed by an algorithm before it can
produce output. It is distinct from the travel times and Pre-delay that form the
reverb itself. TfReverb adds zero implementation latency and therefore reports
zero samples to the host. Direct samples pass through current-sample gain and
pan in the same processing call.

Rack's engine block and the audio interface may buffer audio elsewhere in the
signal path. Those buffers affect the complete system equally and leave the
sample offsets inside TfReverb unchanged. The direct branch remains the
zero-time reference for its early and late arrivals.

Wet timing is the sum of acoustic and artistic delays:

- Pre-delay, from 0 to 250 ms, applies to both early and late paths;
- each ER tap occurs at its calculated excess path delay;
- the late input is aligned to the measured ER/late handoff where possible;
- the FDN then contributes its physical main and velvet delays.

At zero Pre-delay the wet delay reader returns the current live sample, giving
zero added Pre-delay. Increasing Pre-delay shifts both the early FIR input and
the FDN send by the same amount. The direct branch keeps its zero-time position,
so Pre-delay changes the onset gap between direct sound and the complete room
response while preserving the early/late join.

Front-panel parameters are smoothed over approximately 20 ms. Control state is
polled every 64 samples. Pre-delay, main-delay, and main-ratio changes crossfade
fixed read heads over 50 ms; scene, handoff, and velvet-transform/delay
transitions use approximately 100 ms. These transitions retain buffer contents
and pitch while moving between fixed read positions. Room-tail and Pre-delay
history continue through parameter changes.

The dedicated CV destinations are Pre-delay, Decay, Damping, ER / TAIL
balance, Mix, and Listener X/Y. Pre-delay automation uses the read-head
crossfade above. Listener CV changes direct placement immediately through the
smoothed control path, while its new early-reflection geometry is built on the
background worker and transitioned at a maximum of 20 scene updates per second.
The update rate suits gestures and slow LFOs on Listener X/Y. Size is a panel
control because every change rebuilds the complete room geometry and retimes
the late network.

## 15. Factory presets

Each context-menu preset sets all acoustic controls and preserves the user's
master Output Level.

| Control | Medium Hall (factory) | Small Room | Superlush |
|---|---:|---:|---:|
| Size / approximate dimensions | 60%; 9.35 x 12.51 x 5.24 m | 25%; 3.74 x 4.76 x 2.90 m | 78%; 14.72 x 19.02 x 6.89 m |
| Aspect | 50% | 50% | 55% |
| Pre-delay | 12 ms | 0 ms | 28 ms |
| Source XY | 50%, 35% | 38%, 35% | 50%, 32% |
| Listener XY | 50%, 68.2% | 55%, 68% | 50%, 68% |
| Decay | 2.30 s | 0.81 s | 4.59 s |
| Damping | 28% | 40% | 26% |
| Diffusion | 82% | 62% | 100% |
| Modulation | 30% | 15% | 100% |
| Shimmer | 0% | 0% | 45% |
| Wet low/high cut | 20 Hz / 12 kHz | 40 Hz / 9 kHz | 120 Hz / 8 kHz |
| ER / tail balance | Centre (50%) | 40%, early-leaning | 75%, tail-leaning |
| Stereo width | 100% | 86% | 127% |
| Mix | 35% | 22% | 40% |

Medium Hall is the neutral reference: a clear onset, dense tail, and modulation
set below chorus-like audibility. Small Room emphasizes geometric
ambience and a short damped tail, like an instrument heard in a pleasant room.
Superlush preserves the dry attack with Pre-delay and moderate Mix while a
wide, fully modulated, filtered shimmer tail develops behind it.

## 16. Properties checked by tests

The implementation and native tests enforce the following properties:

- the direct path has correct per-source pan, bounded calibrated distance gain,
  and exactly zero delay;
- ER tap times are direct-relative and independent of convolution partition or
  host block size;
- early stereo directions follow the analytical image paths;
- the late encoder and decoder use fixed, position-independent vectors;
- the static velvet operator conserves energy for every Diffusion setting;
- Decay sets RT60 independently of Size and Diffusion;
- Damping shortens high-frequency decay more strongly than low-frequency decay;
- Diffusion expands temporal scattering within one feedback topology;
- Modulation is active throughout its range, remains bounded, and uses
  fixed-head crossfades for parameter transitions;
- Shimmer re-enters the recursive late field before decoding and remains stable
  at maximum settings;
- the stereo decoder rows have equal norm and zero inner product;
- scene and control transitions preserve stored wet history.

## 17. Implementation map

The principal implementation files are:

| File | Responsibility |
|---|---|
| `src/tfdsp/room_reverb.hpp` | Complete signal routing, room mapping, direct path, handoff, wet filters, and width |
| `src/tfdsp/early_reflections.hpp` | Image paths, material response, handoff analysis, FIR construction, and convolution |
| `src/tfdsp/early_reflections_worker.*` | Rate-limited background scene preparation and delivery |
| `src/tfdsp/late_reverb.hpp` | Main FDN, input/output projections, decay, modulation, and shimmer routing |
| `src/tfdsp/velvet_feedback_matrix.hpp` | Paraunitary transforms, velvet delays, and their fractional modulation |
| `src/tfdsp/multiband_decay_filter.hpp` | Complementary per-segment RT60 filtering |
| `src/tfdsp/windowed_pitch_shifter.hpp` | Grain-based octave shifter used inside the feedback loop |
| `src/tfdsp/late_reverb_coefficients.hpp` | Base main-delay ratios, shared velvet coefficients, transforms, Walsh vectors, and calibration constants |
| `src/tfdsp/late_reverb_optimized_coefficients.hpp` | Optimized main-delay ratios used by the factory selection |
| `src/tfdsp/reverb_defaults.hpp` | Factory state and the three acoustic presets |
| `src/tfdsp/reverb_output.hpp` | Constant-power dry/wet output law |
| `src/TfReverb.cpp` | Rack parameters, smoothing, display, persistence, and context menus |
| `python/triggerfish_reverb/velvet.py` | Differentiable PyTorch model of the static late network |
| `python/triggerfish_reverb/objectives.py` | Multi-scale resonance loss used for coefficient optimization |
| `tools/optimize_velvet_reverb.py` | Gradient-descent search, validation, and coefficient-artifact export |

## 18. References

1. J. B. Allen and D. A. Berkley, [“Image method for efficiently simulating
   small-room acoustics”](https://jontalle.web.engr.illinois.edu/Public/AllenBerkley79.pdf),
   *Journal of the Acoustical Society of America*, 65(4), 1979.
2. J.-M. Jot and A. Chaigne, [“Digital Delay Networks for Designing Artificial
   Reverberators”](https://aes2.org/publications/elibrary-page/?id=5663),
   90th Audio Engineering Society Convention, 1991.
3. J. Fagerström, S. J. Schlecht, and V. Välimäki, [“Velvet-Noise Feedback
   Delay Network”](https://www.dafx.de/paper-archive/2020/proceedings/papers/DAFx2020_paper_23.pdf),
   23rd International Conference on Digital Audio Effects, 2020.
4. S. J. Schlecht and E. A. P. Habets, [“Accurate Reverberation Time Control in
   Feedback Delay Networks”](https://www.dafx17.eca.ed.ac.uk/papers/DAFx17_paper_11.pdf),
   20th International Conference on Digital Audio Effects, 2017.
5. S. J. Schlecht and E. A. P. Habets, [“Scattering in Feedback Delay
   Networks”](https://arxiv.org/abs/1912.08888), 2019.
