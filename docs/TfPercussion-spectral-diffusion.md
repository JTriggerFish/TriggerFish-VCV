# Experimental nonlinear spectral-energy diffusion

The September 2026 [beating-packet refinement](TfPercussion-beating-packets.md)
updates oscillator allocation/distributions and adds regional texture checks.
The diffusion control now separates concentration dependence from total-energy
sensitivity. Historical calibration results below are not
acceptance evidence for the changed packet allocator.

Status: the **only metallic transfer law exposed by the workbench** while we
test it. The oscillator bank remains. The previous transport primitives remain
available in C++, but the workbench recipe has no legacy switches, random
neighbour exchange or arrival-phase randomization parameters.

## Why test this?

The old cascade sends energy upward even in a quiet tail. Increasing its speed
can empty a gong's low body before its high spectrum develops. Its passivity is
useful, but does not establish that its spectral evolution is realistic.

Humbert, Josserand, Touzé and Cadot propose a **phenomenological** plate-turbulence
model for spectral energy density:

$$
\partial_t E_\omega = \partial_\omega
  (\omega E_\omega^2\partial_\omega E_\omega).
$$

Transfer follows a density gradient, with energy-dependent conductivity. A flat
density is an equilibrium. This is a statistical spectrum model, not a waveform
or phase generator, and the authors do not claim a formal derivation from the
full kinetic equation. [Paper and equations 3–4](https://arxiv.org/html/1709.09884).

Our experiment borrows that transport law. It does **not** make our random-phase
packets a physical nonlinear plate simulation.

## Signal and energy paths

```text
Contact excitation ───────────────► direct contact observation
        │
        ▼
Modal states ◄── repeated strikes add to existing states
        │
        ├── measure energy per painted packet
        │        │
        │        ▼
        │   spectral-energy diffusion
        │        │
        ◄────────┘ rescale states / seed previously silent packets
        │
        ├── shared frequency-dependent T60 damping
        └── painted prominence → contact/body mix → final EQ → output
```

This diagram separates responsibilities, not individual instruction order.
Turbulence still determines sideband distribution and phase motion. Observation
weights and EQ do not drive the transfer equation. Their output can change as
energy moves: conserving state energy does not imply constant audible loudness.

## Our discrete implementation

The coordinate is $x=f/(1000\,\mathrm{Hz})$. The total squared input-vector norm
defines one reference energy $E_\mathrm{ref}$; it is independent of observation
gain. With packet energy $e_i$ and frequency-cell width $w_i$:

$$
\rho_i = \frac{e_i}{E_\mathrm{ref}w_i},\quad
s=\frac{\sum_i e_i}{E_\mathrm{ref}},\quad p_i=\frac{\rho_i}{s},\qquad
g_{i+1/2}=\kappa\frac{x_{i+1/2}}{x_{i+1}-x_i}
s^b\left[\frac{p_i^2+p_ip_{i+1}+p_{i+1}^2}{3}\right]^a.
$$

Here $\kappa$ is **Diffusion strength**, $a\in[0,1]$ is **Concentration
dependence**, and $b\in[0,2]$ is **Energy sensitivity**. A silent bank bypasses
the solve; no division by zero or self-excitation occurs. Relative densities
are used only to calculate conductance, never to normalize the audio or states.

The old equation is recovered exactly algebraically at $b=2a$. In particular,
$a=1,b=2$ gives the paper-inspired quadratic closure. The independent exponents
are constructive controls, not a physical derivation. Existing workbench fits
import with the explicit value $b=2a$, visible in the UI and saved JSON. There
is no continuing link between the sliders. The historical JSON key
`bloom_energy_acceleration` now names concentration dependence; the new key is
`bloom_energy_sensitivity`.

At $a=0$, spatial conductivity is independent of relative concentration,
including at empty cells. At $b=0$, scaling all stored energies does not alter
conductance. For $b>0$, conductivity vanishes as total energy vanishes. Holding
spectral shape fixed, doubling oscillator amplitudes quadruples energy and
multiplies conductivity by $4^b$. Excitation brightness and contact behaviour
can still change with strike velocity when $b=0$. Transfer also follows the
density gradient: this factorization does not make all perceptual effects
independent, nor does it add a threshold, delay, gain or damping envelope.

Cell faces are arithmetic midpoints between **painted, tuned packet centres**;
the lower face is DC and the upper face extends half a neighbour spacing beyond
the highest centre, capped at Nyquist. Both have zero external flux. This avoids
changing the top cell's density just by raising the sample rate. Sideband positions do
not move these coordinates. Coincident centres share a cell, preserving their
existing energy proportions; newly excited silent cells use input weights.

The initial excitation uses the **same cell measure**: squared anchor input
weights are proportional to the excitation shelf squared times $w_i$, then
normalized once at preparation. Coincident handles split their cell weight.
Previously each handle received equal shelf-weighted energy irrespective of
its cell width: clustered handles invented density spikes and accelerated
nonlinear transfer. Double-precision normalization also prevents an extremely
dark shelf becoming an unintended volume control.

Each audio sample freezes conductivities at the previous energies, then solves:

$$
(w_i+\Delta t(g_L+g_R))\rho_i^{n+1}
-\Delta t g_L\rho_{i-1}^{n+1}-\Delta t g_R\rho_{i+1}^{n+1}
= e_i^n/E_\mathrm{ref}.
$$

The tridiagonal solve is linear in packet count, allocation-free, double
precision, and uses cancellation-free elimination. Nonnegative conductances
give nonnegative solutions; closed boundary fluxes conserve total energy up to
roundoff. There is no post-solve audio normalization or additional decay law.
No Newton solve is required. Stability at a large step is not proof of temporal
accuracy: the semi-implicit approximation still needs convergence checks.

Packet amplitudes are rescaled to their new energies. Arrival-phase diffusion
and random neighbour exchange are disconnected in this recipe. Local noisiness may be zero for a
tonal low packet and high for the upper receiving packets.

## Workbench controls and limitations

- **Diffusion strength:** `bloom_rate`; a coefficient in the stated
  coordinate/energy units, **not octaves per second**. Equal knob values do not
  imply equal audible transfer speed between laws.
- **Concentration dependence:** `bloom_energy_acceleration`; zero removes the
  local-concentration factor, one gives quadratic dependence on relative density.
  The cubic slider mapping gives fine adjustment near zero.
- **Energy sensitivity:** `bloom_energy_sensitivity`; zero removes the overall
  energy-speed dependency, one is proportional to total energy, two is quadratic.
  This affects soft/hard and repeated-strike response, not input velocity gain.
  Hold decay leaves this value fixed, even while compensating other controls.
- **Packet noisiness:** `field_turbulence`; the smooth extended packet response
  is now the only response used by this recipe. It controls sideband share,
  spread and stochastic bandwidth, not another energy-transfer process.

The removed JSON controls are `bloom_spectral_diffusion`,
`field_relaxed_turbulence`, `field_exchange`, and `bloom_phase_diffusion`.
These are removed from the compiled workbench descriptors, not merely hidden.
The fixed recipe definition explicitly selects spectral diffusion and extended
packet broadening, and excludes both extra exchange mechanisms.

All parameters remain serialized in JSON and exposed in the UI. No new hidden
gain, per-mode damping multiplier or onset delay is introduced.

### Bloom timing helper

The **Bloom timing…** button opens an in-page popover. Its earlier/later slider
shows before/current values and highlights the three ordinary sliders it moves.
The surrounding section groups controls into **Energy transfer** and
**Initial excitation**; the helper adds no DSP parameter.

For a captured centre and a gesture $x\in[-1,1]$, diffusion strength is multiplied
by $2^{-x}$, initial excitation tilt changes by $-6x$ dB/octave, and excitation
centre is multiplied by $2^{-x/4}$. Values respect the ordinary control bounds.
Positive $x$ therefore tends to start darker and redistribute more slowly.
Neither concentration dependence nor energy sensitivity is changed. This is a
relative starting-point gesture, not an exact onset-time guarantee or a delay.

**Set centre here** captures the current three values without changing the sound.
**Return to centre** (or double-clicking the slider) restores them. It does not
undo T60 compensation already accepted by **Hold decay**. With Hold decay enabled,
releasing the timing slider still runs the normal visible-parameter compensation;
progress and the result appear beside/below the helper button. Escape or clicking
outside closes the popover without undoing edits.

`tools/check_bloom_timing_meta.py` renders the current gong and crash at earlier,
centre and later, both at their saved energy sensitivity and at sensitivity one,
without Hold decay. All four checked combinations move the 3–12 kHz half-energy
time later. For the saved gong, high-band peaks are about 0.39/0.60/0.87 seconds.
At sensitivity one the gong's initial attack can dominate the peak statistic,
even though half-energy time still moves later; peak timing alone is insufficient.
These are diagnostic checks of these presets, not universal monotonicity claims.
The browser test checks exact parameter changes, live preview, reset, recentring,
keyboard/light dismissal, viewport placement and cleanup on preset switching.

### Split-control verification

`percussion_spectral_diffusion_tests.cpp` checks the old two-cell analytic law
at linked exponents, independent energy scaling, concentration influence,
zero-energy silence, positivity and conservation from tiny to very large states.
The existing frequency-grid and time-step convergence tests remain enabled.
The extra runtime work is a total-energy reduction and one global power;
the allocation-free, linear-time tridiagonal solve is unchanged.

`tools/check_diffusion_control_split.py` compares the archived old gong render
with its explicitly migrated settings and renders 0.25/0.5/1-strength strikes
at energy exponents 0/1/2. The checked six-second gong render was bit-identical;
that is a measured preset result, not a general bit-identity guarantee.
The diagnostic archive is `build/diffusion-control-split`. Snapshot import
tests verify that the new value is saved and never relinked after editing;
browser checks cover the visible controls/help and all four metallic presets.

Changing the painted frequency grid changes the discretization. Cell widths
prevent the trivial serial-stage slowdown, but sparse and refined grids are not
claimed to be identical. Large empty spectral gaps are coarsely represented.
Turbulence still affects frequency-dependent losses and audible colour; this
experiment cannot make every perceptual interaction disappear.

## Fitting and acceptance procedure

The prior diffusion-only presets failed audibly. In particular, gong 3–16 kHz
energy peaked at 30–50 ms instead of approximately 720–960 ms. Turning transfer
off barely changed those highs. Matching decay while immediately exciting the
upper spectrum was not a successful bloom fit.

The current developer workflow uses the **actual workbench Wasm**, not a
surrogate implementation:

1. **Ablate first.** `tools/audit_metal_strength.py` compares ordinary, no
   transfer, no stochastic phase width, and both disabled, across strengths.
   Keep output gain and reference-family gain fixed.
2. **Fit onset, rise and decay.** `tools/fit_metal_bloom_surface.py` uses
   `SpectralBloomLoss`: 24 logarithmic frequency bands, a 4096-sample Hann
   STFT with approximately 10 ms hops, and explicit regions from 0–50 ms
   through 5–6 s. Penalize absolute log-band energy and each early-to-bloom
   contrast separately. The fixed floor is 70 dB below reference peak power.
   This is a diagnostic objective, not a claim of perceptual equivalence.
3. **Search the actual controls.** Joint starts and finite-difference,
   bounded trust-region steps fit diffusion strength, concentration and energy sensitivity,
   excitation tilt/knee, packet width/noisiness and two T60 endpoints.
   Positive wide-range controls use logarithmic solver coordinates.
   Fixed parameters and finite-difference influences are recorded.
4. **Solve observation efficiently.** `SpectralBloomBasis` caches exact
   STFT cross-power matrices from an independently validated affine basis of
   actual renders. It includes interference terms. Analytic derivatives fit
   positive painted amplitudes without rerendering the DSP at every step.
5. **Check spectral quality independently.**
   `tools/polish_metal_bloom_perceptual.py` minimizes published auraloss
   multi-resolution Mel loss with a per-cell/rise guard. It uses analysis
   autograd, not differentiation through an invented instrument.
   `tools/refine_metal_ridges.py` allows small, visible frequency/packet-width
   adjustments when envelope agreement still leaves the wrong ridges.
6. **Check playing response.** `tools/audit_metal_velocity_grid.py` uses the
   actual crash edge layers at velocities 24/48/72/96/120. The optional
   `tools/fit_metal_velocity_colour.py` screens the existing velocity-colour
   control with nominal-tilt compensation, selects on 48/72/96 and reports
   24/120 as holdouts. There are no separate per-velocity gains or compressed
   strength curves.
7. **Review and publish.** Inspect fixed-scale spectrograms, differences and
   band envelopes; replay the complete snapshot; check held-out seeds and
   repeated hits. Publish the chosen JSON in the **main workbench**.
   Spectrogram pictures are local diagnostics, not separate audition pages.

The initial pass fixes recorded gesture, mode frequencies/count, model level
and body excitation gain. Later ridge refinement explicitly releases only
selected frequencies. No per-mode damping multipliers are fitted and both
presets retain two T60 endpoints. A smaller subproblem score does not authorize
a worse spectrum, uncontrolled playing response or a claim of listening approval.

## Evidence and remaining limitations

The updated gong is `build/metal-bloom-perceptual/gong`. On a common six-second
crop at unchanged gain, auraloss Mel improved from 1.60 on the actual previous
WAV to 1.11, alongside much closer delayed high-band growth and decay. Four
seeds measure 1.11–1.48. Eight quarter-note strikes peak at −6.9 dBFS raw;
eight rapid full-strength strikes peak at −1.8 dBFS raw. These finite tests
are not a universal headroom guarantee or listening approval.

Crash refinement also needs ridge placement, not just envelope matching:
the first new envelope-only candidate improved rise timing but worsened Mel
from 1.75 to 2.03. It was **not** accepted on that basis. Guarded spectral
polishing and subsequent ridge/width refinement address that regression.
The selected crash is `build/metal-bloom-velocity-fit/crash`: Mel 1.38 versus
1.75 for the actual previous WAV. The exposed velocity-colour setting is 1,
selected with three source velocities and two additional held-out velocities.
Eight quarter-note strikes peak at −4.3 dBFS raw. Eight rapid full-strength
strikes reach +7.7 dBFS **before browser master attenuation and limiting**;
this is not a unity-headroom guarantee. No limiter is inserted into the voice.

A final shared 400-Hz T60-knot trial (`build/metal-bloom-publish/crash`)
slightly improved the nominal Mel score but regressed two of three additional
seeds against the simpler candidate. It is not published. The two-endpoint
curve remains. The low tail, 3–6-kHz bloom and velocity-dependent peak times
still differ from the source; neither candidate is declared calibrated by ear.

## Harmonic design tools and initial gong pitch

The modal editor's **Generate modes** panel provides **Series**, **Base note**,
octave and an exact **Hz** field, with paired sliders/numeric inputs for count,
falloff, level and local noisiness response. The old **Quick shape** menu and
its ad-hoc layouts have been removed; use either explicit series as a starting point.
The metallic generator/parameter API accept 1 Hz through 15 kHz; the modal
display starts at 20 Hz. Lower modes are retained but not drawn at the edge.
The engine's former hidden 20-Hz modal clamp is removed too. The editor uses
a logarithmic axis so low-frequency anchors are not crowded into a few pixels.
The 40-Hz decay endpoint is retained deliberately: its damping value extends
flat below 40 Hz, without changing the interpretation of saved curves.
Observation high-pass filtering remains an explicit, separate control.

**Upper-mode stretch** is a generator operation, not runtime modulation:

$$
u_n=\max\left(0,\frac{n-h}{h}\right),\qquad
f_n = f_0 r_n\sqrt{1+(s u_n)^2},\qquad 0\le s\le1.
$$

Here $r_n=n$ for the harmonic series, or the tabulated membrane ratio; $n$ is
the one-based mode index. **Protected low modes** $h$ (1–8, default 4) keeps the first
$h$ modes unchanged. Above the core the multiplicative bend starts with zero
slope and increasingly spreads the upper modes. For a membrane template the
protected core retains membrane ratios, not integer harmonics. This is a
constructive design curve, not a calibrated physical gong formula.

At half stretch and a 100-Hz base, the first four modes remain 100, 200, 300
and 400 Hz; mode 8 becomes about 894 Hz and mode 16 about 2884 Hz. The law
does not depend on total mode count: adding modes never retunes existing
frequencies. Falloff remains dB per actual octave. Preview and generation use
the same formula and frequency ceiling. The harmonic guide remains unstretched.

This replaces the earlier power law $f_0r_n^{1+s}$, which detuned even the
second harmonic heavily. Existing JSON mode frequencies are not transformed.
The offline root/stretch fitter uses the shared JS inverse and rejects old
power-law grids rather than interpreting them as the new formula.

The generator is expanded by default for metallic instruments. Custom Hz is
not rounded to a musical note. **Replace modes** writes ordinary editable
frequencies/levels/local noisiness into the patch; there is no runtime harmonic
rule. It replaces the whole mode list, not the damping or diffusion controls.
The **Harmonic guide** shares the generator's base pitch; it only draws a grid
unless snapping is enabled. Generation no longer forces local noisiness to
zero: response 1 follows the global packet controls, response 0 explicitly
disables local broadening. This does not change any existing saved snapshot.

The generator previews count/frequency extent and enforces the formula,
capacity and frequency limits before applying. Both series support up to 32
handles; the membrane table contains the first 32 distinct circular-membrane
Bessel roots, checked against SciPy, with the original first 16 unchanged.
The kick editor still has a 16-handle capacity. Membrane roots are more closely
spaced than integer harmonics: the 32nd is about 6.746 times the base frequency,
not 32 times. Use upper-mode stretch to broaden that range if desired.
Invalid counts are rejected, not truncated.
Clipped observation levels and modes falling below the level floor are
reported in the preview. The modal graph retains a usable minimum height.

A persistent footer banner reports generator validation/application errors,
startup failures, worker failures and uncaught browser/promise errors. Ordinary
Rendering/Ready messages never dismiss it. It is deliberately not a popup.

The gong reference's first 300 ms contains strong ridges around 122, 247 and
375 Hz, but also a strong nonharmonic ridge around 344 Hz. The previous fit
underemphasizes the lowest ridge and the 344-Hz component. A pure harmonic
quantization of every handle would erase useful source structure.
`tools/refine_gong_strike.py` compares harmonic low-frequency starts and fits
existing low frequencies, levels and local noisiness with early (0–400 ms)
and full-duration Mel objectives, guarded by the existing bloom envelope.
Diffusion, excitation, global texture, gain and both T60 endpoints stay fixed.
This is an explicit small voicing subproblem, not a new hidden model parameter.

The preceding follow-up is `build/metal-gong-strike`: lowest centre 123.43 Hz,
with its observation bar raised 3.2 dB and the separate ~340-Hz bar raised
10 dB. Full Mel improves modestly on all four tested seeds (1.10–1.32 versus
1.11–1.34 for their previous renders). Fixed-scale plots retain delayed upper
bloom and matching broad-band decay. The low fundamental remains weaker
relative to the strongest ridge than in the source, and narrow ridge positions
are not all matched. This is a clearer starting point, not a completed fit.
Repeated-hit raw peaks are −7.6 dBFS for quarters and −1.5 dBFS for rapid hard
strikes in the eight-hit tests. User listening approval remains separate.

The current audition candidate uses the protected-core stretched series with
no individual frequency changes and only four broad level-shape coordinates.
See [Gong coarse fitting](TfPercussion-gong-coarse-fit.md) for the exact stages,
validation trade-offs and limitations. The [earlier detailed fitting run](TfPercussion-gong-stretched-fit.md)
is historical. Neither is an approved calibration.

Historical artifacts in `build/spectral-diffusion-v1`,
`build/diffusion-fit-v1`, `build/diffusion-fit-v2` and
`build/metal-bloom-cells` document earlier laws or fits. Their scores must not
be presented as current-renderer results. A true before/after comparison uses
the old saved waveform, cropped to the same duration—not old parameters
rendered through a changed model. The older switch-comparison scripts target
a retired parameter surface and are not the current fitting entry points.

This remains a reduced constructive model. Packet energy coordinates do not
capture every sideband's actual radiated frequency; broad stochastic phase
motion can create spectral skirts outside a packet. Noisiness therefore still
affects audible colour, delivered excitation and effective damping. Sparse
grids approximate large gaps coarsely. Matching one gong recording does not
establish accurate velocity behaviour for all gongs.

## Tests

- Spectral diffusion: conservation, positivity, silence, both directions,
  nonlinear scaling without a leakage floor, coincident cells and grid/time/
  sample-rate convergence.
- Cymbal preparation: excitation quadrature, duplicate-handle energy
  allocation and dark-shelf normalization; sideband changes do not move
  painted transport coordinates.
- Analysis: identity, fixed gain, rejection of premature highs, low-frequency
  leakage rejection, exact cached cross-power/Jacobian checks, and recovery of
  known rise time and T60 from enveloped noise.
- Runtime: native/Wasm parity, prepared/block replay and silent browser loading.

Build/test through `dev.ps1` (MinGW native and the optional Wasm target).
Python analysis and the browser remain optional development dependencies.
