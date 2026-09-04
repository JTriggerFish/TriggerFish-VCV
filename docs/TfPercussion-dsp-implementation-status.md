# Percussion DSP implementation status

This note records what the clean-slate component library implements and what
its tests currently prove. It is an engineering status report, not a claim
that the first crash model is calibrated.

## Implemented foundations

Reusable code lives in `src/tfdsp/percussion/`; no component is hidden in an
instrument or Rack module.

| Family | Components | Current analytic coverage |
| --- | --- | --- |
| Control | fixed-capacity linear/geometric breakpoint trajectory, asymmetric passive-loss controller | exact endpoints, continuous retrigger, sanitization, fast loss/slow release, ordered attenuation-only gains |
| Contact and compact body | half-sine force pulse, tonal chirp, enveloped noise, finite micro-contact burst, gated/finite renewal micro-contact process, 2x correlated FM burst with 4x reference, implement-derived direct/body projection | duration, endpoints, strength-to-energy mapping, deterministic seeds, requested stochastic event rate, coincident dense contacts, gate release, FM spectral convergence, bounds, routing |
| Delay and diffusion | 12-tap/2048-phase moving sinc delay, static Thiran delay, shared cubic Lagrange delay, fractional Schroeder all-pass | tone gain, low-frequency delay, polynomial exactness, integer-boundary continuity, impulse energy, five sample rates |
| Spectral motion | 255-tap antisymmetric FIR Hilbert transformer, phase-continuous signed SSB frequency shifter, and fourth-order translation-band guards | wanted level, image rejection, exact zero shift, through-zero automation, DC/Nyquist translated-content rejection, five sample rates |
| Resonance and loss | arbitrary damped modal bank, deterministic statistical modal cloud, modal passive-loss adapter, complementary three-band delay loss, orthogonal Givens mixer, projected wet-only fractional-comb network, explicit output submix | analytic modal frequency/T60 recurrence, independent projections, cloud determinism/range/normalization, passive attenuation, exact zero-coupling identity, scattering energy, delay T60 recurrence, excitation isolation, group routing |
| Reusable nonlinear dispersion | slow stochastic delay, 2x oversampled bounded self-phase delay with a 4x reference implementation, serial four-all-pass dispersion loop with explicit outer feedback; retained as an optional building block but disconnected from the active cymbal recipe | 1x/2x comparison against 4x, causal onset and declared recurrence delay, zero-drive linear reduction, no hidden feedback sample, long contractive stress |
| Active cymbal energy transport | passive adjacent-packet upward cascade operating on the modal field's stored complex states, with energy-dependent rate and magnitude-preserving destination phase diffusion | energy conservation, one-boundary-per-sample progression, upward centroid motion, stronger-strike acceleration isolated from velocity tilt, maximum-control boundedness |
| Radiation | guarded TDF2 biquad designs and a static high-pass/colour/low-pass chain | centre gain, pass/reject bands, sample-rate sweep, non-finite recovery |
| Observation | zero-capable static fractional delay and per-source gain, polarity, delay, and radiation paths | exact zero and integer delay, source isolation, polarity, non-finite recovery |

All sample loops are allocation-free. Tables and delay storage are created in
`Prepare`; filter and recurrence coefficients are built during preparation,
static configuration, hit creation, or quantized/smoothed frequency-shift
boundary updates. The active chirp evaluates one sine and
the self-phase delay evaluates its intentional nonlinear map at 2x per sample.
The reverb's velvet mixer advances its smoothed rotation algebraically, without
per-sample trigonometric calls. The
release micro-benchmark command, `dev.ps1 benchmark-percussion`, reports
component cost through the MinGW launcher without imposing machine-dependent
CI limits.

Rack calls `rack::system::resetFpuFlags()` on the main engine thread and every
engine worker. On x64 this enables FTZ and DAZ; on ARM64 it enables FTZ, with
DAZ supplied by the architecture. These flags are thread-local, so setting them
once during plugin initialization would not protect worker threads. Reusable
DSP still flushes non-finite and subnormal values at recursive filter states,
feedback delay writes, and public audio boundaries. This keeps the standalone
tests and non-Rack consumers safe without changing a caller's floating-point
environment from inside a component.

## Deliberate distinctions

`EnvelopedNoiseBurst`, `MicroContactBurst`, and `MicroContactProcess` are
separate. The first is a
compact oscillator-plus-noise source for tutorial-style cymbal/hi-hat graphs;
zero tilt is exactly white before its envelope. The second represents many
small implement contacts in one finite hit. The process adds finite clusters
and gated streams driven by a Poisson scheduler and smooth overlapping contact
windows. Density and incident amplitude remain separate controls. None is
mixed into an instrument output implicitly. The scheduler retains fractional
continuous event times and can return several contacts in one audio sample;
the optional contact exactly at gesture onset is a separate parameter.

`CorrelatedFmBurst` is a general compact body/contact source rather than a kick
voice. Amplitude, carrier frequency, and frequency-deviation trajectories are
independent. The modulator continuously blends seeded band-limited irregular
motion with a periodic oscillator, and optional perturbation is explicit.

The reusable dispersion return is stored in its base propagation delay. The signal is
read, passed serially through slow delay, four all-passes, self-phase delay and
loss, then written back with the new body drive. This avoids an accidental
extra sample at the outer feedback sum. The dispersion tap remains an analysis
and body-renderer drive output, not audible dry leakage. This component is not
part of the active metallic-plate recipe.

The coupled resonator is one superset: coupling zero gives independent combs;
nonzero coupling applies an orthogonal matrix before line loss. It is not the
room-reverb velvet FDN or its coefficient set.

## Performance profile

Measurements on 2026-09-04 used an AMD Ryzen 9 PRO 8945HS, GCC 16.2 MinGW
Release build, 48 kHz, and the current default crash graph. They are development
measurements rather than cross-machine budgets.

The reproducible probes are `dev.ps1 benchmark-percussion` for native
components and `dev.ps1 benchmark-workbench` for WebAssembly rendering and
STFT analysis.

| Path | Current |
| --- | ---: |
| Unified 408-mode crash with intrinsic cascade | 1,858 ns/sample |
| Compact FM kick | 275 ns/sample |
| Isolated unified field, cascade, phase and exchange | 1,946 ns/sample |
| Isolated dispersion loop | 143 ns/sample |
| Isolated 512-mode cloud | 529 ns/sample |

The unified crash uses about 8.9% of one core at 48 kHz. Its measured 128-frame
p99, including a hit in the sampled tail, is 284 microseconds inside a
2,667-microsecond deadline (10.6%). Its measured 16-frame p99 is 44 microseconds
inside a 333-microsecond deadline (13.3%). These figures qualify
one native monophonic core, not a complete Rack patch or a 96 kHz budget.

The modal bank sanitizes stable hit projections once, caches safe per-mode
damping multipliers until a damping control changes, and separates independent
state recurrence from ordered output reduction. This lets GCC vectorize the
hot recurrence and lets Emscripten use `simd128` without `-ffast-math` or a
change to finite-value recovery. A separate `-ffast-math` experiment reached
836 ns/sample, so the safe implementation is now within about 9% of that unsafe
upper bound.

The unified field additionally composes its fixed modal and random-phase
rotations during preparation, caches projected drive gains at hit time, and
expands complete PRNG words into branch-free sign arrays before its recurrence
and exchange passes. Static packet compatibility is removed from the audio
loop. Neighbour exchange operates on states already sanitized by propagation,
so it does not repeat finite/subnormal classification inside a bounded
orthogonal transform. A unified hit no longer prepares the inactive legacy
4,096-mode projections. Together these changes reduced the same 408-state
topology from 4,014 to 1,149 ns/sample; they change stochastic realization but
not the fitted modal, decay, projection, diffusion, or exchange parameters.
The intrinsic cascade first measured 3,844 ns/sample when it recomputed octave
gaps and transcendental transfer/rotation functions at every packet boundary.
Frequency gaps are now prepared once; packet energy is measured once per
sample; all adjacent transfers are solved before one combined state update per
packet; and bounded Padé/polynomial evaluations replace the hot `expm1`, `sin`,
and `cos` calls while retaining explicit energy normalization. This reduced the
same crash to 1,858 ns/sample without changing its energy-flow topology.

The remaining optimization order is:

1. Qualify the optimized stochastic statistics and sound against retained
   snapshots before considering mode-count changes.
2. Profile the cascade state-update pass and observation only if the assembled
   Rack module misses its production CPU budget.

For the browser workbench, the optimized Wasm crash graph renders ten seconds
in about 646 ms (6.5% of real time), while the compact kick takes about 57 ms
(0.6% of real time). A ten-second 2048-point, 75%-overlap STFT
takes about 73 ms and stores 3.67 MiB. A 1,089 x 506 heatmap redraw takes about
22 ms. Rendering and analysis remain asynchronous; Canvas replacement remains
lower priority. Reference spectra use an eight-entry LRU cache (about 29 MiB at
the default analysis settings). Snapshots retain audio and controls but
recompute their spectrogram on restore, reducing a ten-second snapshot from
about 5.5 MiB to about 1.8 MiB.

`RadiationFilter` is an object/observation voicing stage. Its final biquad is a
two-pole low-pass, so the asymptotic electrical slope is 12 dB/octave. It is
not, by itself, a physical distance-dependent air-absorption model: atmospheric
loss belongs in an optional propagation/observation stage and must vary with
distance, humidity, and frequency when that distinction matters.

## Current experimental graph and replacement target

`CrashCymbal` composes the tested contact, intrinsic modal energy cascade,
passive-loss, unified modal-field, and observation primitives. The exact C++ graph is
exposed to Python and to the optional browser workbench through WebAssembly.
Graph tests cover repeatability, strength, location, hardness, passive mute,
finiteness, and five sample rates.

The complete implemented signal flow and the unified field's equations are in
the self-contained
[nonlinear resonator architecture](TfPercussion-nonlinear-resonator-architecture.md).

`CompactKick` is the second assembled recipe. It uses two independently
parameterized, overlap-safe correlated-FM bursts plus a short tilted-noise
click, a fixed three-source mixer, and the shared observation/radiation stage.
Eight event voices preserve decaying hits during retriggering. Strength scales
level, initial pitch, and FM deviation; hardness increases secondary-branch FM
and click energy. Its three optional source routes are prepared switches at fixed
call sites. The defaults are a constructive starting voice, not a fitted claim
about a particular acoustic kick.

`SnareDrum` is the fourth compiled recipe, after the general membrane recipe.
It reuses the membrane's direct and body taps and drives `WireRack` only from
high-passed membrane motion. The rack adds a bounded contact envelope,
correlated tilted noise, and 8--48 normalized short wire modes. Static
membrane displacement therefore cannot sustain wire sound, while a rapid
strike can. Direct contact, membrane, and wires remain separate until a
three-source observation stage. The workbench exposes wire engagement,
release, threshold, bandwidth, density, colour, and decay, plus a prominent
body-ring frequency, decay, and level; these are visible fitted starting
controls rather than hidden constants.
Its rejected first numerical start and the lesson from that failure are
recorded in [the compact-snare fit note](TfPercussion-snare-initial-fit.md).

An earlier coupled-comb/frequency-shift graph was rejected during calibration:
its controls could not place persistent ridges independently. The implemented
replacement originally sent direct body drive to a small arbitrary modal bank
and dispersed drive to a deterministic statistical cloud, with a separate
turbulent residual. Those primitive implementations remain available for later
instruments, but the old graph, routing, state, and A/B selector are
disconnected from `CrashCymbal` and the workbench.

The active body is one 408-mode stochastic field derived from
twenty-four editable anchors in a 40 Hz--15 kHz constructive design range. Each anchor expands to a coherent centre mode and
sixteen nearby satellites. Turbulence transfers normalized excitation energy
from the centre to the satellites, widens their ERB-scaled frequency packet,
and enables energy-preserving phase diffusion. At high settings it also turns
on alternating local Givens rotations, allowing passive energy exchange within
and between neighbouring packets. A per-anchor response scaler can keep
selected ridges clean while the global control diffuses the rest. Pole radii
and the external mute controller
remain the only declared loss mechanisms. Bloom moves actual stored energy
upward through adjacent packets; it is not a delayed signal or output layer.
Its rate is explicitly energy-dependent, so stronger strikes progress farther
upward even when the separate velocity-brightness projection is disabled. The
field uses one body observation/radiation path.

## Calibration boundary

These tests are the first analytic quality level only. The following gates are
still open:

- add smoothed, state-preserving live retuning for size and tune;
- extend the existing five-rate, non-finite, denormal, retrigger and automation
  tests to rapid simultaneous control changes across the assembled graph;
- capture and isolate private corpus A, then audition its shared body and
  location/hardness projections against fit repeats;
- add remaining component-ablation views to the implemented residual-isolation,
  modal-resynthesis, branch-solo, and control-sweep report; and
- validate against held-out repeats, independent licensed crashes, and open
  corpora.

The previous aggregate-loss checkpoints were rejected by direct listening.
Automated optimization is now experimental diagnostic tooling; the active path
is a versioned, browser-based ear-fitting workbench using the native C++ graph.
Numerical views explain differences and protect regressions, but do not approve
a fit.

The modular-patch layer serializes that graph as named contact, modal-body,
observation, and mono-output modules with validated
ports and an acyclic routing schedule. Connections are Boolean switches only;
all adjustable levels are named module parameters present in the UI, JSON, and
C++ descriptor surface. The workbench displays the same routing
and links module-role colours to the complete parameter panel. Three optional
metallic-recipe connections are executable as prepared switches at fixed C++ call
sites and can be toggled in the expanded view. The adapter rejects silent or
unsupported structures. A second compact-kick patch exercises the same schema,
shared C/WebAssembly API, snapshots, routing view, contextual controls, live
AudioWorklet path, and fixed-schedule execution with a structurally different
six-node graph. A third membrane patch adds fixed two-source mixers, an
eight-voice contact/FM transient pool, a shared 16-mode body, accumulating
strike-history-driven tension, and selectable bypass/radiation/four-band
presentation. The history envelope is intentionally distinct from measured
modal energy. The modal drive and readout are mode-count normalized, while a
passive stored-state capacity bounds rapid retrigger growth without an output
clipper. Prepared blobs contain the compiled modal coefficients, and the short
observation delays use fixed storage, so worklet installation does not allocate.
Its tom and acoustic-kick starting patches share that compiled topology. The
native membrane graph measured about 216 ns/sample at 48 kHz on the development
machine (about 1.0% of one core at nominal rate). Module replacement and
arbitrary new connection endpoints wait
for another registered, statically ordered recipe.

The fourth patch is the statically scheduled membrane-plus-wires snare. Seven
prepared route switches cover its four exciter paths, body observation,
body-to-wire drive, and wire observation. The four exciter-to-bus levels are
ordinary visible mixer parameters, not connection coefficients. Its browser panels and live
AudioWorklet use the same C++ recipe and prepared-state handoff. The local
reference browser also discovers small allow-listed grids for one acoustic
snare, one kick, a gong, a 21-inch ride, and a 14-inch hi-hat when those development samples are
installed; no sample audio or machine path is copied into the repository. The
toolbar exposes one obvious middle-velocity calibration entry for the crash,
snare, kick, gong, ride, and hi-hat, rather than requiring the user to reconstruct a
reference/recipe pairing from separate selectors. The exact pairings and load
semantics are recorded in
[reference calibration starting points](TfPercussion-reference-calibrations.md).

The Python numerical-analysis implementation and its remaining real-data gates
are recorded in `TfPercussion-analysis-toolkit.md`. The static Plotly report is
a regression fallback. Live WebAssembly audition is implemented as a separate
development-only target and remains absent from normal Rack and release builds.

## Deferred structured-model extensions

`EnergyCoupler`, `AcousticCavity`, and `DistributedContactCoupler` remain
documented future components rather than prerequisites for the compact
video-derived instruments. Existing orthogonal mixing covers passive exchange
inside one resonator bank. Cross-body energy coupling waits for a selected
membrane/body model with energy-normalized ports; distributed collision waits
until the shared metallic hi-hat fit demonstrates a state-dependent interaction
that passive loss and driven stochastic contact cannot reproduce.
