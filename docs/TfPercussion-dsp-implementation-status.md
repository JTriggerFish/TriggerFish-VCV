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
| Contact and compact body | half-sine force pulse, tonal chirp, enveloped noise, finite micro-contact burst, gated/finite renewal micro-contact process, 2x correlated FM burst with 4x reference, explicit direct/body router | duration, endpoints, strength-to-energy mapping, deterministic seeds, requested stochastic event rate, coincident dense contacts, gate release, FM spectral convergence, bounds, routing |
| Delay and diffusion | 12-tap/2048-phase moving sinc delay, static Thiran delay, shared cubic Lagrange delay, fractional Schroeder all-pass | tone gain, low-frequency delay, polynomial exactness, integer-boundary continuity, impulse energy, five sample rates |
| Spectral motion | 255-tap antisymmetric FIR Hilbert transformer, phase-continuous signed SSB frequency shifter, and fourth-order translation-band guards | wanted level, image rejection, exact zero shift, through-zero automation, DC/Nyquist translated-content rejection, five sample rates |
| Resonance and loss | arbitrary damped modal bank, deterministic statistical modal cloud, modal passive-loss adapter, complementary three-band delay loss, orthogonal Givens mixer, projected wet-only fractional-comb network, explicit output submix | analytic modal frequency/T60 recurrence, independent projections, cloud determinism/range/normalization, passive attenuation, exact zero-coupling identity, scattering energy, delay T60 recurrence, excitation isolation, group routing |
| Cymbal bloom | slow stochastic delay, 2x oversampled bounded self-phase delay with a 4x reference implementation, serial four-all-pass dispersion loop with explicit outer feedback | 1x/2x comparison against 4x, causal onset and declared recurrence delay, zero-drive linear reduction, no hidden feedback sample, long contractive stress |
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

The dispersion return is stored in its base propagation delay. The signal is
read, passed serially through slow delay, four all-passes, self-phase delay and
loss, then written back with the new body drive. This avoids an accidental
extra sample at the outer feedback sum. The dispersion tap remains an analysis
and body-renderer drive output, not audible dry leakage.

The coupled resonator is one superset: coupling zero gives independent combs;
nonzero coupling applies an orthogonal matrix before line loss. It is not the
room-reverb velvet FDN or its coefficient set.

## Performance profile

Measurements on 2026-09-02 used an AMD Ryzen 9 PRO 8945HS, GCC 16.2 MinGW
Release build, 48 kHz, and the current default crash graph. They are development
measurements rather than cross-machine budgets.

The reproducible probes are `dev.ps1 benchmark-percussion` for native
components and `node workbench/tests/performance_probe.mjs` after
`dev.ps1 build-workbench` for WebAssembly rendering and STFT analysis.

| Path | Current |
| --- | ---: |
| Crash without modal bodies | 344 ns/sample |
| Crash with 12 sparse modes only | 362 ns/sample |
| Legacy 2,048-mode dense body | 2,497 ns/sample |
| Complete legacy crash | 2,531 ns/sample |
| Legacy 4,096-mode experiment | 4,693 ns/sample |
| Experimental unified 408-mode field | 4,036 ns/sample |
| Isolated dispersion loop | 141 ns/sample |
| Isolated 512-mode cloud | 524 ns/sample |

The complete legacy graph therefore uses about 12.1% of one core at 48 kHz;
the 408-mode unified-field experiment uses about 19.4%. This is the deliberate
quality-first workbench configuration; the anchor/packet allocation should be
profiled before it becomes a production module. Native and WebAssembly render
summaries remain closely matched, so JavaScript/Wasm boundary copying is not
the dominant expense.

The modal bank now sanitizes stable hit projections once, caches safe per-mode
damping multipliers until a damping control changes, and separates independent
state recurrence from ordered output reduction. This lets GCC vectorize the
hot recurrence and lets Emscripten use `simd128` without `-ffast-math` or a
change to finite-value recovery. A separate `-ffast-math` experiment reached
836 ns/sample, so the safe implementation is now within about 9% of that unsafe
upper bound.

The remaining optimization order is:

1. Benchmark aligned modal arrays and vector-width reductions before changing
   the number of modes. Mode-count or residual-substitution changes require
   controlled listening tests and are not the first optimization.
2. Only then inspect the smaller dispersion and observation paths.

For the browser workbench, a 10-second 2048-point, 75%-overlap STFT takes about
70 ms and stores 3.67 MiB. A 1,089 x 506 heatmap redraw takes about 22 ms. A
slider-to-completed-analysis cycle previously measured 1.69 seconds: 1.36
seconds was DSP, 220 ms was the deliberate debounce, and the remainder covered
analysis and UI painting without a main-thread task over 50 ms. The current
modal-field experiment remains compatible with asynchronous rendering. Canvas
replacement remains lower priority. Reference spectra now use an eight-entry LRU cache (about 29 MiB at
the default analysis settings). Snapshots retain audio and controls but
recompute their spectrogram on restore, reducing a ten-second snapshot from
about 5.5 MiB to about 1.8 MiB.

`RadiationFilter` is an object/observation voicing stage. Its final biquad is a
two-pole low-pass, so the asymptotic electrical slope is 12 dB/octave. It is
not, by itself, a physical distance-dependent air-absorption model: atmospheric
loss belongs in an optional propagation/observation stage and must vary with
distance, humidity, and frequency when that distinction matters.

## Current experimental graph and replacement target

`CrashCymbal` composes the tested contact, serial feedback-dispersion,
passive-loss, modal, and observation primitives. The exact C++ graph is
exposed to Python and to the optional browser workbench through WebAssembly.
Graph tests cover repeatability, strength, location, hardness, passive mute,
finiteness, and five sample rates.

The complete implemented signal flow and the unified field's equations are in
the self-contained
[nonlinear resonator architecture](TfPercussion-nonlinear-resonator-architecture.md).

An earlier coupled-comb/frequency-shift graph was rejected during calibration:
its controls could not place persistent ridges independently. The implemented
replacement originally sent direct body drive to a small arbitrary modal bank
and dispersed drive to a deterministic statistical cloud, with a separate
turbulent residual. That graph remains available as the legacy A/B path.

The experimental replacement is one 408-mode stochastic field derived from
twenty-four editable anchors in a 40 Hz--15 kHz constructive design range. Each anchor expands to a coherent centre mode and
sixteen nearby satellites. Turbulence transfers normalized excitation energy
from the centre to the satellites, widens their ERB-scaled frequency packet,
and enables energy-preserving phase diffusion. At high settings it also turns
on alternating local Givens rotations, allowing passive energy exchange within
and between neighbouring packets. A per-anchor response scaler can keep
selected ridges clean while the global control diffuses the rest. Pole radii
and the external mute controller
remain the only declared loss mechanisms. The raw dispersion tap remains
inaudible, and the field uses one body observation/radiation path. A/B controls
also expose bloom all-pass diffusion without altering its nominal delay.

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
until compact hi-hat or snare fits demonstrate a state-dependent interaction
that passive loss and driven stochastic contact cannot reproduce.
