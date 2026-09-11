# Modal motion optimisation

The current gong's two-second WASM render dropped from approximately **647 ms
to 327 ms** in warmed local benchmarks: about twice as fast, with movement still
enabled. In the browser, the complete 7.7-second render dropped from about
2.5 seconds to **1.25 seconds**; the first analysed preview arrived in **88 ms**.
These are local measurements, not cross-platform performance guarantees.

## Implementation

The previous implementation mixed regular oscillator arithmetic with random
branches and per-mode packet lookup inside the modal propagation loop. That
obstructed SIMD vectorisation. The new `BoundedModalMotion`:

1. Stores trajectory phases, endpoints and speeds in separate contiguous arrays.
2. Separates phase advance and quintic-curve evaluation from random knot changes,
   preserving the original random-draw order.
3. Evaluates shared movement once per packet, then prepares per-mode phase
   increments and rotations in regular loops.
4. Applies prepared rotations in modal propagation, removing the indirect packet
   lookup and large-angle branch from that hot loop.

This is portable C++ arranged for compiler vectorisation, shared by native
voices and both WASM targets. The added scratch arrays are fixed-size: no
real-time allocations. There is no control-rate approximation, reduced mode
count, changed seed sequence or relaxed finite/denormal handling. No fit,
damping, gain or user control changed.

The existing Cayley rotation remains algebraically unit-norm. Its small-angle
tangent polynomial and large-angle `tan()` fallback are unchanged. Preparation
selects a polynomial-only loop when a bound proves the fallback unnecessary:

$$
|\Delta\phi| \leq 1.875 \times 2 \times 1.25
\frac{r d}{f_s} = 4.6875\frac{r d}{f_s}.
$$

Here $r$ is knot rate, $d$ is phase depth and $f_s$ is sample rate. Selection
uses 0.29 radians, leaving rounding margin below the existing 0.3-radian
threshold. Low sample rates retain the per-mode fallback. This chooses
equivalent evaluation paths, not a different sound model.

## Verification

- A frozen scalar implementation in `tests/reference_bounded_modal_motion.hpp`
  checks prepared rotations at 1, 8 and 48 kHz, with 9, 17 and 513 modes and
  sharing values 0, 0.5 and 1. Odd counts exercise SIMD remainder loops.
- Existing deterministic reset, bounded-displacement, damping/energy and
  wrapped-sideband tests pass. The native percussion suite has 23 tests.
- All 21 WASM/workbench tests, two native API tests and native/WASM signature
  comparison pass. Builds use `dev.ps1` and MinGW/Emscripten.
- Gong/crash/ride at 44.1, 48 and 96 kHz and strengths 0.3 and 1 produce
  bit-identical two-second renders against the archived pre-change WASM.
  An eight-second, four-hit gong sequence is also bit-identical.

## Reproduction

`tools/compare_modal_optimization.mjs BASELINE_MODULE` compares actual PCM and
writes `build/modal-optimization-comparison.json`. The local baseline module is
under `build/motion-optimization-baseline/`.

`tools/profile_gong_render_components.mjs` supplies warmed median timings;
`tools/profile_workbench_updates.mjs` measures the served UI in a silent,
disposable browser tab. Neither changes presets or the user's playing tab.

Presets without ridge movement showed no meaningful speed change. This pass
does not claim a twofold speedup for every preset, nor that no further
optimisation is possible.
