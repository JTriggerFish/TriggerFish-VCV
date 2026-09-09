# Percussion review checkpoint — 9 September 2026

Fitting is paused. This is an engineering checkpoint, **not certification that
the gong/crash match their references**. The gong still has overly regular
beating and overly persistent straight ridges. No presets were retuned during
this review.

## Current design

- One modal packet body: painted centres, allocated surrounding modes, selectable
  beating layouts, phase blur and independent smooth Hz wander.
- Spectral diffusion redistributes stored energy; the visible T60 curve damps it.
- Bloom timing and Hold decay are design-time helpers. They update visible,
  saved parameters, not additional runtime envelopes or hidden gains.
- Contact/body levels feed one final EQ. Browser master/limiter remain separate.

See [architecture](TfPercussion-nonlinear-resonator-architecture.md),
[current movement/EQ](TfPercussion-modal-wander-and-eq.md), and
[Hold decay](TfPercussion-bloom-decay-coupling.md).

## Review fixes

1. Exposed packet allocation/layout, beating, blur tilt, wander and generic model
   selection flags in the optional Python native wrapper. Previously these new
   fields were available through WASM but absent from the Python binding.
2. EQ response now follows the actual live sample rate; idle uses the rendered
   spectrum's rate. Curves stop at Nyquist instead of folding above it.
3. Coarse perceptual polishing protects and scores against the actual incoming
   audio. Clipping the initial knot amplitudes to optimizer bounds must not
   silently redefine the baseline. Infeasible or worse candidates retain the
   incoming parameters. Added an exact-reference/out-of-bounds regression case.
4. Removed obsolete workbench fitting controls. Snapshot import tooling now uses
   the explicitly converted fit returned by the renderer, not stale archived keys.
5. Optional Torch tests skip cleanly without the perceptual development extras.
   A dependency-blocked collection check exercises the normal CI setup.
6. Silent modulation references produce an explicit error, not a NaN score.
7. Removed the stale one-way-cascade equations from the current architecture
   overview and corrected its diagram. The active bidirectional diffusion law
   is linked as the authoritative implementation description.

Validation: 568 Python tests, 33 native tests, 17 WASM/workbench tests, two native
API tests, native/WASM numerical agreement, silent browser checks, formatting
hooks and the normal MinGW Rack build passed. Optional-dependency-blocked test
collection also passed. These establish functionality, not perceptual fit quality.

The native wrapper retains the generic C++ defaults; the workbench recipe selects
spectral diffusion and relaxed turbulence explicitly. They are not interchangeable
parameter schemas. Use the WASM renderer for current workbench fitting and snapshots.

## Resume fitting

Start from the current workbench snapshots and actual user/reference strike.
Compare reference and candidate at fixed levels, with early/late spectra and
modulation diagnostics. Lower periodicity alone is not acceptance: added noise,
quieter ridges or a changed decay can reduce that metric while sounding worse.

The individual `tools/refine_*`, `screen_*` and dated investigation documents
record experiments, not one authoritative calibration pipeline. Archived build
artifacts belong to the model revision that produced them. Preserve originals;
explicitly import or rerender under a new output directory after schema/model
changes. Do not silently reuse their old scores as current measurements.
