# TfReverb current-architecture differentiable reference

## Purpose

The PyTorch reference searches static coefficients for the exact late-field
architecture used by the Rack module. It does not fit measured impulse
responses and it is never executed by the audio thread.

The hand-authored coefficients in
`src/tfdsp/late_reverb_coefficients.hpp` remain the immutable **Base**
flavour. A fitted artifact is exported separately as
`src/tfdsp/late_reverb_optimized_coefficients.hpp` for the future
**Optimized** flavour. Exporting can therefore never overwrite Base.

No optimization has been run for this generation yet.

## Production operator

With modulation and shimmer disabled, the recursive loop is

```text
six wall-domain input connections
        |
        P                       16 x 6 orthonormal Walsh projection
        |
16 geometry-scaled main delays Zmain
        |
three-band main-delay loss Lmain
        |
U0 -> D1 -> L1 -> U1 -> D2 -> L2 -> U2
 ^                                      |
 +--------------------------------------+
        |
       P^T -> six wall-domain outputs
```

The complete feedback operator is

\[
F(z)=U_2L_2(z)D_2(z)U_1L_1(z)D_1(z)U_0L_{main}(z),
\]

and the delayed-line state is

\[
x=(I-Z_{main}(z)F(z))^{-1}Z_{main}(z)B.
\]

There is one topology for every Decay value. There is no Householder fallback,
feed-forward substitution, or decay-dependent hybrid.

## Exact parity with C++

The optional PyTorch test suite renders the six-by-six wall-domain impulse
response through the actual C++ production loop and compares its complex FFT
against the reference at multiple control points and frequencies. This is a
numerical implementation-parity test, in addition to checking that exported
coefficient literals match. PyTorch remains an opt-in development dependency;
normal local builds and tests do not require it.

The reference mirrors these production details:

- 16 main delays expressed as ratios of the room mean free time;
- two velvet delay banks and three fixed signed-Hadamard transforms;
- the same signed permutations and signs as the Base header;
- room dimensions derived from Size and Aspect using the geometric
  interpolation used by the early-reflection engine;
- mean free time `4V/(cS)`, including the production clamps;
- Diffusion scale `0.20 + 1.30 * smoothstep(diffusion)`;
- integer rounding of every scaled velvet tap;
- the four-point Lagrange response of each fractional main-delay read head;
- complementary one-pole bands split at 220 Hz and
  `min(3500 Hz, 0.20 * sample_rate)`;
- attenuation on every physical main and velvet delay segment;
- Decay mapped to 0.25--8 seconds;
- high-band T60 `decay * 0.20^smoothstep(damping)`;
- low-band T60 `decay * 0.72^smoothstep(damping)`;
- the 0.42 injection gain and normalized six-wall Walsh projection.

Source/listener connection delays and weights sit outside the recursive loop.
They change excitation and observation but not its poles. The objective scores
both the physically reachable six-wall transfer and the complete 16-by-16
resolvent, so a problematic internal mode cannot hide outside the default
spatial subspace.

Runtime crossfades, random modulation, and shimmer are deliberately excluded.
They are time-varying presentation effects and are verified by native
time-domain tests; they must not be used to conceal a deficient static tail.

## Search space

The topology is fixed. Searchable coefficients are:

1. The 16 normalized main-delay ratios.
2. The two sets of 16 base velvet delays.
3. Three signed permutations, one after each normalized Hadamard transform.

The signed-Hadamard factorization is orthogonal by construction. Pure velvet
delays are lossless, so every discrete candidate remains paraunitary before
the explicit decay filters. Stability is consequently not a learned
constraint.

Signed permutations and integer velvet taps are screened as discrete
candidates. Main-delay ratios receive differentiable local refinement around
the selected candidate. Their mean is normalized to one after every
evaluation, preserving the meaning of room mean free time and therefore the
Size control.

Regularization rejects closely spaced main delays, repeated pairwise
path-length differences, and coincident velvet taps. It does not learn output
balance or compensate for an unstable matrix.

## Objective and coverage

The loss is reference-free and measures upward pole prominence rather than
matching a target IR. It removes local spectral means at physical 2, 8, 32,
and 128 Hz scales and measures excess prominence across octave-like bands from
20 Hz to 20 kHz. Deep transfer zeros are not penalized like upward resonances.

Every continuous proposal covers the full audible grid. The default 32-second
grid has 0.03125 Hz spacing. The control grid contains all combinations of:

- Decay: minimum, factory default, maximum;
- Damping: minimum, factory default, maximum;
- Diffusion: minimum, factory default, maximum;

at five room geometries spanning minimum/default/maximum Size and the two
Aspect extremes. This produces 135 control points. Control sub-batches reduce
peak VRAM only. A first pass obtains all per-control values and global
smooth-worst weights; a second pass accumulates their exact aggregate gradient.

Discrete screening uses the same controls over 20 Hz--20 kHz at 1 Hz spacing.
Continuous refinement uses normalized steepest descent with backtracking.
Every proposal is accepted only if the complete objective improves, and the
best checkpoint is retained.

## Artifact and flavour boundary

The optimizer writes
`data/reverb-calibration/velvet-vfm-current-v1.json`. The artifact records:

- an explicit two-stage architecture signature;
- sample rate, frequency range and grid duration;
- complete control-grid size;
- selected discrete seed;
- normalized main-delay ratios;
- two velvet-delay rows;
- three permutations and sign rows;
- optimization history and best checkpoint.

The exporter validates every dimension, permutation, sign, positive delay,
ratio ordering, and unit-mean ratio constraint. Its output namespace is
`late_reverb_optimized_coefficients`, distinct from Base. Runtime flavour
selection will be added only after a fitted artifact passes objective,
time-domain, and listening acceptance.

## Commands

Do not run the optimizer until the baseline diagnostics have been reviewed.
When authorized, the intended GPU command is:

```powershell
$env:PYTHONPATH = "python"
.venv\Scripts\python.exe tools\optimize_velvet_reverb.py `
  --device cuda --dtype float32 --steps 160 --seconds 32 `
  --control-batch 5 --candidates 128 --seed 73021 `
  --lr 0.20
```

Export remains a separate explicit action:

```powershell
$env:PYTHONPATH = "python"
.venv\Scripts\python.exe tools\export_reverb_coefficients.py
```

Compare Base and Optimized without changing the Rack runtime:

```powershell
$env:PYTHONPATH = "python"
.venv\Scripts\python.exe tools\compare_reverb_flavours.py `
  data\reverb-calibration\velvet-vfm-current-v1.json --device cuda
```

## Acceptance

An objective improvement is necessary but insufficient. Before the Optimized
flavour is exposed in the module it must also pass:

- all differentiable architecture/parity tests;
- the existing native paraunitarity and control-response tests;
- long maximum-Size/Decay impulse-response tests at every Diffusion corner;
- boundedness and non-periodicity checks with modulation and shimmer disabled;
- automation, stereo-balance, source/listener, damping, and shimmer tests;
- level-matched listening comparison against Base.

Once the runtime flavour selector is implemented, Base must remain available
so the optimizer's contribution can be auditioned directly and reversed
without changing patches. At present the module runs Base only; no optimized
artifact has passed the acceptance gate yet.
