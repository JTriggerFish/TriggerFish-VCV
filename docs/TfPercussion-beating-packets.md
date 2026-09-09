# Beating packets: oscillator placement, blur and fitting

See [phase coherence refinement](TfPercussion-phase-coherence-refinement.md)
for the current blur balance, doublet controls and revised calibration values.

For the added paired-centre character and its simple Beat rate control, see
[Paired ring](TfPercussion-paired-ring.md). The three original layouts below
remain available.

This experiment starts from the user's `gong test` snapshot dated 2026-09-09.
Its 32 centres form one series: root 120 Hz, four protected harmonics, upper
stretch approximately 0.48, equal -8 dB observation bars. The original snapshot
and its exact old-engine audio are archived in `build/gong-user-beating-baseline`.
Fitting uses the reference's standard gesture; the user's slightly different
gesture is retained separately in that archive. No reference/master gain changes.

## Why beating rather than simply more noise?

[Perrin et al., *The normal modes of cymbals* (2008)](https://www.ioa.org.uk/system/files/proceedings/r_perrin_gm_swallowe_sa_zietlow_tr_moore_the_normal_modes_of_cymbals.pdf)
report many modes with split partners, consistent with small symmetry breaking.
This motivates clusters of stable resonances. It does not establish a harmonic
series or a universal separation in Hz. [Their small-gong study (2014)](https://scholarship.rollins.edu/as_facpub/122/)
also relates measured modes to perturbed axial symmetry and reports nonlinear
components. Results from those instruments are not measurements of our Dresden gong.

[Skare and Abel (DAFx 2019)](https://www.dafx.de/paper-archive/2019/DAFx2019_paper_48.pdf)
use a large bank of high-Q complex resonators and discuss approximate coupling
for delayed high-frequency emergence. This supports distinguishing modal density
from noisy broadening; it does not imply our 512 states are sufficient for every
cymbal or validate our particular energy-transport approximation.

[McDermott and Simoncelli (2011)](https://mcdermottlab.mit.edu/papers/McDermott_Simoncelli_2011_sound_texture_synthesis.pdf)
show the importance of auditory-envelope statistics and modulation structure
beyond a time-averaged spectrum. Their textures are substantially stationary;
we therefore measure short attack/bloom/tail regions separately. Our diagnostic
is inspired by that principle, **not a reproduction of their auditory model or
a validated perceptual loss**.

## Engine and controls

The energy path is unchanged: normalized contact input → stored modal packets
↔ passive spectral-energy diffusion → shared frequency-dependent damping →
painted observation levels → output filters. No new energy source, gain scaling,
per-mode decay or random arrival-phase process has been added.

The old allocator tied state count directly to packet width and overlap. Narrow
packets could not use a large state count, even when stable beating was wanted.
Now **Satellite density** requests a fraction of the available pair budget:

$$
P=\operatorname{round}\left[d\left\lfloor\frac{512-H}{2}\right\rfloor\right],
\qquad q_i=\sqrt{w_i}\,a_i.
$$

$H$ is the active handle count, $w_i$ its ERB spread and $a_i$ its visible local
allocation weight (0–4). Largest-remainder apportionment distributes $P$ pairs
proportionally to $q_i$. Zero width or zero allocation receives no satellites.
Every active handle retains its centre. The pool may leave one slot unused
because satellites are pairs. Reducing one allocation releases budget to others;
uniformly scaling all weights has no effect. The editor reports pool/handle counts.

The three **Sideband layouts** are constructive options, not measured material laws:

- **Scattered:** nested low-discrepancy radii with small deterministic jitter.
- **Even coverage:** the nested radii without jitter. Increasing count fills
  gaps without moving existing oscillator frequencies.
- **Beating doublets:** nested cluster positions, each containing close pairs.
  **Doublet separation** specifies their difference in Hz rather than cents;
  e.g. 6 Hz can produce a 6-Hz beat. It is greyed out for the other layouts.

Frequency support compresses near the represented boundaries instead of clipping
large groups onto them. Doublet separation also compresses near boundaries.
Oscillator phases are deterministic per packet, so reallocating another packet
does not change this packet's phase sequence. Adding/removing handles can still
change packet indexing. Odd pair counts can leave one incomplete doublet group.

**Phase blur** remains the existing passive stochastic phase process, with extra
slider precision near zero. At zero, sidebands are stable damped oscillators;
they can beat without any noise process. Its actual linewidth still follows the
global/local noisiness profile. Width/share and phase blur are not perceptually
orthogonal, but count no longer requires wider packets or more blur.

All new controls are compiled descriptors, serialized in JSON and exposed in the
UI. Generating a fresh series resets local allocation weights to one. Defaults
are explicitly added to factory metallic JSON files; this is an engine experiment
and old metallic presets need auditioning under the changed allocator, not an
assertion of bit-compatible sound. Snapshot loading does not secretly fit anything.

## Why small diffusion exponents are useful

The current conductivity is proportional to

$$
\kappa\left[(\rho_i^2+\rho_i\rho_{i+1}+\rho_{i+1}^2)/3\right]^a.
$$

When the bracket is below one, raising $a$ **reduces** transport. At large density
it can increase transport instead. Lower $a$ therefore helps the weak upper
front travel; it is not an unconditional "more bloom" control. At $a=0$ transport
is energy-independent. This behaviour follows the implemented law, not a latch.

[Humbert et al.'s plate model](https://arxiv.org/html/1709.09884) motivates the
quadratic endpoint $a=1$, not the adjustable exponent or our chosen energy units.
This pass keeps the equation unchanged and maps UI travel $u$ to $a=u^3$.
The displayed/saved number is still the actual exponent. Existing values retain
their meaning; there is no clamp to a narrow fitting range or hidden floor.

## Analysis and structured fitting

`ModalTextureLoss` uses smooth analytic auditory-width band filters, envelope
downsampling, and regional envelope modulation power in 2–8, 8–32 and 32–128 Hz,
plus modulation concentration. A reference-fixed activity mask avoids comparing
silent bands. Its texture features are gain-invariant and must accompany
fixed-level envelope/decay losses. Tests distinguish beat rates and noise while
tolerating a small common carrier shift. This diagnostic does not identify
individual high modes, and cannot by itself prove realism.

`tools/fit_structured_metal_texture.py` screens whole families (crash), then
layout/density/phase blur, and uses bounded Powell search of shared excitation,
diffusion, texture and two T60 endpoints. Four broad prominence coordinates are
only accepted if the combined objective improves. Individual frequencies,
high-frequency bars, local allocation weights and per-mode decays are not free
optimization variables. No added T60 knots in this pass.

The explicit search score is bloom-envelope norm/10 + 0.6 texture distance +
0.15 Mel MRSTFT. These experimental weights are not calibrated just-noticeable
differences. Component scores, untouched reference/gesture/gains and additional
seed comparisons are saved. A scalar win remains a proposal until plots,
restrikes and audition have been checked. The user's better-sounding starting
point is valuable evidence against accepting old spectral-loss wins blindly.

Build and tests use `dev.ps1` (MinGW and optional Wasm). Python fitting and
browser tools remain separate from ordinary Rack/release builds.

## Selected workbench trials (2026-09-09)

The main workbench's gong and crash reference targets now load the selected
trials. These are **not listening-approved calibrations**. The original user
snapshot/audio and every search checkpoint remain archived, not overwritten.

The gong retains all 32 original centre frequencies (120-Hz root, protected
four-harmonic core, stretch 0.48). Root changes and low-frequency allocation
ramps were tested and rejected. It uses beating doublets at 6-Hz separation,
density 0.45 (248 oscillators), spread 3.125 ERB and phase blur 0.00689 ERB.
The diffusion rate is 3.459, exponent 0.04987, and endpoint T60 values are
9.164/2.186 seconds. The crash uses one 24-centre family (120 Hz, four protected
harmonics, stretch 0.7), scattered placement, density 0.85 (438 oscillators),
spread 2.189 ERB and phase blur 0.08422 ERB. Its rate/exponent are 1.661/0.06392
and endpoint T60 values 25.114/0.600 seconds. All local allocation weights are
one, and no individual frequency or per-mode decay was fitted.

`tools/refine_texture_candidate.py` adds six broad observation-amplitude
coordinates at 120, 400, 1000, 2500, 6500 and 15000 Hz. These produce the visible
painted levels, not hidden EQ. Two-endpoint damping was compared with one extra
3000-Hz knot. Neither justified the extra control; both selected trials keep
only the endpoints. Powell searches use the actual Wasm renderer, two seeds
and fixed reference/event/gain settings. `search.json` stores the parameter
vector, objective specification, history and renderer provenance.

Independent checks use three additional seeds and optional WaveSpin JTFS.
This JTFS configuration resamples to 16 kHz and therefore **does not assess
content above 8 kHz**. It was not part of the optimization. Standard-seed
results below compare starting parameters and selected parameters on the same
new engine; they are not a comparison against the archived old-engine audio.
Lower values are better, but the rows have different units and cannot be summed.

| Diagnostic | Gong start → trial | Crash start → trial |
|---|---:|---:|
| Mel MRSTFT | 2.188 → 1.268 | 1.399 → 1.600 |
| Band decay shape error, dB | 8.008 → 4.783 | 3.718 → 3.668 |
| Envelope-modulation texture | 0.290 → 0.309 | 0.528 → 0.283 |
| JTFS below 8 kHz | 0.172 → 0.065 | 0.103 → 0.110 |

The gong's attack/bloom/decay match improves across the extra seeds, while its
texture diagnostic is slightly worse. Its 3–6-kHz rise still starts too early
and the upper tail persists too long. The crash is explicitly a **texture
trial**, not an across-the-board improvement: its 300–700-Hz attack is weak,
700–1500 Hz is excessive, and its low tail is too short. The flat late upper
reference floor is not evidence that extra synthetic noise should be added.
The paired STFT/difference and absolute band-decay plots were inspected; these
mismatches must remain visible rather than be concealed by the combined score.

Quarter-note and eight rapid full-strength restrikes remain finite and build
energy. Raw rapid-strike peaks reach +4.42 dBFS (gong) and +7.87 dBFS (crash),
before browser master volume/limiter. They are **not** unity-headroom renders;
float WAVs preserve those peaks without clipping. No automatic gain matching
was added to hide this, and the browser safety limiter remains enabled.

Artifacts: `build/{gong,crash}-texture-finish/candidate/` contains the selected
JSON, audio, paired/difference plots, band-decay plots and `refit-audit.json`.
These development artifacts and references are not dependencies of Rack builds.
