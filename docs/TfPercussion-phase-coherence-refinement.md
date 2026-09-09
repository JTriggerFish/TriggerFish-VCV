# Phase coherence, low-ring balance and gong beating

Follow-up: [visible controls and slower gong beating](TfPercussion-visible-controls-and-slow-beating.md)
supersedes this note's gong preset values and collapsed Advanced UI placement.

September 2026. Starting crash: user snapshot **crash 2**,
`63433715-ae2a-4095-95d2-687d83c355eb`. No subharmonic generator, additional
damping knots or gain matching was added.

## Why the old low ring was off

The user changed only the first handle: 125.9 → 121.7503 Hz and −20.7114 →
−25.9523 dB. Both the roughly 58-cent pitch change and 5.24 dB reduction in
prominence matter. The saved playing gesture is stronger and at a different
location than the reference cell, so it is audited separately.

Two-second Hann spectra of the reference show:

| Region | Approximately 125.25 Hz | Approximately 126.5 Hz |
|---|---:|---:|
| 0.15–2.15 s | −43.7 dB/Hz | −51.5 dB/Hz |
| 2.15–4.15 s | −63.0 dB/Hz | −53.3 dB/Hz |

The old fixed-balance pair made the **upper partner stronger throughout**.
A long aggregate spectrum missed this reversal. Its midpoint was not the
initially dominant partial. The user's quieter, lower handle is retained as
a perceptual choice, not claimed to identify a physical 121.75-Hz mode.

Weak components also occur around 80.8/82.4 Hz and 40–41 Hz. These are not a
clean octave below 125 Hz; a recording alone does not establish their physical
origin. An additional 81-Hz handle was tested but changed the front/allocation
too much to select. Existing modal editing is sufficient to explore it.
The spectra have approximately 0.5-Hz Fourier resolution and 0.125-Hz
interpolated bins; zero padding is not extra resolving power.

## Research interpretation

[Legge & Fletcher (1989)](https://www.phys.unsw.edu.au/music/people/publications/Leggeetal1989.pdf)
report nearby peaks evolving in frequency and amplitude, and subharmonic and
chaotic regimes at strong excitation. Their experiments largely concern
axisymmetric modes and sinusoidal forcing, not this edge-struck recording.

[Chaigne, Touzé & Thomas (2005)](https://www.jstage.jst.go.jp/article/ast/26/5/26_5_403/_pdf/-char/en)
study nonlinear gong/cymbal vibration and modal coupling. This supports
distinguishing nonlinear motion from a stationary noise bed, not a universal
bass-heavy phase-noise law.

[Skare & Abel (2019), §§4.1–4.3](https://www.dafx.de/paper-archive/2019/DAFx2019_paper_48.pdf)
describe beating from nearby modes with differing decays, and frequency shifts
across excitation levels. A fixed-balance pair with nearly identical damping
cannot reproduce every such evolution.

**Constructive inference:** independent frequency balance of phase coherence
is useful. The tilt below is a design control, not a measured material law or
a substitute for all nonlinear effects.

## One new control: Blur tilt

For local noisiness $I_i$, Phase blur $b$, the existing centre/satellite factor
$q_i$ (0.35/1), and tilt $s$:

$$ B_i=q_i I_i^2 b\,\operatorname{ERB}(f_i)(f_i/1000)^s,
\qquad -2\leq s\leq2. $$

Negative tilt softens lower rings while keeping upper ridges clearer. Zero
retains the previous profile, which already grows strongly with frequency
under a positive noisiness slope. This changes no mode placement, allocation,
excitation/observation gain, transport or damping coefficient. Computation is
at preparation time; no new audio-rate processing. Zero noisiness/blur stays
coherent; the prepared rotation retains its existing sample-rate bandwidth bound.

`field_phase_tilt` appears under Modal body → Advanced trajectory as **Blur tilt
(1 kHz pivot)**. C++, UI and JSON store the same value. The metallic surface has
186 parameters, still API V1. Importing the preceding V1 surface explicitly
adds neutral tilt zero; unrelated missing parameters still fail validation.

## Doublet controls

Beat depth and rate tilt now also apply to **Beating doublets**. Previously that
layout imposed equal partners and a common gap, with those controls inactive.

A complete doublet multiplies its former equal weights by
$d\sqrt{2/(1+d^2)}$ and $\sqrt{2/(1+d^2)}$. Their squared weights still sum to two.
An odd leftover tone retains weight one. Zero depth suppresses one partner
without reallocating the pool; unrelated nearby tones may still beat.
Separation follows the existing $r(f/125)^t$ law, bounded to 0–80 Hz.
This varies speed across packets, not randomly over time. Smooth drift remains
available but is **off in both selected fits**.

Old doublet snapshots import with depth 1 and tilt 0, preserving their formerly
implicit equal weights and constant gap as visible, serializable settings.

## Selected fits and validation

**Crash — user low ring, clearer texture:** preserve the user's two handle
edits; set Blur tilt −1.5 and Packet spread 2.0 (previously 2.6225). Other sound
parameters, including damping, stay fixed. Four-seed reference-floor Mel
distance improves 1.050 → 1.006 and decay error 3.181 → 3.135 dB. Upper texture
distance is essentially unchanged, 0.3723 → 0.3729. At the user's stronger
gesture the first 450 ms changes 1.67 dB RMS; peak −3.10 → −3.34 dBFS.

**Gong — varied gentle beating:** retain series, bloom, T60 and levels; use
depth 0.5, rate tilt +0.25 and blur 0.012. Rate remains 3 Hz at 125 Hz.
Four-seed Mel distance improves 0.796 → 0.780, texture 0.328 → 0.322 and decay
error 4.309 → 4.222 dB. Shared modulation-line concentration falls 0.167 → 0.081.

`modulation_signature.py` separates depth, periodicity, shared spectral lines
and envelope correlation. Correlation alone is not an acceptance test: common
bloom and analysis boundaries also influence it. The shared-line metric above
1 Hz complements it; sub-hertz beating needs the separate long-window diagnostic.
Known-tone tests cover shallow/deep and same/different-rate beating, silence
and invalid inputs. Native tests cover unchanged placement, damping and energy
for blur tilt, plus doublet energy normalization including odd counts.

Actual Wasm checks use four seeds, fixed reference gain, the saved strong
gesture and repeated quarter/rapid-hard hits. Tools: `inspect_metal_texture.py`,
`review_texture_september.py`, `review_metal_refit.py`. Artifacts are under
`build/{crash,gong}-texture-september`; the main workbench is the listening UI.
Fixed-scale plots still show crash ridge-placement errors/excess late mids,
and excess late gong treble. Neither fit is claimed to be perceptually solved.
