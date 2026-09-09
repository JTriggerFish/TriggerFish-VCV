# Paired ring: simple control surface, shared modal engine

Current extension: [phase coherence and doublet refinement](TfPercussion-phase-coherence-refinement.md).
Beat depth/rate tilt now also work on Beating doublets. Historical preset values
and the 185-control count below describe the preceding revision, not the current
186-control surface with Blur tilt.

Implemented after the [low-ring beating investigation](TfPercussion-low-ring-beating.md).
It is another character of the existing modal packets, not an extra synth,
LFO, output tremolo or energy source.

## Using it

Choose **Ring character → Paired ring**, then set **Beat rate**. Around
0.5–2 Hz gives slow breathing; higher values give faster shimmer. Zero removes
the centre split. The slider gives most of its travel to slow rates, with
two-decimal display precision below 10 Hz. Other layouts remain available;
the same rate control applies to surrounding pairs in Beating doublets and
is greyed out for Scattered/Even coverage.

**Beat depth** now independently controls the main pair's amplitude ratio.
Try 0.1–0.2 for subtle movement, 0.3 for a clearer pulse. The original fixed
0.5 ratio could produce 9.54 dB peak-to-trough variation in an isolated pair;
it was not gentle enough for the crash. At 0.3 that theoretical swing is
5.38 dB, and at 0.2 it is 3.52 dB. The full packet's audible depth also depends
on its surrounding tones and phase blur. Zero depth gives an unsplit main
ring, without removing or reallocating the surrounding cloud.

**Advanced → Beat rate tilt** changes speed with frequency, anchored at 125 Hz:

$$
\Delta f(f)=\operatorname{clip}_{[0,80]}\left(r(f/125)^{s}\right).
$$

Here $r$ is Beat rate and $s$ is Beat rate tilt. Zero tilt retains a shared
rate; +0.5 doubles it every two frequency octaves. This is an explicit
constructive control, not a fitted physical law. The selected modal handle
displays its resulting rate. Both new controls are inactive outside Paired ring.

Spread still sets surrounding-tone spacing. Noisiness controls how clear the
main ring is relative to that shimmer. Clear lows and noisier highs therefore
allow a slow low ring beneath a broader crash without extra low/high controls.
Local noisiness zero leaves a clear pair in this character, not a single tone.

`field_distribution` accepts 3; `field_doublet_split` is displayed as Beat
rate and accepts zero. The depth/tilt refinement adds `field_beat_depth` and
`field_beat_rate_tilt`, with matching C++, JSON and UI controls (185 total).
API remains V1. Workbench creation and double-click reset now use Paired ring,
1.25 Hz, depth 0.3 and tilt +0.25, matching the current crash's beating settings.
The depth slider uses a square-law mapping, giving over half its travel to
0–0.3 without reducing the available range or changing the stored parameter.

Ride and hi-hat keep their saved Scattered character; their inactive pairing
values now hold these useful starting settings. Gong retains Beating doublets
and its active 6-Hz satellite split, so its approved texture is not retuned.
Only its inactive paired-depth/tilt values change. Selecting a different
character never silently overwrites a user's rate. Double-click Beat rate to
restore 1.25 Hz when changing from the gong's doublets to a slow paired ring.

## Energy and observation

Each active handle reserves two centre states for nonzero beat rate; remaining
states are allocated to surrounding tones as before, within the same 512-state
budget. Both centre states belong to the **same** transport cell. The painted
frequency is their midpoint; the requested separation compresses only at the
represented frequency boundaries. Surrounding tones retain scattered placement.

The original fixed 2:1 balance has been replaced by the explicit depth $d$.
The weaker/stronger amplitude ratio is $d$, constrained to $[0,1]$.

An orthogonal launch simultaneously preserves excitation energy and the
unsplit ring's initial observation. In units of the original core input:

$$
b_- = \frac{d^2+jd}{1+d^2}, \qquad b_+ = \frac{1-jd}{1+d^2},
$$

$$
b_-+b_+=1, \qquad |b_-|^2+|b_+|^2=1.
$$

Their amplitudes are $d/\sqrt{1+d^2}$ and $1/\sqrt{1+d^2}$, at a quadrature phase
difference. These phases are part of this fixed normalized construction, not
free hidden fit parameters. With no stochastic blur, the zero-separation limit
sums to the original coherent oscillator. At exactly zero, one centre is used.
Beat maxima can exceed the single-ring amplitude, but do not create energy.
The same damping, mute and diffusion processes act on the stored modes.

This remains constructive, not a complete physical model of degenerate mode
orientation. A smooth rate tilt also cannot reproduce every irregular local
split in a real instrument; independent per-packet rate/depth controls are
deliberately not added. Adjacent-mode beating and different decay rates are
discussed in [Skare & Abel, §4.1](https://www.dafx.de/paper-archive/2019/DAFx2019_paper_48.pdf);
their work does not prescribe a common split or this power-law tilt.

## Numerical correction found by testing

The direct audio test exposed pitch error from reconstructing a small sine
from a rounded cosine during oscillator preparation. The cancellation-safe
region now extends to a smaller component below 0.1 rather than 0.01. Away
from axes, the existing unit-norm reconstruction is retained for energy
precision under stochastic phase rotation. Nothing extra runs per sample.
Both the analytic paired-tone test and existing phase-energy tests pass
without relaxing their tolerances.

This correction also affects existing voices slightly: the unchanged gong
snapshot's six-second RMS changes by −0.0011 dB; waveform-difference RMS is
−37.55 dB relative to its original render. Its preset and topology are
unchanged, but the final engine is not bit-identical to the old engine.

## Original paired-ring trial

The first trial was **Crash — paired ring trial**, with 1.25-Hz beating.
The lowest painted centre moves from 120 to 125.9 Hz, the midpoint suggested
by the reference peaks. A smooth low-end clarity edit is stored in existing
local-noisiness controls, fading back to unchanged settings by 900 Hz.
The other frequencies, levels, contact, diffusion, gains and two-point decay
curve remain fixed. There are still 24 painted handles.

`tools/fit_paired_ring_crash.py` screens four rates and four low-clarity scales.
Selection first follows the reference's measured slow pulse, then low-band
texture; a broad-band scalar score is not allowed to silently select a faster
beat. Other spectral/decay scores are reported separately, not claimed to all
improve. Additional phase seeds and repeated strikes are audited with the
same exact Wasm renderer. This is an audition trial, not listening approval.

Across four tested seeds, the 90–180 Hz envelope's strongest pulse changes
from 6.5 Hz to 1.25 Hz, matching the reference's measured pulse. Rapid
modulation's share drops from about 99.5% to 11.7–12.7% (reference: 6.8%).
This is a targeted improvement, not an overall fit victory: the broad-band
spectral score slightly worsens, and decay still needs listening review.

Tests cover pair spacing/midpoint, normalized input energy, zero rate,
boundary frequencies, analytic damped audio, pool limits, UI rate mapping,
serialization, inactive controls and browser loading. Normal Rack builds do
not depend on the web or analysis tools.

## Gentle-ring and body-balance revision

The workbench now uses **Crash — gentle ring and fuller body**. It keeps the
1.25-Hz base rate, uses depth 0.3 and tilt +0.25. A nine-candidate depth/tilt
screen preceded a three-way smooth middle-frequency warp comparison (0.9,
1.0, 1.1); the selected 1.1 warp leaves frequencies below 180 Hz and above
4 kHz unchanged. No independent ridge placement was fitted.

Six broad observation coordinates at 125, 240, 450, 1000, 3500 and 15000 Hz
were fitted using exact affine render bases, validated against actual Wasm
renders. Powell search used two seeds and 160 evaluations per geometry.
The score was reference-floor Mel + 0.3 attack Mel + 0.05 band-envelope error.
The reference floor is 60 dB; attack is 0–300 ms. This is a proposal score,
not a perceptual guarantee. The selected coordinates are baked into the
ordinary painted bars. No new runtime EQ, gain normalization, damping points,
contact, diffusion or velocity changes were introduced.

The standard hit's 90–180 Hz tail rises approximately 6 dB, while 180–320 Hz
falls 6.6 dB. Measured slow low-band modulation power falls from 0.0946 to
0.0400; the reference measures 0.0901. Thus the gentler setting follows the
audition feedback, not a claim that this one diagnostic now matches exactly.
The relative fast-modulation fraction rises as slow modulation is reduced;
absolute fast power does not rise. Across four seeds, reference-floor Mel
falls from 0.933/0.995/0.995/1.002 to 0.858/0.934/0.918/0.939. Individual
ridge differences and excessive late low-band decay remain visible.

Artifacts: `build/crash-gentle-ring`, script: `tools/refine_crash_ring_balance.py`.
The old patch and render are retained under `before`. The gong/ride/hi-hat
are not retuned; inactive pairing defaults are updated as described above.

Verification: 22 native percussion tests, 560 Python tests, 14 Wasm tests
and two native API tests passed, together with native/Wasm parity and a
silent browser check of all four metallic presets and the new controls.
