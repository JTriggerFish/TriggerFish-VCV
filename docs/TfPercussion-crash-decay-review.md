# Crash: fit the shared T60 curve against decay, not overall loudness

The subsequent [texture review](TfPercussion-crash-texture-review.md) retains
this four-point T60 curve and records the next workbench audition.

This follows the [spectral/bloom audition](TfPercussion-crash-dynamics-review.md).
The user judged its spectral profile better but its decay too slow. The previous
full-tail audit already showed late low-mid excess: this should have prevented
presenting the decay as handled.

## What changes, and what does not

Only the existing **shared T60 curve** changes. Its points set damping in the
modal field; the C++ implementation interpolates log T60 on an ERB frequency
axis. This is not a new output envelope. All modal frequencies, painted levels,
per-mode settings, excitation, contact, turbulence, cascade, radiation and
playback gains remain fixed. More curve points add no new DSP mechanism.

Two points are a useful starting simplification, not a technical constraint.
The UI supports eight points total. Here points are added progressively where
the measured decay requires independent shaping; other presets are untouched.

## Fitting procedure

1. Render the actual workbench C++ engine against the same ten-second reference,
   Private crash A, edge, velocity 72, repeat 1, with unchanged onset and gain.
2. Fit only T60 values using `BandDecayShapeLoss` and bounded finite-difference
   least squares. Before each search, probe each control's actual influence.
3. Begin with two endpoints, then add a 600 Hz point, then a 2.5 kHz point.
   New points start on the previous curve. All values are serialized through
   the same UI snapshot path, and snapshot audio is checked against the render.
4. Inspect absolute band envelopes and signed spectrogram differences. Check
   late excess separately; improving a shape score is not a publication gate.
5. Check fresh synthesis seeds and repeated strikes before replacing the
   main workbench preset. No separate audition server or report page is created.

The shape loss uses 4096/512 STFT analysis. Each band's dB curve is anchored by
its mean during 0.2–0.5 s, so a level change cannot conceal an incorrect decay
slope. This changes **analysis coordinates only**, never playback amplitude.
Shape comparisons exclude low-confidence reference bins using a reference-only
40 dB range and a terminal-contamination margin. Absolute late-tail checks remain
necessary because the shape mask alone does not constrain all audible excess.

The upper spectrum is measured separately in 3–4.5, 4.5–6.5, 6.5–8.5,
8.5–11 and 11–15 kHz bands. The lower bands are 100–300, 300–700, 700–1500
and 1500–3000 Hz. The default library bands remain unchanged; explicit bands
are an optional analysis configuration, not extra synthesis parameters.

## Reusable inspection

`tools/plot_decay_comparison.py <candidate-directory>` writes
`decay.plotly.json`. It shows reference and model absolute band-power envelopes
with identical smoothing and a shared reference-derived floor. A gain-difference
test ensures the display does not silently normalize the model. Existing
`tools/plot_spectral_difference.py` provides the complementary signed heatmap.

Private trial artifacts and saved parameter histories are under
`build/crash-decay-v8/`. This analysis remains optional development tooling;
normal VCV builds do not require Python, Plotly or WebAssembly.

## Audition result

The workbench uses **four shared points**: 40 Hz / 27.932 s, 600 Hz / 8.309 s,
2.5 kHz / 5.220 s and 15 kHz / 1.288 s. These are the existing curve's values,
not measured output T60s: passive cascade continues redistributing energy while
the curve damps it. Only `body_decay_*` settings changed from the previous
audition. All modal frequencies and levels, contact, bloom, radiation and gains
remain byte-for-byte equivalent as parameter values.

With the same nine-band measurement, standard-strike decay-shape error changes
from 4.527 to 3.412 dB. Three previously unused synthesis seeds improve from
5.036 / 4.427 / 5.223 to 3.527 / 3.383 / 3.007 dB. All use the same reference
recording; these are not additional reference samples. At the standard strike,
300–700 Hz excess during 6–10 s falls from +7.48 to +0.53 dB.

This does **not** settle the whole calibration. The 700 Hz–3 kHz tail is now
several dB too weak in places. Full-ten-second mel error worsens from 1.706 to
2.001; it also worsens on the three review seeds. This is a decay-focused
audition trade-off, not a claim of universally improved perceptual similarity.
Both signed spectrograms and absolute band curves were inspected.

A fifth point at 8.5 kHz improves shape error only from 3.412 to 3.381 dB and is
not retained. A further trial penalizes excess power during the second after
each reference band drops 40 dB, with a 3 dB allowance. It suppresses too much
mid/upper body energy (up to roughly 10 dB below reference in important tail
regions) and worsens shape error to 3.955 dB, so it is rejected. Late excess
still needs inspection; a one-sided ceiling alone is not a sufficient objective.

The selected artifact directory is `build/crash-decay-v8/four-points/`. Exact
snapshot replay, full-duration rendering and repeated strikes pass. Eight hits
at 250 ms spacing have raw pre-limiter peaks 0.790 at the reference velocity and
2.049 at full strength; browser limiter protection remains necessary. The prior
workbench audition is recoverable from `build/crash-bloom-v7/audition/`.
