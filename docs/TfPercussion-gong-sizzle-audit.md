# Gong upper-layer audit — 10 September

Question: the tuned body is useful, but the metallic layer is too low-pitched
and insufficiently grainy. Can the measurements and fitting objective see it?

**The diagnostics can; the last fitting objective is not sufficient.** This
audit changes neither the published preset nor its engine. It examines the
`gong-layered-final` checkpoint, exact WASM renders and the same saved reference.

## Measurements

Over 0.5–1.5 seconds, measuring only the 3–15 kHz layer:

| Measurement | Reference | Standard model seed | Four model seeds |
| --- | ---: | ---: | ---: |
| Power-weighted spectral centre | 5995 Hz | 4883 Hz | 4004–4883 Hz |
| Fraction of this layer's power above 7 kHz | 26.5% | 17.4% | 1.3–17.4% |

The absolute spectrum also shows too much emphasis around 3–4 kHz and
insufficient energy at the very top. Whole-instrument centroid would be
misleading here because the deliberately stronger bass dominates it.

Existing `ModalTextureLoss` measurements reveal frequency-dependent errors,
not a uniform need for more modulation. For example, the standard seed has
excess fast modulation near 9.5 kHz but deficient modulation near 12 kHz.
At 12 kHz, the four-seed mean normalized modulation-power deficit is about
4 dB in both the 8–32 and 32–128 Hz modulation ranges (0.6–1.6 s).

As an independent timing check, divide the analytic 12 kHz-band envelope by
its Gaussian-smoothed local trend (sigma 100 ms), **for analysis only**.
The deficit persists: 8–32 Hz modulation is 0.4–7.3 dB below reference, and
32–128 Hz modulation is 1.3–5.0 dB below, across four seeds. This supports a
texture difference beyond the late/weak arrival of high-frequency energy.
It does not establish that these two modulation bands fully describe the
listener's word “grainy.”

## Why the previous objective missed it

The latest layered-gong optimization used five broad spectral-power bands and
twelve time regions. It did **not** include the existing auditory-band texture
loss or narrow-band spectral-shape terms. An earlier general analysis contract
describes richer tests; those must not be attributed to this particular fit.

Two controlled counterexamples expose the lost information:

- Equal-amplitude 8 and 13 kHz tones have effectively identical power envelopes
  in its 7–14 kHz band despite their different pitches.
- A 10 kHz carrier with strong 40 Hz AM, RMS-compensated analytically, also has
  effectively identical broad-band envelope measurements to the plain carrier.
  The existing texture diagnostic clearly distinguishes them.

The tests use interior windows, away from onset/end boundaries. These are
counterexamples to broad-power sufficiency, not perceptual validation using
real percussion. A weighted mean over a few stochastic seeds also hid too much
variation in the location and strength of the upper energy.

## Consequence for the next fit

Keep the pitched-body and bloom-envelope constraints, but add regional upper
spectral shape/centre/roll-off and auditory-band envelope modulation. Inspect
individual bands and individual seeds, not only a single aggregate score.
Assess fine fluctuations after separating the slow bloom trend. Envelope
sparsity and cross-band correlation are worthwhile additional diagnostics;
they are missing from our simplified `ModalTextureLoss`.

This direction is supported by [McDermott and Simoncelli's auditory texture
study](https://mcdermottlab.mit.edu/papers/McDermott_Simoncelli_2011_sound_texture_synthesis.pdf),
which found that spectral power and marginal statistics alone were insufficient
and tested modulation and cross-band statistics. That work primarily addresses
stationary textures; applying its ideas regionally to a nonstationary gong is
an inference, not a validated gong-specific loss.

Do not interpret “grainier” as permission to add uncorrelated noise or turn up
global turbulence. First test the existing packet distribution/density, spread
and phase coherence against these diagnostics while protecting the low body.
Only then can a remaining mismatch be attributed confidently to the model.

## Reproduction

`tools/audit_gong_sizzle.py build/gong-layered-final` uses the development
Python environment and `EMSDK_NODE`. It verifies the checkpoint, renders seeds
1675/1982/2586/3276, and writes `sizzle-audit.json` and `sizzle.plotly.json`.
`tools/capture_fit_plot.mjs build/gong-layered-final sizzle` creates the inspected
PNG in a disposable browser tab without audio playback. Plots show the standard
seed except the explicitly labelled four-seed modulation heatmap. No separate
audition page or changed workbench preset is created.

## Follow-up fitting experiment: rejected blur-heavy candidate

The follow-up used actual six-second WASM renders, never a differentiable
surrogate. `refine_gong_sizzle.py` screened individual controls and then three
coordinate grids. `UpperSizzleLoss` compares 18 equal-ERB upper power bands in
five time regions, upper spectral centroids, and 8–32/32–128 Hz auditory-band
envelope modulation. The older screen included 2–8 Hz modulation too; its
aggregate scores must not be compared directly with the revised objective.
The current objective records its windows, weights and normalization in JSON.
Three search seeds were 1675/1982/2586. Each low-body region on each seed must
stay within 2 dB of the accepted preset, including the first 100 ms.

A joint candidate at phase blur 0.25, blur tilt 1.5, packet spread 1.8 and
noisiness slope 0.2 improved the revised search ranking from 26.71 to 15.93.
All frequencies, painted levels, damping, excitation and bloom controls stayed
unchanged. Training-seed low-region changes stayed below 1 dB. Its standard
seed upper centroid increased from 4883 to 5163 Hz, still below the reference's
5995 Hz; independent locally detrended 12 kHz modulation also became closer.

**This candidate was not published.** Inspected STFT and spectral plots showed
an overly smooth upper spectrum. A separate diagnostic subtracts an 80 Hz
smoothed log spectrum after temporal power averaging, measuring fine spectral
contrast in 3–7 and 7–14 kHz, during 0.6–1.6 seconds. Reference contrast was
5.78/5.10 dB, accepted preset 5.71/4.70, candidate 3.93/3.99 on the standard
seed. Seeds 3276/4219/5917 confirmed the smoothing; they were not search seeds.
This statistic is a diagnostic, not a validated perceptual threshold.

`polish_gong_sizzle.py` subsequently tried gentler blur and one broad painted
level adjustment (900/4000/12000 Hz). Rejecting more than 0.75 dB loss of fine
contrast in either band on any training seed left no improvement over the
accepted preset. That tolerance is an explicit design guard, not a hearing
threshold. No individual frequencies, per-mode damping, extra EQ or audio
normalization were used. Artifacts are under `build/gong-sizzle-fit`,
`build/gong-sizzle-joint` and `build/gong-sizzle-polish`.

The next question is thus not simply how much more blur to fit. The present
per-sample independent phase kicks diffuse phase and weaken coherent lines.
`SmoothModalDrift` already supplies bounded, smooth frequency wander, but its
current equal-Hz depth principally addresses slow beating. Test coherent
movement separately from stochastic linewidth before changing the fit again:
in particular bounded phase modulation can retain a carrier while adding
moving sidebands, unlike accumulating independent phase kicks. This is a
constructive synthesis proposal, not an identified physical cause in this gong.
[Fitz and Haken's bandwidth-enhanced oscillator work](https://www.cerlsoundgroup.org/Loris/ICMC95/BandwidthOscillators.html)
is relevant to this carrier/sideband distinction; it does not validate the
current fit or establish suitable control ranges.

The bounded-movement prototype is now implemented and tested; see
[Bounded ridge movement](TfPercussion-bounded-ridge-movement.md) for its controls,
equations, verification and the remaining gong mismatch. Following the user's
request, the gong preset now loads movement 1 rad / 100 changes/s / sharing 0.25
together with its existing blur; the blur-heavy candidate remains unpublished.
