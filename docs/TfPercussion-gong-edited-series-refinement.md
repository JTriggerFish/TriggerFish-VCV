# Refining the user's edited gong without fitting individual ridges

## Published result

Historical result below. The newer user edit and current controlled trials are
documented in [Gong upper-mid study](TfPercussion-gong-upper-mid-study.md).

At this stage the main Gong target loaded **Gong — edited series, refined dynamics**.
The Texture trials menu contains that candidate and the unchanged user edit,
`gong ridge movements etdited`. The rejected cloud/older texture experiments
are no longer in the active picker; their local renders remain archived.

Only two parameters differ from the user's snapshot:

| Control | User edit | Candidate |
|---|---:|---:|
| Bloom rate | 4 | 3.2 |
| High-frequency T60 endpoint | 2.14 s | 1.75 s |

Every modal frequency, prominence, local noisiness and allocation is unchanged.
So are global texture, movement, excitation and gain controls. Phase blur and
output EQ remain off. The low T60 endpoint remains 10 s, with no interior knots.
No new DSP or runtime controls were needed.

The source snapshot is `28bdb36c-56b7-4dad-b0fb-fca74788e6d7`. Its strike is
strength 0.9286317627, location 0.5031711669, mallet, hardness .35, seed 1695.
This differs from the reference cell's nominal strength .76, location .55 and
seed 1675. Rendering the full saved snapshot avoids silently substituting that
factory gesture. The source file in Documents is not modified.

## What was fitted, and how

This is a refinement of an approved **sound-design starting point**, not a
claim that its deliberately strong low body matches the recording. The target
preserves that body's band envelopes below 800 Hz, transitions smoothly over
800–1800 Hz, and uses Dresden Gong03 above 1800 Hz. Reference gain stays −6 dB.

The search uses the actual C++ WASM renderer. A coarse one-control-at-a-time
screen tests shared T60 endpoints, bloom rate/concentration dependence, movement
depth/rate and packet spread. A subsequent 20-point joint grid tests bloom
rate × high T60, retaining the original texture. Each trial is rendered at
seeds 1695 and 1982 separately; their scores are averaged, not their audio.
The original is always a candidate. No optimizer has a per-ridge coordinate.

The objective separates an overall band-level error from its time evolution:

$$e_{bt}=L_{bt}^{\mathrm{model}}-L_{bt}^{\mathrm{target}},\qquad
\bar e_b=\frac{1}{T}\sum_t e_{bt}.$$

Envelope-shape error is the RMS of $e_{bt}-\bar e_b$. The score combines it
in quadrature with the absolute low-body envelope error, preserving the user's
starting body. Bands and time cells are those of `SpectralBloomLoss`: 24
log-spaced bands and 15 regions separating attack, rise and tail over 6 s.
**This is loss decomposition, not audio normalization.** Raw level errors
remain in the audit; all playback/render gains stay fixed. It prevents the
search from shortening decay solely to turn down the user's louder upper end.

Numerical selection is also constrained by the independent finer-frequency
regional spectrum audit: upper decay-shape RMS must not worsen on either seed.
The lowest pooled score, rate 2.8/T60 1.5 s, failed this check and visibly cut
the 9–14 kHz tail too quickly. Rate 3.2/T60 1.75 s is the best eligible grid
point. A small movement-depth reduction was likewise not retained: its tiny
timing improvement moved ridge contrast away from the reference.

## Evidence and limits

At the saved seed, upper regional shape error improves 2.59 → 2.30 dB and
absolute upper error 8.77 → 6.13 dB. Low-body envelope change is about .55 dB
RMS; there is no erased low ridge. These numbers are diagnostics, not a
perceptual acceptance certificate. The inspected band plots show a closer
high tail, but the upper bloom is still early/strong in places.

Holdout seeds 2673/3911 also improve at the saved gesture. At the nominal
strength .76, aggregate timing improves but upper absolute level gets worse
(roughly 3.6–4.3 dB below the reference on those holdouts). Do not label this
a completed multi-velocity calibration. The user's exact saved gesture is
loaded with the candidate so the workbench plays the sound actually refined.

Checks: exact saved-snapshot/WAV replay; all modal controls unchanged; finite,
increasing energy at strengths .3/.5/.76/1; single/four-hit renders through the
actual 48 kHz browser limiter. No limiting occurred at Master −12 or 0 dB in
those tests (four-hit peak about −2.29 dBFS at Master 0). Browser save/load tests
also check the saved gesture, not just parameter values.

## Reproduction

- `tools/refine_edited_gong.py <fit.json> --output <directory>`: coarse screen.
- Add `--fine`: joint dynamics grid with original texture retained.
- `--fine --reuse <previous-fine-audit.json>`: reselect a validated existing
  grid, then render/export/verify the selected snapshot again.
- `tools/audit_edited_gong.py <candidate-directory>`: held-out gestures/seeds.
- `tools/audit_fit_limiter.mjs <candidate.fit.json> <output>`: safety audit.
- `tools/plotly_png.mjs <figure.plotly.json> <image.png>`: silent, isolated
  diagnostic image; the main workbench remains the only listening page.

Use the project's analysis environment and `EMSDK_NODE`; ordinary Rack builds
need none of these tools. Local evidence is in
`build/gong-edited-refinement/{fine,selected}`. Build the served UI through
`dev.ps1 build-workbench`; do not start another server or reload a user's
unsaved browser edits automatically.
