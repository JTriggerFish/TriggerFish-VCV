# EQ-free gong texture comparisons

The current picker and results are documented in
[Gong upper-mid study](TfPercussion-gong-upper-mid-study.md).
The trials below and the subsequent
[edited-series refinement](TfPercussion-gong-edited-series-refinement.md)
are historical.

## Listening correction: structured series, not independently fitted ridges

The user rejected the cloud's resonance/tuning and preferred ridge movement.
Those listening results supersede any favourable numerical ranking below.
The published trial files remain historical candidates, not approved fits.

The user's `gong-ridge-movements-etdited.json` (snapshot
`28bdb36c-56b7-4dad-b0fb-fca74788e6d7`, parent movement trial
`055e2c5f-f9f6-4050-ad18-74408cc356d5`) changes only painted observation
levels in the instrument patch. Frequencies, texture, transport and damping
are unchanged. The edited levels form a much smoother rising upper envelope,
instead of the previous alternating peaks and troughs. This is evidence for
improving the level parameterization before inventing another texture mechanism;
it is not an isolated experiment about frequency spacing.

**Current tuning rule:** do not optimize individual ridge frequencies, gains,
or decay multipliers. Start from the user's edited series. Any subsequent
frequency fitting must use shared series parameters (root/stretch), not a
peak-picking list. Texture and shared damping remain separate search stages.
Do not silently reconstruct the edited series from an older preset.

`polish_gong_texture_trials.py` now has only two observation search coordinates:
common gain and tilt in dB/octave around 1 kHz. Both start at zero. Their bounds
are ±12 dB and ±3 dB/octave, with additional linear constraints keeping all
active bars within their visible range. No individual clipping, gain nullspace,
per-ridge optimizer variable, or runtime parameter is introduced: results bake
into the ordinary editable JSON bar values. The incoming shape is otherwise
preserved; this cannot fix a bad starting envelope. Tests check that an isolated
ridge correction is outside this two-dimensional search space.

The movement family is now the default for the polishing/dynamics tools.
The old independent-gain results below are **not** evidence that these new
constraints have already produced a better fit. No new fit has been published
as part of this correction.

Physics supports structured families, but not a universal equal-spacing or
equal-volume rule: measured cymbal modes can be described by modified Chladni
families ([Wilbur and Rossing, 1996](https://www.auditory.org/asamtgs/asa96haw/5pMU/5pMU8.html)).
Constant *average* plate modal density does not imply an equally spaced ladder.
A smooth prominence envelope is our constructive/perceptual fitting prior,
not a claim that every real mode radiates equally strongly. A single flat-Hz
cloud kept only an average density idea and also changed the transport geometry;
that was insufficient to preserve the desired metallic tuning.

## Original experiment (historical)

The main workbench has a **Texture trials** menu beside the reference targets.
Each entry restores a complete editable patch and the same Dresden Gong03
reference/gesture. Nothing plays automatically. Use **Reference**, **Trigger
synth**, or the strike surface. Snapshot edits before changing trials.

| Trial | Mechanism | Painted handles |
|---|---|---:|
| Gong — phase blur only | Phase blur .025 ERB, tilt 0; ridge movement off | 32 |
| Gong — ridge movement only | Blur off; movement 1.5 rad, 200 changes/s, sharing .15 | 32 |
| Gong — one upper cloud | Blur off; movement 1.5 rad, 120 changes/s; one wide 9 kHz upper packet | 16 |

All three disable pitch wander too, keeping that separate mechanism out of the
comparison. Existing beating remains available. These are **independently
retuned listening candidates**, not identical-parameter on/off comparisons.
The normal saved gong now also bypasses final output EQ. The optional EQ module
remains available for sound design; no EQ setting is fitted in this experiment.

## What the cloud changes

The optional **Plate cloud** ring character samples stable side resonances
approximately uniformly in Hz, instead of approximately uniformly in ERB rate.
It retains the paired central ring and the existing normalized energy budget.
The surrounding resonances are not freshly generated noise grains.

The cloud uses the existing controls: centre frequency, prominence, local
noisiness, global packet spread, satellite density and local sideband allocation.
It introduces **no extra sliders**. The final trial keeps 15 lower handles and
replaces the upper collection with one 9 kHz handle controlling roughly 300
oscillators. Its allocation is 4; other handles are .2. The total pool remains
512 states, not 512 per handle. A unit test also checks that a single handle
can use essentially the entire pool by itself.

For this layout, packet half-width is spread times the local ERB bandwidth,
measured in Hz. Each side compresses at the supported frequency boundaries;
there is no pile-up of clamped frequencies. The editor's width handles use the
same Hz mapping. Its translucent bell remains an illustrative packet display,
not the measured frequency response.

Every side mode retains its own T60 evaluated from the shared curve. The packet
shares one energy-transport coordinate: fewer upper handles therefore change
bloom timing as well as the editing surface. Retuning the cloud's diffusion
strength from 4 to 2.5 improved its initially premature high-frequency arrival.
It still does not reproduce the recording's whole time/frequency evolution.

## Ridge movement range

Ridge movement now spans 0–3 radians; speed remains .1–200 changes/s. Higher
depth creates stronger moving sidebands and can substantially reduce the central
ridge. It is not a free density increase with unchanged tone. The bounded phase
trajectory remains energy-neutral and cannot accumulate a pitch walk.

The 3-radian maximum passes the existing 15 kHz sideband-wrap test at 32, 44.1,
48 and 96 kHz (below −60 dB wrapped energy in that test). This is a measured test
case, not a claim of strict mathematical bandlimiting for every input.

## Fitting and checks

The explicit task is to retain the current low body while refining the upper
texture against the recording. The target uses the incoming **EQ-bypassed**
body below 800 Hz, transitions smoothly in log frequency, and uses the real
reference above 1.8 kHz. This is deliberately a hybrid sound-design target—not
a claim that the stronger low body matches the recording exactly.

We screened 52 texture settings, then varied shared bloom/damping and cloud
placement/allocation. Actual painted amplitudes are optimized independently,
with a soft .15-times-RMS dB-change penalty. There is no hard-interpolated gain
curve, protected bass bar, per-mode damping fit, extra T60 knot, or audio
normalization. The shared two-endpoint T60 curve remained 10 s / 2.14 s after
the tested alternatives. Frequency placement is unchanged except for replacing
the upper collection with a cloud centre.

`RegionalObservationBasis` caches the exact quadratic regional Welch power of
the validated actual-WASM observation basis. Tests cover partial FFT-bin
integration and interference between columns. Each final candidate is rendered
again and verified against its saved patch; the cache is not a substitute DSP.
Reference-fixed regional spectra, band timing and fine ridge contrast remain
separate diagnostics, not a single perceptual acceptance certificate.

The standard seed is 1675; 1982 is a holdout. Holdout upper levels remain several
dB lower and are a limitation of these candidates, not something normalized
away. Movement-only and the cloud retain more fine ridge contrast than blur-only.
This needs listening judgement; no automated score establishes a winner.

Checks include strengths .3/.5/.76/1, four half-second-spaced hits, eight rapid
hits, native DSP tests, WASM/native API tests, and browser patch round-trips.
At 48 kHz, standard single hits and four repeats do not engage the limiter even
at Master 0 dB. Hard/rapid strikes can exceed digital full scale before the
browser limiter at Master 0; the normal −12 dB Master leaves headroom. No limiter
is part of the fitting render or instrument DSP.

Reproduction tools:

- `tools/fit_gong_texture_trials.py`: texture-family screen.
- `tools/polish_gong_texture_trials.py`: now shared gain/tilt only; see correction above.
- `tools/refine_gong_texture_dynamics.py`: shared dynamics/upper-cloud retuning.
- Local results: `build/gong-texture-trials/{blur,movement,cloud}-final`.

## References and limits of the analogy

[Woodhouse's plate modal-density derivation](https://euphonics.org/3-2-4-the-modal-density-of-a-vibrating-plate/)
gives approximately constant average modal density per Hz for thin plates. This
motivates the cloud layout; it does not establish the exact distribution of a
shaped gong's modes or our chosen finite oscillator count.

[Fitz and Haken's bandwidth-enhanced sinusoidal modeling](https://www.cerlsoundgroup.org/Loris/ICMC95/BandwidthOscillators.html)
describes retaining a sinusoidal representation while adding bandwidth through
noise modulation. It supports comparing broadened resonances with coherent
resonances plus moving sidebands, not assuming either must sound right.

[Touzé's cymbal synthesis examples](https://perso.ensta.fr/~touze/tapercymbals.html)
illustrate nonlinear modal synthesis and strike-dependent behaviour. Our cloud
remains a constructive approximation, not a reproduction of that plate model.
