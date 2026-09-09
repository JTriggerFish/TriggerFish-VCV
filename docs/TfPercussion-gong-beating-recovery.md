# Gong: recovering a gentler ring

The September 2026 workbench preset **Gong — gentler beating** is a two-control
re-edit of **Gong — structured beating**, not an engine or controls rollback.

- Doublet separation: **6 → 3 Hz**.
- Noisiness slope: **0.21266 → 0.4 per octave**.

Modal frequencies and prominence, contact, bloom, T60, observation and gain
remain unchanged. The steeper slope cleans the lower packets without replacing
the upper sizzling texture with a blurred noise wash. Phase blur is unchanged.
This preset still uses *Beating doublets*: paired-ring depth and rate-tilt
controls are correctly inactive for that layout.

## How the choice was tested

`tools/recover_gong_beating.py` screens 17 texture-only variants at three seeds:
slower doublets, modest phase blur, steeper noisiness slopes, and the current
paired-ring layout at several depths/rates. Reference-relative low-band envelope
modulation, absolute reference-floor Mel distance and upper-band texture are
measured separately. None normalizes playback gain.

Simply halving separation improves the modulation-frequency match, but barely
reduces low-band pulse depth. Switching to paired rings reduces depth much more,
but changes the attack and upper texture more than wanted. It was **not** chosen
just because it suppressed beating. The selected slope-plus-separation edit
reduces standard-strike 90–180 Hz relative envelope RMS from 0.354 to 0.279
(reference 0.223). These values remove an exponential trend for measurement
only; they are not audio envelope controls or proof of listening equivalence.

`tools/check_gong_recovery.py` checkpoints the selected parameters, source and
reference audio, and checks exact-render provenance. `review_metal_refit.py`
then checks four seeds and quarter-note / rapid-hard repeated hits.

| Four-seed mean metric (lower is better) | Before | Candidate |
|---|---:|---:|
| Mel distance | 1.290 | 1.233 |
| Attack Mel distance | 1.253 | 1.234 |
| Band decay-shape error, dB | 4.391 | 4.309 |
| Upper texture distance | 0.343 | 0.328 |

The initial 450 ms changes by 3.37 dB RMS in active spectral cells versus the
previous synth. This is not an inaudible repair, even though the reference
attack score improves slightly. Fixed-scale band-envelope plots retain the
bloom pattern; remaining errors include excess early 3–6 kHz and too much
late energy above 6 kHz. Listening approval is still needed.

Artifacts are in `build/gong-beating-recovery/{before,candidate}`, with the
screen in `screen.json`. The main workbench preset is the audition surface;
no additional server or report page is created. The old preset's complete
parameter set and WAV are archived in `before`.
