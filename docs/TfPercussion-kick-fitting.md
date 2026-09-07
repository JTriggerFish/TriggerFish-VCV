# Fitting the unified kick

This procedure fits **one centre beater hit, velocity 64, take 1** from the local
kick collection, using the actual `drum.kick.v1` C++/Wasm renderer. It does not
fit a surrogate synthesizer. Reference file identity, SHA-256, onset, source
gain, sample rate and performance inputs are saved with every candidate.

The signal path and controls are described in
[kick architecture](TfPercussion-kick-architecture.md).

## What is fitted

The variables are the visible Contact, Thump, Resonance, strike/tension and
explicit active-mode frequencies and prominence. Their exact numerical search bounds are in
`tools/kick_fit_stages.py`; the renderer rejects values outside UI bounds.
Master attenuation is fixed during joint search, because fitting it alongside
all three source gains creates an exact gain nullspace.
The strongest starting modal bar is also held fixed during each local pass;
the other levels are relative to it, avoiding a redundant common scale.

Velocity, implement, hardness, routing, sample alignment and source
gain are **not optimized**. Resonance uses one T60 and one damping tilt, not
independent mode decays. Thump and contact noise have their own finite source
envelopes; those are distinct from resonator energy loss.

Kick no longer has location or spatial centre/edge coefficients. The original
procedure bypasses output EQ. Current diagnostic trials additionally compare
a fitted low-pass cutoff and one broad radiation colour peak, with a fixed 5 Hz high-pass;
they do not use a fitted multiband correction curve.
The renderer exposes 65 parameters (33 scalar controls and 16 × 2 modal values);
only the subset named in the stage bounds is optimized.

Saved starts must have exactly the current parameter keys and reference identity.
Old fixed-bank/EQ-assisted candidates are rejected, not silently converted.
The old EQ-assisted report has been withdrawn. The current explicit-modal
candidate below is an experimental fit, not a listening-approved calibration.

## Loss

No signal is independently normalized. Four equally weighted representations
compare absolute amplitude:

| Representation | Frequencies | Purpose |
| --- | --- | --- |
| Hann STFT, 512 samples | 250 Hz–16 kHz | Attack/noise spectrum |
| Hann STFT, 2048 samples | 20 Hz–3 kHz | Body structure |
| Hann STFT, 8192 samples | DC–500 Hz | Bass pitch and unwanted subsonics |
| Power envelope, Gaussian 12 ms, sampled every 2 ms | Full signal | Attack/decay timing without following individual bass cycles |

STFT hop is 256 samples. Each representation uses a fixed reference-only
peak-minus-70-dB floor and reference-only salience weights. Weighted squared
dB residuals are summed. Time regions 0–30, 30–100, 100–250, 250–600 and
600–1200 ms have weights .30, .30, .25, .10 and .05.

The long bass window spans adjacent regions; it cannot resolve a 30 ms event
independently. Short-window spectra and the time-domain envelope supply timing.
The late region penalizes excess tail, not a T60 extrapolated from padded
digital silence. The aggregate error is an engineering discrepancy, **not a
perceptual quality score** or listening acceptance.

## Search procedure

1. Preserve the original reference and baseline render.
2. Load an explicitly identified starting parameter vector.
3. Without a resume file, compare explicit reference-spectrum modal proposals
   at several counts/falloffs before local fitting. Log their frequencies and
   scores. Editor formulas can also supply starts; none remain runtime rules.
4. Refine thump/shared damping/tension, active modal frequencies and prominence,
   then contact. Modal search ranges are local to their initial frequencies.
5. Jointly refine all active controls over two noise seeds. Residuals are
   concatenated; stochastic audio is never averaged.
6. Validate three held-out seeds, actual fit-file reload, branch superposition,
   resonance-level linearity and repeated hits.
7. Inspect the shared-scale spectrograms, attack/bass spectra and waveforms;
   retain explicit mismatches even when the aggregate improves.

Bounded least squares uses parameters scaled to their search ranges.
A 2% finite-difference probe records influence and freezes directions below
0.05 dB residual change. Central-difference Jacobians use 0.5% steps normally,
0.1% for fine refinement. Bounds, steps, renderer hash, objective specification,
actual render counts and parameter values are logged per stage.
A step must improve on the actual preceding patch, not a clamped substitute.

## Tests and reproduction

Synthetic recovery through the exact renderer recovers known 49 Hz / 0.7 s
resonance, and known 1.2 / 0.18 s contact-noise level/T60. Python tests cover
gain, pitch, damping, subsonic contamination and phase-insensitive RMS behavior.
These demonstrate known-parameter recovery, not real-instrument accuracy.
Shared safeguards and implementation entry points are documented in
[reusable fitting lessons](TfPercussion-fitting-lessons.md).

```powershell
$env:TF_KICK_SELF_TEST = '1'
.\dev.ps1 fit-kick-start
Remove-Item Env:TF_KICK_SELF_TEST

# Continue a saved candidate, including a finer joint pass.
$env:TF_KICK_FIT_RESUME = 'build/workbench-wasm/site/kick-review/search.json'
$env:TF_KICK_FINE_ONLY = '1'
.\dev.ps1 fit-kick-start
```

With no resume path, fitting starts from the current reference-target preset.
A resume path must contain the current explicit-mode parameter set.
`TF_KICK_AUDIT_ONLY=1` rebuilds the report without fitting.

The selected fit is published directly to the **main workbench**. Choose
**Acoustic kick — medium centre** in the reference-target selector; this loads
both the current audition candidate and its reference. First selecting the Kick recipe
uses the same calibration. Existing edits are preserved when switching recipes.
The versioned source is `workbench/web/kick_calibration.fit.json`; the build
copies that same file into the served site. No manually maintained second vector
or additional report-page audition is required. Audio and plots are pre-limiter;
browser monitoring alone uses the 3 ms safety limiter. Samples stay local.

Publication verifies the saved parameter render and fit-file reload before
replacing the source/served JSON. `TF_KICK_VERIFY_PRESET=1` then checks that the
actual reference-target preset renders sample-identically. Browser tests export
the current sound and compare all current parameters for both entry points. Refresh
an already-open page to load updated modules; never overwrite unsaved browser
state remotely. Debug reports remain optional (`TF_KICK_AUDIT_ONLY=1`).

## 2026-09-05 experiment

The revised factory start scored 9.187 on the fixed two-seed objective. Local
refinement reached 6.137. Raw reference-only modal layouts initially scored
worse, but independently refitting layouts that retained two existing handles
and added reference peaks reached 5.451. Joint and fine refinement reached
**5.032**. The primary-seed score is 5.018; held-out seeds score approximately
5.40, 5.13 and 5.28. These are weighted dB discrepancies, not perceptual scores.

Six handles remain active at approximately 30.38, 54.61, 87.26, 141.55, 586.08
and 816.40 Hz. Thump settles near 26.96 Hz. Shared resonance T60 is 0.307 s at
100 Hz with slope 1; no per-mode decay parameters exist. Model level remains
-12 dB, source monitoring gain +2 dB, and output EQ is bypassed.

I inspected the shared-scale plots. Upper ringing is now represented, and the
early RMS envelope is substantially closer, but the 0–30 ms attack-spectrum
error remains 9.54 dB and the 100–250 ms bass-spectrum error is 6.64 dB. The
low tail also extends too long below about -60 dBFS. The optimizer effectively
mutes direct contact and puts contact colour and damping slope on their bounds.
Those are diagnostic clues for the next contact/trajectory investigation,
not evidence that an EQ stage or extra decay knots should be added.

The exact saved UI fit reloads sample-identically. Held-out-seed, source-sum,
resonance-gain and repeated-hit checks completed; report playback/downloads
were checked in a temporary browser tab. At this stage the workbench preset was
not replaced, which meant the user was hearing the older sound. This workflow
error is corrected below.

## Workbench correction and band-balance refinement

The user's snare-like attack / weak thump feedback concerned the **workbench
preset**, not the separately rendered candidate. Short-noise constrained trials
were started against that wrong baseline and not selected. In particular,
limiting noise T60 to 70 ms improved some band envelopes but worsened the overall
spectral match. Do not present those trials as an accepted improvement or infer
a required exciter redesign from them.

The selected refinement warm-starts the existing six-mode candidate, adds
`DrumBalanceLoss` (see reusable lessons), and retains the ordinary parameter
bounds. It improves the two-seed combined objective from 4.12 to 3.97. Those
numbers are **not comparable** to the earlier spectral-only score of 5.03.
Primary-seed spectral discrepancy changes from 5.02 to 5.17 while band-envelope
discrepancy improves from 2.59 to 2.35; this is a measured tradeoff, not listening
approval. No EQ, per-mode damping or additional DSP was introduced.

Compared with the older workbench preset, direct contact observation is nearly
muted and clean thump gain rises from 1.77 to 3.21 (about 5.2 dB). Its settled
pitch stays near 27 Hz, consistent with the reference's low tail. The body uses
six active modes and shared T60 of 0.309 s at 100 Hz, slope 1. Source noise still
excites the body: its base T60 is 0.221 s, not a newly shortened excitation.
Model gain remains -12 dB and reference gain +2 dB.

```powershell
$env:TF_KICK_FIT_RESUME = 'build/workbench-wasm/site/kick-review/search.json'
$env:TF_KICK_BALANCE_REFINEMENT = '1'
$env:TF_KICK_JOINT_ONLY = '1'
.\dev.ps1 fit-kick-start
# Publishes the verified selected fit to the workbench automatically.
```

The current full editable [workbench calibration](../workbench/web/kick_calibration.fit.json)
contains the reference identity and fixed performance inputs, but no sample audio.

The subsequent listening review rejected this fit. See the
[source-isolation diagnosis](TfPercussion-kick-diagnosis.md): this candidate
fails the new independent shape/decay checks. The diagnostic trials do not
replace it merely for improving an aggregate score; publication now requires
those necessary checks in addition to reproducibility.

The subsequent [matched perceptual-loss experiment](TfPercussion-perceptual-loss-experiment.md)
tests established mel MR-STFT and JTFS objectives on the unchanged C++ voice.
All eight bounded real-reference trials are rejected; known contact-noise
recovery succeeds for each objective. No new workbench preset is published.

## Current audition (2026-09-07)

The standard workbench Kick now uses the pulse-driven-body / noise-only-direct
candidate documented in [body-envelope experiments](TfPercussion-kick-body-envelope.md).
It is explicitly **audition-unapproved**, not a full-match calibration. Its exact
saved-fit reload and repeated-hit checks pass; the independent full-match gate
still fails. This is a manually staged listening experiment, not a weakening of
the automatic calibration publication gate. Older results above are history,
not the parameters currently served.
