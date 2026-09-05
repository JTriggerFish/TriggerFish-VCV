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

Velocity, implement, hardness, location, routing, sample alignment and source
gain are **not optimized**. Resonance uses one T60 and one damping tilt, not
independent mode decays. Thump and contact noise have their own finite source
envelopes; those are distinct from resonator energy loss.

Spatial centre/edge couplings stay fixed for this one-location fit. Output EQ
is not an optimization variable: the default start bypasses it. A manually
chosen radiation response can stay fixed, but a multiband-EQ start is rejected.
The renderer exposes 93 parameters (29 scalar controls and 16 × 4 modal values);
only the subset named in the stage bounds is optimized.

Saved starts must have exactly the current parameter keys and reference identity.
Old fixed-bank/EQ-assisted candidates are rejected, not silently converted.
The old EQ-assisted report has been withdrawn; no new explicit-modal real-sample
fit is claimed yet.

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

The existing server serves `/kick-review/`: one real reference versus one
current candidate, with shared scales, downloadable full fit JSON and protected
playback. Audio and plots are pre-limiter; browser monitoring alone uses the
3 ms safety limiter. Generated sample audio stays outside version control.
