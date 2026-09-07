# Kick body-envelope experiment

**Status: the 20 ms-hold candidate below was rejected in listening.** Its broad
envelope score did not establish that the ringing was fixed. Do not reuse it as
an approved calibration or describe its lower aggregate loss as sound quality.

## Scope

The standard workbench Kick, against its existing medium-velocity oak-kick
reference. Source decoding, onset, reference gain, strike strength and implement
are unchanged. Model level remains −12 dB; no automatic matching. The original
hold trials bypassed EQ; the final audition below uses visible radiation EQ.

## Fixed beater

Kick modes now expose frequency and prominence only. Centre/edge coefficients
and the kick location control were removed, including from serialized patches.
The generic trigger ABI retains location for other instruments, but ignores it
for kick. Shared membrane spatial controls remain available to toms and snares.

The starting preset's old spatial weights were folded into its visible modal
prominences and resonance level before refinement. This preserves individual
transfer magnitudes, not the relative phase of an old negative projection.
No hidden replacement phase or spatial coefficient was introduced.

## Envelope experiment

The thump originally rose in 0.4 ms and immediately decayed exponentially.
`thump_hold_seconds` inserts a constant-amplitude segment after that rise.
Zero preserves the previous trajectory. Range: 0–80 ms; pitch continues moving
during the hold. T60 measures the subsequent decay, excluding the hold.
The envelope is amplitude-continuous, but its slope has corners: this is a
plateau experiment, not yet a rounded-shoulder envelope. It does not clip audio,
gate the modal tail, change contact, or inject energy into the resonators.

## Fit procedure and limitations

`tools/kick_body_refinement.py` renders the actual workbench Wasm, with training
seeds 1449 and 1450. It tests shared modal damping, thump duration/level and modal
prominences with bounded least squares and logged central finite differences.
Frequencies, pitch trajectories, contact parameters and output gain stay fixed.
The existing multi-resolution spectral and band-envelope objective is used;
this experiment does not introduce or claim a perceptual listening score.

Trials used 0, 20 and 40 ms holds. The 20 ms trial was then allowed two existing
upper-bass handles near 167 and 199 Hz. No per-mode decay controls were added.
The chosen intermediate preset has a 20 ms hold, thump T60 about 289 ms,
modal T60 at 100 Hz about 456 ms and damping slope about 0.35 instead of 1.
The two prominent modes near 574 and 808 Hz were reduced.

`tools/plot_kick_body_trials.py` writes offline Plotly band-envelope comparisons;
`tools/capture_local_plot.mjs` can capture them without touching the user's tab.
Plots were visually inspected. The hold gives only a small incremental benefit.
Late 90–350 Hz energy remains insufficient and the initial 350–1400 Hz energy
is too weak. This is an audition candidate, not a completed calibration.
Audio audition stays in the main workbench; no additional server is needed.

Next: compare rounded shoulders against the plateau, and address those remaining
band/time discrepancies without per-mode damping or silently changing contact.

## Follow-up: direct ridge diagnosis

`tools/inspect_kick_ridges.py` plots identical-resolution reference/model STFTs,
region spectra and isolated thump/body contributions. The rejected preset had
about 18 dB excess near 813 Hz at 80–160 ms, plus excess near 576 Hz and a new
late 167 Hz ridge. The first narrower-band refinement suppressed some of this
but made the 87 Hz mode too strong, so it was not published.

`RidgeBalanceLoss` now supplements the existing envelope/spectral objective with
4096-sample spectra at 100–4000 Hz and 8192-sample spectra at 20–250 Hz. Reference
floors are fixed. Quiet bins retain nonzero weights so moving a false ridge into
a previously quiet bin does not make it disappear from the objective. Positive
excess beyond 3 weighted residual units contributes an additional residual:
3 dB in active bins, about 9.5 dB in quiet bins (weight 0.1). Regions are 0–40, 40–80,
80–160, 160–260 and 260–400 ms. Long windows smear time: they are used for ridge
discrimination, not as evidence of an exact acoustic onset or decay knee.

Synthetic tests check zero loss for identical input, increasing penalty for a
false ridge, and sensitivity to an overlong mode. These are measurement checks,
not proof that a synthesized kick is realistic. Every candidate still requires
direct reference comparison; no automatic publication runs in these scripts.

The follow-up freezes one relative modal level to remove the common modal-level
scale ambiguity during optimization, while leaving audible bank gain editable.
Thump pitch trajectories and contact width are also eligible: incorrect source
envelopes must not be compensated solely by longer or louder resonant modes.

## Current standard Kick audition — 2026-09-07

The rejected preset did not merely need a longer envelope. Its long contact
noise also drove the membrane, continually exciting its poles. Fitting that
shared signal made direct noise, modal ringing and decay difficult to separate.
The selected candidate uses two explicit, serialized Contact selectors:

- **Body drive: Pulse only.** A roughly 2.7 ms contact pulse excites the body;
  subsequent source noise no longer drives it.
- **Observation: Noise only.** Direct noise remains audible independently,
  without directly observing the pulse/chirp/grain components.

The two offending handles at 574/808 Hz are disabled at the −72 dB off floor.
Six handles remain near 25, 57, 87, 118, 138 and 164 Hz. Shared modal T60 is
397 ms at 100 Hz, with slope 0.401; there are no per-mode decay multipliers.
The thump holds about 26 ms, then has T60 164 ms. Its fitted curvature is
effectively zero: the rounded-shoulder trial did not establish an advantage.
Radiation EQ uses a 5 Hz high-pass, approximately 833 Hz low-pass and a broad
+7.1 dB colour peak around 964 Hz. No multiband correction is fitted.

Refinement uses `refine_kick_ridges.py`, then `refine_kick_low_body.py`; bounds,
fixed controls, renderer hashes and central-difference steps are recorded in
each local `search.json`. `audit_kick_candidates.py` is a separate band-power
check, not the optimizer. For the primary seed, 100–200 Hz errors at 80–160 /
160–260 ms improve from −5.8 / −8.2 dB to approximately +1.2 / −1.1 dB.
The former narrow upper ridges are removed, but broad 540–650 Hz energy remains
too weak. Do not describe this as an exact spectral or perceptual match.

`verify_kick_audition.py` rerenders with the current engine, verifies exact fit
reload and reference identity, and tests repeated hits at 500/125 ms intervals.
Standard-strength peaks are about 0.44; full-strength repeated hits peak around
0.84 before the browser limiter. Gain-linearity checks exposed cancellation in
the float 5 Hz high-pass; double coefficient/state arithmetic fixes that without
relaxing the check or clipping audio. Current-render provenance is recorded
separately from the earlier fitting renderer.

The single main-workbench preset is now this **unapproved audition candidate**.
The full-match gate remains false and listening approval remains false. No
automatic publication gate was relaxed; this is an explicitly labelled manual
audition update. Source/reference gains and performance inputs are unchanged.

## Verification and audit safeguards

Current-engine rerendering checks the archived reference hash, sample rate,
monitoring gain, onset and performance event **before** rendering or writing.
Recipe and parameter descriptors must match too; only the renderer binary may
change. The new record retains the original metadata, parameters and source-file
hash. Old measurements remain history, not scores from the new engine.

The audit requires explicit archived `search.json` files, never today's preset:

```powershell
$env:EMSDK_NODE = 'D:\dev\emsdk\node\24.19.0_64bit\node.exe'
.\.venv\Scripts\python.exe tools/audit_kick_candidates.py `
  --baseline build/frozen-baseline/search.json `
  --candidate build/kick-final-audition/search.json `
  --output build/kick-comparison-new.json
```

Choose an actual frozen baseline with the current parameter schema; the example
path is not generated automatically. Reports retain both parameter vectors,
source hashes and current-renderer identity and refuse to overwrite an existing
report. These are current-engine comparisons of archived parameters, not replays
of an old binary. Held-out seeds exclude the primary seed and every recorded
training seed. Missing training-seed provenance is an error, not an assumed split.
