# Shared crash, ride and gong fitting

This pass starts from the saved workbench candidates and changes the fitting
method, not the instrument DSP. The kick, snare and hi-hat are left alone. Each
metallic target keeps its own fixed reference, source gain, onset and performance
event. Scores are compared within a target; no weights are adjusted per instrument.

## What was wrong with the previous search

Smoothed log-band energy can improve while narrow ringing becomes diffuse. It
also spreads the importance of a short, loud attack over many quieter frames.
The previous six-second search did not fit the ride's longer decay. Finally, a
small dB perturbation of an already quiet observation bar could mark it inactive,
even though a larger amplitude increase would be useful.

There was also a concrete reference-alignment error in an experimental ride
script. Its early timing scores were discarded; it never replaced the workbench
preset. `aligned_reference` now centralizes onset slicing and zero padding, and
publication checks compare the actual saved reference samples against that result.
A finer subsequent audit found the catalog marker itself was early: source onset
70.104 ms included about 14 ms of stationary noise before the strike. The
authoritative medium-bow marker is now **84.125 ms** (sample 4038 at 48 kHz),
selected from the reference's sharp rise, not candidate correlation. All ride
timing comparisons were rerun after this correction. Crash at 50 ms and gong at
0 ms had no analogous quiet lead-in and were left unchanged. `audit_onset`
reports millisecond RMS and threshold crossings without automatically shifting
anything; a gradual mallet attack must not be mistaken for noise pre-roll.

## One objective across all three

`MetallicBalanceLoss` exposes four separate error components:

| Squared share | Measurement | Failure it makes visible |
| --- | --- | --- |
| 35% | Seven band-power envelopes | Missing attack energy, premature decay, wrong bloom timing |
| 30% | Linear regional spectral magnitude | Strong real components missing despite a reasonable log-average |
| 20% | Local log-spectral contrast | Blurred ridges, filled valleys, overly tonal wash |
| 15% | Short-window attack spectrum | Wrong contact brightness and strength |

Band filters are causal second-order Butterworth filters. Power smoothing is a
**centered** 12 ms Gaussian; it is not a causal onset detector. Identical filtering
and smoothing are applied to reference and candidate. The first region splits
0–30 ms from 30–120 ms; other regions are 120–500 ms, 0.5–1.5 s, 1.5–3 s and
3 s to the chosen end. The fitting duration is recorded explicitly.

Resolved spectra use 16384-sample Hann windows, hop 4096, and only half-bin
smoothing. Contrast subtracts an eight-bin Gaussian-smoothed log spectrum from
the log spectrum. Attack spectra use 512-sample windows and hop 128. Floors and
salience weights depend on the reference only. Waveforms are never normalized.

For regional magnitudes $m$ and reference $r$, the linear contribution is

$$
20\frac{m-r}{\max(\lVert r\rVert_2,\epsilon_r)}.
$$

The denominator scales the residual, not the candidate audio. Consequently this
objective has **mixed units**, not simply dB, and is not a validated perceptual
distance. Separate components, plots, held-out seeds and listening remain needed.

An additional shared trial uses equal-ERB-rate contrast weighting, rather than
equal counts of linear FFT bins. This makes low/mid ridges count appropriately
without creating instrument-specific weights. The choice is recorded as
`contrast_weighting`; `linear` remains available for direct ablation. Set
`TF_FIT_CONTRAST=erb` to reproduce that trial.

The correctly aligned ride then exposed another blind spot: its first 10 ms was
8–12 dB too quiet even with a reasonable 512-window attack score. The optional
`fast_attack=True` variant splits the existing 15% attack share equally between
the STFT and raw power in **0–1, 1–3, 3–10, 10–30 and 30–100 ms** bins. It does
not compare waveform phase or maximize a random sample peak. A known-answer test
removes early energy while preserving the total first-100-ms energy. This must
still be detected. `TF_FIT_FAST_ATTACK=1` enables the same variant for all targets.
Crash supplies a useful opposite case: its candidate had excess early energy but
too little later attack energy, so simply turning every contact up is not valid.

This remains a trade-off objective, not an acceptance certificate. The ride
optimizer could still sacrifice roughly 6 dB in its first millisecond for other
terms. The subsequent bounded experiment therefore keeps the contact fixed and
constrains each raw-power bin to within 3.5 dB of the reference while optimizing
observation amplitudes on both training seeds. Feasibility is checked separately
from the average score. This is a fitting constraint, not a playback envelope or
gain correction.
The [attack-gate reference](TfPercussion-attack-energy-gate.md) documents the
reusable solver, analytic constraints, contact-presentation refinement and tested
commands.

The first constrained ride trial exposed a second loophole: it met the raw-power
gate partly through a false bass impulse (+21 dB in the first 30 ms, 40–250 Hz).
That trial was **not published**. The existing contact high-pass removes it;
subsequent contact-presentation trials must pass both early-power and independent
band-energy checks. The raw-power gate alone is not a timbral acceptance rule.

There are two important measurement limits. The long resolved windows overlap
region boundaries: the first regional spectrum is **not** an isolated attack
measurement. Use the short-window spectrum and raw millisecond bins for that.
Also, spectral terms currently stop at 16 kHz, while raw-power bins include all
frequencies. The ride reference has appreciable energy above 16 kHz; this remains
an explicit blind spot in the spectral comparison, not evidence of a full-band
match.

The combination of linear and log spectral terms follows established practice in
[DDSP](https://arxiv.org/abs/2001.04643) and
[auraloss](https://github.com/csteinmetz1/auraloss/blob/main/auraloss/freq.py).
The particular contrast term, region splits and weights here are engineering
choices, not claims from those sources. The need to distinguish tonal and noisy
energy is also discussed in
[DDSP-SFX](https://www.dafx.de/paper-archive/2024/papers/DAFx24_paper_51.pdf).

## Efficient observation fitting without another synthesizer

With state evolution and active mode membership fixed, output is an affine
combination of the painted observation amplitudes. The basis comes from actual
C++ renders and is checked with mixed-gain renders before use.

The bounded finite-difference fallback now searches **linear amplitude** while
JSON and UI values remain dB. All active observation bars are considered. Solver
coordinates are offset to avoid a tiny trust radius when the starting bank is
quiet. Other nonlinear controls retain their ordinary sensitivity checks.

The optional `polish_observation_autograd` path differentiates only the linear
combination and its analysis using PyTorch. It does **not** implement a second
percussion model. L-BFGS-B jointly updates amplitudes; each proposed final fit is
judged and saved using the actual C++ engine. Runtime checks compare the Torch
objective with canonical NumPy measurements and check an objective gradient by
finite differences before optimizing.

The analysis uses FFT convolution of fixed band-filter impulse responses. Their
one-second truncation is checked through measurement parity, not silently assumed
equivalent. On the saved crash, ride and gong waveforms, the squared objectives
agreed with the canonical implementation to relative errors below $10^{-13}$.
Tests also cover identity, literal gain recovery, missing attack, late-tail
truncation, shifted bloom, and analytic-versus-numerical gradients.

## Search and review policy

Use measured spectral peaks as layout proposals, not as proven physical modes or
ready-made amplitude coefficients. Compare layouts after refitting their levels.
Separate satellite density from phase broadening: reducing both together can
remove necessary wash while fixing an unrelated ringing problem.

Alternate energy/damping blocks with observation updates. Input excitation gain
and output observation gain are not interchangeable when energy-dependent
transport is active. Both are existing visible controls; no hidden correction is
introduced. Two T60 endpoints remain the default, with no per-mode damping fits.
Any eventual extra knot must be justified by measured decay after other errors
are addressed, not used to conceal an incorrect attack or energy trajectory.

For this ride alone, a single additional knot near 700 Hz was retained after
measuring stable ridge decays: approximately 27–30 s around 345/685 Hz, versus
7.2 s at 1319 Hz, 4.5 s at 2442 Hz and 1.7 s at 7775 Hz. It improved the longer
decay comparison across validation seeds. A two-endpoint ERB-interpolated curve
could not describe that bend. This is not a new default, and no per-mode damping
was fitted; neighbouring real ridges still have different decays the shared curve
cannot reproduce exactly.

Seeds repeatedly examined during search are **validation seeds**, not an untouched
test set. Once a candidate is frozen, `tools/verify_metallic_candidate.py` tests
three previously unused seeds, checks exact snapshot/WAV reproduction and freezes
the previous parameter vector inside its audit. Rerunning after publication
reuses that comparator, not the newly published preset. Do not optimize those
fresh audit seeds and then continue calling them held-out tests.
Use `--seed-offsets` to allocate a new final audit after another fitting round;
merely creating a new output directory does not make previously inspected seeds
unseen again.

Nonlinear state controls can produce rough objectives, including different
stochastic phase paths. A tiny finite-difference step is not automatically a
useful gradient. Short population searches are appropriate for those few controls;
the validated affine observation step remains a separate smooth subproblem.

## Reproduction

```powershell
$env:TF_FIT_TARGETS = 'crash,ride,gong'
$env:TF_FIT_WORKERS = '3'
$env:TF_FIT_OBJECTIVE = 'balance-v2'
$env:TF_FIT_CONTRAST = 'erb'
$env:TF_FIT_FAST_ATTACK = '1'
$env:TF_FIT_SECONDS = '12'
$env:TF_FIT_OBSERVATION = 'autograd' # optional; requires development PyTorch
$env:TF_FIT_OUTPUT = 'build/shared-metallic-experiment'
./dev.ps1 fit-instruments
```

The launcher builds once; each target then has an independent renderer/output.
New metallic experiments default to the shared ERB + fast-attack objective.
`TF_FIT_RESUME=1` preserves the recorded duration, objective options and training
seeds when not overridden. A changed objective or duration is rejected: use a
new experiment directory. `polish-instruments` likewise reconstructs the saved
objective rather than silently reverting to the historical trajectory loss.
BLAS defaults to one thread per target to avoid oversubscription. Source presets
are not updated automatically. The ordinary Rack build requires none of these
Python, Torch, browser or Wasm tools. Auditioning stays in the main workbench.

## Selected workbench candidates — 2026-09-07

All three use the unchanged C++ engine and fixed performance inputs/reference
gains. These are improved starting points for listening, **not accepted matches
at the kick's quality**. Common-scale plots, full tails and repeated strikes were
reviewed in addition to the objective. The kick, snare and hi-hat fits are unchanged
from the preceding pass.

| Target | Standard shared score, previous → selected | Changes retained | Remaining mismatch |
| --- | --- | --- | --- |
| Crash | 10.57 → 9.47, 6 s | Four measured midrange centres, observation amplitudes, modestly narrower phase broadening, two-endpoint damping | Too diffuse; 1–4 kHz energy about 5 dB low through part of the early body |
| Ride | 13.04 → 9.33, 12 s | Correct reference onset, stronger and high-passed contact, measured ridge layout, one 700 Hz damping knot, constrained observation fit | Fine ridge/valley detail, upper attack spectrum and late low-mid decay still differ |
| Gong | 12.49 → 8.89, 8 s | Shorter mallet contact, lower excitation centre/tilt, input energy versus observation gain, energy travel, two-endpoint damping | First millisecond too spiky; 140–250 Hz gap; low tail too short |

Scores above use ERB contrast plus the fast-attack term. They are mixed-unit
engineering measurements, not percentages or perceptual grades. Different
durations also prevent comparing scores across rows.

Three fresh, previously unused seeds were tested after each candidate was frozen.
Crash and ride improve their overall score on all three. Gong improves on two;
the third changes 9.10 → 9.73. Its average gain does not erase that regression or
the weak low-mid attack on other seeds. No fresh-seed retuning was performed.
The ride's early-power and first-30-ms band constraints hold on the two fitting
seeds. Its final fresh audit uses offsets 4001/5003/6007; some intervals still
differ by 4.3–4.7 dB, outside the 3.5 dB fitting tolerance. Previously,
first-millisecond deficits were around 15–19 dB. No playback normalization hides
those differences.

The gong supplied a useful control-surface finding. Its old width setting made
the mallet half-sine about 13.9 ms long, putting a spectral zero near its 124 Hz
ring. Shortening it to about 5.5 ms restores that drive. However, the same width
macro also scales noise attack, so slowing the noise rise also lengthens the
coherent contact again. Noise-mixture trials 0/.03/.05 did not resolve that
conflict, and zero would silence the brush path. The nonzero mixture was retained.
This is a documented limitation to revisit explicitly, not hidden compensation.

Private audit outputs (not committed audio) are under:

- `build/instrument-refits-v3/crash/erb-refined`
- `build/instrument-refits-v3/ride/corrected-onset/contact-shape-final`
- `build/instrument-refits-v3/gong/short-erb-frequency`

Each includes exact candidate JSON/WAV provenance, shared measurements, a
`fresh-seed-audit.json` and full-response/restrike checks. Source JSON records the
selected fit and its remaining limitations. The main workbench's **Reference
targets** entries load these exact patches; separate comparison pages are not
published.

Verification for this publication includes 119 fitting-tool tests, native/Wasm
parity and an isolated-browser round trip of all five non-kick Reference targets.
Each browser-saved patch matches the source parameters, event, reference hash,
gain and onset. Repeated maximum-strength hits may exceed unity before the
browser limiter; neither normalization nor a new limiter was added to the core.

The pre-commit review added checks for duplicate or altered snapshot parameters,
reference identity and performance metadata, even when those changes leave the
render unchanged. Attack plots now crop the selected audio region before the
STFT, so later energy cannot leak into a panel labelled "Initial 120 ms". This
plotting correction does not change the recorded fitting objective or the audio.
Torch-dependent optimizer tests skip when the optional package is absent;
ordinary analysis tests continue to run in development CI.
