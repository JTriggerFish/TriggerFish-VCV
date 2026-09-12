# Crash: ringing texture, shared tuning and decay

The current preset is **Crash — beating doubles, refined decay**; see
[Final reviewed result](#final-reviewed-result). The earlier experiments below
are retained as methodology history, not alternative current presets.

The first pass started from the saved **Crash — user low ring, gentle movement**
(`cf357fa8-925a-45d9-afd5-c98908a9e0ee`) at its stored edge strike, velocity 72.
The reference, onset alignment, event, body drive, source mix and output gains
are fixed. Gong is a useful texture starting point, **not** the crash's fitting
target. No DSP or workbench control changes are involved.

## What the starting comparison shows

Fixed-level spectra and band-envelope plots show a low body that fades too
quickly, while the 6–16 kHz band remains too strong through much of the first
two seconds. Making every decay shorter cannot fix both. The dominant low
peak is approximately 122.2 Hz versus 125.6 Hz in the reference, measured over
0.2–2 seconds with a Hann periodogram (about 0.56 Hz bin spacing).

The saved crash uses paired rings, phase blur 0.04, a negative blur tilt,
no shimmer, and slow pitch drift. The current gong uses zero blur and bounded
irregular shimmer. Comparing those mechanisms is preferable to adding another
noise layer. Shimmer does not constitute proof of a better texture on its own.

## What is fitted

`tools/refine_crash_shimmer.py` renders the actual workbench Wasm with
`SavedFitRenderer`, preserving the full saved patch and gesture.

1. Compare the original, EQ bypass alone, quiet paired rings, and nine shimmer
   trials: amounts 0.3/0.7/1.2 at 20/60/160 changes per second. These trials
   bypass EQ, turn blur off, use 512 available oscillators, reduce deliberate
   beat depth to 0.15 and use mostly independent shimmer (sharing 0.15).
   These are combined texture starting points, not a causal one-control audit.
2. From the best EQ-free texture start, optimize nine shared coordinates with
   bounded Powell searches: whole-series pitch scale, a smooth upper-frequency
   stretch, low/high T60 endpoints, diffusion rate, excitation tilt/centre, and
   two broad prominence shelves. No individual frequency, prominence,
   allocation or local noisiness value is searched independently.
3. Inspect separate metrics and fixed-colour plots, then use independent seeds
   and repeated strikes before choosing a workbench candidate. If needed, test
   sparse extra knots in the **shared** T60 curve only after this stage.

The follow-up decay search first tries energy-sensitivity values 0.7/1/1.4
and bloom-rate multipliers 0.6/1. It then screens **one** interior knot at
150/600/3000 Hz, never per-mode decay multipliers. Coarse trials scale the
interpolated T60 by 0.6/1.4/2 and the upper endpoint by 0.5/1; fine trials use
0.9/1/1.1 and 0.7/1 respectively. The unmodified candidate remains eligible,
so an extra knot is not compulsory.

The prominence shelves are baked into the visible bars; series edits are baked
into visible frequencies. They are fitting coordinates, not hidden runtime
coefficients. The patch retains its full parameter surface. Output EQ is not
used as a fitting variable; EQ-free candidates carry its visible bypass state.

## Ranking and limitations

`tools/crash_texture_diagnostics.py` combines existing, tested measurements:

| Term | Weight |
|---|---:|
| Reference-floor Mel MRSTFT, full six seconds | 1 |
| Same spectral comparison over the first 300 ms | 0.3 |
| Fixed-level spectral-envelope and rise error | 0.05 |
| Auditory-band envelope-modulation distance | 0.15 |
| Band-decay shape error, dB | 0.15 |

The Mel floor is fixed from the reference at 60 dB below each resolution's
peak; it does not normalize playback or change with the candidate. The decay
term removes a constant per band **only in that subproblem**, while the spectral
terms retain absolute level. Two seeds are used during search; separate seeds
73519 and 41273 are reserved for audit.

The first 200-evaluation shared-coordinate search used the unmasked spectral
envelope term (V2). Inspection then exposed an incentive to match the very
quiet, nearly flat upper tail. The current **V3** ranking and follow-up decay
search floor both sides of that sub-comparison, only after 1.7 seconds and
only for terminal bands whose 3–6 second variation is below 1.5 dB and whose
level is more than 35 dB below that band's peak. The floor is 10 dB above
the terminal mean. A louder synth tail is still penalized; decaying low
rings are not floored this way. The unmasked measurement is also retained.
This treats a quiet plateau as uncertain evidence, not proof of recording
noise, and changes neither WAV. Before/after scores must use the **same**
version; the original first-pass score is not comparable to a V3 score.

An additional diagnostic measures locally whitened spectral flatness in
1.5–3/3–6/6–12 kHz, with 1-second windows and a 50-Hz smoothing scale. It
distinguishes concentrated ridges from wash without treating broad EQ colour
as texture. It is descriptive, not another target to blindly minimize.

This is an engineering proposal ranking, not a validated perceptual distance.
In particular, neither Mel bins nor modulation statistics guarantee correct
isolated pitches, metallic character, or pleasant beating. We therefore report
the low pitch explicitly, inspect spectra and signed differences, and retain
the user's listening decision. Tests verify that an identical synthetic target
scores near zero, stretched noise decay is penalized, and turning everything
down cannot pass as a match.

## Reproduction and artifacts

Build with `dev.ps1 build-workbench`. Run the project Python environment with
`PYTHONPATH=python` and `EMSDK_NODE` set to the SDK Node executable:

```powershell
.\.venv\Scripts\python.exe tools/refine_crash_shimmer.py --budget 200
.\.venv\Scripts\python.exe tools/refine_crash_shimmer.py --fit build/crash-shimmer-refinement/candidate.fit.json --fine-decay --budget 0 --output build/crash-shimmer-fine-decay
.\.venv\Scripts\python.exe tools/audit_crash_shimmer.py build/crash-shimmer-refinement
```

These commands run a new search with the current V3 ranking; they do not
reproduce the earlier V2 optimizer trajectory. The saved candidate and source
snapshots are the exact inputs for repeatable before/after audits.

The scripts write only to `build/`: exact source snapshot, candidate snapshot,
parameter trace, per-seed metrics, WAVs, and Plotly spectra/envelopes and signed
STFT differences. PNG inspection uses `tools/capture_fit_plot.mjs` in a silent
disposable browser tab. No extra audition server or report page is created.
Publishing to the main workbench is a separate, reviewed file edit.

## Published candidate — 12 September 2026

The first published candidate was **Crash — clearer rings and shimmer**
(`4e29137c-32ab-41b8-bbec-46630d62fcd0`). This is a partial improvement, not a
listening-approved match. Gong, ride, kick and the DSP implementation are unchanged.

- Phase blur 0.04 → 0; shimmer amount 0 → 1.2 radians, speed 60 changes/s,
  sharing 0.15. Slow detuning remains 0.3 Hz at speed 1.
- Beat depth 0.30 → 0.15, density 0.85 → 1; packet spread stays 2.
- Shared series adjustment: first handle 121.75 → 124.30 Hz, progressively
  less shift towards the top; no independently optimized ridge frequencies.
- Bloom rate 6.60 → 5.19; excitation tilt −8.71 → −8.10 dB/oct,
  centre 2130 → 2284 Hz. Energy sensitivity and concentration are unchanged.
- T60 remains **two points**, approximately 29.58 s / 0.911 s at 40 Hz / 15 kHz.
  Neither coarse nor fine one-knot trials improved the overall decay trade-off.
- Broad prominence rebalancing adds about 1.1 dB to the low handles, reduces
  mid prominence, then applies a final smooth upper reduction reaching −2 dB
  at 6 kHz. These values are stored directly in the visible bars. Output EQ
  is bypassed; contact/body mix, strike energy and output gains are unchanged.

The final upper reduction was chosen from `--upper-balance` trials (−2/−4/−6 dB
reached at 4/6 kHz). Its two-seed V3 score is **2.1174**, versus **2.1093**
without it: a deliberate small score trade-off for less excess upper energy.
Larger reductions increasingly worsened the attack spectral comparison.
Independent-seed mean scores are effectively tied between those two variants;
the objective does not justify claiming one is perceptually superior.

Final audit: `build/crash-shimmer-published/` contains original/candidate
snapshots, reference/candidate/before WAVs, fixed-colour plots, signed STFT
differences, three velocity renders, five quarter-note strikes, and `audit.json`.

| Measurement | Original | Candidate | Reference |
|---|---:|---:|---:|
| Dominant 80–180 Hz peak, standard seed | 122.22 Hz | 125.00 Hz | 125.56 Hz |
| Mean 6–16 kHz level error, 0.3–1.5 s | +10.78 dB | +4.22 dB | 0 |
| Mean 100–300 Hz level error, 2–5 s | −8.02 dB | −5.72 dB | 0 |
| Mean 300–700 Hz level error, 2–5 s | +3.86 dB | +5.04 dB | 0 |
| V3 score, unseen seed 73519 | 2.2918 | 2.1676 | — |
| V3 score, unseen seed 41273 | 2.2401 | 2.1437 | — |

The last two band rows are important: the low body still fades too much and
the adjacent low-mid decay is now somewhat worse. Upper-mid ridges also remain
less distinct than the reference late in the sound. The candidate is not an
across-the-board match, and the scalar score must not conceal those limitations.

Five repeated hits peak at −5.63 dBFS before browser master gain, or −17.63 dBFS
with a −12 dB master. No output-level compensation was used. Verification also
checks that loading the published workbench preset produces the exact same PCM
as the reviewed render, rather than testing a different offline approximation.

Validation: 10 focused Python tests; 22 Wasm/workbench tests; two native API
tests and native/Wasm signature comparison; browser load/save checks for crash,
gong and ride. Audition remains the final texture judgement.

## Follow-up: user's beating-doublets starting point

Input: **crash beating test**, `8a53169a-671c-407b-907c-29552eca25c6`.
The attachment and `Documents/crash-beating-test.json` contain identical JSON.
The original is archived, unmodified, in
`build/crash-beating-refinement/user-snapshot.fit.json`.

This is a different texture choice, not just a small amount adjustment:
**Doublets** pairs the surrounding oscillators throughout each packet;
**PairedRing**, used previously, pairs its central ring. Packet spread rises
from 2 to 6, phase blur is 0.069696 with tilt −1.5, and shimmer is slower and
more shared. The user also adds three low-mid handles (27 active in total),
changes low-mid prominence and changes bloom dynamics. These choices are kept
fixed because the user prefers their texture; the previous numerical texture
score is not sufficient grounds to undo that judgement.

### Conditioning and allowed fitting coordinates

The saved performance is a **strong bow strike** (strength 0.954985, location
0.5), not the selected medium edge reference (strength 0.566929, location 1).
`source.fit.json` uses the reference cell's gesture and seed; the original
gesture is tested separately, not treated as a like-for-like reference fit.
The original strong strike peaks at −3.44 dBFS before master/limiter: the
preferred single-hit texture is not dependent on limiter clipping at unity
or the normal −12 dB master setting.

The `--locked-texture` stage fits **four** shared coordinates: low/high T60
endpoints, common pitch scale and smooth upper stretch. Bounds are 12–30 s,
0.25–3 s, 0.96–1.04 and −0.03–0.03 respectively. All prominence bars, texture,
bloom, gain, local allocation and contact parameters are fixed. No extra
decay points or per-mode decay fitting. Two training seeds and separate audit
seeds are retained; the loss is the same declared V3 ranking above.

```powershell
.\.venv\Scripts\python.exe tools/refine_crash_shimmer.py --fit build/crash-beating-refinement/source.fit.json --locked-texture --budget 80 --name 'Crash — beating doubles, refined decay' --output build/crash-beating-fit
.\.venv\Scripts\python.exe tools/audit_crash_shimmer.py build/crash-beating-fit --user-snapshot build/crash-beating-refinement/user-snapshot.fit.json
```

Search stages now live in `tools/crash_refinement_search.py`; rendering,
provenance and CLI orchestration stay in `tools/refine_crash_shimmer.py`.
A regression test checks that texture-locked edits can change only active
modal frequencies and the two shared decay endpoints.

### Final reviewed result

The main workbench now loads **Crash — beating doubles, refined decay**,
`8b5d14e9-93af-4f33-b950-496865caf9cc`. The original Documents snapshot is untouched.

After the four-coordinate search, one interior knot was justified: without
it, lengthening 1.5–6 kHz ringing also prolonged the very top. The same fine
one-knot screen described earlier selected **3 kHz / 4.443 s**, with a shorter
15 kHz endpoint. A final low-endpoint-only screen (multipliers 0.65/0.8/0.9)
then reduced low-mid overhang without changing the new 3 kHz setting.

Final shared curve: **40 Hz / 23.846 s; 3 kHz / 4.443 s; 15 kHz / 0.812 s**.
All other interior points remain inactive. This is a sparse final-stage
refinement, not permission to fit arbitrary per-mode decay multipliers.

The whole modal pattern is gently retuned (first handle 124.296 → 125.319 Hz),
with smooth upper stretching. All **27** active handles and their prominence,
local noisiness and allocation are retained. No texture, contact, bloom, EQ
or output-gain control is changed from the user's test. The reference gesture
is loaded by default; the user's stronger bow strike is separately validated.

| Independent-seed measurement | User test at reference gesture | Refined |
|---|---:|---:|
| Decay-shape error, seed 73519 | 8.55 dB | 4.05 dB |
| Decay-shape error, seed 41273 | 9.10 dB | 4.19 dB |
| V3 score, seed 73519 | 3.3495 | 2.4286 |
| V3 score, seed 41273 | 3.4317 | 2.4267 |

The low peak moves from 124.44 to 125.56 Hz in those audits; the reference
peak is 125.56 Hz at the analysis resolution. The original strong bow gesture
peaks at −3.44 dBFS, the refined one at −2.10 dBFS before master attenuation.
Five quarter-note medium strikes peak at −5.26 dBFS. There is no gain matching.

The 300–700 Hz decay remains too prominent, and detailed ridge placement is
not an exact match. Neither an improved aggregate score nor fixed texture
controls proves perceptual equivalence to the user's preferred sound;
audition is still required. The reviewed plots and before/after renders use
the original conditioned user snapshot, not the intermediate fitting stages.
Artifacts: `build/crash-beating-reviewed/` (including the strong-gesture WAVs).

### Tooling review safeguards

Each fitting run requires a fresh, empty `--output` directory. Existing results
are never silently overwritten or mixed into a new trace. Audit seeds exclude
the training seeds of both compared fits, including unusual user-selected seeds.
The final progress file is written even for a baseline-only run (`--budget 0`).

Texture-locked edits preserve even near-silent prominence bars exactly; broad
attenuation edits cannot accidentally raise a quiet bar. Optimizer starting
coordinates are bounded without rewriting the original candidate.

Repeated-hit audits now pass the complete saved fit to the renderer, including
its routing and strike gesture, instead of inheriting the current calibration.
These review fixes do not change the published crash preset's sound.

`dev.ps1 test-fitting-tools` includes the fitting regression tests. The saved-fit
sequence integration tests are opt-in: set `TF_TEST_WORKBENCH_BRIDGE=1` with the
local reference server and compiled Wasm available. They compare a one-hit
sequence against a normal saved-fit render, including edited routing and gesture.
Normal plugin builds and CI do not require that private sample corpus or server.

```powershell
.\.venv\Scripts\python.exe tools/refine_crash_shimmer.py --fit build/crash-beating-fit/candidate.fit.json --fine-decay --budget 0 --name 'Crash — beating doubles, refined decay' --output build/crash-beating-fine
.\.venv\Scripts\python.exe tools/refine_crash_shimmer.py --fit build/crash-beating-fine/candidate.fit.json --low-decay --budget 0 --name 'Crash — beating doubles, refined decay' --output build/crash-beating-final
```
