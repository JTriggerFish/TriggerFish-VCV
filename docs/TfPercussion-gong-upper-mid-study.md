# Gong upper-mid texture: controlled refinement

## Listening in the workbench

This study used the user's unchanged **gong test**, snapshot
`04af4b2a-10a6-4a72-a98d-0f18f7406056`, from `gong-test (1).json`.
The standard target now loads the [quiet-onset refinement](TfPercussion-gong-early-highs.md#published-refinement).
The **Texture trials** picker preserves the original and two alternatives:

| Trial | Packet spread | Ridge movement depth |
|---|---:|---:|
| gong test | 1.8 | 1.5 |
| Gong — tighter packets | 1.4 | 1.5 |
| Gong — tighter packets + movement | 1.4 | 2.2 |

Everything else is unchanged: all 32 painted frequencies and prominences,
allocation weights, initial excitation, bloom, T60, gains, and saved strike.
Phase blur and output EQ remain off. No DSP mechanism or control was added.
These are listening trials, **not a claim of an overall calibration improvement**.
The user's file in Documents is untouched.

Narrower packets concentrate the surrounding oscillators closer to their
handles; this is not extra damping or a high-pass filter. The second trial
adds bounded ridge movement without adding phase blur. Compare the first trial
before deciding whether the extra movement is useful.

## What the implementation actually produces

The new developer-only `triggerfish_modal_inspect` utility calls the real C++
parameter preparation, rather than approximating the allocation in Python.
For this snapshot at 44.1 kHz it reports **512 active oscillator states**.

| Frequency band | States | Median adjacent spacing | Largest adjacent spacing |
|---|---:|---:|---:|
| 20–800 Hz | 68 | 8.48 Hz | 66.06 Hz |
| 800–2000 Hz | 72 | 15.12 Hz | 56.69 Hz |
| 2000–5000 Hz | 121 | 14.17 Hz | 90.17 Hz |
| 5000–9000 Hz | 118 | 22.49 Hz | 176.08 Hz |
| 9000–15000 Hz | 111 | 35.31 Hz | 280.61 Hz |

The remaining 22 states are above 15 kHz. A large total count does not guarantee
locally dense or perceptually smooth coverage. For example, the 8.57 kHz handle
has 18 states spanning approximately 6.20–11.83 kHz in the original snapshot.

At the 2.37 kHz handle, approximately 10% of its squared input weights belong
to the central doublet; at 3.95 kHz this falls to 3.2%, and at 8.57 kHz to 0.16%.
These are allocation fractions, **not measured radiated energy fractions**.
Thus the dedicated central-doublet control cannot explain all the upper ringing:
the numerous surrounding modes also interfere. Their spacing, movement and
observation matter.

The user's prominence curve is a broad upper-mid scoop followed by a high shelf,
not simply a low-pass rolloff. The render is substantially weaker than the
reference around 2.5–5 kHz. Automatically refilling that scoop to improve a
spectral score would ignore why the user removed it: objectionable timbre.

## Research: implications and limits

Ducceschi and Touzé model gong/cymbal sounds using nonlinear plate modes,
energy-conserving integration and frequency-dependent damping. One illustrative
large circular plate has its thousandth mode near 5.75 kHz. This supports
investigating modal coverage, but its geometry is not established for our
reference and it does not prescribe a universal oscillator count.
[Paper, JSV 2015](https://www.mdphys.org/PDF/jsv_2015.pdf).

Their earlier modal/finite-difference comparison identifies audible consequences
of truncation: a 150-mode example stops transferring energy near 5 kHz and loses
high-frequency richness. It also describes symmetry-dependent nonlinear
coupling. More arbitrary random modes are therefore not automatically equivalent
to a better physical model.
[Ducceschi, Touzé and Bilbao, SMAC 2013](https://www.pure.ed.ac.uk/ws/portalfiles/portal/11221992/smacsmc2013_submission_286.pdf).

Measurements of Balinese gong ageng identify beating from nearby natural modes
and nonlinear harmonic/sum/difference components. That is a different instrument
from the Dresden reference, but demonstrates why independent phase noise need
not reproduce metallic shimmer. Our passive energy redistribution and bounded
phase movement do not reproduce those intermodulation relationships exactly.
[BYU experimental research](https://physics.byu.edu/faculty/gee/gamelan).

Working hypothesis: local modal coverage and correlated spectral motion deserve
attention before adding more noise. The recordings and tests here do **not**
establish which missing physical mechanism is responsible, nor prove that a
larger state budget would solve it. No new physical-model claim is made.

## Experiment and acceptance criteria

Render the actual saved JSON through WASM, preserving its strength 0.950279,
location 0.536998, hardness 0.35, mallet 0.5 and seed 1796. The reference remains
Gong Dresden03 at its fixed −6 dB gain, with verified SHA-256
`a36721cff1f77c22484ac026330aa4da953826ff0e92376d7fb756e27d945147`.
Do not normalize each candidate or substitute the factory gesture.

The sweep checks movement depth/rate, packet width, density, distribution,
central beating, and smooth allocation tilts. It then crosses width with
movement. All painted frequencies and levels stay fixed. Two seeds, 1796 and
1982, reduce dependence on one random realization.

Diagnostics separate two questions:

- **Texture:** RMS error in auditory-band envelope modulation power, expressed
  in dB, in 8–32 and 32–128 Hz modulation bands, over 0.6–1.6 and 1.6–5 s.
  Six analysis bands span 2.5–12.5 kHz. This is not a validated perceptual loss.
- **Timing/body:** existing regional envelope and spectral comparisons, retaining
  the user's low body as a target below 800 Hz, blending towards the reference
  through 1.8 kHz, then using the reference above. Band-mean bias is separated
  from envelope shape; absolute errors are also recorded. Removing bias from a
  diagnostic does not change audio gain.

| Candidate | Timing/body score ↓ | Texture error, dB ↓ | Low-body envelope change, dB RMS |
|---|---:|---:|---:|
| User baseline | 4.34 | 3.21 | 0 |
| Tighter packets | 4.51 | 2.40 | 1.24 |
| Tighter packets + movement | 4.62 | 2.31 | 1.32 |

The trade-off is explicit: modulation statistics improve, but the combined
timing/body score does not. Extra movement gives only a small additional
improvement once packets are narrower. Other tested layouts and upper-biased
allocations did not improve this texture diagnostic; reducing density worsened
it. This does not rule out a redesigned distribution or a larger budget.

A subsequent bloom/T60 search selected a three-control candidate with bloom
2.8 rather than 3.2. Its very small timing benefit did not justify obscuring the
comparison: **that candidate is not published in the picker**. Its holdout
artifacts must not be mistaken for validation of the two published trials.

## Reproduction and artifacts

- Build native inspection/API tests with `./dev.ps1 test-workbench-api`.
- `tools/study_saved_gong_texture.py FIT --output DIRECTORY` runs the scalar
  screen; `--joint` runs width × movement; `--finish` explores dynamics. Run in
  the project Python environment with `PYTHONPATH=python` and `EMSDK_NODE` set.
- `build/gong-upper-mid-study/` contains actual mode lists, fixed-gain WAVs,
  envelope plots, parameter sets, per-seed metrics and renderer provenance.
  Its `texture/`, `joint/` and `final/` directories are experimental artifacts.
- `tools/audit_fit_limiter.mjs FIT_ARRAY OUTPUT INDEX` checks each published
  trial through the actual WASM engine and browser limiter, silently. The
  published audits are under `build/gong-upper-mid-study/published/`.
- `workbench/tests/calibration_ui_probe.mjs gong-standard` checks that the main
  target and trial picker restore all parameters, reference and saved gesture.

No Python, inspection binary or analysis dependencies are required by a normal
VCV build. The existing workbench server is reused.

Validation: 17 focused Python tests pass, including inspection, exact saved-fit
rendering, structured modal levels and upper-sizzle diagnostics. The browser
probe passes for the standard Gong and all three picker entries. At 48 kHz,
both published alternatives render finite single and repeated strikes. Neither
engages the limiter at master −12 dB. Four strikes 0.5 s apart at master 0 dB
produce approximately 0.49 dB and 0.32 dB maximum gain reduction respectively;
single strikes do not. Use −12 dB master for the controlled listening comparison.
