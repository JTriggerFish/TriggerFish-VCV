# Low-ring beating: diagnosis and control design

The subsequent [paired-ring implementation](TfPercussion-paired-ring.md)
implements the recommended character. This document records the preceding
diagnosis and its unchanged-engine experiments.

The reference's slow breathing low ring versus the published crash's rapid
flutter is measurable. This investigation leaves the engine, workbench
presets, gong and saved snapshots unchanged.

## Theory and evidence

Two stable nearby tones beat at their frequency difference; neither an LFO
nor phase blur is required. Unequal levels reduce beat depth; multiple nearby
tones and differing decays produce less regular motion.

[Perrin et al., *The normal modes of cymbals*](https://www.ioa.org.uk/system/files/proceedings/r_perrin_gm_swallowe_sa_zietlow_tr_moore_the_normal_modes_of_cymbals.pdf)
report split modal partners, consistent with imperfect axial symmetry.
[Skare and Abel, DAFx 2019, section 4.1](https://dafx.de/paper-archive/2019/DAFx2019_paper_48.pdf)
describe beating and differing adjacent decay rates in dense high-Q synthesis.
Neither establishes a universal desirable beat rate or proves this recording's
physical modulation mechanism.

The reference's five-second Hann spectrum has prominent peaks near 125.2 and
126.6 Hz (0.2-Hz bins). This separation is consistent with the slow pulse and
motivates testing a close pair, not imposing that separation on every mode.

## Fitting blind spot and new diagnostic

`ModalTextureLoss` starts at 1.5 kHz and excludes modulation below 2 Hz.
It never validated the low-frequency, approximately one-second beating here.

`LowModeBeating` analyzes 90–180, 180–320, 320–550 and 550–900 Hz. Analytic
envelopes have a fitted exponential decay removed before measuring modulation
over 0.5–4.5 seconds, giving 0.25-Hz resolution. It reports power in 0.5–3,
3–8 and 8–30-Hz bands, the strongest modulation bin and absolute band RMS.

Only diagnostic envelopes are normalized, not audio. This is not a standalone
perceptual objective: wrong pitch, decay or perfectly regular tremolo must not
win merely by matching a modulation statistic. Short windows cannot resolve
arbitrarily slow beats. Tests cover known 1.5/12-Hz tone pairs, a single damped
tone, different damping, invalid inputs and 48-kHz resampling.

## Controlled tests

Same reference, event and gains; 90–180-Hz band, standard seed:

| Case | Strongest pulse | Modulation power above 3 Hz |
|---|---:|---:|
| Reference | 1.25 Hz | 6.8% |
| Published crash | 6.5 Hz | 99.4% |
| Phase blur off | 6.5 Hz | 99.7% |
| Diffusion off | 6.5 Hz | 99.8% |
| Both off | 6.5 Hz | 99.9% |
| Explicit close low pair | 1.25 Hz | 1.1% |

Spacing, rather than blur or diffusion, is the leading explanation. Other low
bands have different reference modulation; not all lows should pulse slowly.

The proof uses visible handles at 126.6/125.2 Hz, local noisiness zero, and a
companion observation bar 8 dB lower. Upper controls remain fixed, but adding
a handle reallocates the pool and changes later packet phase seeds: upper DSP
states are not perfectly isolated. The single-seed spectral score worsens
0.924 → 0.963, and the pair is more regular than the reference. This is an
unpublished mechanism test, **not an accepted crash recalibration**.

## Model/control limitation and recommendation

- Wide low packets with relatively few tones have gaps of several Hz, hence
  fast beating even with completely stable phases.
- Reducing local noisiness narrows gaps but also weakens surrounding tones:
  the desired slow beating can disappear instead of becoming controllable.
- **Beating doublets** pairs surrounding tones, not the prominent centre.
  That unpaired centre can dominate beating against distant neighbours.

Recommended next design: an explicit **paired-ring layout**, preserving the
existing layouts. Pair the prominent ring itself and reuse the separation
control as an audible beat rate. Keep broad spread/high-frequency shimmer
separate. Preserve the oscillator budget and normalized excitation energy;
do not add output tremolo, hidden damping or another energy source.

Check zero separation, observation gain and low/high-noisiness behaviour
before adoption. Avoid new per-packet controls unless the shared rate and
existing local controls prove insufficient. The user should choose ring
character and beat rate, rather than balance several width/noise controls.
This recommendation is not implemented by this investigation.

## Reproduction and validation

With the development environment and `EMSDK_NODE` configured:

```powershell
.venv/Scripts/python.exe tools/audit_crash_beating.py
.venv/Scripts/python.exe tools/plot_low_beating.py
```

Artifacts are in `build/crash-low-beating`; the explicit pair snapshot and
independent-seed audit are in its `pair-proof` subdirectory. Standard crash
publication audits now include the low-band diagnostic alongside upper
texture, absolute spectrum, attack and decay, without changing optimization
weights. All 559 Python tests pass. Normal Rack builds need none of this.
