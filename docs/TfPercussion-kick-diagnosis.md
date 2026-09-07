# Kick: source isolation and control reachability

This investigation concerns the actual **main-workbench** calibration, not a
different report candidate. The preset is not changed by the diagnostic command.
Reproducibility tests do not establish a good acoustic fit.

## Fixed experiment

Target: acoustic-kick-oak, medium centre, velocity 64, take 1. Source SHA256:
`72631be1485f64018f990e2c0b6d62969578a374d6a5a7d83059a2f0ba4b50b8`.
Reference gain is +2 dB; onset trim is 1.315 ms. Model gain stays -12 dB.
Strength/location/hardness/implement are 0.5; spread is 0.2. Rendering uses
the workbench's exact C++/Wasm voice at 44.1 kHz, with no limiter or level matching.
The current fit contains six active modal frequencies, all below 808 Hz.

The three observations are:

- Contact direct: pulse + chirp + enveloped noise + microcontacts.
- Thump: the independent clean pitch-swept oscillator.
- Resonance: the modal body, driven by contact even when direct contact is muted.

Isolate observations by zeroing only the other two observation gains. Their
summed render matches the complete voice within 3.8e-8 peak error. Separately
set source noise amplitude to zero: this tests excitation, not observation.
The difference from the full render is saved as the source-noise contribution.
This does not remove the separate microcontact source.

## What is wrong in the published sound

Measured model-minus-reference band power:

| Frequency band | 0–30 ms | 30–100 ms | 100–250 ms |
| --- | ---: | ---: | ---: |
| 20–120 Hz | -0.5 dB | -0.7 dB | +2.2 dB |
| 120–250 Hz | -5.0 dB | -4.4 dB | -10.8 dB |
| 450–1000 Hz | +0.4 dB | +6.9 dB | +3.0 dB |
| 1–2 kHz | -6.8 dB | -3.0 dB | below useful comparison level |
| 2–4 kHz | -8.9 dB | +1.3 dB | below useful comparison level |

The loud sub-1kHz body is sustained too strongly, while the initial upper
spectrum and longer low-mid body are missing. Reducing all noise or increasing
overall bass gain is not an adequate diagnosis. Removing source noise exposes
how strongly the present body relies on noise drive rather than the short impact.

## Why earlier validation missed it

1. An aggregate score was treated as justification to publish despite visible
   spectral holes. Sample-identical UI playback only validates transport.
2. Broad, shallow analysis filters let strong sub-1kHz resonances contribute to
   nominal upper-band measurements. Narrower bands and steeper filters reveal
   a much larger deficit. Integration across wide bands also hides holes inside
   those bands.
3. Local fitting adjusted only active modes. With all active centres below
   808 Hz and direct contact almost muted, the parameter search did not explore
   the absent upper-body coverage. More iterations in that restricted region
   are not an adequate remedy.
4. Contact colour is a complementary one-pole shelf with limited contrast,
   not an independently adjustable noise bandwidth or low-pass cutoff.

## Tests of the existing controls

`dev.ps1 diagnose-kick` performs observation isolation, source-noise removal,
one-control probes, joint direct-contact/noise-duration probes, and trials of
the existing simple output low-pass. It also tests an explicit coverage layout
using the ten previously unused modal slots (200 Hz through 8 kHz).

With `TF_KICK_SPAN_REFIT=1`, two bounded least-squares trials jointly adjust
contact, relative modal levels and shared damping. Frequencies and thump stay
fixed: this is a targeted reachability test, **not an exhaustive fit**. One trial
bypasses EQ; the other permits the existing ordinary output low-pass. No
multiband EQ, extra DSP capacity or per-mode damping is introduced. Two noise
seeds are optimized separately and three additional seeds are evaluated.
Actual bounds, fixed values, finite differences and solver results are logged.
The implementation reuses `Search`; its residuals are band errors / 3 dB and
log-spaced regional spectral errors / 6 dB, without candidate normalization.

Failure of these finite trials does not prove that the architecture cannot fit
the reference. In particular, their added modal frequencies are starting
placements, not recovered acoustic modes.

The generic coverage trial still has 11.1 dB worst evaluated band error. Adding
the existing low-pass reduces this to 3.9 dB for the primary seed, with a cutoff
near 1.94 kHz and audible direct contact restored. However, its worst spectral
interval is still about 25 dB wrong, and held-out seeds have 6.3–8.5 dB worst
band errors. Both trials are rejected. The low-pass result is evidence that
source bandwidth and modal placement matter; it is not a satisfactory fit.

A separate `TF_KICK_MEASURED_LAYOUT=1` trial warm-starts that low-pass result
and replaces generic added frequencies with reference-spectrum peak proposals.
This includes the previously uncovered 425–485 Hz region. Proposals initialize
frequencies only; their measured power is not copied into modal gain. This trial
also permits the full existing direct-contact gain range (0–4), and adjusts
thump gain/decay while keeping pitch fixed. It introduces no new DSP controls.

## Independent checks, not another acceptance average

`BandRegionAudit` measures seven bands over five time regions through 1.2 s.
It uses causal third-order Butterworth bandpass construction for 20–120 Hz and
sixth-order construction above that. There is no envelope smoothing across
region boundaries. Filters still have delay/memory: these are output comparisons,
not estimates of individual physical-mode T60.

Every region/band cell exceeding -80 dBFS in either signal must be within 3 dB.
`RegionSpectrumAudit` additionally checks 24 logarithmic spectral intervals over
attack and early decay. The 90th-percentile discrepancy must be <=6 dB, and no
evaluated interval may exceed 12 dB. These tolerances are explicit engineering
screens, not hearing thresholds or a perceptual quality model. Passing is
necessary for automatic publication, not proof that the sound is convincing.

Kick publication now stops before replacing the preset if these checks fail.
The primary event and three held-out noise seeds must pass.
Results are written to `quality-checks.json`; an improved aggregate score cannot
override the failure. The current published preset itself fails these checks.

The initial diagnostic search exposed another mask error: its residual selected
only reference-active cells. The independent union-of-reference-and-candidate
checks rejected the excess high-frequency tail that this permitted. `RegionFitLoss`
v2 retains a fixed residual for **every** cell: reference-active cells use signed
error; quiet cells penalize only excess energy above a fixed reference-derived
ceiling, with a 3 dB margin below the evaluation threshold. An identity render
still has exactly zero error. A regression test adds a late 4 kHz tone to a
reference with no such tail; it must produce a large fitting penalty as well as
fail validation. This correction is in the fitter, not the drum DSP.

After that correction, the reference-peak trial has 5.3 dB worst band error,
6.1 dB spectral-error P90 and 14.2 dB worst spectral interval. Held-out seeds
have 7.4–10.1 dB worst band errors. It is also rejected: recovering some attack
coverage has not recovered the complete decay or a stable source balance.

## Subsequent agreed experiment

Do not mistake the noise-disabled diagnostic for removing the contact generator
or a reason to remove and reintroduce it during fitting. The subsequent agreed
work keeps the architecture and its existing sources fixed and compares
established perceptually motivated losses under matched starts and budgets.
See [the perceptual-loss experiment](TfPercussion-perceptual-loss-experiment.md)
for the method, exact-model recovery tests and rejected real-reference trials.
Do not publish any of these rejected candidates as an improvement.

All WAVs, parameter vectors, diagnostics and offline Plotly plots are retained
under `build/kick-diagnosis/`. They are local analysis artifacts, not a separate
audition page. The main workbench remains the listening destination.
