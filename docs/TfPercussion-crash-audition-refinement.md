# Crash: broader shimmer audition refinement

This follows the user's listening feedback on the low-blur crash: a little
more phase blur and packet spread, less low-mode prominence, adjusted decay,
and more crash development. It changes the preset, not the DSP. Gong and
user snapshots remain untouched.

## Search and constraints

`tools/refine_crash_audition.py` starts from the verified
`build/crash-linewidth9/candidate` checkpoint. It renders the actual workbench
Wasm against the same medium edge reference, with the same source gain,
gesture and output gains. No loudness matching is performed.

Six normalized coordinates describe shared edits. The first five become
ordinary UI controls; the sixth is a broad attenuation curve baked into the
visible modal bars, not a hidden runtime parameter:

| Coordinate | Search range | Previous → published |
|---|---|---|
| Phase blur | 0.002–0.006 | 0.0015 → 0.00495 |
| Packet spread | 2.4–3.3 | 2.189 → 2.623 |
| Diffusion strength | 2.8–5 | 2.361 → 2.829 |
| Low T60 endpoint | 18–28 s | 25.827 → 27.868 s |
| High T60 endpoint | 0.4–0.85 s | 0.655 → 0.611 s |
| Low prominence edit | −7 to −1 dB | −2.00 dB below 250 Hz |

The prominence cut tapers linearly in log frequency to zero at 450 Hz.
It changes three bars: −2 dB at 120/240 Hz and −0.759 dB at 360 Hz.
Frequencies, contact parameters, all other bars, packet layout/noisiness,
allocation weights and energy dependence are fixed. Two T60 points remain;
there are no extra knots or per-mode decay changes.

The search screens 47 deterministic Latin-hypercube points and the box centre,
then uses bounded Powell refinement for at most 100 further evaluations. It
retains the unchanged starting point as an explicit candidate. Each score is
the mean of two separately rendered phase seeds, not their averaged audio.
The scalar objective and its reference-fixed analysis floor are unchanged
from [the low-blur methodology](TfPercussion-crash-low-blur-fit.md).

Training seeds are 1396978464 and 1396978771. The final audit uses the standard
seed and three additional seeds absent from this search: 1396979375,
1396980065 and 1396980667. These are same-sample robustness checks, not
validation against different cymbals or velocities.

## Result and limitations

The two-seed proposal score falls from 1.83851 to 1.81072. Independent
envelope-modulation texture scores improve on all four audited realizations.
However, reference-floor spectral comparison worsens slightly on three of
four, and band-decay error increases on all four (approximately 0.20–0.37 dB).
This is a listening-directed texture alternative, not a universal metric win.

The inspected absolute band envelopes and paired/difference spectrograms
still show missing early 300–700-Hz energy, mismatched low ridges, and a weak
late low tail. The main upper decay contour remains similar. No claim of
listening approval is made. Stronger diffusion can redistribute energy faster
without necessarily shortening the audible upper tail; the damping and
redistribution controls are not interchangeable.

The published **Crash — broader shimmer trial** reproduces
`build/crash-audition-broader/candidate` exactly. Eight quarter-note hits peak
at −3.25 dBFS before the browser master, and rapid hard strikes at +8.75 dBFS
in floating point (−3.25 dBFS after the normal −12-dB master). The existing
browser safety limiter is unchanged.

## UI and validation

Help text now leads with audible effects, control direction and practical
interactions. Equations remain in architecture documentation. Spread changes
tone spacing, density changes tone count, and phase blur softens stable
ringing; the tooltips distinguish these explicitly. Longer tips use larger
text and flip above controls when they would run off the bottom of the screen.

All 555 Python tests, 14 Wasm tests and two native API tests pass. Browser
checks cover preset loading, help text, keyboard-triggered tooltips and viewport
placement without touching the user's tab or playing audio.
