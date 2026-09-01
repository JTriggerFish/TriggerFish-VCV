# Crash cymbal fitting methodology

## Frozen first object

The first complete metallic instrument is an 18-inch Zildjian K
Constantinople crash from the installed Superior Drummer 3 Core library. Its
container index proves that edge, bow-tip, bow-shank, bell-tip, and bell-shank
are separate recordings with 9--16 velocity layers and multiple round robins.
The machine-readable selection is
`data/crash-calibration/sd3-18-k-constantinople-v1.json`.

The source library is to be rendered locally through SD3. Humanization, internal
processing, velocity randomization, and normalization are disabled. Source
audio and derived isolated samples remain ignored and are never distributed.
Bitwig, Iowa, Salamander, and VCSL are held-out checks, not cells mixed into a
synthetic same-object grid.

## Fixed synthesis graph

The first fitted graph is intentionally compact:

```text
contact direct --------------------------------------> observation
       |
       `-> nonlinear dispersion loop
                       `-> one 32-line coupled wet resonator bank
                              `-> four disjoint observation submixes
                                     `-> small signed frequency shifts -> observation
```

There is no audible dry feed from the dispersion stage. The four submixes are
disjoint views of one body, not four bodies or four independently tuned modal
banks. Small signed shifts reduce persistent harmonic alignment without moving
the shared low-frequency delay structure wholesale. Bow, bell, and edge are
excitation projections into that body. The
contact and body observation paths have separate gains and radiation filters.
Mute changes losses inside both feedback systems and cannot excite them.

The initial fit keeps size, implement family, and neutral hardness fixed. This
does not freeze them out of the architecture: delay frequencies, decay bands,
contact macro parameters, and location projections remain explicit parameters.
It prevents the optimizer from explaining one cymbal by silently changing
several physical/perceptual axes at once.

## Cells and split

The capture plan renders all five articulations and the full 16-value MIDI
velocity probe twice. `build_toontrack_cymbal_sweep.py` creates the MIDI and
nominal-onset manifest; `slice_cymbal_capture.py` detects actual onsets and
creates level-preserving cells with 50 ms of pre-onset noise retained. One
repeat per cell is used for fitting. The second repeat measures
irreducible take variation and supplies the validation tolerance. Exact MIDI
velocities are inputs; discovered SD3 layer identity is descriptive metadata,
not an assumed linear strength scale.

The first acceptance matrix is:

- edge at all sampled velocities for nonlinear onset and bloom;
- bow-tip and bell-tip at low, medium, and high velocity for radial location;
- bow-shank and bell-shank at medium and high velocity for hardness/contact;
- the remaining cells and all second repeats as held-out validation.

## Staged fit

Every stage reports contact (0--15 ms), bloom (15--250 ms), early body
(250 ms--1.5 s), and tail (after 1.5 s) separately.

1. Fit onset alignment, unnormalized early level, crest, and contact spectrum
   with the resonator observation muted.
2. Fit velocity-to-incident-energy and nonlinear dispersion controls from edge
   bloom time, ERB trajectories, spectral flux, flatness, and early level.
3. Fit resonator tuning, coupling, and three-band decay scales from persistent
   modal evidence, peak-density statistics, ERB decay, and late spectral crest.
4. Fit bow/bell/edge excitation projections jointly. A location fit must
   improve held-out velocity cells without changing the shared resonator poles.
5. Fit the two observation/radiation paths and absolute output gain last.
6. Validate mute and then expand hardness, implement, and size on independent
   libraries. These controls are never optimized to rescue a failed neutral fit.

No single scalar loss approves a fit. The first implemented Plotly report shows
aligned waveforms, common-scale spectrograms, and reference/synthesis audio
side by side. It also renders one-second stateful velocity, location, and
hardness sweeps with the other controls held neutral. ERB trajectories, band decays, stable modal evidence, named-loss
tables, and ablations remain report gates before calibration can be accepted.

## Reproducibility and acceptance

The C++ engine is the sole renderer for fitting, tests, Rack, and future
WebAssembly. Fit files record parameter bounds, optimizer seed, source hashes,
code commit, sample rate, segmentation, transform versions, and every named
loss. A candidate is rejected when it improves the aggregate objective while
worsening an audible region outside measured repeat-to-repeat variation.

The currently generated Bitwig edge fit is a pipeline and graph check, not the
primary calibration. The instrument is not called calibrated until the SD3
audio grid has been captured, fitted, and compared with its second repeats,
and until it passes the held-out Bitwig
18-inch crashes and the differently sized Iowa crashes without per-sample
retuning. That check is what prevents a flexible feedback network from merely
memorizing one recording chain.

The reproducible provisional result, including input hashes, optimizer seed,
parameters, and objective, is stored in
`data/crash-calibration/bitwig-a-custom-18-edge-preliminary-fit.json`. It is
kept out of the production defaults so an edge-only result cannot silently
define bell, bow, or hardness behavior.
