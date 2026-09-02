# Percussion ear-fitting workbench

## Purpose and boundary

This workbench is an optional developer application for constructing and
auditioning percussion instruments. It is not part of the Rack plugin and is
not built by `make`, `make dist`, package jobs, or the default standalone CMake
configuration. The browser runs the same C++ DSP as Rack through WebAssembly;
there is no second JavaScript synthesizer.

The first target is the experimental mono `CrashCymbal` renderer. The design is
instrument-neutral so later ride, hi-hat, snare, and kick graphs can supply
their own control descriptors without replacing the workbench.

## Separate kinds of control

The UI must not confuse a performance input with a fitted object parameter.

Performance controls create an event or change live state:

- strength;
- radial strike location;
- implement family (Brush, Mallet, or Stick);
- hardness within that implement family;
- mute; and
- event seed or variation index.

They are always visible beside the strike pad and are stored with an audition
gesture, but they do not silently alter the fitted object. The default pad maps
horizontal position to location and vertical position to strength. Three radio
buttons select Brush, Mallet, or Stick. One adjacent Character slider is
contextual: bristle stiffness, mallet firmness, or tip hardness. The brush
family suppresses the stick chirp/impulse and routes a correlated, smoothly
windowed bristle-contact stream to both direct sound and body drive. Brush
gesture duration remains an internal event property until trigger/gate duration
controls it directly; it is not exposed as an unrelated fitting slider. Clicking
triggers; dragging may
generate a rate-limited sequence only in an explicitly labelled repeat mode.

Fit controls change the instrument definition. The default surface exposes
perceptual macros and compact curve editors. Every raw C++ parameter remains
available in an advanced diagnostic panel, with its units and legal range, but
raw controls are not the normal fitting interface.

Monitor controls change only audition presentation. Reference and synthesis
share the same monitor gain and safety limiter. Monitor gain, level-match mode,
solo, and limiter metering are never serialized as instrument parameters.

## First single-hit fitting surface

The first pass intentionally fits one qualified reference hit. Velocity and
location response curves are locked in this mode because one hit cannot
identify them. The renderer still has a conservative built-in velocity law:
stronger contact is shorter and brighter, high-frequency modal projection
opens progressively, and the nonlinear bloom/turbulent wash grows faster than
the direct strike. This preserves an immediately playable instrument while a
multi-velocity grid is used later to calibrate the exact relationships. The
always-visible fitting controls are:

| Control | UI | Perceptual role | Initial mapping |
| --- | --- | --- | --- |
| Model level | bipolar dB slider | match the stored source level | one pre-monitor output trim; never the limiter |
| Contact | independent sliders | near-field level, ping/noise balance and contact width | direct presentation plus coordinated pulse, chirp, noise and micro-contact gains/durations; hardness and implement remain event inputs |
| Bloom | independent sliders | nonlinear spread amount and development time | dispersion drive/excursion and constrained feedback mapping |
| Body balance | sliders plus branch toggle | resolved/tonal to dense/wash; global wash tilt; explicit resolved-mode ablation | sparse/dense balance and broad dense tilt |
| Resolved modal ridges | 12-bar ERB editor plus enable switch | direct placement of deliberately audible modes | each bar controls one resolved mode only; disabling the branch makes this editor acoustically inert |
| Dense modal wash | smooth ERB curve plus continuous sliders | broad spectral colour, range, density and statistical irregularity | the relative colour curve changes modal energy without moving resonances. Density continuously spans the equivalent of 64--4096 active modes. The implementation activates a deterministic nested subset and crossfades the boundary mode; the original 2048-mode bank remains the `1x` factory case and values above it progressively add a separately seeded extension bank. |
| Body T60 | five-knot log curve | absolute frequency-dependent decay | one curve sampled by sparse modes, dense modes and turbulence |
| Turbulence | three-knot log curve plus amount/persistence | stochastic residual colour and duration | a bloom-driven branch parallel to the dense modal cloud; mean-free colour curve, orthogonal overall level and relative persistence |
| Radiation | per-path EQ controls | observation colour, not stored-body loss | direct, resolved and dense-path static filters |
| Size meta | unlabelled scalar slider | useful broad starting point | expands visibly into the detailed controls; centre is the exact neutral default and saved fits contain only expanded values |

The curve editors use log-frequency or ERB horizontal spacing and raised-cosine
interpolation, so moving one knot has a local, smooth effect without spline
overshoot. Frequencies remain ordered. Resolved-mode and dense-wash curves are
independent and relative because the modal banks are energy-normalized; global
balance and level are controlled elsewhere.

Curve backgrounds are inert by default. A handle drag starts only after a small
movement threshold and remains relative to its starting value. The resolved editor
adds an explicit **Draw bars** mode for a continuous sweep begun anywhere in the
grid; Shift erases, while Alt temporarily enables the same gesture without
latching the mode. Mirror-spectrogram source alignment similarly requires
Shift-drag, and the decorative waveform does not capture background drags.
Turbulence subtracts the mean of its three colour levels before applying its
separate Amount control. The underlying 33-point dense profile is an
implementation detail and is never shown as unrelated sliders.

## Multi-hit and generalization controls

After several independently saved single-hit fits exist, the workbench can
open a grid view and fit relationships:

- a velocity-response curve for level, contact bandwidth, bloom and decay;
- bell/bow/edge projection curves over the shared modal set;
- implement/hardness response for contact duration, tonal/noise balance and
  body coupling; and
- later, size/tune mappings that distinguish sparse-frequency movement, dense
  spectral scaling, contact colour and decay instead of shifting everything
  linearly.

These relationship editors must be derived from named fitted cells. They do not
appear as plausible-looking controls before the dataset contains the evidence
to constrain them.

## Control qualification

A macro is admitted to the default surface only after deterministic sweep tests
establish all of the following:

1. Safety: every endpoint and a dense interior grid stays finite and within the
   declared pre-limiter headroom under isolated and repeated hits.
2. Audibility: its full range creates a material change in the intended time
   region, including a level-matched audition render.
3. Direction: named proxies move consistently—for example, contact extent,
   spectral centroid, ridge/residual ratio, bloom time, or band T60.
4. Local smoothness: adjacent control steps do not create unexplained clicks,
   state resets, seed changes, or large discontinuities.
5. Independence: range-normalized audio/feature derivatives are not nearly
   collinear with another default control. Redundant controls are merged or
   moved to the advanced panel.
6. Usefulness: min/centre/max and slow sweeps are included in an audition page;
   passing numerical checks alone is insufficient.

The qualification output records effect size per contact/bloom/body/tail
region, monotonicity violations, maximum adjacent-step change, derivative
correlation, render time, and the exact macro mapping version.

## Workbench layout

The spectrogram is the central, persistent surface. On a desktop display:

```text
+---------------- reference / fit / snapshot bar ----------------+
| reference browser | transport + A/B | snapshot / compare / save |
+------------------------+-----------------------------------------+
| performance + fitting  | waveform + spectrogram                  |
| controls (collapsible) | reference / synth / mirror / difference |
|                        | synchronized zoom and cursor             |
+------------------------+-----------------------------------------+
| strike pad + sequence  | analysis settings + branch/limiter meter|
+------------------------+-----------------------------------------+
```

The control column remains narrow and scrolls independently; transport and the
currently selected macro stay pinned. The plot never disappears while editing.
Reference/synthesis audio buttons, A/B switching, a short selected-region loop,
and contact/body/tail region shortcuts remain next to the plot rather than at
the bottom of the page.

## Snapshots and fitted states

Every edit changes a working state. A snapshot is immutable and records:

- snapshot UUID, name, timestamp, notes, and optional parent UUID;
- instrument graph ID, parameter schema, macro-mapping version, and code commit;
- complete low-level C++ parameters plus the perceptual macro state;
- reference corpus alias, cell ID, articulation, velocity, repeat, source hash,
  transform/alignment version, sample rate, and pre-onset duration;
- audition gesture or hit sequence and random seeds;
- spectrogram settings and saved viewport; and
- optional cached audio/analysis hashes.

Undo/redo covers the working state. “Snapshot” places an immutable marker in a
comparison tray. Selecting it restores controls and can render old and current
states side by side. “Save fit” writes a versioned JSON document; loading first
validates schema, source identity, and renderer compatibility and never clips
unknown values silently. Browser storage is convenient recovery only. Explicit
JSON export is the durable source of truth.

## Rendering and analysis scheduling

The C++/WebAssembly renderer runs continuously in a dedicated worker and fills
a lock-free shared ring. An `AudioWorklet` consumes 128-frame browser quanta
from a 512-frame target lead; neither DSP rendering nor queue waiting runs on
the DOM thread. The worklet wakes the producer with shared-memory atomics rather
than relying on timer polling. A strike addresses the persistent renderer state
after the already-published 512-frame lead, so repeated hits accumulate in the
same cymbal with deterministic low latency. The worklet owns the read index and
the DSP worker owns the write index; neither side rewrites the other side's
cursor. Structural macro commits happen synchronously in the DSP worker and
become audible after the published lead block. The header shows an elapsed
“Preparing live DSP” timer while a rebuild is pending, then retains the
measured completion time. A later shadow-state swap can remove this short gap.

The current corpus browser exposes one qualified private 320-cell grid as an
instrument: five articulations, sixteen velocities and four variations. Only
the selected WAV is fetched. Its fixed +25.5 dB audition trim comes from the
preserved-level 64--112 velocity range and is applied equally to reference and
synthesis; individual cells are never peak-normalized.

Offline comparison renders use a second worker, independent of the live audio
engine. Slider changes are coalesced while it is busy, so only the newest queued
state is rendered and the DOM remains responsive. Spectrogram calculation is
also asynchronous and generation-tagged: stale jobs are discarded rather than
painted after a newer control state. Rendering a new sound invalidates the
synthesis analysis only, never the reference analysis. The reference STFT peak
fixes one absolute colour ceiling for both mirror halves. Synthesis edits
cannot re-normalize that domain, so a low-band change cannot recolour unrelated
high-frequency content.

## Audition safety

The listening path is:

```text
reference or synthesis -> shared audition trim -> lookahead limiter -> master -> output
```

The initial limiter uses a 5 ms lookahead, immediate downward gain, a lookahead
hold, 100 ms smooth release and a conservative -1 dBFS sample-peak ceiling.
Gain reduction, audio-context state, output level, acknowledged strikes and
queue underruns are visible. It starts enabled and cannot be bypassed from the
ordinary UI. A later true-peak qualification pass may replace its detector
without changing the graph boundary. An explicit offline diagnostic may export
pre-limiter audio to confirm that the model is stable and not being “fitted
into” protection.

Analysis defaults to the pre-limiter renderer output. Reference and synthesis
must use one declared audition-gain policy; independent peak normalization is
never the default because it hides velocity and decay-level errors.

The model-level default is reference-relative rather than a fixed factory
number. On reference selection, the workbench aligns the declared source onset
and matches synthesis and reference RMS over the first 300 ms with the model's
linear output trim, then renders again at that value. This keeps a long wash
from turning down and perceptually distancing the contact. Double-clicking
model level repeats that match; all other sliders, pads, selectors, and monitor
controls reset to their declared defaults.

In the default mirror view, an unmodified wheel pans both time origins together.
Equal forward and backward wheel deltas cancel exactly and are not clamped at
the onset. Ctrl+wheel performs shared time zoom, dragging either half adjusts
only that half's alignment, and double-clicking the plot restores both time
origins and the full viewport.

## Implementation order

Completed foundations include the separate Emscripten target, deterministic
native/Wasm comparison, persistent real-time triggering, A/B reference
transport, safety limiter, corpus browser, asynchronous STFT, comparison views,
the resolved-mode/dense-colour/shared-T60/turbulence editors, and JSON snapshots. Remaining
work is undo/redo and snapshot A/B, rolling live-output capture, true-peak
limiter qualification, and macro qualification before the first retained
listening fit.
