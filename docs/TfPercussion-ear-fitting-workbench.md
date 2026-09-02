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
- hardness/implement coordinate;
- mute; and
- event seed or variation index.

They are always visible beside the strike pad and are stored with an audition
gesture, but they do not silently alter the fitted object. The default pad maps
horizontal position to location and vertical position to strength. Either axis
can instead select hardness/implement from a dropdown. Clicking triggers;
dragging may generate a rate-limited sequence only in an explicitly labelled
repeat mode.

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
identify them. The always-visible fitting controls are:

| Control | UI | Perceptual role | Initial mapping |
| --- | --- | --- | --- |
| Model level | bipolar dB slider | match the stored source level | one pre-monitor output trim; never the limiter |
| Impact character | 2D pad | tonal/ping to diffuse/noisy; compact to broad contact | coordinated pulse, chirp, noise and micro-contact gains/durations; hardness remains an event input |
| Bloom | 2D pad | weak to strong nonlinear spread; immediate to slower development | dispersion drive/excursion and constrained feedback/decay mapping |
| Body character | 2D pad | resolved/tonal to dense/wash; dark to bright | equal-power sparse/dense balance plus broad dense tilt and radiation colour |
| Modal structure | spectral point editor | persistent ridge frequency and salience | sparse-mode frequency and output residue; selected-point detail exposes T60 |
| Wash spectrum | low-order ERB curve | broad unresolved spectral distribution | 8 perceptual nodes resampled into the internal dense gain envelope |
| Decay spectrum | low-order ERB curve | frequency-dependent T60 | one shared curve sampled at sparse modes and into the dense decay envelope |

Three 2D controls are preferable to a dozen coupled sliders only if their axes
remain perceptually separable. Each pad therefore has independent X/Y reset,
numeric readout, shift-drag fine adjustment, and a one-axis sweep audition. If
testing shows curved or interacting useful regions, the background displays a
light response field and the mapping uses a bounded two-dimensional surface,
not a misleading pair of linear labels.

The curve editors use log-frequency/ERB horizontal spacing. Nodes are few and
smooth by default; alt-click adds a node only up to a fixed small limit. A
vertical ribbon shows the editable decay uncertainty/range where applicable.
The underlying 33-point dense profile is an implementation detail and is never
shown as 33 unrelated sliders.

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
cursor. Structural macro commits happen in the DSP worker and become audible
after the published lead block.

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
synthesis analysis only, never the reference analysis.

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
number. On reference selection, the workbench matches synthesis and reference
RMS over their first second with the model's linear output trim, then renders
again at that value. Double-clicking model level repeats that match; all other
sliders, pads, selectors, and monitor controls reset to their declared defaults.

In the default mirror view, an unmodified wheel pans both time origins together.
Equal forward and backward wheel deltas cancel exactly and are not clamped at
the onset. Ctrl+wheel performs shared time zoom, dragging either half adjusts
only that half's alignment, and double-clicking the plot restores both time
origins and the full viewport.

## Implementation order

Completed foundations include the separate Emscripten target, deterministic
native/Wasm comparison, persistent real-time triggering, A/B reference
transport, safety limiter, corpus browser, asynchronous STFT, comparison views,
and JSON snapshots. Remaining work is modal/wash curve editing, undo/redo and
snapshot A/B, rolling live-output capture, true-peak limiter qualification, and
macro qualification before the first retained listening fit.
