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
routes are energy-normalized for their much larger number of contacts, with
more energy sent into the cymbal body than the dry near-field presentation.
Stiffness continuously moves a fixed seeded gesture from dark, overlapping
fine-bristle contact toward a brighter, more articulated coarse-bristle sound;
it does not regenerate stochastic contacts or drive the full stick-strength
nonlinear bloom.
Brush gesture duration remains an internal event property until trigger/gate duration
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
| Bloom | independent sliders | dispersed-body level, nonlinear character, development time, and diffusion | a true body-route gain plus independent dispersion drive/excursion, constrained feedback, and serial-allpass amount |
| Body model | separate focused and legacy-diagnostic views | one anchor-driven stochastic field by default; the older sparse/dense/noise branches remain only for comparison | modal-field turbulence, ERB packet spread, phase bandwidth and passive neighbour exchange |
| Modal field anchors | wide 24-anchor ERB editor from 40 Hz to 15 kHz | direct placement of coherent centre lines and their constructive stochastic neighbourhoods | each anchor controls frequency, packet energy, and a `0..2x` response to global turbulence |
| Dense modal wash | smooth ERB curve plus continuous sliders | broad spectral colour, range, density and statistical irregularity | the relative colour curve changes modal energy without moving resonances. Density continuously spans the equivalent of 64--4096 active modes. The implementation activates a deterministic nested subset and crossfades the boundary mode; the original 2048-mode bank remains the `1x` factory case and values above it progressively add a separately seeded extension bank. |
| Body T60 | editable two-to-eight-knot ERB/log curve | absolute frequency-dependent decay | fixed DC/Nyquist boundary positions with editable finite T60 values; one prepared curve sampled by every body mode |
| Unified turbulence | one primary slider, a per-anchor response curve, and advanced trajectory controls | continuously broadens selected anchors into modal wash and finally noise-like response while others can remain clean | normalized satellite energy, ERB spread, passive phase decorrelation and local orthogonal energy exchange |
| Radiation | compact bandwidth/colour controls with an Advanced section | observation colour, not stored-body loss | direct and body static filters; enable state and Q remain reachable until ablation establishes that fixed values are sufficient |
| Size meta | unlabelled scalar slider | useful broad starting point | expands visibly into the detailed controls; centre is the exact neutral default and saved fits contain only expanded values |

The spectral-colour editors use log-frequency or ERB horizontal spacing and
raised-cosine interpolation. The body-decay editor instead uses straight
segments in ERB versus log-T60 space, exactly matching the DSP preparation
rule; a new fit starts with the two boundary knots, double-click adds or removes
interior knots, and its central `ALL` handle translates every knot in log time.
Frequencies remain ordered. Anchor energies are normalized as a
group, so their dB bars define relative body colour while the exact `-72 dB`
floor deactivates an anchor. Global body and model levels are controlled
elsewhere.

The focused modal editor presents each constructive packet as a narrow centre
line with a translucent bell indicating its turbulence-controlled frequency
neighbourhood. Dragging a centre changes frequency and energy; dragging a wing,
Ctrl-dragging the centre, or using the mouse wheel changes its local turbulence
response. Double-click inserts or removes anchors. Draw mode sweeps a Gaussian
brush across all active anchors; Shift narrows and Ctrl widens the brush.
Presets and Clear alter only anchor frequency, energy, and local turbulence.
The spectrogram height has its own draggable separator.

Mirror-spectrogram source alignment requires Shift-drag, and the decorative
waveform does not capture background drags.
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

The C++/WebAssembly renderer runs directly inside an `AudioWorkletProcessor`.
Each process callback advances the persistent cymbal state for the browser's
render quantum, currently 128 frames in supported browsers, and writes those
samples straight to the worklet output. There is no Worker audio producer,
shared ring, timer polling, or extra TriggerFish render lead. Repeated strikes
therefore address the same state at the next message/worklet boundary rather
than after a pre-rendered queue.

Emscripten emits a dedicated single-file `worklet` build of the same C++ source
as the offline renderer. Initial startup outputs silence until its Wasm module
is ready. A later parameter or routing edit never re-prepares the active DSP
instance. An ordinary Worker expands the high-level controls and compiles the
modal recurrence coefficients into an immutable, version-one prepared blob.
The AudioWorklet then performs only bounded component installation and state
initialization in a muted standby processor; it does not evaluate 408 modes'
trigonometric or exponential coefficient formulas. The old processor keeps
sounding and crossfades to the ready candidate over 3 ms. Rapid edits are
coalesced, stale candidates are destroyed, and retirement is serialized so
each AudioWorklet module's fixed session pool contains at most the active and
one candidate. Prepared blobs are internal same-build messages, not portable
fit files. The unreleased C/WebAssembly ABI remains version one until the first
release. The header reports “Preparing live DSP” until the requested
generation becomes active.

The top-right settings dialog owns workstation concerns rather than instrument
parameters. Its MIDI panel requests Web MIDI access on startup, enumerates and
hot-plugs connected inputs, and permits either all inputs/channels or one
device/channel. Browsers that require a permission gesture retain the Enable
MIDI button as a retry. Every note-on triggers the current location and
implement; MIDI velocity replaces strike strength. Note number is deliberately
unassigned until the later pitch-following design is specified.
The same dialog reports the worklet quantum, limiter lookahead, and browser's
best available device-latency estimate. MIDI preferences are local browser
settings and are not part of a fitted instrument snapshot.

The current corpus browser exposes a deliberately small 25-cell view of the
qualified private crash grid: five articulations, five representative
velocities and one repeat. Only the selected WAV is fetched. Its fixed
+42 dB audition trim comes from the preserved-level velocity grid
and is applied equally to reference and synthesis; individual cells are never
peak-normalized. Small allow-listed snare, kick, gong, ride, and hi-hat corpora provide
one standard calibration cell each. A highlighted **Calibrations** chooser in
the persistent toolbar loads the recipe, starting parameters, reference cell,
strike controls, and shared level match together.

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

The initial limiter uses a 3 ms lookahead, immediate downward gain, a lookahead
hold, 100 ms smooth release and a conservative -1 dBFS sample-peak ceiling.
Gain reduction, audio-context state, output level, acknowledged strikes and
queue underruns are visible. It starts enabled and cannot be bypassed from the
ordinary UI. A later true-peak qualification pass may replace its detector
without changing the graph boundary. An explicit offline diagnostic may export
pre-limiter audio to confirm that the model is stable and not being “fitted
into” protection. This is browser audition protection only: it is not compiled
into, or placed in the signal path of, the percussion DSP library or Rack
modules.

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
the modal-packet/shared-T60/body-colour editors, JSON snapshots, MIDI note-on
audition, and direct Wasm rendering in an AudioWorklet. Remaining work is
undo/redo and snapshot A/B, rolling live-output capture, true-peak limiter
qualification, and macro qualification before the first retained listening
fit.
