# Prog Sequencer development roadmap

This is an internal design note for planned Prog Sequencer work. It preserves
extension points and implementation constraints from the superseded proposal,
syntax-study, grammar, and current-design documents without retaining their
discarded alternatives or development history.

Prog Sequencer is in beta. Nothing in this document is accepted syntax or a
compatibility promise; names, priorities, and exact spellings may change. The
[Prog Sequencer reference](TfProgSequencer-reference.md) is the user-facing
contract for behavior that works now.

## Architectural boundary

New features must continue through the prepared execution path:

```text
cpp-peglib syntax tree
  -> typed syntax records
  -> checked semantic graph
  -> prepared immutable program and bounded workspaces
  -> allocation-free runtime events
```

Parsing, dynamic name lookup, allocation, locks, and user callbacks do not
belong on the audio thread. Features that can emit additional notes or control
events must establish finite preparation bounds before a program can be
published.

The persisted language remains version 1 throughout beta. Beta grammar and
semantics may change without migration support. Once the language leaves beta,
an incompatible change requires an explicit version and migration decision; it
must not silently accumulate parallel compatibility parsers.

## MIDI and live performance

The central plan is to let a compiled sequence respond to live performance
without rewriting the program or giving MIDI ownership of musical time.

The intended order is:

```text
written note or chord
  -> sequence and pitch transforms
  -> live input binding
  -> effective pitch material
  -> instrument interpretation
  -> notes, voicings, or bounded sub-events
```

Input binding and interpretation are separate concepts:

- `midi.root` retargets a written melody or chord progression from its declared
  tonic to a played root note;
- `midi.chord` supplies the currently held notes as harmonic material;
- an interpreter decides how that material is voiced or played, with directions
  such as `piano.simple`, `piano.rootless4`, `bass.two`, `bass.walk`,
  `arp.index`, `arp.up`, and `arp.down`.

The following illustrates the intended model, but is not accepted syntax yet:

```text
rooted = sequence {
  tonic C
  scale dorian
  notes 1 3 5 b7 6 5 3 2
  input midi.root
}

arp1 = sequence {
  notes 1 2 3 4 3 2
  input midi.chord
  interpret arp.index
}
```

`arp.index` treats the written degree as an index into the ordered held-note
set. The sequence therefore continues to own rests, ties, durations,
probability, subdivisions, ratchets, and variation. `arp.up` and `arp.down`
represent a different policy that may generate several notes inside one source
event span.

A source-first shorthand such as `midi.in |> arp arp1` may also be explored.
An ordinary sequence will not become implicitly callable through an ambiguous
form such as `midi.in |> arp1`.

The `input` and `interpret` lanes are intended to follow the Notes pass like
Velocity or Duration and restart at the same boundary. Explicit parenthesized
voicings remain exact and bypass interpretation by default.

MIDI input will be represented as a timestamped performance state containing
held notes, velocities, the most recent note, sustain/pedal state, and relevant
controllers. It must not alter editor text or the compiled pattern. A program
without live input remains seeded and repeatable; a performance can be replayed
when supplied with the same timestamped input events.

Planned supporting work includes Rack MIDI input selection, routing several
interpreted parts from one module, and allowing CV-derived performance input to
use the same binding model.

## Harmony and instrument interpretation

Jazz and Roman chords already retain semantic roots, intervals, bass notes,
scale, and register information. Planned interpreters can use that information
for contextual voicing and voice leading without placing instrument-specific
rules in chord spelling.

An interpreter is intended to receive the transformed current item, neighboring
items, active scale and tonic, register and range, event span and beat, previous
realized voicing, deterministic seed, and the current performance snapshot.
It must consume the effective transformed chord or pitch, not reconstruct
meaning from the original jazz-symbol text.

Intended applications include:

- simple and rootless piano voicings;
- two-feel and walking bass;
- chord-aware arpeggiation;
- continuity-aware voice leading using the previous realized voicing; and
- simultaneous voicings or bounded rhythmic figures inside the current event.

A scale-stacked chord generator is also planned as a distinct operation, with a
possible form such as `chord 2`. It will not change the meaning of `II`, which
explicitly means a major chord rooted on scale degree 2. Keeping these separate
prevents a scale change from silently changing an explicitly written chord
quality.

A `chords` lane may be added as a readability alias for chord-focused material;
it will not introduce a second incompatible timing model.

## Pattern structure and rhythmic vocabularies

The current language executes atomic alternatives such as `[1|3|5]`. Planned
pattern work includes:

- alternatives whose branches contain several events;
- recursively nested random and alternating branches; and
- a drum or rhythm-mask vocabulary in which standalone `x` can mean a hit while
  `x1` remains a ghosted pitched event in note syntax.

The exact drum-lane surface is still exploratory. It should reuse the prepared,
bounded pattern model rather than overload pitch tokens ambiguously.

## CV and signal expressions

The current module has two monophonic sequenced CV outputs. Planned extensions
include:

- physical CV outputs beyond CV1, CV2, and CV3;
- an explicit output safety clamp rather than an undocumented language-level
  voltage ceiling;
- continuous interpolation on `...`-aligned CV lanes;
- prepared control events that traverse every CV knot even when no note event
  occurs at that time;
- typed per-voice expressions for polyphonic CV; and
- top-level expressions and small typed functions over MIDI, CV, event values,
  and previously computed signals.

Possible signal-expression forms include:

```text
cv2 = exp(0.1 * vel) + cv1
cv3 = midi_pitch * cv2
```

These assignments are not accepted today. When introduced, sequence and
arrangement values, event-rate MIDI values, and control/audio-rate signals will
remain distinct types rather than being inferred from textual context.

The parser should produce typed calls, unary/binary operations, references, and
assignments. Semantic preparation must reject dependency cycles, unavailable
inputs, non-finite results, and unbounded state before producing a fixed-storage
numeric evaluator. Signal expressions form a prepared stage before live input
binding and interpretation.

## Pitch and tuning extensions

Longer-term pitch directions retained from the earlier design work are:

- scales containing arbitrary cent or frequency-ratio offsets without changing
  ordinary degree notation; and
- a separate typed literal for voltage or frequency pitch when direct physical
  pitch is needed.

These are exploratory rather than committed syntax. Named pitches will remain
chromatic, while scale-degree notation will continue to describe positions in
the active scale.

## Editing and advanced control

Potential later improvements include preserving playback phase when a sequence
is renamed, an editor view that expands aligned/polymetric lanes without
rewriting their concise source, and evaluating whether an advanced imperative
control layer is useful alongside the declarative sequence language. Expanded
views must retain source provenance for selection, diagnostics, and playback
cursors. None of these additions should complicate ordinary pattern syntax or
weaken safe, predictable realtime playback.

## Boundaries to preserve

As the roadmap evolves, these musical boundaries should remain stable:

- the compiled sequence owns time; MIDI supplies performance state;
- input binding determines what a controller means, while interpretation
  determines how the resulting pitch material is realized;
- live input never edits source text;
- explicit voicings remain explicit;
- Roman chord case remains an explicit quality choice rather than automatic
  scale stacking; and
- future examples in this document remain unavailable until they appear in the
  main reference.
