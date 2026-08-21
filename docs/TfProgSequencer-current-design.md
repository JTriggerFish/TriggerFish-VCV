# TfProgSequencer current design

This is the authoritative implementation note for Prog Sequencer. The syntax
study records alternatives and the v1 proposal retains historical context;
this document records the boundaries the current code must preserve.

The current persisted language version is 3. Version 3 makes pipeline order
musically significant for `modulate_degree` followed or preceded by `scale`,
and extends explicit named-note voicings without changing their old meanings.

## Compiler architecture

The implementation has four explicit representations:

```text
cpp-peglib document AST
  -> typed syntax Document and domain-neutral PatternNode tree
  -> checked SemanticProgram
  -> prepared CompiledProgram
  -> bounded RuntimeEvent stream
```

`PatternNode` owns reusable structure—atom, subdivision, cycle choice, random
choice, Euclidean generation, and repeat—before a lane assigns musical meaning.
cpp-peglib productions emit these structural nodes directly; there is no
secondary delimiter or nesting parser. Voicing tones, shared register suffixes,
and slash-bass relationships are likewise syntax nodes rather than text split
again by the compiler. Lane compilers turn the leaves into typed pitches/chords,
articulations, or scalars. Unsupported nested forms are
therefore retained by the grammar and rejected at an explicit semantic
boundary.

Lane names and transform spellings, aliases, argument shapes, and domains are
central specification tables. The audio runtime never dispatches source names.
Every sequence also has a deterministic stable definition ID, currently
derived from its resolved name, which carries lane counters through content
edits. Runtime state transfer compares both ID and name, so a hash collision
cannot attach one sequence's state to another. Both programs contain a
pre-sorted transfer index, making an audio-thread program swap linear in the
number of definitions. A future editor-owned ID can replace the derivation if
phase should also survive renames.

`SemanticProgram` contains only the resolved musical graph. `CompiledProgram`
contains one privately constructed `SemanticProgram` and adds checked timing
bounds plus preallocated step, schedule, and state workspaces. The compiler
factory is the only construction path, so an unprepared semantic graph cannot
masquerade as an audio-ready program. Only a `CompiledProgram` may cross the
audio-thread handoff.

## Execution and ownership

The editor owns draft text. Ctrl+`.` compiles the complete draft and makes it
the evaluated document. Ctrl+Enter evaluates the top-level statement
containing the selection or current line against the last successful document.
A line inside `sequence { ... }` therefore updates that complete sequence, but
edits to other definitions remain inactive. `play`, `stop`, and `seed` replace
their corresponding evaluated command state. Failed compilation leaves both
the evaluated document and playing program unchanged.

This also applies while loading a preset or saved document into a running
module: invalid or unsupported source reports an error but never stops or
replaces the last valid active program. During initial construction the factory
program is removed before saved source is evaluated, so a document that has
never compiled remains silent rather than playing unrelated factory notes.

The evaluated state is a typed `Document`, not reconstructed source text.
Ctrl+Enter parses the draft once, merges selected statements by definition or
command identity, and compiles the merged tree directly. An unchanged retained
statement adopts current editor spans only when its source text is identical;
an inactive edited statement keeps its old meaning with invalid spans so it
cannot produce a misleading cursor over text that has not run.

Accepted v1 limitation: Ctrl+Enter still requires the complete editor draft to
be syntactically valid because cpp-peglib first identifies the draft's complete
top-level structure. Semantic errors in unselected statements remain isolated,
but an unfinished delimiter or bracket elsewhere in the draft can block the
selected evaluation.

cpp-peglib parses only on the UI thread. The semantic compiler produces an
immutable program plus mutable workspaces whose sizes are prepared before
publication. An atomic pending/retired handoff swaps programs on clock edges;
the audio thread never parses, allocates, locks, logs, or destroys a program.
Saved language versions which cannot be interpreted must not fall back to an
unrelated default program.

Publication is transactional: all potentially throwing editor-state copies and
status preparation finish before the final atomic exchange. Setting lanes such
as `scale`, `tonic`, `cycle`, and `glide` reject inline pipelines instead of
silently ignoring them; whole-sequence transforms belong after the closing
brace.

Lookahead scheduling retains a dedicated checkpoint at the next incoming-clock
boundary. A valid replacement restores that checkpoint, transfers named
sequence state, preserves the first unscheduled logical onset and currently
held outputs, and regenerates only the future queue. Preparing sub-beat events
beyond the boundary therefore cannot advance the state adopted by a hot swap.

## Musical representation

A pitched item explicitly distinguishes:

- one scale degree or absolute note;
- a literal parenthesized voicing; and
- a semantic jazz chord with preserved root and source symbol.

`D@4` is one named note while `D7@3` is a dominant seventh rooted in octave 3.
`(1 b3 5)@4` and `(C E G)@4` are literal simultaneous voicings; individual
registers such as `(C@3 E@4 G@4)` override the shared register. Literal
voicings bypass future interpretation by default.

`@` always denotes an absolute octave, including signed values such as
`C@-1`. Apostrophes and commas are the composable relative register syntax:
`1'` and `Cm7'` are one octave above the active sequence register, while `1,,`
is two octaves below. Language version 4 establishes this distinction; the
former version-3 interpretation of signed `@` as relative is intentionally not
retained because it made negative absolute octaves ambiguous.

Pipelines compose from left to right. In particular:

```text
base |> modulate_degree 3 |> scale major
base |> scale major |> modulate_degree 3
```

are intentionally different when `base` is not already major. A modulation
captures the scale active at that point in the chain. Definition names share
one namespace regardless of whether they name a sequence or an arrangement.

Runtime output is one aligned Rack polyphonic bundle: Pitch, Gate, Trigger,
Velocity, and Accent use the same voice indices, up to Rack's 16 channels.
One musical onset publishes one editor cursor pulse regardless of chord width.

## Timing and prepared scheduling

Incoming clock edges own integer beats. The measured clock period is used only
to interpolate subdivisions and convert millisecond timing offsets. The
scheduler generates far enough ahead to cover the greatest possible combined
`early` displacement in the compiled program. Pending-event capacity covers
lookahead plus combined beat and millisecond `early`/`late` displacement.

Millisecond capacity is prepared conservatively for clocks up to 1 kHz. This
is far above normal musical clock rates and avoids a small arbitrary event
ceiling. Faster clocks may report prepared-workspace exhaustion rather than
allocate in the audio callback. Event-density arithmetic is overflow-checked;
allocation failure is a compile diagnostic. There is no separate fixed
64-step scheduler guard.

The first external clock interval cannot predict a millisecond-relative early
event because no period has yet been measured. Beat-relative lookahead remains
available immediately; millisecond lookahead becomes exact after the second
clock edge.

`swing RATIO` uses each event's subdivision as its grid. The optional
`swing RATIO SUBDIVISION` form selects an explicit grid in incoming-clock beats,
so `swing .6 1/8` and `swing .6 1/16` can groove independently of the duration
lane. Swing is a beat-grid transform and therefore continues across evaluator
step boundaries rather than being confined to one bracket group.

Pattern repetition, Euclidean step counts, ratchets, timing displacement, and
glide have no small musical constants. The compiler expands and validates them
on the UI thread using checked addressability and event-density arithmetic.
Rack's 16-channel polyphony remains a real output-format boundary.

## Effective pitch and future interpreters

Future input binding and interpretation happen after pattern selection and
active pitch transforms, but before final event realization:

```text
compiled pitched item
  -> lane/conditional transforms
  -> live input binding
  -> effective pitched item
  -> instrument interpreter
  -> prepared note/voicing/sub-event schedule
```

An interpreter must not consume the raw stored `jazzSymbol` as though it were
already transposed. It receives an effective view containing the transformed
root/pitches, chord intervals and slash bass, active scale/tonic, register and
range, previous/current/next item, prior realized voicing, event span,
absolute beat, deterministic seed, and immutable live-performance snapshot.

The planned `input` lane chooses what a controller means. `midi.root` can
transpose either an ordinary note sequence or chord progression from its
declared tonic; `midi.chord` supplies held notes as harmonic material. The
independently cycling `interpret` lane chooses realization, for example
`piano.simple`, `piano.rootless4`, `bass.two`, `bass.walk`, or `arp.up`.
Interpreters may emit simultaneous voices or bounded sub-events within the
current item span, but user callbacks and dynamic allocation remain forbidden
on the audio thread.

MIDI input is normalized into a fixed-size timestamped performance state:
held notes, velocities, most recent note, and pedal/controller state. It never
mutates editor text or the immutable compiled graph. With no live input the
program stays deterministic; with an input event log it is replayable.

### Future signal expressions

The language must later accept small user functions and expressions over MIDI,
CV, and previously computed signals, for example:

```text
cv2 = exp(0.1 * vel) + cv1
cv3 = midi_pitch * cv2
```

cpp-peglib will emit typed call, unary, binary, name-reference, and assignment
nodes. Name resolution classifies sequence/arrangement values, event-rate MIDI
values, and audio/control-rate scalar signals explicitly; it must not infer a
domain by reparsing source text. The semantic compiler checks types, dependency
cycles, available inputs, and bounded state, then lowers the expression graph
to a prepared immutable evaluator representation. Runtime evaluation uses
fixed storage and direct numeric operations only—no parsing, allocation,
dynamic name lookup, or user callback from `process()`.

The input snapshot/binding boundary above remains the source of MIDI values.
CV expressions form a prepared signal stage before live input binding and
instrument interpretation, and their outputs may be referenced by sequence
lanes or other signal expressions. Adding this requires a language-version
bump and new expression AST alternatives, but does not require changing the
current `SemanticProgram`/`CompiledProgram` publication boundary.

## Deferred implementation

The following are designed but not wired yet:

- Rack MIDI input UI and timestamped performance-state adapter;
- executable `input`, `interpret`, and `chords` lane vocabulary;
- contextual voicing/voice-leading and bass/arp interpreters;
- routing for several interpreted parts from one module; and
- CV-derived performance input using the same binding boundary;
- typed MIDI/CV function expressions and their prepared evaluator.

When these arrive, add them as prepared semantic/runtime stages. Do not teach
the chord parser instrument-specific voicings and do not invoke an embedded
language or arranger callback from `process()`.
