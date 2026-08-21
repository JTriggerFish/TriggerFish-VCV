# TfProgSequencer current design

This is the authoritative implementation note for Prog Sequencer. The syntax
study records alternatives and the v1 proposal retains historical context;
this document records the boundaries the current code must preserve.

The authoritative executable language is specified in
[TfProgSequencer-language-grammar.md](TfProgSequencer-language-grammar.md).
The current persisted language version is 5. It replaces the earlier prototype
outright: this build contains one cpp-peglib grammar and one semantic compiler,
with no v3/v4 parser or compatibility branch. The saved version is an explicit
incompatibility marker and a future migration boundary, not runtime dispatch.

## Compiler architecture

The implementation has four explicit representations:

```text
cpp-peglib document AST
  -> typed syntax Document and domain-neutral PatternNode tree
  -> checked SemanticProgram
  -> prepared CompiledProgram
  -> bounded RuntimeEvent stream
```

`PatternNode` owns reusable structure—complete note event, rest, tie,
subdivision, cycle choice, random choice, voicing, slash bass, and typed
suffixes—before semantic lowering assigns runtime meaning.
cpp-peglib productions emit these structural nodes directly; there is no
secondary delimiter or nesting parser. Voicing tones, shared register suffixes,
and slash-bass relationships are likewise syntax nodes rather than text split
again by the compiler. Lane compilers turn the validated leaves into typed
pitches/chords, note-entry semantics, scalars, timing offsets, or CV knots.
Unsupported multi-event nested choice branches are retained by the grammar and
rejected at an explicit semantic boundary rather than flattened ambiguously.

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
Ctrl+Enter first parses the complete draft so unchanged statements can adopt
current source spans. If unrelated draft syntax is incomplete, a line-bounded
fallback searches for the smallest complete statement containing the selection
and parses it with the same cpp-peglib grammar; there is no second syntax
parser. The selected statements are merged by definition or command identity
and the typed tree is compiled directly. An unchanged retained statement adopts
current editor spans only when its source text is identical; an inactive edited
statement keeps its old meaning with invalid spans so it cannot produce a
misleading cursor over text that has not run. A continuation line beginning
with `|>` is always retained with its containing sequence.

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

## Editor rendering

The text display renderer expresses every visual state as one normalized
intensity. Background, ordinary text, comments, selection, status text,
execution flashes, and cursor persistence all pass through the same active
heatmap; none owns a semantic RGB colour. The current map uses uniformly spaced
samples from Matplotlib's canonical magma table. Replacing the single active
table changes the complete display scheme without modifying drawing or decay
logic. Cursor histories remain fixed-size UI-thread state and never affect the
audio thread.

For short moves, the cursor itself is a fractional-pixel phosphor beam; the
destination has only a faint resting marker until the beam arrives. Each lane
retains four timestamped motions, so a new clock does not erase an older trail.
Every beam deposits up to fourteen overlapping additive glow samples across
glyphs and whitespace. Each sample represents an earlier beam timestamp and
decays according to the time since the beam passed it, giving the completed
path spatially coherent exponential persistence.
These are fixed-capacity, UI-thread-only arrays rather than per-frame heap
allocations.

Events without local travel--especially repeated locations--also create a
subtle bloom around their glyph. Eighteen antialiased, sub-pixel-spaced fills
diffuse outward from the rectangular caret's footprint. The tight inner field
retains that shape; successively wider layers become fainter and more rounded
until the outside is only a soft haze, with no stroked rings or contours. This
fixed layer count avoids heap allocation and the unreliable tiny-glyph corner
case in NanoVG's box-gradient shader. Four independent blooms are retained per
lane, so retriggers add energy rather than resetting the previous halo. First
events and deliberate long jumps use the same stationary feedback.

Successful `Ctrl+Enter` and `Ctrl+.` evaluations reuse the diffusion field at
lower opacity. The editor derives one source rectangle per executed visual text
row, so a selection, line, or complete program blooms from its actual occupied
area instead of flashing a panel-sized rectangle.

Envelope intensity is calculated from event timestamps rather than accumulated
per-frame state, so dropped UI frames cannot change brightness. Travel uses 75%
of the observed interval between that lane's pulses, bounded to 35--120 ms; its
afterglow also follows the interval, bounded to 75--300 ms. Thus motion remains
visible across several display frames without turning very fast patterns into
a permanent smear. Repeated positions pulse in place. The saved `Cursor travel`
menu setting switches only spatial progress between constant-speed `Linear` and
eased `Smoothstep`; both use the identical temporal envelope and heatmap.

The sequence/arrangement cursor observes the outer played expression. A named
term is a feedback abstraction boundary: in `song = acid + iv + v`, playback
highlights those three tokens on the `song` row and never leaks an internal
`acid` span from the definitions of `iv` or `v`. Explicit parenthesized groups
retain their inner term spans. With the default duration of one and no explicit
subdivision, the active arrangement term pulses once per incoming-clock beat.

## Musical representation

The `notes` lane is event-first. A pitched token is an ordinary onset; `x`,
`^`/`^^`, and `>` prefix ghost, accent, and slide entry. `_` is a standalone
tie cell, `~` a rest, and duration/Euclidean/ratchet/probability/replication
suffixes belong to the event or group they follow. There is no source-level
articulation lane and no legacy articulation parser. Optional scalar lanes are
overrides or independent polymeters rather than the primary event model.

A pitched item explicitly distinguishes:

- one scale degree or absolute note;
- a literal parenthesized voicing; and
- a semantic jazz chord with preserved root and source symbol.

Roman symbols `I` through `VII` retain a scale-relative root plus explicit
chromatic chord intervals. Uppercase implies major and lowercase minor, so the
same prepared representation can stay in key under `shift_degree`, modulate as
a complete event under `modulate_degree`, and later be handed to an instrument
interpreter without losing harmonic intent.

`D@4` is one named note while `D7@3` is a dominant seventh rooted in octave 3.
`(1 b3 5)@4` and `(C E G)@4` are literal simultaneous voicings; individual
registers such as `(C@3 E@4 G@4)` override the shared register. Literal
voicings bypass future interpretation by default.

`@` always denotes an absolute octave, including signed values such as
`C@-1`. Apostrophes and commas are the composable relative register syntax:
`1'` and `Cm7'` are one octave above the active sequence register, while `1,,`
is two octaves below. Version 5 retains this unambiguous distinction while
moving articulation, duration, replication, probability, ratchet, and
Euclidean structure onto complete note events.

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
CV1 and CV2 are monophonic sequence controls. One musical onset publishes one
editor cursor pulse regardless of chord width.

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

The `offset` lane is the patternable form of early/late displacement. Negative
values are early, positive values late, and both beat fractions and `ms` are
typed values. Free offset and CV lanes derive phase from absolute incoming-clock
score time; `rate R` scales only that lane. Thus subdivisions sample the same
lane more often without accelerating it. Edge-aligned ellipsis and numeric
`rate` are intentionally distinct modes and cannot be combined.

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

Pattern-driven arpeggiation is a first-class instance of this binding rather
than a MIDI special case. Conceptually, `midi.in |> arp arp1_seq` combines a
timestamped performance source with a compiled event template. The canonical
sequence form keeps timing ownership visible:

```text
arp1 = sequence {
  notes 1 2 3 4 3 2
  input midi.chord
  interpret arp.index
}
```

`arp.index` treats the written degree as an index into the ordered held-note
set, so rests, ties, durations, probability, subdivisions, and ratchets remain
properties of `arp1`. `arp.up`/`arp.down` are different interpreters: they may
emit several bounded sub-events inside one source event. A future source-first
expression may be accepted as typed sugar, but bare `midi.in |> arp1_seq` is
not used because it would make an ordinary sequence value implicitly callable.
Both surfaces lower to a node typed as `(PerformanceSource, EventPattern) ->
PitchedEventStream` and read the immutable MIDI snapshot at scheduled event
time.

MIDI input is normalized into a fixed-size timestamped performance state:
held notes, velocities, most recent note, and pedal/controller state. It never
mutates editor text or the immutable compiled graph. With no live input the
program stays deterministic; with an input event log it is replayable.

### Sequenced CV lanes

The grammar recognizes positive `cvN` lane names and the module implements
`cv1` and `cv2`; unavailable output indices are semantic diagnostics. These are
time-indexed scalar patterns, not pitch events and not the future top-level
signal assignments below. They use the same polymetric and edge-aligned forms
as other numerical lanes. Step output is the default; `interp linear`,
`interp smooth`, and `interp power P` create curve segments between explicit
score-time points. Every logical event, including a rest or tie, samples the
lane; ratchets share that sample. Edge-aligned CV is stepped in the current
direct interpreter, because interpolating it correctly requires prepared times
for every following structural boundary.

Compilation lowers every curve to typed scalar knots with source spans and
enough cyclic lookahead to identify the next deterministic point. The audio
thread evaluates only the active segment into preallocated output state. Lane
phase derives from absolute incoming-clock score time; subdivisions may sample
more often but do not advance it, and numeric `rate` is explicit. If an event
span crosses several CV knots, the current interpreter reaches its pending
target and holds until the next logical event; a future prepared control-event
stream will traverse the intermediate segments independently. See the
normative interpolation rules in `TfProgSequencer-language-grammar.md`.

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

- prepared execution of multi-event and recursively nested random/alternate
  branches (atomic choices already execute);
- Rack MIDI input UI and timestamped performance-state adapter;
- executable `input`, `interpret`, and `chords` lane vocabulary;
- contextual voicing/voice-leading and bass/arp interpreters;
- routing for several interpreted parts from one module; and
- CV-derived performance input using the same binding boundary;
- typed MIDI/CV function expressions and their prepared evaluator.

When these arrive, add them as prepared semantic/runtime stages. Do not teach
the chord parser instrument-specific voicings and do not invoke an embedded
language or arranger callback from `process()`.
