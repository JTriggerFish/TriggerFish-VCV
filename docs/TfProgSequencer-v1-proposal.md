# TriggerFish Prog Sequencer: v1 proposal

> **Historical event syntax:** The executable v5 event and numerical
> lane grammar superseding the event examples below is recorded in
> [TfProgSequencer-language-grammar.md](TfProgSequencer-language-grammar.md).

> **Design status:** This is the original architecture proposal. Option A in
> [TfProgSequencer-syntax-options.md](TfProgSequencer-syntax-options.md) led to
> the now-replaced first prototype. The implementation uses a
> version-pinned, vendored cpp-peglib PEG front end, a separate semantic
> compiler, and a separate allocation-free playback runtime. The Lua material
> below is historical, not the current direction. The implementation has also
> moved beyond this proposal's monophonic scope with `note@octave` registers,
> an independent octave lane, explicit and jazz/slash chords, and 16-channel
> Rack polyphonic outputs.
> The authoritative implemented boundaries and deferred interpreter/MIDI work
> are maintained in
> [TfProgSequencer-current-design.md](TfProgSequencer-current-design.md).

## Current architecture addendum: pitch and live interpretation

Every compiled pitched item distinguishes a single note, an explicit voicing,
or a semantic jazz chord. Jazz symbols retain their root and harmonic identity;
their current close-position voices are a default realization, not a loss of
meaning. This supports future sequenceable interpreters such as `bass.two`,
`bass.walk`, `piano.simple`, `piano.rootless4`, and arpeggiators without
changing the note or chord language.

Live controller input belongs on the real-time side of the same boundary. A
MIDI adapter should normalize timestamped messages into a fixed-size
performance state containing held notes, velocities, sustain/pedal state, and
the most recent note. At each scheduled pitched event, the interpreter is
given an immutable view containing:

- previous, current, and next pitched items;
- event span and absolute beat phase;
- active scale, tonic, register/range, and deterministic seed;
- current live performance state; and
- the prior realization when voice leading needs continuity.

Input-binding policy and realization policy remain separate, sequenceable
lanes. For example, `midi.root` can transpose an ordinary melody or chord
progression from its declared tonic to the played root. Another binding can
offer the complete held MIDI chord. An interpreter may then produce direct
notes, a block or rootless voicing, a two-feel or walking line, or an arpeggio.
It can schedule several prepared sub-events within the current event span, but
no parser, dynamic allocation, lock, or user callback is allowed in the audio
callback.

MIDI activity never mutates editor text or the immutable compiled graph. The
same source therefore remains deterministic when no live input is present and
replayable given the same timestamped input event stream. Direct MIDI UI/port
wiring and the final policy vocabulary are intentionally deferred; the graph,
scheduler, and semantic-pitch boundary must not assume they are absent.

## Original recommendation (superseded)

Build a monophonic, externally clocked sequencer with two language layers:

1. A small pattern notation for the musical material. It makes notes, rests,
   ties, durations, repetitions, ratchets, slides, accents, subdivisions, and
   variations quick to type.
2. An embedded Lua host for naming patterns and combining or transforming them.
   Lua is used to construct an immutable pattern graph; it is never run on the
   audio thread.

This keeps the common case close to Tidal's mini-notation while retaining a
real language when a piece grows beyond a single pattern. It also gives the DSP
side a bounded, deterministic representation that can be exchanged safely
while Rack is processing audio.

The first editor should be intentionally plain. It needs selection, multiline
editing, horizontal and vertical scrolling, execution shortcuts, and a status
area. It does not need syntax highlighting. Bracket matching can follow later,
but is not required for the first usable build.

## Product character

The module should feel like an instrument, not an IDE:

- A four-note idea should take one short line.
- Editing valid code should never interrupt the clock.
- Invalid code should never stop the last valid pattern.
- Randomness should be easy to invite but repeatable when it produces something
  worth keeping.
- Scale-degree notation should make transposition and mutation musically useful
  by default.
- Explicit pitch, velocity, duration, gate, and slide values should always be
  available as an escape hatch.
- The language should be compositional: a pattern made by an operation can be
  passed to any other operation.

## Proposed module

### Panel and I/O

A 22 HP panel is a reasonable starting point: a 16 HP text display derived from
Rack's Notes module, with a 6 HP strip for ports and minimal transport feedback.
This preserves useful editor height instead of placing every jack below it.

Inputs:

- `CLOCK`: rising edges define beat boundaries.
- `RESET`: cancels pending sub-events, drops Gate and Trig, and makes the next
  Clock edge beat zero. Clock is ignored for 1 ms after Reset, following Rack's
  recommendation for cable-delay tolerance.
- `RUN`: optional for v1. When low it stops scheduling events and drops Gate;
  pitch and velocity retain their last values.

Outputs:

- `V/OCT`: up to 16 Rack polyphonic pitch channels, with C@4 = 0 V.
- `GATE`: 10 V while the current note is held.
- `TRIG`: 10 V for 1 ms at each articulated note onset.
- `VEL`: 0-10 V velocity. It is held at the last note value during a rest.
- `ACC`: 0-10 V sparse accent, high only while the accented Gate is high.

The pitch output should also retain its last note through a rest. This avoids
unnecessary pitch jumps while Gate is low.

No pitch knob, scale knob, or step controls are needed in v1. The text is the
source of truth. A context menu can provide trigger duration, live-update mode,
and light/dark editor appearance without consuming panel space.

### Editor layout

Use one plain monospaced editor viewport and reserve its bottom 2-3 lines for a
status panel. The status panel can show:

- `READY`, `QUEUED`, `PLAYING`, `WAITING FOR CLOCK`, or `STOPPED`;
- current beat and pattern length;
- the last successful evaluation time; and
- one concise diagnostic, including line and column.

For example:

```text
ERROR 8:17  expected ')' after every(...)
last good program is still playing
```

The first cut should have no syntax highlighting. An error underline and a
matching-bracket box are diagnostics/navigation rather than full highlighting,
but both can wait. Scrolling the caret into view is not optional: Rack's Notes
field clips text and does not provide the viewport behavior needed for code.

### Execution model

- Primary modifier + `.` evaluates the entire document in a fresh sandbox.
- Primary modifier + `Enter` evaluates the top-level statement containing the
  selection in the context of the last successful document. With no selection,
  it evaluates the statement containing the current line. A lane line updates
  its containing sequence; unrelated draft edits remain inactive.
- Primary modifier + `Shift` + `Enter` evaluates and explicitly restarts the
  pattern at the next Clock edge.
- A successful evaluation is queued atomically for the next Clock edge.
- Normal evaluation preserves the absolute beat phase. Reset or explicit
  restart begins at beat zero.
- A failed evaluation reports an error and leaves the last valid program
  running unchanged.
- `stop` is a transport command; like `play`, it takes effect when the queued
  program activates on the next Clock edge.

"Primary modifier" means Ctrl on Windows/Linux and the Rack-equivalent command
modifier on macOS.

Whole-document evaluation deliberately starts from a clean environment so a
saved patch is reproducible. Selection evaluation uses the current successful
environment, which makes incremental live coding practical.

## Language proposal

### Original host-language hypothesis

Use vendored, version-pinned Lua and expose only a small pattern-building API.
Lua gives the advanced user variables, functions, loops, tables, and ordinary
arithmetic without requiring TriggerFish to invent and maintain a complete
general-purpose language.

The sandbox should omit file, process, package, network, and debug access. Apply
an instruction budget, memory budget, source-size limit, AST-depth limit, and
maximum event-density limit. Patterns are immutable values.

Lua is preferable here to the main alternatives:

| Option | Strength | Reason not to choose it for v1 |
| --- | --- | --- |
| Mini notation only | smallest and safest runtime | creates a low ceiling and eventually grows into an accidental language |
| Entirely custom language/VM | complete control of surface syntax | parser, runtime, tooling, and compatibility all become TriggerFish responsibilities |
| Embedded JavaScript | closest conceptual match to Gibber/Strudel | a larger and more complex host than this monophonic module needs |
| Embedded Lua plus mini notation | small C embedding API and a genuine language behind the shorthand | chosen; the two layers and sandbox need a clear boundary |

A complete program could begin like this:

```lua
set {
  tonic = "D@4",
  scale = "dorian",
  velocity = .72,
  accent = .88,
  strong_accent = 1.0,
  gate = .80,
  slide = .25,
  seed = 23,
}

riff = p [[1 2 b3+ [4 5] 5@2 ~ <b7 8>]]

lead = riff
  :every(4, rev)
  :sometimes(.25, shift_degree(2))

play(lead)
```

Lua's `p [[...]]` form avoids quote escaping, while `p "..."` remains convenient
for short patterns.

### Pitch atoms

Scale degrees are one-based and can continue across octaves:

| Text | Meaning in the active scale |
| --- | --- |
| `1 2 3 4 5 6 7` | degrees in the tonic octave |
| `8 9 10` | degrees one octave above |
| `0 -1` | degrees below the tonic octave |
| `b3 #4` | chromatic alteration of a degree |
| `C@4 F#@5 Bb@3` | absolute equal-tempered note names |

Continuing degree numbers avoid adding a second octave punctuation system and
make transformations predictable. `shift_degree(7)` moves a pattern by one
scale octave; `transpose_semitone(12)` moves by one chromatic octave.

The active tonic is translated to Rack pitch with C@4 = 0 V. Named notes bypass
the scale. A later tuning extension can let a scale contain arbitrary cent or
ratio offsets without changing the event notation.

### First-class event syntax

These constructs are part of the mini parser, not library functions:

| Syntax | Meaning |
| --- | --- |
| `1 3 5` | three notes, one input-clock beat each |
| `~` | one beat of silence; Gate falls, pitch is retained |
| `_` | tie the prior note for one more beat; no Gate or Trig retrigger |
| `1@3` | hold for three beats; equivalent to `1 _ _` |
| `1!3` | repeat over three beats with three separate articulations |
| `1*3` | ratchet three articulations inside one beat |
| `[1 2 3]` | subdivide one beat into three equal events |
| `1 -> 3` | slide into `3`; Gate remains high and `3` does not retrigger |
| `1 -> 3{slide=.5}` | use a half-beat pitch glide for this transition |
| `3+` / `3++` | accent / strong accent velocity |
| `3{vel=.63}` | exact velocity, normalized to the 0-10 V output |
| `3{gate=.4}` | exact gate fraction of this event span |
| `3{dur=3/2}` | exact duration in incoming-clock beats |

The important distinctions are:

- A tie extends the same event and keeps Gate high.
- A long duration is the compact form of one or more ties.
- A repetition creates new note onsets and therefore new Gate/Trig
  articulations.
- A ratchet preserves the event's overall duration but puts repeated onsets
  inside it.
- A slide connects two different pitched events, keeps Gate high across the
  boundary, suppresses the target Trig, and slews V/OCT to the target.
- A rest consumes time and explicitly makes Gate low.

The slide amount is a duration in beats, capped by the target event span. A
zero slide is an immediate legato pitch change. Velocity changes at the target
boundary even when the target is reached by a slide.

Durations may also be applied to rests (`~@4`). A tie without a preceding note
is a compile error, and a slide whose source or destination is a rest is a
compile error. These diagnostics are preferable to silently surprising output.

### Pattern structure and variation

The following subset gives v1 a large musical range without requiring chords or
drum lanes:

| Syntax/API | Meaning |
| --- | --- |
| `[a b]` | subdivision inside one beat |
| `<a b c>` | use the next alternative on each pattern pass |
| `[a|b|c]` | deterministic random choice |
| `a?` / `a?.2` | play with probability .5 / .2 |
| `a(3,8)` / `a(3,8,1)` | Euclidean 3-in-8 rhythm, optionally rotated; consumes 8 beats |
| `cat(a, b)` | concatenate patterns |
| `repeatn(4, a)` | concatenate four copies |
| `a:rev()` | reverse event order |
| `a:rotate(2)` | rotate by two beats |
| `a:early(1/4)` / `a:late(1/4)` | time offset |
| `a:fast(2)` / `a:slow(2)` | time scaling |
| `a:shift_degree(2)` | diatonic in-key remapping |
| `a:transpose_semitone(1)` | chromatic transposition |
| `a:every(4, rev)` | apply a transform every fourth pass |
| `a:sometimes(.25, rev)` | probabilistically apply a transform |
| `a:quantize("minor pentatonic")` | quantize absolute/generated pitch |

Random decisions should be derived from `(seed, absolute beat, event path)`.
They then survive edits and transport jumps without depending on how many
random calls happened previously. Changing `seed` creates a new variation;
keeping it preserves a happy accident in a saved patch.

Top-level commas should be reserved for future simultaneous patterns/chords.
In v1 they should produce a clear "polyphony is not supported yet" diagnostic
rather than acquiring another meaning.

### More examples

A straightforward sequence:

```lua
set { tonic="C@4", scale="minor", gate=.75 }
play(p [[1 3 5 8]])
```

Tie, rest, duration, repeated articulation, ratchet, and slide:

```lua
play(p [[1 _ ~ 3@2 5!3 8*4 5 -> 1]])
```

Human-readable emphasis with one exact exception:

```lua
play(p [[1 2+ b3 4++ 5{vel=.57} ~ 8]])
```

A variation that remains deterministic after saving:

```lua
motif = p [[1 [2 b3] 5 <4 6> ~ 8]]
play(motif
  :every(4, rotate(1))
  :sometimes(.20, shift_degree(2)))
```

Reusable ordinary Lua removes the ceiling without complicating the common
notation:

```lua
function rise(notes, degrees)
  local result = notes
  for i = 1, degrees do
    result = cat(result, notes:shift_degree(i))
  end
  return result
end

cell = p [[1 2 b3 5]]
play(rise(cell, 3))
```

## Clock and scheduling semantics

Each `CLOCK` rising edge is one beat and is authoritative. Plain top-level
events therefore behave exactly like a conventional sequencer: `1 2 3 4`
produces one note on each of four successive edges.

Subdivisions and ratchets require events between external edges. The engine
tracks the recent clock interval and schedules their rational positions at
sample accuracy. It should:

- use a short robust average/median rather than chase single-edge jitter;
- re-anchor phase on every external edge so it cannot drift away;
- cancel stale sub-events when Reset arrives;
- wait for a usable period estimate before emitting off-edge events;
- report `WAITING FOR CLOCK` or `CLOCK UNSTABLE` in the status area; and
- stop sub-events and lower Gate after a conservative clock-loss timeout.

This model deliberately differs from Tidal's fixed-length cycle: the outer
sequence is measured in incoming-clock beats, while brackets compress their
contents into one beat. It fits modular clock expectations better. Pattern
passes still provide a stable unit for `<...>`, `every`, and seeded variation.

Gate length is based on the event's scheduled span, not merely on the last
external period. Ties and slides override the normal gate gap. Trigger length is
independent and defaults to Rack's recommended 1 ms.

## Runtime architecture

The implementation should have four distinct layers:

```text
plain editor text
       |
       v
Lua sandbox + mini-notation parser        UI/worker side only
       |
       v
immutable typed pattern graph
       |
       v
lock-free pending/retired graph mailboxes
       |
       v
clock scheduler -> V/OCT, GATE, TRIG, VEL, ACC audio thread
```

### Compilation side

Lua evaluation builds C++ pattern nodes through a small direct C API binding.
The mini parser produces source-located nodes for events and transforms. Before
publication, validation computes bounds such as maximum nesting, event density,
and scheduled work per audio block.

Do not expose a Lua callback that is invoked from `process()`. Advanced Lua code
may build graphs, loop over data, and define graph-building functions, but the
published result must be the typed declarative graph.

### Audio side

The audio thread owns transport, clock estimation, gate/trigger timers, slide
state, and the currently active graph. Graph queries must be bounded and must
not allocate, lock, log, parse, or call Lua.

The implemented scheduler prepares queue capacity from event density plus the
largest combined beat/millisecond early and late offsets. Its lookahead uses
the compiled maximum early displacement rather than assuming one beat. The
millisecond bound is conservatively prepared for clocks up to 1 kHz; faster
clocks report workspace exhaustion instead of allocating in `process()`.

A small single-producer/single-consumer mailbox publishes a candidate graph.
The engine swaps it on a Clock edge and returns the retired graph through a
second mailbox so destruction occurs away from the audio thread. This avoids a
mutex and avoids the last smart-pointer destructor running in `process()`.

### Persistence

Save at least:

- source text;
- editor scroll/cursor state if useful;
- current seed and context-menu preferences; and
- a language/schema version.

On patch load, evaluate the full document in a fresh sandbox. The version field
allows a future parser to preserve or migrate old semantics. Do not silently
reinterpret a saved source file after a breaking language change.

## Relationship to existing work

### Retain

- From Tidal/Strudel: mini-notation, nested subdivision, alternation, Euclidean
  and probability-friendly patterns, immutable transforms, and patterns as
  composable values rather than flattened arrays.
- From Gibber: immediate selection evaluation, named reusable objects, and
  method-style transformation chains.
- From Teletype: external events as the centre of the instrument, short commands,
  deterministic state, and clear live/edit feedback.
- From crow: Lua as a small embedded host paired with a purpose-built musical
  vocabulary.
- From VCV Notes: the native Rack LED display appearance and basic text editing
  behavior.
- From Spellbook: on-module plain-text sequencing and explicit modular voltage
  concepts.

### Avoid in v1

- Tidal's full polyphonic/control-pattern surface;
- a general DSP scripting API like VCV Prototype;
- a tracker table with one row per clock;
- a home-grown general-purpose interpreter;
- executing user language code in the audio callback;
- syntax highlighting, autocomplete, a minimap, or other IDE features;
- chords, drums, multiple lanes, or polyphonic output; and
- filesystem/network/package access from scripts.

## Fit with this repository

The plugin currently targets Rack SDK 2.6.6 and C++17, uses the Rack SDK
Makefile for distributable builds, and separately builds testable DSP code with
CMake. Lua would be the repository's first embedded language dependency, so it
should be vendored and built statically for Windows x64, Linux x64, macOS x64,
and macOS arm64.

The implemented code is split at the parser/compiler/real-time boundaries:

```text
src/TfProgSequencer.cpp                 Rack module, persistence, and editor UI
src/tfseq.hpp                           immutable graph and public runtime API
src/tfseq_parser.*                      PEG grammar and typed syntax records
src/tfseq_compiler.cpp                  musical validation and graph lowering
src/tfseq_runtime.cpp                   clock-time, allocation-free playback
tests/text_sequencer_tests.cpp          parser, compiler, transform, and timing tests
vendor/cpp-peglib/                      pinned single-header parser and MIT license
res-src/TfProgSequencer*.svg            editable panels
res/TfProgSequencer*.svg                outlined Rack runtime panels
```

The Rack Makefile picks up the top-level translation units. The standalone
CMake target lists the parser, compiler, and runtime explicitly. cpp-peglib is
included by only the parser translation unit and never reaches the Rack module
or audio-runtime interfaces.

## Delivery plan

### Phase 0: risk spikes

1. Vendor and build a pinned Lua release on every CI architecture.
2. Prove that a custom `LedDisplayTextField` subclass can scroll, retain
   selection behavior, intercept the two execution shortcuts, and report a
   line/column diagnostic.
3. Prototype the graph mailbox and clock-subdivision scheduler with no Lua or UI.

These are throwaway-quality spikes intended to settle build, editor, and timing
risk before the language surface hardens.

### Phase 1: headless musical core

Implement and test pitch parsing, rests, ties, duration, repeats, ratchets,
slides, accents, exact attributes, grouping, alternation, probability, and the
core transformations. Include golden event-timeline tests at multiple sample
rates and clock tempos.

### Phase 2: module and editor

Add Clock/Reset and the five outputs, text persistence, the plain scrolling
editor, status feedback, whole-buffer/selection evaluation, and transactional
hot swap. Supply a generated smoke-test patch.

### Phase 3: hardening

Add sandbox limits, fuzz the mini parser, test hostile and very dense programs,
measure worst-case `process()` time, test clock loss/jitter/reset races, and run
the plugin build on all four supported architectures.

### Phase 4: restrained polish

Add bracket matching only if it proves useful. Improve diagnostics, add a small
built-in example/preset set, and document the language. Do not expand into
chords or drums until the monophonic timing and live-edit experience are proven.

## Acceptance criteria for v1

- A new user can type and run `play(p [[1 3 5 8]])` without opening a manual.
- A normal note, rest, tie, long duration, repeated note, ratchet, and slide have
  observably distinct Gate/Trig behavior.
- Every successful edit becomes active atomically on a Clock edge.
- A syntax/runtime error leaves the prior sequence running and identifies its
  line and column.
- The audio thread allocates no memory, takes no lock, and executes no Lua.
- Reset/Clock behavior follows Rack's 1 ms reset-ignore guidance.
- V/OCT, Gate, Trig, Velocity, and Accent follow Rack voltage conventions.
- The patch restores its source and behavior after save/reload.
- The module builds on Windows x64, Linux x64, macOS x64, and macOS arm64.
- Parser, pattern, and scheduler behavior is covered by deterministic headless
  tests.

## Decisions to validate before implementation

The proposal makes defaults so work can start, but three choices deserve a
short playable prototype rather than prolonged syntax debate:

1. Whether `@`, `!`, and `*` feel sufficiently distinct for duration,
   repetition, and ratcheting.
2. Whether slide is clearest as infix `1 -> 3` or as a modifier on the target
   note.
3. Whether selection evaluation should preserve absolute transport phase by
   default or restart the changed pattern.

Everything else can evolve compatibly behind the language-version field.

## References

- Tidal Cycles mini-notation:
  https://tidalcycles.org/docs/reference/mini_notation/
- Tidal's cyclical time model:
  https://tidalcycles.org/docs/reference/cycles/
- Strudel time modifiers and pattern effects:
  https://strudel.cc/learn/time-modifiers/
  https://strudel.cc/workshop/pattern-effects/
- Gibber live-coding environment:
  https://github.com/gibber-cc/gibber
- monome Teletype:
  https://monome.org/docs/teletype/
- monome crow and its Lua environment:
  https://monome.org/docs/crow/
- VCV Prototype scripting host:
  https://vcvrack.com/Prototype
- VCV Spellbook/RhythML:
  https://library.vcvrack.com/TMT/Spellbook
- VCV voltage standards:
  https://vcvrack.com/manual/VoltageStandards
- Lua reference manuals and C API:
  https://www.lua.org/manual/
