# TriggerFish Prog Sequencer: syntax options

> **Current syntax decision:** The executable v5 prefix/suffix event
> grammar and edge-aligned numerical-lane ellipsis are recorded in
> [TfProgSequencer-language-grammar.md](TfProgSequencer-language-grammar.md).
> The alternatives below remain design history rather than the current syntax
> decision.

> **Historical implementation note:** Option A became the first playable
> prototype, including a separate `articulation` lane. Version 5 has replaced
> that source language outright; use the authoritative grammar above for all
> executable syntax. The remainder of this file records the alternatives that
> led to it.

## Why revisit the first proposal

Attaching `+`, `++`, velocity, gate, and slide metadata directly to notes is
useful for one-off exceptions, but it should not be the only model. Musical
parameters often have independent phrase lengths:

```text
notes:    1 2 b3 4 5 6 b7 8       (8 steps)
accent:   + . .                    (3 steps)
```

If the accent stream advances once per note, the combined phrase returns to its
starting state after 24 note attacks. That displacement is a central feature,
not an edge case.

The first decision should therefore be the event/pattern model and the shape of
ordinary sequences. The embedded implementation language can be chosen after
that. Lua is one possible backend, but it should not dictate the notation.

## What established systems suggest

### FoxDot and SuperCollider: independently cycling event attributes

FoxDot players accept patterns for pitch, duration, amplitude, sustain, and
other attributes. A short attribute pattern cycles against the others. A
typical FoxDot expression has the conceptual shape:

```python
p1 >> pluck([0, 2, 4, 6, 7], dur=[1, 1/2], amp=[1, .6, .6])
```

SuperCollider's `Pbind` has the same important idea: an event is assembled from
named streams such as `degree`, `dur`, `legato`, and `amp`:

```supercollider
Pbind(
  \degree, Pseq([0, 1, 2, 3, 4, 5, 6, 7], inf),
  \dur,    Pseq([1, 0.5, 0.5], inf),
  \amp,    Pseq([1, 0.6, 0.6], inf)
)
```

This is the closest precedent for the desired displacement behavior. The
TriggerFish version should be much less syntactically heavy than `Pbind`.

### Gibber and Sardine: a normal imperative language plus pattern arguments

Gibber associates note and duration sequences with object properties and allows
their pattern data to be transformed live. Sardine similarly accepts pattern
strings in ordinary Python function arguments, including parameter arguments.
Both demonstrate that a conventional, imperative host can coexist with a small
musical pattern language.

### Tidal and Strudel: everything is a pattern

Tidal control values such as gain are themselves patterns and are combined with
note/sound patterns. Strudel expresses the same idea with method chains such as
`.velocity("...")`. Their transformation vocabulary and mini-notation are worth
borrowing.

Their normal mental model is nevertheless cycle-based: several values are
distributed over a common cycle and combination operators decide which
pattern's event structure wins. Tidal polymeters and Strudel's newer stepwise
functions can express differing metres, but a simple independent per-note
cursor should be more direct for this module.

### Sonic Pi: explicit independent cursors

Sonic Pi rings and named `tick` counters make independence completely explicit.
One cursor may advance notes while another advances durations or timbre. This is
excellent as an advanced imperative model, although too verbose for every basic
sequence.

### Alda: readable event notation

Alda is a useful reference for readable notes, rests, lengths, and ties. Its
score-oriented syntax is less suited to indefinitely cycling control streams,
but it is a reminder that common musical actions should not require a nest of
functions.

## Shared semantics for all options

The following behavior is proposed independently of surface syntax.

### Evaluation granularity

Ctrl+`.` evaluates the complete draft. Ctrl+Enter evaluates the top-level
statement containing the selection or caret against the last successful
document. Selecting a lane inside `sequence { ... }` therefore updates that
complete sequence, while edits to unrelated definitions remain inactive.
This statement boundary keeps incremental evaluation deterministic without
turning individual lane fragments into a second incomplete grammar.

### Streams

A sequence may have these named streams:

| Stream | Advances | Default |
| --- | --- | --- |
| `notes` | on a new pitched event | required |
| `octave` | on a new pitched event or chord | tonic octave |
| `duration` | on every scheduled slot | `1` clock beat |
| `accent` | on a new pitched event | normal velocity |
| `velocity` | on a new pitched event | configured base velocity |
| `gate` | on a new pitched event | configured base gate fraction |
| `slide` | on a new pitched event | no slide |
| `ratchet` | on a new pitched event | one articulation |

Every stream loops independently. With eight notes and three accent values, the
accent phase moves against the melody and the combined state repeats after 24
note attacks. Adding a five-value duration stream makes the total state repeat
after 120 scheduled attacks.

### Sparse modifier notation

Modifier streams should use `.` for "use the default/no articulation":

```text
accent    + . .       accent, normal, normal
accent    ++ . + .    strong, normal, accent, normal
slide     . . > .     slide into every third note of four
ratchet   . . 3       ratchet every third note three times
gate      . . .4      use a 40% gate every third note
velocity  . . .63     use exact velocity .63 every third note
```

This allows articulation to be written only where it matters. Inline note
modifiers such as `5+` can remain as a local override, but a separate stream is
the natural choice for repeating articulation.

`>` is shown as a slide *into* the current note. This works better in a separate
stream than marking the source note, because the stream describes how each new
event is entered.

### Rests and ties

The note stream may still contain first-class time events:

```text
notes  1 2 _ ~ b3 >4 5
```

- `_` ties the previous note and does not advance note-onset modifier streams.
- `~` creates silence and does not advance note-onset modifier streams.
- `>4` consumes the next pitched event as a slide target and suppresses its
  trigger.
- The duration stream advances for pitches, ties, and rests because all consume
  scheduled time.

The implemented `articulation` stream keeps pitch data free of `_`, `~`, and
`>` when desired, and supports independently cycling ties, rests, slides,
ratchets, probabilities, subdivisions, and Euclidean patterns.

### Registers and chords

The absolute register marker follows the note: `D@4` is the individual note
D4, while `D7` is the conventional dominant-seventh symbol and `D7@3` places
that chord's root in octave 3. `@` is always absolute, including a signed
register such as `C@-1`. Trailing apostrophes and commas are relative to the
active `octave` lane: `1'` and `1''` raise by one and two octaves; `1,` and
`1,,` lower by one and two. They also apply to jazz chords and explicit
voicings. Parentheses mean simultaneity, not choice:

```text
octave 3 3 4
notes (1 b3 5) Cm7 D7 Bbm7b5 / D
notes (C E G)@4 (C@3 E@4 G@4) Cm7' (1 b3 5),
```

Inside an explicit parenthesized voicing, bare note names are individual
pitches, so one trailing register can be shared by the chord. Outside
parentheses, the octave marker retains the crucial `D@4` note versus `D7`
chord distinction.

`[a|b]` remains deterministic random choice and `<a b>` remains pass-by-pass
alternation, so neither construct clashes with chord voicing. Pitch, Gate,
Trigger, Velocity, and Accent use aligned Rack polyphonic channels; slash bass
is channel 1. Rack's 16-channel cable capacity is the only chord-width limit.

Jazz symbols remain semantic chord values in the compiled graph even though
the current default realization is close position. This reserves a small,
orthogonal `interpret` lane and transform rather than baking arranging policy
into chord parsing:

```text
changes = sequence {
  chords     Cm7 F7 Bbmaj7 G7
  duration   4 4 4 4
  interpret  piano.rootless4 piano.simple
}

bass = changes |> interpret bass.walk
keys = changes |> interpret piano.rootless4
```

`interpret` and the `chords` readability alias are reserved directions, not
yet executable syntax. Like every other lane, `interpret` will loop
independently and advance once per pitched event. An interpreter receives the
previous, current, and next pitched items, their duration, register/range, and
deterministic seed. A pitched item may be a single note, an explicit voicing,
or a semantic jazz chord. The interpreter may emit either one simultaneous
voicing or a timed pattern inside the event span. This covers simple and
rootless piano voicings, as well as two-feel or walking bass, without making
core pitch notation instrument-specific. Parenthesized chords remain explicit
voicings and bypass interpretation by default.

#### Live performance input

The same interpreter context is designed to accept a timestamped snapshot of
live controller input. MIDI does not modify source text or the immutable
compiled graph. The audio side maintains held notes, note velocities,
pedal state, and the most recent note, then supplies that state alongside the
currently playing note or chord whenever an event is interpreted.

One possible compact surface is an independently cycling `input` lane:

```text
line = sequence {
  tonic      C
  scale      dorian
  notes      1 3 5 b7 6 5 3 2
  input      midi.root
}

live = sequence {
  chords     Cm7 F7 Bbmaj7 G7
  duration   4 4 4 4
  input      midi.root midi.chord
  interpret  piano.rootless4 arp.up
}
```

`midi.root` transposes the complete written pitch material from its declared
tonic to the played root. It therefore works for ordinary note sequences as
well as chords, including absolute note literals when that policy is selected.
`midi.chord` offers the held note set to the interpreter, which may combine it
with or replace the written item according to its named policy. `arp.up` may
emit a timed note pattern inside the current event span. Both lanes cycle
independently, just like `octave`, `accent`, and `duration`.

These names are reserved design examples, not yet executable syntax. The
important boundary is fixed: input binding chooses what the controller means;
interpretation chooses how the resulting pitch material is realized. This
also leaves room for CV-derived input later without coupling MIDI handling to
the parser.

The interpreter must receive the effective pitched item after active
transposition/modulation and input binding. The raw stored jazz symbol records
source intent; it is not itself the final transformed chord. The complete
ordering and real-time boundary are recorded in
[TfProgSequencer-current-design.md](TfProgSequencer-current-design.md).

### Line comments and quick truncation

`//` ignores the rest of the current line. It works on an otherwise empty line,
after a lane, and after top-level commands. This makes it quick to audition a
shorter pattern without deleting the remainder:

```text
notes 1 2 3 // 4 5 6 7 8
// accent + . .
play riff // optional reminder
```

Here the note lane has three active entries and therefore wraps after `3`.
`#` is deliberately not a comment marker because it already denotes sharp
scale degrees such as `#4`.

### Finite sections

Independently cycling streams are naturally infinite. Arrangement operations
therefore need an explicit finite view:

```text
riff.take(16)       first 16 scheduled slots
riff.once()         one traversal of the primary note stream
```

Concatenation accepts finite sections. `play(riff)` may play an unbounded riff,
while `play(song)` may play and loop a finite arrangement.

## Option A: musical lane block

This is a small purpose-built surface syntax inspired by FoxDot players and
SuperCollider event bindings, but stripped of host-language ceremony.

```text
riff = sequence {
  notes     1 2 b3 4 5 6 b7 8
  accent    + . .
  duration  1 1/2 1/2 1 2
  slide     . . . >
  ratchet   . . 3
}

play riff
```

Only the contents of a lane use mini-notation. The surrounding language remains
ordinary statements, assignments, calls, and blocks, while musical
transformations use concise functional pipelines. A pipeline on a lane affects
only that lane:

```text
riff = sequence {
  notes     1 2 b3 4 5 6 b7 8
            |> rotate 2
            |> every 4 rev
            |> sometimes .20 (shift_degree 2)

  accent    + . . |> rotate 1
  duration  1 1/2 1/2 1 2 |> rev
}

play riff
```

A pipeline after the closing brace transforms the complete sequence:

```text
variation = sequence {
  notes   1 2 b3 4 5 6 b7 8
  accent  + . .
}
|> every 4 rev
|> sometimes .10 (fast 2)
```

Stitching sections:

```text
intro = riff.take(8)
verse = riff.take(16)

fill = sequence {
  notes    5 6 b7 8 5 >1
  accent   . + ++
  ratchet  . . 3 4
}

song = intro + verse.repeat(2) + fill
play song
```

Advantages:

- Most readable representation of several independent streams.
- No quotes or commas in the musical rows.
- Easy to edit one lane without touching the others.
- Imperative outer program with compact, local functional transformation
  chains.
- Does not commit the implementation to Lua, JavaScript, Wren, or a custom VM.

Costs:

- Requires a small top-level parser or a preprocessing/transpilation stage.
- Operator and mutation semantics must be specified carefully.
- Less immediately familiar than adopting an existing language unchanged.

## Option B: imperative object/property syntax

This follows the object-oriented side of Gibber and is compatible in spirit
with JavaScript, Wren, Python, and Lua, although the exact punctuation varies.

```text
riff = sequence()
riff.notes    = pattern("1 2 b3 4 5 6 b7 8")
riff.accent   = pattern("+ . .")
riff.duration = pattern("1 1/2 1/2 1 2")
riff.slide    = pattern(". . . >")
riff.ratchet  = pattern(". . 3")

play(riff)
```

Patterns can then be changed in place during a performance:

```text
riff.notes.rotate(2)
riff.accent.rotate(1)
riff.notes.shift_degree(2)
```

Stitching remains conventional:

```text
intro = riff.take(8)
verse = riff.take(16)
song = concat(intro, verse.repeat(2), fill)
play(song)
```

A constructor shorthand could reduce the setup:

```text
riff = sequence(
  notes    = "1 2 b3 4 5 6 b7 8",
  accent   = "+ . .",
  duration = "1 1/2 1/2 1 2",
  slide    = ". . . >"
)
```

Advantages:

- Imperative and unsurprising to general programmers.
- Re-evaluating one assignment naturally updates one aspect of a live object.
- Maps readily onto several embeddable languages.
- Requires little custom syntax outside the pattern strings.

Costs:

- Quotes, `pattern(...)`, punctuation, and repeated `riff.` add visual noise.
- The constructor form differs across candidate host languages.
- In-place mutation needs transactional semantics so a failed edit cannot leave
  half an update active.

## Option C: chained player syntax

This is closest to Strudel and the fluent parts of Gibber/Sardine.

```text
riff = notes("1 2 b3 4 5 6 b7 8")
  .accent("+ . .")
  .duration("1 1/2 1/2 1 2")
  .slide(". . . >")
  .ratchet(". . 3")

play(riff)
```

Manipulation and stitching:

```text
variation = riff
  .rotate_notes(2)
  .shift_degree(2)
  .rotate_accent(1)
  .every(4, reverse)

song = concat(
  intro,
  variation.take(16).repeat(2),
  fill
)

play(song)
```

Advantages:

- Compact and strongly compositional.
- Very familiar to Strudel and modern JavaScript users.
- Each call can return a new immutable pattern, simplifying rollback.
- Maps especially naturally to JavaScript/QuickJS or Wren.

Costs:

- Long chains become punctuation-heavy on a small panel.
- Targeting a single lane requires extra method names or a lane selector.
- More functional in feel than the requested imperative preference.
- Lua would require colon method syntax unless wrapped or preprocessed.

## Option D: Tidal-style control overlay

This retains the concise "everything is a pattern" character of Tidal while
changing the semantics so modifier streams advance per note rather than being
normalized automatically to a common cycle.

```text
riff = notes "1 2 b3 4 5 6 b7 8"
     # accent "+ . ."
     # duration "1 1/2 1/2 1 2"
     # slide ". . . >"
     # ratchet ". . 3"

play riff
```

Manipulation and stitching:

```text
variation = riff
  |> rotate notes 2
  |> rotate accent 1
  |> shift_degree 2
  |> every 4 reverse

song = concat intro (repeat 2 variation) fill
play song
```

Advantages:

- Densest option and close to established algorave notation.
- Uniform treatment of notes, velocity, duration, gates, and later CV streams.
- Pipelines make transformations easy to insert and remove live.

Costs:

- `#`, `|>`, and prefix calls have a learning curve.
- It is the most functional-looking option.
- Reusing Tidal punctuation with different timing semantics could mislead Tidal
  users.
- Requires a custom parser or extensive preprocessing.

## Option E: explicit imperative cursors

This takes Sonic Pi's independent named ticks to its logical conclusion. It is
best considered an advanced escape hatch rather than the everyday notation.

```text
notes     = cycle(1, 2, b3, 4, 5, 6, b7, 8)
accents   = cycle(+, ., .)
durations = cycle(1, 1/2, 1/2, 1, 2)
slides    = cycle(., ., ., >)

riff = sequence {
  loop {
    emit(
      note     = notes.next(),
      accent   = accents.next(),
      duration = durations.next(),
      slide    = slides.next()
    )
  }
}

play(riff)
```

It supports arbitrary conditionals and state naturally:

```text
if beat % 7 == 0 {
  notes.rotate(1)
}
```

Advantages:

- Most explicit and most naturally imperative.
- Independent phase and conditional advancement are easy to understand.
- Almost no compositional ceiling.

Costs:

- Far too verbose for the common case.
- Harder to analyze, bound, hot-swap, and execute safely in real time.
- Stitching declarative sections becomes less elegant.
- Encourages runtime callbacks when the audio engine should receive a compiled
  pattern graph.

## Comparison

| Option | Basic readability | Live mutation | Composition | Host coupling | Imperative feel |
| --- | --- | --- | --- | --- | --- |
| A. Lane block | excellent | excellent | good | very low | high |
| B. Object/property | good | excellent | good | medium | highest |
| C. Chained player | good | fair | excellent | medium | medium |
| D. Control overlay | compact but symbolic | fair | excellent | low/custom | low |
| E. Explicit cursors | verbose | excellent | fair | high | highest |

## Current recommendation

Prototype **A and B**, using exactly the same engine underneath:

- Option A is the strongest musical front end. It shows independent streams
  clearly, wastes little screen space, and keeps ordinary code imperative.
- Option B is the best baseline for judging whether a custom top-level syntax is
  worth its parser. It also maps most directly to embeddable imperative
  languages.

Option C is a viable fallback if fluent chaining feels especially good in a
playable test. Option D contributes useful operators but is not the best overall
fit for the stated imperative preference. Option E should inform a later
advanced API, not define v1.

The backend is now selected: cpp-peglib owns the declarative PEG grammar and
source locations, while TriggerFish owns musical validation and graph lowering.
The resulting bounded declarative pattern graph is the only representation the
audio runtime sees. A future advanced imperative layer can be evaluated on its
own merits without making the concise Option-A notation depend on that VM.

## References

- Tidal control patterns and combination:
  https://tidalcycles.org/docs/reference/controls/
  https://tidalcycles.org/docs/reference/pattern_structure/
- Strudel stepwise patterns, constructors, and velocity:
  https://strudel.cc/learn/stepwise/
  https://strudel.cc/learn/factories/
  https://strudel.cc/learn/effects/
- FoxDot tutorials and pattern functions:
  https://foxdot-community.codeberg.page/tutorials/
  https://foxdot681713046.wordpress.com/docs/pattern-functions/
- SuperCollider `Pbind` and pattern guide:
  https://doc.sccode.org/Classes/Pbind.html
  https://doc.sccode.org/Tutorials/A-Practical-Guide/PG_03_What_Is_Pbind.html
- Gibber:
  https://github.com/gibber-cc/gibber
- Sonic Pi rings and named ticks:
  https://sonic-pi.net/tutorial.html
- Sardine pattern language:
  https://sardine.raphaelforment.fr/pattern_languages/sardine.html
- Alda tutorial:
  https://alda.io/tutorial/
- Wren embedding and syntax:
  https://wren.io/embedding/
  https://wren.io/syntax.html
- QuickJS:
  https://github.com/bellard/quickjs/blob/master/doc/quickjs.texi
