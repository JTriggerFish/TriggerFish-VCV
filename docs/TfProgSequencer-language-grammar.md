# TfProgSequencer language grammar

This is the authoritative syntax and semantic contract for persisted language
version 5. Version 5 replaces the earlier prototype grammar outright. Saved
source still carries the version as an incompatibility marker, but this build
contains exactly one parser and compiler implementation.

This document is a compiler contract, not only a syntax sketch. Sections
labelled as future work reserve typed extension points but are not accepted by
the current grammar. Everything else is normative: the PEG creates typed nodes,
the semantic compiler applies the rules below, and unsupported combinations
are diagnostics rather than guesses. Whitespace separates adjacent pattern
elements but is otherwise insignificant outside comments and line endings.

The central decision is that `notes` contains complete musical events. Pitch,
rests, articulation, register, duration, repetition, ratchets, and local exact
values normally live together. Separate numerical lanes remain optional for
free-running modulation and sparse aligned overrides. A normal phrase should
not require an `articulation` lane.

## Event vocabulary

Articulation prefixes describe how an event is entered. Register and rhythmic
suffixes describe what happens after its pitch or chord is selected.

| Form | Meaning |
| --- | --- |
| `1` | ordinary attack |
| `x1` | ghost attack |
| `^1`, `^^1` | accent, strong accent |
| `>1`, `>^1` | slide into the target, optionally accented |
| `~` | silence |
| `_` | extend the preceding event by one base slot |
| `1_` | two-slot event |
| `1__`, `1_3` | three-slot event |
| `1.`, `1..` | dotted and double-dotted event |
| `1!3` | three consecutive logical events |
| `1*3` | three ratchets inside one logical event span |
| `1?.25` | event with probability `.25` |
| `[1 2 3]` | three events subdividing one parent slot |
| `[1\|2]` | deterministic random choice |
| `<1 2 3>` | pass-by-pass alternation |
| `1(3,8,1)` | rotated Euclidean distribution |

The only accepted reading order is:

```text
[slide] [dynamic] [pitch or chord] [register] [duration]
  [Euclidean] [ratchet] [probability] [replication] [exact attributes]
```

For example, `>^1'_3*2?.8!4` means an accented slide into the upper
tonic, with a three-slot span, two ratchets, 80 percent probability, replicated
as four logical events. Dense combinations are legal but ordinary source
should remain sparse. A suffix written out of order is a syntax error; the
compiler never silently reorders source operators.

Euclidean notation is attached to the event or group it distributes:

```text
1(3,8)       // three attacks in eight cells: x..x..x.
1(3,8,1)     // rotate one cell later/right: .x..x..x
1(3,8,-1)    // rotate one cell earlier/left
^1(5,8)      // accents on every active cell
[1 3](3,8)   // the complete subdivision occupies each active cell
```

The arguments are `(pulses, steps, rotation?)`; rotation defaults to zero and
wraps modulo `steps`. Inactive cells are true rests and still consume score
time and aligned-lane positions. Ratchet, probability, and replication remain
outside the Euclidean expansion according to the canonical suffix order, so
`1(3,8)*2?.7!4` means four copies of an eight-cell rhythm whose active cells
each make one probability decision and, when sounding, emit two ratchets.

`x` is a ghost only as a prefix attached to a pitched note or chord. Standalone
`x` remains available to a future drum or rhythm-mask grammar, where it can
naturally mean a hit. There is no separate `articulation` lane.

`^` is emphasis instead of overloading `!` at both ends of an event. `!` is
therefore always replication when it follows an event. A future general signal
expression may still use `^` differently because expressions and note
mini-notation have separate PEG contexts.

## Registers

`@` remains absolute and follows the pitched value. Apostrophes and commas are
relative to the active sequence register:

```text
D@4      // absolute D4
Cm7@3    // chord rooted absolutely in octave 3
1'  1'' // one and two octaves above the active register
1,  1,, // one and two octaves below it
```

Registers precede duration and structural suffixes, so `^1'_3`, `x3,`,
`Cm7,_2`, and `(1 b3 5)'` are unambiguous. A relative suffix contains either
apostrophes or commas, never a mixture such as `1',`. An absolute register may
receive one relative offset, although `D@5` is clearer than `D@4'`.

Continuing degrees retain their existing scale-cardinality meaning. In a
seven-note scale, `8` is the next tonic; in a five-note scale, `6` is the next
tonic. This is orthogonal to register punctuation.

## Scale degrees and named pitches

`notes` accepts scale degrees, chromatic named pitches, explicit voicings, and
jazz chord symbols in the same pattern:

```text
notes 1 2 b3 Eb F#@3 5' D (C E G)@4 Cm7
```

Degree accidentals modify the degree selected from the active scale, so `b3`
means one semitone below that scale's third degree. Named pitches use uppercase
letters with `b` or `#` accidentals and bypass scale lookup. A bare named pitch
uses the active `octave` lane; `F#@3` supplies an absolute register. `C4` is not
an octave spelling: `@` remains mandatory so jazz extensions and note registers
cannot collide. Thus `D@4` is the note D4 and `D7` is a dominant-seventh chord.

Degree arithmetic is exact for any signed integer. For a scale with `K` ordered
intervals, degree `D` uses `q,r = floor_divmod(D-1, K)` and resolves to
`12q + scale[r]` semitones from the tonic before written accidentals. Thus `1`
is tonic, `K+1` is the next tonic, `0` is the last scale member below tonic,
and negative values continue downward without a C++ remainder-sign ambiguity.

The transformation domains are deliberate:

- `scale` and `shift_degree` remap degree events only;
- named notes and named chord roots retain their chromatic identity under a
  scale change;
- `transpose_semitone` and `transpose_octave` affect degrees, named notes, and
  chords after their pitch material is resolved;
- `modulate_degree` derives a chromatic interval from the active scale and then
  transposes the complete event, including named notes and chords, preserving
  its chromatic shape;
- a future `midi.root` binding likewise applies after event selection so it may
  transpose either a degree sequence or a named-note sequence coherently.

“Absolute” in `F#@3` therefore fixes the source register; it does not make that
note immune to explicit later transposition. A future literal-voltage or
frequency event, if needed, should use a separate typed form.

## Roman-degree and explicit-root chords

Chord events accept both Roman scale-degree roots and explicit A--G roots:

```text
notes I i ii7 iim7 V7 bIImaj7 #ivdim7
      Cm7 D7 Bbm7b5@3 / D@2
```

A numeric degree remains a single note, while a Roman numeral is a semantic
chord. The active scale supplies a Roman chord's root pitch; case supplies its
default triad quality independently of the scale:

- `I` is a major triad on scale degree 1;
- `i` is a minor triad on scale degree 1;
- `II` and `ii` are major and minor triads on degree 2;
- `ii7` is a minor-seventh chord;
- `iim7` is an accepted explicit spelling of the same minor-seventh quality;
- `V7` is a dominant-seventh chord;
- `Imaj7`, `idim7`, `iim7b5`, `Isus4`, and alterations use the same quality
  vocabulary as explicit-root jazz symbols.

The initial chord vocabulary is deliberately finite and tokenized rather than
captured as arbitrary text:

- explicit roots are uppercase `A` through `G` with at most one `b` or `#`;
- triad qualities are `maj`, `m`, `min`, `dim`, `aug`, `sus2`, and `sus4`;
- extensions are `5`, `6`, `7`, `9`, `11`, and `13`;
- alterations are `b` or `#` followed by `1`, `2`, `3`, `4`, `5`, `6`, `7`,
  `9`, `11`, or `13`;
- additions use `add` followed by the same degree vocabulary.

At least one quality, extension, alteration, or addition must follow an
explicit root for it to be a chord: `D` is a named note, while `Dmaj`, `Dm`, and
`D7` are chords. `maj` selects a major triad and a major seventh when an
extension of seven or greater is present; an unqualified `7`, `9`, `11`, or
`13` uses a minor seventh. Unsupported, contradictory, or duplicate quality
parts are semantic diagnostics with their own source spans.

Leading `b` and `#` alter the resolved Roman root chromatically, so `bII` is a
major triad one semitone below the scale's second degree. Roman numeral case
must be uniform; mixed forms such as `Iv` are diagnostics. A Roman root names
one of the seven conventional harmonic degrees `I`--`VII`; that syntactic range
does not grow with scale cardinality. It selects that ordinal member without
wrapping: a Roman degree
greater than the active scale's cardinality is a diagnostic. A pentatonic scale
therefore accepts `I`--`V`, while an octatonic scale still deliberately stops at
`VII`; use an explicit root or a future scale-stacked chord generator for other
harmony. Register marks express the next scale octave, so `II'` is the upper
version of degree 2 rather than `IX`. Numeric note degrees continue across
scale octaves and may address additional degrees of octatonic or future larger
scales; explicit-root chords and a future typed chord generator cover harmony
outside the seven Roman functions.

For Roman roots, case supplies the triad only when no explicit triad quality is
present. A repeated compatible spelling such as `iim7` is accepted. A
contradiction such as `Im7` or `imaj7` is rejected instead of silently
overriding Roman case; `dim`, `aug`, and `sus` are explicit non-major/minor
qualities and may follow either case. This keeps `I`, `i`, and jazz quality
suffixes readable without maintaining two chord parsers.

Restricting Roman roots to seven degrees also keeps lowercase `x`
unambiguously available as the ghost prefix: `xii7` is a ghost `ii7`, never a
Roman degree twelve chord.

Registers and event operators compose normally:

```text
^V7'      // accented dominant, one relative octave up
xii7,     // ghost minor seventh, one relative octave down
>bII@3_2  // slide to an absolute-register flat-II chord for two slots
```

A slash bass is a pitch, not another chord. It may be a scale degree or named
pitch:

```text
I / 3
V7 / F#@2
Bbm7b5 / D@2
```

Slash syntax is valid only after a semantic chord or explicit voicing, never a
single note. The canonical explicit-voicing form is `(C E G) / B`; all tones
inside the parentheses are pitches, not chord symbols. A slash bass changes
inversion metadata but does not replace the chord root.

Roman chords preserve semantic root degree, accidental, stated quality,
extensions, alterations, slash bass, and original source span in the compiled
graph. A chord interpreter can therefore realize `ii7` as a rootless voicing,
walking bass source, arpeggio, or another policy without reconstructing intent
from emitted pitches.

Transform behavior follows the note rules with one extra distinction:

- `scale` changes the pitch of a Roman root but does not silently change its
  case-stated chord quality;
- `shift_degree` moves the Roman root within the active scale while preserving
  the stated quality;
- `modulate_degree`, `transpose_semitone`, and `transpose_octave` move the
  complete chord while preserving its interval structure;
- explicit-root jazz chords retain their root under `scale` and
  `shift_degree`, but respond to modulation and chromatic/octave transposition.

Roman case is therefore an explicit quality choice, not shorthand for
automatically stacking thirds from the current scale. A future scale-stacked
chord generator should be a separate typed operation such as `chord 2`, so a
scale edit cannot unexpectedly turn an explicitly written `II` major chord
minor.

## Duration, gate, and repetition

Duration is rhythmic span: the time from an onset to the next logical event.
Gate is how long the Gate output stays high inside that span. Slide is the
pitch slew time. They remain separate values. `gate=1` fills a span but does
not turn the following attack into a tie: only standalone `_` suppresses the
next onset while preserving pitch, and only `>` creates a legato pitch change.
Fractional Gate is multiplied by event span; millisecond Gate is capped at that
span. Slide time is likewise capped at the target span. Ratchet sub-triggers may
occur while a still-high Gate overlaps their subspan, but neither Gate nor
slide leaks past the logical event unless a tie/slide rule explicitly holds it.

The base duration comes from the independently cycling `duration` lane or its
one-beat default. Inline duration multiplies that base:

- `1_` has total weight 2;
- `1__` and `1_3` have total weight 3;
- `1.` has weight `3/2`;
- `1..` has weight `7/4`;
- `1_3.` has weight `9/2`.

Numeric `_N` states total weight, not the number of added slots. Standalone `_`
is deliberately different: it is a semantic tie cell. It extends the preceding
pitched event by one copy of that event's base duration, produces no onset, and
forces Gate to remain continuous across the boundary. `_!3` creates three tie
cells. A duration suffix lengthens a note's scheduling span but retains its
normal Gate setting, so `1_3{gate=.5}` is a three-slot event with a half-span
Gate whereas `1 _!2` is held continuously for all three slots. This distinction
is required because note duration and Gate duration are independently useful.
Exact nonstandard spans use `1{len=5/4}`.

The `duration` lane advances once for every selected note or rest cell,
including a note later suppressed by probability. It does not advance for a
tie cell or a ratchet. A tie reuses the preceding event's already resolved base
duration; editing another lane cannot change the length of a tie after its
source event has begun.

Replication and ratcheting differ deliberately:

- `1!3` creates three logical events; onset-property lanes advance for every
  copy.
- `1*3` resolves one logical event and creates three sub-onsets inside its
  span; all ratchets inherit the resolved event properties.
- `1_3*4` is one three-slot event containing four ratchets.

A ratchet divides its resolved span into equal subspans. Each sub-onset emits a
Trigger and applies Gate as a fraction of that subspan. An absolute millisecond
Gate may overlap the next ratchet; this keeps Gate high while Trigger still
marks every sub-onset. Velocity, Accent, pitch, probability result, inline
attributes, and CV cell position are shared by all ratchets.

Probability turns a pitched event into a rest without removing its time. Bare
`?` means `.5`. It is decided once per logical event, before ratchets are
emitted, so either all of a note's ratchets sound or none do. Replicated copies
receive independent deterministic decisions because their event paths differ.
Probability on a literal rest is rejected as a no-op.

## Pattern structure and operator scope

Square brackets divide one parent span among their selected children. With no
weights, `[1 2 3]` produces three equal thirds. Child duration values are
relative weights inside the bracket, so `[1_ 2]` divides the parent span `2:1`;
they do not expand the bracket to three beats. A duration suffix on the bracket
changes its parent span: `[1 2]_2` occupies two base slots. This normalization
applies recursively.

`[a | b | c]` chooses one complete branch deterministically from the program
seed and event path. A branch may contain several elements, for example
`[1 2 | 3 [4 5]]`. `<a b c>` chooses one element on each pass through its
containing pattern; use a bracketed subdivision when one alternate must contain
several events, as in `<[1 2] [3 4]>`.

The PEG retains those complete branches in the typed tree. The current direct
interpreter executes atomic alternatives such as `<1 3 5>` and `[1|3|5]`.
Multi-event and recursively nested choice branches are intentionally a compile
diagnostic until the prepared branch graph is wired; they are specified here
so that implementation does not require another syntax change. Edge-aligned
lanes additionally require equal structural cell counts across every branch.

Duration, Euclidean, probability, and replication suffixes may follow a square
or angle group and apply to the group as a unit:

```text
[1 2]_2       // the complete subdivision occupies two base slots
[1 2](3,8)    // place the complete subdivision on three of eight slots
<1 3>?.5      // choose the pass alternate, then probabilistically sound it
[1 2]!3       // three consecutive copies of the subdivision
```

Ratchets and exact event attributes apply only to pitched atomic events. Group
ratchets and group attributes are rejected because distributing them over
nested choices would hide scope. A Euclidean suffix creates its stated number
of logical cells; active cells contain the suffixed event or group and inactive
cells are rests. The hit count must be from zero through the positive step
count; rotation is any integer and is reduced modulo the step count. Positive
rotation moves hits later/right by that many cells and negative rotation moves
them earlier/left.

Every suffix has one decision level:

- group probability gates the selected group as a unit while preserving its
  complete span;
- event probability gates one logical event;
- Euclidean expansion occurs before probability at each generated hit;
- replication is outermost and gives each copy an independent event path; and
- a random or alternate branch is selected before suffixes on its containing
  group are applied.

## Articulation semantics

An accent is a normal triggered onset with an accent value and velocity floor.
Initial built-ins remain `.88` for `^` and `1` for `^^`.

A ghost remains a triggered onset. After numerical lanes and exact inline
values resolve, it attenuates velocity by `.5`, suppresses Accent output, and
caps Gate to a short value. The built-in cap is the lesser of ten percent
of event span and 20 ms once an input-clock period is known; before the first
measured period, ten percent of span is used. A quiet sustained note should use
explicit velocity and gate without the `x` prefix.

A slide target suppresses Trigger, keeps the preceding Gate continuous, and
slews pitch for the resolved `slide` time. `>^1` is valid. Ghost plus accent is
impossible because dynamic is one choice. Ghost plus slide is rejected
initially because short-gate ghost semantics contradict legato slide semantics.

A slide requires a preceding pitched candidate in the same statically ordered
branch. A literal rest breaks that relationship. A failed-probability slide is
a rest and does not hold the source Gate. Slide plus ratchet is rejected because
one requests trigger suppression while the other requests sub-triggers.
Replication copies articulation, so `>3!2` makes two legato copies; use `>3 3`
when only the first target should slide.

The compiler marks an ordinary source event before `>` as legato-to-next, which
holds its Gate through the boundary regardless of its normal Gate fraction. A
ghost source may not feed a tie or slide because either would contradict the
ghost Gate cap. If an otherwise valid source was silent because its probability
failed,
the slide target degrades to an ordinary triggered attack rather than
disappearing. A tie after a probability miss remains a rest until a later
pitched attack establishes a new source.

Slides between single pitches are fully defined. A chord interpreter receives
the same semantic slide marker and owns its voicing policy. The initial direct
polyphonic interpreter may slide chord-to-chord only when source and target have
the same voice count, matching source order; other chord slides are a compile
diagnostic rather than an implicit and unstable voice-leading algorithm.

Rests consume duration but not pitch-onset properties such as velocity, gate,
accent, or slide. Ties advance time-indexed CV lanes but no onset-property lane.
An extension without a preceding pitched event is a diagnostic. Initially a
tie or slide may not reach backward across a random/alternate branch boundary,
an arrangement boundary, or the start/end wrap of a sequence; those cases need
explicit branch-aware dataflow before they can be admitted safely.

## Exact inline attributes

Rare precision has a named, order-independent escape hatch:

```text
notes 1 2{vel=.63, gate=12ms}
      >3{slide=80ms}
      Cm7_3{gate=1}
      5{len=5/4}
```

Initial keys are `len`, `gate`, `vel`, and `slide`. Bare `gate` values from zero
to one are fractions of event span; `ms` supplies absolute time. Bare `len` and
`slide` values use incoming-clock beats, while slide also accepts `ms`.
Supplying both a duration suffix and `len` is a duplicate specification. Every
attribute requires `=`, keys may occur only once, and unknown keys are
diagnostics. A rest accepts `len` only. `vel`, `gate`, and `slide` apply only to
pitched events, and `slide` on an event without the `>` prefix is rejected as an
unused specification.

Resolution order is:

1. sequence default;
2. free-running or aligned numerical lane;
3. exact inline attribute;
4. articulation constraint such as ghost attenuation/gate cap, accent floor,
   or slide trigger suppression.

Inline attributes travel with notes under note transformations. Aligned lane
overrides belong to playback positions instead.

## Numerical lanes

Without alignment syntax, numerical lanes keep their polymetric behavior:

```text
notes    1 2 3 4 5 6 7 8
velocity .72 .55
gate     .8 .4 .6
slide    . 80ms
duration 1 1 3/2
offset   -10ms!2 +8ms
```

Each short lane loops independently. `.` is the typed no-op meaning “inherit
the normal value for this cell”; `.1` is numeric one tenth. Event prefixes no
longer use `.`, so there is no `notes .1` versus `gate .1` ambiguity.

The accepted scalar lanes and domains are:

| Lane | Values | Default/no-op | Advancement |
| --- | --- | --- | --- |
| `octave` | signed integer octave offset | zero | pitched logical events |
| `velocity`, `vel` | finite `0..1` | `.72` | pitched logical events |
| `accent` | finite `0..1` | no accent | pitched logical events |
| `duration`, `dur` | positive beat value or fraction | one beat | notes and rests |
| `gate` | fraction `0..1`, or non-negative `ms` | `.8` of span | pitched logical events |
| `slide` | non-negative beats/fraction or `ms` | `glide` setting | pitched logical events |
| `ratchet` | positive integer | one onset | pitched logical events |
| `offset` | signed beat/fraction or `ms` displacement | on-grid | incoming-clock score time |

An accent-lane value raises Velocity to at least that value and drives Accent
output for the effective Gate; zero explicitly produces no accent. A `slide`
value is sampled on every pitched logical event so polymetric phase remains
stable, but it is audible only when the event has `>`. A `ratchet` lane value
and an inline `*N` on the same event are duplicate specifications.

`offset` is the patternable form of early/late timing. Negative values schedule
early and positive values schedule late, so `offset -10ms!2 +8ms` is early for
two incoming-clock beats and late for one. It repeats as a three-beat pattern
independently of note density: subdivisions sample the score-time lane more
often but do not make it run faster. A positive finite `rate R` scales a
free-running lane's phase:

```text
offset -10ms!2 +8ms |> rate 1/2
```

This runs the offset phrase at half speed. Whole-sequence `early A` and
`late A` remain constant additive conveniences and never alter lane phase.
Rate-scaled and edge-aligned lanes are distinct modes; combining `rate` with
`...` is initially a diagnostic.

“Pitched logical event” includes an event suppressed by probability: lane phase
is a property of the score, not of the random result. It excludes rests and tie
cells. Ratchets do not advance any scalar lane. The `duration` lane is
time-indexed and therefore also advances on rests. Aligned lanes use structural
position instead of these free-running advancement rules; ignored values at a
rest or tie still retain their aligned positions.

Bare numeric values use normalized units only where the table says so. `ms` is
case-sensitive. Scientific notation, `NaN`, and infinities are initially
rejected; a later general expression grammar may produce checked finite values
without changing mini-notation tokens.

### Edge-aligned ellipsis

One edge ellipsis changes a numerical lane from a free-running pattern into an
aligned override map for one structural cycle of `notes`:

```text
notes 1 2 3 4 5 6 7 8 9
gate  ... .1!2
```

Leading `...` right-aligns the explicit cells, producing defaults in the first
seven note cells and `.1` in the final two. Trailing `...` left-aligns:

```text
gate .1!2 ...
```

An explicit `.` consumes one aligned cell without changing its normal value:

```text
gate ... .1 . .2
```

The last three note cells therefore receive `.1`, no override, and `.2`.
Ellipsis is permitted only once and only as the first or final term. Middle
forms such as `.2 ... .1` and two-ended forms such as `... .1 ...` are rejected.

Rules:

- terms other than `...` expand repetitions such as `.1!2` before alignment;
- ellipsis pads with exactly enough no-op cells to match the statically
  expanded note-cell count;
- too many explicit cells is a compile diagnostic;
- replication creates cells; ratcheting and elongation remain properties of
  one cell;
- a Euclidean step count creates that many cells;
- rests occupy aligned cells even though onset-only properties are not read at
  those cells;
- tie cells occupy aligned positions; only time-indexed CV lanes read them;
- subdivisions flatten to their ordered child cells;
- random and alternate branches must have equal structural cell counts when
  referenced by an aligned lane; unequal branches are initially rejected;
- note transformations move inline properties with notes, while an aligned
  lane resolves against resulting playback positions.

This explicit mode preserves `gate .8 .4` as a free-running two-cycle pattern.
Ellipsis never changes an ordinary lane merely because its size matches notes.

A future editor view may render the expanded matrix without rewriting concise
source. Displayed cells must retain PEG/semantic provenance so selection,
diagnostics, and play cursors remain meaningful.

## General CV lanes

`cv1` and `cv2` are time-indexed control lanes belonging to the same sequence:

```text
notes 1 2 3 4
cv1   0 5 -2 3
cv2   0 . . 5 |> interp linear
```

Bare values are finite volts. There is no arbitrary language-level voltage
ceiling; the eventual output stage may offer an explicit safety clamp. The
default interpolation is `step`. Every logical note cell samples CV, including
rests, probability misses, and standalone tie cells; ratchets share their
parent cell. A free lane's phase is nevertheless score-time based, so an
elongated event can cross lane cells without advancing an event cursor.

`cvN` uses the same free-running, repetition, `.` no-op, edge-ellipsis, and
numeric `rate R` forms as other scalar lanes. Free-running CV phase derives
from absolute score time so independent rational rates do not accumulate drift.
The implementation provides physical outputs
for `cv1` and `cv2`; the grammar recognizes positive `cvN` names so later output
counts do not require another syntax change. An unavailable output number is a
semantic diagnostic.

For a stepped CV, `.` consumes a lane cell and holds the preceding explicit
point. For an interpolated CV, it consumes a lane cell but creates no control
point, so the curve continues to the next explicit value. Leading no-ops look
backward through the cyclic lane to its final explicit point; a lane containing
no explicit value outputs 0 V.

Interpolation is a lane pipeline:

```text
cv1 0 5 0             // stepped sample-and-hold
cv1 0 5 0 |> interp step
cv1 0 5 0 |> interp linear
cv1 0 5 0 |> interp smooth
cv1 0 5 0 |> interp power 2
```

`linear` uses constant slope, `smooth` uses cubic smoothstep, and `power P`
uses `a + (b-a)t^P` with finite `P > 0`; it is shape interpolation rather than
an exponential voltage domain and therefore works through zero and negative
values. At each logical event, interpolation starts from the current output and
targets the next explicit point on the lane's incoming-clock score-time grid,
including across the lane loop. If no event occurs at a later knot, the output
reaches the pending target and holds until the next event; a future prepared
control-event stream can traverse every intervening segment without changing
this syntax. Subdivisions may sample the curve more frequently and ratchets
share one sample, but neither advances its phase; `rate` is the explicit way to
change that phase. Final event timing transforms move the emitted sample and
its target together. Hot replacement starts from the current output value and
joins the newly compiled curve at its next knot, avoiding a forced
discontinuity.

Edge-aligned CV currently supports `step` only. Continuous interpolation on an
ellipsis-aligned lane is a diagnostic until prepared structural-boundary times
can supply exact knot durations for unequal note lengths.

CV values are monophonic sequence controls initially. When a note/chord emits
polyphony, the value is broadcast to all output channels. A later typed
per-voice expression can specialize it without changing the lane syntax.

Top-level signal expressions remain a separate future grammar:

```text
cv2 = exp(0.1 * vel) + cv1
```

The `=` makes this unambiguous from a `cv2 ...` lane inside a sequence. Such
expressions are reserved and are not accepted until the signal-expression
compiler and bounded runtime evaluator exist.

Live MIDI binding is likewise a typed future stage, not an overloaded note
token. The canonical pattern-driven arpeggio keeps the compiled sequence on the
left because it owns musical time:

```text
arp_voice = arp1_seq
  |> input midi.chord
  |> interpret arp.index
```

The conceptual source-first form `midi.in |> arp arp1_seq` may be future sugar
for the same `(PerformanceSource, EventPattern) -> PitchedEventStream` node.
Bare `midi.in |> arp1_seq` remains invalid because ordinary sequence values are
not implicitly callable. `arp.index` selects held MIDI notes with the written
event indices; interpreters such as `arp.up` instead generate bounded
sub-events inside an event span. These productions remain reserved until MIDI
input, `input`, and `interpret` have typed AST/runtime support.

## Document, arrangement, and transform syntax

The language uses this document shape:

```text
acid = sequence {
  cycle 16
  tonic D#@2
  scale minor
  notes 1 ^1 >5 ~ 1 _
  velocity .72 .55
  cv1 0 5 0 |> interp smooth
}

iv   = acid |> modulate_degree 4
song = acid * 8 + iv * 4
play song
```

`//` starts a line comment. A sequence body accepts one `notes` lane, optional
scalar/CV lanes, and the settings `cycle`, `tonic`, `scale`, and `glide`.
Aliases are limited to `vel` and `dur`; adding aliases casually makes future
identifiers harder to reserve. A lane or sequence may continue pipelines on
following indented lines. Duplicate settings or canonical lane names are
diagnostics even when different aliases are used.

| Setting | Domain and default |
| --- | --- |
| `cycle` | positive arrangement span in incoming-clock beats; `8` |
| `tonic` | named pitch with optional register; `C@4` |
| `scale` | registered scale identifier; `major` |
| `glide` | non-negative beats/fraction or `ms`; `.25` beat |

`cycle` controls when an arrangement advances; it does not rescale the note
pattern. A shorter pattern loops inside the cycle and a longer one continues
from its preserved phase when that sequence is visited again. Pipeline-only
continuation lines attach to the immediately preceding lane, sequence closing
brace, or assignment. A blank or another statement ends that attachment; an
orphan continuation is a diagnostic.

At an arrangement boundary, pending ratchets/subdivision events belonging to
the old part are cancelled and an event that would cross the boundary is
clipped there. The new part owns the boundary attack. Cross-part ties and
slides are deliberately invalid in this version, so Gate does not acquire an
implicit legato rule from adjacency between two named sections.

The initial scale registry contains `major`/`ionian`, `minor`/`aeolian`,
`harmonic_minor`, `dorian`, `phrygian`, `lydian`, `mixolydian`, `locrian`,
`major_pentatonic`/`pentatonic`, `minor_pentatonic`,
`octatonic_whole_half`/`whole_half`, and
`octatonic_half_whole`/`half_whole`. Aliases resolve to one canonical scale ID
in the semantic graph. Adding a scale is data, not a grammar change.

At top level, `+` concatenates named sequences/arrangements, `*` repeats its
immediately preceding term by a positive integer, and parentheses group those
operators. `*` binds more tightly than `+`; a pipeline applies after the full
concatenation to its left unless parentheses establish a smaller assignment.
Assignments are immutable graph definitions. Recursive references and duplicate
definitions are diagnostics. `play NAME`, `stop`, and `seed UNSIGNED_INTEGER`
are the only transport/state statements in this version.

Keywords, setting/lane names, transform names, and names matching `cv` followed
by a positive integer are reserved and cannot be user definition names. In
particular, the compiler rejects `cv2 = acid` rather than treating it as an
arrangement; this preserves `cv2 = signal_expression` for the typed future
grammar without later reinterpretation.

Pipelines are left-to-right and every transform requires its documented numeric
argument; `fast` and `slow` are never bare qualitative commands:

| Transform | Meaning |
| --- | --- |
| `rev` | reverse structural event order |
| `rotate N` | rotate structural cells by signed integer `N` |
| `shift_degree N` | remap every degree by signed diatonic index offset `N` |
| `modulate_degree N` | transpose the tonal centre by scale degree `N`, preserving chromatic shapes |
| `transpose_semitone N` | transpose all resolved pitch by signed semitones |
| `octave N`, `transpose_octave N` | transpose all pitch by signed octaves |
| `scale NAME` | change the active scale at this point in the pipeline |
| `fast F` | divide score time by positive factor `F` |
| `slow F` | multiply score time by positive factor `F` |
| `rate R` | multiply one free-running lane's independent phase by positive ratio `R` |
| `swing R` | swing the event's own subdivision with `.5 <= R < 1` |
| `swing R GRID` | swing a positive incoming-clock beat/fraction grid |
| `early A`, `late A` | fixed non-negative beat/fraction or `ms` displacement |
| `early random A`, `late random A` | deterministic displacement from zero through `A` |
| `interp MODE ...` | CV-lane-only interpolation described above |

For swing, `.5` is straight and `R` is the fraction of each adjacent grid pair
given to its first member. Swing and early/late displacement change scheduled
positions but not score span. Events that land at the same sample are ordered
by stable source path. `fast` and `slow` scale event time inside the containing
sequence cycle; the cycle/arrangement window remains on the authoritative
external-clock grid. Thus `fast 2` fits two passes of an eight-beat source into
one eight-beat cycle rather than ending an arrangement term early. Random
early/late amounts use a deterministic uniform value in the stated closed
interval.

Degree arguments are index offsets: `shift_degree 1` moves tonic to the second
scale degree. Modulation arguments are one-based destinations:
`modulate_degree 1` is unison, `4` moves the centre to the fourth, and negative
values select downward scale degrees. These intentionally different domains
must use separate compiler types so an off-by-one conversion cannot leak from
one operation into the other.

More precisely, `modulate_degree N` uses the continuing positive degree
interval for `abs(N)` and applies the sign of `N`; both `0` and magnitude `1`
are unison. Therefore `-4` is a downward fourth while `shift_degree -4` is four
scale indices downward. This preserves the established modulation behavior
without calling either operation the ambiguous `transpose_degree`.

`every N TRANSFORM` applies its transform on every positive `N`th pass;
`sometimes P TRANSFORM` applies it with deterministic probability `0..1`.
Their final argument is one transform call, optionally parenthesized when it has
arguments: `|> every 4 (rotate 1)`. Conditions cannot directly contain another
condition in this version. Randomness is derived from seed, stable definition
identity, pass, and event path, never mutable call order.

Pitch transforms carry inline event attributes with their events. An aligned
scalar or CV lane is applied to the resulting playback positions after
structural transforms. Free-running lanes retain their own phase and compatible
lane-local transforms. Compatibility is closed rather than inferred:

- `notes` accepts `rev`, `rotate`, pitch transforms, and conditions around
  those transforms;
- scalar lanes accept `rev`, `rotate`, `rate`, and conditions around structural
  transforms;
- CV lanes accept `rev`, `rotate`, `rate`, `interp`, and conditions around
  `rev` or `rotate` only;
- a containing sequence or arrangement accepts structural, pitch, scale, and
  timing transforms, which move all note and CV boundaries coherently; and
- `rate` and `interp` are never sequence transforms, while `fast`, `slow`,
  `swing`, `early`, and `late` are never lane-local in this version.

This matrix can later grow by adding typed operations; accepting an operation
in an undefined domain and hoping it sounds sensible is not permitted.

## PEG-shaped grammar

The production names specify structure; cpp-peglib may factor lexical details
differently. It must emit these typed nodes directly. The semantic compiler
validates chord qualities, duplicates, articulation combinations, units, and
ranges; it must not recover event structure by splitting an opaque token.

```ebnf
Document        <- (SequenceDefinition / Play / Stop / Seed / Assignment /
                    BlankLine)* Trailing End
SequenceDefinition
                <- H Identifier H '=' H 'sequence' H '{' H Comment? Newline
                   SequenceBodyLine*
                   H '}' H Comment? LineEnd PipelineContinuation*
SequenceBodyLine
                <- NotesLane / ScalarLane / CvLane / SettingLine /
                   PipelineContinuation / BlankLine
NotesLane       <- H 'notes' S NotePattern (H Pipeline)* H Comment? LineEnd
ScalarLane      <- H ScalarLaneName S ScalarPattern
                   (H Pipeline)* H Comment? LineEnd
CvLane          <- H CvLaneName S ScalarPattern
                   (H Pipeline)* H Comment? LineEnd
SettingLine     <- H SettingName S SettingValue H Comment? LineEnd
PipelineContinuation
                <- H Pipeline (H Pipeline)* H Comment? LineEnd

Play            <- H 'play' S Identifier H Comment? LineEnd
Stop            <- H 'stop' H Comment? LineEnd
Seed            <- H 'seed' S UnsignedInteger H Comment? LineEnd
Assignment      <- H Identifier H '=' H Expression H Comment? LineEnd
                   PipelineContinuation*
Expression      <- Concatenation (H Pipeline)*
Concatenation   <- Term (H '+' H Term)*
Term            <- Primary (H '*' H PositiveInteger)?
Primary         <- Identifier / '(' H Expression H ')'

Pipeline        <- '|>' H (ConditionalTransform / TransformCall)
ConditionalTransform
                <- 'every' S PositiveInteger S TransformCall /
                   'sometimes' S Probability S TransformCall
TransformCall   <- '(' H Operation (S TransformArgument)* H ')' /
                   Operation (S TransformArgument)*
Operation       <- Identifier
TransformArgument
                <- ScalarValue / Identifier

NotePattern     <- NoteElement (S NoteElement)*
NoteElement     <- GroupElement / NoteEvent / RestEvent / TieExtension
GroupElement    <- GroupPrimary DurationSuffix? EuclideanSuffix?
                   ProbabilitySuffix? ReplicationSuffix?
GroupPrimary    <- BracketGroup / Alternate
BracketGroup    <- '[' H (RandomChoice / NotePattern) H ']'
RandomChoice    <- NotePattern H '|' H NotePattern
                   (H '|' H NotePattern)*
Alternate       <- '<' H NoteElement (S NoteElement)* H '>'

NoteEvent       <- OnsetPrefix? PitchedValue DurationSuffix?
                   EuclideanSuffix? RatchetSuffix? ProbabilitySuffix?
                   ReplicationSuffix? Attributes?
RestEvent       <- '~' DurationSuffix? ReplicationSuffix? Attributes?
TieExtension    <- '_' ReplicationSuffix?

OnsetPrefix     <- SlidePrefix? DynamicPrefix?
SlidePrefix     <- '>'
DynamicPrefix   <- '^^' / '^' / 'x'

PitchedValue    <- ChordValue (H '/' H SlashBass)? /
                   PitchValue
ChordValue      <- (ExplicitVoicing / RomanChord / JazzChord) RegisterSuffix?
PitchValue      <- (NamedPitch / ScaleDegree) RegisterSuffix?
ExplicitVoicing <- '(' H PitchValue (S PitchValue)+ H ')'
SlashBass       <- PitchValue

RomanChord      <- RootAccidental* (UpperRoman / LowerRoman) ChordSuffix?
JazzChord       <- NamedRoot ChordSuffix
ChordSuffix     <- TriadQuality ChordExtension? ChordModification* /
                   ChordExtension ChordModification* /
                   ChordModification+
TriadQuality    <- 'sus2' / 'sus4' / 'maj' / 'min' / 'dim' / 'aug' / 'm'
ChordExtension  <- '13' / '11' / '9' / '7' / '6' / '5'
ChordModification
                <- ('b' / '#') ChordDegree / 'add' ChordDegree
ChordDegree     <- '13' / '11' / '9' / '7' / '6' / '5' / '4' / '3' / '2' / '1'
RootAccidental  <- 'b' / '#'
UpperRoman      <- 'VII' / 'VI' / 'V' / 'IV' / 'III' / 'II' / 'I'
LowerRoman      <- 'vii' / 'vi' / 'v' / 'iv' / 'iii' / 'ii' / 'i'
NamedPitch      <- NamedRoot
NamedRoot       <- [A-G] ('b' / '#')?
ScaleDegree     <- RootAccidental* SignedInteger

RegisterSuffix  <- AbsoluteRegister RelativeRegister? / RelativeRegister
AbsoluteRegister
                <- '@' SignedInteger
RelativeRegister
                <- Apostrophe+ / Comma+
Apostrophe      <- "'"
Comma           <- ','

DurationSuffix  <- Elongation Dots? / Dots
Elongation      <- '_' PositiveInteger / '_'+
Dots            <- '..' / '.'
EuclideanSuffix <- '(' H EuclideanPulses H ',' H EuclideanSteps
                   (H ',' H EuclideanRotation)? H ')'
EuclideanPulses <- UnsignedInteger
EuclideanSteps  <- PositiveInteger
EuclideanRotation
                <- SignedInteger
RatchetSuffix   <- '*' PositiveInteger
ProbabilitySuffix
                <- '?' Probability?
ReplicationSuffix
                <- '!' PositiveInteger

Attributes      <- '{' H Attribute (H ',' H Attribute)* H '}'
Attribute       <- Identifier H '=' H ScalarValue

ScalarPattern   <- RightAligned / LeftAligned / FreePattern
RightAligned    <- Ellipsis S ScalarTerm (S ScalarTerm)*
LeftAligned     <- ScalarTerm (S ScalarTerm)* S Ellipsis
FreePattern     <- ScalarTerm (S ScalarTerm)*
ScalarTerm      <- (ScalarValue / Default) ReplicationSuffix?
Ellipsis        <- '...'
Default         <- '.'

ScalarLaneName  <- 'octave' / 'velocity' / 'vel' / 'accent' /
                   'duration' / 'dur' / 'gate' / 'slide' / 'ratchet' /
                   'offset'
CvLaneName      <- 'cv' PositiveInteger
SettingName     <- 'cycle' / 'tonic' / 'scale' / 'glide'
SettingValue    <- ScalarValue / PitchValue / Identifier

ScalarValue     <- Number Unit?
Unit            <- 'ms'
Probability     <- Number
Number          <- SignedInteger '/' PositiveInteger / Decimal
Decimal         <- Sign? ([0-9]+ ('.' [0-9]+)? / '.' [0-9]+)
SignedInteger   <- Sign? UnsignedInteger
PositiveInteger <- [1-9] [0-9]*
UnsignedInteger <- [0-9]+
Sign            <- '+' / '-'
Identifier      <- [A-Za-z_] [A-Za-z0-9_]*

BlankLine       <- H Comment? Newline
Trailing        <- H Comment?
Comment         <- '//' (!Newline .)*
H               <- [ \t]*
S               <- [ \t]+
LineEnd         <- Newline / End
Newline         <- '\r\n' / '\n' / '\r'
End             <- !.
```

The PEG chooses `...` before decimal/default productions, `^^` before `^`, and
the longest duration/chord/Roman token before its prefix. Ordered choice also
tries `ChordValue` before `PitchValue`, `RandomChoice` before a plain
subdivision, and a decimal before the scalar `.` default. Keyword productions
require the shown whitespace or delimiter so identifiers such as `playhead` and
`cv1_shape` are not split. Lane-specific productions make `x1` a typed event,
`.1` a typed scalar, and `cv2 = ...` a future top-level assignment without
semantic token splitting.

## Prototype replacement and compiler boundary

This grammar replaces the earlier prototype language outright. In particular,
the `articulation`/`art` lane no longer exists: articulation belongs to note
events. There is one parser, one typed tree, and one semantic lowering path.
Prototype patches and fixtures in this repository are updated with the grammar;
there is deliberately no legacy parser or source-version dispatch while the
language remains experimental.

Saved source still carries a language-version number. This replacement bumps
that number, but the version is currently only an incompatibility marker: the
binary contains one grammar and one compiler. A future incompatible change may
introduce an explicit migration at this boundary; it must not silently revive
or accumulate old parser implementations.

cpp-peglib produces document, expression, group, event, chord, suffix,
attribute, scalar, CV, and transform nodes with source spans. Semantic lowering
resolves names, types, units, branch cardinality, operator compatibility, and
bounded event density. Prepared immutable graphs alone cross to the audio
thread. Parsing or semantic failure leaves the currently valid program playing;
`stop` takes effect only when it belongs to a successfully compiled program.

No second delimiter or expression parser follows cpp-peglib. It owns statement,
lane, group, choice, chord/slash, suffix, attribute, unit, and pipeline
boundaries. Semantic lowering converts the already validated captured leaf
lexemes to numbers, registers, chord intervals, and enum values; it never tries
to rediscover nesting or operator scope from source text. This is the
implementability criterion for the cpp-peglib pass.

The lowering order is fixed:

1. resolve document names, settings, lane aliases, and transform signatures;
2. type note/group/chord/CV nodes while retaining nested structure and source
   paths;
3. apply structural and pitch transformations to that typed graph;
4. compute branch cardinalities and bind edge-aligned lanes to resulting
   positions;
5. attach free-running lane readers and perform predecessor dataflow for ties,
   slides, and CV knots;
6. derive checked event-density/lookahead requirements and allocate exact
   workspaces on the compilation thread; and
7. publish only the immutable prepared graph.

Replication and Euclidean counts use overflow-checked size arithmetic. The
compiler may allocate according to the successfully checked graph while off the
audio thread; it does not impose small musical maxima merely to avoid runtime
allocation. Choices remain graph branches rather than a Cartesian flattening.
Allocation failure, unrepresentable counts, or an event schedule too dense for
addressable storage is a compile diagnostic, leaving the old program active.

## Minimum conformance cases

Parser/compiler tests must cover at least these valid distinctions:

```text
notes 1 x2 ^3 ^^4 >5 >^6 ~ 7_3 8. 9!2 10*3 11?.2
notes 1 _!2 [2_ 3]_2 [4 5 | 6 [7 8]] <1 [2 3]>
notes D@4 Eb' I i ii7 Bbm7b5@3 / D@2 (C E G) / B
gate ... .1 . .2
cv1 0 . . 5 |> interp smooth
```

They must also reject each of these for the stated reason:

| Source fragment | Diagnostic class |
| --- | --- |
| `^>1`, `1!3?.5` | noncanonical operator order |
| `>x1`, `x1 _`, `x1 >2` | incompatible ghost/legato articulation |
| `1*3{slide=1}`, `>1*3` | unused or contradictory slide/ratchet |
| `~?.5`, `_` at a pattern start | meaningless probability or missing tie source |
| `VIII`, `Iv`, `Im7` | invalid Roman degree, case, or contradictory quality |
| `C4`, `D / F#` | ambiguous octave spelling or slash bass on a note |
| `.2 ... .1`, `... .1 ...` | illegal ellipsis position/count |
| `[1 | 2 3]` with an aligned lane | unequal branch cardinality |
| `cv0 1`, `cv3 1` with two outputs | invalid or unavailable CV output |

Golden runtime timelines must separately demonstrate attack, rest, long span,
tie, slide, replication, ratchet, subdivision weights, probability miss,
arrangement clipping, CV step, and all three continuous CV curves. These are
semantic tests, not snapshots of parser internals.

## Clash audit

- `@` means absolute register only.
- `'` and `,` are relative-register suffixes only in note mini-notation.
- `^` is articulation only as a note prefix.
- `x` is ghost only when attached as a note prefix; standalone rhythm/drum `x`
  can remain a hit.
- `>` is a slide only as a note prefix; it is not a scalar-lane value.
- `_` is a tie as a complete note element and duration only as a suffix attached
  to a note or group. Both cases are decided by the PEG position, not token
  post-processing.
- `!` is replication only as a suffix.
- `*` is note ratcheting in mini-notation and arrangement repetition at the
  expression level; context distinguishes them.
- `/` remains fractions in numerical values and slash bass after chords; it is
  not a duration suffix.
- parentheses start voicings before a pitch and Euclidean arguments after a
  complete event.
- braces after an event are attributes; sequence-body braces belong to the
  document grammar.
- `.` is no-op/decimal in scalar lanes and dotted duration only after an event.
- `cvN` without `=` is a sequence lane; `cvN =` is reserved for a future
  top-level signal definition.

This leaves a coherent path to chords, drums, interpreted MIDI, CV expressions,
and an aligned editor view without changing core event punctuation.
