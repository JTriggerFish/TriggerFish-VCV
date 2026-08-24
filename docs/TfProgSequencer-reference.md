# Prog Sequencer reference

Prog Sequencer is an externally clocked, polyphonic live-coding sequencer for
VCV Rack. A program describes one or more reusable sequences, combines them
into an arrangement, and selects what to play. Pitch, Gate, Trigger, Velocity,
and Accent share up to 16 polyphonic channels; CV1, CV2, and CV3 are
monophonic.

Prog Sequencer is currently in beta. Its language and feature set may change
as musical workflows are refined; check this reference when upgrading.

## Quick start

```text
riff = sequence {
  subdiv 8
  tonic D@3
  scale dorian
  notes 1 2 ^3{stacc} 4 5{quiet} >6 7{ten} _
  velocity .72 .58 .82
}

answer = riff |> modulate_degree 4 |> octave 1
song = riff * 2 + answer

seed 42
play song
```

Each incoming Clock advances one beat. `subdiv 8` makes an unsuffixed note an
eighth note (`1/2` beat), so the eight-step Notes pass lasts four incoming-clock
beats. Short auxiliary lanes loop within that pass and restart with it.

Use Ctrl+`.` to evaluate the complete document. Ctrl+Enter evaluates the
top-level statement containing the selection or cursor, leaving other edited
statements inactive. A successful change takes effect on the next Clock while
preserving the phase of sequences whose names still exist. The active
arrangement term follows that named sequence across insertions and reordering;
if the active term was removed, the nearest remaining structural position
starts cleanly on the activation clock. If an edit has an error, the last valid
program keeps playing.

`//` starts a comment. `stop` stops the arrangement without deleting the
program.

## Program structure

A sequence has one `notes` lane, optional settings, and optional numerical or
CV lanes:

```text
name = sequence {
  subdiv 16
  tonic C@4
  scale minor
  glide 1/4
  notes 1 3 5 7
  octave 0 1
  velocity .7 .55 .85
  gate .8
  cv1 0 5 2
}
```

Settings do not cycle:

| Setting | Meaning | Default |
| --- | --- | --- |
| `subdiv` | denominator of an unsuffixed note (`4/N` beats) | `4` |
| `tonic` | tonal centre and optional default octave | `C@4` |
| `scale` | scale used by degrees and Roman chords | `major` |
| `glide` | default slide time in beats | `1/4` |

One complete top-level Notes pass is the sequence boundary. Its actual elapsed
length follows the default subdivision, explicit duration suffixes, the
Duration lane, and timing transforms. At the boundary every auxiliary lane
restarts, and an arrangement advances after the requested number of complete
passes. A bracket group such as `[1 2]` subdivides one step and therefore still
occupies one top-level position.

## Notes and pitch

### Scale degrees

Numbers are degrees in the active scale:

```text
scale minor
notes 1 2 3 4 5 6 7 8 b3 #4
```

In a minor scale, `3` is already the minor third; `b3` lowers it by another
semitone. Degrees continue across scale octaves. In a seven-note scale, `8` is
the next tonic; in a pentatonic scale it is `6`; in an octatonic scale it is
`9`. `0` is the scale member immediately below the tonic and negative degrees
continue downward.

Available scales are:

- `major`, `ionian`
- `minor`, `aeolian`, `harmonic_minor`
- `dorian`, `phrygian`, `lydian`, `mixolydian`, `locrian`
- `major_pentatonic`, `pentatonic`, `minor_pentatonic`
- `octatonic_whole_half`, `whole_half`
- `octatonic_half_whole`, `half_whole`

### Named pitches and registers

Named pitches are chromatic and do not snap to the scale:

```text
notes C D# Eb F#
```

`@` sets an absolute octave. Apostrophes and commas shift relative to the
sequence octave:

```text
notes D@4 1' 1'' 1, Cm7' C@-1
```

`D@4` is the single note D4. `D7` is a dominant-seventh chord; `D7@3` is that
chord rooted in octave 3.

### Chords and voicings

Parentheses create a simultaneous voicing. Jazz and Roman-degree symbols are
also accepted:

```text
notes (1 b3 5) (C E G)@4 (C@3 E@4 G@4)
notes Cm7 D7 Cmaj9 C7#9 Bbm7b5@3 / D@2
notes I i iim7 bVII V7
```

Jazz chords support major, minor, diminished, augmented, suspended, power,
sixth, seventh, ninth, eleventh, and thirteenth forms, plus `add` and `b`/`#`
alterations. Roman roots are limited to `I` through `VII`; uppercase implies a
major triad and lowercase a minor triad. Roman chords follow the active scale.

A slash bass is emitted first and placed below the chord unless it has an
explicit octave. Chords and voicings use matching Pitch, Gate, Trigger,
Velocity, and Accent channels, up to Rack's 16-channel limit.

## Rhythm and articulation

Common note forms are:

| Form | Meaning |
| --- | --- |
| `1` | ordinary attack |
| `x1` | short, quiet ghost attack |
| `^1`, `^^1` | accent and strong accent |
| `1{quiet}` | quiet attack with normal Gate |
| `1{stacc}` | staccato attack with a quarter-span Gate |
| `1{ten}` | tenuto attack with a 95 percent Gate |
| `>1`, `>^1` | slide to the target, optionally accented |
| `~` | rest |
| `_` | extend the preceding event by one base slot |
| `1_`, `1_3` | event lasting two or three base slots |
| `1.`, `1..` | dotted or double-dotted duration |
| `1!3` | repeat as three logical events |
| `1*3` | three ratchets within one event |
| `1?.25` | play with 25 percent probability |
| `1(3,8,1)` | three hits in eight cells, rotated one cell later |

Velocity is a continuous `0..1` performance intensity. Accent is separate
note intent which can engage a dedicated instrument response. `^` and `^^`
assert Accent and raise Velocity to `.88` and `1` respectively, so the common
case still needs only one mark. There is no numerical Accent lane. A destination
that does not use the Accent output can simply use Velocity.

`{quiet}` halves the resolved Velocity without changing Gate or Accent, so
`^1{quiet}` is a quieter note which retains its timbral Accent output. A ghost
is deliberately a stronger compound gesture: `x1` halves Velocity, suppresses
Accent, uses at most 10 percent of the event span and 20 ms of Gate, and cannot
feed a tie or slide. `quiet`, `stacc`, and `ten` are therefore rejected on a
ghost when redundant or contradictory.

A bare probability suffix, as in `1?`, means `.5`. Probability is decided
once per logical event, so all of its ratchets either sound or rest together.
A missed event still occupies time and advances its numerical lanes.

Suffixes have one accepted order:

```text
[slide] [dynamic] pitch [register] [duration]
  [Euclidean] [ratchet] [probability] [replication] [attributes]
```

For example, `>^1'_3*2?.8!4` creates four three-slot upper-tonic slides; each
has two ratchets and an 80 percent chance of sounding.

### Groups and choices

Square brackets subdivide one parent span:

```text
notes [1 2] [3 4 5] [1_ 5]
```

The first group plays two halves, the second three thirds, and `[1_ 5]` divides
its span in a 2:1 ratio. A suffix on a group applies to the complete group:

```text
notes [1 3 5]_2 [2 4](3,8) [6 7]?.5 [1 8]!2
```

Angle brackets alternate once per sequence pass; a vertical bar makes a seeded
random choice:

```text
notes <1 3 5 7> [1|3|5] [(1 3 5)|(2 4 6)]
```

Atomic pitch and chord choices are available now. Multi-event branches and
recursively nested alternatives are listed under planned features below.

### Ties, slides, and gates

`_` holds the previous pitch and Gate without retriggering. `>1` glides from
the preceding sounding pitch; after a rest or probability miss it becomes a
normal attack. Ties and slides do not cross a random/alternate branch,
arrangement boundary, or sequence wrap. Chord slides require the same known
voice count on both sides.

`gate .5` holds Gate for half an event. `{stacc}` and `{ten}` provide concise
Gate presets of `.25` and `.95` without changing the event span. They combine
with ordinary notes, accents, quiet notes, slides, and ratchets. On a slide they
control the target's eventual release; the source still remains connected into
the slide. Ratchets apply the selected Gate fraction to every subspan.

Only `_` is a true tie. `{stacc,ten}` is contradictory and rejected.

### Event attributes

Attributes travel with a note under transforms. Bare names are expressive
flags; attributes with values provide exact local controls:

```text
notes 1{quiet} ^3{stacc} 5{ten} >6{slide=80ms,ten} 7{len=5/4}
notes 1{vel=.95} 3{gate=12ms}
```

| Attribute | Meaning |
| --- | --- |
| `quiet` | halve resolved Velocity; Gate and Accent are unchanged |
| `stacc` | Gate preset `.25` |
| `ten` | Gate preset `.95` |
| `vel=VALUE` | note-local base Velocity from `0` to `1`; dynamic modifiers still apply |
| `gate=VALUE` | Gate fraction `0..1` or non-negative `ms` |
| `slide=VALUE` | slide time in beats or `ms`; valid only on `>` |
| `len=VALUE` | positive exact event span in beats |

Gate resolution follows musical specificity: the normal `.8` or named
`stacc`/`ten` preset is the fallback; an explicit Gate-lane value overrides it;
`.` in that lane retains it; and inline `gate=VALUE` wins last.

## Random notes and values

Random results are repeatable. `seed NUMBER` changes the variation while the
same seed and program produce the same event sequence.

Pitch random forms are:

| Form | Meaning in `notes` |
| --- | --- |
| `$` | uniform choice from one octave of the active scale |
| `${1,8}` | uniform scale degree from inclusive bounds |
| `$u{1,8}` | explicit spelling of `${1,8}` |
| `$n{4,1.5}` | normal distribution in scale-degree space |
| `$c{0,11}` | uniform chromatic semitone offset from the tonic |
| `$cn{6,2}` | normal chromatic offset from the tonic |

Scale forms always resolve to a note in the active scale. Normal pitch values
are rounded to the nearest degree. Chromatic forms bypass scale quantization
and round to the nearest semitone. Uniform pitch bounds must be ordered
integers and are inclusive; a normal standard deviation must be positive.

Use `[1|3|5]` when the desired choices are an explicit set rather than a range.

Numerical and CV lanes use continuous uniform or normal values:

```text
velocity $u{.45,.9}
velocity $n{.7,.12}
gate     $u{10ms,25ms}
octave   $u{-1,1}
ratchet  $n{3,.6}
cv1      $u{-5,5}
```

Arguments must have matching units. Octave and Ratchet draws are integers;
normal draws are rounded. Normal results are kept inside the lane's legal
range, such as `0..1` for Velocity and fractional Gate. Normal draws
are limited to four standard deviations from the mean, keeping realtime
scheduling and event capacity finite and predictable.

One pitch and one set of pitched-event controls are drawn per logical event;
ratchets share those results. A random CV token is one repeatable control
point, not audio-rate noise.

## Numerical lanes

Numerical lanes loop independently unless aligned with `...`:

```text
notes    1 2 3 4 5 6 7 8
octave   0 1 0
velocity .72 .55
duration 1 1 3/2
gate     .8 .4 .6
slide    . 80ms
ratchet  1 1 2
offset   -10ms!2 +8ms
```

| Lane | Accepted values | Normal default | Advances on |
| --- | --- | --- | --- |
| `octave` | signed integer octave offset | `0` | pitched events |
| `velocity`, `vel` | `0..1` | `.72` | pitched events |
| `duration`, `dur` | positive beats or fraction | `1` | notes and rests |
| `gate` | `0..1` or non-negative `ms` | `.8` | pitched events |
| `slide` | non-negative beats/fraction or `ms` | `glide` | pitched events |
| `ratchet` | positive integer | `1` | pitched events |
| `offset` | signed beats/fraction or `ms` | on-grid | score time |

`.` consumes a lane position while inheriting the normal value. `.1` is the
number one tenth. `.!3` repeats the default across three positions.

An accent prefix raises Velocity to its built-in accent floor and drives the
separate Accent output for the effective Gate. A Slide value is read for every
pitched event but matters only when the note uses `>`.

### Aligning a lane to Notes

A leading `...` right-aligns explicit values with one structural cycle of the
Notes lane; a trailing `...` left-aligns them:

```text
notes 1 2 3 4 5 6 7 8
gate  ... .1!2          // short gates on notes 7 and 8

notes 1 2 3 4 5 6 7 8
gate  .1!2 ...          // short gates on notes 1 and 2
```

`...` must be the first or last term and may appear only once. Replication and
Euclidean cells count as positions; ratchets and elongation do not. Random and
alternate branches must have equal cell counts when another lane is aligned to
them. An aligned lane cannot also use `rate`.

## CV1, CV2, and CV3

CV lanes use volts and follow score time, including rests, ties, and
probability misses. Ratchets share their parent's CV sample.

```text
notes 1 2 3 4
cv1   0 5 0 |> interp smooth
cv2   0 . . 5 |> interp linear
cv3   5 2 0 |> interp power 2
```

Interpolation modes are:

- `step` — sample and hold; the default;
- `linear` — constant slope;
- `smooth` — smoothstep curve; and
- `power P` — curve `a + (b-a)t^P`, where `P > 0`.

On a stepped CV lane, `.` holds the preceding value. On an interpolated lane,
`.` leaves no control point, so the curve continues to the next explicit
value. Leading defaults look backward through the repeating lane.

Free CV lanes accept `|> rate R`. Edge-aligned CV supports `step` only.

## Transforms

Pipelines transform a lane, a sequence, or an arrangement:

```text
notes 1 2 3 4 |> rotate 1 |> every 4 rev
offset -10ms!2 +8ms |> rate 1/2

variation = riff
  |> shift_degree 2
  |> sometimes .25 (rotate 1)
```

Transform scope is explicit:

| Target | Available transforms |
| --- | --- |
| Notes lane | `rev`, `rotate`, and pitch transforms |
| Numerical lane | `rev`, `rotate`, `rate` |
| CV lane | numerical-lane transforms plus `interp` |
| Complete sequence/arrangement | `rev`, `rotate`, pitch, scale, and timing transforms |

`rate` and `interp` are lane-local. `fast`, `slow`, `swing`, `early`, and
`late` apply to a complete sequence or arrangement.

### Structural and pitch transforms

| Transform | Effect |
| --- | --- |
| `rev`, `reverse` | reverse order |
| `rotate N` | rotate by signed positions |
| `shift_degree N` | move notes by scale indices, staying in key |
| `modulate_degree N` | move the tonal centre by a scale degree |
| `transpose_semitone N` | transpose chromatically |
| `octave N`, `transpose_octave N` | transpose by octaves |
| `scale NAME` | change the active scale |

`shift_degree 1` moves the tonic to degree 2. `modulate_degree 1` is unison,
`modulate_degree 4` moves the centre to the fourth, and negative values move
downward. The operations are deliberately different: shifting stays in the
current key, while modulation preserves the riff's chromatic shape around a
new tonal centre.

Pipelines run left to right, so these need not sound the same:

```text
riff |> modulate_degree 3 |> scale major
riff |> scale major |> modulate_degree 3
```

### Timing transforms

| Transform | Effect |
| --- | --- |
| `fast F` | divide event time by positive factor `F` |
| `slow F` | multiply event time by positive factor `F` |
| `swing R` | swing each event's own subdivision |
| `swing R GRID` | swing an explicit positive clock grid |
| `early A`, `late A` | fixed beat/fraction or `ms` displacement |
| `early random A`, `late random A` | seeded displacement from zero to `A` |

Swing `.5` is straight; larger values lengthen the first member of each grid
pair. Timing transforms move events without changing the arrangement window.

`rate R` is lane-local: it changes the phase speed of Offset or another
free-running numerical/CV lane without speeding up the Notes lane.

### Conditional transforms

```text
riff |> every 4 rev
riff |> sometimes .3 (rotate 1)
```

`every N` applies a transform on each positive Nth sequence pass.
`sometimes P` applies it with seeded probability `0..1`.

## Arrangements

`+` concatenates sections and `*` repeats them:

```text
acid = sequence {
  tonic D#@2
  scale minor
  notes ^1 1!3 ^5 7 ^1 ~ 1 ~ ^8 >1, ~ 1 >^1, >1
  velocity .5
  gate .5
  glide .8
}

iv = acid |> modulate_degree 4 |> octave -1
v  = acid |> modulate_degree 5 |> octave -1 |> scale major

song = acid * 8 + iv * 4 + v * 4
play song
```

Parentheses group arrangement expressions:

```text
song = (verse * 2 + fill) * 2 + ending
```

A section boundary occurs only after the final Notes step, including its
ratchets or subdivisions, has completed. The new section owns the next attack;
ties and slides do not cross between different sections.

## Musical examples

### Seeded generative melody

```text
melody = sequence {
  subdiv 16
  tonic A@3
  scale minor_pentatonic
  notes $u{1,10}(5,8) $n{5,1.25}(3,8,2)
  octave 0 1 0
  velocity $n{.72,.12}
  gate $u{.35,.8}
  cv1 0 2 5 3 |> interp smooth
}

seed 2026
play melody
```

The two Euclidean figures fill sixteen cells. Every generated pitch belongs to
A minor pentatonic, while the three-value Octave lane and independent controls
create a longer variation.

### Polymetric chord movement

```text
changes = sequence {
  tonic C@3
  scale major
  notes iim7 V7 Imaj7 VI7
  duration 2
  velocity .62 .78 .7
  gate .9
  cv1 0 3 5 2 |> interp power 2
}

lift = changes |> modulate_degree 4
song = changes * 2 + lift * 2
play song
```

Each pass contains four two-beat chords. The three-step Velocity lane moves
against the harmony within a pass and restarts when the four-chord Notes phrase
restarts; the second half reuses the phrase around the fourth.

### Humanized timing

```text
groove = sequence {
  notes [1 5] 3 [4 6 8] 5
  duration 1 1 2
  offset -7ms 4ms 0 |> rate 1/2
  ratchet 1 1 2
}
|> swing .58 1/8
|> early random 3ms

seed 99
play groove
```

Subdivision, Duration, Offset, and Ratchet each have independent phrases.
Swing supplies the main groove while the seeded early displacement adds a
repeatable small variation.

## Current limits

The following are not available yet:

- multi-event or recursively nested random/alternate branches;
- continuous interpolation on `...`-aligned CV lanes;
- CV outputs beyond CV1, CV2, and CV3;
- MIDI input bindings and performance interpreters for bass, piano, and
  arpeggiation;
- scale-stacked chord generation; and
- expressions that derive CV signals from MIDI, other CVs, or functions.

Until those features arrive, use atomic alternatives, free-running CV for
continuous curves, explicit chord symbols/voicings, and the three physical CV
outputs.
