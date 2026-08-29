# Prog Sequencer reference

Prog Sequencer is an externally clocked, polyphonic live-coding sequencer for
VCV Rack. A program describes one or more reusable sequences, combines them
into an arrangement, and selects what to play. Pitch, Gate, Trigger, Velocity,
and Accent share up to 16 polyphonic channels; CV1, CV2, and CV3 are
monophonic.

The language takes inspiration from live-coding systems including TidalCycles
and Gibber, together with pattern languages, trackers, and step sequencers. Its
goal is fast musical iteration: a short line should produce a useful phrase,
and structural variations should be possible without copying or manually
editing every step. The syntax is designed for improvising, sequencing, and
reshaping material while it plays.

Prog Sequencer is currently in beta. Its syntax, semantics, and feature set are
subject to change as musical workflows are refined; check this reference when
upgrading.

## Contents

- [Quick start](#quick-start)
- [Editing and execution](#editing-and-execution)
  - [Live activation and phase](#live-activation-and-phase)
- [Clock and shared transport](#clock-and-shared-transport)
- [Musical time values](#musical-time-values)
- [Program structure](#program-structure)
- [Notes and pitch](#notes-and-pitch)
  - [Scale degrees](#scale-degrees)
  - [Named pitches and registers](#named-pitches-and-registers)
  - [Chords and voicings](#chords-and-voicings)
- [Rhythm and articulation](#rhythm-and-articulation)
  - [Direct events](#direct-events)
  - [Reusable rhythmic gestures](#reusable-rhythmic-gestures)
  - [Slides with a gesture](#slides-with-a-gesture)
  - [Probability: silence and omission](#probability-silence-and-omission)
  - [Groups and choices](#groups-and-choices)
  - [Ties, slides, and gates](#ties-slides-and-gates)
  - [Event attributes](#event-attributes)
- [Random notes and values](#random-notes-and-values)
- [Numerical lanes](#numerical-lanes)
  - [Aligning a lane to events](#aligning-a-lane-to-events)
- [CV1, CV2, and CV3](#cv1-cv2-and-cv3)
  - [CV envelopes](#cv-envelopes)
- [Transforms](#transforms)
  - [Structural and pitch transforms](#structural-and-pitch-transforms)
  - [Timing transforms](#timing-transforms)
  - [Conditional transforms](#conditional-transforms)
- [Arrangements](#arrangements)
- [Musical examples](#musical-examples)
  - [Expressive bass line](#expressive-bass-line)
  - [Song with distinct sections](#song-with-distinct-sections)
  - [Seeded generative melody](#seeded-generative-melody)
  - [One intentional CV curve](#one-intentional-cv-curve)
  - [Humanized timing](#humanized-timing)
- [Current limits](#current-limits)
- [Formal grammar reference](#formal-grammar-reference)

## Quick start

```text
riff = sequence {
  subdiv 8n
  tonic D@3
  scale dorian
  notes 1 2 ^3{stacc} 4 5{quiet} >6 7{ten} _
}

answer = riff |> shift_degree 3 |> octave 1
song = riff * 2 + answer

seed 42
play song
```

The Clock input uses a fixed 24 PPQN transport signal. Twenty-four pulses
advance one quarter-note beat. `subdiv 8n` makes an unsuffixed note an eighth
note (`1/2` beat), so the eight-step Notes pass lasts four quarter-note beats.
Ordinary auxiliary lanes keep their own cycle lengths across repeated Notes
passes. Use `...` to attach values to the beginning or end, or both edges of
each Notes pass.

The essential live controls are:

- Ctrl+Enter evaluates the selected source or the top-level statement
  containing the cursor;
- Ctrl+Shift+Enter evaluates the same source on the next scheduler step when
  its playback representation is compatible;
- Ctrl+`.` evaluates the complete document while preserving playback phase;
- Ctrl+Shift+`.` evaluates the complete document and restarts the arrangement
  on the next quarter beat;
- Ctrl+Space mutes the module immediately while its score keeps running, or
  queues its unmute for the next shared quarter beat; and
- Ctrl+R restarts only this module on the next shared quarter beat without
  recompiling its draft.

When RUN is connected directly to a TriggerFish Transport, Ctrl+Shift+Space
pauses or plays that Transport, Ctrl+Shift+R restarts it from beat zero, and
Ctrl+Shift+Backspace stops it at beat zero. These commands control every
sequencer receiving the same Transport outputs.

Ctrl+Z undoes program edits. Ctrl+Shift+Z and Ctrl+Y redo them.

A successful ordinary evaluation takes effect on the next quarter beat and
preserves playback phase whenever the old and new representations are
compatible. Structural changes have a deliberately stricter restart policy,
described under [Live activation and phase](#live-activation-and-phase). If an
edit has an error, the last valid program keeps playing.

During playback, the active note, rest, or tie remains softly illuminated and
an underline advances across it for its complete duration. The RUN light is dim
while the transport is enabled and flashes brightly on each quarter beat,
so sustained or slow phrases still show that time is advancing.

`//` starts a comment.

## Editing and execution

Editing changes the draft shown in the module. Playback changes only after an
evaluation shortcut succeeds; commenting, duplication, and ordinary typing do
not alter the active program.

| Shortcut | Action |
| --- | --- |
| Ctrl+Enter | Evaluate the selection, or the containing top-level statement when there is no selection, on the next quarter beat |
| Ctrl+Shift+Enter | Evaluate the same source on the next scheduler step when phase can be preserved |
| Ctrl+`.` | Evaluate the complete document and preserve the phase of sequences whose names remain active |
| Ctrl+Shift+`.` | Evaluate the complete document and restart its arrangement from the beginning on the next quarter beat |
| Ctrl+Space | Mute this module immediately while its score keeps running, or unmute it on the next master quarter beat |
| Ctrl+R | Restart only this module from its active arrangement on the next master quarter beat |
| Ctrl+Shift+Space | Pause or play the TriggerFish Transport connected directly to RUN |
| Ctrl+Shift+R | Restart the connected TriggerFish Transport from beat zero |
| Ctrl+Shift+Backspace | Stop the connected TriggerFish Transport and return it to beat zero |
| Ctrl+`/` | Toggle line comments on the selected lines or the line containing the cursor |
| Ctrl+D | Duplicate the selection or, when there is no selection, the current line |
| Ctrl+Z | Undo the last program edit |
| Ctrl+Shift+Z or Ctrl+Y | Redo the last undone program edit |

### Live activation and phase

Evaluation has two independent concerns: the boundary at which compiled source
becomes active, and the playback state carried into it. Compilation itself is
immediate, but a successful program remains pending until its musical
activation boundary. Typing and compiling therefore never replace code halfway
through an audio sample or an existing scheduler step.

Ctrl+Enter is the normal live-evaluation command. It activates on the strictly
following quarter-note beat and preserves the active arrangement term,
repetition count, sequence cursor, independent lane phases, random cycle, and
transform phase. Ctrl+`.` applies the same quarter-beat policy to the complete
document.

Ctrl+Shift+Enter is the faster alternative for a selected statement or source
region. It activates on the next scheduler-step boundary while preserving the
same state. A scheduler step is the next boundary produced by the currently
active event scheduler, not merely a division of the clock grid. It therefore
respects duration suffixes, the Duration lane, omissions, and—in separated
mode—the rhythm gesture's own subdivision.

Direct pitched events and a separated held-pitch/rhythm gesture use different
playback representations. There is no reliable one-to-one cursor conversion
between them. If an edit changes the active sequence from one representation
to the other, its activation is promoted to the next quarter beat regardless
of shortcut. At that beat:

- the current arrangement term and its repetition count are retained;
- the current sequence pass restarts from its beginning;
- its pitch, rhythm, numerical, and CV lane states restart together;
- stale Gates and Slides are cleared; and
- the shared Transport and all sibling sequencers continue untouched.

Removing the active term is also incompatible: the nearest remaining
arrangement position begins cleanly on the activation beat. Ctrl+Shift+`.` is
the explicit stronger operation; it evaluates the complete document and
restarts the entire arrangement on the next quarter beat. If playback is
paused or stopped, a pending change waits until score time advances to its
requested boundary.

The local transport starts enabled. It is combined with the RUN input, so both
must permit playback. The RUN jack is normalled high, and an unpatched module
therefore follows the local transport. Ctrl+Space mutes locally and
immediately: Gate, Trigger, and Accent go low, but the score, note and lane
cursors, random state, slides, CV envelopes, and graphical beat indicators keep
following the shared Clock. Pressing Ctrl+Space again queues unmute on the next
true master quarter-note boundary; another press before the boundary cancels
the queued resume. This keeps the module in musical phase with its siblings
throughout the mute. The status strip turns red and reads `MUTED`; a queued
unmute turns it amber and reads `MUTED / UNMUTE NEXT BEAT`, while the sequence
cursor continues to show the current phase. A connected RUN voltage below 1 V
remains a genuine score pause, and raising it continues from the retained
position. RESET is the separate operation that returns the arrangement and all
lanes to their beginnings. A shared RESET also clears a local Ctrl+Space mute or queued resume,
so a Transport Restart reliably brings every connected sequencer back into
playback. Ctrl+R instead leaves the shared Transport and the other sequencers
untouched; it uses the sequencer's current musical clock phase to apply that
local restart on the next quarter-note boundary.

The transport shortcuts leave the program source unchanged. A connected
Transport command is handed to the Transport's audio processing through an
atomic mailbox. Timing then reaches every sequencer through the ordinary
Clock, Run, and Reset cables. If RUN has no direct connection from a
TriggerFish Transport, a connected-Transport shortcut reports this in the
status area and leaves the patch unchanged. The shortcut requires a direct RUN
cable, which identifies one unambiguous Transport even in a larger patch.

Double-click selects a word and triple-click selects a complete row. The usual
text-field selection, clipboard, and navigation controls remain available.
The module context menu links directly to this reference. Its Examples submenu
can load an Acid bassline, Slow bassline, or Descending arpeggio program into
the editor. Loading replaces the current draft and participates in Rack undo;
the active program continues until the example is evaluated. The Editor
submenu contains the width and heatmap controls. Heatmaps include the
perceptually uniform Magma, Inferno, Plasma, Viridis, and Cividis ramps, plus
black-based CRT Green, CRT Blue, CRT Yellow, and CRT Red phosphor-style ramps.
The heatmap choice is stored with the patch.

Each active `cv1`, `cv2`, or `cv3` line displays a six-second scrolling trace
of that lane's output in the unused space to the right of its text. The fixed
vertical range is -5 V to +5 V. Positive-only traces place 0 V at the bottom,
negative-only traces place it at the top, and bipolar traces place it at the
centre. The line and its subtle phosphor glow use the selected heatmap. If the
program text leaves less than 42 pixels on that row, the trace is hidden rather
than drawn over the source; widen the module or shorten the line to reveal it.
The trace is capped at the width of 12 editor characters so it remains a
compact lane indicator.

## Clock and shared transport

Prog Sequencer interprets CLOCK as a fixed 24-pulse-per-quarter-note transport
stream. Score time is still written in quarter-note beats: `subdiv 4n` produces
quarter notes, `subdiv 8n` eighth notes, and `subdiv 16n` sixteenth notes. The
additional clock edges provide timing resolution and rapid period acquisition;
program speed continues to follow quarter-note score time. An external clock
source must therefore be configured for 24 PPQN.

The receiver uses direct edge locking. Each incoming pulse is an authoritative
phase point. The interval between adjacent pulses supplies a lightly smoothed
quarter-note period estimate for interpolation until the next pulse, and every
twenty-fourth pulse is the exact beat boundary. Phase is acquired on the first
pulse, the period is known on the second, and every later pulse applies an exact
phase correction.

For several sequencers, connect one TriggerFish Transport in parallel:

| Transport output | Prog Sequencer input | Purpose |
| --- | --- | --- |
| CLOCK | CLOCK | Shared fixed-rate timing and phase |
| RUN | RUN | Common pause and resume state |
| RESET | RESET | Common return to beat zero |

RESTART sends Reset before the first Clock and Run signals, giving every
sequencer time to clear its arrangement before beat zero arrives. PAUSE freezes
the master-clock phase and lowers Run. PLAY continues that phase. STOP lowers
Run, resets every connected sequencer, and leaves the master phase at zero.
A Clock cable carries timing edges. Run carries the immediate play/pause state,
and Reset marks the common beginning; all three connections form the complete
shared transport.

With the editor focused, Ctrl+Shift+Space controls Play/Pause on the Transport
connected directly to RUN. Ctrl+Shift+R requests Restart and
Ctrl+Shift+Backspace requests Stop. The keyboard request changes only the
connected Transport; its Clock, Run, and Reset outputs remain the source of
synchronization for the sequencers.

Changing Rack's audio driver or sample rate preserves both the learned clock
period and the elapsed fraction of the current pulse. This avoids a false
short or long interval on the first edge after an audio-device change.

## Musical time values

Tempo-relative durations have an explicit note-value form:

| Literal | Duration in quarter-note beats |
| --- | ---: |
| `1n` | 4; whole note |
| `2n` | 2; half note |
| `4n` | 1; quarter note |
| `8n` | 1/2; eighth note |
| `16n` | 1/4; sixteenth note |
| `32n` | 1/8; thirty-second note |

Append `d` for a dotted value or `t` for a triplet value. For example, `8nd`
is a dotted eighth lasting 3/4 beat, while `8nt` is an eighth-note triplet
lasting 1/3 beat. Denominators are powers of two from 1 through 128.

`subdiv` always takes a note value. Straight, dotted, and triplet steps are
therefore written directly:

```text
subdiv 16n   // sixteenth-note steps
subdiv 8nd   // dotted-eighth steps
subdiv 8nt   // eighth-note triplet steps
```

Other tempo-relative controls accept either note values or quarter-note beats.
Thus `16n` and `1/4` have the same duration, while the explicit note form often
reads more naturally in a musical program. Bare values retain beat arithmetic:
`1` is one quarter-note beat and `3/2` is one and a half beats.

Gate, Slide, Offset, `early`, and `late` also accept `ms`. Envelope segments
additionally accept seconds with `s`. A signed note value such as `-16n` is
valid where signed timing is meaningful, including Offset.

## Program structure

A sequence has one pitched lane, spelled `notes` or `chords`, plus optional
settings and numerical or CV lanes. The two spellings have identical syntax;
`chords` simply makes harmonic intent easier to see. Normally each pitched
event owns its attack. A reusable `rhythm` definition can instead supply the
attacks for a slower held-pitch timeline.

```text
name = sequence {
  subdiv 16n
  tonic C@4
  scale minor
  glide 16n
  notes 1 ^3{stacc} 5 >7 8{ten}
  cv1 0 5 2 |> interp smooth
}
```

A document is assembled from a small number of building blocks:

- a `sequence` definition creates an independently playable musical section;
- a `rhythm` definition creates a reusable unpitched gesture;
- an assignment such as `song = verse * 2 + chorus` names an arrangement or a
  transformed variation;
- `play NAME` selects the sequence or arrangement that should loop;
- `seed NUMBER` selects a repeatable set of random choices.

Definitions and assignments may appear in any order, but every name must be
unique. Names use letters, digits, and underscores and cannot begin with a
digit. Names such as `cv1`, `cv2`, and other `cv` followed by a non-zero number
are reserved for control-voltage lanes and future signal values.
Names which are themselves complete inline rhythm patterns are also reserved
for that syntax. In particular, a reusable rhythm cannot be named `x`, `_`,
`x_2`, or `x__`; descriptive names such as `comp`, `xbeat`, and `pulse_2` are
unambiguous.

Blank lines and `//` comments may appear between definitions or inside a
sequence, and a comment may follow any complete line. If an evaluated document
contains several `seed` commands, the last seed takes effect. The final `play`
command selects the active arrangement. Playback state is controlled locally
with Ctrl+Space or by the connected Transport.

Every sequence requires exactly one pitched lane: either `notes` or `chords`.
Settings and auxiliary lanes are optional, but the same lane cannot be declared
twice. `notes` and `chords` therefore cannot both occur in one sequence; a
single lane may freely mix individual notes and voicings. `vel` is a shorter
spelling of `velocity`, and `dur` is a shorter spelling of `duration`.

Settings do not cycle:

| Setting | Meaning | Default |
| --- | --- | --- |
| `subdiv` | note value of an unsuffixed note | `4n` |
| `tonic` | tonal centre and optional default octave | `C@4` |
| `scale` | scale used by degrees and Roman chords | `major` |
| `key` | written key used by `transpose_key` | none |
| `voicing` | automatic Jazz-chord voicing recipe | `basic` |
| `glide` | default slide time in beats or note values | `16n` |

Without a separate rhythm, one complete top-level pitched pass is the sequence
boundary. Its elapsed length follows `subdiv`, duration suffixes, the Duration
lane, `??` presence decisions, and timing transforms. With a separate rhythm,
the held pitched timeline defines the boundary and the shorter rhythm gesture
loops inside it. At the boundary the pitched pattern and gesture restart;
independently cycling numerical and CV lanes retain their normal phases while
the same sequence continues. An arrangement advances after the requested
number of complete pitched passes.

A bracket group such as `[1 2]` subdivides one direct event and therefore still
occupies one top-level position. When an arrangement moves to a different named
sequence, that sequence starts its lane cursors and phases from the beginning.

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
- `major_pentatonic`, `minor_pentatonic`
- `octatonic_whole_half`, `whole_half`
- `octatonic_half_whole`, `half_whole`

### Named pitches and registers

Named pitches are chromatic and do not snap to the scale:

```text
notes C D# Eb F#
```

`@` sets an absolute octave. Unmarked pitches use the octave from `tonic`, or
the current absolute value from an `octave` lane when that lane is present.
Apostrophes and commas shift relative to that octave:

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
In a `chords` lane only, a bare named root such as `C` or `F#` is its standard
major triad. The same spelling in `notes` remains one pitch; use `Cmaj:(1)` for
an explicit one-note factor inside a chord lane.

For a chord progression, prefer an explicit `key` and named chord roots:

```text
changes = sequence {
  key C
  voicing basic
  chords Dm9 G13 Cm9 F13 Bbm9 Eb13 Abm9 Db13
}
```

`key` records the progression's written transposition anchor; it does not
infer harmony, constrain chord roots, or make a chord diatonic. `Dm9` always
has D as its written root. `transpose_key D` on a sequence written with
`key C` shifts every realized pitch up two semitones. A target key chooses the
nearest chromatic interval (an exact tritone goes upward); octave transforms
remain explicit. `transpose_key` is a compile error without `key`. `tonic` and
`scale` remain the pitch system for scale degrees and Roman chords, so a fast
sequence of ii-V changes never needs a second synchronized centre lane.

#### Automatic voicing recipes

`voicing` accepts three recipes:

| Recipe | Behaviour |
| --- | --- |
| `basic` | ordinary stacked chord formula, followed by smooth inversions |
| `rootless_3notes` | three-note rootless guide-tone voicings |
| `rootless_4notes` | four-note Bill Evans-style A/B voicings |

The first `basic` chord uses the obvious root-position stack. Later automatic
chords choose inversions by shortest ordered semitone motion, then smallest
individual leap, stationary voices, register, and a deterministic recipe
order. Automatic voices never cross.

The engine tries to retain the previous note count when the next formula has
enough useful tones. It protects thirds, sevenths, written upper extensions,
and altered tones; the
`basic` recipe omits a perfect fifth first and then the root. For example, a
four-note `C7` followed by `F13` realizes the latter as `3 7 9 13`, rather than
abruptly growing to six voices. Rootless recipes omit the root by definition
unless an explicit factor override requests it, usually omit the fifth, and
fall back to a normal triad when no useful guide pair exists. Suspended seventh
chords use the suspension and seventh as their guides; sixth chords such as
`C6add9` use the third and sixth. A three-note `Dm11`, for example, retains its
third, seventh, and defining eleventh rather than silently becoming `Dm9`.

Recipes are contextual but deterministic. Repeated rhythm hits reuse the same
realized chord. Rests retain the last harmonic context, slash basses do not
enter the upper-voice calculation, and a continuing loop voice-leads across
its seam.

#### Explicit factors

Append a factor list to retain exact chord content while allowing the recipe
to choose its inversion and register:

```text
chords Cm9:(1 3 5 7 9) Cm9:(3)
chords G13:(3 7 9 13) C6add9:(1 3 6 9)
```

Factors inherit the written chord quality. Thus `3` in `Cm9:(3)` is the minor
third; the redundant spelling `b3` means the same thing. Omitted factors are
not emitted, so that example really is one note. A contradictory factor such
as `Cmaj7:(b3)`, or the duplicate `Cm9:(3 b3)`, is a compile error. The older
pitch voicing `(C@3 E@4 G@4)` remains the form for fixing exact pitches and
registers instead of chord factors.

`C7alt` creates separate altered candidates. Each contains exactly one of
`b9` or `#9` and exactly one of `#11` or `b13`; the closest legal candidate is
selected from the preceding voicing. It never automatically produces a
simultaneous `b9/#9` cluster. An explicit override can deliberately request an
unusual combination:

```text
chords C7alt C7alt:(3 7 b9 b13) C7alt:(3 b7 #9 #5)
```

`13` chords include 9 and 13 but not the normally avoided natural 11.
Conventional post-extension suspension spellings such as `G9sus`, `G7sus4`,
and `G7sus2` are accepted.

The recipe split follows the practical distinction between chord vocabulary
and voice-leading selection. The motion policy is informed by Dmitri
Tymoczko's discussion of short, non-crossing voice leading in
[*The Geometry of Musical Chords*](https://www.brainmusic.org/EducationalActivities/Tymoczko_chords2006.pdf);
crossing and overlap are kept as explicit constraints, as in the
[music21 voice-leading model](https://music21.org/music21docs/moduleReference/moduleVoiceLeading.html).
The rootless A/B factor sets follow the common 3-5-7-9 / 7-9-3-5 and
3-13-7-9 / 7-9-3-13 practice summarized in
[Piano Encyclopaedia's rootless-voicing reference](https://piano.org/theory/rootless-voicings/).

A slash bass is emitted first and placed below the chord unless it has an
explicit octave. Chords and voicings use matching Pitch, Gate, Trigger,
Velocity, and Accent channels, up to Rack's 16-channel limit.

`chords` is an optional readability alias for `notes`. It does not select a
different scheduler, so chords may be sequenced directly with all ordinary
event syntax:

```text
changes = sequence {
  subdiv 8n
  chords Cmaj9_2 [Dm11 G13] ^Am11 ~ Fmaj9*3
}
```

Conversely, `notes` may contain chords, and either spelling may mix single
pitches with voicings.

## Rhythm and articulation

### Direct events

In the compact form, pitch and rhythm stay together. This remains the best
form for a melody, bass line, or progression whose attacks do not repeat as a
separate idea:

```text
lead = sequence {
  subdiv 8n
  notes 1_2 ^3 [4 5] >6{ten} ~ 8*3
}
```

Each pitched token attacks at the beginning of its own span. All the familiar
duration, articulation, probability, subdivision, and ratchet syntax applies.

### Reusable rhythmic gestures

Define a gesture once when the same attack pattern should play different notes
or chords. A rhythm has its own subdivision and uses `x` for an unpitched hit:

```text
comp = rhythm {
  subdiv 16n
  events x_3 x ~!6 x_3 x ~!2
}
```

The gesture above lasts four beats. Applying it to a held chord timeline keeps
the harmonic durations separate and visible:

```text
keys = sequence {
  subdiv 2n
  chords (C F A)_2 (C D G)_2 (B, D G)_2 (A, C F) (G, C E)
  rhythm comp
  velocity .72 .8
  gate .95 .8
}
```

The sequence subdivision is a half note. The first three `_2` chords therefore
last four beats each, and the final two unsuffixed chords last two beats each:
`4 + 4 + 4 + 2 + 2 = 16` beats. The four-beat `comp` gesture loops four times
inside that phrase. At every `x`, the sequencer samples whichever pitched value
is active at that written musical time. Hits never advance the chord cursor,
so editing a rest, probability, or ratchet cannot move a later chord boundary.
A random or alternate pitched value is resolved once for its held span rather
than being redrawn by every hit. These explicit three-note voicings retain
their written notes and registers at every change.

The same gesture works without polyphony:

```text
bass = sequence {
  subdiv 2n
  notes 1_2 ; b7_2 ; 4_2 ; 5 7
  rhythm comp
}
```

Single notes, chords, and mixtures of both are the same kind of pitched value.
`notes` and `chords` are interchangeable readability choices, but only one may
appear in a sequence.

An inline rhythm is still useful when it shares the sequence subdivision:

```text
notes 1_2 4_2
rhythm x ~ x x
```

With a rhythm gesture, the pitched lane owns values and their held durations;
the rhythm owns attacks, rests, attack probability, ratchets, and attack-local
attributes. The `duration` lane consequently affects rhythm cells. Gate,
Velocity, Slide, Ratchet, and other event-indexed lanes advance on rhythm hits.
Ambiguous attack modifiers in the held pitch lane are rejected instead of
being silently overridden.

Standalone `x` is a rhythm hit only in a rhythm pattern. `x1` in a direct
`notes` or `chords` lane remains a ghosted pitched event.

### Slides with a gesture

A slide can express either a particular pitch transition or a particular
gesture position. Mark a target value to slide on the first successful hit in
its span:

```text
notes 1_2 >4_2 5_2 >b7_2
rhythm comp
```

Or put the slide in a reusable gesture:

```text
legato = rhythm {
  subdiv 8n
  events x >x x ~
}
```

`>x` connects that hit from the preceding sounding value. `>pitch` is consumed
by the first successful hit after entering its pitch span; later hits in the
same span are ordinary attacks. With no valid predecessor, either form becomes
a normal attack. Chord slides retain the existing rule that source and target
must have the same statically known voice count. An explicit slide time may be
placed on the target or the rhythm hit, but not both; the `slide` lane and then
the sequence `glide` value provide the usual fallbacks.

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
| `~_`, `~{len=2}` | rest for two base slots |
| `~!3` | repeat a rest three times |
| `_` | extend the preceding event by one base slot |
| `_!3` | extend the preceding event by three base slots |
| `1_`, `1__`, `1_3` | event lasting two or three base slots |
| `1.`, `1..` | dotted or double-dotted duration |
| `1!3` | repeat as three logical events |
| `1*3` | three ratchets within one event |
| `1?0.25` | sound with 25 percent probability; otherwise rest |
| `1??0.25` | exist with 25 percent probability; otherwise omit the span |
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

### Probability: silence and omission

The two probability suffixes answer different musical questions:

| Suffix | Question | Result on a miss |
| --- | --- | --- |
| `?P` | Should this event sound? | A rest occupying the complete event span |
| `??P` | Should this event exist in this pass? | No event and no elapsed span |

`P` is a probability from zero to one. For example, `1?0.35` sounds degree 1
on 35 percent of passes and rests on the others. `1??0.35` omits the event on a
miss, so the next event occurs earlier and the complete Notes pass is shorter.
A bare `?` or `??` means `0.5`.

Both decisions are made once per logical event, so all ratchets share the same
result. Replication creates separate logical events: the two copies in
`1??0.5!2` receive independent, seeded decisions. The seed and written event
identity make the results repeatable; omitting one event does not change the
random pitch chosen by a later written event.

A `?` miss becomes an ordinary rest. It closes Gate, consumes Duration, and
advances the numerical lanes which apply to that event. A `??` miss is resolved
before Duration or any event-based lane is read. Independent numerical lanes
therefore advance only for surviving events and retain their own cycle lengths.
Offset and CV lanes follow the resulting shorter score time. Every lane written
with `...` is remapped to the shortened Notes pass; ordinary lanes continue
their independent cycles across its boundary.

Presence probability is accepted on a top-level pitched event, a top-level
rest, or a complete top-level group. An omitted rest does not interrupt the
currently sounding predecessor, so a later slide can continue from it. At
least one top-level element must have guaranteed presence, preventing a pass
from collapsing to zero time. Presence probability inside a subdividing group
and its combination with a Euclidean suffix are rejected because the parent
span would otherwise have ambiguous subdivision semantics.

Repeated underscores and the numbered form express the same elongation, so
`1__` and `1_3` both occupy three base slots. Dots may follow either form:
`1_.` starts from two slots and dots that duration to three. In a Euclidean
suffix `(PULSES,STEPS,ROTATION)`, pulses may be zero but cannot exceed the
positive step count; rotation is an optional signed number of cells.

When several marks belong to one note, build it from left to right in this
order:

1. put `>` first for a slide;
2. add `x`, `^`, or `^^` for the dynamic gesture;
3. write the pitch or chord and any register mark;
4. add duration, Euclidean, and ratchet marks;
5. add probability and repetition; and
6. finish with `{attributes}`.

For example, `>^1'_3*2?0.8!4` creates four three-slot upper-tonic slides; each
has two ratchets and an 80 percent chance of sounding.

### Groups and choices

Square brackets subdivide one parent span:

```text
notes [1 2] [3 4 5] [1_ 5]
```

The first group plays two halves, the second three thirds, and `[1_ 5]` divides
its span in a 2:1 ratio. A suffix on a group applies to the complete group:

```text
notes [1 3 5]_2 [2 4](3,8) [6 7]?0.5 [1 8]!2
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
the preceding sounding pitch; after a rest or `?` miss it becomes a normal
attack. A `??` miss contributes no event and therefore does not break that
connection. Ties and slides do not cross a random/alternate branch,
arrangement boundary, or sequence wrap. Chord slides require the same known
voice count on both sides. If a presence omission exposes a predecessor with a
different voice count, the target becomes a normal attack.

`gate .5` holds Gate for half an event. `{stacc}` and `{ten}` provide concise
Gate presets of `.25` and `.95` without changing the event span. They combine
with ordinary notes, accents, quiet notes, slides, and ratchets. On a slide they
control the target's eventual release; the source still remains connected into
the slide. Ratchets apply the selected Gate fraction to every subspan.

The sequence `glide` value supplies the default slide time. An explicit zero in
the Slide lane or `{slide=0}` requests an immediate pitch transition rather
than inheriting that default; `.` in the Slide lane inherits `glide`.

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
| `gate=VALUE` | Gate fraction `0..1`, note value, or non-negative `ms` |
| `slide=VALUE` | slide time in beats, note values, or `ms`; valid only on `>` |
| `len=VALUE` | positive base-slot multiplier; the explicit form of a duration suffix |

`len` can also set a rest's span, as in `~{len=3/2}`. The other attributes
describe a pitched event and are not accepted on rests or ties.

Gate resolution follows musical specificity: the normal `.8` or named
`stacc`/`ten` preset is the fallback; an explicit Gate-lane value overrides it;
`.` in that lane retains it; and inline `gate=VALUE` wins last.

For expressive parts, prefer `stacc`, `ten`, `quiet`, accents, ghosts, and
inline overrides on the notes that need them. A dense Gate lane replaces the
named Gate articulation at every position containing a number. Use a Gate or
Velocity lane when its independent cycle is the musical idea, and use `.` where
the note's own articulation should remain in force.

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
Probability can be applied to either kind of random pitch: `[1|3|5]?0.5` chooses
one of those notes or silence, while `$?0.5` chooses a random scale note or
silence. Both forms retain their complete score-time slot. Use
`[1|3|5]??0.5` or `$??0.5` when a miss should remove the slot and shorten the
Notes pass.

Numerical and CV lanes use continuous uniform or normal values:

```text
velocity $u{.45,.9}
velocity $n{.7,.12}
gate     $u{10ms,25ms}
octave   $u{3,5}
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
octave   3 4 3
velocity .72 .55
duration 4n 4n 4nd
gate     .8 .4 16n
slide    . 80ms 32n
ratchet  1 1 2
offset   -10ms!2 +8ms
```

This block is a syntax inventory, not a recommended sequence template. Most
phrases need few or none of these lanes because notes already carry dynamics,
Gate articulation, duration, slides, ratchets, and octave marks.

| Lane | Accepted values | Normal default | Advances on |
| --- | --- | --- | --- |
| `octave` | absolute integer octave | octave declared by `tonic` | pitched events |
| `velocity`, `vel` | `0..1` | `.72` | pitched events |
| `duration`, `dur` | positive beats or note values | subdivision step | notes and rests |
| `gate` | `0..1`, note values, or non-negative `ms` | `.8` | pitched events |
| `slide` | non-negative beats, note values, or `ms` | `glide` | pitched events |
| `ratchet` | positive integer | `1` | pitched events |
| `offset` | signed beats, note values, or `ms` | on-grid | score time |

`.` consumes a lane position while inheriting the normal value. `.1` is the
number one tenth. `.!3` repeats the default across three positions.

Free numerical lanes advance independently according to the final column of
the table. A `?` miss is a rest and advances the applicable lanes. A `??` miss
does not exist as an event, so it advances none of them. A free lane may cycle
several times during one Notes pass or leave values unused. It keeps its cursor
when the Notes pattern repeats, wraps only at the end of its own values, and
therefore need not share the Notes cycle length. A `rate` pipeline likewise
keeps score-time phase across repeated passes of the same sequence.

A semicolon between pattern terms is a visual separator. It adds no event and
no duration, so these two lanes are identical:

```text
notes 1 2 3 4 5 6 7 8
notes 1 2 3 4 ; 5 6 7 8
```

This is especially useful with `subdiv 16n`, where groups of four terms make
the beats visible. Semicolons work in pitched, rhythm, numerical, and CV
patterns. They must occur between terms; they are neither barlines nor
statement terminators.

An accent prefix raises Velocity to its built-in accent floor and drives the
separate Accent output for the effective Gate. A Slide value is read for every
pitched event but matters only when the note uses `>`.

### Aligning a lane to events

A leading `...` right-aligns explicit values with one realized event-pattern
pass; a trailing `...` left-aligns them. That is the pitched pattern in direct
mode and one rhythm-gesture pass in separated mode. An ellipsis between two
groups applies values at both edges and leaves the middle positions at the
lane's normal default:

```text
notes 1 2 3 4 5 6 7 8
gate  ... .1!2          // short gates on notes 7 and 8

notes 1 2 3 4 5 6 7 8
gate  .1!2 ...          // short gates on notes 1 and 2

notes 1 2 3 4 5 6 7 8
gate  .5 .1 ... .1 .5   // shape notes 1, 2, 7, and 8 only
```

`...` may appear once, at the beginning, at the end, or between at least one
value on either side. `!N` repetitions count as positions; ratchets and note
elongation do not. Random and alternate Notes branches must have equal cell
counts when another lane is aligned to them. An aligned lane cannot also use
`rate`.

Presence omissions are resolved before alignment. A right-aligned lane remains
attached to the final surviving positions, and a left-aligned lane remains
attached to the first surviving positions. If omissions leave fewer positions
than explicit lane values, right alignment uses the final values and left
alignment uses the first values. Values on both sides of a middle ellipsis stay
attached to their respective edges. If omissions make the pass shorter than
the two edge groups together, each cell uses the value from the nearer outer
edge; a central tie uses the left group. Aligned lanes are recalculated from
the start of every realized event-pattern pass; free lanes continue
independently.

## CV1, CV2, and CV3

CV lanes use volts and follow score time, including rests, ties, and `?`
misses. A `??` miss contributes no score time. Ratchets share their parent's CV
sample.

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

### CV envelopes

A CV output can generate an AD, AR, or ADSR envelope from the Notes events.
The short form uses the envelope as the complete CV signal:

```text
notes 1 3 5 7
cv1 env ad
cv2 env ar 10ms 400ms depth 3
cv3 env adsr 5ms 200ms .4 600ms depth 2 curve .25
```

An envelope can also be added to a constant or sequenced CV lane. This is
useful for placing a transient sweep over a slower cutoff or timbre contour:

```text
cv1 0 .5 1 .5 |> interp smooth
  |> add env ad 5ms 300ms depth .8

cv2 1 |> add env ad . 16n depth 2
```

The parameter order follows the envelope name:

| Form | Positional parameters |
| --- | --- |
| `env ad` | attack, decay |
| `env ar` | attack, release |
| `env adsr` | attack, decay, sustain, release |

Trailing parameters may be omitted. A `.` keeps the default at an earlier
position. Note values such as `16n`, `8nd`, and `8nt` follow standard musical
durations. Times ending in `ms` or `s` use elapsed time. Bare numbers and
fractions use quarter-note beats, so `1/4` and `16n` have equal duration.
Sustain is a proportion from 0 to 1.

| Option | Meaning | Default |
| --- | --- | --- |
| `depth V` | signed peak voltage; a negative value inverts the envelope | `5` |
| `curve C` | segment shape from linear at `-1` through increasingly analogue-like curvature at `0` and `1` | `0` |
| `follow velocity` | scale the peak by the resolved note Velocity | disabled |
| `follow vel` | concise alias for `follow velocity` | disabled |
| `accent M` | additionally multiply accented attacks by `M` | `1` |

The time defaults are 5 ms attack, 300 ms AD decay, 300 ms AR release, and
5 ms / 250 ms / 0.5 / 300 ms for ADSR. `depth`, `curve`, `follow`, and
`accent` may follow the positional parameters in any order. One envelope is
available on each CV lane, and `add env` is the final operation in that lane's
pipeline. The result is left unclamped, allowing deliberate voltage ranges and
signed modulation.

AD begins on every attack, including ratchets, and completes independently of
the note gate. AR begins when the gate rises, holds its peak while the gate is
high, and releases when the gate falls. ADSR begins on an attack, moves to its
sustain level, and releases with the gate. A new AD or ADSR attack starts from
the current voltage, which keeps rapid retriggers continuous. Ties preserve the
current gate and slides introduce no new attack. Rests and sounded-note misses
release gate-controlled envelopes; a `??` omission creates no event.

Chords drive one envelope from their first logical voice. The peak captured at
each attack is

```text
depth * (resolved Velocity when follow is enabled, otherwise 1)
      * (accent multiplier on an accented attack, otherwise 1).
```

The resolved Velocity already includes the built-in effect of `^`, `^^`,
`quiet`, and the Velocity lane. `accent M` supplies an additional timbral lift
when desired. Pausing freezes the envelope at its current voltage. Stop and
reset clear it. If the clock disappears while transport remains active, a
held gate is released and the envelope completes against the most recently
measured tempo or elapsed time.

The inline CV sparkline displays the final output voltage. Envelope-only lanes
therefore show the envelope, while additive lanes show the base CV and envelope
combined. Editing the row preserves its existing trace history.

## Transforms

Pipelines transform a lane, a sequence, or an arrangement:

```text
notes 1 2 3 4 |> rotate 1 |> every 4 rev
offset -10ms!2 +8ms |> rate 1/2

variation = riff
  |> shift_degree 2
  |> sometimes .25 (rotate 1)
```

A pipeline written on a lane line transforms that lane. Further indented
pipeline lines inside the braces continue the preceding lane. Pipelines after
the closing brace transform the complete sequence, while pipelines after an
assignment transform every sequence referenced by that arrangement.

Transform scope is explicit:

| Target | Available transforms |
| --- | --- |
| Notes lane | `rev`, `rotate`, and pitch transforms |
| Rhythm lane or definition | `rev` and `rotate` |
| Numerical lane | `rev`, `rotate`, `rate` |
| CV lane | numerical-lane transforms plus `interp` and a final `add env` |
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
| `transpose_key KEY` | transpose from the sequence's written `key` to the nearest occurrence of `KEY` |
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
free-running numerical/CV lane without speeding up the pitched lane.

### Conditional transforms

```text
riff |> every 4 rev
riff |> sometimes .3 (rotate 1)
```

`every N` applies a transform on each positive Nth sequence pass.
`sometimes P` applies it with seeded probability `0..1`.
Conditional structural transforms applied to a complete sequence share one
decision across Notes and every auxiliary lane, preserving their relationship.
`scale` is an unconditional whole-sequence transform, while `interp` is an
unconditional CV-lane operation.

## Arrangements

An arrangement may combine independently defined sequences, transformed
versions of existing material, or both. Each named sequence keeps its own
settings, lanes, and length; arrangement playback moves to the next term after
the current sequence completes its pitched pass. In gesture mode that means
the complete held-value timeline, not one repetition of the gesture. See
[Song with distinct sections](#song-with-distinct-sections) for a complete
multi-section example.

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

The fixed `.5` Velocity and Gate lanes set the base intensity and half-step
Gate of this particular 303 phrase. They are sound settings shared by every
note; the accent marks still raise the selected steps above that base.

Parentheses group arrangement expressions:

```text
song = (verse * 2 + fill) * 2 + ending
```

A section boundary occurs only after the final direct event or held pitch span
has completed. The new section owns the next attack; ties and slides do not
cross between different sections.

## Musical examples

### Expressive bass line

```text
bass = sequence {
  subdiv 16n
  tonic E@2
  scale minor
  glide 32n
  notes ^1 1!2 x1 [5 6] >b7{stacc} ~ 1{ten} 3?0.35
}

fill = bass |> rotate 2 |> every 2 rev
song = bass * 3 + fill

seed 73
play song
```

The Notes pattern contains the performance intent: accent, repetition, ghost,
subdivision, slide, staccato, rest, and tenuto. Its final slot plays degree 3
with 35 percent probability and is otherwise silent. That slot always occupies
one sixteenth-note span, so probability changes the performance without
changing the bass line's length. The fill reuses the articulations while
changing structure. No auxiliary lane is needed.

### Song with distinct sections

```text
verse = sequence {
  subdiv 8n
  tonic D@3
  scale minor
  notes 1 3 4 ^5{stacc} 1' 7 5 4
}

chorus = sequence {
  subdiv 4n
  tonic Bb@3
  scale major
  notes I_2{ten} V_2{ten} vi_2{ten} IV_2{ten}
}

fill = sequence {
  subdiv 16n
  tonic D@3
  scale harmonic_minor
  notes [5 6 7] ^8*2 ~ V7
}

song = verse * 2 + chorus + verse + fill + chorus * 2
play song
```

`verse`, `chorus`, and `fill` are independent sequences. Each owns its
subdivision, tonal settings, pitched pass, and optional lanes. The arrangement
advances after each sequence completes its own pitched pass, so sections can have
different lengths without padding them to a shared grid.

### Seeded generative melody

```text
melody = sequence {
  subdiv 16n
  tonic A@3
  scale minor_pentatonic
  notes $u{1,10}(5,8) $n{5,1.25}(3,8,2)
}
|> every 4 (rotate 1)

seed 2026
play melody
```

The two Euclidean figures fill sixteen cells. Every generated pitch belongs to
A minor pentatonic. The seed makes the result repeatable, while every fourth
pass rotates the prepared structure.

### One intentional CV curve

```text
texture = sequence {
  subdiv 8n
  tonic D@3
  scale dorian
  notes 1{ten} [3 5] 7{quiet} <8 6>
  cv1 0 . 5 . |> interp smooth
}

play texture
```

Here CV1 has a separate job, such as opening a filter across the phrase. The
`.` entries leave room for the smooth curve to continue. Note articulation
still controls the polyphonic Gate and Velocity outputs.

### Humanized timing

```text
groove = sequence {
  notes [1 5] 3 [4 6 8] 5
  offset -7ms 4ms 0 |> rate 1/2
}
|> swing .58 32n
|> early random 3ms

seed 99
play groove
```

The Offset lane is present because its slower three-value cycle is the point of
the example. Swing supplies the main groove while the seeded early displacement
adds a repeatable small variation.

## Current limits

The following are not available yet:

- multi-event or recursively nested random/alternate branches;
- grouped subdivisions and attack-local dynamics inside a held `notes` or
  `chords` timeline; put those attacks in its rhythm gesture instead;
- `??` presence probability inside a subdividing group or beside a Euclidean
  suffix;
- continuous interpolation on `...`-aligned CV lanes;
- CV outputs beyond CV1, CV2, and CV3;
- MIDI input bindings and performance interpreters for bass, piano, and
  arpeggiation;
- scale-stacked chord generation; and
- expressions that derive CV signals from MIDI, other CVs, or functions.

Until those features arrive, use atomic alternatives, free-running CV for
continuous curves, explicit chord symbols/voicings, and the three physical CV
outputs.

## Formal grammar reference

This appendix is the exact syntax map for implementers and for resolving edge
cases. The preceding chapters introduce the same language through musical
examples. In the productions below, square brackets mean an optional part,
braces mean repetition, and `A | B` means a choice. Quoted words and
punctuation are literal source text.

### Top-level statements and names

```text
program     ::= { sequence | rhythm-definition | assignment | play | seed |
                  comment }
sequence    ::= NAME "=" "sequence" "{" sequence-line+ "}" pipeline*
rhythm-definition ::= NAME "=" "rhythm" "{" rhythm-definition-line+ "}"
rhythm-definition-line ::= "subdiv" positive-note-value |
                           "events" rhythm-pattern pipeline*
assignment  ::= NAME "=" expression pipeline*
play        ::= "play" NAME
seed        ::= "seed" UNSIGNED_INTEGER

expression  ::= term { "+" term }
term        ::= primary [ "*" POSITIVE_INTEGER ]
primary     ::= NAME | "(" expression ")"

NAME        ::= (LETTER | "_") { LETTER | DIGIT | "_" }
comment     ::= "//" text-to-end-of-line
```

A sequence, rhythm, or assignment name may be defined once. `+` concatenates
arrangement parts, `*` repeats its immediately preceding name or parenthesized
expression, and parentheses group an expression. Repetition therefore binds
more tightly than concatenation. Pipelines apply after the complete expression.
An assignment names an arrangement expression; `play` accepts a sequence or an
arrangement, not a rhythm definition. The selected arrangement loops until
another executed `play` command replaces it. The final `play` command selects
the arrangement, and the final `seed` command selects the random seed.

Names beginning with `cv` followed only by a non-zero decimal integer are
reserved for CV lanes and future typed signal values.
An identifier which is also a complete inline rhythm pattern is reserved for
that pattern and cannot name a reusable rhythm.

Blank lines and comments may appear between statements and sequence lines. A
comment may also follow any complete line. Statements are line-oriented; the
closing brace and a pipeline continuation may begin subsequent lines.

### Sequence lines and pipeline attachment

```text
sequence-line ::= setting-line | notes-line | chords-line | rhythm-line |
                  scalar-line | cv-line
setting-line  ::= "subdiv" positive-note-value |
                  "tonic" named-pitch ["@" SIGNED_INTEGER] |
                  "scale" SCALE_NAME |
                  "key" named-pitch |
                  "voicing" VOICING_RECIPE |
                  "glide" nonnegative-time
notes-line    ::= "notes" note-pattern pipeline*
chords-line   ::= "chords" note-pattern pipeline*
rhythm-line   ::= "rhythm" (NAME | rhythm-pattern) pipeline*
scalar-line   ::= scalar-lane scalar-pattern pipeline*
cv-line       ::= cv-name envelope-source |
                  cv-name scalar-pattern cv-pipeline*
cv-name       ::= "cv1" | "cv2" | "cv3"
cv-pipeline   ::= pipeline | "|>" "add" envelope-source

envelope-source ::= "env" envelope-shape envelope-option*
envelope-shape  ::= "ad" [env-time-or-default [env-time-or-default]] |
                    "ar" [env-time-or-default [env-time-or-default]] |
                    "adsr" [env-time-or-default
                      [env-time-or-default
                        [sustain-or-default [env-time-or-default]]]]
env-time-or-default ::= envelope-time | "."
sustain-or-default  ::= PROBABILITY | "."
envelope-time       ::= positive-note-value |
                        NONNEGATIVE_NUMBER ["ms" | "s"]
envelope-option     ::= "depth" NUMBER |
                        "curve" NUMBER |
                        "follow" ("velocity" | "vel") |
                        "accent" NONNEGATIVE_NUMBER

scalar-lane   ::= "octave" | "velocity" | "vel" |
                  "duration" | "dur" | "gate" | "slide" |
                  "ratchet" | "offset"
```

`VOICING_RECIPE` is `basic`, `rootless_3notes`, or `rootless_4notes`. `vel` and
`dur` are scalar aliases. Settings accept one value and cannot have an
inline pipeline. `SCALE_NAME` is one of the names listed under
[Scale degrees](#scale-degrees). A sequence requires exactly one `notes` or
`chords` lane. A `rhythm` lane is optional and may contain an inline pattern or
name a reusable rhythm definition.
Each setting and canonical lane may appear at most once, so `velocity` and its
`vel` alias cannot both occur in the same sequence. A pipeline on the same line
as a lane transforms that lane. A pipeline on a following line inside the
braces continues the preceding lane:

```text
notes 1 2 3 4
  |> rotate 1
  |> every 4 rev
```

A pipeline following the closing brace transforms the complete sequence:

```text
riff = sequence {
  notes 1 2 3 4
}
|> swing .58 32n
|> sometimes .25 (rotate 1)
```

### Note-pattern syntax

```text
note-pattern ::= note-element { pattern-separator note-element }
note-element ::= group | note-event | rest | tie

note-event   ::= [">"] ["x" | "^" | "^^"] pitched-value
                 [duration] [euclidean] [ratchet]
                 [event-probability] [replication] [attributes]
rest         ::= "~" [duration] [presence-probability] [replication]
                 ["{" "len=" value "}"]
tie          ::= "_" [replication]

rhythm-pattern ::= rhythm-event { pattern-separator rhythm-event }
rhythm-event   ::= rhythm-hit | rest | tie
rhythm-hit     ::= [">"] "x" [duration] [euclidean] [ratchet]
                   [event-probability] [replication] [attributes]

group        ::= group-primary [duration] [euclidean]
                 [event-probability] [replication]
group-primary ::= "[" note-pattern "]" |
                  "[" note-pattern ("|" note-pattern)+ "]" |
                  "<" note-element+ ">"

pitched-value ::= random-pitch | pitch | chord ["/" pitch]
pitch         ::= (named-pitch | scale-degree) [register]
chord         ::= (explicit-voicing | roman-chord |
                   jazz-chord [factor-voicing]) [register]
named-pitch   ::= ("A" | "B" | "C" | "D" | "E" | "F" | "G")
                  ["b" | "#"]
scale-degree  ::= {"b" | "#"} SIGNED_INTEGER
explicit-voicing ::= "(" pitch SPACE pitch { SPACE pitch } ")"
factor-voicing ::= ":" "(" chord-factor { SPACE chord-factor } ")"
chord-factor  ::= {"b" | "#"} chord-degree
roman-chord   ::= {"b" | "#"} roman-root [chord-suffix]
roman-root    ::= "I" | "II" | "III" | "IV" | "V" | "VI" | "VII" |
                  "i" | "ii" | "iii" | "iv" | "v" | "vi" | "vii"
jazz-chord    ::= named-pitch chord-suffix
chord-suffix  ::= triad-quality [extension] {modification} ["alt"] |
                  extension [post-extension-quality] {modification} ["alt"] |
                  modification {modification} ["alt"] |
                  "alt"
triad-quality ::= "maj" | "min" | "m" | "dim" | "aug" |
                  "sus2" | "sus4"
post-extension-quality ::= "sus" | "sus2" | "sus4"
extension     ::= "5" | "6" | "7" | "9" | "11" | "13"
modification  ::= ("b" | "#" | "add") chord-degree
chord-degree  ::= "1" | "2" | "3" | "4" | "5" | "6" | "7" |
                  "9" | "11" | "13"
random-pitch  ::= "$" | "${" SIGNED_INTEGER "," SIGNED_INTEGER "}" |
                  "$u{" SIGNED_INTEGER "," SIGNED_INTEGER "}" |
                  "$n{" NUMBER "," POSITIVE_NUMBER "}" |
                  "$c{" SIGNED_INTEGER "," SIGNED_INTEGER "}" |
                  "$cn{" NUMBER "," POSITIVE_NUMBER "}"
register      ::= "@" SIGNED_INTEGER [relative-register] |
                  relative-register
relative-register ::= "'" {"'"} | "," {","}
```

The order shown above is the accepted suffix order. Atomic note events support
ratchets and attributes; groups support their own duration, Euclidean,
probability, and replication suffixes. Rests accept duration, presence
probability, replication, and `len`; ties accept replication. A standalone `_`
is a tie, while `_` attached to an event is a duration suffix.

Angle brackets choose one plain pitched alternative on each sequence pass.
Brackets containing `|` make a seeded random choice. In the current version,
each alternate or random branch must contain one plain pitch or chord; its
duration and probability suffixes belong to the group as a whole.

### Suffixes, attributes, and scalar values

```text
duration     ::= (underscore-run | "_" POSITIVE_INTEGER) [dots] | dots
underscore-run ::= "_" {"_"}
dots         ::= "." | ".."
euclidean    ::= "(" UNSIGNED_INTEGER "," POSITIVE_INTEGER
                 ["," SIGNED_INTEGER] ")"
ratchet      ::= "*" POSITIVE_INTEGER
event-probability ::= sound-probability | presence-probability
sound-probability ::= "?" [NUMBER]
presence-probability ::= "??" [NUMBER]
replication  ::= "!" POSITIVE_INTEGER
attributes   ::= "{" attribute {"," attribute} "}"
attribute    ::= "quiet" | "stacc" | "ten" |
                 "vel=" value | "gate=" value |
                 "slide=" value | "len=" value

scalar-pattern ::= scalar-term { pattern-separator scalar-term } |
                   "..." SPACE scalar-term
                     { pattern-separator scalar-term } |
                   scalar-term { pattern-separator scalar-term } SPACE "..." |
                   scalar-term { pattern-separator scalar-term } SPACE "..."
                     SPACE scalar-term { pattern-separator scalar-term }
scalar-term    ::= (value | "." | random-scalar) [replication]
pattern-separator ::= SPACE | [SPACE] ";" [SPACE]
random-scalar  ::= "$u{" value "," value "}" |
                   "$n{" value "," value "}"

value               ::= NUMBER ["ms"] | note-value
note-value          ::= ["+" | "-"] NOTE_DENOMINATOR ("nd" | "nt" | "n")
positive-note-value ::= ["+"] NOTE_DENOMINATOR ("nd" | "nt" | "n")
nonnegative-time    ::= positive-note-value | NONNEGATIVE_NUMBER
NUMBER              ::= SIGNED_INTEGER "/" POSITIVE_INTEGER | DECIMAL
DECIMAL             ::= ["+" | "-"]
                        (DIGITS ["." DIGITS] | "." DIGITS)
SIGNED_INTEGER      ::= ["+" | "-"] DIGITS
UNSIGNED_INTEGER    ::= DIGITS
POSITIVE_INTEGER    ::= NONZERO_DIGIT {DIGIT}
POSITIVE_NUMBER     ::= NUMBER whose value is greater than zero
NONNEGATIVE_NUMBER  ::= NUMBER whose value is at least zero
NOTE_DENOMINATOR    ::= a power of two from 1 to 128
PROBABILITY         ::= NUMBER whose value is from zero to one
DIGITS              ::= DIGIT {DIGIT}
DIGIT               ::= "0" | NONZERO_DIGIT
NONZERO_DIGIT       ::= "1" | "2" | "3" | "4" | "5" |
                        "6" | "7" | "8" | "9"
LETTER              ::= "A".."Z" | "a".."z"
SPACE               ::= one or more spaces or tabs
```

Euclidean pulses may be zero and must not exceed the positive step count.
Presence probability is restricted to top-level events, rests, and complete
top-level groups, and cannot accompany a Euclidean suffix. A direct pitched pass must
contain at least one element whose presence is guaranteed. Bare `?` and `??`
suffixes mean probability `0.5`. `NUMBER` may be an integer, decimal, or
fraction where the surrounding feature permits it. A note value uses a
power-of-two denominator from 1 through 128 and the suffix `n`, `nd`, or `nt`.
`subdiv` requires a positive note value. Millisecond units are accepted by
Gate, Slide, Offset, timing transforms, and the corresponding random forms. CV
values are volts and omit a unit. Settings and duration spans use score beats,
note values, or dimensionless values as described above.

An envelope time without a unit is measured in quarter-note beats.
Envelope curve values range from -1 to 1, and ADSR sustain values range from
zero to one. An envelope-only CV line accepts no pipeline. On a scalar CV line,
`add env` is unconditional, appears once, and closes the pipeline.

Scalar patterns contain flat values rather than note groups. `.` consumes a
position and inherits that lane's normal value. One ellipsis may occur at an
edge, or between two sets of values, to align the supplied values with the
beginning, end, or both edges of the active event structure. A semicolon is accepted
where whitespace would separate two pattern terms and has no semantic effect.

### Pipelines

```text
pipeline  ::= "|>" [condition] transform-call
condition ::= "every" POSITIVE_INTEGER |
              "sometimes" PROBABILITY
transform-call ::= OPERATION { argument } |
                   "(" OPERATION { argument } ")"
```

Parentheses delimit a transform call used after a condition; they do not make a
new sequence. Pipelines compose from left to right. The Transforms section lists
every accepted operation, its arguments, and the targets to which it applies.
