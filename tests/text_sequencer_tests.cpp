#include "tfseq.hpp"
#include "tfseq_parser.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

int failures = 0;

void checkImpl(bool condition, const std::string &message) {
  if (condition)
    return;
  std::cerr << "FAIL: " << message << '\n';
  ++failures;
}

#define check(condition, message)                                              \
  checkImpl((condition), std::string(__func__) + ": " + (message))

bool close(float left, float right, float tolerance = 1e-5f) {
  return std::fabs(left - right) <= tolerance;
}

void pegFrontendBuildsTypedSyntax() {
  const std::string source = R"(acid = sequence {
  notes 1 <2 b3> [4|5] 1'!2 |> rotate 1 |> every 4 rev
  articulation x [x x] x(3,8,1)
}
|> sometimes .25 rev
song = acid * 8 + acid
play song
)";
  const auto parsed = tfseq::syntax::Parse(source);
  check(static_cast<bool>(parsed), parsed.diagnostic.message);
  if (!parsed)
    return;
  check(parsed.document.statements.size() == 3,
        "PEG front end emits three typed statements");
  const auto *sequence = std::get_if<tfseq::syntax::SequenceDefinition>(
      &parsed.document.statements[0]);
  check(sequence && sequence->lanes.size() == 2,
        "PEG front end retains sequence lanes");
  if (sequence) {
    const auto &pattern = sequence->lanes[0].pattern;
    check(
        pattern.steps.size() == 4 &&
            pattern.steps[1].kind == tfseq::syntax::PatternKind::CycleChoice &&
            pattern.steps[2].kind == tfseq::syntax::PatternKind::RandomChoice &&
            pattern.steps[3].kind == tfseq::syntax::PatternKind::Repeat,
        "PEG front end produces typed choice and repetition nodes");
    check(sequence->lanes[0].pipelines.size() == 2,
          "PEG pattern grammar separates inline pipelines");
    check(sequence->lanes[1].pattern.steps[2].kind ==
                  tfseq::syntax::PatternKind::Euclidean &&
              sequence->lanes[1].pattern.steps[2].arguments.size() == 3,
          "Euclidean structure is typed before articulation semantics");
    check(sequence->pipelines.size() == 1,
          "PEG document grammar retains sequence pipelines");
  }
  const auto *assignment =
      std::get_if<tfseq::syntax::Assignment>(&parsed.document.statements[1]);
  check(assignment && assignment->expression.terms.size() == 2 &&
            assignment->expression.terms[0].repeats == 8,
        "PEG assignment grammar parses concatenation and repetition");
  if (sequence) {
    check(sequence->name.span.line == 1 && sequence->name.span.column == 1,
          "PEG AST locations are retained for editor cursors");
    check(sequence->lanes[0].pattern.steps[1].span.line == 2,
          "grouped pattern tokens retain their source line");
  }
}

void pegFrontendRejectsMalformedStructure() {
  const auto unclosed =
      tfseq::syntax::Parse("riff = sequence {\n  notes 1 [2 3\n}\nplay riff\n");
  check(!unclosed, "PEG front end rejects an unclosed pattern group");
  check(unclosed.diagnostic.line == 2,
        "PEG syntax errors identify the malformed source line");

  const auto unexpected =
      tfseq::syntax::Parse("riff = sequence {\n  notes 1 2)\n}\n");
  check(!unexpected, "PEG front end rejects an unexpected closing bracket");
  check(unexpected.diagnostic.line == 2,
        "unexpected bracket diagnostics retain their source line");

  const auto unclosedVoicing =
      tfseq::syntax::Parse("riff = sequence {\n  notes (1 b3 5\n}\n");
  check(!unclosedVoicing,
        "PEG front end owns and rejects an unclosed voicing delimiter");

  const auto comments = tfseq::Compile(
      "// header\r\nriff = sequence {\r\n  notes 1 2 // notes\r\n}\r\n"
      "play riff // transport\r\n");
  check(static_cast<bool>(comments),
        "PEG document grammar accepts comments and CRLF input");
}

void typedPatternTreeOwnsReusableStructure() {
  const auto nested = tfseq::syntax::Parse(R"(a = sequence {
  notes <1 [2|3]>!2
  articulation [x x!2]!2
}
play a
)");
  check(static_cast<bool>(nested), nested.diagnostic.message);
  if (!nested)
    return;
  const auto *sequence = std::get_if<tfseq::syntax::SequenceDefinition>(
      &nested.document.statements.front());
  check(sequence &&
            sequence->lanes[0].pattern.steps[0].kind ==
                tfseq::syntax::PatternKind::Repeat &&
            sequence->lanes[0].pattern.steps[0].children[0].kind ==
                tfseq::syntax::PatternKind::CycleChoice &&
            sequence->lanes[0].pattern.steps[0].children[0].children[1].kind ==
                tfseq::syntax::PatternKind::RandomChoice,
        "nested grouping has one domain-neutral syntax tree");

  const auto articulation = tfseq::Compile(R"(a = sequence {
  notes 1
  articulation [x x!2]!2
}
play a
)");
  check(articulation &&
            articulation.program->semantic().sequences[0].articulation.size() ==
                2 &&
            articulation.program->semantic()
                    .sequences[0]
                    .articulation[0]
                    .atoms.size() == 3,
        "generic repeat nodes compose at group and child levels");

  const auto before = tfseq::Compile("a = sequence {\n notes 1\n}\nplay a\n");
  const auto after = tfseq::Compile("a = sequence {\n notes 2\n}\nplay a\n");
  check(before && after &&
            before.program->semantic().sequences[0].stableId ==
                after.program->semantic().sequences[0].stableId,
        "definition identity remains stable across content edits");
}

void lineCommentsCanTruncatePatterns() {
  const auto compiled = tfseq::Compile(R"(riff = sequence {
  // Keep alternatives close by while auditioning a shorter version.
  notes 1 2 // 8 7 5 4
  articulation x // x ~ x
}

play riff // comments are also valid on commands
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;

  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const auto first = runtime.next(0.0);
  const auto firstCount = first.count;
  const auto firstPitch = first.count ? first.events[0].pitchVolts : 0.f;
  const auto second = runtime.next(1.0);
  const auto secondCount = second.count;
  const auto secondPitch = second.count ? second.events[0].pitchVolts : 0.f;
  const auto wrapped = runtime.next(2.0);
  check(firstCount == 1 && close(firstPitch, 0.f),
        "a truncated pattern keeps notes before //");
  check(secondCount == 1 && close(secondPitch, 2.f / 12.f),
        "the final uncommented note remains active");
  check(wrapped.count == 1 && close(wrapped.events[0].pitchVolts, 0.f),
        "notes after // are excluded from the compiled pattern");
}

void selectionEvaluationReplacesOnlyContainingStatements() {
  const std::string evaluated = R"(a = sequence {
  notes 1 2
}
b = sequence {
  notes 3 4
}
play a
)";
  const std::string draft = R"(a = sequence {
  notes 5 6
}
b = sequence {
  notes 7 8
}
play a
)";
  const auto evaluatedDocument = tfseq::syntax::Parse(evaluated);
  const auto draftDocument = tfseq::syntax::Parse(draft);
  check(evaluatedDocument && draftDocument,
        "selection fixtures parse before typed merging");
  if (!evaluatedDocument || !draftDocument)
    return;
  const auto selectionBegin = static_cast<int>(draft.find("notes 5 6"));
  const auto contextual = tfseq::syntax::MergeSelectionDocuments(
      evaluatedDocument.document, evaluated, draftDocument.document, draft,
      selectionBegin, selectionBegin + 9);
  check(static_cast<bool>(contextual),
        "selection merge: " + contextual.diagnostic.message);
  if (!contextual)
    return;
  const auto compiled = tfseq::Compile(contextual.document);
  check(static_cast<bool>(compiled),
        "selection compile: " + compiled.diagnostic.message);
  if (!compiled)
    return;
  const tfseq::Sequence *a = nullptr;
  const tfseq::Sequence *b = nullptr;
  for (const auto &sequence : compiled.program->semantic().sequences) {
    if (sequence.name == "a")
      a = &sequence;
    else if (sequence.name == "b")
      b = &sequence;
  }
  check(a && b, "selection context retains unselected definitions");
  if (a && b) {
    check(a->notes[0].values[0].voices[0].degree == 5,
          "the containing selected sequence is updated");
    check(b->notes[0].values[0].voices[0].degree == 3,
          "an unrelated draft edit remains unevaluated");
    check(!b->notes[0].span.valid(),
          "inactive draft text cannot receive a stale playback cursor");
    check(a->notes[0].span.begin == selectionBegin + 6,
          "selected statement source spans remain aligned to the editor");
  }

  const auto stopBegin = static_cast<int>(draft.find("play a"));
  std::string stoppedDraft = draft;
  stoppedDraft.replace(static_cast<std::size_t>(stopBegin), 6, "stop  ");
  const auto stoppedDocument = tfseq::syntax::Parse(stoppedDraft);
  const auto stoppedContext =
      stoppedDocument
          ? tfseq::syntax::MergeSelectionDocuments(
                evaluatedDocument.document, evaluated, stoppedDocument.document,
                stoppedDraft, stopBegin, stopBegin + 4)
          : tfseq::syntax::SelectionDocumentResult{};
  const auto stopped = stoppedContext ? tfseq::Compile(stoppedContext.document)
                                      : tfseq::Compile("");
  check(stoppedContext && stopped && stopped.program->semantic().stopped,
        "a selected transport statement replaces the active transport only");
}

void explicitSingleRepeatRemainsAnArrangement() {
  const auto compiled = tfseq::Compile(R"(a = sequence {
  notes 1
}
song = a * 1
play song
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (compiled) {
    check(compiled.program->semantic().arrangement.size() == 1 &&
              compiled.program->semantic().arrangement[0].cycles == 1,
          "an explicit * 1 remains a finite arrangement section");
  }

  const auto largeRepeat = tfseq::Compile(R"(a = sequence {
  notes 1
}
song = a * 257
play song
)");
  check(largeRepeat &&
            largeRepeat.program->semantic().arrangement.size() == 1 &&
            largeRepeat.program->semantic().arrangement[0].cycles == 257,
        "arrangement repeats are not constrained by a small musical ceiling");
}

void compileAndCycleIndependentLanes() {
  const std::string source = R"(riff = sequence {
  cycle 8
  tonic C@4
  scale minor
  notes 1 2 b3 4 5 6 b7 8 |> rotate 1
  articulation x x x x x x x x
  velocity .61
  accent + . .
  duration 1
}
play riff
)";
  auto compiled = tfseq::Compile(source);
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;

  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const float expectedPitch[] = {2.f / 12.f, 2.f / 12.f, 5.f / 12.f,
                                 7.f / 12.f};
  const float expectedVelocity[] = {.88f, .61f, .61f, .88f};
  for (int beat = 0; beat < 4; ++beat) {
    const auto events = runtime.next(beat);
    check(events.count == 1, "one event should be emitted per plain step");
    if (events.count != 1)
      continue;
    check(close(events.events[0].pitchVolts, expectedPitch[beat]),
          "note lane rotates independently");
    check(close(events.events[0].velocity, expectedVelocity[beat]),
          "three-item accent lane displaces against notes");
    check(events.events[0]
              .cursors[static_cast<std::size_t>(tfseq::CursorLane::Notes)]
              .valid(),
          "emitted note carries a source cursor");
  }
}

void parseFirstClassArticulation() {
  const std::string source = R"(voice = sequence {
  notes 1 2 3 4 5 6 7 8
  articulation x _ ~ > [x x] x*3 x!2 x(3,8,1)
  duration 1
}
play voice
)";
  auto compiled = tfseq::Compile(source);
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  const auto &art = compiled.program->semantic().sequences[0].articulation;
  check(art.size() == 16,
        "articulation expands replication and Euclidean steps");
  check(art[0].atoms[0].kind == tfseq::ArticulationKind::Attack,
        "x is an attack");
  check(art[1].atoms[0].kind == tfseq::ArticulationKind::Tie,
        "underscore is a tie");
  check(art[2].atoms[0].kind == tfseq::ArticulationKind::Rest,
        "tilde is a rest");
  check(art[3].atoms[0].kind == tfseq::ArticulationKind::Slide,
        "greater-than is a slide");
  check(art[4].atoms.size() == 2, "brackets subdivide one slot");
  check(art[5].atoms[0].ratchets == 3, "asterisk ratchets one slot");

  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  check(runtime.next(0).events[0].kind == tfseq::EventKind::Attack,
        "attack reaches runtime");
  check(runtime.next(1).events[0].kind == tfseq::EventKind::Tie,
        "tie reaches runtime without consuming a note");
  check(runtime.next(2).events[0].kind == tfseq::EventKind::Rest,
        "rest reaches runtime without consuming a note");
  check(runtime.next(3).events[0].kind == tfseq::EventKind::Slide,
        "slide reaches runtime");
  const auto subdivision = runtime.next(4);
  check(subdivision.count == 2 &&
            close(static_cast<float>(subdivision.events[1].beat), 4.5f),
        "bracketed attacks are evenly subdivided");
  const auto ratchet = runtime.next(5);
  check(ratchet.count == 3, "ratchet emits three attacks in one slot");
}

void articulationModifiersComposeInEitherOrder() {
  const auto compiled = tfseq::Compile(R"(a = sequence {
  notes 1
  articulation x*3?1 x?1*2
}
play a
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const auto first = runtime.next(0.0);
  check(first.count == 3,
        "ratchet and probability suffixes compose in either order");
  const auto second = runtime.next(1.0);
  check(second.count == 2,
        "articulation suffix parsing is not order-sensitive");
}

void transformByCycleAndArrange() {
  const std::string source = R"(a = sequence {
  cycle 2
  notes 1 2 |> every 2 rev
}
b = sequence {
  cycle 1
  notes 5
}
song = a * 2 + b
play song
)";
  auto compiled = tfseq::Compile(source);
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  check(compiled.program->semantic().arrangement.size() == 2,
        "arrangement stitches named sections");
  check(compiled.program->semantic().arrangement[0].cycles == 2,
        "arrangement repetition is retained");

  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  check(close(runtime.next(0).events[0].pitchVolts, 0.f),
        "first cycle starts in source order");
  check(close(runtime.next(1).events[0].pitchVolts, 2.f / 12.f),
        "first cycle continues in source order");
  check(close(runtime.next(2).events[0].pitchVolts, 2.f / 12.f),
        "every 2 reverses the second cycle");
  check(close(runtime.next(4).events[0].pitchVolts, 7.f / 12.f),
        "arrangement advances to the stitched section");
}

void rejectInvalidInput() {
  auto compiled =
      tfseq::Compile("bad = sequence {\n notes 1 nope\n}\nplay bad\n");
  check(!compiled, "invalid degree is rejected");
  check(compiled.diagnostic.line == 2 && compiled.diagnostic.column > 1,
        "diagnostic includes a useful source location");
}

void degreesContinueAcrossOctaves() {
  auto compiled = tfseq::Compile(R"(range = sequence {
  notes 0 -1 1 8
}
play range
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  check(close(runtime.next(0).events[0].pitchVolts, -1.f / 12.f),
        "degree zero is the scale step below C@4");
  check(close(runtime.next(1).events[0].pitchVolts, -3.f / 12.f),
        "negative degrees continue down the scale");
  check(close(runtime.next(2).events[0].pitchVolts, 0.f),
        "degree one is the tonic");
  check(close(runtime.next(3).events[0].pitchVolts, 1.f),
        "degree eight is the next tonic");
}

void scaleCardinalityAndExplicitOctaves() {
  auto compiled = tfseq::Compile(R"(p = sequence {
  cycle 3
  tonic C@4
  scale major_pentatonic
  notes 1 6 1'
}

o = sequence {
  cycle 3
  tonic C@4
  scale octatonic_half_whole
  notes 1 9 1'
}
song = p + o
play song
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  check(close(runtime.next(0).events[0].pitchVolts, 0.f),
        "pentatonic degree one is the tonic");
  check(close(runtime.next(1).events[0].pitchVolts, 1.f),
        "pentatonic degree six is the next tonic");
  check(close(runtime.next(2).events[0].pitchVolts, 1.f),
        "apostrophe raises a degree by one octave in any scale");
  check(close(runtime.next(3).events[0].pitchVolts, 0.f),
        "octatonic degree one is the tonic");
  check(close(runtime.next(4).events[0].pitchVolts, 1.f),
        "octatonic degree nine is the next tonic");
  check(close(runtime.next(5).events[0].pitchVolts, 1.f),
        "explicit octave remains scale-cardinality independent");
}

void absoluteAndRelativeRegistersRemainDistinct() {
  const auto compiled = tfseq::Compile(R"(registers = sequence {
  tonic C@4
  scale major
  notes 1, 1'' 1@-1 1@+1 Cm7' (C E G)' (C@3 E@4 G@4),
}
play registers
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;

  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  check(close(runtime.next(0).events[0].pitchVolts, -1.f),
        "comma lowers a scale degree by one relative octave");
  check(close(runtime.next(1).events[0].pitchVolts, 2.f),
        "relative octave marks compose after a scale degree");
  check(close(runtime.next(2).events[0].pitchVolts, -5.f),
        "a negative @ register is an absolute octave");
  check(close(runtime.next(3).events[0].pitchVolts, -3.f),
        "an explicitly signed positive @ register is also absolute");

  const auto jazz = runtime.next(4);
  check(jazz.count == 4 && close(jazz.events[0].pitchVolts, 1.f),
        "relative marks transpose an interpreted jazz chord");
  const auto shared = runtime.next(5);
  check(shared.count == 3 && close(shared.events[0].pitchVolts, 1.f) &&
            close(shared.events[1].pitchVolts, 4.f / 12.f + 1.f),
        "relative marks transpose a parenthesized shared-register voicing");
  const auto spread = runtime.next(6);
  check(spread.count == 3 && close(spread.events[0].pitchVolts, -2.f) &&
            close(spread.events[1].pitchVolts, -1.f + 4.f / 12.f),
        "relative marks compose with absolute registers in a spread voicing");
}

void octaveSuffixesAndChordsAreUnambiguous() {
  const auto compiled = tfseq::Compile(R"(harmony = sequence {
  tonic C
  scale major
  octave 3 4 3
  notes D@4 D7@3 Cm7@3 Bbm7b5@3 / D@2 (1 b3 5)@4 1'
}
play harmony
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  const auto &semanticChord =
      compiled.program->semantic().sequences[0].notes[1].values[0];
  check(semanticChord.meaning == tfseq::ChordValue::Meaning::JazzSymbol &&
            semanticChord.jazzSymbol == "D7@3" &&
            semanticChord.rootPitchClass == 2,
        "jazz harmonic intent survives compilation for future interpreters");
  const auto &explicitChord =
      compiled.program->semantic().sequences[0].notes[4].values[0];
  check(explicitChord.meaning == tfseq::ChordValue::Meaning::ExplicitVoicing,
        "parenthesized voicings remain explicitly voiced");
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());

  const auto single = runtime.next(0.0);
  check(single.count == 1 && close(single.events[0].pitchVolts, 2.f / 12.f),
        "D@4 is an individual absolute note, not a chord");

  const auto dominant = runtime.next(1.0);
  const float dominantExpected[] = {-10.f / 12.f, -6.f / 12.f, -3.f / 12.f,
                                    0.f};
  check(dominant.count == 4, "D7 expands to a four-voice dominant seventh");
  for (std::size_t voice = 0; voice < dominant.count; ++voice)
    check(dominant.events[voice].voice == voice &&
              dominant.events[voice].voiceCount == 4 &&
              close(dominant.events[voice].pitchVolts, dominantExpected[voice]),
          "jazz chord voice " + std::to_string(voice) + " expected " +
              std::to_string(dominantExpected[voice]) + " V, got " +
              std::to_string(dominant.events[voice].pitchVolts) + " V");

  const auto minor = runtime.next(2.0);
  const float minorExpected[] = {-1.f, -9.f / 12.f, -5.f / 12.f, -2.f / 12.f};
  check(minor.count == 4, "Cm7 expands to a minor seventh chord");
  for (std::size_t voice = 0; voice < minor.count; ++voice)
    check(close(minor.events[voice].pitchVolts, minorExpected[voice]),
          "minor seventh chord intervals are voiced in close position");

  const auto slash = runtime.next(3.0);
  check(slash.count == 5 && close(slash.events[0].pitchVolts, -22.f / 12.f),
        "an explicit slash bass becomes polyphonic channel zero");

  const auto degrees = runtime.next(4.0);
  check(degrees.count == 3 && close(degrees.events[0].pitchVolts, 0.f) &&
            close(degrees.events[1].pitchVolts, 3.f / 12.f) &&
            close(degrees.events[2].pitchVolts, 7.f / 12.f),
        "parentheses specify simultaneous scale-degree voices");

  const auto relative = runtime.next(5.0);
  check(relative.count == 1 && close(relative.events[0].pitchVolts, 0.f),
        "apostrophe shifts relative to the independently cycling octave lane");

  const auto choices = tfseq::Compile(R"(choices = sequence {
  cycle 1
  notes <(1 3 5) (2 4 6)> [(1 b3 5)|(2 4 6)]
}
play choices
)");
  check(static_cast<bool>(choices), choices.diagnostic.message);
  if (choices) {
    tfseq::Runtime choiceRuntime;
    choiceRuntime.setProgram(choices.program.get());
    check(choiceRuntime.next(0.0).count == 3 &&
              choiceRuntime.next(1.0).count == 3,
          "alternation and random choice retain nested chord atoms");
  }

  const auto inversion = tfseq::Compile(R"(inversion = sequence {
  tonic C
  octave 4
  notes C/E
  articulation x _
  gate .5
}
play inversion
)");
  check(static_cast<bool>(inversion), inversion.diagnostic.message);
  if (inversion) {
    tfseq::Runtime inversionRuntime;
    inversionRuntime.setProgram(inversion.program.get());
    const auto attack = inversionRuntime.next(0.0);
    check(attack.count == 3 &&
              close(attack.events[0].pitchVolts, -8.f / 12.f) &&
              close(attack.events[1].pitchVolts, 0.f) &&
              close(attack.events[2].pitchVolts, 7.f / 12.f),
          "an implicit slash bass inverts below the close-position chord");
    bool everyVoiceLegato = attack.count == 3;
    for (std::size_t voice = 0; voice < attack.count; ++voice)
      everyVoiceLegato = everyVoiceLegato && attack.events[voice].legatoToNext;
    check(everyVoiceLegato,
          "a following tie extends every voice of the current chord");
  }

  check(!tfseq::Compile(R"(too_wide = sequence {
  notes (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 / C)
}
play too_wide
)"),
        "slash chords respect Rack's 16-channel cable limit");

  const auto namedVoicing = tfseq::Compile(R"(named = sequence {
  notes (C E G)@4
}
play named
)");
  check(static_cast<bool>(namedVoicing),
        "named voicing: " + namedVoicing.diagnostic.message);
  if (namedVoicing) {
    tfseq::Runtime namedRuntime;
    namedRuntime.setProgram(namedVoicing.program.get());
    const auto chord = namedRuntime.next(0.0);
    check(chord.count == 3 && close(chord.events[0].pitchVolts, 0.f) &&
              close(chord.events[1].pitchVolts, 4.f / 12.f) &&
              close(chord.events[2].pitchVolts, 7.f / 12.f),
          "bare named tones share an explicit parenthesized chord register");
    bool oneCursorPerChord = chord.count == 3;
    for (std::size_t voice = 1; voice < chord.count; ++voice) {
      const auto &first =
          chord.events[0]
              .cursors[static_cast<std::size_t>(tfseq::CursorLane::Notes)];
      const auto &current =
          chord.events[voice]
              .cursors[static_cast<std::size_t>(tfseq::CursorLane::Notes)];
      oneCursorPerChord = oneCursorPerChord && first.begin == current.begin &&
                          first.end == current.end;
    }
    check(oneCursorPerChord,
          "all chord voices report one coherent musical cursor span");
  }
}

void distinguishInKeyShiftsFromModulation() {
  auto compiled = tfseq::Compile(R"(base = sequence {
  cycle 1
  tonic C@4
  scale major
  notes 7
}

in_key = base |> shift_degree 4
modulated = base |> modulate_degree 5
song = in_key + modulated
play song
)");
  check(static_cast<bool>(compiled),
        "ordered modulation: " + compiled.diagnostic.message);
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  check(close(runtime.next(0).events[0].pitchVolts, 17.f / 12.f),
        "shift_degree remaps the note inside the current key");
  check(close(runtime.next(1).events[0].pitchVolts, 18.f / 12.f),
        "modulate_degree preserves the riff shape at a new tonal centre");
}

void scaleAndModulationPipelinesComposeLeftToRight() {
  const auto compiled = tfseq::Compile(R"(base = sequence {
  cycle 1
  tonic C@4
  scale minor
  notes 1
}
minor_third = base |> modulate_degree 3 |> scale major
major_third = base |> scale major |> modulate_degree 3
song = minor_third + major_third
play song
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  check(close(runtime.next(0.0).events[0].pitchVolts, 3.f / 12.f),
        "modulation captures the minor scale active before scale replacement");
  check(close(runtime.next(1.0).events[0].pitchVolts, 4.f / 12.f),
        "a preceding scale replacement affects subsequent modulation");
}

void duplicateNamesAreRejectedInEitherOrder() {
  check(!tfseq::Compile(R"(same = sequence {
  notes 1
}
same = same
play same
)"),
        "an assignment cannot silently shadow an earlier sequence");
  check(!tfseq::Compile(R"(same = other
same = sequence {
  notes 1
}
play same
)"),
        "a sequence cannot shadow an earlier assignment");
}

void settingLanePipelinesAreNeverIgnored() {
  const auto cycle = tfseq::Compile(R"(a = sequence {
  cycle 8 |> fast 2
  notes 1
}
play a
)");
  check(!cycle && cycle.diagnostic.message.find("does not accept") !=
                      std::string::npos,
        "a pipeline attached to cycle is rejected instead of discarded");
  const auto scale = tfseq::Compile(R"(a = sequence {
  scale minor |> rev
  notes 1
}
play a
)");
  check(!scale &&
            scale.diagnostic.message.find("closing brace") != std::string::npos,
        "setting diagnostics explain where sequence transforms belong");
}

void harmonicMinorBuildsMajorDominantInKey() {
  auto compiled = tfseq::Compile(R"(tonic_triad = sequence {
  cycle 3
  tonic C@4
  scale harmonic_minor
  notes 1 3 5
}
v = tonic_triad |> shift_degree 4
play v
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  check(close(runtime.next(0).events[0].pitchVolts, 7.f / 12.f),
        "harmonic-minor in-key shift places the root on degree five");
  check(close(runtime.next(1).events[0].pitchVolts, 11.f / 12.f),
        "harmonic-minor raised seventh makes the dominant third major");
  check(close(runtime.next(2).events[0].pitchVolts, 14.f / 12.f),
        "harmonic-minor in-key shift completes the major dominant triad");
}

void stopIsAFirstClassTransportCommand() {
  auto stopped = tfseq::Compile("stop\n");
  check(static_cast<bool>(stopped), stopped.diagnostic.message);
  if (!stopped)
    return;
  check(stopped.program->semantic().stopped,
        "stop compiles without a sequence definition");
  check(stopped.program->semantic().arrangement.empty(),
        "stop schedules no arrangement");

  auto lastCommandWins = tfseq::Compile(R"(riff = sequence {
  notes 1
}
play riff
stop
)");
  check(static_cast<bool>(lastCommandWins), lastCommandWins.diagnostic.message);
  check(lastCommandWins && lastCommandWins.program->semantic().stopped,
        "a trailing stop overrides play");
}

void conciseAcidSyntax() {
  auto compiled = tfseq::Compile(R"(acid = sequence {
  cycle 16
  tonic D#@2
  scale minor
  notes 1!4 5 7 1!2 8 1, 1 1, 1
  articulation x!7 ~ > ~ >!2 ~ >!3
  velocity .5
  accent + .!3 + . + . + .!2 + .
  gate .5
  glide .8
}
iv = acid |> modulate_degree 4 |> octave -1
v = acid |> modulate_degree 5 |> octave -1 |> scale major
song = acid * 8 + iv * 4 + v * 4
play song
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  check(compiled.program->semantic().sequences.size() == 3,
        "derived acid sections reuse one concise sequence");
  check(compiled.program->semantic().arrangement.size() == 3,
        "acid sections form a three-part song");
  check(compiled.program->semantic().sequences[0].notes.size() == 13,
        "note repetition expands only the pitched events");
  check(compiled.program->semantic().sequences[0].accent.size() == 13,
        "sparse scalar repetition follows pitched events");

  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const auto first = runtime.next(0);
  check(first.count == 1 && close(first.events[0].pitchVolts, -21.f / 12.f),
        "acid starts on D-sharp 2");
  check(close(first.events[0].velocity, .88f), "first acid step is accented");
  check(close(first.events[0].accent, .88f),
        "accent remains independently available from velocity");
  for (int beat = 1; beat < 7; ++beat)
    runtime.next(beat);
  check(runtime.next(7).events[0].kind == tfseq::EventKind::Rest,
        "compact replicated articulation retains the first rest");
  const auto slide = runtime.next(8);
  check(slide.events[0].kind == tfseq::EventKind::Slide &&
            close(slide.events[0].slideBeats, .8f),
        "greater-than uses the sequence glide time");

  runtime.reset();
  int acidCursor = -1;
  int ivCursor = -1;
  const int minor[] = {-21, -21, -21, -21, -14, -11, -21, 0,
                       -21, 0,   -9,  -33, 0,   -21, -33, -21};
  const int major[] = {-26, -26, -26, -26, -19, -15, -26, 0,
                       -26, 0,   -14, -38, 0,   -26, -38, -26};
  for (int beat = 0; beat < 256; ++beat) {
    const int position = beat % 16;
    const bool rest = position == 7 || position == 9 || position == 12;
    const int transpose = beat < 128 ? 0 : beat < 192 ? -7 : 0;
    const int expected =
        beat < 192 ? minor[position] + transpose : major[position];
    const auto events = runtime.next(beat);
    check(events.count == 1, "each acid slot produces one bounded event");
    if (events.count != 1)
      continue;
    const int sequenceCursor =
        events.events[0]
            .cursors[static_cast<std::size_t>(tfseq::CursorLane::Sequence)]
            .begin;
    if (beat == 0)
      acidCursor = sequenceCursor;
    else if (beat == 128)
      ivCursor = sequenceCursor;
    if (rest) {
      check(events.events[0].kind == tfseq::EventKind::Rest,
            "acid rests survive all derived sections");
    } else {
      check(close(events.events[0].pitchVolts,
                  static_cast<float>(expected) / 12.f),
            "derived acid sections reproduce the original pitch song");
      check(close(events.events[0].gateFraction, .5f),
            "acid smoke pattern preserves its half-step gate");
    }
  }
  check(acidCursor >= 0 && ivCursor >= 0 && acidCursor != ivCursor,
        "the play cursor follows each named arrangement section");
}

void hotSwapPreservesNamedSequencePhase() {
  auto original = tfseq::Compile(R"(riff = sequence {
  cycle 4
  notes 1 2 3 4
}
play riff
)");
  auto edited = tfseq::Compile(R"(riff = sequence {
  cycle 4
  notes 1 2 7 4
}
play riff
)");
  check(static_cast<bool>(original), original.diagnostic.message);
  check(static_cast<bool>(edited), edited.diagnostic.message);
  if (!original || !edited)
    return;

  tfseq::Runtime runtime;
  runtime.setProgram(original.program.get());
  runtime.next(0.0);
  runtime.next(1.0);
  runtime.replaceProgram(edited.program.get(), 2.0);
  const auto next = runtime.next(2.0);
  check(next.count == 1 && close(next.events[0].pitchVolts, 11.f / 12.f),
        "hot swap continues at the current named-sequence cursor");
}

void preparedWorkspaceHasNoSmallEventCeiling() {
  std::string subdivisions = "[";
  for (int index = 0; index < 40; ++index)
    subdivisions += index == 0 ? "x" : " x";
  subdivisions += "]";
  const auto compiled =
      tfseq::Compile("dense = sequence {\n  notes 1 2 3 4\n  articulation " +
                     subdivisions + "\n}\nplay dense\n");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  check(compiled.program->maximumEventsPerStep >= 40,
        "the compiler prepares workspace from actual event density");
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const auto events = runtime.next(0.0);
  check(events.count == 40 && !events.overflowed,
        "a legal step above the former 32-event ceiling is not truncated");
}

void playbackStateHasNoFixedSequenceCeiling() {
  std::string source;
  for (int index = 0; index < 40; ++index) {
    source += "s" + std::to_string(index) + " = sequence {\n";
    source += "  cycle 1\n";
    source += "  notes " + std::to_string(index + 1) + "\n}\n";
  }
  source += "song = ";
  for (int index = 0; index < 40; ++index) {
    if (index != 0)
      source += " + ";
    source += "s" + std::to_string(index);
  }
  source += "\nplay song\n";
  const auto compiled = tfseq::Compile(source);
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  check(compiled.program->stateWorkspace.size() == 40,
        "playback state is prepared for every compiled sequence");
  const auto &sequences = compiled.program->semantic().sequences;
  const auto &order = compiled.program->stateTransferOrder;
  bool transferOrderIsSorted = order.size() == sequences.size();
  for (std::size_t position = 1;
       transferOrderIsSorted && position < order.size(); ++position) {
    const auto &left = sequences[order[position - 1]];
    const auto &right = sequences[order[position]];
    transferOrderIsSorted =
        left.stableId < right.stableId ||
        (left.stableId == right.stableId && left.name < right.name);
  }
  check(transferOrderIsSorted,
        "the compiler prepares a collision-safe linear state-transfer index");
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  for (int beat = 0; beat < 40; ++beat) {
    const auto events = runtime.next(static_cast<double>(beat));
    check(events.count == 1,
          "arrangements beyond the old sequence ceiling keep playing");
  }
}

void longDurationLegatoDoesNotDependOnSchedulerLookahead() {
  const auto compiled = tfseq::Compile(R"(legato = sequence {
  notes 1 3
  articulation x _ >
  duration 2
  gate .5
}
play legato
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const auto attack = runtime.next(0.0);
  check(attack.count == 1 && attack.events[0].legatoToNext,
        "an attack knows that a multi-beat tie follows it");
  const auto tie = runtime.next(2.0);
  check(tie.count == 1 && tie.events[0].kind == tfseq::EventKind::Tie,
        "the following tie remains scheduled at its true duration boundary");
}

void nestedArrangementsAndCompositeTransformsAreTyped() {
  const auto compiled = tfseq::Compile(R"(a = sequence {
  cycle 1
  notes 1
}
b = sequence {
  cycle 1
  notes 3
}
phrase = a + b
song = (phrase * 2 + a) |> transpose_semitone 12
play song
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  check(compiled.program->semantic().arrangement.size() == 5,
        "named and grouped arrangements compose recursively");
  check(compiled.program->semantic().sequences.size() == 4,
        "a transformed repeated phrase reuses one clone per source sequence");
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const float expected[] = {1.f, 1.f + 4.f / 12.f, 1.f, 1.f + 4.f / 12.f, 1.f};
  for (int beat = 0; beat < 5; ++beat) {
    const auto events = runtime.next(beat);
    check(events.count == 1 &&
              close(events.events[0].pitchVolts, expected[beat]),
          "a pipeline on a composition transforms each component");
  }

  const auto trailing = tfseq::Compile(R"(a = sequence {
  notes 1 |> rev ignored
}
play a
)");
  check(!trailing,
        "typed transform parsing rejects ignored trailing arguments");
  const auto wrongLane = tfseq::Compile(R"(a = sequence {
  notes 1
  accent + . |> transpose_semitone 2
}
play a
)");
  check(!wrongLane,
        "a pitch transform on a scalar lane is a diagnostic, not a no-op");
}

void sequenceConditionsStayCoherentAcrossLanes() {
  const auto compiled = tfseq::Compile(R"(coherent = sequence {
  cycle 1
  notes 1 3
  velocity .1 .9
}
|> sometimes .5 rev
play coherent
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  for (int beat = 0; beat < 16; ++beat) {
    const auto events = runtime.next(beat);
    if (events.count != 1) {
      check(false, "coherence fixture emits one event per beat");
      continue;
    }
    const bool lowPitch = close(events.events[0].pitchVolts, 0.f);
    check(close(events.events[0].velocity, lowPitch ? .1f : .9f),
          "one sometimes decision is shared by notes and velocity");
  }
}

void quantifiedTimingTransforms() {
  const auto compiled = tfseq::Compile(R"(timed = sequence {
  notes 1 2
  articulation [x x]
  duration 1
}
|> fast 2
|> swing .6
|> late 1/16
|> early random 8ms
play timed
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const auto events = runtime.next(0.0);
  check(close(static_cast<float>(events.durationBeats), .5f),
        "fast requires and applies an explicit numeric factor");
  check(events.count == 2 &&
            close(static_cast<float>(events.events[0].beat), 1.f / 16.f) &&
            close(static_cast<float>(events.events[0].spanBeats), .3f) &&
            close(static_cast<float>(events.events[1].beat), .3625f) &&
            close(static_cast<float>(events.events[1].spanBeats), .2f),
        "swing and beat-relative late offsets retain the step duration");
  check(events.events[0].timingOffsetMilliseconds <= 0.0 &&
            events.events[0].timingOffsetMilliseconds >= -8.0,
        "random millisecond offsets remain in their declared range");

  const auto slow = tfseq::Compile(R"(a = sequence {
  notes 1
}
|> slow 3/2
|> late 12ms
play a
)");
  check(static_cast<bool>(slow), slow.diagnostic.message);
  if (slow) {
    tfseq::Runtime slowRuntime;
    slowRuntime.setProgram(slow.program.get());
    const auto event = slowRuntime.next(0.0);
    check(
        close(static_cast<float>(event.durationBeats), 1.5f) &&
            close(static_cast<float>(event.events[0].timingOffsetMilliseconds),
                  12.f),
        "slow fractions and absolute millisecond offsets are explicit");
  }
}

void timingPreparationCoversEarlyLookaheadAndMilliseconds() {
  const auto timed = tfseq::Compile(R"(a = sequence {
  notes 1
  duration 3/2
}
|> early 1
|> early 8ms
|> late 100ms
play a
)");
  check(static_cast<bool>(timed),
        "timing preparation: " + timed.diagnostic.message);
  if (!timed)
    return;
  check(close(static_cast<float>(timed.program->maximumEarlyBeats), 1.f) &&
            close(static_cast<float>(timed.program->maximumEarlyMilliseconds),
                  8.f),
        "the compiler retains worst-case early timing bounds");
  const auto lookahead =
      tfseq::SchedulingLookaheadBeats(*timed.program, true, 24000.0, 48000.0);
  check(close(static_cast<float>(lookahead), 2.016f),
        "scheduler lookahead includes beat and millisecond early offsets");

  const auto plain = tfseq::Compile(R"(a = sequence {
  notes 1
  duration 3/2
}
play a
)");
  check(plain &&
            timed.program->scheduleCapacity > plain.program->scheduleCapacity,
        "millisecond lateness contributes to prepared queue capacity");
}

void pegFrontendOwnsVoicingAndSlashStructure() {
  const auto parsed = tfseq::syntax::Parse(R"(harmony = sequence {
  notes Bbm7b5@3 / D@2 (1 b3 5)@4 (C E G / B)@3
}
play harmony
)");
  check(static_cast<bool>(parsed), parsed.diagnostic.message);
  if (!parsed)
    return;
  const auto *sequence = std::get_if<tfseq::syntax::SequenceDefinition>(
      &parsed.document.statements.front());
  check(sequence && sequence->lanes.front().pattern.steps.size() == 3,
        "spaced slash chords remain one PEG pattern element");
  if (!sequence)
    return;
  const auto &steps = sequence->lanes.front().pattern.steps;
  check(steps[0].kind == tfseq::syntax::PatternKind::Slash &&
            steps[0].children.size() == 2 &&
            steps[0].children[0].atom.text == "Bbm7b5@3" &&
            steps[0].children[1].atom.text == "D@2",
        "PEG emits the jazz chord and slash bass relationship directly");
  check(steps[1].kind == tfseq::syntax::PatternKind::Voicing &&
            steps[1].children.size() == 3 && steps[1].suffix.text == "@4",
        "PEG emits explicit tones and their shared register directly");
  check(steps[2].kind == tfseq::syntax::PatternKind::Slash &&
            steps[2].children[0].kind == tfseq::syntax::PatternKind::Voicing,
        "a slash inside an explicit voicing remains structural AST");
}

void hotSwapCheckpointSurvivesLookahead() {
  auto original = tfseq::Compile(R"(riff = sequence {
  notes 1 2 3 4 5 6 7 8
  duration 1/4
}
play riff
)");
  auto edited = tfseq::Compile(R"(riff = sequence {
  notes 1 2 3 4 8 6 7 8
  duration 1/4
}
play riff
)");
  check(original && edited, "sub-beat checkpoint fixtures compile");
  if (!original || !edited)
    return;
  tfseq::Runtime runtime;
  tfseq::Runtime activation;
  runtime.setProgram(original.program.get());
  for (double beat : {0.0, 0.25, 0.5, 0.75})
    runtime.next(beat);
  runtime.snapshot(activation);
  runtime.next(1.0);
  runtime.next(1.25);
  runtime = activation;
  runtime.replaceProgram(edited.program.get(), 1.0);
  const auto next = runtime.next(1.0);
  check(next.count == 1 && close(next.events[0].pitchVolts, 1.f),
        "a clock-boundary snapshot ignores counters advanced by lookahead");
}

void swingCrossesStepBoundariesAtExplicitSubdivisions() {
  const auto eighths = tfseq::Compile(R"(a = sequence {
  notes 1 2 3 4
  duration 1/8
}
|> swing .6 1/8
play a
)");
  check(static_cast<bool>(eighths), eighths.diagnostic.message);
  if (eighths) {
    tfseq::Runtime runtime;
    runtime.setProgram(eighths.program.get());
    const auto first = runtime.next(0.0);
    const auto firstSpan = first.events[0].spanBeats;
    const auto second = runtime.next(0.125);
    check(close(static_cast<float>(firstSpan), .15f) &&
              close(static_cast<float>(second.events[0].beat), .15f) &&
              close(static_cast<float>(second.events[0].spanBeats), .10f),
          "eighth-beat swing pairs events returned by separate runtime steps");
  }

  const auto sixteenths = tfseq::Compile(R"(a = sequence {
  notes 1 2
  duration 1/16
}
|> swing .625 1/16
play a
)");
  check(static_cast<bool>(sixteenths), sixteenths.diagnostic.message);
  if (sixteenths) {
    tfseq::Runtime runtime;
    runtime.setProgram(sixteenths.program.get());
    runtime.next(0.0);
    const auto second = runtime.next(0.0625);
    check(close(static_cast<float>(second.events[0].beat), .078125f),
          "the swing grid accepts a sixteenth-beat subdivision");
  }
}

void patternExpansionUsesAddressableLimits() {
  const auto repeated = tfseq::Compile(R"(a = sequence {
  notes 1!65
  articulation x(3,65)
  glide 65
  ratchet 17
}
play a
)");
  check(static_cast<bool>(repeated), repeated.diagnostic.message);
  if (!repeated)
    return;
  check(repeated.program->semantic().sequences[0].notes.size() == 65 &&
            repeated.program->semantic().sequences[0].articulation.size() ==
                65 &&
            repeated.program->maximumEventsPerStep >= 17,
        "former fixed expansion ceilings are derived from allocation");
}

void unsafeNumericInputsAreDiagnostics() {
  check(!tfseq::Compile("seed .5\nstop\n"),
        "fractional seeds are rejected rather than narrowed");
  check(!tfseq::Compile("a = sequence {\n notes 1\n duration 0\n}\nplay a\n"),
        "zero duration is rejected before playback");
  check(!tfseq::Compile("a = sequence {\n notes 1\n gate 2\n}\nplay a\n"),
        "out-of-range gate values are rejected rather than clamped");
  check(!tfseq::Compile(
            "a = sequence {\n notes 999999999999999999999\n}\nplay a\n"),
        "overflowing scale degrees are rejected safely");
}

} // namespace

int main() {
  pegFrontendBuildsTypedSyntax();
  pegFrontendRejectsMalformedStructure();
  typedPatternTreeOwnsReusableStructure();
  pegFrontendOwnsVoicingAndSlashStructure();
  lineCommentsCanTruncatePatterns();
  selectionEvaluationReplacesOnlyContainingStatements();
  explicitSingleRepeatRemainsAnArrangement();
  compileAndCycleIndependentLanes();
  parseFirstClassArticulation();
  articulationModifiersComposeInEitherOrder();
  transformByCycleAndArrange();
  rejectInvalidInput();
  degreesContinueAcrossOctaves();
  scaleCardinalityAndExplicitOctaves();
  absoluteAndRelativeRegistersRemainDistinct();
  octaveSuffixesAndChordsAreUnambiguous();
  distinguishInKeyShiftsFromModulation();
  scaleAndModulationPipelinesComposeLeftToRight();
  duplicateNamesAreRejectedInEitherOrder();
  settingLanePipelinesAreNeverIgnored();
  harmonicMinorBuildsMajorDominantInKey();
  stopIsAFirstClassTransportCommand();
  conciseAcidSyntax();
  hotSwapPreservesNamedSequencePhase();
  hotSwapCheckpointSurvivesLookahead();
  preparedWorkspaceHasNoSmallEventCeiling();
  playbackStateHasNoFixedSequenceCeiling();
  longDurationLegatoDoesNotDependOnSchedulerLookahead();
  nestedArrangementsAndCompositeTransformsAreTyped();
  sequenceConditionsStayCoherentAcrossLanes();
  quantifiedTimingTransforms();
  swingCrossesStepBoundariesAtExplicitSubdivisions();
  patternExpansionUsesAddressableLimits();
  timingPreparationCoversEarlyLookaheadAndMilliseconds();
  unsafeNumericInputsAreDiagnostics();
  if (failures != 0)
    std::cerr << failures << " text sequencer test(s) failed\n";
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
