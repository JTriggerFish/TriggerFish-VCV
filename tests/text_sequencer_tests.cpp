#include "tfseq.hpp"
#include "tfseq_editor.hpp"
#include "tfseq_parser.hpp"
#include "tfui_animation.hpp"
#include "tfui_colormap.hpp"

#include <algorithm>
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

void editorShortcutTextOperationsAreStructural() {
  const std::string lanes = "  notes 1\n  gate .5\n";
  const auto commented = tfseq::editor::ToggleLineComments(
      lanes, 0, static_cast<int>(lanes.size()));
  check(commented.changed && commented.text == "  // notes 1\n  // gate .5\n",
        "comment toggle preserves indentation across selected lines");
  const auto restored = tfseq::editor::ToggleLineComments(
      commented.text, commented.cursor, commented.selection);
  check(restored.changed && restored.text == lanes,
        "a second comment toggle restores the selected lines");

  const std::string windowsLines = "notes 1\r\n\r\ngate .5\r\n";
  const auto windowsCommented = tfseq::editor::ToggleLineComments(
      windowsLines, 0, static_cast<int>(windowsLines.size()));
  check(windowsCommented.changed &&
            windowsCommented.text == "// notes 1\r\n\r\n// gate .5\r\n",
        "comment toggle preserves CRLF and leaves blank rows blank");
  const auto windowsRestored = tfseq::editor::ToggleLineComments(
      windowsCommented.text, windowsCommented.cursor,
      windowsCommented.selection);
  check(windowsRestored.changed && windowsRestored.text == windowsLines,
        "comment toggle restores CRLF source exactly");

  const std::string carriageReturnLines = "notes 1\rgate .5";
  const auto carriageReturnCommented = tfseq::editor::ToggleLineComments(
      carriageReturnLines, 0, static_cast<int>(carriageReturnLines.size()));
  check(carriageReturnCommented.changed &&
            carriageReturnCommented.text == "// notes 1\r// gate .5",
        "comment toggle recognizes lone carriage-return line endings");

  const std::string program = "notes 1\nplay a\n";
  const auto duplicatedLine = tfseq::editor::Duplicate(program, 3, 3);
  check(duplicatedLine.changed &&
            duplicatedLine.text == "notes 1\nnotes 1\nplay a\n",
        "duplicate repeats the complete current line");
  const auto duplicatedSelection = tfseq::editor::Duplicate("abc", 0, 2);
  check(duplicatedSelection.changed && duplicatedSelection.text == "ababc" &&
            duplicatedSelection.selection == 2 &&
            duplicatedSelection.cursor == 4,
        "duplicate repeats and selects the explicit selection");
  const auto duplicatedWindowsLine =
      tfseq::editor::Duplicate("notes 1\r\nplay a", 12, 12);
  check(duplicatedWindowsLine.changed &&
            duplicatedWindowsLine.text == "notes 1\r\nplay a\r\nplay a",
        "duplicate preserves CRLF when the final line has no terminator");
  const auto duplicatedCarriageReturnLine =
      tfseq::editor::Duplicate("notes 1\rplay a", 2, 2);
  check(duplicatedCarriageReturnLine.changed &&
            duplicatedCarriageReturnLine.text == "notes 1\rnotes 1\rplay a",
        "duplicate preserves lone carriage-return line endings");
}

void heatmapMapsScalarIntensity() {
  const auto dark = tfui::sampleHeatmap(tfui::ProgramEditorHeatmap, -1.f);
  const auto bright = tfui::sampleHeatmap(tfui::ProgramEditorHeatmap, 2.f);
  check(close(dark.red, tfui::MagmaHeatmap.front().red) &&
            close(dark.green, tfui::MagmaHeatmap.front().green) &&
            close(dark.blue, tfui::MagmaHeatmap.front().blue),
        "heatmap clamps intensity below zero to its dark endpoint");
  check(close(bright.red, tfui::MagmaHeatmap.back().red) &&
            close(bright.green, tfui::MagmaHeatmap.back().green) &&
            close(bright.blue, tfui::MagmaHeatmap.back().blue),
        "heatmap clamps intensity above one to its bright endpoint");

  for (const auto palette : {
           tfui::HeatmapPalette::Magma, tfui::HeatmapPalette::Inferno,
           tfui::HeatmapPalette::Plasma, tfui::HeatmapPalette::Viridis,
           tfui::HeatmapPalette::Cividis, tfui::HeatmapPalette::CrtGreen,
           tfui::HeatmapPalette::CrtBlue, tfui::HeatmapPalette::CrtYellow,
           tfui::HeatmapPalette::CrtRed}) {
    float previousLuminance = -1.f;
    for (int step = 0; step <= 100; ++step) {
      const auto color = tfui::sampleHeatmap(
          palette, static_cast<float>(step) / 100.f);
      const float luminance =
          0.2126f * color.red + 0.7152f * color.green + 0.0722f * color.blue;
      check(luminance + 1e-6f >= previousLuminance,
            "increasing editor intensity never lowers heatmap luminance");
      previousLuminance = luminance;
    }
  }

  const auto viridis = tfui::sampleHeatmap(tfui::HeatmapPalette::Viridis, 0.f);
  check(close(viridis.red, tfui::ViridisHeatmap.front().red) &&
            close(viridis.green, tfui::ViridisHeatmap.front().green) &&
            close(viridis.blue, tfui::ViridisHeatmap.front().blue),
        "named heatmaps select their own color table");
  check(tfui::heatmapPaletteFromInt(3) == tfui::HeatmapPalette::Viridis &&
            tfui::heatmapPaletteFromInt(8) == tfui::HeatmapPalette::CrtRed &&
            tfui::heatmapPaletteFromInt(99) == tfui::HeatmapPalette::Magma,
        "saved heatmap values restore valid choices and default safely");

  for (const auto palette : {
           tfui::HeatmapPalette::CrtGreen, tfui::HeatmapPalette::CrtBlue,
           tfui::HeatmapPalette::CrtYellow, tfui::HeatmapPalette::CrtRed}) {
    const auto unlit = tfui::sampleHeatmap(palette, 0.f);
    check(close(unlit.red, 0.f) && close(unlit.green, 0.f) &&
              close(unlit.blue, 0.f),
          "CRT heatmaps begin at an unlit black screen");
  }
  const auto green = tfui::sampleHeatmap(tfui::HeatmapPalette::CrtGreen, 0.7f);
  const auto blue = tfui::sampleHeatmap(tfui::HeatmapPalette::CrtBlue, 0.7f);
  const auto yellow =
      tfui::sampleHeatmap(tfui::HeatmapPalette::CrtYellow, 0.7f);
  const auto red = tfui::sampleHeatmap(tfui::HeatmapPalette::CrtRed, 0.7f);
  check(green.green > green.red && green.green > green.blue &&
            blue.blue > blue.red && blue.blue > blue.green &&
            yellow.red > yellow.blue && yellow.green > yellow.blue &&
            red.red > red.green && red.red > red.blue,
        "CRT heatmaps retain their named phosphor hue below the bright bloom");
}

void cursorAnimationIsFrameIndependentAndTempoBounded() {
  check(tfui::arrangementCursorGroup(0.0, 4) == 0 &&
            tfui::arrangementCursorGroup(3.999, 4) == 0 &&
            tfui::arrangementCursorGroup(4.0, 4) == 1,
        "arrangement feedback groups four incoming clocks by default");
  check(close(static_cast<float>(tfui::cursorTravelDuration(1.0)), 0.12f),
        "slow pulses use the absolute cursor-travel ceiling");
  check(close(static_cast<float>(tfui::cursorTravelDuration(0.1)), 0.075f),
        "cursor travel uses 75 percent of a fast pulse interval");
  check(close(static_cast<float>(tfui::cursorTravelDuration(0.001)), 0.035f),
        "cursor travel remains visible for multiple display frames");
  check(close(static_cast<float>(tfui::cursorMotionTailDuration(0.1)), 0.11f),
        "motion persistence scales with its lane's pulse interval");
  check(close(static_cast<float>(tfui::cursorMotionTailDuration(1.0)), 0.30f),
        "motion persistence has an absolute ceiling");
  check(close(static_cast<float>(tfui::cursorBloomExpansionDuration(0.1)),
              0.085f),
        "stationary bloom expansion follows the pulse interval");
  check(close(static_cast<float>(tfui::cursorBloomExpansionDuration(0.001)),
              0.070f),
        "stationary blooms remain visible for multiple display frames");
  check(
      close(static_cast<float>(tfui::cursorBloomTailDuration(0.1)), 0.15f) &&
          close(static_cast<float>(tfui::cursorBloomTailDuration(1.0)), 0.42f),
      "stationary bloom persistence is tempo-scaled and bounded");
  check(close(tfui::cursorBloomExpansion(0.0, 0.1), 0.f) &&
            close(tfui::cursorBloomExpansion(0.05, 0.1), 0.5f) &&
            close(tfui::cursorBloomExpansion(0.1, 0.1), 1.f),
        "stationary bloom radius expands smoothly from its source");
  check(
      close(tfui::cursorTravelPosition(tfui::CursorTravelCurve::Linear, 0.25f),
            0.25f),
      "linear travel has constant spatial progress");
  check(close(tfui::cursorTravelPosition(tfui::CursorTravelCurve::Smoothstep,
                                         0.25f),
              0.15625f),
        "smoothstep travel eases spatial progress only");

  const double duration = tfui::cursorTravelDuration(0.2);
  const double tailDuration = tfui::cursorMotionTailDuration(0.2);
  check(close(tfui::cursorMotionIntensity(0.0, duration, tailDuration), 0.65f),
        "the travelling beam is visible immediately at its source");
  check(tfui::cursorMotionIntensity(duration * 0.5, duration, tailDuration) >
            0.99f,
        "the travelling beam reaches a bright midpoint");
  check(close(tfui::cursorMotionIntensity(duration, duration, tailDuration),
              0.65f),
        "motion joins its exponential persistence continuously at arrival");
  check(tfui::cursorMotionIntensity(duration + tailDuration, duration,
                                    tailDuration) < 0.25f,
        "motion persistence is derived directly from timestamp age");
  check(close(tfui::activeStepProgress(3.0, 2.0, 6.0), .25f) &&
            close(tfui::activeStepProgress(1.0, 2.0, 6.0), 0.f) &&
            close(tfui::activeStepProgress(7.0, 2.0, 6.0), 1.f) &&
            close(tfui::activeStepProgress(3.0, 2.0, 2.0), 0.f),
        "active-step progress is duration-aware, clamped, and safe");
  check(close(tfui::transportLightBrightness(true, false), .16f) &&
            close(tfui::transportLightBrightness(true, true), 1.f) &&
            close(tfui::transportLightBrightness(false, true), 0.f) &&
            close(tfui::transportLightBrightness(false, false), 0.f),
        "the transport light distinguishes an armed transport from a clock "
        "heartbeat");
}

void pegFrontendBuildsTypedSyntax() {
  const std::string source = R"(acid = sequence {
  notes 1 <2 b3> [4|5] 1'!2 6(3,8,1) |> rotate 1 |> every 4 rev
}
|> sometimes .25 rev
song = acid * 8 + acid
  |> fast 2
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
  check(sequence && sequence->lanes.size() == 1,
        "PEG front end retains the note-first sequence lane");
  if (sequence) {
    const auto &pattern = sequence->lanes[0].pattern;
    check(
        pattern.steps.size() == 5 &&
            pattern.steps[1].kind == tfseq::syntax::PatternKind::CycleChoice &&
            pattern.steps[2].kind == tfseq::syntax::PatternKind::RandomChoice &&
            pattern.steps[3].kind == tfseq::syntax::PatternKind::Event &&
            pattern.steps[3].repeatCount.text == "2",
        "PEG front end produces typed choice and repetition nodes");
    check(sequence->lanes[0].pipelines.size() == 2,
          "PEG pattern grammar separates inline pipelines");
    check(pattern.steps[4].kind == tfseq::syntax::PatternKind::Event &&
              pattern.steps[4].arguments.size() == 3,
          "Euclidean suffixes are typed on complete note events");
    check(sequence->pipelines.size() == 1,
          "PEG document grammar retains sequence pipelines");
  }
  const auto *assignment =
      std::get_if<tfseq::syntax::Assignment>(&parsed.document.statements[1]);
  check(assignment && assignment->expression.terms.size() == 2 &&
            assignment->expression.terms[0].repeats == 8 &&
            assignment->expression.pipelines.size() == 1,
        "PEG assignment grammar parses repetition and pipeline continuation");
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
  check(!tfseq::syntax::Parse(
            "a = sequence {\n  notes 1\n}\nsong = a * 0\nplay song\n"),
        "arrangement repetition is positive in the PEG itself");
}

void typedPatternTreeOwnsReusableStructure() {
  const auto nested = tfseq::syntax::Parse(R"(a = sequence {
  notes <1 [2|3]>!2
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
                tfseq::syntax::PatternKind::CycleChoice &&
            sequence->lanes[0].pattern.steps[0].repeatCount.text == "2" &&
            sequence->lanes[0].pattern.steps[0].children[1].kind ==
                tfseq::syntax::PatternKind::RandomChoice,
        "nested grouping has one domain-neutral syntax tree");

  const auto articulation = tfseq::Compile(R"(a = sequence {
  notes [1 1!2]!2
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

  const auto playBegin = static_cast<int>(draft.find("play a"));
  std::string alternateDraft = draft;
  alternateDraft.replace(static_cast<std::size_t>(playBegin), 6, "play b");
  const auto alternateDocument = tfseq::syntax::Parse(alternateDraft);
  const auto alternateContext =
      alternateDocument
          ? tfseq::syntax::MergeSelectionDocuments(
                evaluatedDocument.document, evaluated,
                alternateDocument.document, alternateDraft, playBegin,
                playBegin + 6)
          : tfseq::syntax::SelectionDocumentResult{};
  const auto alternate =
      alternateContext ? tfseq::Compile(alternateContext.document)
                       : tfseq::Compile("");
  const bool playsB = alternate &&
                      !alternate.program->semantic().arrangement.empty() &&
                      alternate.program->semantic()
                              .sequences[alternate.program->semantic()
                                             .arrangement.front()
                                             .sequence]
                              .name == "b";
  check(alternateContext && playsB,
        "a selected play statement replaces the active arrangement selection");

  const std::string brokenElsewhere = R"(a = sequence {
  notes 9 10
}
|> octave 1
b = sequence {
  notes [7 8
}
play a
)";
  const auto selectedBegin =
      static_cast<int>(brokenElsewhere.find("notes 9 10"));
  check(!tfseq::syntax::Parse(brokenElsewhere),
        "the complete unrelated draft is intentionally invalid");
  const auto selectedDocument = tfseq::syntax::ParseStatementsContaining(
      brokenElsewhere, selectedBegin, selectedBegin);
  check(selectedDocument && selectedDocument.document.statements.size() == 1,
        "the shared PEG isolates a valid containing statement");
  if (selectedDocument) {
    const auto *selectedSequence =
        std::get_if<tfseq::syntax::SequenceDefinition>(
            &selectedDocument.document.statements.front());
    check(selectedSequence && selectedSequence->pipelines.size() == 1,
          "statement isolation retains following sequence continuations");
    const auto isolated = tfseq::syntax::MergeSelectionDocuments(
        evaluatedDocument.document, evaluated, selectedDocument.document,
        brokenElsewhere, selectedBegin, selectedBegin);
    const auto isolatedProgram =
        isolated ? tfseq::Compile(isolated.document) : tfseq::Compile("");
    check(isolated && isolatedProgram,
          "a valid selected edit compiles despite unrelated malformed text");
    if (isolatedProgram) {
      const auto &selected = isolatedProgram.program->semantic().sequences[0];
      check(selected.notes[0].values[0].voices[0].degree == 9,
            "only the selected definition is replaced by isolated evaluation");
    }
  }
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
  tonic C@4
  scale minor
  notes 1 2 b3 4 5 6 b7 8 |> rotate 1
  velocity .88 .61 .61
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
          "three-item velocity lane displaces against notes");
    check(events.events[0]
              .cursors[static_cast<std::size_t>(tfseq::CursorLane::Notes)]
              .valid(),
          "emitted note carries a source cursor");
  }
}

void parseFirstClassArticulation() {
  const std::string source = R"(voice = sequence {
  notes 1 _ >2 ~ 3 [4 5] 6*3 7!2 8(3,8,1)
  duration 1
}
play voice
)";
  auto compiled = tfseq::Compile(source);
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  const auto &art = compiled.program->semantic().sequences[0].articulation;
  check(art.size() == 17, "note events expand replication and Euclidean cells");
  check(art[0].atoms[0].kind == tfseq::ArticulationKind::Attack,
        "a plain note is an attack");
  check(art[1].atoms[0].kind == tfseq::ArticulationKind::Tie,
        "underscore is a tie");
  check(art[2].atoms[0].kind == tfseq::ArticulationKind::Slide,
        "greater-than is a slide");
  check(art[3].atoms[0].kind == tfseq::ArticulationKind::Rest,
        "tilde is a rest");
  check(art[5].atoms.size() == 2, "brackets subdivide one slot");
  check(art[6].atoms[0].ratchets == 3, "asterisk ratchets one slot");

  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  check(runtime.next(0).events[0].kind == tfseq::EventKind::Attack,
        "attack reaches runtime");
  check(runtime.next(1).events[0].kind == tfseq::EventKind::Tie,
        "tie reaches runtime without consuming a note");
  check(runtime.next(2).events[0].kind == tfseq::EventKind::Slide,
        "slide reaches runtime");
  check(runtime.next(3).events[0].kind == tfseq::EventKind::Rest,
        "rest reaches runtime without consuming a note");
  runtime.next(4);
  const auto subdivision = runtime.next(5);
  check(subdivision.count == 2 &&
            close(static_cast<float>(subdivision.events[1].beat), 5.5f),
        "bracketed attacks are evenly subdivided");
  const auto ratchet = runtime.next(6);
  check(ratchet.count == 3, "ratchet emits three attacks in one slot");

  const auto nestedTie = tfseq::Compile(R"(a = sequence {
  notes [1 _] [2 >3]
  gate .25
}
play a
)");
  check(static_cast<bool>(nestedTie), nestedTie.diagnostic.message);
  if (nestedTie) {
    tfseq::Runtime tieRuntime;
    tieRuntime.setProgram(nestedTie.program.get());
    const auto tied = tieRuntime.next(0.0);
    check(tied.count == 2 && tied.events[0].legatoToNext &&
              tied.events[1].kind == tfseq::EventKind::Tie,
          "a source holds Gate into a tie inside one subdivision");
    const auto slid = tieRuntime.next(1.0);
    check(slid.count == 2 && slid.events[0].legatoToNext &&
              slid.events[1].kind == tfseq::EventKind::Slide,
          "a source holds Gate into a slide inside one subdivision");
  }

  const auto ghostMilliseconds = tfseq::Compile(R"(a = sequence {
  notes x1{gate=12ms}
}
play a
)");
  check(static_cast<bool>(ghostMilliseconds),
        ghostMilliseconds.diagnostic.message);
  if (ghostMilliseconds) {
    tfseq::Runtime ghostRuntime;
    ghostRuntime.setProgram(ghostMilliseconds.program.get());
    const auto ghost = ghostRuntime.next(0.0);
    check(ghost.count == 1 && close(ghost.events[0].gateMilliseconds, 12.f) &&
              close(ghost.events[0].gateFraction, .1f) &&
              close(ghost.events[0].gateCapMilliseconds, 20.f),
          "ghost caps preserve millisecond gate units and a short fallback");
  }

  const auto slideMilliseconds = tfseq::Compile(R"(a = sequence {
  notes 1 >2{slide=80ms}
}
play a
)");
  check(static_cast<bool>(slideMilliseconds),
        slideMilliseconds.diagnostic.message);
  if (slideMilliseconds) {
    tfseq::Runtime slideRuntime;
    slideRuntime.setProgram(slideMilliseconds.program.get());
    slideRuntime.next(0.0);
    const auto slide = slideRuntime.next(1.0);
    check(slide.count == 1 && close(slide.events[0].slideMilliseconds, 80.f),
          "exact millisecond slide values override the numerical lane");
  }
}

void namedGateAndDynamicsArticulations() {
  const auto compiled = tfseq::Compile(R"(a = sequence {
  notes 1{quiet} ^2{stacc} ^3{quiet,ten} 4{ten} 5{stacc,gate=.4}
}
play a
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;

  const auto &steps = compiled.program->semantic().sequences[0].articulation;
  check(steps.size() == 5 && steps[0].atoms[0].quiet &&
            steps[1].atoms[0].gateArticulation ==
                tfseq::GateArticulation::Staccato &&
            steps[2].atoms[0].quiet &&
            steps[2].atoms[0].gateArticulation ==
                tfseq::GateArticulation::Tenuto,
        "quiet, stacc, and ten are typed note-owned articulations");

  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const auto quiet = runtime.next(0.0).events[0];
  const auto accentedStaccato = runtime.next(1.0).events[0];
  const auto quietAccentedTenuto = runtime.next(2.0).events[0];
  const auto tenuto = runtime.next(3.0).events[0];
  const auto exactGate = runtime.next(4.0).events[0];
  check(close(quiet.velocity, .36f) && close(quiet.gateFraction, .8f),
        "quiet halves Velocity without shortening Gate");
  check(close(accentedStaccato.velocity, .88f) &&
            close(accentedStaccato.accent, .88f) &&
            close(accentedStaccato.gateFraction, .25f),
        "staccato is orthogonal to accent and Velocity");
  check(close(quietAccentedTenuto.velocity, .44f) &&
            close(quietAccentedTenuto.accent, .88f) &&
            close(quietAccentedTenuto.gateFraction, .95f),
        "quiet can retain a timbral accent while tenuto controls Gate");
  check(close(tenuto.velocity, .72f) && close(tenuto.gateFraction, .95f),
        "tenuto changes Gate without changing Velocity");
  check(close(exactGate.gateFraction, .4f),
        "an inline gate value overrides named gate articulation");

  const auto laneOverride = tfseq::Compile(R"(a = sequence {
  notes 1{stacc} 2{ten}
  gate . .6
}
play a
)");
  check(static_cast<bool>(laneOverride), laneOverride.diagnostic.message);
  if (laneOverride) {
    tfseq::Runtime laneRuntime;
    laneRuntime.setProgram(laneOverride.program.get());
    const auto inherited = laneRuntime.next(0.0).events[0];
    const auto overridden = laneRuntime.next(1.0).events[0];
    check(close(inherited.gateFraction, .25f) &&
              close(overridden.gateFraction, .6f),
          "a lane no-op retains articulation while a value overrides it");
  }

  const auto ratcheted = tfseq::Compile(R"(a = sequence {
  notes 1*2{stacc}
}
play a
)");
  check(static_cast<bool>(ratcheted), ratcheted.diagnostic.message);
  if (ratcheted) {
    tfseq::Runtime ratchetRuntime;
    ratchetRuntime.setProgram(ratcheted.program.get());
    const auto events = ratchetRuntime.next(0.0);
    check(events.count == 2 && close(events.events[0].gateFraction, .25f) &&
              close(events.events[1].gateFraction, .25f),
          "named gate articulation applies to every ratchet subspan");
  }

  const auto slide = tfseq::Compile(R"(a = sequence {
  notes 1 >2{stacc}
}
play a
)");
  check(static_cast<bool>(slide), slide.diagnostic.message);
  if (slide) {
    tfseq::Runtime slideRuntime;
    slideRuntime.setProgram(slide.program.get());
    const bool sourceLegato = slideRuntime.next(0.0).events[0].legatoToNext;
    const auto target = slideRuntime.next(1.0);
    check(sourceLegato && target.count == 1 &&
              target.events[0].kind == tfseq::EventKind::Slide &&
              close(target.events[0].gateFraction, .25f),
          "a slide stays connected at onset and uses the target release");
  }

  check(!tfseq::Compile("a = sequence {\n notes x1{quiet}\n}\nplay a\n"),
        "quiet is rejected when ghost already supplies it");
  check(!tfseq::Compile("a = sequence {\n notes x1{stacc}\n}\nplay a\n"),
        "stacc is rejected when ghost already supplies it");
  check(!tfseq::Compile("a = sequence {\n notes x1{ten}\n}\nplay a\n"),
        "tenuto is rejected as contradictory to ghost");
  check(!tfseq::Compile("a = sequence {\n notes 1{stacc,ten}\n}\nplay a\n"),
        "staccato and tenuto are mutually exclusive");
  check(!tfseq::Compile("a = sequence {\n notes 1{quiet=.5}\n}\nplay a\n"),
        "articulation flags reject numerical values");
  check(!tfseq::Compile("a = sequence {\n notes 1{gate}\n}\nplay a\n"),
        "numerical event attributes still require values");
}

void zeroSlideOverridesGlideFallback() {
  const auto compiled = tfseq::Compile(R"(a = sequence {
  notes 1 >2 >3 >4{slide=0} >5{slide=40ms}
  glide 3/4
  slide . 0 . . .
}
play a
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;

  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  runtime.next(0.0);
  const auto zeroFromLane = runtime.next(1.0);
  check(zeroFromLane.count == 1 &&
            close(zeroFromLane.events[0].slideBeats, 0.f) &&
            zeroFromLane.events[0].slideMilliseconds < 0.f,
        "an explicit zero-beat Slide lane value overrides Glide");

  const auto glideFallback = runtime.next(2.0);
  check(glideFallback.count == 1 &&
            close(glideFallback.events[0].slideBeats, .75f),
        "a Slide lane default marker inherits Glide");

  const auto zeroFromAttribute = runtime.next(3.0);
  check(zeroFromAttribute.count == 1 &&
            close(zeroFromAttribute.events[0].slideBeats, 0.f) &&
            zeroFromAttribute.events[0].slideMilliseconds < 0.f,
        "an explicit zero-beat slide attribute overrides Glide");

  const auto milliseconds = runtime.next(4.0);
  check(milliseconds.count == 1 &&
            close(milliseconds.events[0].slideMilliseconds, 40.f),
        "a millisecond slide attribute remains an exact override");
}

void articulationModifiersComposeInEitherOrder() {
  const auto compiled = tfseq::Compile(R"(a = sequence {
  notes 1*3?1 1*2?1
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
  check(second.count == 2, "canonical ratchet/probability suffixes compose");
  check(!tfseq::Compile(R"(a = sequence {
  notes 1?1*2
}
play a
)"),
        "noncanonical suffix order is rejected rather than reordered");

  const auto missed = tfseq::Compile(R"(a = sequence {
  notes 1?0 2
  velocity .1 .9
}
play a
)");
  check(static_cast<bool>(missed), missed.diagnostic.message);
  if (missed) {
    tfseq::Runtime missedRuntime;
    missedRuntime.setProgram(missed.program.get());
    const auto silent = missedRuntime.next(0.0);
    const bool firstWasRest =
        silent.count == 1 && silent.events[0].kind == tfseq::EventKind::Rest;
    const auto sounded = missedRuntime.next(1.0);
    check(firstWasRest && sounded.count == 1 &&
              close(sounded.events[0].velocity, .9f),
          "a probability miss still advances pitched numerical lanes");
  }

  const auto orphanedSlide = tfseq::Compile(R"(a = sequence {
  notes 1?0 >2
}
play a
)");
  check(static_cast<bool>(orphanedSlide), orphanedSlide.diagnostic.message);
  if (orphanedSlide) {
    tfseq::Runtime slideRuntime;
    slideRuntime.setProgram(orphanedSlide.program.get());
    slideRuntime.next(0.0);
    const auto target = slideRuntime.next(1.0);
    check(target.count == 1 &&
              target.events[0].kind == tfseq::EventKind::Attack,
          "a slide after a probability miss degrades to an attack");
  }
}

void structuralPresenceShortensRealizedPasses() {
  const auto independent = tfseq::Compile(R"(a = sequence {
  notes    1 2??0 3
  velocity 0.9 0.4 0.2
  duration 0.5 1
}
play a
)");
  check(static_cast<bool>(independent), independent.diagnostic.message);
  if (independent) {
    tfseq::Runtime runtime;
    runtime.setProgram(independent.program.get());
    const auto first = runtime.next(0.0);
    const auto firstCount = first.count;
    const auto firstVelocity = first.count ? first.events[0].velocity : 0.f;
    const auto firstDuration = first.durationBeats;
    const auto third = runtime.next(0.5);
    const auto thirdCount = third.count;
    const auto thirdPitch = third.count ? third.events[0].pitchVolts : 0.f;
    const auto thirdVelocity = third.count ? third.events[0].velocity : 0.f;
    const auto thirdDuration = third.durationBeats;
    const auto wrapped = runtime.next(1.5);
    check(firstCount == 1 && close(firstVelocity, 0.9f) &&
              close(static_cast<float>(firstDuration), 0.5f) &&
              thirdCount == 1 && close(thirdPitch, 4.f / 12.f) &&
              close(thirdVelocity, 0.4f) &&
              close(static_cast<float>(thirdDuration), 1.f),
          "an omitted event consumes neither time nor independent lane values");
    check(wrapped.count == 1 && close(wrapped.events[0].pitchVolts, 0.f) &&
              close(wrapped.events[0].velocity, 0.2f),
          "independent lanes continue across the shortened Notes boundary");
  }

  const auto rightAligned = tfseq::Compile(R"(a = sequence {
  notes    1 2??0 3
  velocity ... 0.1 0.2 0.3
  cv1      ... 2 4
}
play a
)");
  check(static_cast<bool>(rightAligned), rightAligned.diagnostic.message);
  if (rightAligned) {
    tfseq::Runtime runtime;
    runtime.setProgram(rightAligned.program.get());
    const auto first = runtime.next(0.0);
    const auto firstVelocity = first.count ? first.events[0].velocity : 0.f;
    const auto firstCv = first.count ? first.events[0].cvValue[0] : 0.f;
    const auto third = runtime.next(1.0);
    check(first.count == 1 && close(firstVelocity, 0.2f) &&
              close(firstCv, 2.f) && third.count == 1 &&
              close(third.events[0].velocity, 0.3f) &&
              close(third.events[0].cvValue[0], 4.f),
          "right alignment is recalculated against the surviving Notes pass");
  }

  const auto edgeAligned = tfseq::Compile(R"(a = sequence {
  notes    1 2 3 4 ; 5 6 7 8
  velocity 0.5 ; 0.1 ... 0.1 ; 0.5
  gate     0.25!2 ... 0.5!2
  cv1      4 ; 3 ... 2 ; 1
}
play a
)");
  check(static_cast<bool>(edgeAligned), edgeAligned.diagnostic.message);
  if (edgeAligned) {
    tfseq::Runtime runtime;
    runtime.setProgram(edgeAligned.program.get());
    std::array<float, 8> velocity{};
    std::array<float, 8> gate{};
    std::array<float, 8> cv{};
    for (std::size_t step = 0; step < velocity.size(); ++step) {
      const auto events = runtime.next(static_cast<double>(step));
      if (events.count > 0) {
        velocity[step] = events.events[0].velocity;
        gate[step] = events.events[0].gateFraction;
        cv[step] = events.events[0].cvValue[0];
      }
    }
    check(close(velocity[0], .5f) && close(velocity[1], .1f) &&
              close(velocity[2], .72f) && close(velocity[5], .72f) &&
              close(velocity[6], .1f) && close(velocity[7], .5f),
          "a middle ellipsis overrides both edges and leaves the gap at the "
          "lane default");
    check(close(gate[0], .25f) && close(gate[1], .25f) &&
              close(gate[2], .8f) && close(gate[5], .8f) &&
              close(gate[6], .5f) && close(gate[7], .5f),
          "repetition counts expand independently on both sides of an "
          "ellipsis");
    check(close(cv[0], 4.f) && close(cv[1], 3.f) && close(cv[2], 0.f) &&
              close(cv[5], 0.f) && close(cv[6], 2.f) && close(cv[7], 1.f),
          "sparse edge alignment leaves uncovered CV cells at zero volts");
  }

  const auto overlappingEdges = tfseq::Compile(R"(a = sequence {
  notes    1 2??0 3??0 4
  velocity .1 .2 ... .8 .9
}
play a
)");
  check(static_cast<bool>(overlappingEdges),
        overlappingEdges.diagnostic.message);
  if (overlappingEdges) {
    tfseq::Runtime runtime;
    runtime.setProgram(overlappingEdges.program.get());
    const auto first = runtime.next(0.0);
    const auto firstVelocity = first.count ? first.events[0].velocity : 0.f;
    const auto last = runtime.next(1.0);
    check(first.count == 1 && last.count == 1 &&
              close(firstVelocity, .1f) && close(last.events[0].velocity, .9f),
          "overlapping edge groups retain their outermost values after "
          "presence omissions");
  }

  const bool rejectsTrailingSeparator = !tfseq::Compile(R"(a = sequence {
  notes 1 2 ;
}
play a
)");
  const bool rejectsSecondEllipsis = !tfseq::Compile(R"(a = sequence {
  notes 1 2
  gate 0.5 ... 0.5 ...
}
play a
)");
  check(rejectsTrailingSeparator && rejectsSecondEllipsis,
        "visual separators cannot trail and an aligned lane has one ellipsis");

  const auto alignmentBeforeEdit = tfseq::Compile(R"(a = sequence {
  notes 1 2 3
  velocity ... 0.2 0.4
}
play a
)");
  const auto alignmentAfterEdit = tfseq::Compile(R"(a = sequence {
  notes 1 2??0 3
  velocity ... 0.2 0.4
}
play a
)");
  check(static_cast<bool>(alignmentBeforeEdit) &&
            static_cast<bool>(alignmentAfterEdit),
        alignmentBeforeEdit ? alignmentAfterEdit.diagnostic.message
                            : alignmentBeforeEdit.diagnostic.message);
  if (alignmentBeforeEdit && alignmentAfterEdit) {
    tfseq::Runtime runtime;
    runtime.setProgram(alignmentBeforeEdit.program.get());
    runtime.next(0.0);
    runtime.replaceProgram(alignmentAfterEdit.program.get(), 1.0);
    const auto next = runtime.next(1.0);
    check(
        next.count == 1 && close(next.events[0].pitchVolts, 4.f / 12.f) &&
            close(next.events[0].velocity, 0.4f),
        "hot replacement recalculates realized alignment for the new program");
  }

  const auto leftAligned = tfseq::Compile(R"(a = sequence {
  notes    1 2??0 3
  velocity .1 .2!2 ...
}
play a
)");
  check(static_cast<bool>(leftAligned), leftAligned.diagnostic.message);
  if (leftAligned) {
    tfseq::Runtime runtime;
    runtime.setProgram(leftAligned.program.get());
    const auto first = runtime.next(0.0);
    const auto firstVelocity = first.count ? first.events[0].velocity : 0.f;
    const auto third = runtime.next(1.0);
    check(first.count == 1 && close(firstVelocity, 0.1f) && third.count == 1 &&
              close(third.events[0].velocity, 0.2f),
          "left alignment retains the first values of the surviving pass");
  }

  const auto omittedStructures = tfseq::Compile(R"(a = sequence {
  notes 1 [2 3]??0 ~??0 4
}
play a
)");
  check(static_cast<bool>(omittedStructures),
        omittedStructures.diagnostic.message);
  if (omittedStructures) {
    tfseq::Runtime runtime;
    runtime.setProgram(omittedStructures.program.get());
    const auto first = runtime.next(0.0);
    const auto last = runtime.next(1.0);
    check(first.count == 1 && last.count == 1 &&
              close(last.events[0].pitchVolts, 5.f / 12.f),
          "complete groups and rests can be omitted as top-level spans");
  }

  const auto slideAcrossOmission = tfseq::Compile(R"(a = sequence {
  notes 1 ~??0 >2
}
play a
)");
  check(static_cast<bool>(slideAcrossOmission),
        slideAcrossOmission.diagnostic.message);
  if (slideAcrossOmission) {
    tfseq::Runtime runtime;
    runtime.setProgram(slideAcrossOmission.program.get());
    const auto source = runtime.next(0.0);
    const bool sourceLegato =
        source.count == 1 && source.events[0].legatoToNext;
    const auto target = runtime.next(1.0);
    check(sourceLegato && target.count == 1 &&
              target.events[0].kind == tfseq::EventKind::Slide,
          "an omitted span preserves the sounding predecessor for a slide");
  }

  const auto changedVoiceCount = tfseq::Compile(R"(a = sequence {
  notes (1 3 5) [1]??0 >2
}
play a
)");
  check(static_cast<bool>(changedVoiceCount),
        changedVoiceCount.diagnostic.message);
  if (changedVoiceCount) {
    tfseq::Runtime runtime;
    runtime.setProgram(changedVoiceCount.program.get());
    const auto source = runtime.next(0.0);
    const bool sourceLegato = source.count > 0 && source.events[0].legatoToNext;
    const auto target = runtime.next(1.0);
    check(!sourceLegato && target.count == 1 &&
              target.events[0].kind == tfseq::EventKind::Attack,
          "voice-count changes neither hold the source nor emit an invalid "
          "slide");
  }

  const auto arrangement = tfseq::Compile(R"(a = sequence {
  notes 1 2??0
}
b = sequence {
  notes 5
}
song = a + b
play song
)");
  check(static_cast<bool>(arrangement), arrangement.diagnostic.message);
  if (arrangement) {
    tfseq::Runtime runtime;
    runtime.setProgram(arrangement.program.get());
    const auto first = runtime.next(0.0);
    const auto nextPart = runtime.next(1.0);
    check(first.count == 1 && nextPart.count == 1 &&
              close(nextPart.events[0].pitchVolts, 7.f / 12.f),
          "an omitted final event advances the arrangement at the earlier "
          "boundary");
  }

  const auto replicated = tfseq::Compile(R"(a = sequence {
  notes 1 2??0.5!2 3
}
play a
)");
  check(static_cast<bool>(replicated), replicated.diagnostic.message);
  if (replicated) {
    const auto &steps =
        replicated.program->semantic().sequences.front().articulation;
    check(steps.size() == 4 && steps[1].presenceProbability == 0.5f &&
              steps[1].presenceIdentity != steps[2].presenceIdentity,
          "replicated optional events receive independent stable decisions");
  }

  const auto omissionBeforeRandom = tfseq::Compile(R"(a = sequence {
  notes 1??0 $
}
seed 17
play a
)");
  const auto silenceBeforeRandom = tfseq::Compile(R"(a = sequence {
  notes 1?0 $
}
seed 17
play a
)");
  check(static_cast<bool>(omissionBeforeRandom) &&
            static_cast<bool>(silenceBeforeRandom),
        omissionBeforeRandom ? silenceBeforeRandom.diagnostic.message
                             : omissionBeforeRandom.diagnostic.message);
  if (omissionBeforeRandom && silenceBeforeRandom) {
    tfseq::Runtime omissionRuntime;
    tfseq::Runtime silenceRuntime;
    omissionRuntime.setProgram(omissionBeforeRandom.program.get());
    silenceRuntime.setProgram(silenceBeforeRandom.program.get());
    const auto omittedResult = omissionRuntime.next(0.0);
    silenceRuntime.next(0.0);
    const auto silentResult = silenceRuntime.next(1.0);
    check(omittedResult.count == 1 && silentResult.count == 1 &&
              close(omittedResult.events[0].pitchVolts,
                    silentResult.events[0].pitchVolts),
          "an earlier omission does not perturb a later written random pitch");
  }

  check(!tfseq::Compile(R"(a = sequence {
  notes 1??0.5
}
play a
)"),
        "a Notes pass cannot consist entirely of optional-presence events");
  check(!tfseq::Compile(R"(a = sequence {
  notes 1 [2 3??0.5]
}
play a
)"),
        "presence probability inside a subdividing group is rejected");
  check(!tfseq::Compile(R"(a = sequence {
  notes 1 2(3,8)??0.5
}
play a
)"),
        "presence probability combined with Euclidean expansion is rejected");
}

void nestedGroupReplicationAndProbabilityKeepPreparedIdentity() {
  const auto compiled = tfseq::Compile(R"(a = sequence {
  notes [[1 2]!3 3]
}
play a
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (compiled) {
    const auto &step =
        compiled.program->semantic().sequences[0].articulation.front();
    check(step.atoms.size() == 7,
          "nested group replication emits consecutive copies");
    check(step.atoms.size() < 6 ||
              (close(static_cast<float>(step.atoms[0].offsetFraction), 0.f) &&
               close(static_cast<float>(step.atoms[2].offsetFraction), .25f) &&
               close(static_cast<float>(step.atoms[4].offsetFraction), .5f)),
          "nested replicas retain their own event positions");
  }

  const auto probabilities = tfseq::Compile(R"(a = sequence {
  notes [[[1 2]?.5 3]?.5!2]
}
play a
)");
  check(static_cast<bool>(probabilities), probabilities.diagnostic.message);
  if (!probabilities)
    return;
  const auto &atoms =
      probabilities.program->semantic().sequences[0].articulation.front().atoms;
  check(atoms.size() == 6, "nested probability example expands six events");
  if (atoms.size() != 6)
    return;
  const auto outerFirst = atoms[2].enclosingProbabilityGates.front().group;
  const auto outerSecond = atoms[5].enclosingProbabilityGates.front().group;
  check(atoms[0].enclosingProbabilityGates.size() == 2 &&
            atoms[1].enclosingProbabilityGates.size() == 2 &&
            atoms[2].enclosingProbabilityGates.size() == 1 &&
            atoms[0].enclosingProbabilityGates.back().group == outerFirst &&
            atoms[1].enclosingProbabilityGates.back().group == outerFirst,
        "one outer group probability is shared by every member of its copy");
  check(outerFirst != outerSecond &&
            atoms[0].enclosingProbabilityGates.front().group !=
                atoms[3].enclosingProbabilityGates.front().group,
        "replicated nested groups receive independent prepared decisions");
}

void transformByCycleAndArrange() {
  const std::string source = R"(a = sequence {
  subdiv 8
  notes 1 2 |> every 2 rev
}
b = sequence {
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
  const auto first = runtime.next(0);
  check(close(first.events[0].pitchVolts, 0.f) &&
            close(static_cast<float>(first.durationBeats), .5f),
        "subdiv sets the default note value without setting pass length");
  check(close(runtime.next(1).events[0].pitchVolts, 2.f / 12.f),
        "first cycle continues in source order");
  check(close(runtime.next(2).events[0].pitchVolts, 2.f / 12.f),
        "every 2 reverses the second cycle");
  check(close(runtime.next(3).events[0].pitchVolts, 0.f),
        "the complete Notes pass defines the second cycle boundary");
  check(close(runtime.next(4).events[0].pitchVolts, 7.f / 12.f),
        "arrangement advances to the stitched section");
}

void rejectInvalidInput() {
  auto compiled =
      tfseq::Compile("bad = sequence {\n notes 1 nope\n}\nplay bad\n");
  check(!compiled, "invalid degree is rejected");
  check(compiled.diagnostic.line == 2 && compiled.diagnostic.column > 1,
        "diagnostic includes a useful source location");
  check(!tfseq::Compile(R"(bad = sequence {
  notes 1
  accent .88
}
play bad
)"),
        "accent is note articulation, not a numerical lane");
}

void deterministicRandomPitchAndScalarValues() {
  const std::string source = R"(randoms = sequence {
  tonic D@4
  scale minor_pentatonic
  notes $ ${2,4} $n{3,.75} $c{0,0} $cn{6,1.5}
  velocity $u{.2,.4}
  gate $n{.6,.15}
  cv1 $u{-2,2}
}
seed 1729
play randoms
)";
  const auto compiled = tfseq::Compile(source);
  const auto replayCompiled = tfseq::Compile(source);
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  check(static_cast<bool>(replayCompiled), replayCompiled.diagnostic.message);
  if (!compiled || !replayCompiled)
    return;

  const auto &sequence = compiled.program->semantic().sequences[0];
  check(sequence.notes.size() == 5 && sequence.notes[0].randomDefaultRange &&
            sequence.notes[1].randomDomain ==
                tfseq::PitchItem::RandomDomain::ScaleDegree &&
            sequence.notes[3].randomDomain ==
                tfseq::PitchItem::RandomDomain::ChromaticSemitone,
        "random pitch forms lower to typed semantic values");
  check(sequence.velocity[0].randomDistribution ==
                tfseq::ScalarItem::RandomDistribution::Uniform &&
            sequence.gate[0].randomDistribution ==
                tfseq::ScalarItem::RandomDistribution::Normal,
        "random scalar forms retain their distributions");

  tfseq::Runtime first;
  tfseq::Runtime replay;
  first.setProgram(compiled.program.get());
  replay.setProgram(replayCompiled.program.get());
  constexpr int pentatonicSemitones[] = {0, 3, 5, 7, 10};
  for (int step = 0; step < 100; ++step) {
    const auto one = first.next(static_cast<double>(step));
    const auto two = replay.next(static_cast<double>(step));
    check(one.count == 1 && two.count == 1,
          "each random pitch remains one prepared voice");
    if (one.count != 1 || two.count != 1)
      continue;
    check(close(one.events[0].pitchVolts, two.events[0].pitchVolts) &&
              close(one.events[0].velocity, two.events[0].velocity) &&
              close(one.events[0].gateFraction, two.events[0].gateFraction) &&
              close(one.events[0].cvValue[0], two.events[0].cvValue[0]),
          "the seed and event path reproduce every random draw");
    check(one.events[0].velocity >= .2f && one.events[0].velocity <= .4f &&
              one.events[0].gateFraction >= 0.f &&
              one.events[0].gateFraction <= 1.f &&
              one.events[0].cvValue[0] >= -2.f &&
              one.events[0].cvValue[0] <= 2.f,
          "random scalar draws remain within lane domains");

    const int patternPosition = step % 5;
    const int semitone =
        static_cast<int>(std::lround(one.events[0].pitchVolts * 12.f)) - 2;
    if (patternPosition == 0 || patternPosition == 2) {
      const int pitchClass = ((semitone % 12) + 12) % 12;
      check(std::find(std::begin(pentatonicSemitones),
                      std::end(pentatonicSemitones),
                      pitchClass) != std::end(pentatonicSemitones),
            "scale-domain random notes are quantized to the active scale");
    } else if (patternPosition == 1) {
      check(semitone == 3 || semitone == 5 || semitone == 7,
            "uniform scale-degree bounds are inclusive and quantized");
    } else if (patternPosition == 3) {
      check(semitone == 0,
            "chromatic uniform zero is the unquantized tonic offset");
    } else {
      check(close(one.events[0].pitchVolts * 12.f,
                  std::round(one.events[0].pitchVolts * 12.f)),
            "chromatic normal pitches round to semitone offsets");
    }
  }

  const auto milliseconds = tfseq::Compile(R"(a = sequence {
  notes 1
  gate $u{10ms,20ms}
}
play a
)");
  check(static_cast<bool>(milliseconds), milliseconds.diagnostic.message);
  if (milliseconds) {
    tfseq::Runtime runtime;
    runtime.setProgram(milliseconds.program.get());
    const auto event = runtime.next(0.0);
    check(event.count == 1 && event.events[0].gateMilliseconds >= 10.f &&
              event.events[0].gateMilliseconds <= 20.f,
          "uniform scalar arguments preserve matching millisecond units");
  }

  check(!tfseq::Compile("a = sequence {\n notes $u{1.5,4}\n}\nplay a\n"),
        "uniform pitch ranges require integer endpoints");
  check(!tfseq::Compile("a = sequence {\n notes $n{3,0}\n}\nplay a\n"),
        "normal pitch standard deviation must be positive");
  check(!tfseq::Compile(
            "a = sequence {\n notes 1\n velocity $u{.8,.2}\n}\nplay a\n"),
        "uniform scalar bounds must be ordered");
  check(
      !tfseq::Compile("a = sequence {\n notes 1\n cv1 $u{1ms,2}\n}\nplay a\n"),
      "random scalar units must match and suit the lane");
}

void randomPreparationUsesBoundedRanges() {
  const auto bounded = tfseq::Compile(R"(a = sequence {
  notes 1
  duration $u{1/4,1/4}
  ratchet $u{2,2}
  offset $u{-2,-2}
}
play a
)");
  check(static_cast<bool>(bounded), bounded.diagnostic.message);
  if (bounded) {
    check(bounded.program->maximumEventsPerStep >= 2,
          "uniform random ratchets contribute their upper bound to capacity");
    check(close(static_cast<float>(bounded.program->maximumEarlyBeats), 2.f),
          "uniform random offsets contribute their early bound to lookahead");
    check(bounded.program->scheduleCapacity >= 28,
          "uniform random durations contribute their lower bound to capacity");
    tfseq::Runtime runtime;
    runtime.setProgram(bounded.program.get());
    const auto step = runtime.next(0.0);
    check(step.count == 2 && !step.overflowed &&
              close(static_cast<float>(step.durationBeats), .25f) &&
              close(static_cast<float>(step.events[0].beat), -2.f),
          "bounded random timing and ratchets fit their prepared workspace");
  }

  const auto normal = tfseq::Compile(R"(a = sequence {
  notes 1
  ratchet $n{3,.6}
  offset $n{0,.5}
}
play a
)");
  check(static_cast<bool>(normal), normal.diagnostic.message);
  if (normal) {
    check(normal.program->maximumEventsPerStep >= 5 &&
              close(static_cast<float>(normal.program->maximumEarlyBeats), 2.f),
          "four-sigma normal bounds prepare ratchet and timing capacity");
    tfseq::Runtime runtime;
    runtime.setProgram(normal.program.get());
    for (int stepIndex = 0; stepIndex < 64; ++stepIndex) {
      const auto step = runtime.next(static_cast<double>(stepIndex));
      check(step.count >= 1 && step.count <= 5 && !step.overflowed,
            "normal random ratchets stay inside their prepared bound");
    }
  }
}

void selectiveEvaluationPreservesRandomIdentity() {
  const std::string evaluated = R"(a = sequence {
  notes 1
}
b = sequence {
  notes ${1,97}
}
play b
)";
  const std::string draft = R"(a = sequence {
  notes 2
}
b = sequence {
  notes ${2,98}
}
play b
)";
  const auto evaluatedDocument = tfseq::syntax::Parse(evaluated);
  const auto draftDocument = tfseq::syntax::Parse(draft);
  check(evaluatedDocument && draftDocument,
        "random selective-evaluation fixtures parse");
  if (!evaluatedDocument || !draftDocument)
    return;
  const auto selectionBegin = static_cast<int>(draft.find("notes 2"));
  const auto merged = tfseq::syntax::MergeSelectionDocuments(
      evaluatedDocument.document, evaluated, draftDocument.document, draft,
      selectionBegin, selectionBegin + 7);
  const auto original = tfseq::Compile(evaluatedDocument.document);
  const auto selected =
      merged ? tfseq::Compile(merged.document) : tfseq::Compile("");
  check(merged && original && selected,
        "random selective-evaluation programs compile");
  if (!merged || !original || !selected)
    return;

  tfseq::Runtime before;
  tfseq::Runtime after;
  before.setProgram(original.program.get());
  after.setProgram(selected.program.get());
  for (int stepIndex = 0; stepIndex < 64; ++stepIndex) {
    const auto expected = before.next(static_cast<double>(stepIndex));
    const auto actual = after.next(static_cast<double>(stepIndex));
    check(
        expected.count == 1 && actual.count == 1 &&
            close(expected.events[0].pitchVolts, actual.events[0].pitchVolts),
        "an inactive random edit does not change the evaluated random stream");
  }
}

void documentedMusicalExamplesCompile() {
  const auto rhythmForms = tfseq::Compile(R"(forms = sequence {
  notes 1 _!3 ~_ ~{len=2} ~!3 1_ 1__ 1_3 1. 1.. 1_. 3?0.35 [1|3|5]?0.5 $?0.5
}
play forms
)");
  check(static_cast<bool>(rhythmForms),
        "documented rhythm forms: " + rhythmForms.diagnostic.message);

  const auto quickStart = tfseq::Compile(R"(riff = sequence {
  subdiv 8
  tonic D@3
  scale dorian
  notes 1 2 ^3{stacc} 4 5{quiet} >6 7{ten} _
}
answer = riff |> shift_degree 3 |> octave 1
song = riff * 2 + answer
seed 42
play song
)");
  check(static_cast<bool>(quickStart),
        "quick-start reference example: " + quickStart.diagnostic.message);

  const auto readme = tfseq::Compile(R"(riff = sequence {
  subdiv 16
  tonic D@3
  scale dorian
  glide 1/8
  notes ^1 1!2 x1 [3 4] >5{stacc} ~ 6{quiet} 7{ten}
}
fill = sequence {
  subdiv 16
  tonic D@3
  scale harmonic_minor
  notes [5 6 7] ^8*2 ~ V7
}
song = riff * 3 + fill
seed 42
play song
)");
  check(static_cast<bool>(readme),
        "README reference example: " + readme.diagnostic.message);

  const auto bass = tfseq::Compile(R"(bass = sequence {
  subdiv 16
  tonic E@2
  scale minor
  glide 1/8
  notes ^1 1!2 x1 [5 6] >b7{stacc} ~ 1{ten} 3?0.35
}
fill = bass |> rotate 2 |> every 2 rev
song = bass * 3 + fill
seed 73
play song
)");
  check(static_cast<bool>(bass),
        "bass reference example: " + bass.diagnostic.message);

  const auto sections = tfseq::Compile(R"(verse = sequence {
  subdiv 8
  tonic D@3
  scale minor
  notes 1 3 4 ^5{stacc} 1' 7 5 4
}
chorus = sequence {
  subdiv 4
  tonic Bb@3
  scale major
  notes I_2{ten} V_2{ten} vi_2{ten} IV_2{ten}
}
fill = sequence {
  subdiv 16
  tonic D@3
  scale harmonic_minor
  notes [5 6 7] ^8*2 ~ V7
}
song = verse * 2 + chorus + verse + fill + chorus * 2
play song
)");
  check(static_cast<bool>(sections),
        "multi-section reference example: " + sections.diagnostic.message);

  const auto generative = tfseq::Compile(R"(melody = sequence {
  subdiv 16
  tonic A@3
  scale minor_pentatonic
  notes $u{1,10}(5,8) $n{5,1.25}(3,8,2)
}
|> every 4 (rotate 1)
seed 2026
play melody
)");
  check(static_cast<bool>(generative),
        "generative reference example: " + generative.diagnostic.message);

  const auto cv = tfseq::Compile(R"(texture = sequence {
  subdiv 8
  tonic D@3
  scale dorian
  notes 1{ten} [3 5] 7{quiet} <8 6>
  cv1 0 . 5 . |> interp smooth
}
play texture
)");
  check(static_cast<bool>(cv),
        "CV reference example: " + cv.diagnostic.message);

  const auto groove = tfseq::Compile(R"(groove = sequence {
  notes [1 5] 3 [4 6 8] 5
  offset -7ms 4ms 0 |> rate 1/2
}
|> swing .58 1/8
|> early random 3ms
seed 99
play groove
)");
  check(static_cast<bool>(groove),
        "timing reference example: " + groove.diagnostic.message);
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
  check(!tfseq::Compile(R"(ambiguous = sequence {
  scale pentatonic
  notes 1
}
play ambiguous
)"),
        "bare pentatonic is rejected instead of silently selecting major");

  auto compiled = tfseq::Compile(R"(p = sequence {
  tonic C@4
  scale major_pentatonic
  notes 1 6 1'
}

o = sequence {
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

void octaveLaneUsesAbsoluteValuesAndTonicFallback() {
  const auto compiled = tfseq::Compile(R"(octaves = sequence {
  tonic D@4
  scale major
  notes 1 1 1
  octave . 3 .
}
play octaves
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;

  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  check(close(runtime.next(0.0).events[0].pitchVolts, 2.f / 12.f),
        "an Octave lane default marker inherits the tonic octave");
  check(close(runtime.next(1.0).events[0].pitchVolts, -10.f / 12.f),
        "an Octave lane number selects an absolute octave");
  check(close(runtime.next(2.0).events[0].pitchVolts, 2.f / 12.f),
        "the inherited tonic octave remains cyclically available");
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
  notes Cmaj/E _
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

void romanChordsRetainDegreeAndQualitySemantics() {
  const auto compiled = tfseq::Compile(R"(harmony = sequence {
  tonic C@4
  scale major
  notes I i iim7 bVII
}
play harmony
)");
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const auto major = runtime.next(0.0);
  check(major.count == 3 && close(major.events[0].pitchVolts, 0.f) &&
            close(major.events[1].pitchVolts, 4.f / 12.f),
        "uppercase Roman degrees imply a major triad");
  const auto minor = runtime.next(1.0);
  check(minor.count == 3 && close(minor.events[1].pitchVolts, 3.f / 12.f),
        "lowercase Roman degrees imply a minor triad");
  const auto supertonic = runtime.next(2.0);
  check(supertonic.count == 4 &&
            close(supertonic.events[0].pitchVolts, 2.f / 12.f) &&
            close(supertonic.events[3].pitchVolts, 1.f),
        "Roman chord extensions preserve the scale-relative root");
  const auto flatSeven = runtime.next(3.0);
  check(flatSeven.count == 3 &&
            close(flatSeven.events[0].pitchVolts, 10.f / 12.f),
        "accidentals apply to a Roman chord root");
  check(!tfseq::Compile(R"(a = sequence {
  notes Im7
}
play a
 )"),
        "contradictory Roman case and quality are rejected");
  check(!tfseq::Compile(R"(a = sequence {
  scale major_pentatonic
  notes VI
}
play a
 )"),
        "Roman degree cannot exceed the active scale cardinality");
  check(!tfseq::Compile(R"(a = sequence {
  scale major
  notes VII
}
reduced = a |> scale major_pentatonic
play reduced
 )"),
        "derived scale changes revalidate Roman cardinality");
}

void scaleAndModulationPipelinesComposeLeftToRight() {
  const auto compiled = tfseq::Compile(R"(base = sequence {
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
  subdiv 8 |> fast 2
  notes 1
}
play a
)");
  check(!cycle,
        "a pipeline attached to subdiv is rejected instead of discarded");
  const auto scale = tfseq::Compile(R"(a = sequence {
  scale minor |> rev
  notes 1
}
play a
)");
  check(!scale, "setting diagnostics explain where sequence transforms belong");
}

void harmonicMinorBuildsMajorDominantInKey() {
  auto compiled = tfseq::Compile(R"(tonic_triad = sequence {
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

void conciseAcidSyntax() {
  const std::string source = R"(acid = sequence {
  tonic D#@2
  scale minor
  notes ^1 1!3 ^5 7 ^1 ~ 1 ~ ^8 >1, ~ 1 >^1, >1
  velocity .5
  gate .5
  glide .8
}
iv = acid |> modulate_degree 4 |> octave -1
v = acid |> modulate_degree 5 |> octave -1 |> scale major
song = acid * 8 + iv * 4 + v * 4
play song
)";
  auto compiled = tfseq::Compile(source);
  check(static_cast<bool>(compiled), compiled.diagnostic.message);
  if (!compiled)
    return;
  check(compiled.program->semantic().sequences.size() == 3,
        "derived acid sections reuse one concise sequence");
  check(compiled.program->semantic().arrangement.size() == 3,
        "acid sections form a three-part song");
  check(compiled.program->semantic().sequences[0].notes.size() == 13,
        "note repetition expands only the pitched events");
  std::size_t accented = 0;
  for (const auto &step :
       compiled.program->semantic().sequences[0].articulation) {
    for (const auto &atom : step.atoms)
      accented += atom.hasAccent ? 1u : 0u;
  }
  check(accented == 5, "acid accents are attached to their note events");

  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const auto first = runtime.next(0);
  check(first.count == 1 && close(first.events[0].pitchVolts, -21.f / 12.f),
        "acid starts on D-sharp 2");
  check(close(first.events[0].velocity, .88f), "first acid step is accented");
  check(close(first.events[0].accent, .88f),
        "note articulation drives the dedicated accent output");
  for (int beat = 1; beat < 7; ++beat)
    runtime.next(beat);
  check(runtime.next(7).events[0].kind == tfseq::EventKind::Rest,
        "compact replicated articulation retains the first rest");
  const auto slide = runtime.next(8);
  check(slide.events[0].kind == tfseq::EventKind::Attack &&
            close(slide.events[0].slideBeats, .8f),
        "the first event after a rest retriggers while retaining glide setup");

  runtime.reset();
  const auto songLine = source.find("song =");
  const int songAcid = static_cast<int>(source.find("acid", songLine));
  const int songIv = static_cast<int>(source.find("iv", songLine));
  const int songV = static_cast<int>(source.find("v", songIv + 2));
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
    const int expectedCursor = beat < 128   ? songAcid
                               : beat < 192 ? songIv
                                            : songV;
    check(sequenceCursor == expectedCursor,
          "the beat cursor stays on the active term of the song row");
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
}

void hotSwapPreservesNamedSequencePhase() {
  auto original = tfseq::Compile(R"(riff = sequence {
  notes 1 2 3 4
}
play riff
)");
  auto edited = tfseq::Compile(R"(riff = sequence {
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

void hotSwapPreservesActiveArrangementTermIdentity() {
  auto original = tfseq::Compile(R"(a = sequence {
  notes 1
}
b = sequence {
  notes 5
}
c = sequence {
  notes 7
}
song = a + b + c
play song
)");
  auto edited = tfseq::Compile(R"(a = sequence {
  notes 1
}
b = sequence {
  notes 5
}
c = sequence {
  notes 7
}
d = sequence {
  notes 2
}
song = b + d + a + c
play song
)");
  check(static_cast<bool>(original), original.diagnostic.message);
  check(static_cast<bool>(edited), edited.diagnostic.message);
  if (!original || !edited)
    return;

  tfseq::Runtime runtime;
  runtime.setProgram(original.program.get());
  runtime.next(0.0); // Complete a; b is now the active arrangement term.
  runtime.replaceProgram(edited.program.get(), 1.0);
  const auto next = runtime.next(1.0);
  check(next.count == 1 && close(next.events[0].pitchVolts, 7.f / 12.f),
        "hot swap follows the active named arrangement term across reordering");
}

void preparedWorkspaceHasNoSmallEventCeiling() {
  std::string subdivisions = "[";
  for (int index = 0; index < 40; ++index)
    subdivisions += index == 0 ? "1" : " 1";
  subdivisions += "]";
  const auto compiled = tfseq::Compile("dense = sequence {\n  notes " +
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

  const auto combinedDensity = tfseq::Compile(R"(dense = sequence {
  notes 1{len=1/10}
  duration 1/10
}
|> slow 1/2
|> early 1
play dense
)");
  check(static_cast<bool>(combinedDensity),
        "combined density: " + combinedDensity.diagnostic.message);
  if (combinedDensity) {
    check(combinedDensity.program->scheduleCapacity >= 400,
          "prepared capacity uses combined lane, event, and time factors");
  }
}

void playbackStateHasNoFixedSequenceCeiling() {
  std::string source;
  for (int index = 0; index < 40; ++index) {
    source += "s" + std::to_string(index) + " = sequence {\n";
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
  notes 1 _ >3
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
  notes 1
}
b = sequence {
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
  velocity .5 |> transpose_semitone 2
}
play a
)");
  check(!wrongLane,
        "a pitch transform on a scalar lane is a diagnostic, not a no-op");
}

void sequenceConditionsStayCoherentAcrossLanes() {
  const auto compiled = tfseq::Compile(R"(coherent = sequence {
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
  notes [1 2]
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
  offset -1/2 -4ms
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
  check(close(static_cast<float>(timed.program->maximumEarlyBeats), 1.5f) &&
            close(static_cast<float>(timed.program->maximumEarlyMilliseconds),
                  12.f),
        "whole and lane timing bounds combine for prepared scheduling");
  const auto lookahead =
      tfseq::SchedulingLookaheadBeats(*timed.program, true, 24000.0, 48000.0);
  check(close(static_cast<float>(lookahead), 2.524f),
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
  notes Bbm7b5@3 / D@2 (1 b3 5)@4 (C E G)@3 / B@3
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
  check(steps[0].kind == tfseq::syntax::PatternKind::Event &&
            steps[0].children[0].kind == tfseq::syntax::PatternKind::Slash &&
            steps[0].children[0].children[0].atom.text == "Bbm7b5@3" &&
            steps[0].children[0].children[1].atom.text == "D@2",
        "PEG emits the jazz chord and slash bass relationship directly");
  check(steps[1].kind == tfseq::syntax::PatternKind::Event &&
            steps[1].children[0].kind == tfseq::syntax::PatternKind::Voicing &&
            steps[1].children[0].children.size() == 3 &&
            steps[1].children[0].suffix.text == "@4",
        "PEG emits explicit tones and their shared register directly");
  check(steps[2].kind == tfseq::syntax::PatternKind::Event &&
            steps[2].children[0].kind == tfseq::syntax::PatternKind::Slash &&
            steps[2].children[0].children[0].kind ==
                tfseq::syntax::PatternKind::Voicing,
        "a slash after an explicit voicing remains structural AST");
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
  notes 1(3,65)*17
  glide 65
}
play a
)");
  check(static_cast<bool>(repeated), repeated.diagnostic.message);
  if (!repeated)
    return;
  check(repeated.program->semantic().sequences[0].notes.size() == 3 &&
            repeated.program->semantic().sequences[0].articulation.size() ==
                65 &&
            repeated.program->maximumEventsPerStep >= 17,
        "former fixed expansion ceilings are derived from allocation");
}

void euclideanTimingAndCvLanesHaveScoreTimeSemantics() {
  const auto euclidean = tfseq::Compile(R"(a = sequence {
  notes 1(3,8,1)
}
play a
)");
  check(static_cast<bool>(euclidean), euclidean.diagnostic.message);
  if (euclidean) {
    const auto &steps = euclidean.program->semantic().sequences[0].articulation;
    check(steps.size() == 8, "Euclidean syntax creates its stated cell count");
    std::string rhythm;
    for (const auto &step : steps)
      rhythm += step.atoms.front().kind == tfseq::ArticulationKind::Attack
                    ? 'x'
                    : '.';
    check(rhythm == ".x..x..x",
          "positive Euclidean rotation moves hits later/right");
  }

  const auto timed = tfseq::Compile(R"(a = sequence {
  notes [1 2] 3 4
  offset -1/8!2 1/8 |> rate 1/2
  cv1 0 . . 6 |> interp linear
  cv2 0 4 |> interp power 2
  cv3 -1 1 |> interp smooth
}
play a
)");
  check(static_cast<bool>(timed), timed.diagnostic.message);
  if (timed) {
    tfseq::Runtime runtime;
    runtime.setProgram(timed.program.get());
    const auto first = runtime.next(0.0);
    check(first.count == 2 &&
              close(static_cast<float>(first.events[0].beat), -.125f) &&
              close(static_cast<float>(first.events[1].beat), .375f),
          "subdivisions sample one score-time offset without accelerating it");
    check(
        close(first.events[0].cvValue[0], 0.f) &&
            close(first.events[0].cvTarget[0], 6.f) &&
            close(static_cast<float>(first.events[0].cvTargetBeat[0]), 2.875f),
        "linear CV skips no-op cells and targets the next explicit knot");

    const auto halfBeat = runtime.next(.5);
    check(halfBeat.count == 1 && close(halfBeat.events[0].cvValue[0], 1.f) &&
              close(halfBeat.events[0].cvValue[1], 1.f) &&
              close(halfBeat.events[0].cvValue[2], 0.f),
          "all three CV lanes interpolate within Notes-pass score time");
    runtime.next(2.0);
    runtime.next(3.0);
    const auto later = runtime.next(4.0);
    check(later.count > 0 &&
              close(static_cast<float>(later.events[0].beat), 4.125f),
          "a rate-controlled lane preserves phase across Notes passes");
  }

  const auto polyrhythmic = tfseq::Compile(R"(arpeggio = sequence {
  subdiv 16
  tonic D@4
  scale dorian
  notes 3' 1' 5 3 ; 1' 5 3 1 ; 2' 7 5 2 ; 7 5 2 1
  gate .52
  cv1 4 0 |> interp linear |> rate 1/5
}
play arpeggio
)");
  check(static_cast<bool>(polyrhythmic), polyrhythmic.diagnostic.message);
  if (polyrhythmic) {
    tfseq::Runtime runtime;
    runtime.setProgram(polyrhythmic.program.get());
    const auto first = runtime.next(0.0);
    const auto firstCv = first.events[0].cvValue[0];
    for (int step = 1; step < 16; ++step)
      runtime.next(step * .25);
    const auto secondPass = runtime.next(4.0);
    const auto secondPassCv = secondPass.events[0].cvValue[0];
    for (int step = 17; step < 20; ++step)
      runtime.next(step * .25);
    const auto low = runtime.next(5.0);
    check(first.count == 1 && secondPass.count == 1 && low.count == 1,
          "a semicolon-grouped 16-step phrase keeps one event per sixteenth");
    check(close(firstCv, 4.f) && close(secondPassCv, .8f) &&
              close(low.events[0].cvValue[0], 0.f),
          "a ten-beat CV triangle continues through a four-beat Notes pass");
  }

  check(!tfseq::Compile(R"(a = sequence {
  notes 1 2
  offset ... -1/8 |> rate 1/2
}
play a
 )"),
        "rate and edge alignment are distinct lane modes");
  check(!tfseq::Compile(R"(a = sequence {
  notes 1 2
  cv1 0 5 ... |> interp linear
}
play a
 )"),
        "continuous aligned CV waits for exact structural knot timing");
  check(!tfseq::Compile(R"(a = sequence {
  notes 1
  cv4 1
}
play a
 )"),
        "unavailable CV output numbers are semantic diagnostics");

  const auto cyclic = tfseq::Compile(R"(a = sequence {
  notes 1
  cv1 . . 5
}
play a
)");
  check(static_cast<bool>(cyclic), cyclic.diagnostic.message);
  if (cyclic) {
    tfseq::Runtime runtime;
    runtime.setProgram(cyclic.program.get());
    const auto event = runtime.next(0.0);
    check(event.count == 1 && close(event.events[0].cvValue[0], 5.f),
          "leading CV no-ops look backward through the cyclic lane");
  }
}

void unsafeNumericInputsAreDiagnostics() {
  check(!tfseq::Compile("a = sequence {\n subdiv 0\n notes 1\n}\nplay a\n"),
        "zero subdivision is rejected");
  check(!tfseq::Compile("a = sequence {\n subdiv 7/2\n notes 1\n}\nplay a\n"),
        "fractional subdivision denominators are rejected");
  check(!tfseq::Compile("a = sequence {\n cycle 8\n notes 1\n}\nplay a\n"),
        "the fixed wall-clock cycle setting is not part of the language");
  check(!tfseq::Compile("seed .5\n"),
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
  editorShortcutTextOperationsAreStructural();
  heatmapMapsScalarIntensity();
  cursorAnimationIsFrameIndependentAndTempoBounded();
  pegFrontendBuildsTypedSyntax();
  pegFrontendRejectsMalformedStructure();
  typedPatternTreeOwnsReusableStructure();
  pegFrontendOwnsVoicingAndSlashStructure();
  lineCommentsCanTruncatePatterns();
  selectionEvaluationReplacesOnlyContainingStatements();
  explicitSingleRepeatRemainsAnArrangement();
  compileAndCycleIndependentLanes();
  parseFirstClassArticulation();
  namedGateAndDynamicsArticulations();
  zeroSlideOverridesGlideFallback();
  articulationModifiersComposeInEitherOrder();
  structuralPresenceShortensRealizedPasses();
  nestedGroupReplicationAndProbabilityKeepPreparedIdentity();
  transformByCycleAndArrange();
  rejectInvalidInput();
  deterministicRandomPitchAndScalarValues();
  randomPreparationUsesBoundedRanges();
  selectiveEvaluationPreservesRandomIdentity();
  documentedMusicalExamplesCompile();
  degreesContinueAcrossOctaves();
  scaleCardinalityAndExplicitOctaves();
  absoluteAndRelativeRegistersRemainDistinct();
  octaveLaneUsesAbsoluteValuesAndTonicFallback();
  octaveSuffixesAndChordsAreUnambiguous();
  distinguishInKeyShiftsFromModulation();
  romanChordsRetainDegreeAndQualitySemantics();
  scaleAndModulationPipelinesComposeLeftToRight();
  duplicateNamesAreRejectedInEitherOrder();
  settingLanePipelinesAreNeverIgnored();
  harmonicMinorBuildsMajorDominantInKey();
  conciseAcidSyntax();
  hotSwapPreservesNamedSequencePhase();
  hotSwapPreservesActiveArrangementTermIdentity();
  hotSwapCheckpointSurvivesLookahead();
  preparedWorkspaceHasNoSmallEventCeiling();
  playbackStateHasNoFixedSequenceCeiling();
  longDurationLegatoDoesNotDependOnSchedulerLookahead();
  nestedArrangementsAndCompositeTransformsAreTyped();
  sequenceConditionsStayCoherentAcrossLanes();
  quantifiedTimingTransforms();
  swingCrossesStepBoundariesAtExplicitSubdivisions();
  patternExpansionUsesAddressableLimits();
  euclideanTimingAndCvLanesHaveScoreTimeSemantics();
  timingPreparationCoversEarlyLookaheadAndMilliseconds();
  unsafeNumericInputsAreDiagnostics();
  if (failures != 0)
    std::cerr << failures << " text sequencer test(s) failed\n";
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
