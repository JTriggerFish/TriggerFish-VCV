#include "plugin.hpp"
#include "tfseq.hpp"
#include "tfseq_editor.hpp"
#include "tfseq_envelope.hpp"
#include "tfseq_parser.hpp"
#include "tftransport.hpp"
#include "tfui_animation.hpp"
#include "tfui_colormap.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace EditorIntensity {
constexpr float Background = 0.f;
constexpr float Selection = 0.42f;
constexpr float Comment = 0.60f;
constexpr float Status = 0.76f;
constexpr float Text = 0.86f;
} // namespace EditorIntensity

NVGcolor editorColor(tfui::HeatmapPalette palette, float intensity) noexcept {
  const auto rgb = tfui::sampleHeatmap(palette, intensity);
  return nvgRGBf(rgb.red, rgb.green, rgb.blue);
}

constexpr const char *DefaultSource = R"(riff = sequence {
  subdiv 8n
  tonic C@4
  scale minor
  notes 1 x2 3{quiet} _ [>4 ^5{stacc}] 6*3 ~ 8{ten}
  offset -6ms 0 +6ms |> rate 1/2
  cv1 0 5 0 |> interp smooth
}

play riff
)";

constexpr int ArrangementCursorPulseBeats = 4;
constexpr const char *ProgSequencerDocumentationUrl =
    "https://github.com/JTriggerFish/TriggerFish-VCV/blob/master/docs/"
    "TfProgSequencer-reference.md";

constexpr int LanguageVersion = 1;

std::uint64_t packSpan(const tfseq::SourceSpan &span) noexcept {
  if (!span.valid())
    return 0;
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(span.begin + 1))
          << 32) |
         static_cast<std::uint32_t>(span.end + 1);
}

tfseq::SourceSpan unpackSpan(std::uint64_t packed) noexcept {
  if (packed == 0)
    return {};
  tfseq::SourceSpan span;
  span.begin = static_cast<int>(packed >> 32) - 1;
  span.end = static_cast<int>(packed & 0xffffffffU) - 1;
  return span;
}

int validPanelWidth(int widthHp) noexcept {
  return widthHp == 30 || widthHp == 38 ? widthHp : 22;
}

} // namespace

struct TfProgSequencer : Module {
  enum ParamIds { NUM_PARAMS };
  enum InputIds { CLOCK_INPUT, RESET_INPUT, RUN_INPUT, NUM_INPUTS };
  enum OutputIds {
    PITCH_OUTPUT,
    GATE_OUTPUT,
    TRIGGER_OUTPUT,
    VELOCITY_OUTPUT,
    ACCENT_OUTPUT,
    CV1_OUTPUT,
    CV2_OUTPUT,
    CV3_OUTPUT,
    NUM_OUTPUTS
  };
  enum LightIds { RUN_LIGHT, NUM_LIGHTS };

  enum class TransportStatus { Waiting, Playing, Paused, Stopped };

  std::string source = DefaultSource;
  std::string evaluatedSource = DefaultSource;
  tfseq::syntax::Document evaluatedDocument;
  std::string compileMessage =
      "READY - Ctrl+. runs all, Ctrl+Enter runs statement/selection";
  std::atomic<std::uint64_t> sourceRevision{1};
  std::atomic<int> restoredEditorCursor{-1};
  std::atomic<int> restoredEditorSelection{-1};
  // CompiledProgram is at least two-byte aligned, leaving the low pointer bit
  // available to carry the restart-on-activation request in the same atomic
  // publication as the program itself.
  static_assert(alignof(tfseq::CompiledProgram) >= 2);
  static constexpr std::uintptr_t PendingRestartBit = 1;
  std::atomic<std::uintptr_t> pendingProgram{0};
  std::atomic<tfseq::CompiledProgram *> retiredProgram{nullptr};
  tfseq::CompiledProgram *activeProgram = nullptr;
  tfseq::Runtime runtime;
  tfseq::Runtime activationRuntime;
  double activationCheckpointBeat = 0.0;
  double activationNextStepBeat = 0.0;
  bool activationCheckpointValid = false;

  std::array<std::atomic<std::uint64_t>,
             static_cast<std::size_t>(tfseq::CursorLane::Count)>
      cursorSpans{};
  std::array<std::atomic<std::uint64_t>,
             static_cast<std::size_t>(tfseq::CursorLane::Count)>
      cursorPulses{};
  std::atomic<std::uint64_t> executionSpan{0};
  std::atomic<std::uint64_t> executionPulse{0};
  std::atomic<TransportStatus> transportStatus{TransportStatus::Waiting};
  std::atomic<bool> workspaceOverflow{false};
  std::atomic<double> visibleBeat{0.0};
  std::atomic<std::uint64_t> activeStepSpan{0};
  std::atomic<double> activeStepBeginBeat{0.0};
  std::atomic<double> activeStepEndBeat{0.0};
  std::array<std::atomic<std::uint64_t>, tfseq::CvLaneCount>
      activeCvLineSpans{};
  std::array<std::atomic<float>, tfseq::CvLaneCount> visibleCvValues{};
  std::atomic<int> panelWidthHp{30};
  std::atomic<tfui::HeatmapPalette> editorHeatmap{tfui::HeatmapPalette::Magma};
  std::atomic<bool> editorRunEnabled{true};

  dsp::SchmittTrigger clockTrigger;
  dsp::SchmittTrigger resetTrigger;
  dsp::PulseGenerator clockLightPulse;
  struct OutputVoiceState {
    dsp::PulseGenerator triggerPulse;
    double gateOffBeat = -1.0;
    bool gateHigh = false;
    float pitch = 0.f;
    float slideFrom = 0.f;
    float slideTo = 0.f;
    double slideBeginBeat = 0.0;
    double slideEndBeat = 0.0;
    bool sliding = false;
    float velocity = 7.2f;
    float accent = 0.f;
  };

  std::array<OutputVoiceState, tfseq::MaximumPolyphony> outputVoices{};
  struct CvOutputState {
    float value = 0.f;
    float output = 0.f;
    float from = 0.f;
    float target = 0.f;
    double beginBeat = 0.0;
    double endBeat = 0.0;
    float power = 1.f;
    bool initialized = false;
    tfseq::CvInterpolation interpolation = tfseq::CvInterpolation::Step;
    tfseq::CvEnvelopeSpec envelopeSpec;
    tfseq::CvEnvelopeEngine envelope;
    float envelopePeak = 0.f;
    bool envelopeTriggerPending = false;
  };
  std::array<CvOutputState, tfseq::CvLaneCount> cvOutputs{};
  std::size_t outputVoiceCount = 1;
  std::size_t scheduledCount = 0;
  std::uint64_t nextScheduleOrder = 0;
  bool clockSeen = false;
  bool periodKnown = false;
  bool clockTimedOut = false;
  std::int64_t samplesSinceClock = 0;
  std::int64_t transportPulseCount = 0;
  std::int64_t resetIgnoreSamples = 0;
  double periodSamples = 0.0;
  double pulseBeat = 0.0;
  double programStartBeat = 0.0;
  double nextStepBeat = 0.0;
  bool wasRunning = true;
  double sampleRateHz = 44100.0;
  std::int64_t lastArrangementCursorGroup =
      std::numeric_limits<std::int64_t>::min();
  std::uint64_t lastArrangementCursorSpan = 0;

  TfProgSequencer() {
    config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
    configInput(CLOCK_INPUT, "24 PPQN master clock");
    configInput(RESET_INPUT, "Reset");
    configInput(RUN_INPUT, "Run gate (low pauses; normalled high)");
    configOutput(PITCH_OUTPUT, "Polyphonic 1 V/oct pitch");
    configOutput(GATE_OUTPUT, "Polyphonic gate");
    configOutput(TRIGGER_OUTPUT, "Polyphonic trigger");
    configOutput(VELOCITY_OUTPUT, "Polyphonic velocity");
    configOutput(ACCENT_OUTPUT, "Polyphonic accent");
    configOutput(CV1_OUTPUT, "CV 1");
    configOutput(CV2_OUTPUT, "CV 2");
    configOutput(CV3_OUTPUT, "CV 3");
    configLight(RUN_LIGHT, "Running");
    publishSource(source);
  }

  static tfseq::CompiledProgram *
  pendingProgramPointer(std::uintptr_t tagged) noexcept {
    return reinterpret_cast<tfseq::CompiledProgram *>(tagged &
                                                      ~PendingRestartBit);
  }

  static std::uintptr_t tagPendingProgram(tfseq::CompiledProgram *program,
                                          bool restart) noexcept {
    return reinterpret_cast<std::uintptr_t>(program) |
           (restart ? PendingRestartBit : 0);
  }

  ~TfProgSequencer() override {
    delete pendingProgramPointer(pendingProgram.exchange(0));
    delete retiredProgram.exchange(nullptr);
    delete activeProgram;
  }

  void reportDiagnostic(const tfseq::Diagnostic &diagnostic) {
    compileMessage = "ERROR " + std::to_string(diagnostic.line) + ":" +
                     std::to_string(diagnostic.column) + "  " +
                     diagnostic.message;
  }

  bool publishDocument(const tfseq::syntax::Document &document,
                       const std::string &editorSource,
                       bool restartOnActivation = false) {
    tfseq::CompileResult compiled;
    try {
      compiled = tfseq::Compile(document);
    } catch (const std::bad_alloc &) {
      compileMessage = "ERROR not enough memory to compile this program";
      return false;
    } catch (const std::length_error &) {
      compileMessage = "ERROR program exceeds addressable memory";
      return false;
    }
    if (!compiled) {
      reportDiagnostic(compiled.diagnostic);
      return false;
    }
    // Prepare every potentially-throwing UI-side update before the atomic
    // publication. Once exchanged, success is guaranteed and the audio thread
    // sees a program matching the retained evaluated document.
    tfseq::syntax::Document nextDocument = document;
    std::string nextSource = editorSource;
    std::string nextMessage = "QUEUED - activates on the next quarter beat";
    evaluatedDocument.statements.swap(nextDocument.statements);
    evaluatedSource.swap(nextSource);
    compileMessage.swap(nextMessage);
    auto *candidate = compiled.program.release();
    delete pendingProgramPointer(pendingProgram.exchange(
        tagPendingProgram(candidate, restartOnActivation),
        std::memory_order_acq_rel));
    return true;
  }

  bool publishSource(const std::string &text,
                     bool restartOnActivation = false) {
    try {
      const auto parsed = tfseq::syntax::Parse(text);
      if (!parsed) {
        reportDiagnostic(parsed.diagnostic);
        return false;
      }
      return publishDocument(parsed.document, text, restartOnActivation);
    } catch (const std::bad_alloc &) {
      compileMessage = "ERROR not enough memory to parse this program";
      return false;
    } catch (const std::length_error &) {
      compileMessage = "ERROR program exceeds addressable memory";
      return false;
    }
  }

  void flashExecution(int begin, int end) {
    if (end <= begin)
      return;
    executionSpan.store(packSpan({begin, end, 1, 1}),
                        std::memory_order_relaxed);
    executionPulse.fetch_add(1, std::memory_order_release);
  }

  bool publishSelection(const std::string &selection, int begin) {
    try {
      auto draft = tfseq::syntax::Parse(source);
      if (!draft) {
        draft = tfseq::syntax::ParseStatementsContaining(
            source, begin, begin + static_cast<int>(selection.size()));
        if (!draft) {
          reportDiagnostic(draft.diagnostic);
          return false;
        }
      }
      const auto contextual = tfseq::syntax::MergeSelectionDocuments(
          evaluatedDocument, evaluatedSource, draft.document, source, begin,
          begin + static_cast<int>(selection.size()));
      if (!contextual) {
        reportDiagnostic(contextual.diagnostic);
        return false;
      }
      const bool accepted = publishDocument(contextual.document, source);
      if (accepted)
        flashExecution(begin, begin + static_cast<int>(selection.size()));
      return accepted;
    } catch (const std::bad_alloc &) {
      compileMessage = "ERROR not enough memory to evaluate this selection";
      return false;
    } catch (const std::length_error &) {
      compileMessage = "ERROR selection exceeds addressable memory";
      return false;
    }
  }

  void collectRetired() {
    delete retiredProgram.exchange(nullptr, std::memory_order_acq_rel);
  }

  json_t *dataToJson() override {
    json_t *root = json_object();
    json_object_set_new(root, "source",
                        json_stringn(source.c_str(), source.size()));
    json_object_set_new(root, "languageVersion", json_integer(LanguageVersion));
    json_object_set_new(
        root, "panelWidthHp",
        json_integer(panelWidthHp.load(std::memory_order_relaxed)));
    json_object_set_new(root, "editorHeatmap",
                        json_integer(static_cast<int>(
                            editorHeatmap.load(std::memory_order_relaxed))));
    return root;
  }

  void dataFromJson(json_t *root) override {
    int storedLanguageVersion = LanguageVersion;
    if (json_t *versionJson = json_object_get(root, "languageVersion")) {
      if (json_is_integer(versionJson))
        storedLanguageVersion =
            static_cast<int>(json_integer_value(versionJson));
    }
    if (json_t *widthJson = json_object_get(root, "panelWidthHp")) {
      if (json_is_integer(widthJson))
        panelWidthHp.store(
            validPanelWidth(static_cast<int>(json_integer_value(widthJson))),
            std::memory_order_relaxed);
    }
    if (json_t *heatmapJson = json_object_get(root, "editorHeatmap")) {
      if (json_is_integer(heatmapJson)) {
        editorHeatmap.store(tfui::heatmapPaletteFromInt(static_cast<int>(
                                json_integer_value(heatmapJson))),
                            std::memory_order_relaxed);
      }
    }
    if (json_t *sourceJson = json_object_get(root, "source")) {
      if (json_is_string(sourceJson)) {
        source = json_string_value(sourceJson);
        sourceRevision.fetch_add(1, std::memory_order_release);
        // Remove only the constructor/pending candidate before evaluating
        // saved source. A compile failure preserves any already-active valid
        // program; on initial construction there is no active program, so the
        // module remains silent instead of falling back to factory notes.
        delete pendingProgramPointer(
            pendingProgram.exchange(0, std::memory_order_acq_rel));
        evaluatedSource.clear();
        evaluatedDocument.statements.clear();
        if (storedLanguageVersion == LanguageVersion) {
          publishSource(source);
        } else {
          compileMessage = "ERROR saved language version " +
                           std::to_string(storedLanguageVersion) +
                           " is not supported by this build";
        }
      }
    }
  }

  void restoreEditedSource(const std::string &text, int cursor,
                           int selection) {
    if (source == text)
      return;
    source = text;
    restoredEditorCursor.store(cursor, std::memory_order_relaxed);
    restoredEditorSelection.store(selection, std::memory_order_relaxed);
    sourceRevision.fetch_add(1, std::memory_order_release);
  }

  void resetTransport() noexcept {
    runtime.reset();
    scheduledCount = 0;
    nextScheduleOrder = 0;
    clockSeen = false;
    periodKnown = false;
    clockTimedOut = false;
    samplesSinceClock = 0;
    transportPulseCount = 0;
    periodSamples = 0.0;
    pulseBeat = 0.0;
    programStartBeat = 0.0;
    nextStepBeat = 0.0;
    outputVoiceCount = 1;
    for (auto &voice : outputVoices) {
      voice.gateOffBeat = -1.0;
      voice.gateHigh = false;
      voice.sliding = false;
      voice.triggerPulse.reset();
    }
    cvOutputs = {};
    clockLightPulse.reset();
    activeStepSpan.store(0, std::memory_order_relaxed);
    activeStepBeginBeat.store(0.0, std::memory_order_relaxed);
    activeStepEndBeat.store(0.0, std::memory_order_relaxed);
    for (auto &span : activeCvLineSpans)
      span.store(0, std::memory_order_relaxed);
    for (auto &value : visibleCvValues)
      value.store(0.f, std::memory_order_relaxed);
    activationCheckpointValid = false;
    workspaceOverflow.store(false, std::memory_order_relaxed);
    for (auto &span : cursorSpans)
      span.store(0, std::memory_order_relaxed);
    lastArrangementCursorGroup = std::numeric_limits<std::int64_t>::min();
    lastArrangementCursorSpan = 0;
    transportStatus.store(TransportStatus::Stopped, std::memory_order_relaxed);
  }

  void swapProgramAtClock(double beat) noexcept {
    if (retiredProgram.load(std::memory_order_acquire) != nullptr)
      return;
    const auto pending = pendingProgram.exchange(0, std::memory_order_acq_rel);
    auto *candidate = pendingProgramPointer(pending);
    if (!candidate)
      return;
    const bool restartOnActivation = (pending & PendingRestartBit) != 0;
    auto *previousProgram = activeProgram;
    const bool preservePhase = !restartOnActivation && previousProgram &&
                               activationCheckpointValid &&
                               std::abs(activationCheckpointBeat - beat) < 1e-7;
    activeProgram = candidate;
    if (preservePhase) {
      runtime = activationRuntime;
      runtime.replaceProgram(activeProgram, beat - programStartBeat);
      nextStepBeat = activationNextStepBeat;
    } else {
      runtime.setProgram(activeProgram);
      programStartBeat = beat;
      nextStepBeat = beat;
    }
    scheduledCount = 0;
    nextScheduleOrder = 0;
    workspaceOverflow.store(false, std::memory_order_relaxed);
    activationCheckpointValid = false;
    if (!preservePhase) {
      outputVoiceCount = 1;
      for (auto &voice : outputVoices) {
        voice.gateHigh = false;
        voice.sliding = false;
        voice.triggerPulse.reset();
      }
    }
    for (auto &span : cursorSpans)
      span.store(0, std::memory_order_relaxed);
    activeStepSpan.store(0, std::memory_order_release);
    for (auto &span : activeCvLineSpans)
      span.store(0, std::memory_order_release);
    if (previousProgram)
      retiredProgram.store(previousProgram, std::memory_order_release);
  }

  void captureActivationCheckpoint() noexcept {
    if (!activeProgram || activationCheckpointValid ||
        nextStepBeat + 1e-9 < activationCheckpointBeat)
      return;
    runtime.snapshot(activationRuntime);
    activationNextStepBeat = nextStepBeat;
    activationCheckpointValid = true;
  }

  bool enqueue(const tfseq::RuntimeEvent &sourceEvent,
               double absoluteBeat) noexcept {
    if (!activeProgram ||
        scheduledCount >= activeProgram->scheduleWorkspace.size()) {
      workspaceOverflow.store(true, std::memory_order_relaxed);
      return false;
    }
    auto &target = activeProgram->scheduleWorkspace[scheduledCount++];
    target = sourceEvent;
    target.scheduleOrder = nextScheduleOrder++;
    const double absoluteShift = absoluteBeat - sourceEvent.beat;
    target.beat = absoluteBeat;
    for (auto &targetBeat : target.cvTargetBeat)
      targetBeat += absoluteShift;
    if (periodKnown && periodSamples > 0.0 &&
        target.timingOffsetMilliseconds != 0.0) {
      const double millisecondShift = target.timingOffsetMilliseconds * 0.001 *
                                      sampleRateHz / periodSamples;
      target.beat += millisecondShift;
      for (auto &targetBeat : target.cvTargetBeat)
        targetBeat += millisecondShift;
    }
    auto &schedule = activeProgram->scheduleWorkspace;
    auto later = [](const tfseq::RuntimeEvent &left,
                    const tfseq::RuntimeEvent &right) {
      return left.beat > right.beat ||
             (left.beat == right.beat &&
              left.scheduleOrder > right.scheduleOrder);
    };
    std::size_t child = scheduledCount - 1;
    while (child > 0) {
      const std::size_t parent = (child - 1) / 2;
      if (!later(schedule[parent], schedule[child]))
        break;
      std::swap(schedule[parent], schedule[child]);
      child = parent;
    }
    return true;
  }

  void scheduleDueSteps(double horizon) noexcept {
    if (!activeProgram || activeProgram->semantic().arrangement.empty())
      return;
    std::size_t preparedSteps = 0;
    const auto maximumPreparedSteps =
        std::max<std::size_t>(1, activeProgram->scheduleCapacity);
    captureActivationCheckpoint();
    while (nextStepBeat <= horizon + 1e-9 &&
           preparedSteps++ < maximumPreparedSteps) {
      captureActivationCheckpoint();
      const double localBeat = nextStepBeat - programStartBeat;
      const auto step = runtime.next(localBeat);
      if (step.overflowed) {
        workspaceOverflow.store(true, std::memory_order_relaxed);
        break;
      }
      if (step.durationBeats <= 0.0)
        break;
      for (std::size_t index = 0; index < step.count; ++index) {
        const auto &event = step.events[index];
        if (!enqueue(event, programStartBeat + event.beat))
          return;
      }
      const double followingStep = nextStepBeat + step.durationBeats;
      if (!(followingStep > nextStepBeat)) {
        workspaceOverflow.store(true, std::memory_order_relaxed);
        return;
      }
      nextStepBeat = followingStep;
    }
    captureActivationCheckpoint();
    if (nextStepBeat <= horizon + 1e-9)
      workspaceOverflow.store(true, std::memory_order_relaxed);
  }

  void showCursors(const tfseq::RuntimeEvent &event) noexcept {
    const auto notesLane = static_cast<std::size_t>(tfseq::CursorLane::Notes);
    const auto noteSpan = packSpan(event.cursors[notesLane]);
    activeStepBeginBeat.store(event.beat - programStartBeat,
                              std::memory_order_relaxed);
    activeStepEndBeat.store(event.beat - programStartBeat + event.spanBeats,
                            std::memory_order_relaxed);
    activeStepSpan.store(noteSpan, std::memory_order_release);
    for (std::size_t index = 0; index < activeCvLineSpans.size(); ++index) {
      const auto lane = static_cast<std::size_t>(tfseq::CvCursorLane(index));
      activeCvLineSpans[index].store(packSpan(event.cursors[lane]),
                                     std::memory_order_release);
    }
    for (std::size_t lane = 0; lane < cursorSpans.size(); ++lane) {
      if (!event.cursors[lane].valid())
        continue;
      const auto packed = packSpan(event.cursors[lane]);
      cursorSpans[lane].store(packed, std::memory_order_relaxed);
      bool pulse = true;
      if (lane == static_cast<std::size_t>(tfseq::CursorLane::Sequence)) {
        const auto group = tfui::arrangementCursorGroup(
            event.beat - programStartBeat, ArrangementCursorPulseBeats);
        pulse = group != lastArrangementCursorGroup ||
                packed != lastArrangementCursorSpan;
        lastArrangementCursorGroup = group;
        lastArrangementCursorSpan = packed;
      }
      if (pulse)
        cursorPulses[lane].fetch_add(1, std::memory_order_release);
    }
  }

  void applyEvent(const tfseq::RuntimeEvent &event) noexcept {
    if (event.kind == tfseq::EventKind::Tie ||
        event.kind == tfseq::EventKind::Rest || event.voice == 0)
      showCursors(event);
    if (event.voice == 0) {
      for (std::size_t index = 0; index < cvOutputs.size(); ++index) {
        auto &state = cvOutputs[index];
        state.interpolation = event.cvInterpolation[index];
        state.power = event.cvPower[index];
        state.envelopeSpec = event.cvEnvelope[index];
        if (tfseq::CvEnvelopeTriggers(event.kind)) {
          state.envelopePeak = tfseq::CvEnvelopePeak(
              state.envelopeSpec, event.velocity, event.accent > 0.f);
          state.envelopeTriggerPending = true;
        }
        if (state.interpolation == tfseq::CvInterpolation::Step) {
          state.value = event.cvValue[index];
          state.from = state.value;
          state.target = state.value;
          state.beginBeat = event.beat;
          state.endBeat = event.beat;
        } else {
          if (!state.initialized)
            state.value = event.cvValue[index];
          state.from = state.value;
          state.target = event.cvTarget[index];
          state.beginBeat = event.beat;
          state.endBeat = std::max(event.beat, event.cvTargetBeat[index]);
        }
        state.initialized = true;
      }
    }
    if (event.kind == tfseq::EventKind::Rest) {
      for (std::size_t voice = 0; voice < outputVoiceCount; ++voice) {
        outputVoices[voice].gateHigh = false;
        outputVoices[voice].sliding = false;
        outputVoices[voice].triggerPulse.reset();
      }
      return;
    }
    if (event.kind == tfseq::EventKind::Tie) {
      for (std::size_t voice = 0; voice < outputVoiceCount; ++voice) {
        auto &state = outputVoices[voice];
        if (state.gateHigh)
          state.gateOffBeat =
              std::max(state.gateOffBeat, event.beat + event.spanBeats);
      }
      return;
    }

    const auto requestedVoices =
        std::clamp<std::size_t>(event.voiceCount, 1, tfseq::MaximumPolyphony);
    if (requestedVoices < outputVoiceCount) {
      for (std::size_t voice = requestedVoices; voice < outputVoiceCount;
           ++voice) {
        outputVoices[voice].gateHigh = false;
        outputVoices[voice].sliding = false;
        outputVoices[voice].triggerPulse.reset();
      }
    }
    outputVoiceCount = requestedVoices;
    const auto voiceIndex =
        std::min<std::size_t>(event.voice, tfseq::MaximumPolyphony - 1);
    auto &voice = outputVoices[voiceIndex];
    voice.velocity = event.velocity * 10.f;
    voice.accent = event.accent * 10.f;
    double gateDuration =
        event.gateBeats >= 0.f
            ? event.gateBeats
        : event.gateMilliseconds >= 0.f && periodKnown && periodSamples > 0.0
            ? event.gateMilliseconds * 0.001 * sampleRateHz / periodSamples
            : event.spanBeats * event.gateFraction;
    voice.gateHigh = gateDuration > 0.0;
    if (event.gateCapMilliseconds >= 0.f && periodKnown &&
        periodSamples > 0.0) {
      gateDuration = std::min(gateDuration, event.gateCapMilliseconds * 0.001 *
                                                sampleRateHz / periodSamples);
      gateDuration =
          std::min(gateDuration, event.spanBeats * event.gateFraction);
    }
    voice.gateOffBeat = event.beat + std::min(event.spanBeats, gateDuration);
    const double boundary = event.beat + event.spanBeats;
    if (event.legatoToNext)
      voice.gateOffBeat = std::max(voice.gateOffBeat, boundary);
    if (event.kind == tfseq::EventKind::Slide) {
      voice.slideFrom = voice.pitch;
      voice.slideTo = event.pitchVolts;
      voice.slideBeginBeat = event.beat;
      voice.slideEndBeat =
          event.beat + std::min<double>(event.spanBeats,
                                        event.slideMilliseconds >= 0.f &&
                                                periodKnown &&
                                                periodSamples > 0.0
                                            ? event.slideMilliseconds * 0.001 *
                                                  sampleRateHz / periodSamples
                                            : event.slideBeats);
      voice.sliding = voice.slideEndBeat > voice.slideBeginBeat;
      if (!voice.sliding)
        voice.pitch = voice.slideTo;
    } else {
      voice.pitch = event.pitchVolts;
      voice.sliding = false;
      voice.triggerPulse.trigger(0.001f);
    }
  }

  void processScheduled(double phase) noexcept {
    if (!activeProgram)
      return;
    auto &schedule = activeProgram->scheduleWorkspace;
    auto earlier = [](const tfseq::RuntimeEvent &left,
                      const tfseq::RuntimeEvent &right) {
      return left.beat < right.beat ||
             (left.beat == right.beat &&
              left.scheduleOrder < right.scheduleOrder);
    };
    while (scheduledCount > 0 && schedule[0].beat <= phase + 1e-9) {
      const auto event = schedule[0];
      --scheduledCount;
      if (scheduledCount > 0) {
        schedule[0] = schedule[scheduledCount];
        std::size_t parent = 0;
        for (;;) {
          const std::size_t left = parent * 2 + 1;
          if (left >= scheduledCount)
            break;
          const std::size_t right = left + 1;
          const std::size_t child =
              right < scheduledCount && earlier(schedule[right], schedule[left])
                  ? right
                  : left;
          if (!earlier(schedule[child], schedule[parent]))
            break;
          std::swap(schedule[parent], schedule[child]);
          parent = child;
        }
      }
      applyEvent(event);
    }
  }

  void process(const ProcessArgs &args) override {
    if (periodKnown && sampleRateHz > 0.0 && args.sampleRate > 0.0 &&
        args.sampleRate != sampleRateHz) {
      periodSamples *= args.sampleRate / sampleRateHz;
    }
    sampleRateHz = args.sampleRate;
    if (resetTrigger.process(inputs[RESET_INPUT].getVoltage(), 0.1f, 2.f)) {
      resetTransport();
      resetIgnoreSamples =
          std::max<std::int64_t>(1, std::llround(args.sampleRate * 0.001));
    }
    if (resetIgnoreSamples > 0)
      --resetIgnoreSamples;

    const bool running = editorRunEnabled.load(std::memory_order_relaxed) &&
                         (!inputs[RUN_INPUT].isConnected() ||
                          inputs[RUN_INPUT].getVoltage() >= 1.f);
    bool clockEdge =
        clockTrigger.process(inputs[CLOCK_INPUT].getVoltage(), 0.1f, 2.f);
    if (resetIgnoreSamples > 0)
      clockEdge = false;

    if (!running) {
      if (wasRunning) {
        clockLightPulse.reset();
        for (auto &voice : outputVoices)
          voice.triggerPulse.reset();
      }
      if (transportStatus.load(std::memory_order_relaxed) !=
          TransportStatus::Stopped) {
        transportStatus.store(TransportStatus::Paused,
                              std::memory_order_relaxed);
      }
    } else if (clockEdge) {
      bool beatBoundary = false;
      if (clockSeen && samplesSinceClock > 0) {
        const double measured =
            tftransport::BeatPeriodFromPulseSamples(samplesSinceClock);
        // A stopped external clock produces one very long apparent interval.
        // Retain the last valid period on the returning edge; folding the gap
        // into the smoother would make subdivisions crawl for several beats.
        if (!clockTimedOut) {
          periodSamples =
              periodKnown ? 0.75 * periodSamples + 0.25 * measured : measured;
        }
        periodKnown = true;
        clockTimedOut = false;
        ++transportPulseCount;
        pulseBeat = tftransport::BeatAtPulse(transportPulseCount);
        beatBoundary = tftransport::IsQuarterNotePulse(transportPulseCount);
      } else {
        clockSeen = true;
        clockTimedOut = false;
        transportPulseCount = 0;
        pulseBeat = 0.0;
        beatBoundary = true;
      }
      samplesSinceClock = 0;
      if (beatBoundary) {
        clockLightPulse.trigger(0.075f);
        swapProgramAtClock(pulseBeat);
        activationCheckpointBeat = pulseBeat + 1.0;
        activationCheckpointValid = false;
      }
      // A subdivision exactly on this clock, and any first-beat subdivisions
      // held while the clock period was still unknown, already occupy the
      // prepared queue. Apply them before filling the next lookahead window so
      // two adjacent windows never need to coexist transiently.
      processScheduled(pulseBeat);
      if (activeProgram)
        scheduleDueSteps(pulseBeat + tfseq::SchedulingLookaheadBeats(
                                         *activeProgram, periodKnown,
                                         periodSamples, sampleRateHz));
      processScheduled(pulseBeat);
      transportStatus.store(activeProgram ? TransportStatus::Playing
                                           : TransportStatus::Waiting,
                            std::memory_order_relaxed);
    } else if (!wasRunning) {
      transportStatus.store(activeProgram ? TransportStatus::Playing
                                           : TransportStatus::Waiting,
                            std::memory_order_relaxed);
    }

    double phase = pulseBeat;
    if (clockSeen && periodKnown && periodSamples > 0.0) {
      // Every transport pulse owns its exact phase point. Interpolation covers
      // only the interval before the following pulse arrives.
      const double withinBeat =
          std::min(tftransport::BeatsPerPulse - 1e-9,
                   static_cast<double>(samplesSinceClock) / periodSamples);
      phase += withinBeat;
    }
    if (running && clockSeen) {
      if (activeProgram)
        scheduleDueSteps(pulseBeat + tfseq::SchedulingLookaheadBeats(
                                         *activeProgram, periodKnown,
                                         periodSamples, sampleRateHz));
      processScheduled(phase);
      if (periodKnown && samplesSinceClock > periodSamples * 4.0) {
        if (!clockTimedOut) {
          for (auto &voice : outputVoices) {
            voice.gateHigh = false;
            voice.sliding = false;
          }
        }
        clockTimedOut = true;
        transportStatus.store(TransportStatus::Waiting,
                              std::memory_order_relaxed);
      }
    }

    for (std::size_t index = 0; index < outputVoiceCount; ++index) {
      auto &voice = outputVoices[index];
      if (voice.sliding) {
        const double denominator = voice.slideEndBeat - voice.slideBeginBeat;
        const float amount =
            denominator > 0.0
                ? static_cast<float>((phase - voice.slideBeginBeat) /
                                     denominator)
                : 1.f;
        voice.pitch = voice.slideFrom + (voice.slideTo - voice.slideFrom) *
                                            std::clamp(amount, 0.f, 1.f);
        if (phase >= voice.slideEndBeat) {
          voice.pitch = voice.slideTo;
          voice.sliding = false;
        }
      }
      if (voice.gateHigh && phase >= voice.gateOffBeat)
        voice.gateHigh = false;
    }

    for (auto &state : cvOutputs) {
      if (state.interpolation == tfseq::CvInterpolation::Step ||
          state.endBeat <= state.beginBeat) {
        state.value = state.target;
      } else {
        float amount = static_cast<float>((phase - state.beginBeat) /
                                          (state.endBeat - state.beginBeat));
        amount = std::clamp(amount, 0.f, 1.f);
        if (state.interpolation == tfseq::CvInterpolation::Smooth)
          amount = amount * amount * (3.f - 2.f * amount);
        else if (state.interpolation == tfseq::CvInterpolation::Power)
          amount = std::pow(amount, state.power);
        state.value = state.from + (state.target - state.from) * amount;
      }
      if (running) {
        const double beatDelta =
            periodKnown && periodSamples > 0.0 ? 1.0 / periodSamples : 0.0;
        const float envelope = state.envelope.process(
            outputVoices[0].gateHigh, state.envelopeTriggerPending,
            state.envelopePeak, state.envelopeSpec, args.sampleTime, beatDelta);
        state.envelopeTriggerPending = false;
        state.output =
            tfseq::CvEnvelopeOutput(state.value, envelope, state.envelopeSpec);
      }
    }

    for (const auto output : {PITCH_OUTPUT, GATE_OUTPUT, TRIGGER_OUTPUT,
                              VELOCITY_OUTPUT, ACCENT_OUTPUT})
      outputs[output].setChannels(static_cast<int>(outputVoiceCount));
    for (std::size_t index = 0; index < outputVoiceCount; ++index) {
      auto &voice = outputVoices[index];
      outputs[PITCH_OUTPUT].setVoltage(voice.pitch, static_cast<int>(index));
      outputs[GATE_OUTPUT].setVoltage(running && voice.gateHigh ? 10.f : 0.f,
                                      static_cast<int>(index));
      outputs[TRIGGER_OUTPUT].setVoltage(
          running && voice.triggerPulse.process(args.sampleTime) ? 10.f : 0.f,
          static_cast<int>(index));
      outputs[VELOCITY_OUTPUT].setVoltage(voice.velocity,
                                          static_cast<int>(index));
      outputs[ACCENT_OUTPUT].setVoltage(running && voice.gateHigh ? voice.accent
                                                                  : 0.f,
                                        static_cast<int>(index));
    }
    constexpr std::array<int, tfseq::CvLaneCount> cvOutputIds{
        CV1_OUTPUT, CV2_OUTPUT, CV3_OUTPUT};
    for (std::size_t index = 0; index < cvOutputIds.size(); ++index) {
      outputs[cvOutputIds[index]].setChannels(1);
      outputs[cvOutputIds[index]].setVoltage(cvOutputs[index].output);
      visibleCvValues[index].store(cvOutputs[index].output,
                                   std::memory_order_relaxed);
    }
    const bool transportActive = running && activeProgram;
    const bool clockFlash = running && clockLightPulse.process(args.sampleTime);
    lights[RUN_LIGHT].setBrightness(
        tfui::transportLightBrightness(transportActive, clockFlash));
    visibleBeat.store(clockSeen ? phase - programStartBeat : 0.0,
                      std::memory_order_relaxed);
    if (running && clockSeen)
      ++samplesSinceClock;
    wasRunning = running;
  }
};

struct TfProgSequencerSourceChange : history::ModuleAction {
  std::string oldSource;
  std::string newSource;
  int oldCursor = 0;
  int oldSelection = 0;
  int newCursor = 0;
  int newSelection = 0;

  TfProgSequencerSourceChange() { name = "edit sequencer program"; }

  void apply(const std::string &source, int cursor, int selection) {
    if (!APP || !APP->scene || !APP->scene->rack)
      return;
    auto *moduleWidget = APP->scene->rack->getModule(moduleId);
    auto *module =
        moduleWidget ? moduleWidget->getModule<TfProgSequencer>() : nullptr;
    if (module)
      module->restoreEditedSource(source, cursor, selection);
  }

  void undo() override { apply(oldSource, oldCursor, oldSelection); }
  void redo() override { apply(newSource, newCursor, newSelection); }
};

struct TfSequenceEditor : app::LedDisplayTextField {
  static constexpr std::size_t CursorLaneCount =
      static_cast<std::size_t>(tfseq::CursorLane::Count);
  static constexpr std::size_t CursorTrailCapacity = 8;
  static constexpr std::size_t CursorMotionCapacity = 4;
  static constexpr std::size_t CursorBloomCapacity = 4;
  static constexpr int CursorGlowSamples = 14;
  // Enough fixed storage for the full six-second window at the capped 90 Hz
  // sampling rate, with a small margin for timer jitter.
  static constexpr std::size_t CvTraceCapacity = 576;
  static constexpr double CvTraceWindowSeconds = 6.0;
  static constexpr double CvTraceSampleIntervalSeconds = 1.0 / 90.0;
  static constexpr float CvTraceVoltageRange = 5.f;
  static constexpr float MinimumCvTraceWidth = 42.f;
  static constexpr int MaximumCvTraceCharacters = 12;

  struct CursorTrailPoint {
    tfseq::SourceSpan span;
    int drawBegin = -1;
    int drawEnd = -1;
    double pulsedAt = 0.0;
    float tailEnergy = 0.f;
  };

  struct CursorMotion {
    tfseq::SourceSpan from;
    tfseq::SourceSpan to;
    double startedAt = 0.0;
    double duration = 0.0;
    double tailDuration = 0.0;

    bool valid() const noexcept {
      return from.valid() && to.valid() && duration > 0.0 && tailDuration > 0.0;
    }
  };

  struct CursorBloom {
    tfseq::SourceSpan span;
    double startedAt = 0.0;
    double expansionDuration = 0.0;
    double tailDuration = 0.0;
    float strength = 0.f;

    bool valid() const noexcept {
      return span.valid() && expansionDuration > 0.0 && tailDuration > 0.0 &&
             strength > 0.f;
    }
  };

  struct GlyphBox {
    float x = 0.f;
    float y = 0.f;
    float width = 0.f;
    float height = 0.f;
    bool valid = false;
  };

  struct CvTraceSample {
    double sampledAt = 0.0;
    float value = 0.f;
  };

  TfProgSequencer *module = nullptr;
  std::uint64_t loadedRevision = 0;
  std::array<std::uint64_t, CursorLaneCount> observedPulses{};
  std::array<tfseq::SourceSpan, CursorLaneCount> observedSpans{};
  std::array<std::array<CursorTrailPoint, CursorTrailCapacity>, CursorLaneCount>
      cursorTrails{};
  std::array<std::array<CursorMotion, CursorMotionCapacity>, CursorLaneCount>
      cursorMotions{};
  std::array<std::array<CursorBloom, CursorBloomCapacity>, CursorLaneCount>
      cursorBlooms{};
  std::array<double, CursorLaneCount> lastCursorPulseTimes{};
  std::array<std::array<CvTraceSample, CvTraceCapacity>, tfseq::CvLaneCount>
      cvTraces{};
  std::array<std::size_t, tfseq::CvLaneCount> cvTraceWrite{};
  std::array<std::size_t, tfseq::CvLaneCount> cvTraceCount{};
  std::array<int, tfseq::CvLaneCount> cvTraceLineOrdinals{};
  double lastCvTraceSampleAt = -INFINITY;
  std::uint64_t observedExecutionPulse = 0;
  double executionPulsedAt = 0.0;
  float executionTailEnergy = 0.f;
  double lastClickTime = -INFINITY;
  Vec lastClickPosition;
  int clickCount = 0;
  int synchronizedCursor = 0;
  int synchronizedSelection = 0;
  int cursorBeforeChange = 0;
  int selectionBeforeChange = 0;
  bool changeSnapshotValid = false;
  bool suppressShortcutSpace = false;

  TfSequenceEditor() {
    multiline = true;
    cvTraceLineOrdinals.fill(-1);
    color = editorColor(tfui::HeatmapPalette::Magma, EditorIntensity::Text);
    bgColor =
        editorColor(tfui::HeatmapPalette::Magma, EditorIntensity::Background);
  }

  tfui::HeatmapPalette heatmap() const noexcept {
    return module ? module->editorHeatmap.load(std::memory_order_relaxed)
                  : tfui::HeatmapPalette::Magma;
  }

  NVGcolor heatmapColor(float intensity) const noexcept {
    return editorColor(heatmap(), intensity);
  }

  bool sourcePositionsMatchActiveProgram() const noexcept {
    return module &&
           module->pendingProgram.load(std::memory_order_acquire) == 0 &&
           module->source == module->evaluatedSource;
  }

  void draw(const DrawArgs &args) override {
    nvgBeginPath(args.vg);
    nvgRect(args.vg, 0.f, 0.f, box.size.x, box.size.y);
    nvgFillColor(args.vg, bgColor);
    nvgFill(args.vg);
    app::LedDisplayTextField::draw(args);
  }

  void step() override {
    app::LedDisplayTextField::step();
    if (!module)
      return;
    color = heatmapColor(EditorIntensity::Text);
    bgColor = heatmapColor(EditorIntensity::Background);
    module->collectRetired();
    const auto revision =
        module->sourceRevision.load(std::memory_order_acquire);
    if (revision != loadedRevision) {
      setText(module->source);
      const int restoredCursor =
          module->restoredEditorCursor.exchange(-1, std::memory_order_relaxed);
      const int restoredSelection = module->restoredEditorSelection.exchange(
          -1, std::memory_order_relaxed);
      if (restoredCursor >= 0 && restoredSelection >= 0) {
        cursor = std::clamp(restoredCursor, 0, static_cast<int>(text.size()));
        selection =
            std::clamp(restoredSelection, 0, static_cast<int>(text.size()));
      }
      synchronizedCursor = cursor;
      synchronizedSelection = selection;
      clearSourcePositionVisualization();
      loadedRevision = revision;
    }
    sampleCvTraces(system::getTime());
  }

  void clearSourcePositionVisualization() noexcept {
    cursorTrails = {};
    cursorMotions = {};
    cursorBlooms = {};
    lastCursorPulseTimes = {};
    observedSpans = {};
  }

  void sampleCvTraces(double now) noexcept {
    if (!module || !std::isfinite(now))
      return;
    if (lastCvTraceSampleAt > 0.0 &&
        now - lastCvTraceSampleAt < CvTraceSampleIntervalSeconds)
      return;
    if (lastCvTraceSampleAt > 0.0 &&
        now - lastCvTraceSampleAt > 2.0 * CvTraceWindowSeconds) {
      cvTraces = {};
      cvTraceWrite = {};
      cvTraceCount = {};
    }
    lastCvTraceSampleAt = now;
    for (std::size_t lane = 0; lane < cvTraces.size(); ++lane) {
      auto &write = cvTraceWrite[lane];
      cvTraces[lane][write] = {
          now, module->visibleCvValues[lane].load(std::memory_order_relaxed)};
      write = (write + 1) % CvTraceCapacity;
      cvTraceCount[lane] =
          std::min<std::size_t>(cvTraceCount[lane] + 1, CvTraceCapacity);
    }
  }

  void synchronizeEditedSource(int previousCursor, int previousSelection) {
    if (module) {
      const std::string editedSource = getText();
      if (editedSource != module->source) {
        auto *change = new TfProgSequencerSourceChange;
        change->moduleId = module->id;
        change->oldSource = module->source;
        change->newSource = editedSource;
        change->oldCursor = previousCursor;
        change->oldSelection = previousSelection;
        change->newCursor = cursor;
        change->newSelection = selection;
        module->source = editedSource;
        if (APP && APP->history && module->id >= 0)
          APP->history->push(change);
        else
          delete change;
      }
      // Compiled spans refer to the previous source layout. Let the next clock
      // event establish fresh cursor positions rather than drawing stale ones.
      // CV samples are independent of character positions, so editing must not
      // erase the scrolling trace history.
      clearSourcePositionVisualization();
    }
    synchronizedCursor = cursor;
    synchronizedSelection = selection;
  }

  void onChange(const ChangeEvent &) override {
    synchronizeEditedSource(changeSnapshotValid ? cursorBeforeChange
                                                : synchronizedCursor,
                            changeSnapshotValid ? selectionBeforeChange
                                                : synchronizedSelection);
  }

  void applyEdit(tfseq::editor::EditResult edit) {
    if (!edit.changed)
      return;
    const int previousCursor = cursor;
    const int previousSelection = selection;
    text = std::move(edit.text);
    selection = edit.selection;
    cursor = edit.cursor;
    synchronizeEditedSource(previousCursor, previousSelection);
  }

  static bool sameSpan(const tfseq::SourceSpan &left,
                       const tfseq::SourceSpan &right) {
    return left.begin == right.begin && left.end == right.end;
  }

  static float trailPersistence(const CursorTrailPoint &point,
                                double now) noexcept {
    if (!(point.pulsedAt > 0.0))
      return 0.f;
    const double age = std::max(0.0, now - point.pulsedAt);
    const float head = tfui::exponentialDecay(age, tfui::CursorHeadSeconds);
    const float tail =
        point.tailEnergy * tfui::exponentialDecay(age, tfui::CursorTailSeconds);
    return 0.72f * head + 0.28f * tail;
  }

  static void pulseTrailPoint(CursorTrailPoint &point, double now,
                              float impulse) noexcept {
    const double age = point.pulsedAt > 0.0 ? now - point.pulsedAt : 1e9;
    point.tailEnergy = tfui::accumulatedTail(point.tailEnergy, age, impulse);
    point.pulsedAt = now;
  }

  CursorTrailPoint &trailPointFor(std::size_t lane,
                                  const tfseq::SourceSpan &span, double now) {
    auto &trail = cursorTrails[lane];
    for (auto &point : trail) {
      if (point.span.valid() && sameSpan(point.span, span))
        return point;
    }
    for (auto &point : trail) {
      if (!point.span.valid()) {
        point.span = span;
        point.drawBegin = span.begin;
        point.drawEnd = span.begin + 1;
        return point;
      }
    }
    auto *faintest = &trail.front();
    for (auto &point : trail) {
      if (trailPersistence(point, now) < trailPersistence(*faintest, now))
        faintest = &point;
    }
    *faintest = {};
    faintest->span = span;
    faintest->drawBegin = span.begin;
    faintest->drawEnd = span.begin + 1;
    return *faintest;
  }

  CursorMotion &motionSlotFor(std::size_t lane, double now) {
    auto &motions = cursorMotions[lane];
    for (auto &motion : motions) {
      if (!motion.valid())
        return motion;
    }
    auto *faintest = &motions.front();
    for (auto &motion : motions) {
      const auto intensity =
          tfui::cursorMotionIntensity(std::max(0.0, now - motion.startedAt),
                                      motion.duration, motion.tailDuration);
      const auto faintestIntensity = tfui::cursorMotionIntensity(
          std::max(0.0, now - faintest->startedAt), faintest->duration,
          faintest->tailDuration);
      if (intensity < faintestIntensity)
        faintest = &motion;
    }
    return *faintest;
  }

  static float bloomIntensity(const CursorBloom &bloom, double now) noexcept {
    if (!bloom.valid())
      return 0.f;
    return bloom.strength *
           tfui::exponentialDecay(now - bloom.startedAt, bloom.tailDuration);
  }

  CursorBloom &bloomSlotFor(std::size_t lane, double now) {
    auto &blooms = cursorBlooms[lane];
    for (auto &bloom : blooms) {
      if (!bloom.valid())
        return bloom;
    }
    auto *faintest = &blooms.front();
    for (auto &bloom : blooms) {
      if (bloomIntensity(bloom, now) < bloomIntensity(*faintest, now))
        faintest = &bloom;
    }
    return *faintest;
  }

  bool isShortSameLineMove(const tfseq::SourceSpan &from,
                           const tfseq::SourceSpan &to) const {
    if (!from.valid() || !to.valid())
      return false;
    const int begin = std::min(from.begin, to.begin);
    const int end = std::max(from.begin, to.begin);
    if (end - begin > 8 || end > static_cast<int>(text.size()))
      return false;
    return text.find('\n', static_cast<std::size_t>(begin)) >=
           static_cast<std::size_t>(end);
  }

  GlyphBox glyphBoxAt(const DrawArgs &args, int position) const {
    if (position < 0 || position >= static_cast<int>(text.size()))
      return {};
    const char *textBegin = text.c_str();
    const char *textEnd = textBegin + text.size();
    const char *target = textBegin + position;
    const char *remaining = textBegin;
    const float x = textOffset.x + BND_TEXT_RADIUS;
    float baseline = textOffset.y + BND_WIDGET_HEIGHT - BND_TEXT_PAD_DOWN;
    const float rowWidth = box.size.x - 2 * textOffset.x - 2 * BND_TEXT_RADIUS;
    float ascender = 0.f;
    float descender = 0.f;
    float lineHeight = 0.f;
    nvgTextMetrics(args.vg, &ascender, &descender, &lineHeight);

    while (remaining < textEnd) {
      NVGtextRow rows[BND_MAX_ROWS];
      const int rowCount = nvgTextBreakLines(args.vg, remaining, textEnd,
                                             rowWidth, rows, BND_MAX_ROWS);
      if (rowCount <= 0)
        break;
      for (int rowIndex = 0; rowIndex < rowCount; ++rowIndex) {
        const auto &row = rows[rowIndex];
        if (target >= row.start && target < row.end) {
          const float advance =
              nvgTextBounds(args.vg, x, baseline, row.start, target, nullptr);
          const char *next = std::min(target + 1, row.end);
          float glyphWidth = nvgTextBounds(args.vg, x + advance, baseline,
                                           target, next, nullptr);
          if (!(glyphWidth > 0.f))
            glyphWidth = 6.f;
          return {x + advance, baseline - ascender, glyphWidth,
                  std::max(1.f, ascender - descender), true};
        }
        baseline += lineHeight;
      }
      const char *next = rows[rowCount - 1].next;
      if (next <= remaining)
        break;
      remaining = next;
    }
    return {};
  }

  template <typename Callback>
  void forEachSpanBox(const DrawArgs &args, const tfseq::SourceSpan &span,
                      Callback callback) const {
    if (!span.valid() || text.empty())
      return;
    const int begin = std::clamp(span.begin, 0, static_cast<int>(text.size()));
    const int end = std::clamp(span.end, begin, static_cast<int>(text.size()));
    if (begin == end)
      return;

    const char *textBegin = text.c_str();
    const char *textEnd = textBegin + text.size();
    const char *spanBegin = textBegin + begin;
    const char *spanEnd = textBegin + end;
    const char *remaining = textBegin;
    const float x = textOffset.x + BND_TEXT_RADIUS;
    float baseline = textOffset.y + BND_WIDGET_HEIGHT - BND_TEXT_PAD_DOWN;
    const float rowWidth = box.size.x - 2 * textOffset.x - 2 * BND_TEXT_RADIUS;
    float ascender = 0.f;
    float descender = 0.f;
    float lineHeight = 0.f;
    nvgTextMetrics(args.vg, &ascender, &descender, &lineHeight);

    while (remaining < textEnd) {
      NVGtextRow rows[BND_MAX_ROWS];
      const int rowCount = nvgTextBreakLines(args.vg, remaining, textEnd,
                                             rowWidth, rows, BND_MAX_ROWS);
      if (rowCount <= 0)
        break;
      for (int rowIndex = 0; rowIndex < rowCount; ++rowIndex) {
        const auto &row = rows[rowIndex];
        const char *runBegin = std::max(row.start, spanBegin);
        const char *runEnd = std::min(row.end, spanEnd);
        if (runBegin < runEnd) {
          const float left =
              nvgTextBounds(args.vg, x, baseline, row.start, runBegin, nullptr);
          const float right =
              nvgTextBounds(args.vg, x, baseline, row.start, runEnd, nullptr);
          callback(GlyphBox{x + left, baseline - ascender,
                            std::max(2.f, right - left),
                            std::max(1.f, ascender - descender), true});
        }
        baseline += lineHeight;
      }
      const char *next = rows[rowCount - 1].next;
      if (next <= remaining)
        break;
      remaining = next;
    }
  }

  void drawActiveStep(const DrawArgs &args, const tfseq::SourceSpan &span,
                      float progress) const {
    if (!span.valid())
      return;
    forEachSpanBox(args, span, [&](const GlyphBox &glyph) {
      nvgSave(args.vg);
      nvgGlobalCompositeOperation(args.vg, NVG_LIGHTER);
      NVGcolor activeFill = heatmapColor(0.30f);
      activeFill.a = 0.13f;
      nvgBeginPath(args.vg);
      nvgRoundedRect(args.vg, glyph.x - 1.f, glyph.y - 0.5f, glyph.width + 2.f,
                     glyph.height + 1.f, 1.5f);
      nvgFillColor(args.vg, activeFill);
      nvgFill(args.vg);

      const float underlineY = glyph.y + glyph.height - 0.6f;
      NVGcolor track = heatmapColor(0.22f);
      track.a = 0.30f;
      nvgBeginPath(args.vg);
      nvgRect(args.vg, glyph.x, underlineY, glyph.width, 0.7f);
      nvgFillColor(args.vg, track);
      nvgFill(args.vg);

      const float completed = glyph.width * std::clamp(progress, 0.f, 1.f);
      NVGcolor beam = heatmapColor(0.88f);
      beam.a = 0.76f;
      nvgBeginPath(args.vg);
      nvgRoundedRect(args.vg, glyph.x, underlineY - 0.3f,
                     std::max(0.8f, completed), 1.3f, 0.65f);
      nvgFillColor(args.vg, beam);
      nvgFill(args.vg);

      const float headX = glyph.x + completed;
      NVGcolor headCore = heatmapColor(0.96f);
      NVGcolor headEdge = headCore;
      headCore.a = 0.72f;
      headEdge.a = 0.f;
      const auto headPaint = nvgRadialGradient(args.vg, headX, underlineY, 0.4f,
                                               4.f, headCore, headEdge);
      nvgBeginPath(args.vg);
      nvgCircle(args.vg, headX, underlineY, 4.f);
      nvgFillPaint(args.vg, headPaint);
      nvgFill(args.vg);
      nvgRestore(args.vg);
    });
  }

  std::vector<tfseq::SourceSpan> cvLineSpans(std::size_t lane) const {
    std::vector<tfseq::SourceSpan> spans;
    if (lane >= tfseq::CvLaneCount)
      return spans;
    const std::string label = "cv" + std::to_string(lane + 1);
    std::size_t lineBegin = 0;
    while (lineBegin <= text.size()) {
      const auto newline = text.find('\n', lineBegin);
      const auto lineEnd = newline == std::string::npos ? text.size() : newline;
      auto tokenBegin = lineBegin;
      while (tokenBegin < lineEnd &&
             (text[tokenBegin] == ' ' || text[tokenBegin] == '\t'))
        ++tokenBegin;
      const auto afterLabel = tokenBegin + label.size();
      if (afterLabel <= lineEnd &&
          text.compare(tokenBegin, label.size(), label) == 0 &&
          (afterLabel == lineEnd || text[afterLabel] == ' ' ||
           text[afterLabel] == '\t')) {
        spans.push_back(
            {static_cast<int>(lineBegin), static_cast<int>(lineEnd), 1, 1});
      }
      if (newline == std::string::npos)
        break;
      lineBegin = newline + 1;
    }
    return spans;
  }

  tfseq::SourceSpan cvTraceLineSpan(std::size_t lane,
                                    const tfseq::SourceSpan &activeSpan,
                                    bool sourcePositionsCurrent) {
    const auto spans = cvLineSpans(lane);
    if (spans.empty())
      return {};

    if (sourcePositionsCurrent && activeSpan.valid()) {
      for (std::size_t index = 0; index < spans.size(); ++index) {
        if (activeSpan.begin < spans[index].begin ||
            activeSpan.begin > spans[index].end)
          continue;
        cvTraceLineOrdinals[lane] = static_cast<int>(index);
        return spans[index];
      }
    }

    const int ordinal = cvTraceLineOrdinals[lane];
    if (ordinal >= 0 && ordinal < static_cast<int>(spans.size()))
      return spans[static_cast<std::size_t>(ordinal)];

    if (!activeSpan.valid())
      return {};
    const auto nearest = std::min_element(
        spans.begin(), spans.end(),
        [&](const tfseq::SourceSpan &left, const tfseq::SourceSpan &right) {
          return std::abs(left.begin - activeSpan.begin) <
                 std::abs(right.begin - activeSpan.begin);
        });
    cvTraceLineOrdinals[lane] =
        static_cast<int>(std::distance(spans.begin(), nearest));
    return *nearest;
  }

  GlyphBox inlineTraceBox(const DrawArgs &args,
                          const tfseq::SourceSpan &lineSpan) const {
    if (!lineSpan.valid())
      return {};
    GlyphBox finalRow;
    forEachSpanBox(args, lineSpan,
                   [&](const GlyphBox &row) { finalRow = row; });
    if (!finalRow.valid)
      return {};
    const float right = box.size.x - textOffset.x - BND_TEXT_RADIUS;
    const float left = finalRow.x + finalRow.width + 7.f;
    const float availableWidth = right - left;
    if (availableWidth < MinimumCvTraceWidth)
      return {};
    static constexpr char WidthReference[] = "000000000000";
    static_assert(sizeof(WidthReference) - 1 == MaximumCvTraceCharacters,
                  "CV trace width reference must match its character cap");
    const float maximumWidth =
        nvgTextBounds(args.vg, 0.f, 0.f, WidthReference, nullptr, nullptr);
    return {left, finalRow.y + 0.5f, std::min(availableWidth, maximumWidth),
            std::max(6.f, finalRow.height - 1.f), true};
  }

  void drawCvTrace(const DrawArgs &args, std::size_t lane, double now,
                   const tfseq::SourceSpan &activeSpan,
                   bool sourcePositionsCurrent) {
    if (lane >= cvTraces.size() || cvTraceCount[lane] < 2)
      return;
    const auto lineSpan =
        cvTraceLineSpan(lane, activeSpan, sourcePositionsCurrent);
    const auto plot = inlineTraceBox(args, lineSpan);
    if (!plot.valid)
      return;

    const auto count = cvTraceCount[lane];
    const auto first =
        (cvTraceWrite[lane] + CvTraceCapacity - count) % CvTraceCapacity;
    float minimum = std::numeric_limits<float>::infinity();
    float maximum = -std::numeric_limits<float>::infinity();
    for (std::size_t offset = 0; offset < count; ++offset) {
      const auto &sample = cvTraces[lane][(first + offset) % CvTraceCapacity];
      const double age = now - sample.sampledAt;
      if (age < 0.0 || age > CvTraceWindowSeconds ||
          !std::isfinite(sample.value))
        continue;
      minimum = std::min(minimum, sample.value);
      maximum = std::max(maximum, sample.value);
    }
    if (!std::isfinite(minimum) || !std::isfinite(maximum))
      return;
    const auto polarity = tfui::cvTracePolarity(minimum, maximum);

    auto coordinates = [&](const CvTraceSample &sample) {
      const double age = std::max(0.0, now - sample.sampledAt);
      const float x =
          plot.x +
          plot.width * static_cast<float>(1.0 - age / CvTraceWindowSeconds);
      const float y = plot.y + plot.height * tfui::cvTraceVerticalFraction(
                                                 sample.value, polarity,
                                                 CvTraceVoltageRange);
      return Vec(x, y);
    };

    nvgSave(args.vg);
    nvgIntersectScissor(args.vg, plot.x - 4.f, plot.y - 4.f, plot.width + 8.f,
                        plot.height + 8.f);
    NVGcolor zero = heatmapColor(0.22f);
    zero.a = 0.20f;
    nvgBeginPath(args.vg);
    const float zeroY =
        plot.y + plot.height * tfui::cvTraceZeroFraction(polarity);
    nvgMoveTo(args.vg, plot.x, zeroY);
    nvgLineTo(args.vg, plot.x + plot.width, zeroY);
    nvgStrokeWidth(args.vg, 0.6f);
    nvgStrokeColor(args.vg, zero);
    nvgStroke(args.vg);

    auto strokeTrace = [&](float width, float intensity, float alpha) {
      bool started = false;
      nvgBeginPath(args.vg);
      for (std::size_t offset = 0; offset < count; ++offset) {
        const auto &sample = cvTraces[lane][(first + offset) % CvTraceCapacity];
        const double age = now - sample.sampledAt;
        if (age < 0.0 || age > CvTraceWindowSeconds)
          continue;
        const auto point = coordinates(sample);
        if (!started) {
          nvgMoveTo(args.vg, point.x, point.y);
          started = true;
        } else {
          nvgLineTo(args.vg, point.x, point.y);
        }
      }
      if (!started)
        return;
      NVGcolor color = heatmapColor(intensity);
      color.a = alpha;
      nvgStrokeWidth(args.vg, width);
      nvgLineCap(args.vg, NVG_ROUND);
      nvgLineJoin(args.vg, NVG_ROUND);
      nvgStrokeColor(args.vg, color);
      nvgStroke(args.vg);
    };

    nvgGlobalCompositeOperation(args.vg, NVG_LIGHTER);
    strokeTrace(5.0f, 0.58f, 0.055f);
    strokeTrace(2.4f, 0.76f, 0.16f);
    strokeTrace(0.9f, 0.96f, 0.78f);

    const auto newestIndex =
        (cvTraceWrite[lane] + CvTraceCapacity - 1) % CvTraceCapacity;
    const auto head = coordinates(cvTraces[lane][newestIndex]);
    NVGcolor headCore = heatmapColor(1.f);
    NVGcolor headEdge = headCore;
    headCore.a = 0.62f;
    headEdge.a = 0.f;
    const auto headPaint = nvgRadialGradient(args.vg, head.x, head.y, 0.5f,
                                             4.5f, headCore, headEdge);
    nvgBeginPath(args.vg);
    nvgCircle(args.vg, head.x, head.y, 4.5f);
    nvgFillPaint(args.vg, headPaint);
    nvgFill(args.vg);
    nvgRestore(args.vg);
  }

  bool drawCursorMotion(const DrawArgs &args, CursorMotion &motion, double now,
                        tfui::CursorTravelCurve curve) const {
    if (!motion.valid())
      return false;
    const double age = std::max(0.0, now - motion.startedAt);
    const float linearProgress =
        static_cast<float>(std::clamp(age / motion.duration, 0.0, 1.0));
    const float positionProgress =
        tfui::cursorTravelPosition(curve, linearProgress);
    const float intensity =
        tfui::cursorMotionIntensity(age, motion.duration, motion.tailDuration);
    if (intensity < 0.004f) {
      motion = {};
      return false;
    }

    const auto from = glyphBoxAt(args, motion.from.begin);
    const auto to = glyphBoxAt(args, motion.to.begin);
    if (!from.valid || !to.valid ||
        std::abs(from.y - to.y) > std::max(from.height, to.height) * 0.25f) {
      motion = {};
      return false;
    }

    const float fromCenter = from.x + 0.5f * from.width;
    const float toCenter = to.x + 0.5f * to.width;
    const float headCenter =
        fromCenter + (toCenter - fromCenter) * positionProgress;
    const float headWidth =
        from.width + (to.width - from.width) * positionProgress;
    const float top = from.y + (to.y - from.y) * positionProgress;
    const float height =
        from.height + (to.height - from.height) * positionProgress;

    // Render a bounded set of overlapping sub-pixel phosphor deposits. Their
    // additive overlap makes whitespace illuminate continuously; keeping each
    // motion alive after arrival gives successive clock events independent,
    // naturally overlapping trails without any frame-to-frame integration.
    nvgSave(args.vg);
    nvgGlobalCompositeOperation(args.vg, NVG_LIGHTER);
    const float travelledPixels = std::abs(headCenter - fromCenter);
    const int sampleCount =
        std::clamp(static_cast<int>(std::ceil(travelledPixels / 2.5f)), 0,
                   CursorGlowSamples);
    for (int sample = 0; sample < sampleCount; ++sample) {
      const float fraction =
          static_cast<float>(sample + 1) / static_cast<float>(sampleCount + 1);
      // Sample the beam's actual earlier timestamps. Smoothstep therefore
      // changes the distribution of deposits along the path, while their
      // persistence remains a real exponential decay in wall-clock time.
      const float sampleTimeProgress = linearProgress * (1.f - fraction);
      const float sampleProgress =
          tfui::cursorTravelPosition(curve, sampleTimeProgress);
      const float sampleCenter =
          fromCenter + (toCenter - fromCenter) * sampleProgress;
      const float sampleWidth =
          from.width + (to.width - from.width) * sampleProgress;
      const float sampleTop = from.y + (to.y - from.y) * sampleProgress;
      const float sampleHeight =
          from.height + (to.height - from.height) * sampleProgress;
      const double depositedAt =
          motion.startedAt + motion.duration * sampleTimeProgress;
      const float depositedIntensity =
          tfui::cursorMotionIntensity(motion.duration * sampleTimeProgress,
                                      motion.duration, motion.tailDuration);
      const float sampleIntensity = std::clamp(
          depositedIntensity *
              tfui::exponentialDecay(now - depositedAt, motion.tailDuration),
          0.f, 1.f);
      NVGcolor sampleCore = heatmapColor(sampleIntensity);
      NVGcolor sampleEdge = heatmapColor(0.12f * sampleIntensity);
      sampleCore.a = 0.16f;
      sampleEdge.a = 0.f;
      const float sampleLeft = sampleCenter - 0.5f * sampleWidth;
      const auto samplePaint =
          nvgBoxGradient(args.vg, sampleLeft, sampleTop, sampleWidth,
                         sampleHeight, 1.5f, 5.5f, sampleCore, sampleEdge);
      nvgBeginPath(args.vg);
      nvgRect(args.vg, sampleLeft - 5.5f, sampleTop - 5.5f, sampleWidth + 11.f,
              sampleHeight + 11.f);
      nvgFillPaint(args.vg, samplePaint);
      nvgFill(args.vg);
    }

    // The beam head, rather than an already-bright destination caret, is the
    // dominant cursor. It therefore moves visibly at fractional glyph offsets.
    const float headLeft = headCenter - 0.5f * headWidth;
    NVGcolor headCore = heatmapColor(intensity);
    NVGcolor headEdge = heatmapColor(0.14f * intensity);
    headCore.a = 0.68f;
    headEdge.a = 0.f;
    const auto headPaint =
        nvgBoxGradient(args.vg, headLeft, top, headWidth, height, 1.5f, 7.f,
                       headCore, headEdge);
    nvgBeginPath(args.vg);
    nvgRect(args.vg, headLeft - 7.f, top - 7.f, headWidth + 14.f,
            height + 14.f);
    nvgFillPaint(args.vg, headPaint);
    nvgFill(args.vg);
    NVGcolor headFill = heatmapColor(std::min(1.f, intensity + 0.08f));
    headFill.a = 0.34f;
    nvgBeginPath(args.vg);
    nvgRoundedRect(args.vg, headLeft, top, headWidth, height, 1.5f);
    nvgFillColor(args.vg, headFill);
    nvgFill(args.vg);
    nvgRestore(args.vg);
    return true;
  }

  void drawDiffusionRect(const DrawArgs &args, const GlyphBox &source,
                         float expansion, float intensity,
                         float opacityScale = 1.f) const {
    if (!source.valid || intensity < 0.004f || opacityScale <= 0.f)
      return;
    const float maximumSpread = 2.25f + 12.25f * expansion;
    // Diffuse energy directly from the caret's rectangular footprint. NanoVG's
    // box-gradient shader becomes unreliable when its corner radius approaches
    // a glyph's tiny width, so use a bounded stack of sub-pixel rounded fills.
    // The layers are closer than one pixel at full spread and antialias into a
    // smooth field: rectangular and dense near the source, increasingly round,
    // faint, and diffuse outward. There are no stroked contours.
    nvgSave(args.vg);
    nvgGlobalCompositeOperation(args.vg, NVG_LIGHTER);
    constexpr int DiffusionLayers = 18;
    for (int layer = DiffusionLayers; layer >= 1; --layer) {
      const float distance =
          static_cast<float>(layer) / static_cast<float>(DiffusionLayers);
      const float spread = maximumSpread * distance;
      const float falloff = std::exp(-1.65f * distance * distance);
      const float level = (0.52f + 0.40f * falloff) * intensity;
      NVGcolor layerColor = heatmapColor(level);
      layerColor.a = opacityScale * (0.0155f + 0.0085f * (1.f - distance));
      const float cornerRadius = 1.2f + spread * (0.42f + 0.38f * distance);
      nvgBeginPath(args.vg);
      nvgRoundedRect(args.vg, source.x - spread, source.y - spread,
                     source.width + 2.f * spread, source.height + 2.f * spread,
                     cornerRadius);
      nvgFillColor(args.vg, layerColor);
      nvgFill(args.vg);
    }
    nvgRestore(args.vg);
  }

  bool drawCursorBloom(const DrawArgs &args, CursorBloom &bloom,
                       double now) const {
    if (!bloom.valid())
      return false;
    const double age = std::max(0.0, now - bloom.startedAt);
    const float intensity = bloomIntensity(bloom, now);
    if (intensity < 0.004f) {
      bloom = {};
      return false;
    }
    const auto glyph = glyphBoxAt(args, bloom.span.begin);
    if (!glyph.valid) {
      bloom = {};
      return false;
    }
    drawDiffusionRect(args, glyph,
                      tfui::cursorBloomExpansion(age, bloom.expansionDuration),
                      intensity);
    return true;
  }

  std::size_t commentStartFor(std::size_t rowStart) const {
    const auto previousNewline =
        rowStart > 0 ? text.rfind('\n', rowStart - 1) : std::string::npos;
    const auto lineBegin =
        previousNewline == std::string::npos ? 0 : previousNewline + 1;
    const auto lineEnd = text.find('\n', lineBegin);
    const auto comment = text.find("//", lineBegin);
    if (comment == std::string::npos ||
        (lineEnd != std::string::npos && comment >= lineEnd))
      return text.size();
    return comment;
  }

  void drawTextWithDimComments(const DrawArgs &args) {
    const char *textBegin = text.c_str();
    const char *textEnd = textBegin + text.size();
    const char *remaining = textBegin;
    const float x = textOffset.x + BND_TEXT_RADIUS;
    float baseline = textOffset.y + BND_WIDGET_HEIGHT - BND_TEXT_PAD_DOWN;
    const float width = box.size.x - 2 * textOffset.x - 2 * BND_TEXT_RADIUS;
    nvgFontSize(args.vg, 12.f);
    nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_BASELINE);
    float lineHeight = 0.f;
    nvgTextMetrics(args.vg, nullptr, nullptr, &lineHeight);
    const NVGcolor commentColor = heatmapColor(EditorIntensity::Comment);

    while (remaining < textEnd) {
      NVGtextRow rows[BND_MAX_ROWS];
      const int rowCount = nvgTextBreakLines(args.vg, remaining, textEnd, width,
                                             rows, BND_MAX_ROWS);
      if (rowCount <= 0)
        break;
      for (int rowIndex = 0; rowIndex < rowCount; ++rowIndex) {
        const auto &row = rows[rowIndex];
        const auto rowStart = static_cast<std::size_t>(row.start - textBegin);
        const auto rowEnd = static_cast<std::size_t>(row.end - textBegin);
        const auto commentStart = commentStartFor(rowStart);
        if (commentStart >= rowEnd) {
          nvgFillColor(args.vg, color);
          nvgText(args.vg, x, baseline, row.start, row.end);
        } else {
          float runX = x;
          if (commentStart > rowStart) {
            nvgFillColor(args.vg, color);
            runX = nvgText(args.vg, runX, baseline, row.start,
                           textBegin + commentStart);
          }
          nvgFillColor(args.vg, commentColor);
          nvgText(args.vg, runX, baseline,
                  textBegin + std::max(rowStart, commentStart), row.end);
        }
        baseline += lineHeight;
      }
      const char *next = rows[rowCount - 1].next;
      if (next <= remaining)
        break;
      remaining = next;
    }
  }

  static bool isWordCharacter(char character) {
    const auto value = static_cast<unsigned char>(character);
    return std::isalnum(value) || character == '_';
  }

  void selectWordAt(int position) {
    if (text.empty())
      return;
    position = std::clamp(position, 0, static_cast<int>(text.size()) - 1);
    if (!isWordCharacter(text[static_cast<std::size_t>(position)]) &&
        position > 0 &&
        isWordCharacter(text[static_cast<std::size_t>(position - 1)]))
      --position;
    int begin = position;
    int end = position;
    while (begin > 0 &&
           isWordCharacter(text[static_cast<std::size_t>(begin - 1)]))
      --begin;
    while (end < static_cast<int>(text.size()) &&
           isWordCharacter(text[static_cast<std::size_t>(end)]))
      ++end;
    if (begin == end)
      end = std::min(begin + 1, static_cast<int>(text.size()));
    selection = begin;
    cursor = end;
  }

  void selectLineAt(int position) {
    const auto span = lineSpanAt(position);
    selection = span.first;
    cursor = span.second;
  }

  std::pair<int, int> lineSpanAt(int position) const {
    position = std::clamp(position, 0, static_cast<int>(text.size()));
    const auto begin =
        position > 0 ? text.rfind('\n', static_cast<std::size_t>(position - 1))
                     : std::string::npos;
    const auto end = text.find('\n', static_cast<std::size_t>(position));
    return {begin == std::string::npos ? 0 : static_cast<int>(begin + 1),
            end == std::string::npos ? static_cast<int>(text.size())
                                     : static_cast<int>(end)};
  }

  void onSelectText(const SelectTextEvent &event) override {
    if (suppressShortcutSpace) {
      suppressShortcutSpace = false;
      if (event.codepoint == static_cast<std::uint32_t>(' ')) {
        event.consume(this);
        return;
      }
    }
    cursorBeforeChange = cursor;
    selectionBeforeChange = selection;
    changeSnapshotValid = true;
    app::LedDisplayTextField::onSelectText(event);
    changeSnapshotValid = false;
    synchronizedCursor = cursor;
    synchronizedSelection = selection;
  }

  bool requestConnectedTransport(tftransport::Command command,
                                 const char *action) {
    if (!module || !APP || !APP->scene || !APP->scene->rack)
      return false;
    auto *moduleWidget = APP->scene->rack->getModule(module->id);
    auto *runPort = moduleWidget
                        ? moduleWidget->getInput(TfProgSequencer::RUN_INPUT)
                        : nullptr;
    if (runPort) {
      for (auto *cableWidget :
           APP->scene->rack->getCompleteCablesOnPort(runPort)) {
        auto *cable = cableWidget ? cableWidget->getCable() : nullptr;
        if (cable && cable->inputModule == module &&
            cable->inputId == TfProgSequencer::RUN_INPUT &&
            tftransport::RequestModuleCommand(cable->outputModule,
                                              cable->outputId, command)) {
          module->compileMessage = std::string("TRANSPORT - ") + action;
          return true;
        }
      }
    }
    module->compileMessage =
        "TRANSPORT - connect RUN directly to TriggerFish Transport";
    return false;
  }

  void onButton(const ButtonEvent &event) override {
    app::LedDisplayTextField::onButton(event);
    synchronizedCursor = cursor;
    synchronizedSelection = selection;
    if (event.action != GLFW_PRESS || event.button != GLFW_MOUSE_BUTTON_LEFT)
      return;
    const double now = system::getTime();
    const bool nearby = event.pos.minus(lastClickPosition).norm() <= 5.f;
    if (now - lastClickTime <= 0.42 && nearby)
      ++clickCount;
    else
      clickCount = 1;
    lastClickTime = now;
    lastClickPosition = event.pos;
    const int position = cursor;
    if (clickCount == 2)
      selectWordAt(position);
    else if (clickCount >= 3) {
      selectLineAt(position);
      clickCount = 0;
    }
    synchronizedCursor = cursor;
    synchronizedSelection = selection;
  }

  void onDragHover(const DragHoverEvent &event) override {
    app::LedDisplayTextField::onDragHover(event);
    synchronizedCursor = cursor;
    synchronizedSelection = selection;
  }

  void onSelectKey(const SelectKeyEvent &event) override {
    if (event.key == GLFW_KEY_SPACE && event.action == GLFW_RELEASE)
      suppressShortcutSpace = false;
    if (event.action == GLFW_PRESS && APP && APP->history &&
        (event.isKeyCommand(GLFW_KEY_Z, RACK_MOD_CTRL | RACK_MOD_SHIFT) ||
         event.isKeyCommand(GLFW_KEY_Y, RACK_MOD_CTRL))) {
      APP->history->redo();
      event.consume(this);
      return;
    }
    if (event.action == GLFW_PRESS && APP && APP->history &&
        event.isKeyCommand(GLFW_KEY_Z, RACK_MOD_CTRL)) {
      APP->history->undo();
      event.consume(this);
      return;
    }
    if (module && event.action == GLFW_PRESS &&
        event.isKeyCommand(GLFW_KEY_PERIOD, RACK_MOD_CTRL | RACK_MOD_SHIFT)) {
      module->source = getText();
      if (module->publishSource(module->source, true))
        module->flashExecution(0, static_cast<int>(module->source.size()));
      event.consume(this);
      return;
    }
    if (module && event.action == GLFW_PRESS &&
        event.isKeyCommand(GLFW_KEY_PERIOD, RACK_MOD_CTRL)) {
      module->source = getText();
      if (module->publishSource(module->source))
        module->flashExecution(0, static_cast<int>(module->source.size()));
      event.consume(this);
      return;
    }
    if (module && event.action == GLFW_PRESS &&
        (event.isKeyCommand(GLFW_KEY_ENTER, RACK_MOD_CTRL) ||
         event.isKeyCommand(GLFW_KEY_KP_ENTER, RACK_MOD_CTRL))) {
      module->source = getText();
      if (cursor != selection) {
        const int begin = std::min(cursor, selection);
        module->publishSelection(getSelectedText(), begin);
      } else {
        const auto span = lineSpanAt(cursor);
        module->publishSelection(
            text.substr(static_cast<std::size_t>(span.first),
                        static_cast<std::size_t>(span.second - span.first)),
            span.first);
      }
      event.consume(this);
      return;
    }
    if (module && event.action == GLFW_PRESS &&
        event.isKeyCommand(GLFW_KEY_SLASH, RACK_MOD_CTRL)) {
      applyEdit(tfseq::editor::ToggleLineComments(text, cursor, selection));
      event.consume(this);
      return;
    }
    if (module &&
        (event.action == GLFW_PRESS || event.action == GLFW_REPEAT) &&
        event.isKeyCommand(GLFW_KEY_SPACE,
                           RACK_MOD_CTRL | RACK_MOD_SHIFT)) {
      suppressShortcutSpace = true;
      if (event.action == GLFW_PRESS)
        requestConnectedTransport(tftransport::Command::TogglePlayPause,
                                  "PLAY/PAUSE");
      event.consume(this);
      return;
    }
    if (module && event.action == GLFW_PRESS &&
        event.isKeyCommand(GLFW_KEY_R, RACK_MOD_CTRL | RACK_MOD_SHIFT)) {
      requestConnectedTransport(tftransport::Command::PlayFromBeginning,
                                "RESTART");
      event.consume(this);
      return;
    }
    if (module && event.action == GLFW_PRESS &&
        event.isKeyCommand(GLFW_KEY_BACKSPACE,
                           RACK_MOD_CTRL | RACK_MOD_SHIFT)) {
      requestConnectedTransport(tftransport::Command::Stop, "STOP");
      event.consume(this);
      return;
    }
    if (module &&
        (event.action == GLFW_PRESS || event.action == GLFW_REPEAT) &&
        event.isKeyCommand(GLFW_KEY_SPACE, RACK_MOD_CTRL)) {
      suppressShortcutSpace = true;
      if (event.action == GLFW_PRESS) {
        const bool enabled =
            module->editorRunEnabled.load(std::memory_order_relaxed);
        module->editorRunEnabled.store(!enabled, std::memory_order_relaxed);
      }
      event.consume(this);
      return;
    }
    if (module && event.action == GLFW_PRESS &&
        event.isKeyCommand(GLFW_KEY_D, RACK_MOD_CTRL)) {
      applyEdit(tfseq::editor::Duplicate(text, cursor, selection));
      event.consume(this);
      return;
    }
    cursorBeforeChange = cursor;
    selectionBeforeChange = selection;
    changeSnapshotValid = true;
    app::LedDisplayTextField::onSelectKey(event);
    changeSnapshotValid = false;
    synchronizedCursor = cursor;
    synchronizedSelection = selection;
  }

  void drawLayer(const DrawArgs &args, int layer) override {
    if (layer == 1 && module) {
      nvgScissor(args.vg, RECT_ARGS(args.clipBox));
      NVGcolor invisibleText = color;
      invisibleText.a = 0.f;
      const double now = system::getTime();
      const bool sourcePositionsCurrent = sourcePositionsMatchActiveProgram();
      if (!sourcePositionsCurrent) {
        clearSourcePositionVisualization();
        // Ignore cursor pulses produced by the still-playing evaluated source
        // while the editor displays a different draft. The first event from a
        // newly activated program will publish a later pulse and fresh spans.
        for (std::size_t lane = 0; lane < observedPulses.size(); ++lane) {
          observedPulses[lane] =
              module->cursorPulses[lane].load(std::memory_order_acquire);
        }
      }
      // Each event gets an immediate bright head and a longer phosphor-like
      // tail. Intensities are derived from event timestamps, so a late or
      // skipped UI frame cannot alter the envelope or introduce brightness
      // jitter. All state remains fixed-size and UI-thread-only.
      const auto executionPulse =
          module->executionPulse.load(std::memory_order_acquire);
      if (executionPulse != observedExecutionPulse) {
        observedExecutionPulse = executionPulse;
        const double age =
            executionPulsedAt > 0.0 ? now - executionPulsedAt : 1e9;
        executionTailEnergy =
            tfui::accumulatedTail(executionTailEnergy, age, 0.45f);
        executionPulsedAt = now;
      }
      auto font = APP->window->loadFont(fontPath);
      if (font && font->handle >= 0) {
        bndSetFont(font->handle);
        nvgFontFaceId(args.vg, font->handle);
        nvgFontSize(args.vg, 12.f);
        nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_BASELINE);

        // CV history uses only the otherwise empty space to the right of the
        // active lane. Draw it first so the source text and cursor effects
        // remain crisp above the scope-like beam.
        for (std::size_t lane = 0; lane < tfseq::CvLaneCount; ++lane) {
          const auto span = unpackSpan(
              module->activeCvLineSpans[lane].load(std::memory_order_acquire));
          drawCvTrace(args, lane, now, span, sourcePositionsCurrent);
        }

        const auto executed =
            unpackSpan(module->executionSpan.load(std::memory_order_relaxed));
        if (executed.valid() &&
            executed.begin < static_cast<int>(text.size())) {
          const double executionAge =
              executionPulsedAt > 0.0 ? now - executionPulsedAt : 1e9;
          const float executionIntensity =
              std::clamp(0.72f * tfui::exponentialDecay(
                                     executionAge, tfui::CursorHeadSeconds) +
                             0.28f * executionTailEnergy *
                                 tfui::exponentialDecay(
                                     executionAge, tfui::CursorTailSeconds),
                         0.f, 1.f);
          const NVGcolor executionFill = heatmapColor(executionIntensity);
          bndIconLabelCaret(
              args.vg, textOffset.x, textOffset.y,
              box.size.x - 2 * textOffset.x, box.size.y - 2 * textOffset.y, -1,
              invisibleText, 12.f, text.c_str(), executionFill, executed.begin,
              std::min(executed.end, static_cast<int>(text.size())));
          // Successful Ctrl+Enter / Ctrl+. evaluation diffuses from the actual
          // executed row rectangles. Multi-line programs therefore glow as
          // text rows, not as one unrelated panel-sized box.
          nvgFontFaceId(args.vg, font->handle);
          nvgFontSize(args.vg, 12.f);
          nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_BASELINE);
          const float executionExpansion = tfui::cursorBloomExpansion(
              executionAge, tfui::ExecutionBloomExpansionSeconds);
          forEachSpanBox(args, executed, [&](const GlyphBox &rowBox) {
            drawDiffusionRect(args, rowBox, executionExpansion,
                              0.68f * executionIntensity, 0.62f);
          });
        }
        for (std::size_t lane = 0; lane < module->cursorSpans.size(); ++lane) {
          if (!sourcePositionsCurrent)
            continue;
          const auto pulse =
              module->cursorPulses[lane].load(std::memory_order_acquire);
          const auto span = unpackSpan(
              module->cursorSpans[lane].load(std::memory_order_relaxed));
          if (pulse != observedPulses[lane]) {
            const auto events = std::min<std::uint64_t>(
                4, pulse > observedPulses[lane] ? pulse - observedPulses[lane]
                                                : 1);
            observedPulses[lane] = pulse;
            const double previousPulseTime = lastCursorPulseTimes[lane];
            const double pulseInterval =
                previousPulseTime > 0.0
                    ? (now - previousPulseTime) / static_cast<double>(events)
                    : 0.18;
            lastCursorPulseTimes[lane] = now;
            if (span.valid()) {
              const auto previous = observedSpans[lane];
              bool moving = false;
              if (!sameSpan(previous, span)) {
                if (isShortSameLineMove(previous, span)) {
                  const double travelDuration =
                      tfui::cursorTravelDuration(pulseInterval);
                  auto &motion = motionSlotFor(lane, now);
                  motion = {previous, span, now, travelDuration,
                            tfui::cursorMotionTailDuration(pulseInterval)};
                  moving = true;
                }
              }
              auto &point = trailPointFor(lane, span, now);
              point.drawBegin = span.begin;
              point.drawEnd = span.begin + 1;
              // A moving event is represented by its fractional beam. Do not
              // flash the destination ahead of it. Repeats, first events, and
              // non-local jumps still pulse in place.
              if (moving) {
                point.pulsedAt = 0.0;
                point.tailEnergy = 0.f;
              } else {
                pulseTrailPoint(point, now, 0.34f * static_cast<float>(events));
                auto &bloom = bloomSlotFor(lane, now);
                bloom = {span, now,
                         tfui::cursorBloomExpansionDuration(pulseInterval),
                         tfui::cursorBloomTailDuration(pulseInterval),
                         std::min(1.f, 0.88f + 0.04f * static_cast<float>(
                                                           events - 1))};
              }
              observedSpans[lane] = span;
            }
          }
          nvgFontFaceId(args.vg, font->handle);
          nvgFontSize(args.vg, 12.f);
          nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_BASELINE);
          for (auto &point : cursorTrails[lane]) {
            if (!point.span.valid())
              continue;
            const bool current = span.valid() && sameSpan(point.span, span);
            const float persistence = trailPersistence(point, now);
            if (!current && persistence < 0.004f) {
              point = {};
              continue;
            }
            if (point.drawBegin < 0 ||
                point.drawBegin >= static_cast<int>(text.size()))
              continue;
            const float cursorIntensity =
                std::clamp((current ? 0.025f : 0.f) + persistence, 0.f, 1.f);
            const NVGcolor cursorFill = heatmapColor(cursorIntensity);
            // This is trail persistence plus a deliberately faint resting
            // marker. Local movement is drawn by the fractional beam above.
            const int begin = point.drawBegin;
            const int end =
                std::min(point.drawEnd, static_cast<int>(text.size()));
            bndIconLabelCaret(args.vg, textOffset.x, textOffset.y,
                              box.size.x - 2 * textOffset.x,
                              box.size.y - 2 * textOffset.y, -1, invisibleText,
                              12.f, text.c_str(), cursorFill, begin, end);
          }
        }
        // The onset cursor deliberately decays quickly. This quieter layer
        // remains for the complete note/rest/tie span, and its beam makes a
        // slow step's continuing progress explicit between attacks.
        if (sourcePositionsCurrent &&
            module->transportStatus.load(std::memory_order_relaxed) ==
                TfProgSequencer::TransportStatus::Playing) {
          const auto span = unpackSpan(
              module->activeStepSpan.load(std::memory_order_acquire));
          const float progress = tfui::activeStepProgress(
              module->visibleBeat.load(std::memory_order_relaxed),
              module->activeStepBeginBeat.load(std::memory_order_relaxed),
              module->activeStepEndBeat.load(std::memory_order_relaxed));
          drawActiveStep(args, span, progress);
        }
        // Add stationary blooms and moving beams after every resting/history
        // marker. This prevents one lane's marker from masking another lane's
        // glow when their source spans overlap.
        for (auto &blooms : cursorBlooms) {
          for (auto &bloom : blooms)
            drawCursorBloom(args, bloom, now);
        }
        for (auto &motions : cursorMotions) {
          for (auto &motion : motions)
            drawCursorMotion(args, motion, now,
                             tfui::CursorTravelCurve::Linear);
        }
        bndSetFont(APP->window->uiFont->handle);
      }
      nvgResetScissor(args.vg);
    }
    if (layer == 1) {
      nvgScissor(args.vg, RECT_ARGS(args.clipBox));
      auto font = APP->window->loadFont(fontPath);
      if (font && font->handle >= 0) {
        bndSetFont(font->handle);
        NVGcolor invisibleText = color;
        invisibleText.a = 0.f;
        const NVGcolor highlightColor =
            heatmapColor(EditorIntensity::Selection);
        const int begin = std::min(cursor, selection);
        const int end = this == APP->event->selectedWidget
                            ? std::max(cursor, selection)
                            : -1;
        // Draw Rack's normal caret/selection behind the text, then render the
        // glyphs ourselves so comments can be dimmed without touching pixels
        // in the display background.
        bndIconLabelCaret(args.vg, textOffset.x, textOffset.y,
                          box.size.x - 2 * textOffset.x,
                          box.size.y - 2 * textOffset.y, -1, invisibleText,
                          12.f, text.c_str(), highlightColor, begin, end);
        nvgFontFaceId(args.vg, font->handle);
        drawTextWithDimComments(args);
        bndSetFont(APP->window->uiFont->handle);
      }
      nvgResetScissor(args.vg);
    }
    Widget::drawLayer(args, layer);
  }
};

struct TfSequenceStatus : Widget {
  static constexpr float MinimumHeight = 16.f;
  static constexpr float MaximumHeight = 120.f;
  TfProgSequencer *module = nullptr;
  float requiredHeight = MinimumHeight;

  NVGcolor heatmapColor(float intensity) const noexcept {
    const auto palette =
        module ? module->editorHeatmap.load(std::memory_order_relaxed)
               : tfui::HeatmapPalette::Magma;
    return editorColor(palette, intensity);
  }

  std::string statusText() const {
    std::string status =
        module ? module->compileMessage : "PROG SEQUENCER BETA";
    if (!module)
      return status;
    const auto transport =
        module->transportStatus.load(std::memory_order_relaxed);
    const char *name =
        transport == TfProgSequencer::TransportStatus::Playing   ? "PLAY"
        : transport == TfProgSequencer::TransportStatus::Paused  ? "PAUSE"
        : transport == TfProgSequencer::TransportStatus::Stopped ? "STOP"
                                                                 : "WAIT";
    if (status.rfind("READY", 0) == 0)
      status = "READY";
    else if (status.rfind("QUEUED", 0) == 0)
      status = module->pendingProgram.load(std::memory_order_acquire)
                   ? "QUEUED"
                   : "ACTIVE";
    if (module->workspaceOverflow.load(std::memory_order_relaxed))
      status = "INTERNAL ERROR - prepared event workspace exhausted";
    return std::string(name) + " " +
           rack::string::f(
               "%.2f", module->visibleBeat.load(std::memory_order_relaxed)) +
           "  " + status;
  }

  static float wrappedTextHeight(NVGcontext *vg, const std::string &text,
                                 float width) {
    float ascender = 0.f;
    float descender = 0.f;
    float lineHeight = 0.f;
    nvgTextMetrics(vg, &ascender, &descender, &lineHeight);
    if (!(lineHeight > 0.f))
      lineHeight = 10.f;

    const char *remaining = text.c_str();
    const char *end = remaining + text.size();
    std::size_t lineCount = 0;
    do {
      NVGtextRow rows[BND_MAX_ROWS];
      const int rowCount =
          nvgTextBreakLines(vg, remaining, end, width, rows, BND_MAX_ROWS);
      if (rowCount <= 0)
        break;
      lineCount += static_cast<std::size_t>(rowCount);
      const char *next = rows[rowCount - 1].next;
      if (next <= remaining)
        break;
      remaining = next;
    } while (remaining < end);
    return std::max(lineHeight, lineHeight * static_cast<float>(lineCount));
  }

  void drawLayer(const DrawArgs &args, int layer) override {
    if (layer == 1) {
      nvgScissor(args.vg, RECT_ARGS(args.clipBox));
      nvgBeginPath(args.vg);
      nvgRect(args.vg, 0.f, 0.f, box.size.x, box.size.y);
      nvgFillColor(args.vg, heatmapColor(EditorIntensity::Background));
      nvgFill(args.vg);
      auto font = APP->window->loadFont(
          asset::system("res/fonts/ShareTechMono-Regular.ttf"));
      if (font && font->handle >= 0) {
        nvgFontFaceId(args.vg, font->handle);
        nvgFontSize(args.vg, 10.f);
        nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
        nvgTextLineHeight(args.vg, 1.05f);
        nvgFillColor(args.vg, heatmapColor(EditorIntensity::Status));
        const std::string status = statusText();
        const float textWidth = std::max(1.f, box.size.x - 8.f);
        requiredHeight = std::clamp(
            std::ceil(wrappedTextHeight(args.vg, status, textWidth)) + 6.f,
            MinimumHeight, MaximumHeight);
        nvgTextBox(args.vg, 4.f, 3.f, textWidth, status.c_str(), nullptr);
      }
      nvgResetScissor(args.vg);
    }
    Widget::drawLayer(args, layer);
  }
};

struct TfProgSequencerWidget : ModuleWidget {
  static constexpr float EditorTop = 23.f;
  static constexpr float StatusBottom = 375.f;
  static constexpr float DisplayGap = 2.f;
  TfSequenceEditor *editor = nullptr;
  TfSequenceStatus *status = nullptr;
  std::array<PortWidget *, 3> inputs{};
  std::array<PortWidget *, 8> outputs{};
  Widget *runLight = nullptr;
  int displayedWidthHp = 22;

  TfProgSequencerWidget(TfProgSequencer *module) {
    setModule(module);
    setPanel(APP->window->loadSvg(
        asset::plugin(pluginInstance, "res/TfProgSequencer.svg")));

    editor = createWidget<TfSequenceEditor>(Vec(5, EditorTop));
    editor->box.size = Vec(242, StatusBottom - TfSequenceStatus::MinimumHeight -
                                    DisplayGap - EditorTop);
    editor->module = module;
    addChild(editor);

    status = createWidget<TfSequenceStatus>(
        Vec(5, StatusBottom - TfSequenceStatus::MinimumHeight));
    status->box.size = Vec(242, TfSequenceStatus::MinimumHeight);
    status->module = module;
    addChild(status);

    inputs[0] = createInputCentered<PJ301MPort>(Vec(289, 65), module,
                                                TfProgSequencer::CLOCK_INPUT);
    inputs[1] = createInputCentered<PJ301MPort>(Vec(317, 65), module,
                                                TfProgSequencer::RESET_INPUT);
    inputs[2] = createInputCentered<PJ301MPort>(Vec(303, 121), module,
                                                TfProgSequencer::RUN_INPUT);
    for (auto *input : inputs)
      addInput(input);
    outputs[0] = createOutputCentered<PJ301MPort>(
        Vec(289, 188), module, TfProgSequencer::PITCH_OUTPUT);
    outputs[1] = createOutputCentered<PJ301MPort>(Vec(317, 188), module,
                                                  TfProgSequencer::GATE_OUTPUT);
    outputs[2] = createOutputCentered<PJ301MPort>(
        Vec(289, 244), module, TfProgSequencer::TRIGGER_OUTPUT);
    outputs[3] = createOutputCentered<PJ301MPort>(
        Vec(317, 244), module, TfProgSequencer::VELOCITY_OUTPUT);
    outputs[4] = createOutputCentered<PJ301MPort>(
        Vec(289, 300), module, TfProgSequencer::ACCENT_OUTPUT);
    outputs[5] = createOutputCentered<PJ301MPort>(Vec(317, 300), module,
                                                  TfProgSequencer::CV1_OUTPUT);
    outputs[6] = createOutputCentered<PJ301MPort>(Vec(289, 356), module,
                                                  TfProgSequencer::CV2_OUTPUT);
    outputs[7] = createOutputCentered<PJ301MPort>(Vec(317, 356), module,
                                                  TfProgSequencer::CV3_OUTPUT);
    for (auto *output : outputs)
      addOutput(output);
    runLight = createLightCentered<SmallLight<GreenLight>>(
        Vec(284, 121), module, TfProgSequencer::RUN_LIGHT);
    addChild(runLight);
    displayedWidthHp = 0;
    applyPanelWidth(
        module ? module->panelWidthHp.load(std::memory_order_relaxed) : 30);
  }

  void applyPanelWidth(int widthHp) {
    widthHp = validPanelWidth(widthHp);
    if (widthHp == displayedWidthHp)
      return;
    const Vec previousPosition = box.pos;
    const std::string filename = widthHp == 30   ? "TfProgSequencer-30.svg"
                                 : widthHp == 38 ? "TfProgSequencer-38.svg"
                                                 : "TfProgSequencer.svg";
    setPanel(
        APP->window->loadSvg(asset::plugin(pluginInstance, "res/" + filename)));
    displayedWidthHp = widthHp;
    const float right = box.size.x - 27.f;
    const float leftColumn = right - 14.f;
    const float rightColumn = right + 14.f;
    editor->box.size.x = box.size.x - 65.f;
    status->box.size.x = box.size.x - 65.f;
    const Vec inputPositions[] = {Vec(leftColumn, 65.f), Vec(rightColumn, 65.f),
                                  Vec(right, 121.f)};
    for (std::size_t index = 0; index < inputs.size(); ++index)
      inputs[index]->box.pos =
          inputPositions[index].minus(inputs[index]->box.size.div(2));
    const Vec outputPositions[] = {
        Vec(leftColumn, 188.f), Vec(rightColumn, 188.f),
        Vec(leftColumn, 244.f), Vec(rightColumn, 244.f),
        Vec(leftColumn, 300.f), Vec(rightColumn, 300.f),
        Vec(leftColumn, 356.f), Vec(rightColumn, 356.f)};
    for (std::size_t index = 0; index < outputs.size(); ++index)
      outputs[index]->box.pos =
          outputPositions[index].minus(outputs[index]->box.size.div(2));
    runLight->box.pos =
        Vec(right - 19.f, 121.f).minus(runLight->box.size.div(2));
    if (parent && APP && APP->scene && APP->scene->rack)
      APP->scene->rack->setModulePosNearest(this, previousPosition);
  }

  void step() override {
    ModuleWidget::step();
    const float statusHeight =
        status ? status->requiredHeight : TfSequenceStatus::MinimumHeight;
    if (status) {
      status->box.pos.y = StatusBottom - statusHeight;
      status->box.size.y = statusHeight;
    }
    if (editor)
      editor->box.size.y = StatusBottom - statusHeight - DisplayGap - EditorTop;
    auto *prog = getModule<TfProgSequencer>();
    if (prog)
      applyPanelWidth(prog->panelWidthHp.load(std::memory_order_relaxed));
  }

  void appendContextMenu(Menu *menu) override {
    auto *prog = getModule<TfProgSequencer>();
    if (!prog)
      return;
    menu->addChild(new MenuSeparator);
    menu->addChild(createMenuItem("Prog Sequencer documentation", "", []() {
      system::openBrowser(ProgSequencerDocumentationUrl);
    }));
    menu->addChild(new MenuSeparator);
    menu->addChild(createMenuLabel("Editor width"));
    for (const int width : {22, 30, 38}) {
      menu->addChild(createCheckMenuItem(
          rack::string::f("%d HP", width), "",
          [=]() {
            return prog->panelWidthHp.load(std::memory_order_relaxed) == width;
          },
          [=]() {
            prog->panelWidthHp.store(width, std::memory_order_relaxed);
          }));
    }
    menu->addChild(new MenuSeparator);
    menu->addChild(createMenuLabel("Editor heatmap"));
    for (const auto palette :
         {tfui::HeatmapPalette::Magma, tfui::HeatmapPalette::Inferno,
          tfui::HeatmapPalette::Plasma, tfui::HeatmapPalette::Viridis,
          tfui::HeatmapPalette::Cividis, tfui::HeatmapPalette::CrtGreen,
          tfui::HeatmapPalette::CrtBlue, tfui::HeatmapPalette::CrtYellow,
          tfui::HeatmapPalette::CrtRed}) {
      menu->addChild(createCheckMenuItem(
          tfui::heatmapPaletteName(palette), "",
          [=]() {
            return prog->editorHeatmap.load(std::memory_order_relaxed) ==
                   palette;
          },
          [=]() {
            prog->editorHeatmap.store(palette, std::memory_order_relaxed);
          }));
    }
  }
};

Model *modelTfProgSequencer =
    createModel<TfProgSequencer, TfProgSequencerWidget>("TfProgSequencer");
