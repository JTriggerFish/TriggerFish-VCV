#include "plugin.hpp"
#include "tfseq.hpp"
#include "tfseq_editor.hpp"
#include "tfseq_parser.hpp"
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

namespace {

namespace EditorIntensity {
constexpr float Background = 0.f;
constexpr float Selection = 0.42f;
constexpr float Comment = 0.60f;
constexpr float Status = 0.76f;
constexpr float Text = 0.86f;
} // namespace EditorIntensity

NVGcolor editorColor(float intensity) noexcept {
  const auto rgb = tfui::sampleHeatmap(tfui::ProgramEditorHeatmap, intensity);
  return nvgRGBf(rgb.red, rgb.green, rgb.blue);
}

constexpr const char *DefaultSource = R"(riff = sequence {
  subdiv 8
  tonic C@4
  scale minor
  notes 1 x2 3{quiet} _ [>4 ^5{stacc}] 6*3 ~ 8{ten}
  offset -6ms 0 +6ms |> rate 1/2
  cv1 0 5 0 |> interp smooth
}

play riff
)";

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

  enum class TransportStatus { Waiting, Playing, Stopped };

  std::string source = DefaultSource;
  std::string evaluatedSource = DefaultSource;
  tfseq::syntax::Document evaluatedDocument;
  std::string compileMessage =
      "READY - Ctrl+. runs all, Ctrl+Enter runs statement/selection";
  std::atomic<std::uint64_t> sourceRevision{1};
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
  std::atomic<int> panelWidthHp{30};
  std::atomic<tfui::CursorTravelCurve> cursorMotionMode{
      tfui::CursorTravelCurve::Linear};
  std::atomic<int> arrangementCursorClocksPerPulse{4};
  std::atomic<bool> editorRunEnabled{true};

  dsp::SchmittTrigger clockTrigger;
  dsp::SchmittTrigger resetTrigger;
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
    float from = 0.f;
    float target = 0.f;
    double beginBeat = 0.0;
    double endBeat = 0.0;
    float power = 1.f;
    bool initialized = false;
    tfseq::CvInterpolation interpolation = tfseq::CvInterpolation::Step;
  };
  std::array<CvOutputState, tfseq::CvLaneCount> cvOutputs{};
  std::size_t outputVoiceCount = 1;
  std::size_t scheduledCount = 0;
  std::uint64_t nextScheduleOrder = 0;
  bool clockSeen = false;
  bool periodKnown = false;
  std::int64_t samplesSinceClock = 0;
  std::int64_t resetIgnoreSamples = 0;
  double periodSamples = 0.0;
  double clockBeat = 0.0;
  double programStartBeat = 0.0;
  double nextStepBeat = 0.0;
  bool wasRunning = true;
  double sampleRateHz = 44100.0;
  std::int64_t lastArrangementCursorGroup =
      std::numeric_limits<std::int64_t>::min();
  std::uint64_t lastArrangementCursorSpan = 0;

  TfProgSequencer() {
    config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
    configInput(CLOCK_INPUT, "Clock");
    configInput(RESET_INPUT, "Reset");
    configInput(RUN_INPUT, "Run gate (normalled high)");
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

  static tfseq::CompiledProgram *pendingProgramPointer(
      std::uintptr_t tagged) noexcept {
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
    std::string nextMessage = "QUEUED - activates on the next clock";
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
    json_object_set_new(root, "cursorMotionMode",
                        json_integer(static_cast<int>(
                            cursorMotionMode.load(std::memory_order_relaxed))));
    json_object_set_new(
        root, "arrangementCursorClocksPerPulse",
        json_integer(arrangementCursorClocksPerPulse.load(
            std::memory_order_relaxed)));
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
    if (json_t *motionJson = json_object_get(root, "cursorMotionMode")) {
      if (json_is_integer(motionJson)) {
        const auto value = json_integer_value(motionJson);
        cursorMotionMode.store(
            value == static_cast<int>(tfui::CursorTravelCurve::Smoothstep)
                ? tfui::CursorTravelCurve::Smoothstep
                : tfui::CursorTravelCurve::Linear,
            std::memory_order_relaxed);
      }
    }
    if (json_t *divisionJson =
            json_object_get(root, "arrangementCursorClocksPerPulse")) {
      if (json_is_integer(divisionJson)) {
        const int value =
            static_cast<int>(json_integer_value(divisionJson));
        arrangementCursorClocksPerPulse.store(
            value == 1 || value == 2 || value == 4 || value == 8 ? value : 4,
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

  void resetTransport() noexcept {
    runtime.reset();
    scheduledCount = 0;
    nextScheduleOrder = 0;
    clockSeen = false;
    periodKnown = false;
    samplesSinceClock = 0;
    clockBeat = 0.0;
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
    activationCheckpointValid = false;
    workspaceOverflow.store(false, std::memory_order_relaxed);
    for (auto &span : cursorSpans)
      span.store(0, std::memory_order_relaxed);
    lastArrangementCursorGroup = std::numeric_limits<std::int64_t>::min();
    lastArrangementCursorSpan = 0;
    transportStatus.store(TransportStatus::Waiting, std::memory_order_relaxed);
  }

  void swapProgramAtClock(double beat) noexcept {
    if (retiredProgram.load(std::memory_order_acquire) != nullptr)
      return;
    const auto pending =
        pendingProgram.exchange(0, std::memory_order_acq_rel);
    auto *candidate = pendingProgramPointer(pending);
    if (!candidate)
      return;
    const bool restartOnActivation = (pending & PendingRestartBit) != 0;
    auto *previousProgram = activeProgram;
    const bool preservePhase =
        !restartOnActivation && previousProgram &&
        !previousProgram->semantic().stopped &&
        !candidate->semantic().stopped && activationCheckpointValid &&
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
    if (activeProgram->semantic().stopped)
      transportStatus.store(TransportStatus::Stopped,
                            std::memory_order_relaxed);
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
      const double millisecondShift =
          target.timingOffsetMilliseconds * 0.001 * sampleRateHz /
          periodSamples;
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
    if (!activeProgram || activeProgram->semantic().stopped ||
        activeProgram->semantic().arrangement.empty())
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
    for (std::size_t lane = 0; lane < cursorSpans.size(); ++lane) {
      if (!event.cursors[lane].valid())
        continue;
      const auto packed = packSpan(event.cursors[lane]);
      cursorSpans[lane].store(packed, std::memory_order_relaxed);
      bool pulse = true;
      if (lane == static_cast<std::size_t>(tfseq::CursorLane::Sequence)) {
        const int clocksPerPulse =
            arrangementCursorClocksPerPulse.load(std::memory_order_relaxed);
        const auto group = tfui::arrangementCursorGroup(
            event.beat - programStartBeat, clocksPerPulse);
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
    voice.gateHigh = event.gateFraction > 0.f;
    double gateDuration =
        event.gateMilliseconds >= 0.f && periodKnown && periodSamples > 0.0
            ? event.gateMilliseconds * 0.001 * sampleRateHz / periodSamples
            : event.spanBeats * event.gateFraction;
    if (event.gateCapMilliseconds >= 0.f && periodKnown &&
        periodSamples > 0.0) {
      gateDuration = std::min(
          gateDuration, event.gateCapMilliseconds * 0.001 * sampleRateHz /
                            periodSamples);
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
          event.beat +
          std::min<double>(
              event.spanBeats,
              event.slideMilliseconds >= 0.f && periodKnown && periodSamples > 0.0
                  ? event.slideMilliseconds * 0.001 * sampleRateHz /
                        periodSamples
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
    sampleRateHz = args.sampleRate;
    if (resetTrigger.process(inputs[RESET_INPUT].getVoltage(), 0.1f, 2.f)) {
      resetTransport();
      resetIgnoreSamples =
          std::max<std::int64_t>(1, std::llround(args.sampleRate * 0.001));
    }
    if (resetIgnoreSamples > 0)
      --resetIgnoreSamples;

    const bool running =
        editorRunEnabled.load(std::memory_order_relaxed) &&
        (!inputs[RUN_INPUT].isConnected() ||
         inputs[RUN_INPUT].getVoltage() >= 1.f);
    bool clockEdge =
        clockTrigger.process(inputs[CLOCK_INPUT].getVoltage(), 0.1f, 2.f);
    if (resetIgnoreSamples > 0)
      clockEdge = false;

    if (!running) {
      if (wasRunning)
        resetTransport();
      for (auto &voice : outputVoices) {
        voice.gateHigh = false;
        voice.sliding = false;
        voice.triggerPulse.reset();
      }
      transportStatus.store(TransportStatus::Stopped,
                            std::memory_order_relaxed);
    } else if (clockEdge) {
      if (clockSeen && samplesSinceClock > 0) {
        const double measured = static_cast<double>(samplesSinceClock);
        periodSamples =
            periodKnown ? 0.75 * periodSamples + 0.25 * measured : measured;
        periodKnown = true;
        clockBeat += 1.0;
      } else {
        clockSeen = true;
        clockBeat = 0.0;
      }
      samplesSinceClock = 0;
      swapProgramAtClock(clockBeat);
      activationCheckpointBeat = clockBeat + 1.0;
      activationCheckpointValid = false;
      if (activeProgram)
        scheduleDueSteps(clockBeat + tfseq::SchedulingLookaheadBeats(
                                         *activeProgram, periodKnown,
                                         periodSamples, sampleRateHz));
      processScheduled(clockBeat);
      transportStatus.store(activeProgram && !activeProgram->semantic().stopped
                                ? TransportStatus::Playing
                            : activeProgram ? TransportStatus::Stopped
                                            : TransportStatus::Waiting,
                            std::memory_order_relaxed);
    }

    double phase = clockBeat;
    if (clockSeen && periodKnown && periodSamples > 0.0) {
      // Integer beat boundaries remain owned by CLOCK. Extrapolation is used
      // only for subdivisions, and cannot advance into the next beat before
      // its external edge actually arrives.
      const double withinBeat = std::min(
          0.999999, static_cast<double>(samplesSinceClock) / periodSamples);
      phase += withinBeat;
    }
    if (running && clockSeen) {
      if (activeProgram)
        scheduleDueSteps(clockBeat + tfseq::SchedulingLookaheadBeats(
                                         *activeProgram, periodKnown,
                                         periodSamples, sampleRateHz));
      processScheduled(phase);
      if (periodKnown && samplesSinceClock > periodSamples * 4.0) {
        for (auto &voice : outputVoices) {
          voice.gateHigh = false;
          voice.sliding = false;
        }
        transportStatus.store(TransportStatus::Waiting,
                              std::memory_order_relaxed);
      }
    }

    bool anyGateHigh = false;
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
      anyGateHigh = anyGateHigh || voice.gateHigh;
    }

    for (auto &state : cvOutputs) {
      if (state.interpolation == tfseq::CvInterpolation::Step ||
          state.endBeat <= state.beginBeat) {
        state.value = state.target;
        continue;
      }
      float amount = static_cast<float>((phase - state.beginBeat) /
                                        (state.endBeat - state.beginBeat));
      amount = std::clamp(amount, 0.f, 1.f);
      if (state.interpolation == tfseq::CvInterpolation::Smooth)
        amount = amount * amount * (3.f - 2.f * amount);
      else if (state.interpolation == tfseq::CvInterpolation::Power)
        amount = std::pow(amount, state.power);
      state.value = state.from + (state.target - state.from) * amount;
    }

    for (const auto output : {PITCH_OUTPUT, GATE_OUTPUT, TRIGGER_OUTPUT,
                              VELOCITY_OUTPUT, ACCENT_OUTPUT})
      outputs[output].setChannels(static_cast<int>(outputVoiceCount));
    for (std::size_t index = 0; index < outputVoiceCount; ++index) {
      auto &voice = outputVoices[index];
      outputs[PITCH_OUTPUT].setVoltage(voice.pitch, static_cast<int>(index));
      outputs[GATE_OUTPUT].setVoltage(voice.gateHigh ? 10.f : 0.f,
                                      static_cast<int>(index));
      outputs[TRIGGER_OUTPUT].setVoltage(
          voice.triggerPulse.process(args.sampleTime) ? 10.f : 0.f,
          static_cast<int>(index));
      outputs[VELOCITY_OUTPUT].setVoltage(voice.velocity,
                                          static_cast<int>(index));
      outputs[ACCENT_OUTPUT].setVoltage(voice.gateHigh ? voice.accent : 0.f,
                                        static_cast<int>(index));
    }
    constexpr std::array<int, tfseq::CvLaneCount> cvOutputIds{
        CV1_OUTPUT, CV2_OUTPUT, CV3_OUTPUT};
    for (std::size_t index = 0; index < cvOutputIds.size(); ++index) {
      outputs[cvOutputIds[index]].setChannels(1);
      outputs[cvOutputIds[index]].setVoltage(cvOutputs[index].value);
    }
    lights[RUN_LIGHT].setBrightness(anyGateHigh ? 1.f : 0.f);
    visibleBeat.store(clockSeen ? phase - programStartBeat : 0.0,
                      std::memory_order_relaxed);
    if (clockSeen)
      ++samplesSinceClock;
    wasRunning = running;
  }
};

struct TfSequenceEditor : app::LedDisplayTextField {
  static constexpr std::size_t CursorLaneCount =
      static_cast<std::size_t>(tfseq::CursorLane::Count);
  static constexpr std::size_t CursorTrailCapacity = 8;
  static constexpr std::size_t CursorMotionCapacity = 4;
  static constexpr std::size_t CursorBloomCapacity = 4;
  static constexpr int CursorGlowSamples = 14;

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
  std::uint64_t observedExecutionPulse = 0;
  double executionPulsedAt = 0.0;
  float executionTailEnergy = 0.f;
  double lastClickTime = -INFINITY;
  Vec lastClickPosition;
  int clickCount = 0;

  TfSequenceEditor() {
    multiline = true;
    color = editorColor(EditorIntensity::Text);
    bgColor = editorColor(EditorIntensity::Background);
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
    module->collectRetired();
    const auto revision =
        module->sourceRevision.load(std::memory_order_acquire);
    if (revision != loadedRevision) {
      setText(module->source);
      loadedRevision = revision;
    }
  }

  void synchronizeEditedSource() {
    if (module) {
      module->source = getText();
      // Compiled spans refer to the previous source layout. Let the next clock
      // event establish fresh cursor positions rather than drawing stale ones.
      cursorTrails = {};
      cursorMotions = {};
      cursorBlooms = {};
      lastCursorPulseTimes = {};
      observedSpans = {};
    }
  }

  void onChange(const ChangeEvent &) override { synchronizeEditedSource(); }

  void applyEdit(tfseq::editor::EditResult edit) {
    if (!edit.changed)
      return;
    text = std::move(edit.text);
    selection = edit.selection;
    cursor = edit.cursor;
    synchronizeEditedSource();
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
      NVGcolor sampleCore = editorColor(sampleIntensity);
      NVGcolor sampleEdge = editorColor(0.12f * sampleIntensity);
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
    NVGcolor headCore = editorColor(intensity);
    NVGcolor headEdge = editorColor(0.14f * intensity);
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
    NVGcolor headFill = editorColor(std::min(1.f, intensity + 0.08f));
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
      NVGcolor layerColor = editorColor(level);
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
    const NVGcolor commentColor = editorColor(EditorIntensity::Comment);

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

  void onButton(const ButtonEvent &event) override {
    app::LedDisplayTextField::onButton(event);
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
  }

  void onSelectKey(const SelectKeyEvent &event) override {
    if (module && event.action == GLFW_PRESS &&
        event.isKeyCommand(GLFW_KEY_PERIOD,
                           RACK_MOD_CTRL | RACK_MOD_SHIFT)) {
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
    if (module && event.action == GLFW_PRESS &&
        event.isKeyCommand(GLFW_KEY_SPACE, RACK_MOD_CTRL)) {
      const bool enabled =
          module->editorRunEnabled.load(std::memory_order_relaxed);
      module->editorRunEnabled.store(!enabled, std::memory_order_relaxed);
      event.consume(this);
      return;
    }
    if (module && event.action == GLFW_PRESS &&
        event.isKeyCommand(GLFW_KEY_D, RACK_MOD_CTRL)) {
      applyEdit(tfseq::editor::Duplicate(text, cursor, selection));
      event.consume(this);
      return;
    }
    app::LedDisplayTextField::onSelectKey(event);
  }

  void drawLayer(const DrawArgs &args, int layer) override {
    if (layer == 1 && module) {
      nvgScissor(args.vg, RECT_ARGS(args.clipBox));
      NVGcolor invisibleText = color;
      invisibleText.a = 0.f;
      const double now = system::getTime();
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
          const NVGcolor executionFill = editorColor(executionIntensity);
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
            const NVGcolor cursorFill = editorColor(cursorIntensity);
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
        // Add stationary blooms and moving beams after every resting/history
        // marker. This prevents one lane's marker from masking another lane's
        // glow when their source spans overlap.
        for (auto &blooms : cursorBlooms) {
          for (auto &bloom : blooms)
            drawCursorBloom(args, bloom, now);
        }
        const auto travelCurve =
            module->cursorMotionMode.load(std::memory_order_relaxed);
        for (auto &motions : cursorMotions) {
          for (auto &motion : motions)
            drawCursorMotion(args, motion, now, travelCurve);
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
        const NVGcolor highlightColor = editorColor(EditorIntensity::Selection);
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

  std::string statusText() const {
    std::string status = module ? module->compileMessage : "PROG SEQUENCER";
    if (!module)
      return status;
    const auto transport =
        module->transportStatus.load(std::memory_order_relaxed);
    const char *name =
        transport == TfProgSequencer::TransportStatus::Playing   ? "PLAY"
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

  void drawLayer(const DrawArgs &args, int layer) override {
    if (layer == 1) {
      nvgScissor(args.vg, RECT_ARGS(args.clipBox));
      nvgBeginPath(args.vg);
      nvgRect(args.vg, 0.f, 0.f, box.size.x, box.size.y);
      nvgFillColor(args.vg, editorColor(EditorIntensity::Background));
      nvgFill(args.vg);
      auto font = APP->window->loadFont(
          asset::system("res/fonts/ShareTechMono-Regular.ttf"));
      if (font && font->handle >= 0) {
        nvgFontFaceId(args.vg, font->handle);
        nvgFontSize(args.vg, 10.f);
        nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
        nvgTextLineHeight(args.vg, 1.05f);
        nvgFillColor(args.vg, editorColor(EditorIntensity::Status));
        const std::string status = statusText();
        float bounds[4]{};
        nvgTextBoxBounds(args.vg, 4.f, 3.f, box.size.x - 8.f, status.c_str(),
                         nullptr, bounds);
        requiredHeight = std::clamp(
            std::ceil(std::max(10.f, bounds[3] - bounds[1])) + 6.f,
            MinimumHeight, MaximumHeight);
        nvgTextBox(args.vg, 4.f, 3.f, box.size.x - 8.f, status.c_str(),
                   nullptr);
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
    editor->box.size =
        Vec(242, StatusBottom - TfSequenceStatus::MinimumHeight - DisplayGap -
                     EditorTop);
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
    outputs[5] = createOutputCentered<PJ301MPort>(
        Vec(317, 300), module, TfProgSequencer::CV1_OUTPUT);
    outputs[6] = createOutputCentered<PJ301MPort>(
        Vec(289, 356), module, TfProgSequencer::CV2_OUTPUT);
    outputs[7] = createOutputCentered<PJ301MPort>(
        Vec(317, 356), module, TfProgSequencer::CV3_OUTPUT);
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
    const Vec inputPositions[] = {Vec(leftColumn, 65.f),
                                  Vec(rightColumn, 65.f), Vec(right, 121.f)};
    for (std::size_t index = 0; index < inputs.size(); ++index)
      inputs[index]->box.pos =
          inputPositions[index].minus(inputs[index]->box.size.div(2));
    const Vec outputPositions[] = {
        Vec(leftColumn, 188.f),  Vec(rightColumn, 188.f),
        Vec(leftColumn, 244.f),  Vec(rightColumn, 244.f),
        Vec(leftColumn, 300.f),  Vec(rightColumn, 300.f),
        Vec(leftColumn, 356.f),  Vec(rightColumn, 356.f)};
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
    const float statusHeight = status ? status->requiredHeight
                                      : TfSequenceStatus::MinimumHeight;
    if (status) {
      status->box.pos.y = StatusBottom - statusHeight;
      status->box.size.y = statusHeight;
    }
    if (editor)
      editor->box.size.y =
          StatusBottom - statusHeight - DisplayGap - EditorTop;
    auto *prog = getModule<TfProgSequencer>();
    if (prog)
      applyPanelWidth(prog->panelWidthHp.load(std::memory_order_relaxed));
  }

  void appendContextMenu(Menu *menu) override {
    auto *prog = getModule<TfProgSequencer>();
    if (!prog)
      return;
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
    menu->addChild(createMenuLabel("Cursor travel"));
    for (const auto mode : {tfui::CursorTravelCurve::Linear,
                            tfui::CursorTravelCurve::Smoothstep}) {
      const char *label =
          mode == tfui::CursorTravelCurve::Linear ? "Linear" : "Smoothstep";
      menu->addChild(createCheckMenuItem(
          label, "",
          [=]() {
            return prog->cursorMotionMode.load(std::memory_order_relaxed) ==
                   mode;
          },
          [=]() {
            prog->cursorMotionMode.store(mode, std::memory_order_relaxed);
          }));
    }
    menu->addChild(new MenuSeparator);
    menu->addChild(createMenuLabel("Arrangement cursor pulse"));
    for (const int clocks : {1, 2, 4, 8}) {
      menu->addChild(createCheckMenuItem(
          clocks == 1 ? "Every clock"
                      : rack::string::f("Every %d clocks", clocks),
          "",
          [=]() {
            return prog->arrangementCursorClocksPerPulse.load(
                       std::memory_order_relaxed) == clocks;
          },
          [=]() {
            prog->arrangementCursorClocksPerPulse.store(
                clocks, std::memory_order_relaxed);
          }));
    }
  }
};

Model *modelTfProgSequencer =
    createModel<TfProgSequencer, TfProgSequencerWidget>("TfProgSequencer");
