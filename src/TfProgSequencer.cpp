#include "plugin.hpp"
#include "tfseq.hpp"
#include "tfseq_parser.hpp"

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

constexpr const char *DefaultSource = R"(riff = sequence {
  cycle 8
  tonic C@4
  scale minor
  notes 1 2 3 4 5 6 7 8
  articulation x x _ > [x x] x*3 ~ x
  accent + . .
  duration 1
}

play riff
)";

constexpr int LanguageVersion = 4;

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
    NUM_OUTPUTS
  };
  enum LightIds { RUN_LIGHT, NUM_LIGHTS };

  enum class TransportStatus { Waiting, Playing, Stopped };

  std::string source = DefaultSource;
  std::string evaluatedSource = DefaultSource;
  tfseq::syntax::Document evaluatedDocument;
  std::string compileMessage =
      "READY - Ctrl+. runs all, Ctrl+Enter runs line/selection";
  std::atomic<std::uint64_t> sourceRevision{1};
  std::atomic<tfseq::CompiledProgram *> pendingProgram{nullptr};
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
  std::size_t outputVoiceCount = 1;
  std::size_t scheduledCount = 0;
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
    configLight(RUN_LIGHT, "Running");
    publishSource(source);
  }

  ~TfProgSequencer() override {
    delete pendingProgram.exchange(nullptr);
    delete retiredProgram.exchange(nullptr);
    delete activeProgram;
  }

  void reportDiagnostic(const tfseq::Diagnostic &diagnostic) {
    compileMessage = "ERROR " + std::to_string(diagnostic.line) + ":" +
                     std::to_string(diagnostic.column) + "  " +
                     diagnostic.message;
  }

  bool publishDocument(const tfseq::syntax::Document &document,
                       const std::string &editorSource) {
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
    delete pendingProgram.exchange(candidate, std::memory_order_acq_rel);
    return true;
  }

  bool publishSource(const std::string &text) {
    try {
      const auto parsed = tfseq::syntax::Parse(text);
      if (!parsed) {
        reportDiagnostic(parsed.diagnostic);
        return false;
      }
      return publishDocument(parsed.document, text);
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
      const auto draft = tfseq::syntax::Parse(source);
      if (!draft) {
        reportDiagnostic(draft.diagnostic);
        return false;
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
    if (json_t *sourceJson = json_object_get(root, "source")) {
      if (json_is_string(sourceJson)) {
        source = json_string_value(sourceJson);
        sourceRevision.fetch_add(1, std::memory_order_release);
        // Remove only the constructor/pending candidate before evaluating
        // saved source. A compile failure preserves any already-active valid
        // program; on initial construction there is no active program, so the
        // module remains silent instead of falling back to factory notes.
        delete pendingProgram.exchange(nullptr, std::memory_order_acq_rel);
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
    activationCheckpointValid = false;
    workspaceOverflow.store(false, std::memory_order_relaxed);
    for (auto &span : cursorSpans)
      span.store(0, std::memory_order_relaxed);
    transportStatus.store(TransportStatus::Waiting, std::memory_order_relaxed);
  }

  void swapProgramAtClock(double beat) noexcept {
    if (retiredProgram.load(std::memory_order_acquire) != nullptr)
      return;
    auto *candidate =
        pendingProgram.exchange(nullptr, std::memory_order_acq_rel);
    if (!candidate)
      return;
    auto *previousProgram = activeProgram;
    const bool preservePhase =
        previousProgram && !previousProgram->semantic().stopped &&
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
    target.beat = absoluteBeat;
    if (periodKnown && periodSamples > 0.0 &&
        target.timingOffsetMilliseconds != 0.0)
      target.beat += target.timingOffsetMilliseconds * 0.001 * sampleRateHz /
                     periodSamples;
    auto &schedule = activeProgram->scheduleWorkspace;
    std::size_t child = scheduledCount - 1;
    while (child > 0) {
      const std::size_t parent = (child - 1) / 2;
      if (schedule[parent].beat <= schedule[child].beat)
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
      if (event.cursors[lane].valid())
        cursorSpans[lane].store(packSpan(event.cursors[lane]),
                                std::memory_order_relaxed);
      if (event.cursors[lane].valid())
        cursorPulses[lane].fetch_add(1, std::memory_order_release);
    }
  }

  void applyEvent(const tfseq::RuntimeEvent &event) noexcept {
    if (event.kind == tfseq::EventKind::Tie ||
        event.kind == tfseq::EventKind::Rest || event.voice == 0)
      showCursors(event);
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
    voice.gateOffBeat = event.beat + event.spanBeats * event.gateFraction;
    const double boundary = event.beat + event.spanBeats;
    if (event.legatoToNext)
      voice.gateOffBeat = std::max(voice.gateOffBeat, boundary);
    if (event.kind == tfseq::EventKind::Slide) {
      voice.slideFrom = voice.pitch;
      voice.slideTo = event.pitchVolts;
      voice.slideBeginBeat = event.beat;
      voice.slideEndBeat =
          event.beat + std::min<double>(event.spanBeats, event.slideBeats);
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
              right < scheduledCount &&
                      schedule[right].beat < schedule[left].beat
                  ? right
                  : left;
          if (schedule[parent].beat <= schedule[child].beat)
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

    const bool running = !inputs[RUN_INPUT].isConnected() ||
                         inputs[RUN_INPUT].getVoltage() >= 1.f;
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

  struct CursorTrailPoint {
    tfseq::SourceSpan span;
    int drawBegin = -1;
    int drawEnd = -1;
    float tail = 0.f;
    float head = 0.f;
  };

  TfProgSequencer *module = nullptr;
  std::uint64_t loadedRevision = 0;
  std::array<std::uint64_t, CursorLaneCount> observedPulses{};
  std::array<tfseq::SourceSpan, CursorLaneCount> observedSpans{};
  std::array<std::array<CursorTrailPoint, CursorTrailCapacity>, CursorLaneCount>
      cursorTrails{};
  std::uint64_t observedExecutionPulse = 0;
  float executionTail = 0.f;
  float executionHead = 0.f;
  double lastGlowFrame = 0.0;
  double lastClickTime = -INFINITY;
  Vec lastClickPosition;
  int clickCount = 0;

  TfSequenceEditor() { multiline = true; }

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

  void onChange(const ChangeEvent &) override {
    if (module) {
      module->source = getText();
      // Compiled spans refer to the previous source layout. Let the next clock
      // event establish fresh cursor positions rather than drawing stale ones.
      cursorTrails = {};
      observedSpans = {};
    }
  }

  static bool sameSpan(const tfseq::SourceSpan &left,
                       const tfseq::SourceSpan &right) {
    return left.begin == right.begin && left.end == right.end;
  }

  CursorTrailPoint &trailPointFor(std::size_t lane,
                                  const tfseq::SourceSpan &span) {
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
      if (point.head + point.tail < faintest->head + faintest->tail)
        faintest = &point;
    }
    *faintest = {};
    faintest->span = span;
    faintest->drawBegin = span.begin;
    faintest->drawEnd = span.begin + 1;
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
    NVGcolor commentColor = color;
    commentColor.a *= 0.58f;

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
    app::LedDisplayTextField::onSelectKey(event);
  }

  void drawLayer(const DrawArgs &args, int layer) override {
    if (layer == 1 && module) {
      nvgScissor(args.vg, RECT_ARGS(args.clipBox));
      NVGcolor invisibleText = color;
      invisibleText.a = 0.f;
      const double now = system::getTime();
      const double elapsed = lastGlowFrame > 0.0
                                 ? std::clamp(now - lastGlowFrame, 0.0, 0.25)
                                 : 0.0;
      lastGlowFrame = now;
      // Each event gets an immediate bright head and a longer phosphor-like
      // tail. A small fixed history lets several former positions decay at
      // once without allocating or doing work on the audio thread.
      // A slightly slower head avoids an abrupt brightness step, while the
      // shorter tail clears old positions sooner. The two components together
      // feel smoother without leaving a long-lived haze over dense patterns.
      const float tailDecay = static_cast<float>(std::exp(-elapsed / 0.32));
      const float headDecay = static_cast<float>(std::exp(-elapsed / 0.055));
      executionTail *= tailDecay;
      executionHead *= headDecay;
      const auto executionPulse =
          module->executionPulse.load(std::memory_order_acquire);
      if (executionPulse != observedExecutionPulse) {
        observedExecutionPulse = executionPulse;
        executionHead = 1.f;
        executionTail = std::min(1.f, executionTail + 0.45f);
      }
      auto font = APP->window->loadFont(fontPath);
      if (font && font->handle >= 0) {
        bndSetFont(font->handle);
        const auto executed =
            unpackSpan(module->executionSpan.load(std::memory_order_relaxed));
        if (executed.valid() &&
            executed.begin < static_cast<int>(text.size())) {
          NVGcolor executionFill = color;
          executionFill.a =
              std::min(0.92f, 0.58f * executionHead + 0.30f * executionTail);
          bndIconLabelCaret(
              args.vg, textOffset.x, textOffset.y,
              box.size.x - 2 * textOffset.x, box.size.y - 2 * textOffset.y, -1,
              invisibleText, 12.f, text.c_str(), executionFill, executed.begin,
              std::min(executed.end, static_cast<int>(text.size())));
        }
        for (std::size_t lane = 0; lane < module->cursorSpans.size(); ++lane) {
          for (auto &point : cursorTrails[lane]) {
            point.tail *= tailDecay;
            point.head *= headDecay;
          }
          const auto pulse =
              module->cursorPulses[lane].load(std::memory_order_acquire);
          const auto span = unpackSpan(
              module->cursorSpans[lane].load(std::memory_order_relaxed));
          if (pulse != observedPulses[lane]) {
            const auto events = std::min<std::uint64_t>(
                4, pulse > observedPulses[lane] ? pulse - observedPulses[lane]
                                                : 1);
            observedPulses[lane] = pulse;
            if (span.valid()) {
              const auto previous = observedSpans[lane];
              if (!sameSpan(previous, span) &&
                  isShortSameLineMove(previous, span)) {
                auto &bridge = trailPointFor(lane, previous);
                if (previous.begin < span.begin) {
                  bridge.drawBegin = previous.begin;
                  bridge.drawEnd = span.begin;
                } else {
                  bridge.drawBegin = span.begin + 1;
                  bridge.drawEnd = previous.begin + 1;
                }
                bridge.head = 1.f;
                bridge.tail = std::min(
                    1.f, bridge.tail + 0.34f * static_cast<float>(events));
              }
              auto &point = trailPointFor(lane, span);
              point.drawBegin = span.begin;
              point.drawEnd = span.begin + 1;
              point.head = 1.f;
              point.tail = std::min(
                  1.f, point.tail + 0.34f * static_cast<float>(events));
              observedSpans[lane] = span;
            }
          }
          for (auto &point : cursorTrails[lane]) {
            if (!point.span.valid())
              continue;
            const bool current = span.valid() && sameSpan(point.span, span);
            const float persistence = 0.50f * point.head + 0.28f * point.tail;
            if (!current && persistence < 0.004f) {
              point = {};
              continue;
            }
            if (point.drawBegin < 0 ||
                point.drawBegin >= static_cast<int>(text.size()))
              continue;
            NVGcolor cursorFill = color;
            cursorFill.a =
                std::min(0.92f, (current ? 0.14f : 0.f) + persistence);
            // The live head remains one character wide. A recent short move
            // can extend an older fading sample through the intervening text
            // and whitespace, so the phosphor trail reads as continuous.
            const int begin = point.drawBegin;
            const int end =
                std::min(point.drawEnd, static_cast<int>(text.size()));
            bndIconLabelCaret(args.vg, textOffset.x, textOffset.y,
                              box.size.x - 2 * textOffset.x,
                              box.size.y - 2 * textOffset.y, -1, invisibleText,
                              12.f, text.c_str(), cursorFill, begin, end);
          }
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
        NVGcolor highlightColor = color;
        highlightColor.a = 0.5f;
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
  TfProgSequencer *module = nullptr;

  void drawLayer(const DrawArgs &args, int layer) override {
    if (layer == 1) {
      nvgScissor(args.vg, RECT_ARGS(args.clipBox));
      nvgBeginPath(args.vg);
      nvgRect(args.vg, 0.f, 0.f, box.size.x, box.size.y);
      nvgFillColor(args.vg, nvgRGB(0x08, 0x0a, 0x0b));
      nvgFill(args.vg);
      auto font = APP->window->loadFont(
          asset::system("res/fonts/ShareTechMono-Regular.ttf"));
      if (font && font->handle >= 0) {
        nvgFontFaceId(args.vg, font->handle);
        nvgFontSize(args.vg, 10.f);
        nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
        nvgTextLineHeight(args.vg, 1.05f);
        nvgFillColor(args.vg, nvgRGB(0xc9, 0xd1, 0xd3));
        std::string status = module ? module->compileMessage : "PROG SEQUENCER";
        if (module) {
          const auto transport =
              module->transportStatus.load(std::memory_order_relaxed);
          const char *name =
              transport == TfProgSequencer::TransportStatus::Playing   ? "PLAY"
              : transport == TfProgSequencer::TransportStatus::Stopped ? "STOP"
                                                                       : "WAIT";
          if (status.rfind("QUEUED", 0) == 0 &&
              module->pendingProgram.load(std::memory_order_acquire) == nullptr)
            status = "ACTIVE";
          if (module->workspaceOverflow.load(std::memory_order_relaxed))
            status = "INTERNAL ERROR - prepared event workspace exhausted";
          status =
              std::string(name) + " " +
              string::f("%.2f",
                        module->visibleBeat.load(std::memory_order_relaxed)) +
              "  " + status;
        }
        nvgTextBox(args.vg, 4.f, 3.f, box.size.x - 8.f, status.c_str(),
                   nullptr);
      }
      nvgResetScissor(args.vg);
    }
    Widget::drawLayer(args, layer);
  }
};

struct TfProgSequencerWidget : ModuleWidget {
  TfSequenceEditor *editor = nullptr;
  TfSequenceStatus *status = nullptr;
  std::array<PortWidget *, 3> inputs{};
  std::array<PortWidget *, 5> outputs{};
  Widget *runLight = nullptr;
  int displayedWidthHp = 22;

  TfProgSequencerWidget(TfProgSequencer *module) {
    setModule(module);
    setPanel(APP->window->loadSvg(
        asset::plugin(pluginInstance, "res/TfProgSequencer.svg")));

    editor = createWidget<TfSequenceEditor>(Vec(5, 23));
    editor->box.size = Vec(242, 292);
    editor->module = module;
    addChild(editor);

    status = createWidget<TfSequenceStatus>(Vec(5, 317));
    status->box.size = Vec(242, 58);
    status->module = module;
    addChild(status);

    inputs[0] = createInputCentered<PJ301MPort>(Vec(290, 68), module,
                                                TfProgSequencer::CLOCK_INPUT);
    inputs[1] = createInputCentered<PJ301MPort>(Vec(290, 116), module,
                                                TfProgSequencer::RESET_INPUT);
    inputs[2] = createInputCentered<PJ301MPort>(Vec(290, 164), module,
                                                TfProgSequencer::RUN_INPUT);
    for (auto *input : inputs)
      addInput(input);
    outputs[0] = createOutputCentered<PJ301MPort>(
        Vec(290, 214), module, TfProgSequencer::PITCH_OUTPUT);
    outputs[1] = createOutputCentered<PJ301MPort>(Vec(290, 248), module,
                                                  TfProgSequencer::GATE_OUTPUT);
    outputs[2] = createOutputCentered<PJ301MPort>(
        Vec(290, 282), module, TfProgSequencer::TRIGGER_OUTPUT);
    outputs[3] = createOutputCentered<PJ301MPort>(
        Vec(290, 316), module, TfProgSequencer::VELOCITY_OUTPUT);
    outputs[4] = createOutputCentered<PJ301MPort>(
        Vec(290, 350), module, TfProgSequencer::ACCENT_OUTPUT);
    for (auto *output : outputs)
      addOutput(output);
    runLight = createLightCentered<SmallLight<GreenLight>>(
        Vec(267, 191), module, TfProgSequencer::RUN_LIGHT);
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
    editor->box.size.x = box.size.x - 65.f;
    status->box.size.x = box.size.x - 65.f;
    const float inputY[] = {68.f, 116.f, 164.f};
    for (std::size_t index = 0; index < inputs.size(); ++index)
      inputs[index]->box.pos =
          Vec(right, inputY[index]).minus(inputs[index]->box.size.div(2));
    const float outputY[] = {214.f, 248.f, 282.f, 316.f, 350.f};
    for (std::size_t index = 0; index < outputs.size(); ++index)
      outputs[index]->box.pos =
          Vec(right, outputY[index]).minus(outputs[index]->box.size.div(2));
    runLight->box.pos =
        Vec(right - 23.f, 191.f).minus(runLight->box.size.div(2));
    if (parent && APP && APP->scene && APP->scene->rack)
      APP->scene->rack->setModulePosNearest(this, previousPosition);
  }

  void step() override {
    ModuleWidget::step();
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
          string::f("%d HP", width), "",
          [=]() {
            return prog->panelWidthHp.load(std::memory_order_relaxed) == width;
          },
          [=]() {
            prog->panelWidthHp.store(width, std::memory_order_relaxed);
          }));
    }
  }
};

Model *modelTfProgSequencer =
    createModel<TfProgSequencer, TfProgSequencerWidget>("TfProgSequencer");
