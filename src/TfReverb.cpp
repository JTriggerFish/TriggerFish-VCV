#include "plugin.hpp"

#include "components.hpp"
#include "tfdsp/reverb_output.hpp"
#include "tfdsp/room_reverb.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>

namespace {
float Smoothstep(const float value) {
  const float limited = std::clamp(value, 0.f, 1.f);
  return limited * limited * (3.f - 2.f * limited);
}

float InverseSmoothstep(const float targetValue) {
  const float target = std::clamp(targetValue, 0.f, 1.f);
  float lower = 0.f;
  float upper = 1.f;
  for (int iteration = 0; iteration < 20; ++iteration) {
    const float middle = 0.5f * (lower + upper);
    if (Smoothstep(middle) < target)
      lower = middle;
    else
      upper = middle;
  }
  return 0.5f * (lower + upper);
}

struct AspectQuantity : engine::ParamQuantity {
  float getDisplayValue() override {
    return std::exp((2.f * getValue() - 1.f) * std::log(1.8f));
  }
  void setDisplayValue(const float ratio) override {
    setValue(0.5f * (1.f + std::log(std::clamp(ratio, 1.f / 1.8f, 1.8f)) /
                               std::log(1.8f)));
  }
};

struct RoomSizeQuantity : engine::ParamQuantity {
  static constexpr float MinimumFloorScale = 3.1304951685f;
  static constexpr float MaximumFloorScale = 21.2132034356f;
  float getDisplayValue() override {
    return MinimumFloorScale *
           std::exp(Smoothstep(getValue()) *
                    std::log(MaximumFloorScale / MinimumFloorScale));
  }
  void setDisplayValue(const float metres) override {
    const float coordinate =
        std::log(std::clamp(metres, MinimumFloorScale, MaximumFloorScale) /
                 MinimumFloorScale) /
        std::log(MaximumFloorScale / MinimumFloorScale);
    setValue(InverseSmoothstep(coordinate));
  }
};

struct DecayQuantity : engine::ParamQuantity {
  float getDisplayValue() override {
    return 0.25f * std::exp(getValue() * std::log(32.f));
  }
  void setDisplayValue(const float seconds) override {
    setValue(std::log(std::clamp(seconds, 0.25f, 8.f) / 0.25f) /
             std::log(32.f));
  }
};

struct LowCutQuantity : engine::ParamQuantity {
  float getDisplayValue() override {
    return 20.f * std::exp(getValue() * std::log(50.f));
  }
  void setDisplayValue(const float hertz) override {
    setValue(std::log(std::clamp(hertz, 20.f, 1'000.f) / 20.f) /
             std::log(50.f));
  }
};

struct HighCutQuantity : engine::ParamQuantity {
  float getDisplayValue() override {
    return 1'000.f * std::exp(getValue() * std::log(20.f));
  }
  void setDisplayValue(const float hertz) override {
    setValue(std::log(std::clamp(hertz, 1'000.f, 20'000.f) / 1'000.f) /
             std::log(20.f));
  }
};

struct WidthQuantity : engine::ParamQuantity {
  float getDisplayValue() override { return 150.f * Smoothstep(getValue()); }
  void setDisplayValue(const float percent) override {
    setValue(InverseSmoothstep(percent / 150.f));
  }
};
} // namespace

struct TfReverb : Module {
  enum ParamIds {
    SPACE,
    ASPECT,
    LEGACY_HEIGHT,
    PRE_DELAY,
    SOURCE_X,
    SOURCE_Y,
    LISTENER_X,
    LISTENER_Y,
    DECAY,
    DAMPING,
    DIFFUSION,
    MODULATION,
    WIDTH,
    EARLY_LEVEL,
    TAIL_LEVEL,
    LOW_CUT,
    HIGH_CUT,
    MIX,
    LEVEL,
    SHIMMER,
    NUM_PARAMS
  };
  enum InputIds {
    AUDIO_INPUT,
    X_POSITION_INPUT,
    Y_POSITION_INPUT,
    Z_POSITION_INPUT,
    SPACE_CV_INPUT,
    PRE_DELAY_CV_INPUT,
    DECAY_CV_INPUT,
    DAMPING_CV_INPUT,
    MIX_CV_INPUT,
    NUM_INPUTS
  };
  enum OutputIds { LEFT_OUTPUT, RIGHT_OUTPUT, NUM_OUTPUTS };
  enum LightIds { NUM_LIGHTS };

  TfReverb() {
    config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
    configParam<RoomSizeQuantity>(SPACE, 0.f, 1.f,
                                  tfdsp::reverb_defaults::Space,
                                  "Room size", " m");
    configParam<AspectQuantity>(ASPECT, 0.f, 1.f,
                                tfdsp::reverb_defaults::Aspect,
                                "Room width/depth ratio", " W/D");
    // Preserve parameter slot 2 so existing patches keep every subsequent
    // control at the same serialized ID. The value is intentionally ignored:
    // Size now derives a coherent ceiling height for the complete room.
    configParam(LEGACY_HEIGHT, 0.f, 1.f, tfdsp::reverb_defaults::Height,
                "Legacy room height (unused)");
    configParam(PRE_DELAY, 0.f, 1.f, tfdsp::reverb_defaults::PreDelay,
                "Wet pre-delay", " ms", 0.f, 250.f);
    configParam(SOURCE_X, 0.f, 1.f, tfdsp::reverb_defaults::Source[0],
                "Default source X", "%", 0.f, 100.f);
    configParam(SOURCE_Y, 0.f, 1.f, tfdsp::reverb_defaults::Source[1],
                "Default source Y", "%", 0.f, 100.f);
    configParam(LISTENER_X, 0.f, 1.f, tfdsp::reverb_defaults::Listener[0],
                "Listener X", "%", 0.f, 100.f);
    configParam(LISTENER_Y, 0.f, 1.f, tfdsp::reverb_defaults::Listener[1],
                "Listener Y", "%", 0.f, 100.f);
    configParam<DecayQuantity>(DECAY, 0.f, 1.f, tfdsp::reverb_defaults::Decay,
                               "Late decay", " s");
    configParam(DAMPING, 0.f, 1.f, tfdsp::reverb_defaults::Damping,
                "Room damping", "%", 0.f, 100.f);
    configParam(DIFFUSION, 0.f, 1.f, tfdsp::reverb_defaults::Diffusion,
                "Diffusion", "%", 0.f, 100.f);
    configParam(MODULATION, 0.f, 1.f, tfdsp::reverb_defaults::Modulation,
                "Late modulation", "%", 0.f, 100.f);
    configParam<WidthQuantity>(WIDTH, 0.f, 1.f, tfdsp::reverb_defaults::Width,
                               "Stereo width", "%");
    configParam(EARLY_LEVEL, -60.f, 6.f, tfdsp::reverb_defaults::EarlyLevelDb,
                "Early-reflection trim", " dB");
    configParam(TAIL_LEVEL, -60.f, 6.f, tfdsp::reverb_defaults::TailLevelDb,
                "Late-tail trim", " dB");
    configParam<LowCutQuantity>(LOW_CUT, 0.f, 1.f,
                                tfdsp::reverb_defaults::LowCut, "Wet low cut",
                                " Hz");
    configParam<HighCutQuantity>(HIGH_CUT, 0.f, 1.f,
                                 tfdsp::reverb_defaults::HighCut,
                                 "Wet high cut", " Hz");
    configParam(MIX, 0.f, 1.f, tfdsp::reverb_defaults::Mix, "Dry / wet mix",
                "%", 0.f, 100.f);
    configParam(LEVEL, -60.f, 6.f, tfdsp::reverb_defaults::OutputLevelDb,
                "Output level", " dB");
    configParam(SHIMMER, 0.f, 1.f, tfdsp::reverb_defaults::Shimmer,
                "Octave shimmer", "%", 0.f, 100.f);
    getParamQuantity(ASPECT)->description =
        "Reshapes the floor without changing its area. Below centre is "
        "narrower and deeper; above centre is wider and shallower, changing "
        "the reflection pattern while preserving overall room size.";
    getParamQuantity(SPACE)->description =
        "Sets the room's overall dimensions. Larger values give later, more "
        "widely spaced reflections and a slower sense of buildup. Decay sets "
        "how long the tail lasts independently.";
    getParamQuantity(PRE_DELAY)->description =
        "Delays the entire wet response behind the dry signal. Increase it to "
        "separate the source from the room, preserve attack, or create a "
        "rhythmic gap before the reflections begin.";
    getParamQuantity(SOURCE_X)->description =
        "Moves the default source from the left wall to the right wall. Source "
        "position changes early-reflection timing and its relationship to the "
        "listener.";
    getParamQuantity(SOURCE_Y)->description =
        "Moves the default source from the front wall to the back wall. Its "
        "distance from the listener also influences the automatic early/tail "
        "balance.";
    getParamQuantity(LISTENER_X)->description =
        "Moves the listening point from the left wall to the right wall, "
        "changing reflection timing and the stereo perspective.";
    getParamQuantity(LISTENER_Y)->description =
        "Moves the listening point from the front wall to the back wall. "
        "Moving it relative to the source changes perceived distance and the "
        "early/tail balance.";
    getParamQuantity(DECAY)->description =
        "Sets how long the midrange late tail takes to fall by 60 dB. Use "
        "Damping to make low and high frequencies decay at different rates.";
    getParamQuantity(WIDTH)->description =
        "Sets the spread of the wet signal. 0% is mono, 100% keeps the natural "
        "room image, and 150% exaggerates the stereo width.";
    getParamQuantity(MODULATION)->description =
        "Adds slow random movement to the late tail, reducing stationary "
        "ringing and adding animation. The first 35% is static; higher values "
        "progress from subtle motion to audible chorusing.";
    getParamQuantity(DIFFUSION)->description =
        "Controls how quickly individual echoes merge into a smooth tail. Low "
        "values keep more separated, textured reflections; high values create "
        "a denser, softer wash. It does not change the decay time.";
    getParamQuantity(DAMPING)->description =
        "Makes low and high frequencies die away sooner than the midrange. "
        "Increase it for a darker, tighter and more absorbent room; reduce it "
        "for a brighter, more persistent tail.";
    getParamQuantity(EARLY_LEVEL)->description =
        "Adjusts the first discrete wall reflections around their automatic "
        "position-based level. Raise it for more room shape and definition; "
        "lower it for a smoother onset.";
    getParamQuantity(TAIL_LEVEL)->description =
        "Adjusts the diffuse sustain around its automatic position-based "
        "level. Raise it for a larger wash or lower it for a shorter, more "
        "reflection-focused impression.";
    getParamQuantity(SHIMMER)->description =
        "Feeds octave-up energy back through the late tail. Higher values "
        "create a stronger rising harmonic bloom; increase Damping or lower "
        "Wet high cut if the result is too bright.";
    getParamQuantity(LOW_CUT)->description =
        "Removes bass from the wet output without thinning the dry signal. "
        "Raise it to prevent low-frequency mud or leave more space for bass "
        "and kick.";
    getParamQuantity(HIGH_CUT)->description =
        "Softens the wet output above the selected frequency. Lower it for a "
        "darker or more distant room; use Damping when the high frequencies "
        "should also decay faster.";
    getParamQuantity(MIX)->description =
        "Balances the original signal against the complete room response. 0% "
        "is dry and 100% is wet.";
    getParamQuantity(LEVEL)->description =
        "Sets the final output level after the dry/wet mix.";
    getParamQuantity(LOW_CUT)->displayPrecision = 5;
    getParamQuantity(HIGH_CUT)->displayPrecision = 5;
    getParamQuantity(DECAY)->displayPrecision = 4;

    configInput(AUDIO_INPUT, "Mono or positioned polyphonic audio");
    configInput(X_POSITION_INPUT, "Source X positions");
    configInput(Y_POSITION_INPUT, "Source Y positions");
    configInput(Z_POSITION_INPUT, "Source Z positions");
    configInput(SPACE_CV_INPUT, "Room size CV");
    configInput(PRE_DELAY_CV_INPUT, "Pre-delay CV");
    configInput(DECAY_CV_INPUT, "Decay CV");
    configInput(DAMPING_CV_INPUT, "Damping CV");
    configInput(MIX_CV_INPUT, "Dry / wet mix CV");
    configOutput(LEFT_OUTPUT, "Left audio");
    configOutput(RIGHT_OUTPUT, "Right audio");
    configBypass(AUDIO_INPUT, LEFT_OUTPUT);
    configBypass(AUDIO_INPUT, RIGHT_OUTPUT);
    for (auto &position : roomPlanPositions_) {
      position[0].store(tfdsp::reverb_defaults::Source[0],
                        std::memory_order_relaxed);
      position[1].store(tfdsp::reverb_defaults::Source[1],
                        std::memory_order_relaxed);
    }
  }

  void onSampleRateChange(const SampleRateChangeEvent &event) override {
    reverb_.SetSampleRate(event.sampleRate);
    smoothingCoefficient_ = 1.f - std::exp(-1.f / (0.020f * event.sampleRate));
    smoothInitialized_ = false;
    outputGainInitialized_ = false;
    outputGainCountdown_ = 0;
  }

  void onReset() override {
    setLateReverbFlavour(tfdsp::DefaultLateReverbFlavour);
    reverb_.SetLateReverbFlavourImmediate(tfdsp::DefaultLateReverbFlavour);
    reverb_.Reset();
    hasProcessed_ = false;
    smoothInitialized_ = false;
    outputGainInitialized_ = false;
    outputGainCountdown_ = 0;
  }

  void process(const ProcessArgs &args) override {
    const auto flavour = lateReverbFlavour();
    if (!hasProcessed_) {
      reverb_.SetLateReverbFlavourImmediate(flavour);
      hasProcessed_ = true;
    } else {
      reverb_.SetLateReverbFlavour(flavour);
    }
    if (smoothingCoefficient_ <= 0.f)
      smoothingCoefficient_ = 1.f - std::exp(-args.sampleTime / 0.020f);

    std::array<float, NUM_PARAMS> targets{};
    for (int param = 0; param < NUM_PARAMS; ++param)
      targets[static_cast<std::size_t>(param)] = params[param].getValue();
    targets[SPACE] = controlWithCv(SPACE, SPACE_CV_INPUT);
    targets[PRE_DELAY] = controlWithCv(PRE_DELAY, PRE_DELAY_CV_INPUT);
    targets[DECAY] = controlWithCv(DECAY, DECAY_CV_INPUT);
    targets[DAMPING] = controlWithCv(DAMPING, DAMPING_CV_INPUT);
    targets[MIX] = controlWithCv(MIX, MIX_CV_INPUT);
    if (!smoothInitialized_) {
      smoothed_ = targets;
      smoothInitialized_ = true;
    } else {
      for (std::size_t index = 0; index < smoothed_.size(); ++index)
        smoothed_[index] +=
            smoothingCoefficient_ * (targets[index] - smoothed_[index]);
    }

    const std::size_t sourceCount =
        inputs[AUDIO_INPUT].isConnected()
            ? std::min<std::size_t>(tfdsp::RoomReverb::MaximumSources,
                                    inputs[AUDIO_INPUT].getChannels())
            : 0;
    tfdsp::RoomReverb::InputFrame sourceAudio{};
    tfdsp::RoomReverb::SourcePositions positions{};
    float dry = 0.f;
    for (std::size_t source = 0; source < sourceCount; ++source) {
      sourceAudio[source] = inputs[AUDIO_INPUT].getVoltage(source);
      dry += sourceAudio[source];
      positions[source][0] = sourcePosition(
          X_POSITION_INPUT, source, spreadDefaultX(source, sourceCount));
      positions[source][1] =
          sourcePosition(Y_POSITION_INPUT, source, smoothed_[SOURCE_Y]);
      positions[source][2] = sourcePosition(Z_POSITION_INPUT, source,
                                            tfdsp::reverb_defaults::Source[2]);
      roomPlanPositions_[source][0].store(positions[source][0],
                                          std::memory_order_relaxed);
      roomPlanPositions_[source][1].store(positions[source][1],
                                          std::memory_order_relaxed);
    }
    const bool positioned = inputs[X_POSITION_INPUT].isConnected() ||
                            inputs[Y_POSITION_INPUT].isConnected() ||
                            sourceCount > 1;
    roomPlanPositioned_.store(positioned, std::memory_order_relaxed);
    roomPlanSourceCount_.store(sourceCount, std::memory_order_release);

    tfdsp::RoomReverbControls controls;
    controls.space = smoothed_[SPACE];
    controls.aspect = smoothed_[ASPECT];
    controls.listener = {smoothed_[LISTENER_X], smoothed_[LISTENER_Y],
                         tfdsp::reverb_defaults::Listener[2]};
    controls.preDelay = smoothed_[PRE_DELAY];
    controls.decay = smoothed_[DECAY];
    controls.damping = smoothed_[DAMPING];
    controls.diffusion = smoothed_[DIFFUSION];
    controls.modulation = smoothed_[MODULATION];
    controls.shimmer = smoothed_[SHIMMER];
    controls.width = smoothed_[WIDTH];
    controls.earlyLevelDb = smoothed_[EARLY_LEVEL];
    controls.tailLevelDb = smoothed_[TAIL_LEVEL];
    controls.lowCut = smoothed_[LOW_CUT];
    controls.highCut = smoothed_[HIGH_CUT];
    const auto wet =
        reverb_.Process(sourceAudio, positions, sourceCount, controls);

    if (outputGainCountdown_ == 0) {
      const auto target =
          tfdsp::CalculateReverbOutputGains(smoothed_[MIX], smoothed_[LEVEL]);
      if (!outputGainInitialized_) {
        outputGains_ = target;
        outputGainStep_ = {};
        outputGainInitialized_ = true;
      } else {
        outputGainStep_.dry =
            (target.dry - outputGains_.dry) / OutputGainUpdateInterval;
        outputGainStep_.wet =
            (target.wet - outputGains_.wet) / OutputGainUpdateInterval;
      }
      outputGainCountdown_ = OutputGainUpdateInterval;
    }
    const auto mixed = tfdsp::MixReverbOutput(dry, wet, outputGains_);
    outputGains_.dry += outputGainStep_.dry;
    outputGains_.wet += outputGainStep_.wet;
    --outputGainCountdown_;
    outputs[LEFT_OUTPUT].setVoltage(mixed[0]);
    outputs[RIGHT_OUTPUT].setVoltage(mixed[1]);
  }

  json_t *dataToJson() override {
    json_t *root = json_object();
    json_object_set_new(
        root, "lateReverbFlavour",
        json_integer(static_cast<int>(lateReverbFlavour())));
    return root;
  }

  void dataFromJson(json_t *root) override {
    if (json_t *flavourJson = json_object_get(root, "lateReverbFlavour")) {
      if (json_is_integer(flavourJson)) {
        setLateReverbFlavour(
            json_integer_value(flavourJson) ==
                    static_cast<int>(tfdsp::LateReverbFlavour::Optimized)
                ? tfdsp::LateReverbFlavour::Optimized
                : tfdsp::LateReverbFlavour::Base);
      }
    }
  }

  tfdsp::LateReverbFlavour lateReverbFlavour() const noexcept {
    return static_cast<tfdsp::LateReverbFlavour>(
        reverbFlavour_.load(std::memory_order_relaxed));
  }

  void setLateReverbFlavour(const tfdsp::LateReverbFlavour flavour) noexcept {
    reverbFlavour_.store(static_cast<int>(flavour),
                         std::memory_order_relaxed);
  }

  std::size_t roomPlanSourceCount() const noexcept {
    return roomPlanSourceCount_.load(std::memory_order_acquire);
  }

  bool roomPlanPositioned() const noexcept {
    return roomPlanPositioned_.load(std::memory_order_relaxed);
  }

  std::array<float, 2>
  roomPlanSourcePosition(const std::size_t source) const noexcept {
    return {roomPlanPositions_[source][0].load(std::memory_order_relaxed),
            roomPlanPositions_[source][1].load(std::memory_order_relaxed)};
  }

private:
  static constexpr std::size_t OutputGainUpdateInterval = 64;
  tfdsp::RoomReverb reverb_{};
  std::array<std::array<std::atomic<float>, 2>,
             tfdsp::RoomReverb::MaximumSources>
      roomPlanPositions_{};
  std::atomic<std::size_t> roomPlanSourceCount_{0};
  std::atomic<bool> roomPlanPositioned_{false};
  std::atomic<int> reverbFlavour_{
      static_cast<int>(tfdsp::DefaultLateReverbFlavour)};
  std::array<float, NUM_PARAMS> smoothed_{};
  float smoothingCoefficient_{};
  bool smoothInitialized_{};
  bool hasProcessed_{};
  tfdsp::ReverbOutputGains outputGains_{};
  tfdsp::ReverbOutputGains outputGainStep_{};
  std::size_t outputGainCountdown_{};
  bool outputGainInitialized_{};

  float controlWithCv(const int paramId, const int inputId) noexcept {
    const float cv = inputs[inputId].isConnected()
                         ? 0.1f * inputs[inputId].getVoltage(0)
                         : 0.f;
    return std::clamp(params[paramId].getValue() + cv, 0.f, 1.f);
  }

  float sourcePosition(const int inputId, const std::size_t source,
                       const float fallback) noexcept {
    if (!inputs[inputId].isConnected())
      return std::clamp(fallback, 0.f, 1.f);
    return std::clamp(0.1f * inputs[inputId].getPolyVoltage(source), 0.f, 1.f);
  }

  float spreadDefaultX(const std::size_t source,
                       const std::size_t sourceCount) const noexcept {
    if (sourceCount <= 1)
      return smoothed_[SOURCE_X];
    const float coordinate =
        static_cast<float>(source) / static_cast<float>(sourceCount - 1);
    return std::clamp(smoothed_[SOURCE_X] + 0.4f * (coordinate - 0.5f), 0.f,
                      1.f);
  }
};

struct TfRoomPlanWidget : widget::OpaqueWidget {
  enum class Marker { None, Source, Listener };

  struct PlanTooltip : ui::Tooltip {
    TfRoomPlanWidget *plan{};
    void step() override {
      if (plan && plan->module) {
        auto &params = plan->module->params;
        text = string::f(
            "Room plan\nSource: %.1f%%, %.1f%%\nListener: %.1f%%, %.1f%%\n"
            "Distance automatically balances early reflections against the "
            "tail; off-centre placement changes reflection timing and stereo "
            "perspective.\nDrag amber or blue marker; double-click to reset.",
            100.f * params[TfReverb::SOURCE_X].getValue(),
            100.f * params[TfReverb::SOURCE_Y].getValue(),
            100.f * params[TfReverb::LISTENER_X].getValue(),
            100.f * params[TfReverb::LISTENER_Y].getValue());
      } else {
        text = "Room plan\nDrag amber source or blue listener marker.";
      }
      Tooltip::step();
      box.pos = plan->getAbsoluteOffset(Vec(plan->box.size.x, 0.f)).round();
      if (parent)
        box = box.nudge(parent->box.zeroPos());
    }
  };

  TfReverb *module{};
  Marker activeMarker{Marker::None};
  Vec dragPosition{};
  std::array<float, 2> oldValues{};
  ui::Tooltip *tooltip{};

  static constexpr float Margin = 5.f;

  ~TfRoomPlanWidget() override { destroyTooltip(); }

  Vec markerPosition(const int xParam, const int yParam,
                     const Vec fallback) const {
    if (!module)
      return normalizedToScreen(fallback);
    return normalizedToScreen(Vec(module->params[xParam].getValue(),
                                  module->params[yParam].getValue()));
  }

  Vec normalizedToScreen(const Vec normalized) const {
    return Vec(Margin + normalized.x * (box.size.x - 2.f * Margin),
               Margin + (1.f - normalized.y) * (box.size.y - 2.f * Margin));
  }

  Vec screenToNormalized(const Vec screen) const {
    return Vec(
        std::clamp((screen.x - Margin) / (box.size.x - 2.f * Margin), 0.f, 1.f),
        std::clamp(1.f - (screen.y - Margin) / (box.size.y - 2.f * Margin), 0.f,
                   1.f));
  }

  void draw(const DrawArgs &args) override {
    const float left = Margin;
    const float top = Margin;
    const float width = box.size.x - 2.f * Margin;
    const float height = box.size.y - 2.f * Margin;

    nvgBeginPath(args.vg);
    nvgRoundedRect(args.vg, left, top, width, height, 2.f);
    nvgFillColor(args.vg, nvgRGBA(0x1b, 0x20, 0x24, 0xff));
    nvgFill(args.vg);
    nvgStrokeColor(args.vg, nvgRGBA(0x10, 0x10, 0x10, 0xff));
    nvgStrokeWidth(args.vg, 1.5f);
    nvgStroke(args.vg);

    nvgBeginPath(args.vg);
    for (int division = 1; division < 4; ++division) {
      const float fraction = 0.25f * static_cast<float>(division);
      nvgMoveTo(args.vg, left + fraction * width, top);
      nvgLineTo(args.vg, left + fraction * width, top + height);
      nvgMoveTo(args.vg, left, top + fraction * height);
      nvgLineTo(args.vg, left + width, top + fraction * height);
    }
    nvgStrokeColor(args.vg, nvgRGBA(0xff, 0xff, 0xff, 0x18));
    nvgStrokeWidth(args.vg, 0.7f);
    nvgStroke(args.vg);

    // Show live bus/polyphonic positions as smaller, read-only source dots.
    if (module && module->roomPlanPositioned()) {
      const std::size_t sourceCount = module->roomPlanSourceCount();
      for (std::size_t source = 0; source < sourceCount; ++source) {
        const auto position = module->roomPlanSourcePosition(source);
        drawSource(args.vg,
                   normalizedToScreen(Vec(position[0], position[1])), 2.2f,
                   0x78);
      }
    }

    drawSource(args.vg,
               markerPosition(TfReverb::SOURCE_X, TfReverb::SOURCE_Y,
                              Vec(tfdsp::reverb_defaults::Source[0],
                                  tfdsp::reverb_defaults::Source[1])),
               4.f, 0xff);
    drawListener(args.vg,
                 markerPosition(TfReverb::LISTENER_X, TfReverb::LISTENER_Y,
                                Vec(tfdsp::reverb_defaults::Listener[0],
                                    tfdsp::reverb_defaults::Listener[1])));
  }

  void onButton(const ButtonEvent &event) override {
    if (event.action != GLFW_PRESS || event.button != GLFW_MOUSE_BUTTON_LEFT ||
        (event.mods & RACK_MOD_MASK) != 0 || !module)
      return;
    const Vec source =
        markerPosition(TfReverb::SOURCE_X, TfReverb::SOURCE_Y, {});
    const Vec listener =
        markerPosition(TfReverb::LISTENER_X, TfReverb::LISTENER_Y, {});
    const float sourceDistance = event.pos.minus(source).norm();
    const float listenerDistance = event.pos.minus(listener).norm();
    if (std::min(sourceDistance, listenerDistance) > 9.f)
      return;
    activeMarker =
        sourceDistance <= listenerDistance ? Marker::Source : Marker::Listener;
    dragPosition = activeMarker == Marker::Source ? source : listener;
    const auto ids = activeParamIds();
    oldValues = {module->params[ids[0]].getValue(),
                 module->params[ids[1]].getValue()};
    event.consume(this);
  }

  void onDragMove(const DragMoveEvent &event) override {
    if (event.button != GLFW_MOUSE_BUTTON_LEFT || !module ||
        activeMarker == Marker::None)
      return;
    dragPosition = dragPosition.plus(event.mouseDelta.div(getAbsoluteZoom()));
    const Vec position = screenToNormalized(dragPosition);
    const auto ids = activeParamIds();
    module->getParamQuantity(ids[0])->setValue(position.x);
    module->getParamQuantity(ids[1])->setValue(position.y);
  }

  void onDoubleClick(const DoubleClickEvent &event) override {
    if (!module || activeMarker == Marker::None)
      return;
    const auto ids = activeParamIds();
    const std::array<float, 2> defaults =
        activeMarker == Marker::Source
            ? std::array<float, 2>{tfdsp::reverb_defaults::Source[0],
                                   tfdsp::reverb_defaults::Source[1]}
            : std::array<float, 2>{tfdsp::reverb_defaults::Listener[0],
                                   tfdsp::reverb_defaults::Listener[1]};
    module->getParamQuantity(ids[0])->setValue(defaults[0]);
    module->getParamQuantity(ids[1])->setValue(defaults[1]);
  }

  void onDragEnd(const DragEndEvent &event) override {
    if (event.button != GLFW_MOUSE_BUTTON_LEFT || !module ||
        activeMarker == Marker::None)
      return;
    const auto ids = activeParamIds();
    auto *changes = new history::ComplexAction;
    changes->name = activeMarker == Marker::Source ? "move reverb source"
                                                   : "move reverb listener";
    for (std::size_t axis = 0; axis < ids.size(); ++axis) {
      const float newValue = module->params[ids[axis]].getValue();
      if (newValue == oldValues[axis])
        continue;
      auto *change = new history::ParamChange;
      change->name = changes->name;
      change->moduleId = module->id;
      change->paramId = ids[axis];
      change->oldValue = oldValues[axis];
      change->newValue = newValue;
      changes->push(change);
    }
    if (changes->isEmpty())
      delete changes;
    else
      APP->history->push(changes);
    activeMarker = Marker::None;
  }

  void onEnter(const EnterEvent &event) override {
    if (!settings::tooltips || tooltip)
      return;
    auto *created = new PlanTooltip;
    created->plan = this;
    APP->scene->addChild(created);
    tooltip = created;
  }

  void onLeave(const LeaveEvent &event) override { destroyTooltip(); }

private:
  void destroyTooltip() {
    if (!tooltip)
      return;
    APP->scene->removeChild(tooltip);
    delete tooltip;
    tooltip = nullptr;
  }

  std::array<int, 2> activeParamIds() const {
    return activeMarker == Marker::Source
               ? std::array<int, 2>{TfReverb::SOURCE_X, TfReverb::SOURCE_Y}
               : std::array<int, 2>{TfReverb::LISTENER_X, TfReverb::LISTENER_Y};
  }

  static void drawSource(NVGcontext *vg, const Vec position, const float radius,
                         const unsigned char alpha) {
    nvgBeginPath(vg);
    nvgCircle(vg, position.x, position.y, radius);
    nvgFillColor(vg, nvgRGBA(0xff, 0xb0, 0x32, alpha));
    nvgFill(vg);
    nvgStrokeColor(vg, nvgRGBA(0x18, 0x18, 0x18, alpha));
    nvgStrokeWidth(vg, 1.f);
    nvgStroke(vg);
  }

  static void drawListener(NVGcontext *vg, const Vec position) {
    nvgBeginPath(vg);
    nvgCircle(vg, position.x, position.y, 4.6f);
    nvgFillColor(vg, nvgRGBA(0x36, 0xc8, 0xeb, 0xff));
    nvgFill(vg);
    nvgBeginPath(vg);
    nvgCircle(vg, position.x, position.y, 2.f);
    nvgFillColor(vg, nvgRGBA(0x1b, 0x20, 0x24, 0xff));
    nvgFill(vg);
  }
};

struct TfReverbWidget : ModuleWidget {
  TfReverbWidget(TfReverb *module) {
    setModule(module);
    setPanel(APP->window->loadSvg(
        asset::plugin(pluginInstance, "res/TfReverb.svg")));

    addChild(createWidget<ScrewSilver>(Vec(15, 0)));
    addChild(createWidget<ScrewSilver>(Vec(box.size.x - 30, 0)));
    addChild(createWidget<ScrewSilver>(Vec(15, RACK_GRID_HEIGHT - 15)));
    addChild(
        createWidget<ScrewSilver>(Vec(box.size.x - 30, RACK_GRID_HEIGHT - 15)));

    addParam(createParam<TfAudioKob>(Vec(32, 43), module, TfReverb::SPACE));
    addParam(createParam<TfAudioKob>(Vec(132, 43), module, TfReverb::ASPECT));
    addParam(
        createParam<TfAudioKob>(Vec(232, 43), module, TfReverb::PRE_DELAY));

    auto *roomPlan = createWidget<TfRoomPlanWidget>(Vec(54, 101));
    roomPlan->box.size = Vec(192, 61);
    roomPlan->module = module;
    addChild(roomPlan);

    addParam(createParam<TfAudioKob>(Vec(1, 174), module, TfReverb::DECAY));
    addParam(createParam<TfAudioKob>(Vec(51, 174), module, TfReverb::DAMPING));
    addParam(
        createParam<TfAudioKob>(Vec(101, 174), module, TfReverb::DIFFUSION));
    addParam(
        createParam<TfAudioKob>(Vec(151, 174), module, TfReverb::MODULATION));
    addParam(createParam<TfAudioKob>(Vec(201, 174), module, TfReverb::SHIMMER));
    addParam(createParam<TfAudioKob>(Vec(251, 174), module, TfReverb::WIDTH));

    addParam(createParam<TfCvKnob>(Vec(6, 236), module, TfReverb::EARLY_LEVEL));
    addParam(createParam<TfCvKnob>(Vec(56, 236), module, TfReverb::TAIL_LEVEL));
    addParam(createParam<TfCvKnob>(Vec(106, 236), module, TfReverb::LOW_CUT));
    addParam(createParam<TfCvKnob>(Vec(156, 236), module, TfReverb::HIGH_CUT));
    addParam(createParam<TfCvKnob>(Vec(206, 236), module, TfReverb::MIX));
    addParam(createParam<TfCvKnob>(Vec(256, 236), module, TfReverb::LEVEL));

    addInput(
        createInput<PJ301MPort>(Vec(4, 286), module, TfReverb::AUDIO_INPUT));
    addInput(createInput<PJ301MPort>(Vec(37, 286), module,
                                     TfReverb::X_POSITION_INPUT));
    addInput(createInput<PJ301MPort>(Vec(70, 286), module,
                                     TfReverb::Y_POSITION_INPUT));
    addInput(createInput<PJ301MPort>(Vec(103, 286), module,
                                     TfReverb::Z_POSITION_INPUT));
    addInput(createInput<PJ301MPort>(Vec(136, 286), module,
                                     TfReverb::SPACE_CV_INPUT));
    addInput(createInput<PJ301MPort>(Vec(169, 286), module,
                                     TfReverb::PRE_DELAY_CV_INPUT));
    addInput(createInput<PJ301MPort>(Vec(202, 286), module,
                                     TfReverb::DECAY_CV_INPUT));
    addInput(createInput<PJ301MPort>(Vec(235, 286), module,
                                     TfReverb::DAMPING_CV_INPUT));
    addInput(
        createInput<PJ301MPort>(Vec(268, 286), module, TfReverb::MIX_CV_INPUT));

    addOutput(
        createOutput<PJ301MPort>(Vec(101, 335), module, TfReverb::LEFT_OUTPUT));
    addOutput(createOutput<PJ301MPort>(Vec(155, 335), module,
                                       TfReverb::RIGHT_OUTPUT));
  }

  void appendContextMenu(Menu *menu) override {
    auto *reverb = getModule<TfReverb>();
    if (!reverb)
      return;
    menu->addChild(new MenuSeparator);
    menu->addChild(createMenuLabel("Late-tail FDN"));
    for (const auto flavour : {tfdsp::LateReverbFlavour::Base,
                               tfdsp::LateReverbFlavour::Optimized}) {
      const char *label = flavour == tfdsp::LateReverbFlavour::Base
                              ? "Base FDN"
                              : "Optimized FDN";
      menu->addChild(createCheckMenuItem(
          label, "",
          [=]() { return reverb->lateReverbFlavour() == flavour; },
          [=]() { reverb->setLateReverbFlavour(flavour); }));
    }
  }
};

Model *modelTfReverb = createModel<TfReverb, TfReverbWidget>("TfReverb");
