#include "plugin.hpp"

#include "components.hpp"
#include "tfdsp/reverb_output.hpp"
#include "tfdsp/room_reverb.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <string>

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
    BALANCE,
    LOW_CUT,
    HIGH_CUT,
    MIX,
    LEVEL,
    SHIMMER,
    SOURCE_2_X,
    SOURCE_2_Y,
    SOURCE_3_X,
    SOURCE_3_Y,
    SOURCE_4_X,
    SOURCE_4_Y,
    SOURCE_5_X,
    SOURCE_5_Y,
    SOURCE_6_X,
    SOURCE_6_Y,
    SOURCE_7_X,
    SOURCE_7_Y,
    SOURCE_8_X,
    SOURCE_8_Y,
    SOURCE_1_AUTO_X,
    SOURCE_2_AUTO_X,
    SOURCE_3_AUTO_X,
    SOURCE_4_AUTO_X,
    SOURCE_5_AUTO_X,
    SOURCE_6_AUTO_X,
    SOURCE_7_AUTO_X,
    SOURCE_8_AUTO_X,
    NUM_PARAMS
  };
  enum InputIds {
    AUDIO_INPUT,
    LISTENER_X_CV_INPUT,
    DECAY_CV_INPUT,
    DAMPING_CV_INPUT,
    PRE_DELAY_CV_INPUT,
    LISTENER_Y_CV_INPUT,
    BALANCE_CV_INPUT,
    MIX_CV_INPUT,
    NUM_INPUTS
  };
  enum OutputIds { LEFT_OUTPUT, RIGHT_OUTPUT, NUM_OUTPUTS };
  enum LightIds { NUM_LIGHTS };

  inline static constexpr auto SourcePlanDefaults = [] {
    std::array<std::array<float, 2>, tfdsp::RoomReverb::MaximumSources>
        defaults{};
    for (auto &position : defaults)
      position = {tfdsp::reverb_defaults::Source[0],
                  tfdsp::reverb_defaults::Source[1]};
    return defaults;
  }();

  static int sourceXParamId(const std::size_t source) noexcept {
    return source == 0 ? SOURCE_X
                       : SOURCE_2_X + 2 * static_cast<int>(source - 1);
  }

  static int sourceYParamId(const std::size_t source) noexcept {
    return sourceXParamId(source) + 1;
  }

  static int sourceAutoXParamId(const std::size_t source) noexcept {
    return SOURCE_1_AUTO_X + static_cast<int>(source);
  }

  TfReverb() {
    config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
    configParam<RoomSizeQuantity>(
        SPACE, 0.f, 1.f, tfdsp::reverb_defaults::Space, "Room size", " m");
    configParam<AspectQuantity>(ASPECT, 0.f, 1.f,
                                tfdsp::reverb_defaults::Aspect,
                                "Room width/depth ratio", " W/D");
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
    configParam(BALANCE, 0.f, 1.f, tfdsp::reverb_defaults::Balance,
                "Early / late balance", "%", 0.f, 200.f, -100.f);
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
    for (std::size_t source = 1; source < SourcePlanDefaults.size(); ++source) {
      const std::string name = "Source " + std::to_string(source + 1);
      configParam(sourceXParamId(source), 0.f, 1.f,
                  SourcePlanDefaults[source][0], name + " X", "%", 0.f, 100.f);
      configParam(sourceYParamId(source), 0.f, 1.f,
                  SourcePlanDefaults[source][1], name + " Y", "%", 0.f, 100.f);
      getParamQuantity(sourceXParamId(source))->description =
          "Sets this source's left/right room position.";
      getParamQuantity(sourceYParamId(source))->description =
          "Sets this source's front/back room position.";
    }
    for (std::size_t source = 0; source < SourcePlanDefaults.size(); ++source) {
      const std::string name = "Source " + std::to_string(source + 1);
      configSwitch(sourceAutoXParamId(source), 0.f, 1.f, 1.f,
                   name + " horizontal placement", {"Manual", "Automatic"});
    }
    getParamQuantity(ASPECT)->description =
        "Changes room width and depth while preserving floor area.";
    getParamQuantity(SPACE)->description =
        "Sets room dimensions and reflection spacing.";
    getParamQuantity(PRE_DELAY)->description =
        "Delays the wet response behind the dry attack.";
    getParamQuantity(SOURCE_X)->description =
        "Sets default source left/right position and direct panning.";
    getParamQuantity(SOURCE_Y)->description =
        "Sets default source front/back position and distance.";
    getParamQuantity(LISTENER_X)->description =
        "Sets listener left/right position and stereo perspective.";
    getParamQuantity(LISTENER_Y)->description =
        "Sets listener front/back position and source distance.";
    getParamQuantity(DECAY)->description =
        "Sets the midrange late-tail RT60.";
    getParamQuantity(WIDTH)->description =
        "Sets wet width: 0% mono, 100% natural, 150% extra-wide.";
    getParamQuantity(MODULATION)->description =
        "Adds slow late-tail movement, from subtle motion to chorusing.";
    getParamQuantity(DIFFUSION)->description =
        "Sets how quickly echoes merge into a dense, smooth tail.";
    getParamQuantity(DAMPING)->description =
        "Shortens low- and high-frequency decay relative to the midrange.";
    getParamQuantity(BALANCE)->description =
        "Turns from early reflections through the inferred mix to late tail.";
    getParamQuantity(SHIMMER)->description =
        "Feeds octave-up energy through the late tail.";
    getParamQuantity(LOW_CUT)->description =
        "Removes bass from the wet signal only.";
    getParamQuantity(HIGH_CUT)->description =
        "Softens high frequencies in the wet signal only.";
    getParamQuantity(MIX)->description =
        "Balances dry input and complete room response.";
    getParamQuantity(LEVEL)->description =
        "Sets output level after the dry/wet mix.";
    getParamQuantity(LOW_CUT)->displayPrecision = 5;
    getParamQuantity(HIGH_CUT)->displayPrecision = 5;
    getParamQuantity(DECAY)->displayPrecision = 4;

    configInput(AUDIO_INPUT, "Mono or polyphonic source audio");
    configInput(LISTENER_X_CV_INPUT, "Listener X slow CV");
    configInput(DECAY_CV_INPUT, "Decay CV");
    configInput(DAMPING_CV_INPUT, "Damping CV");
    configInput(PRE_DELAY_CV_INPUT, "Pre-delay CV");
    configInput(LISTENER_Y_CV_INPUT, "Listener Y slow CV");
    configInput(BALANCE_CV_INPUT, "Early / late balance CV");
    configInput(MIX_CV_INPUT, "Dry / wet mix CV");
    configOutput(LEFT_OUTPUT, "Left audio");
    configOutput(RIGHT_OUTPUT, "Right audio");
    configBypass(AUDIO_INPUT, LEFT_OUTPUT);
    configBypass(AUDIO_INPUT, RIGHT_OUTPUT);
    roomPlanListenerPosition_[0].store(tfdsp::reverb_defaults::Listener[0],
                                       std::memory_order_relaxed);
    roomPlanListenerPosition_[1].store(tfdsp::reverb_defaults::Listener[1],
                                       std::memory_order_relaxed);
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

    const std::size_t sourceCount =
        inputs[AUDIO_INPUT].isConnected()
            ? std::min<std::size_t>(tfdsp::RoomReverb::MaximumSources,
                                    inputs[AUDIO_INPUT].getChannels())
            : 0;

    std::array<float, NUM_PARAMS> targets{};
    for (int param = 0; param < NUM_PARAMS; ++param)
      targets[static_cast<std::size_t>(param)] = params[param].getValue();
    targets[LISTENER_X] = controlWithCv(LISTENER_X, LISTENER_X_CV_INPUT);
    targets[LISTENER_Y] = controlWithCv(LISTENER_Y, LISTENER_Y_CV_INPUT);
    for (std::size_t axis = 0; axis < 2; ++axis) {
      const int paramId = axis == 0 ? LISTENER_X : LISTENER_Y;
      roomPlanListenerPosition_[axis].store(
          targets[static_cast<std::size_t>(paramId)],
          std::memory_order_relaxed);
      roomPlanListenerCvOffset_[axis].store(
          targets[static_cast<std::size_t>(paramId)] -
              params[paramId].getValue(),
          std::memory_order_relaxed);
    }
    targets[PRE_DELAY] = controlWithCv(PRE_DELAY, PRE_DELAY_CV_INPUT);
    targets[DECAY] = controlWithCv(DECAY, DECAY_CV_INPUT);
    targets[DAMPING] = controlWithCv(DAMPING, DAMPING_CV_INPUT);
    targets[BALANCE] = controlWithCv(BALANCE, BALANCE_CV_INPUT);
    targets[MIX] = controlWithCv(MIX, MIX_CV_INPUT);
    for (std::size_t source = 0; source < sourceCount; ++source) {
      const auto xParam = static_cast<std::size_t>(sourceXParamId(source));
      if (sourceXAutomatic(source))
        targets[xParam] =
            tfdsp::reverb_defaults::ProgressiveSourceX(source, sourceCount);
    }
    if (!smoothInitialized_) {
      smoothed_ = targets;
      smoothInitialized_ = true;
    } else {
      for (std::size_t index = 0; index < smoothed_.size(); ++index)
        smoothed_[index] +=
            smoothingCoefficient_ * (targets[index] - smoothed_[index]);
    }

    tfdsp::RoomReverb::InputFrame sourceAudio{};
    tfdsp::RoomReverb::SourcePositions positions{};
    for (std::size_t source = 0; source < sourceCount; ++source) {
      sourceAudio[source] = inputs[AUDIO_INPUT].getVoltage(source);
      positions[source][0] =
          smoothed_[static_cast<std::size_t>(sourceXParamId(source))];
      positions[source][1] =
          smoothed_[static_cast<std::size_t>(sourceYParamId(source))];
      positions[source][2] = tfdsp::reverb_defaults::Source[2];
    }
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
    controls.balance = smoothed_[BALANCE];
    controls.lowCut = smoothed_[LOW_CUT];
    controls.highCut = smoothed_[HIGH_CUT];
    const auto frame =
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
    const auto mixed =
        tfdsp::MixReverbOutput(frame.direct, frame.wet, outputGains_);
    outputGains_.dry += outputGainStep_.dry;
    outputGains_.wet += outputGainStep_.wet;
    --outputGainCountdown_;
    outputs[LEFT_OUTPUT].setVoltage(mixed[0]);
    outputs[RIGHT_OUTPUT].setVoltage(mixed[1]);
  }

  json_t *dataToJson() override {
    json_t *root = json_object();
    json_object_set_new(root, "lateReverbFlavour",
                        json_integer(static_cast<int>(lateReverbFlavour())));
    json_object_set_new(root, "sourceXPlacementVersion", json_integer(1));
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

  void fromJson(json_t *root) override {
    Module::fromJson(root);
    const json_t *data = json_object_get(root, "data");
    const json_t *placementVersion =
        json_is_object(data) ? json_object_get(data, "sourceXPlacementVersion")
                             : nullptr;
    if (!json_is_integer(placementVersion) ||
        json_integer_value(placementVersion) < 1) {
      // Patches saved before placement mode became explicit used X=0.5 as the
      // automatic-spread sentinel. Preserve their effective positions once.
      for (std::size_t source = 0; source < SourcePlanDefaults.size();
           ++source) {
        const float x = params[sourceXParamId(source)].getValue();
        params[sourceAutoXParamId(source)].setValue(
            std::abs(x - SourcePlanDefaults[source][0]) < 1.e-6f ? 1.f : 0.f);
      }
    }
  }

  tfdsp::LateReverbFlavour lateReverbFlavour() const noexcept {
    return static_cast<tfdsp::LateReverbFlavour>(
        reverbFlavour_.load(std::memory_order_relaxed));
  }

  void setLateReverbFlavour(const tfdsp::LateReverbFlavour flavour) noexcept {
    reverbFlavour_.store(static_cast<int>(flavour), std::memory_order_relaxed);
  }

  std::size_t roomPlanSourceCount() const noexcept {
    return roomPlanSourceCount_.load(std::memory_order_acquire);
  }

  std::array<float, 2>
  roomPlanSourcePosition(const std::size_t source) noexcept {
    const float configuredX = params[sourceXParamId(source)].getValue();
    const float x = sourceXAutomatic(source)
                        ? tfdsp::reverb_defaults::ProgressiveSourceX(
                              source, roomPlanSourceCount())
                        : configuredX;
    return {x, params[sourceYParamId(source)].getValue()};
  }

  bool sourceXAutomatic(const std::size_t source) noexcept {
    return params[sourceAutoXParamId(source)].getValue() >= 0.5f;
  }

  std::array<float, 2> roomPlanListenerPosition() const noexcept {
    return {roomPlanListenerPosition_[0].load(std::memory_order_relaxed),
            roomPlanListenerPosition_[1].load(std::memory_order_relaxed)};
  }

  float roomPlanListenerCvOffset(const std::size_t axis) const noexcept {
    return roomPlanListenerCvOffset_[axis].load(std::memory_order_relaxed);
  }

  void applyPreset(const tfdsp::reverb_defaults::ReverbPreset &preset,
                   const char *presetName) {
    auto *changes = new history::ComplexAction;
    changes->name = std::string("set reverb preset: ") + presetName;
    const auto setParameter = [&](const int id, const float value) {
      const float oldValue = params[id].getValue();
      getParamQuantity(id)->setValue(value);
      const float newValue = params[id].getValue();
      if (oldValue == newValue)
        return;
      auto *change = new history::ParamChange;
      change->name = changes->name;
      change->moduleId = this->id;
      change->paramId = id;
      change->oldValue = oldValue;
      change->newValue = newValue;
      changes->push(change);
    };
    setParameter(SPACE, preset.space);
    setParameter(ASPECT, preset.aspect);
    setParameter(PRE_DELAY, preset.preDelay);
    setParameter(SOURCE_X, preset.source[0]);
    setParameter(SOURCE_Y, preset.source[1]);
    // A preset is a complete room placement, not only a set of visible knob
    // values. Return every channel to automatic horizontal spreading and give
    // every source the preset's calibrated depth.
    for (std::size_t source = 1; source < SourcePlanDefaults.size(); ++source) {
      setParameter(sourceXParamId(source), SourcePlanDefaults[source][0]);
      setParameter(sourceYParamId(source), preset.source[1]);
    }
    for (std::size_t source = 0; source < SourcePlanDefaults.size(); ++source)
      setParameter(sourceAutoXParamId(source), 1.f);
    setParameter(LISTENER_X, preset.listener[0]);
    setParameter(LISTENER_Y, preset.listener[1]);
    setParameter(DECAY, preset.decay);
    setParameter(DAMPING, preset.damping);
    setParameter(DIFFUSION, preset.diffusion);
    setParameter(MODULATION, preset.modulation);
    setParameter(SHIMMER, preset.shimmer);
    setParameter(WIDTH, preset.width);
    setParameter(BALANCE, preset.balance);
    setParameter(LOW_CUT, preset.lowCut);
    setParameter(HIGH_CUT, preset.highCut);
    setParameter(MIX, preset.mix);
    if (changes->isEmpty())
      delete changes;
    else
      APP->history->push(changes);
  }

private:
  static constexpr std::size_t OutputGainUpdateInterval = 64;
  tfdsp::RoomReverb reverb_{};
  std::atomic<std::size_t> roomPlanSourceCount_{0};
  std::array<std::atomic<float>, 2> roomPlanListenerPosition_{};
  std::array<std::atomic<float>, 2> roomPlanListenerCvOffset_{};
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
};

struct TfRoomPlanWidget : widget::OpaqueWidget {
  enum class Marker { None, Source, Listener };

  struct PlanTooltip : ui::Tooltip {
    TfRoomPlanWidget *plan{};
    void step() override {
      if (plan && plan->module) {
        text = rack::string::f(
            "Room plan - %zu source%s\n"
            "Drag markers; double-click to reset.\n"
            "Affects direct sound and early reflections.",
            plan->module->roomPlanSourceCount(),
            plan->module->roomPlanSourceCount() == 1 ? "" : "s");
      } else {
        text = "Room plan\nDrag a source or the listener marker.";
      }
      Tooltip::step();
      if (plan)
        box.pos =
            plan->getAbsoluteOffset(Vec(plan->box.size.x, 0.f)).round();
      if (parent)
        box = box.nudge(parent->box.zeroPos());
    }
  };

  TfReverb *module{};
  Marker activeMarker{Marker::None};
  std::size_t activeSource{};
  Vec dragPosition{};
  std::array<float, 2> oldValues{};
  float oldAutomaticX{1.f};
  ui::Tooltip *tooltip{};

  static constexpr float Margin = 5.f;

  ~TfRoomPlanWidget() override { destroyTooltip(); }

  Vec listenerMarkerPosition() const {
    if (!module)
      return normalizedToScreen(Vec(tfdsp::reverb_defaults::Listener[0],
                                    tfdsp::reverb_defaults::Listener[1]));
    const auto position = module->roomPlanListenerPosition();
    return normalizedToScreen(Vec(position[0], position[1]));
  }

  Vec normalizedToScreen(const Vec normalized) const {
    return Vec(Margin + normalized.x * (box.size.x - 2.f * Margin),
               Margin + normalized.y * (box.size.y - 2.f * Margin));
  }

  Vec screenToNormalized(const Vec screen) const {
    return Vec(
        std::clamp((screen.x - Margin) / (box.size.x - 2.f * Margin), 0.f, 1.f),
        std::clamp((screen.y - Margin) / (box.size.y - 2.f * Margin), 0.f,
                   1.f));
  }

  Vec sourceMarkerPosition(const std::size_t source) const {
    if (!module)
      return normalizedToScreen(Vec(TfReverb::SourcePlanDefaults[source][0],
                                    TfReverb::SourcePlanDefaults[source][1]));
    const auto position = module->roomPlanSourcePosition(source);
    return normalizedToScreen(Vec(position[0], position[1]));
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

    nvgFontFaceId(args.vg, APP->window->uiFont->handle);
    nvgFontSize(args.vg, 5.2f);
    nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
    nvgFillColor(args.vg, nvgRGBA(0xff, 0xb0, 0x32, 0xff));
    nvgText(args.vg, left + 5.f, top + 4.f, "SOURCES", nullptr);
    nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_TOP);
    nvgFillColor(args.vg, nvgRGBA(0x36, 0xc8, 0xeb, 0xff));
    nvgText(args.vg, left + width - 5.f, top + 4.f, "LISTENER", nullptr);

    if (module) {
      const std::size_t sourceCount = module->roomPlanSourceCount();
      for (std::size_t source = 0; source < sourceCount; ++source) {
        drawSource(args.vg, sourceMarkerPosition(source), source);
      }
    }
    drawListener(args.vg, listenerMarkerPosition());
  }

  void onButton(const ButtonEvent &event) override {
    if (event.action != GLFW_PRESS || event.button != GLFW_MOUSE_BUTTON_LEFT ||
        (event.mods & RACK_MOD_MASK) != 0 || !module)
      return;
    const Vec listener = listenerMarkerPosition();
    float nearestDistance = event.pos.minus(listener).norm();
    activeMarker = Marker::Listener;
    dragPosition = listener;
    const std::size_t sourceCount = module->roomPlanSourceCount();
    for (std::size_t source = 0; source < sourceCount; ++source) {
      const Vec position = sourceMarkerPosition(source);
      const float distance = event.pos.minus(position).norm();
      if (distance >= nearestDistance)
        continue;
      nearestDistance = distance;
      activeMarker = Marker::Source;
      activeSource = source;
      dragPosition = position;
    }
    if (nearestDistance > 9.f) {
      activeMarker = Marker::None;
      return;
    }
    const auto ids = activeParamIds();
    oldValues = {module->params[ids[0]].getValue(),
                 module->params[ids[1]].getValue()};
    if (activeMarker == Marker::Source) {
      oldAutomaticX =
          module->params[TfReverb::sourceAutoXParamId(activeSource)].getValue();
    }
    event.consume(this);
  }

  void onDragMove(const DragMoveEvent &event) override {
    if (event.button != GLFW_MOUSE_BUTTON_LEFT || !module ||
        activeMarker == Marker::None)
      return;
    dragPosition = dragPosition.plus(event.mouseDelta.div(getAbsoluteZoom()));
    const Vec position = screenToNormalized(dragPosition);
    const auto ids = activeParamIds();
    if (activeMarker == Marker::Source) {
      module->getParamQuantity(TfReverb::sourceAutoXParamId(activeSource))
          ->setValue(0.f);
    }
    const std::array<float, 2> coordinates{position.x, position.y};
    for (std::size_t axis = 0; axis < ids.size(); ++axis) {
      const float cvOffset = activeMarker == Marker::Listener
                                 ? module->roomPlanListenerCvOffset(axis)
                                 : 0.f;
      module->getParamQuantity(ids[axis])->setValue(coordinates[axis] -
                                                    cvOffset);
    }
  }

  void onDoubleClick(const DoubleClickEvent &event) override {
    if (!module || activeMarker == Marker::None)
      return;
    const auto ids = activeParamIds();
    const std::array<float, 2> defaults =
        activeMarker == Marker::Source
            ? TfReverb::SourcePlanDefaults[activeSource]
            : std::array<float, 2>{tfdsp::reverb_defaults::Listener[0],
                                   tfdsp::reverb_defaults::Listener[1]};
    module->getParamQuantity(ids[0])->setValue(defaults[0]);
    module->getParamQuantity(ids[1])->setValue(defaults[1]);
    if (activeMarker == Marker::Source) {
      module->getParamQuantity(TfReverb::sourceAutoXParamId(activeSource))
          ->setValue(1.f);
    }
    commitPositionHistory();
    activeMarker = Marker::None;
  }

  void onDragEnd(const DragEndEvent &event) override {
    if (event.button != GLFW_MOUSE_BUTTON_LEFT || !module ||
        activeMarker == Marker::None)
      return;
    commitPositionHistory();
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
  void commitPositionHistory() {
    if (!module || activeMarker == Marker::None)
      return;
    const auto ids = activeParamIds();
    auto *changes = new history::ComplexAction;
    changes->name =
        activeMarker == Marker::Source
            ? rack::string::f("move reverb source %zu", activeSource + 1)
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
    if (activeMarker == Marker::Source) {
      const int automaticId = TfReverb::sourceAutoXParamId(activeSource);
      const float newAutomaticX = module->params[automaticId].getValue();
      if (newAutomaticX != oldAutomaticX) {
        auto *change = new history::ParamChange;
        change->name = changes->name;
        change->moduleId = module->id;
        change->paramId = automaticId;
        change->oldValue = oldAutomaticX;
        change->newValue = newAutomaticX;
        changes->push(change);
      }
    }
    if (changes->isEmpty())
      delete changes;
    else
      APP->history->push(changes);
  }

  void destroyTooltip() {
    if (!tooltip)
      return;
    APP->scene->removeChild(tooltip);
    delete tooltip;
    tooltip = nullptr;
  }

  std::array<int, 2> activeParamIds() const {
    return activeMarker == Marker::Source
               ? std::array<int, 2>{TfReverb::sourceXParamId(activeSource),
                                    TfReverb::sourceYParamId(activeSource)}
               : std::array<int, 2>{TfReverb::LISTENER_X, TfReverb::LISTENER_Y};
  }

  static void drawSource(NVGcontext *vg, const Vec position,
                         const std::size_t source) {
    nvgBeginPath(vg);
    nvgCircle(vg, position.x, position.y, 5.2f);
    nvgFillColor(vg, nvgRGBA(0xff, 0xb0, 0x32, 0xff));
    nvgFill(vg);
    nvgStrokeColor(vg, nvgRGBA(0x18, 0x18, 0x18, 0xff));
    nvgStrokeWidth(vg, 1.f);
    nvgStroke(vg);
    nvgFontFaceId(vg, APP->window->uiFont->handle);
    nvgFontSize(vg, 5.8f);
    nvgTextAlign(vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
    nvgFillColor(vg, nvgRGBA(0x18, 0x18, 0x18, 0xff));
    const std::string label = std::to_string(source + 1);
    nvgText(vg, position.x, position.y + 0.2f, label.c_str(), nullptr);
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

    auto *roomPlan = createWidget<TfRoomPlanWidget>(Vec(17, 27));
    roomPlan->box.size = Vec(206, 92);
    roomPlan->module = module;
    addChild(roomPlan);

    addParam(createParam<TfAudioKob>(Vec(6, 130), module, TfReverb::SPACE));
    addParam(createParam<TfAudioKob>(Vec(54, 130), module, TfReverb::ASPECT));
    addParam(
        createParam<TfAudioKob>(Vec(102, 130), module, TfReverb::PRE_DELAY));
    addParam(createParam<TfAudioKob>(Vec(150, 130), module, TfReverb::DECAY));
    addParam(createParam<TfAudioKob>(Vec(198, 130), module, TfReverb::DAMPING));

    addParam(createParam<TfAudioKob>(Vec(6, 184), module, TfReverb::DIFFUSION));
    addParam(
        createParam<TfAudioKob>(Vec(54, 184), module, TfReverb::MODULATION));
    addParam(createParam<TfAudioKob>(Vec(102, 184), module, TfReverb::SHIMMER));
    addParam(createParam<TfAudioKob>(Vec(150, 184), module, TfReverb::WIDTH));
    addParam(createParam<TfAudioKob>(Vec(198, 184), module, TfReverb::BALANCE));

    addParam(
        createParam<TfCvKnob>(Vec(33.826, 242), module, TfReverb::LOW_CUT));
    addParam(
        createParam<TfCvKnob>(Vec(81.826, 242), module, TfReverb::HIGH_CUT));
    addParam(createParam<TfCvKnob>(Vec(129.826, 242), module, TfReverb::MIX));
    addParam(createParam<TfCvKnob>(Vec(177.826, 242), module, TfReverb::LEVEL));

    addInput(createInput<PJ301MPort>(Vec(12.15, 292), module,
                                     TfReverb::AUDIO_INPUT));
    addInput(createInput<PJ301MPort>(Vec(60.15, 292), module,
                                     TfReverb::LISTENER_X_CV_INPUT));
    addInput(createInput<PJ301MPort>(Vec(108.15, 292), module,
                                     TfReverb::DECAY_CV_INPUT));
    addInput(createInput<PJ301MPort>(Vec(156.15, 292), module,
                                     TfReverb::DAMPING_CV_INPUT));
    addInput(createInput<PJ301MPort>(Vec(12.15, 339), module,
                                     TfReverb::PRE_DELAY_CV_INPUT));
    addInput(createInput<PJ301MPort>(Vec(60.15, 339), module,
                                     TfReverb::LISTENER_Y_CV_INPUT));
    addInput(createInput<PJ301MPort>(Vec(108.15, 339), module,
                                     TfReverb::BALANCE_CV_INPUT));
    addInput(createInput<PJ301MPort>(Vec(156.15, 339), module,
                                     TfReverb::MIX_CV_INPUT));

    addOutput(createOutput<PJ301MPort>(Vec(204.15, 292), module,
                                       TfReverb::LEFT_OUTPUT));
    addOutput(createOutput<PJ301MPort>(Vec(204.15, 339), module,
                                       TfReverb::RIGHT_OUTPUT));
  }

  void appendContextMenu(Menu *menu) override {
    auto *reverb = getModule<TfReverb>();
    if (!reverb)
      return;
    menu->addChild(new MenuSeparator);
    menu->addChild(createMenuLabel("Presets"));
    const auto addPreset =
        [=](const char *label,
            const tfdsp::reverb_defaults::ReverbPreset &preset) {
          menu->addChild(createMenuItem(
              label, "", [=]() { reverb->applyPreset(preset, label); }));
        };
    addPreset("Medium Hall (default)", tfdsp::reverb_defaults::MediumHall);
    addPreset("Small Room", tfdsp::reverb_defaults::SmallRoom);
    addPreset("Superlush", tfdsp::reverb_defaults::Superlush);
    menu->addChild(new MenuSeparator);
    menu->addChild(createMenuLabel("Late-tail FDN"));
    for (const auto flavour : {tfdsp::LateReverbFlavour::Base,
                               tfdsp::LateReverbFlavour::Optimized}) {
      const char *label = flavour == tfdsp::LateReverbFlavour::Base
                              ? "Base FDN"
                              : "Optimized FDN";
      menu->addChild(createCheckMenuItem(
          label, "", [=]() { return reverb->lateReverbFlavour() == flavour; },
          [=]() { reverb->setLateReverbFlavour(flavour); }));
    }
  }
};

Model *modelTfReverb = createModel<TfReverb, TfReverbWidget>("TfReverb");
