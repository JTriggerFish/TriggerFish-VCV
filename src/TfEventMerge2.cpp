#include "plugin.hpp"

#include "components.hpp"
#include "tfdsp/event_merge.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>

struct TfEventMerge2 : Module {
  enum ParamIds { NUM_PARAMS };
  enum InputIds {
    A_PITCH_INPUT,
    A_GATE_INPUT,
    A_TRIGGER_INPUT,
    A_VELOCITY_INPUT,
    A_ACCENT_INPUT,
    B_PITCH_INPUT,
    B_GATE_INPUT,
    B_TRIGGER_INPUT,
    B_VELOCITY_INPUT,
    B_ACCENT_INPUT,
    NUM_INPUTS
  };
  enum OutputIds {
    PITCH_OUTPUT,
    GATE_OUTPUT,
    TRIGGER_OUTPUT,
    VELOCITY_OUTPUT,
    ACCENT_OUTPUT,
    NUM_OUTPUTS
  };
  enum LightIds { NUM_LIGHTS };

  static constexpr std::array<const char *, tfdsp::EventSignalCount>
      SignalNames{"Pitch", "Gate", "Trigger", "Velocity", "Accent"};

  TfEventMerge2() {
    config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
    for (std::size_t signal = 0; signal < SignalNames.size(); ++signal) {
      configInput(A_PITCH_INPUT + static_cast<int>(signal),
                  std::string("A ") + SignalNames[signal]);
      configInput(B_PITCH_INPUT + static_cast<int>(signal),
                  std::string("B ") + SignalNames[signal]);
      configOutput(PITCH_OUTPUT + static_cast<int>(signal), SignalNames[signal]);
    }
  }

  void process(const ProcessArgs &) override {
    tfdsp::EventBundle bundles[2];
    for (int source = 0; source < 2; ++source) {
      const int firstInput = source == 0 ? A_PITCH_INPUT : B_PITCH_INPUT;
      int channels = 0;
      for (std::size_t signal = 0; signal < tfdsp::EventSignalCount; ++signal)
        channels = std::max(channels,
                            inputs[firstInput + static_cast<int>(signal)]
                                .getChannels());
      bundles[source].voiceCount = static_cast<std::size_t>(
          std::clamp(channels, 0, static_cast<int>(tfdsp::EventVoiceLimit)));
      for (std::size_t signal = 0; signal < tfdsp::EventSignalCount; ++signal)
        for (std::size_t voice = 0; voice < bundles[source].voiceCount;
             ++voice)
          bundles[source].signals[signal][voice] =
              inputs[firstInput + static_cast<int>(signal)].getPolyVoltage(
                  static_cast<int>(voice));
    }

    const auto merged = tfdsp::MergeEventBundles(bundles[0], bundles[1]);
    for (std::size_t signal = 0; signal < tfdsp::EventSignalCount; ++signal) {
      auto &output = outputs[PITCH_OUTPUT + static_cast<int>(signal)];
      for (std::size_t voice = 0; voice < merged.voiceCount; ++voice)
        output.setVoltage(merged.signals[signal][voice],
                          static_cast<int>(voice));
      output.setChannels(static_cast<int>(merged.voiceCount));
    }
  }
};

struct TfEventMerge2Widget : ModuleWidget {
  TfEventMerge2Widget(TfEventMerge2 *module) {
    setModule(module);
    setPanel(APP->window->loadSvg(
        asset::plugin(pluginInstance, "res/TfEventMerge2.svg")));

    addChild(createWidget<ScrewSilver>(Vec(15, 0)));
    addChild(createWidget<ScrewSilver>(Vec(box.size.x - 30, 0)));
    addChild(createWidget<ScrewSilver>(Vec(15, RACK_GRID_HEIGHT - 15)));
    addChild(
        createWidget<ScrewSilver>(Vec(box.size.x - 30, RACK_GRID_HEIGHT - 15)));

    constexpr std::array<float, tfdsp::EventSignalCount> rows{
        70.f, 124.f, 178.f, 232.f, 286.f};
    for (std::size_t signal = 0; signal < rows.size(); ++signal) {
      addInput(createInput<PJ301MPort>(Vec(18.f, rows[signal]), module,
                                      TfEventMerge2::A_PITCH_INPUT +
                                          static_cast<int>(signal)));
      addInput(createInput<PJ301MPort>(Vec(78.f, rows[signal]), module,
                                      TfEventMerge2::B_PITCH_INPUT +
                                          static_cast<int>(signal)));
      addOutput(createOutput<PJ301MPort>(Vec(138.f, rows[signal]), module,
                                        TfEventMerge2::PITCH_OUTPUT +
                                            static_cast<int>(signal)));
    }
  }
};

Model *modelTfEventMerge2 =
    createModel<TfEventMerge2, TfEventMerge2Widget>("TfEventMerge2");
