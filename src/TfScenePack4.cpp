#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/scene_pack.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>

struct TfScenePack4 : Module {
  static constexpr int LaneCount = 4;
  static constexpr int MaximumSources = 8;

  enum ParamIds {
    X_1,
    Y_1,
    Z_1,
    X_2,
    Y_2,
    Z_2,
    X_3,
    Y_3,
    Z_3,
    X_4,
    Y_4,
    Z_4,
    NUM_PARAMS
  };
  enum InputIds {
    BUS_AUDIO_INPUT,
    BUS_X_INPUT,
    BUS_Y_INPUT,
    BUS_Z_INPUT,
    SOURCE_1_INPUT,
    SOURCE_2_INPUT,
    SOURCE_3_INPUT,
    SOURCE_4_INPUT,
    NUM_INPUTS
  };
  enum OutputIds { AUDIO_OUTPUT, X_OUTPUT, Y_OUTPUT, Z_OUTPUT, NUM_OUTPUTS };
  enum LightIds { NUM_LIGHTS };

  TfScenePack4() {
    config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
    constexpr std::array<float, LaneCount> defaultX{3.5f, 4.5f, 5.5f, 6.5f};
    for (int lane = 0; lane < LaneCount; ++lane) {
      const int first = X_1 + 3 * lane;
      const std::string laneName = "Source " + std::to_string(lane + 1) + " ";
      configParam(first, 0.f, 10.f, defaultX[lane], laneName + "X position");
      configParam(first + 1, 0.f, 10.f, 3.5f, laneName + "Y position");
      configParam(first + 2, 0.f, 10.f, 4.5f, laneName + "Z position");
      configInput(SOURCE_1_INPUT + lane, laneName + "mono audio");
    }
    configInput(BUS_AUDIO_INPUT, "Existing scene audio");
    configInput(BUS_X_INPUT, "Existing scene X positions");
    configInput(BUS_Y_INPUT, "Existing scene Y positions");
    configInput(BUS_Z_INPUT, "Existing scene Z positions");
    configOutput(AUDIO_OUTPUT, "Scene audio");
    configOutput(X_OUTPUT, "Scene X positions");
    configOutput(Y_OUTPUT, "Scene Y positions");
    configOutput(Z_OUTPUT, "Scene Z positions");
  }

  void process(const ProcessArgs &) override {
    tfdsp::ScenePackInput packInput;
    if (inputs[BUS_AUDIO_INPUT].isConnected()) {
      const int busChannels =
          std::min(MaximumSources, inputs[BUS_AUDIO_INPUT].getChannels());
      packInput.busCount = static_cast<std::size_t>(busChannels);
      for (int channel = 0; channel < busChannels; ++channel) {
        auto &source = packInput.bus[static_cast<std::size_t>(channel)];
        source.audio = inputs[BUS_AUDIO_INPUT].getVoltage(channel);
        source.x = positionVoltage(BUS_X_INPUT, channel);
        source.y = positionVoltage(BUS_Y_INPUT, channel);
        source.z = positionVoltage(BUS_Z_INPUT, channel);
      }
    }

    for (int lane = 0; lane < LaneCount; ++lane) {
      auto &input = inputs[SOURCE_1_INPUT + lane];
      packInput.localConnected[static_cast<std::size_t>(lane)] =
          input.isConnected();
      auto &source = packInput.local[static_cast<std::size_t>(lane)];
      source.audio = input.getVoltageSum();
      const int first = X_1 + 3 * lane;
      source.x = params[first].getValue();
      source.y = params[first + 1].getValue();
      source.z = params[first + 2].getValue();
    }

    const auto packed = tfdsp::PackSceneSources(packInput);
    for (std::size_t channel = 0; channel < packed.sourceCount; ++channel) {
      const auto &source = packed.sources[channel];
      outputs[AUDIO_OUTPUT].setVoltage(source.audio, channel);
      outputs[X_OUTPUT].setVoltage(source.x, channel);
      outputs[Y_OUTPUT].setVoltage(source.y, channel);
      outputs[Z_OUTPUT].setVoltage(source.z, channel);
    }
    for (int output = AUDIO_OUTPUT; output < NUM_OUTPUTS; ++output)
      outputs[output].setChannels(static_cast<int>(packed.sourceCount));
  }

private:
  float positionVoltage(const int inputId, const int channel) noexcept {
    auto &input = inputs[inputId];
    if (!input.isConnected())
      return 5.f;
    return tfdsp::SanitizeScenePosition(input.getPolyVoltage(channel));
  }
};

struct TfScenePack4Widget : ModuleWidget {
  TfScenePack4Widget(TfScenePack4 *module) {
    setModule(module);
    setPanel(APP->window->loadSvg(
        asset::plugin(pluginInstance, "res/TfScenePack4.svg")));

    addChild(createWidget<ScrewSilver>(Vec(15, 0)));
    addChild(createWidget<ScrewSilver>(Vec(box.size.x - 30, 0)));
    addChild(createWidget<ScrewSilver>(Vec(15, RACK_GRID_HEIGHT - 15)));
    addChild(
        createWidget<ScrewSilver>(Vec(box.size.x - 30, RACK_GRID_HEIGHT - 15)));

    addInput(createInput<PJ301MPort>(Vec(18, 49), module,
                                     TfScenePack4::SOURCE_1_INPUT));
    addParam(createParam<TfTrimpot>(Vec(20, 94), module, TfScenePack4::X_1));
    addParam(createParam<TfTrimpot>(Vec(20, 128), module, TfScenePack4::Y_1));
    addParam(createParam<TfTrimpot>(Vec(20, 162), module, TfScenePack4::Z_1));

    addInput(createInput<PJ301MPort>(Vec(65, 49), module,
                                     TfScenePack4::SOURCE_2_INPUT));
    addParam(createParam<TfTrimpot>(Vec(67, 94), module, TfScenePack4::X_2));
    addParam(createParam<TfTrimpot>(Vec(67, 128), module, TfScenePack4::Y_2));
    addParam(createParam<TfTrimpot>(Vec(67, 162), module, TfScenePack4::Z_2));

    addInput(createInput<PJ301MPort>(Vec(112, 49), module,
                                     TfScenePack4::SOURCE_3_INPUT));
    addParam(createParam<TfTrimpot>(Vec(114, 94), module, TfScenePack4::X_3));
    addParam(createParam<TfTrimpot>(Vec(114, 128), module, TfScenePack4::Y_3));
    addParam(createParam<TfTrimpot>(Vec(114, 162), module, TfScenePack4::Z_3));

    addInput(createInput<PJ301MPort>(Vec(159, 49), module,
                                     TfScenePack4::SOURCE_4_INPUT));
    addParam(createParam<TfTrimpot>(Vec(161, 94), module, TfScenePack4::X_4));
    addParam(createParam<TfTrimpot>(Vec(161, 128), module, TfScenePack4::Y_4));
    addParam(createParam<TfTrimpot>(Vec(161, 162), module, TfScenePack4::Z_4));

    addInput(createInput<PJ301MPort>(Vec(18, 223), module,
                                     TfScenePack4::BUS_AUDIO_INPUT));
    addInput(createInput<PJ301MPort>(Vec(65, 223), module,
                                     TfScenePack4::BUS_X_INPUT));
    addInput(createInput<PJ301MPort>(Vec(112, 223), module,
                                     TfScenePack4::BUS_Y_INPUT));
    addInput(createInput<PJ301MPort>(Vec(159, 223), module,
                                     TfScenePack4::BUS_Z_INPUT));

    addOutput(createOutput<PJ301MPort>(Vec(18, 306), module,
                                       TfScenePack4::AUDIO_OUTPUT));
    addOutput(createOutput<PJ301MPort>(Vec(65, 306), module,
                                       TfScenePack4::X_OUTPUT));
    addOutput(createOutput<PJ301MPort>(Vec(112, 306), module,
                                       TfScenePack4::Y_OUTPUT));
    addOutput(createOutput<PJ301MPort>(Vec(159, 306), module,
                                       TfScenePack4::Z_OUTPUT));
  }
};

Model *modelTfScenePack4 =
    createModel<TfScenePack4, TfScenePack4Widget>("TfScenePack4");
