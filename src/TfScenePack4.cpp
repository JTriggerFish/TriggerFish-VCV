#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/scene_pack.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

struct TfScenePack4 : Module {
  static constexpr int InputCount =
      static_cast<int>(tfdsp::ScenePackInputCount);
  static constexpr int MaximumSources =
      static_cast<int>(tfdsp::ScenePackMaximumSourceCount);

  enum ParamIds { NUM_PARAMS };
  enum InputIds {
    SOURCE_1_INPUT,
    SOURCE_2_INPUT,
    SOURCE_3_INPUT,
    SOURCE_4_INPUT,
    NUM_INPUTS
  };
  enum OutputIds { AUDIO_OUTPUT, NUM_OUTPUTS };
  enum LightIds { NUM_LIGHTS };

  TfScenePack4() {
    config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
    for (int port = 0; port < InputCount; ++port)
      configInput(SOURCE_1_INPUT + port,
                  "Source bundle " + std::to_string(port + 1));
    configOutput(AUDIO_OUTPUT, "Packed source audio");
  }

  void process(const ProcessArgs &) override {
    tfdsp::ScenePackInput packInput;
    for (int port = 0; port < InputCount; ++port) {
      auto &input = inputs[SOURCE_1_INPUT + port];
      if (!input.isConnected())
        continue;
      // A connected upstream Scene Pack can legitimately expose zero channels
      // while all of its inputs are unused. Preserve that empty bundle rather
      // than manufacturing a silent source and a phantom reverb marker.
      const int channels = std::clamp(input.getChannels(), 0, MaximumSources);
      packInput.channelCounts[static_cast<std::size_t>(port)] =
          static_cast<std::size_t>(channels);
      for (int channel = 0; channel < channels; ++channel)
        packInput.inputs[static_cast<std::size_t>(port)]
                        [static_cast<std::size_t>(channel)] =
            input.getVoltage(channel);
    }

    const auto packed = tfdsp::PackSceneSources(packInput);
    for (std::size_t source = 0; source < packed.sourceCount; ++source)
      outputs[AUDIO_OUTPUT].setVoltage(packed.audio[source], source);
    outputs[AUDIO_OUTPUT].setChannels(static_cast<int>(packed.sourceCount));
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

    addInput(createInput<PJ301MPort>(Vec(13.15, 87), module,
                                     TfScenePack4::SOURCE_1_INPUT));
    addInput(createInput<PJ301MPort>(Vec(53.15, 87), module,
                                     TfScenePack4::SOURCE_2_INPUT));
    addInput(createInput<PJ301MPort>(Vec(13.15, 149), module,
                                     TfScenePack4::SOURCE_3_INPUT));
    addInput(createInput<PJ301MPort>(Vec(53.15, 149), module,
                                     TfScenePack4::SOURCE_4_INPUT));

    addOutput(createOutput<PJ301MPort>(Vec(33.15, 287), module,
                                       TfScenePack4::AUDIO_OUTPUT));
  }
};

Model *modelTfScenePack4 =
    createModel<TfScenePack4, TfScenePack4Widget>("TfScenePack4");
