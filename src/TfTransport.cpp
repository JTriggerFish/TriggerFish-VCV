#include "plugin.hpp"

#include "components.hpp"
#include "tftransport.hpp"

#include <array>

struct TfTransport : Module {
  enum ParamIds { TEMPO, PLAY_FROM_BEGINNING, PAUSE, PLAY, STOP, NUM_PARAMS };
  enum InputIds { NUM_INPUTS };
  enum OutputIds { CLOCK_OUTPUT, RUN_OUTPUT, RESET_OUTPUT, NUM_OUTPUTS };
  enum LightIds { PLAYING_LIGHT, PAUSED_LIGHT, STOPPED_LIGHT, NUM_LIGHTS };

  tftransport::Engine transport;
  tftransport::CommandMailbox remoteCommands;
  std::array<dsp::SchmittTrigger, 4> commandTriggers{};
  dsp::PulseGenerator clockPulse;
  dsp::PulseGenerator resetPulse;

  TfTransport() {
    config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
    configParam(TEMPO, 30.f, 300.f, 120.f, "Tempo", " BPM");
    getParamQuantity(TEMPO)->displayPrecision = 5;
    configButton(PLAY_FROM_BEGINNING, "Play from beginning");
    configButton(PAUSE, "Pause at the current position");
    configButton(PLAY, "Play from the current position");
    configButton(STOP, "Stop and return to the beginning");
    configOutput(CLOCK_OUTPUT, "24 PPQN master clock");
    configOutput(RUN_OUTPUT, "Run gate");
    configOutput(RESET_OUTPUT, "Reset");
    configLight(PLAYING_LIGHT, "Playing");
    configLight(PAUSED_LIGHT, "Paused");
    configLight(STOPPED_LIGHT, "Stopped");
  }

  json_t *dataToJson() override {
    json_t *root = json_object();
    json_object_set_new(root, "state",
                        json_integer(static_cast<int>(transport.state())));
    return root;
  }

  void dataFromJson(json_t *root) override {
    if (json_t *state = json_object_get(root, "state")) {
      const int value = static_cast<int>(json_integer_value(state));
      transport.loadState(value == static_cast<int>(tftransport::State::Playing)
                              ? tftransport::State::Playing
                          : value ==
                                  static_cast<int>(tftransport::State::Paused)
                              ? tftransport::State::Paused
                              : tftransport::State::Stopped);
    }
  }

  void onReset() override {
    transport.command(tftransport::Command::Stop);
    clockPulse.reset();
    resetPulse.reset();
  }

  void process(const ProcessArgs &args) override {
    constexpr std::array<tftransport::Command, 4> Commands{
        tftransport::Command::PlayFromBeginning, tftransport::Command::Pause,
        tftransport::Command::Play, tftransport::Command::Stop};
    constexpr std::array<int, 4> Parameters{PLAY_FROM_BEGINNING, PAUSE, PLAY,
                                            STOP};
    for (std::size_t index = 0; index < commandTriggers.size(); ++index) {
      if (commandTriggers[index].process(params[Parameters[index]].getValue()))
        transport.command(Commands[index]);
    }
    tftransport::Command remoteCommand;
    if (remoteCommands.consume(remoteCommand))
      transport.command(remoteCommand);

    const auto frame =
        transport.process(args.sampleTime, params[TEMPO].getValue());
    if (frame.clock)
      clockPulse.trigger(0.001f);
    if (frame.reset)
      resetPulse.trigger(0.001f);

    outputs[CLOCK_OUTPUT].setVoltage(clockPulse.process(args.sampleTime) ? 10.f
                                                                         : 0.f);
    outputs[RUN_OUTPUT].setVoltage(frame.run ? 10.f : 0.f);
    outputs[RESET_OUTPUT].setVoltage(resetPulse.process(args.sampleTime) ? 10.f
                                                                         : 0.f);

    const auto state = transport.state();
    lights[PLAYING_LIGHT].setBrightness(
        state == tftransport::State::Playing ? 1.f : 0.f);
    lights[PAUSED_LIGHT].setBrightness(
        state == tftransport::State::Paused ? 1.f : 0.f);
    lights[STOPPED_LIGHT].setBrightness(
        state == tftransport::State::Stopped ? 1.f : 0.f);
  }
};

struct TfTempoDisplay : app::LedDisplayChoice {
  TfTransport *module = nullptr;
  app::ParamWidget *tempoParam = nullptr;

  void step() override {
    text = module ? rack::string::f(
                        "%.2f BPM",
                        module->params[TfTransport::TEMPO].getValue())
                  : "120.00 BPM";
    app::LedDisplayChoice::step();
  }

  void drawLayer(const DrawArgs &args, int layer) override {
    nvgSave(args.vg);
    nvgScissor(args.vg, RECT_ARGS(args.clipBox));
    if (layer == 1) {
      auto font = APP->window->loadFont(fontPath);
      if (font && font->handle >= 0) {
        nvgFillColor(args.vg, color);
        nvgFontFaceId(args.vg, font->handle);
        nvgFontSize(args.vg, 9.f);
        nvgTextLetterSpacing(args.vg, 0.f);
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
        nvgText(args.vg, box.size.x * 0.5f, box.size.y * 0.5f + 0.5f,
                text.c_str(), nullptr);
      }
    }
    Widget::drawLayer(args, layer);
    nvgRestore(args.vg);
  }

  void onButton(const ButtonEvent &event) override {
    if (event.action == GLFW_PRESS && event.button == GLFW_MOUSE_BUTTON_LEFT &&
        tempoParam) {
      tempoParam->createContextMenu();
      event.consume(this);
      return;
    }
    app::LedDisplayChoice::onButton(event);
  }
};

struct TfTransportWidget : ModuleWidget {
  TfTransportWidget(TfTransport *module) {
    setModule(module);
    setPanel(APP->window->loadSvg(
        asset::plugin(pluginInstance, "res/TfTransport.svg")));

    addChild(createWidget<ScrewSilver>(Vec(0, 0)));
    addChild(createWidget<ScrewSilver>(Vec(box.size.x - 15, 0)));
    addChild(createWidget<ScrewSilver>(Vec(0, RACK_GRID_HEIGHT - 15)));
    addChild(
        createWidget<ScrewSilver>(Vec(box.size.x - 15, RACK_GRID_HEIGHT - 15)));

    addParam(
        createParam<TfLargeAudioKnob>(Vec(33, 43), module, TfTransport::TEMPO));
    auto *tempoDisplay = createWidget<TfTempoDisplay>(Vec(23, 101));
    tempoDisplay->box.size = Vec(74, 18);
    tempoDisplay->module = module;
    tempoDisplay->tempoParam = getParam(TfTransport::TEMPO);
    addChild(tempoDisplay);

    addParam(createParam<LEDButton>(Vec(25.38, 155), module,
                                    TfTransport::PLAY_FROM_BEGINNING));
    addParam(
        createParam<LEDButton>(Vec(79.38, 155), module, TfTransport::PAUSE));
    addParam(
        createParam<LEDButton>(Vec(25.38, 215), module, TfTransport::PLAY));
    addParam(
        createParam<LEDButton>(Vec(79.38, 215), module, TfTransport::STOP));

    addChild(createLight<MediumLight<GreenLight>>(Vec(29.331, 272), module,
                                                  TfTransport::PLAYING_LIGHT));
    addChild(createLight<MediumLight<YellowLight>>(Vec(54.331, 272), module,
                                                   TfTransport::PAUSED_LIGHT));
    addChild(createLight<MediumLight<RedLight>>(Vec(79.331, 272), module,
                                                TfTransport::STOPPED_LIGHT));

    addOutput(createOutput<PJ301MPort>(Vec(8.15, 325), module,
                                       TfTransport::CLOCK_OUTPUT));
    addOutput(createOutput<PJ301MPort>(Vec(48.15, 325), module,
                                       TfTransport::RUN_OUTPUT));
    addOutput(createOutput<PJ301MPort>(Vec(88.15, 325), module,
                                       TfTransport::RESET_OUTPUT));
  }
};

bool tftransport::RequestModuleCommand(engine::Module *target,
                                       int sourceOutputId,
                                       Command command) noexcept {
  if (!target || target->model != modelTfTransport ||
      sourceOutputId != TfTransport::RUN_OUTPUT)
    return false;
  auto *transport = dynamic_cast<TfTransport *>(target);
  if (!transport)
    return false;
  transport->remoteCommands.post(command);
  return true;
}

Model *modelTfTransport =
    createModel<TfTransport, TfTransportWidget>("TfTransport");
