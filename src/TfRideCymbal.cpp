#include "plugin.hpp"

#include "components.hpp"

struct TfRideCymbal : Module {
	enum ParamIds {
		LOCATION_PARAM,
		HARDNESS_PARAM,
		IMPLEMENT_PARAM,
		MUTE_PARAM,
		DECAY_PARAM,
		WASH_PARAM,
		LEVEL_PARAM,
		TUNE_PARAM,
		NUM_PARAMS
	};
	enum InputIds {
		HIT_INPUT,
		STRENGTH_INPUT,
		LOCATION_INPUT,
		HARDNESS_INPUT,
		IMPLEMENT_INPUT,
		MUTE_INPUT,
		TUNE_INPUT,
		NUM_INPUTS
	};
	enum OutputIds { LEFT_OUTPUT, RIGHT_OUTPUT, NUM_OUTPUTS };
	enum LightIds { NUM_LIGHTS };

	TfRideCymbal() {
		config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
		configParam(LOCATION_PARAM, 0.f, 1.f, 0.56f, "Hit location", "%",
			0.f, 100.f);
		configParam(HARDNESS_PARAM, 0.f, 1.f, 0.55f, "Stick hardness", "%",
			0.f, 100.f);
		configParam(IMPLEMENT_PARAM, 0.f, 1.f, 0.f,
			"Implement (tip to mallet to brush)", "%", 0.f, 100.f);
		configParam(MUTE_PARAM, 0.f, 1.f, 0.f, "Hand mute", "%", 0.f, 100.f);
		configParam(DECAY_PARAM, 0.f, 1.f, 0.64f, "Natural decay", "%",
			0.f, 100.f);
		configParam(WASH_PARAM, 0.f, 1.f, 0.68f, "Dense wash", "%",
			0.f, 100.f);
		configParam(LEVEL_PARAM, 0.f, 1.f, 0.72f, "Output level", "%",
			0.f, 100.f);
		configParam(TUNE_PARAM, 0.f, 1.f, 0.5f, "Tune", " semitones",
			0.f, 24.f, -12.f);
		configInput(HIT_INPUT, "Hit trigger");
		configInput(STRENGTH_INPUT,
			"Strike strength/velocity (0V to 10V; 8V normalled)");
		configInput(LOCATION_INPUT,
			"Hit location override (0V bell, 10V edge)");
		configInput(HARDNESS_INPUT, "Stick hardness override (0V to 10V)");
		configInput(IMPLEMENT_INPUT,
			"Implement override (0V tip, 5V mallet, 10V brush)");
		configInput(MUTE_INPUT, "Hand mute override (0V to 10V)");
		configInput(TUNE_INPUT,
			"Tune override (0V -12 semitones, 5V unison, 10V +12 semitones)");
		configOutput(LEFT_OUTPUT, "Left audio");
		configOutput(RIGHT_OUTPUT, "Right audio");
	}

	void process(const ProcessArgs &) override {
		// Intentional development shell. The previous cymbal engine was removed
		// so new audited percussion components are not calibrated around it.
		outputs[LEFT_OUTPUT].setVoltage(0.f);
		outputs[RIGHT_OUTPUT].setVoltage(0.f);
	}
};

struct TfRideCymbalWidget : ModuleWidget {
	TfRideCymbalWidget(TfRideCymbal *module) {
		setModule(module);
		setPanel(APP->window->loadSvg(
			asset::plugin(pluginInstance, "res/TfRideCymbal.svg")));
		addChild(createWidget<ScrewSilver>(Vec(0, 0)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 15, 0)));
		addChild(createWidget<ScrewSilver>(Vec(0, RACK_GRID_HEIGHT - 15)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 15,
			RACK_GRID_HEIGHT - 15)));

		addParam(createParam<TfAudioKob>(Vec(2, 45), module,
			TfRideCymbal::LOCATION_PARAM));
		addParam(createParam<TfAudioKob>(Vec(39, 45), module,
			TfRideCymbal::HARDNESS_PARAM));
		addParam(createParam<TfAudioKob>(Vec(76, 45), module,
			TfRideCymbal::IMPLEMENT_PARAM));
		addParam(createParam<TfAudioKob>(Vec(113, 45), module,
			TfRideCymbal::MUTE_PARAM));
		addParam(createParam<TfCvKnob>(Vec(2, 113), module,
			TfRideCymbal::DECAY_PARAM));
		addParam(createParam<TfCvKnob>(Vec(39, 113), module,
			TfRideCymbal::WASH_PARAM));
		addParam(createParam<TfCvKnob>(Vec(76, 113), module,
			TfRideCymbal::TUNE_PARAM));
		addParam(createParam<TfCvKnob>(Vec(113, 113), module,
			TfRideCymbal::LEVEL_PARAM));

		addInput(createInput<PJ301MPort>(Vec(1, 194), module,
			TfRideCymbal::LOCATION_INPUT));
		addInput(createInput<PJ301MPort>(Vec(31, 194), module,
			TfRideCymbal::HARDNESS_INPUT));
		addInput(createInput<PJ301MPort>(Vec(61, 194), module,
			TfRideCymbal::IMPLEMENT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(91, 194), module,
			TfRideCymbal::MUTE_INPUT));
		addInput(createInput<PJ301MPort>(Vec(121, 194), module,
			TfRideCymbal::TUNE_INPUT));
		addInput(createInput<PJ301MPort>(Vec(38, 269), module,
			TfRideCymbal::HIT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(86, 269), module,
			TfRideCymbal::STRENGTH_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(38, 333), module,
			TfRideCymbal::LEFT_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(86, 333), module,
			TfRideCymbal::RIGHT_OUTPUT));
	}
};

Model *modelTfRideCymbal =
	createModel<TfRideCymbal, TfRideCymbalWidget>("TfRideCymbal");
