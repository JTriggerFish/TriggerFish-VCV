#include "plugin.hpp"

#include "components.hpp"

struct TfHiHat : Module {
	enum ParamIds {
		LOCATION_PARAM,
		HARDNESS_PARAM,
		IMPLEMENT_PARAM,
		PEDAL_PARAM,
		DECAY_PARAM,
		CONTACT_PARAM,
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
		PEDAL_INPUT,
		TUNE_INPUT,
		NUM_INPUTS
	};
	enum OutputIds { LEFT_OUTPUT, RIGHT_OUTPUT, NUM_OUTPUTS };
	enum LightIds { NUM_LIGHTS };

	TfHiHat() {
		config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
		configParam(LOCATION_PARAM, 0.f, 1.f, 0.68f, "Hit location", "%",
			0.f, 100.f);
		configParam(HARDNESS_PARAM, 0.f, 1.f, 0.52f, "Stick hardness", "%",
			0.f, 100.f);
		configParam(IMPLEMENT_PARAM, 0.f, 1.f, 0.f,
			"Implement (tip to mallet to brush)", "%", 0.f, 100.f);
		configParam(PEDAL_PARAM, 0.f, 1.f, 0.f, "Pedal closure", "%",
			0.f, 100.f);
		configParam(DECAY_PARAM, 0.f, 1.f, 0.52f, "Open decay", "%",
			0.f, 100.f);
		configParam(CONTACT_PARAM, 0.f, 1.f, 0.62f, "Cymbal contact", "%",
			0.f, 100.f);
		configParam(LEVEL_PARAM, 0.f, 1.f, 0.72f, "Output level", "%",
			0.f, 100.f);
		configParam(TUNE_PARAM, 0.f, 1.f, 0.5f, "Tune", " semitones",
			0.f, 24.f, -12.f);
		configInput(HIT_INPUT, "Hit trigger");
		configInput(STRENGTH_INPUT,
			"Strike strength/velocity (0V to 10V; 8V normalled)");
		configInput(LOCATION_INPUT,
			"Hit location override (0V centre, 10V edge)");
		configInput(HARDNESS_INPUT, "Stick hardness override (0V to 10V)");
		configInput(IMPLEMENT_INPUT,
			"Implement override (0V tip, 5V mallet, 10V brush)");
		configInput(PEDAL_INPUT,
			"Pedal override (0V open, 10V closed; gates produce a foot chick)");
		configInput(TUNE_INPUT,
			"Tune override (0V -12 semitones, 5V unison, 10V +12 semitones)");
		configOutput(LEFT_OUTPUT, "Left audio");
		configOutput(RIGHT_OUTPUT, "Right audio");
	}

	void process(const ProcessArgs &) override {
		// Intentional development shell. Hi-hat synthesis resumes only after the
		// reusable ride components and their quality tests are established.
		outputs[LEFT_OUTPUT].setVoltage(0.f);
		outputs[RIGHT_OUTPUT].setVoltage(0.f);
	}
};

struct TfHiHatWidget : ModuleWidget {
	TfHiHatWidget(TfHiHat *module) {
		setModule(module);
		setPanel(APP->window->loadSvg(
			asset::plugin(pluginInstance, "res/TfHiHat.svg")));
		addChild(createWidget<ScrewSilver>(Vec(0, 0)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 15, 0)));
		addChild(createWidget<ScrewSilver>(Vec(0, RACK_GRID_HEIGHT - 15)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 15,
			RACK_GRID_HEIGHT - 15)));

		addParam(createParam<TfAudioKob>(Vec(2, 45), module,
			TfHiHat::LOCATION_PARAM));
		addParam(createParam<TfAudioKob>(Vec(39, 45), module,
			TfHiHat::HARDNESS_PARAM));
		addParam(createParam<TfAudioKob>(Vec(76, 45), module,
			TfHiHat::IMPLEMENT_PARAM));
		addParam(createParam<TfAudioKob>(Vec(113, 45), module,
			TfHiHat::PEDAL_PARAM));
		addParam(createParam<TfCvKnob>(Vec(2, 113), module,
			TfHiHat::DECAY_PARAM));
		addParam(createParam<TfCvKnob>(Vec(39, 113), module,
			TfHiHat::CONTACT_PARAM));
		addParam(createParam<TfCvKnob>(Vec(76, 113), module,
			TfHiHat::TUNE_PARAM));
		addParam(createParam<TfCvKnob>(Vec(113, 113), module,
			TfHiHat::LEVEL_PARAM));

		addInput(createInput<PJ301MPort>(Vec(1, 194), module,
			TfHiHat::LOCATION_INPUT));
		addInput(createInput<PJ301MPort>(Vec(31, 194), module,
			TfHiHat::HARDNESS_INPUT));
		addInput(createInput<PJ301MPort>(Vec(61, 194), module,
			TfHiHat::IMPLEMENT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(91, 194), module,
			TfHiHat::PEDAL_INPUT));
		addInput(createInput<PJ301MPort>(Vec(121, 194), module,
			TfHiHat::TUNE_INPUT));
		addInput(createInput<PJ301MPort>(Vec(38, 269), module,
			TfHiHat::HIT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(86, 269), module,
			TfHiHat::STRENGTH_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(38, 333), module,
			TfHiHat::LEFT_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(86, 333), module,
			TfHiHat::RIGHT_OUTPUT));
	}
};

Model *modelTfHiHat = createModel<TfHiHat, TfHiHatWidget>("TfHiHat");
