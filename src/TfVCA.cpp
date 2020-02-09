#include "models/VCAcore.hpp"

#include <memory>
#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/filters.hpp"
#include "tfdsp/sampleRate.hpp"


// Analog modelled VCA with 2x oversampling
struct TfVCA : Module 
{
	enum ParamIds {
		INPUT_GAIN,
		LIN_INPUT_LEVEL,
		EXP_INPUT_LEVEL,
		CV_BLEED,
		EXP_CV_BASE,
		OUTPUT_LEVEL,
		NUM_PARAMS
	};
	enum InputIds {
		AUDIO_INPUT,
		LIN_CV_INPUT,
		EXP_CV_INPUT,
		NUM_INPUTS
	};
	enum OutputIds {
		MAIN_OUTPUT,
		NUM_OUTPUTS
	};
	enum LightIds {
		CV_LIGHT,
		NUM_LIGHTS
	};

	static const float _maxCvBleed;
	static constexpr float _cvBleedHighPassF = 10.f;
	static constexpr float _audioHighPassF = 5.0f;

	float _normalisedHighPassCv;
	float _normalisedHighPassAudio;

	std::unique_ptr<::VCA_TransistorCore<tfdsp::X2Resampler_Order7> >_vcaTransi;

	tfdsp::FirstOrderHighPassZdf<float> _cvHighPass{};
	tfdsp::FirstOrderHighPassZdf<float> _audioHighPass{};

	//----------------------------------------------------------------

	TfVCA() : _vcaTransi(new ::VCA_TransistorCore<tfdsp::X2Resampler_Order7>(tfdsp::CreateX2Resampler_Chebychev7))
	{
    config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
    configParam(TfVCA::LIN_INPUT_LEVEL, 0.0f, 1.0f, 1.0f, "");
    configParam(TfVCA::EXP_INPUT_LEVEL, 0.0f, 1.0f, 0.0f, "");
    configParam(TfVCA::INPUT_GAIN, 0.0f, 2.0f, 0.5f, "");
    configParam(TfVCA::OUTPUT_LEVEL, 0.0f, 2.0f, 1.0f, "");
    configParam(TfVCA::EXP_CV_BASE, 2.0f, 50.0f, 50.0f, "");
    configParam(TfVCA::CV_BLEED, 0.0f, 1.0f, 0.5f, "");

		//_resampler = tfdsp::CreateX2Resampler_Butterworth5();
    float gSampleRate = APP->engine->getSampleRate();
		init(gSampleRate);
	}

	void process(const ProcessArgs& args) override;
	void init(float sampleRate);
	void onSampleRateChange() override;


	// For more advanced Module features, read Rack's engine.hpp header file
	// - dataToJson, dataFromJson: serialization of internal data
	// - onReset, onRandomize, onCreate, onDelete: implements special behavior when user clicks these from the context menu
};

const float TfVCA::_maxCvBleed = 1.41f * std::pow(10.f, -20.f / 20);

void TfVCA::init(float sampleRate)
{
	_vcaTransi->SetSampleRate(sampleRate);

	 _normalisedHighPassCv = _cvBleedHighPassF / (0.5f * sampleRate);
	 _normalisedHighPassAudio = _audioHighPassF / (0.5f * sampleRate);
}

void TfVCA::process(const ProcessArgs& args) {
	//float deltaTime = args.sampleTime;

	float inputGain = params[INPUT_GAIN].getValue();
	constexpr float audioRenorm = 5.0f;
	inputGain /= audioRenorm;

	float audio = inputs[AUDIO_INPUT].getVoltage() * inputGain;
	//VCA cv should be unipolar between 0 and 10, normalise to 0 to 1. 
	//If no input plugged in then pass zero.
	float linCv = inputs[LIN_CV_INPUT].getNormalVoltage(0.f) / 10.f * params[LIN_INPUT_LEVEL].getValue();
	float expCv = inputs[EXP_CV_INPUT].getNormalVoltage(0.f) / 10.f * params[EXP_INPUT_LEVEL].getValue();

	float expBase = params[EXP_CV_BASE].getValue();
	expCv = (powf(expBase, expCv) - 1.f) / (expBase - 1.f);

	float cv = linCv + expCv;

	//cv bleed
	outputs[MAIN_OUTPUT].setVoltage(_cvHighPass(cv, _normalisedHighPassCv) * params[CV_BLEED].getValue() * _maxCvBleed);


	//Renormalise the audio so that the output level and the input gain
	//are more orthogonal, where the input gain is mostly used for distortion
	auto finalGain = std::min(100.0f, (1.0f + inputGain) / (0.00001f + inputGain));
	finalGain *= params[OUTPUT_LEVEL].getValue();

	//VCA core
	audio = _vcaTransi->Step(audio, cv, finalGain);

	//DC rejection in case there is some DC offset due to aliasing
	audio = _audioHighPass(audio, _normalisedHighPassAudio);

	outputs[MAIN_OUTPUT].value += audio;

	//Deal with input monitoring lights
	lights[CV_LIGHT].setSmoothBrightness(std::max(0.f, cv), args.sampleTime);

}
void TfVCA::onSampleRateChange()
{
	float gSampleRate = APP->engine->getSampleRate();
	init(gSampleRate);
}


struct TfVCAWidget : ModuleWidget {
	TfVCAWidget(TfVCA *module) {
		setModule(module);
		setPanel(APP->window->loadSvg(asset::plugin(pluginInstance, "res/TfVCA.svg")));

		//Panel screws
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		//KnobsAudio
		addParam(createParam<TfCvKnob>(Vec(26,45.5), module, TfVCA::LIN_INPUT_LEVEL));
		addParam(createParam<TfCvKnob>(Vec(26, 104), module, TfVCA::EXP_INPUT_LEVEL));
		addParam(createParam<TfLargeAudioKnob>(Vec(108,79), module, TfVCA::INPUT_GAIN));
		addParam(createParam<TfAudioKob>(Vec(72,154), module, TfVCA::OUTPUT_LEVEL));

		addParam(createParam<TfTrimpot>(Vec(38,245), module, TfVCA::EXP_CV_BASE));
		addParam(createParam<TfTrimpot>(Vec(121,245), module, TfVCA::CV_BLEED));

		//Activity led
		addChild(createLight<MediumLight<BlueLight>>(Vec(85, 250), module, TfVCA::CV_LIGHT));

		//Jacks at the bottom
		constexpr float offset = 15.0f;
		constexpr float spacing = 42.0f;
		addInput(createInput<PJ301MPort>(Vec(offset, 313), module, TfVCA::LIN_CV_INPUT));
		addInput(createInput<PJ301MPort>(Vec(offset + spacing, 313), module, TfVCA::EXP_CV_INPUT));
		addInput(createInput<PJ301MPort>(Vec(offset + 2*spacing, 313 ), module, TfVCA::AUDIO_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(offset + 3*spacing, 313), module, TfVCA::MAIN_OUTPUT));

	}
};


// Specify the Module and ModuleWidget subclass, human-readable
// author name for categorization per pluginInstance, module slug (should never
// change), human-readable module name, and any number of tags
// (found in `include/tags.hpp`) separated by commas.
Model *modelTfVCA = createModel<TfVCA, TfVCAWidget>("TfVCA");
