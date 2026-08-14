#include <memory>
#include <array>
#include <cmath>
#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/control.hpp"
#include "tfdsp/noise.hpp"

// Like TfSlop but with 4 outputs and a common drift on top of the idiosyncratic drifts
struct TfSlop4 : Module
{
	enum ParamIds
	{
		TRACK_SCALING1,
		TRACK_SCALING2,
		TRACK_SCALING3,
		TRACK_SCALING4,
		HUM_LEVEL,
		COMMON_DRIFT_LEVEL,
		INDIVIDUAL_DRIFT_LEVEL,
		NUM_PARAMS
	};
	enum InputIds
	{
		VOCT_INPUT1,
		VOCT_INPUT2,
		VOCT_INPUT3,
		VOCT_INPUT4,
		NUM_INPUTS
	};
	enum OutputIds
	{
		VOCT_OUTPUT1,
		VOCT_OUTPUT2,
		VOCT_OUTPUT3,
		VOCT_OUTPUT4,
		NUM_OUTPUTS
	};
	enum LightIds
	{
		NUM_LIGHTS
	};

	std::random_device _seed{};
	std::minstd_rand _rng;
	// 0.01 V/oct is a peak pitch deviation of 12 cents, common to all outputs.
	static constexpr float _maxHum{12.0f / 1200.0f};
	static constexpr float _humFreq{60.0f};
	tfdsp::RecursiveSineOscillator _humOscillator{};

	// Temperature drift is modeled as an exact sampled OU process.
	static constexpr double _tau{60.0}; //Time constant ( average decay time) in seconds
	static constexpr double _sigmaCents{0.1 / 12};
	static constexpr double _sigmaHz{1.5};
	tfdsp::InterpolatedOrnsteinUhlenbeck _ouCommon{};
	std::array<tfdsp::InterpolatedOrnsteinUhlenbeck, 4> _ouIndividual{};
	float _sampleRate{44100.0f};

	//----------------------------------------------------------------

	TfSlop4() : _rng(_seed())
	{
		config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
		configParam(TfSlop4::HUM_LEVEL, 0.0f, 1.0f, 0.10f, "Common 60 Hz hum depth", " cents peak", 0.0f, 12.0f);
		configParam(TfSlop4::COMMON_DRIFT_LEVEL, 0.0f, 1.0f, 0.05f, "Common drift", "%", 0.0f, 100.0f);
		configParam(TfSlop4::INDIVIDUAL_DRIFT_LEVEL, 0.0f, 1.0f, 0.05f, "Individual drift", "%", 0.0f, 100.0f);
		configParam(TfSlop4::TRACK_SCALING1, 1.0f - 0.2f / 12, 1.0f, 1.0f, "Channel 1 tracking", "%", 0.0f, 100.0f);
		configParam(TfSlop4::TRACK_SCALING2, 1.0f - 0.2f / 12, 1.0f, 1.0f, "Channel 2 tracking", "%", 0.0f, 100.0f);
		configParam(TfSlop4::TRACK_SCALING3, 1.0f - 0.2f / 12, 1.0f, 1.0f, "Channel 3 tracking", "%", 0.0f, 100.0f);
		configParam(TfSlop4::TRACK_SCALING4, 1.0f - 0.2f / 12, 1.0f, 1.0f, "Channel 4 tracking", "%", 0.0f, 100.0f);
		for (int i = 0; i < 4; ++i)
		{
			configInput(VOCT_INPUT1 + i, "Channel " + std::to_string(i + 1) + " 1V/octave pitch");
			configOutput(VOCT_OUTPUT1 + i, "Channel " + std::to_string(i + 1) + " detuned 1V/octave pitch");
			configBypass(VOCT_INPUT1 + i, VOCT_OUTPUT1 + i);
		}
		//_resampler = tfdsp::CreateX2Resampler_Butterworth5();
		float gSampleRate = APP->engine->getSampleRate();
		init(gSampleRate);
	}

	void process(const ProcessArgs &args) override;
	void init(float sampleRate);
	void onSampleRateChange(const SampleRateChangeEvent& event) override;
	void onReset(const ResetEvent& event) override;

	// For more advanced Module features, read Rack's engine.hpp header file
	// - dataToJson, dataFromJson: serialization of internal data
	// - onReset, onRandomize, onCreate, onDelete: implements special behavior when user clicks these from the context menu
};

void TfSlop4::init(float sampleRate)
{
	_sampleRate = sampleRate;
	_humOscillator.SetFrequency(_humFreq, sampleRate);
	_ouCommon.Configure(sampleRate, _tau, _sigmaCents);
	for (auto& process : _ouIndividual)
		process.Configure(sampleRate, _tau, _sigmaHz);
}

void TfSlop4::process(const ProcessArgs &args)
{
	std::array<float, 4> voct;
	for (int i = 0; i < 4; ++i)
	{
		//NOTE! : the parameters that are replicated for each input are put at the beginning to make life easier in loops
		//careful not to add parameters before these in the enum !
		const float input = inputs[i].getVoltage();
		voct[i] = (std::isfinite(input) ? input : 0.0f) * params[i].getValue();
	}

	float hum = _maxHum * params[HUM_LEVEL].getValue() * _humOscillator.Step();

	//The common drift operates in cents
	float driftCommon = params[COMMON_DRIFT_LEVEL].getValue() * _ouCommon.Step(_rng);

	for (int i = 0; i < 4; ++i)
	{
		//The individual drifts operate in hz for linear detuning
		double v = voct[i] + hum + driftCommon;
		double drift = params[INDIVIDUAL_DRIFT_LEVEL].getValue() * _ouIndividual[i].Step(_rng);
		outputs[i].setVoltage(tfdsp::detune::linear(v, drift));
	}
}
void TfSlop4::onSampleRateChange(const SampleRateChangeEvent& event)
{
	init(event.sampleRate);
}
void TfSlop4::onReset(const ResetEvent& event)
{
	Module::onReset(event);
	_humOscillator.Reset();
	_ouCommon.Reset();
	for (auto& process : _ouIndividual)
		process.Reset();
	init(_sampleRate);
}

struct TfSlop4Widget : ModuleWidget
{
	TfSlop4Widget(TfSlop4 *module)
	{
		setModule(module);
		setPanel(APP->window->loadSvg(asset::plugin(pluginInstance, "res/TfSlop4.svg")));

		//Panel screws
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		//Knobs
		addParam(createParam<TfCvKnob>(Vec(61, 66), module, TfSlop4::HUM_LEVEL));
		addParam(createParam<TfCvKnob>(Vec(16, 133), module, TfSlop4::COMMON_DRIFT_LEVEL));
		addParam(createParam<TfCvKnob>(Vec(105, 133), module, TfSlop4::INDIVIDUAL_DRIFT_LEVEL));

		//Tracking trimmers
		auto leftMargin = 13;
		auto spacing = 35;
		addParam(createParam<TfTrimpot>(Vec(leftMargin, 223), module, TfSlop4::TRACK_SCALING1));
		addParam(createParam<TfTrimpot>(Vec(leftMargin + spacing, 223), module, TfSlop4::TRACK_SCALING2));
		addParam(createParam<TfTrimpot>(Vec(leftMargin + 2 * spacing, 223), module, TfSlop4::TRACK_SCALING3));
		addParam(createParam<TfTrimpot>(Vec(leftMargin + 3 * spacing, 223), module, TfSlop4::TRACK_SCALING4));

		//Input jacks
		leftMargin = 10;
		addInput(createInput<PJ301MPort>(Vec(leftMargin, 283), module, TfSlop4::VOCT_INPUT1));
		addInput(createInput<PJ301MPort>(Vec(leftMargin + spacing, 283), module, TfSlop4::VOCT_INPUT2));
		addInput(createInput<PJ301MPort>(Vec(leftMargin + 2 * spacing, 283), module, TfSlop4::VOCT_INPUT3));
		addInput(createInput<PJ301MPort>(Vec(leftMargin + 3 * spacing, 283), module, TfSlop4::VOCT_INPUT4));

		//Output jacks
		addOutput(createOutput<PJ301MPort>(Vec(leftMargin, 319), module, TfSlop4::VOCT_OUTPUT1));
		addOutput(createOutput<PJ301MPort>(Vec(leftMargin + spacing, 319), module, TfSlop4::VOCT_OUTPUT2));
		addOutput(createOutput<PJ301MPort>(Vec(leftMargin + 2 * spacing, 319), module, TfSlop4::VOCT_OUTPUT3));
		addOutput(createOutput<PJ301MPort>(Vec(leftMargin + 3 * spacing, 319), module, TfSlop4::VOCT_OUTPUT4));
	}
};

// Specify the Module and ModuleWidget subclass, human-readable
// author name for categorization per pluginInstance, module slug (should never
// change), human-readable module name, and any number of tags
// (found in `include/tags.hpp`) separated by commas.
Model *modelTfSlop4 = createModel<TfSlop4, TfSlop4Widget>("TfSlop4");
