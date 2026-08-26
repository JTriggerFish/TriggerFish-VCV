#include <memory>
#include <array>
#include <algorithm>
#include <cmath>
#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/control.hpp"
#include "tfdsp/noise.hpp"

// Analog style modulation of pitch for VCOs and filter cutoffs
struct TfSlop : Module
{
	enum ParamIds
	{
		HUM_LEVEL,
		DRIFT_LEVEL,
		TRACK_SCALING,
		DETUNE_MODE,
		NUM_PARAMS
	};
	enum InputIds
	{
		VOCT_INPUT,
		NUM_INPUTS
	};
	enum OutputIds
	{
		VOCT_OUTPUT,
		NUM_OUTPUTS
	};
	enum LightIds
	{
		NUM_LIGHTS
	};

	std::random_device _seed{};
	std::minstd_rand _rng;
	// 0.01 V/oct is a peak pitch deviation of 12 cents.
	static constexpr float _maxHum{12.0f / 1200.0f};
	static constexpr float _humFreq{60.0f};
	tfdsp::RecursiveSineOscillator _humOscillator{};

	// Temperature drift is modeled as an exact sampled OU process.
	static constexpr double _tau{60.0}; //Time constant ( average decay time) in seconds
	static constexpr double _sigmaCents{0.2 / 12};
	static constexpr double _sigmaHz{2};
	std::array<tfdsp::InterpolatedOrnsteinUhlenbeck,
		PORT_MAX_CHANNELS> _ou{};
	float _prevDetuneMode{};
	float _sampleRate{44100.0f};
	int _activeChannels{};

	//----------------------------------------------------------------

	TfSlop() : _rng(_seed())
	{
		config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
		configParam(TfSlop::HUM_LEVEL, 0.0f, 1.0f, 0.25f, "60 Hz hum depth", " cents peak", 0.0f, 12.0f);
		configParam(TfSlop::DRIFT_LEVEL, 0.0f, 1.0f, 0.05f, "Drift", "%", 0.0f, 100.0f);
		configParam(TfSlop::TRACK_SCALING, 1.0f - 0.2f / 12, 1.0f, 1.0f, "Tracking", "%", 0.0f, 100.0f);
		configSwitch(TfSlop::DETUNE_MODE, -1.0f, 1.0f, -1.0f, "Drift mode", {"Linear Hz", "Proportional cents", "Proportional cents"});
		configInput(VOCT_INPUT, "1V/octave pitch");
		configOutput(VOCT_OUTPUT, "Detuned 1V/octave pitch");
		configBypass(VOCT_INPUT, VOCT_OUTPUT);

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

void TfSlop::init(float sampleRate)
{
	_sampleRate = sampleRate;
	_humOscillator.SetFrequency(_humFreq, sampleRate);
	const double sigma = params[DETUNE_MODE].getValue() < 0 ? _sigmaHz : _sigmaCents;
	for (auto& process : _ou)
		process.Configure(sampleRate, _tau, sigma);
	_prevDetuneMode = params[DETUNE_MODE].getValue();
}

void TfSlop::process(const ProcessArgs &args)
{
	if (_prevDetuneMode != params[DETUNE_MODE].getValue())
	{
		const double sigma = params[DETUNE_MODE].getValue() < 0 ? _sigmaHz : _sigmaCents;
		for (auto& process : _ou)
		{
			process.Reset();
			process.Configure(args.sampleRate, _tau, sigma);
		}
		_prevDetuneMode = params[DETUNE_MODE].getValue();
	}
	const int channels = std::clamp(std::max(
		inputs[VOCT_INPUT].getChannels(), 1), 1, PORT_MAX_CHANNELS);
	for (int channel = channels; channel < _activeChannels; ++channel)
		_ou[channel].Reset();
	_activeChannels = channels;
	outputs[VOCT_OUTPUT].setChannels(channels);

	const float tracking = params[TRACK_SCALING].getValue();
	const float hum = _maxHum * params[HUM_LEVEL].getValue() *
		_humOscillator.Step();
	const float driftDepth = params[DRIFT_LEVEL].getValue();
	const bool linearDrift = params[DETUNE_MODE].getValue() < 0;
	for (int channel = 0; channel < channels; ++channel)
	{
		const float input = inputs[VOCT_INPUT].getPolyVoltage(channel);
		const double voct = (std::isfinite(input) ? input : 0.0f) * tracking + hum;
		const double drift = driftDepth * _ou[channel].Step(_rng);
		const double output = linearDrift ?
			tfdsp::detune::linear(voct, drift) : voct + drift;
		outputs[VOCT_OUTPUT].setVoltage(std::isfinite(output) ?
			static_cast<float>(output) : 0.0f, channel);
	}
}
void TfSlop::onSampleRateChange(const SampleRateChangeEvent& event)
{
	init(event.sampleRate);
}
void TfSlop::onReset(const ResetEvent& event)
{
	Module::onReset(event);
	_humOscillator.Reset();
	for (auto& process : _ou)
		process.Reset();
	_activeChannels = 0;
	init(_sampleRate);
}

struct TfSlopWidget : ModuleWidget
{
	TfSlopWidget(TfSlop *module)
	{
		setModule(module);
		setPanel(APP->window->loadSvg(asset::plugin(pluginInstance, "res/TfSlop.svg")));

		//Panel screws
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		//Knobs
		addParam(createParam<TfCvKnob>(Vec(30, 55), module, TfSlop::HUM_LEVEL));
		addParam(createParam<TfCvKnob>(Vec(10, 127), module, TfSlop::DRIFT_LEVEL));
		addParam(createParam<TfCvKnob>(Vec(30, 190), module, TfSlop::TRACK_SCALING));

		//Drift mode switch
		addParam(createParam<CKSS>(Vec(65, 135), module, TfSlop::DETUNE_MODE));

		//Jacks at the bottom
		addInput(createInput<PJ301MPort>(Vec(13.5, 317), module, TfSlop::VOCT_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(55, 317), module, TfSlop::VOCT_OUTPUT));
	}
};

// Specify the Module and ModuleWidget subclass, human-readable
// author name for categorization per pluginInstance, module slug (should never
// change), human-readable module name, and any number of tags
// (found in `include/tags.hpp`) separated by commas.
Model *modelTfSlop = createModel<TfSlop, TfSlopWidget>("TfSlop");
