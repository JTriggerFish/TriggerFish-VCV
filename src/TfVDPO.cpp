#include "models/VdpSplitOscillator.hpp"
#include <memory>
#include <array>
#include <algorithm>
#include <cmath>
#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/noise.hpp"

// Analog style modulation of pitch for VCOs and filter cutoffs
struct TfVDPO : Module
{
	enum ParamIds
	{
		DAMPING,
		FREQ,
		INPUT_GAIN,
		LEVEL,
		VOCT_SCALING,
		DAMPING_ATTENUVERT,
		NUM_PARAMS
	};
	enum InputIds
	{
		VOCT_INPUT,
		AUDIO_INPUT,
		DAMPING_INPUT,
		NUM_INPUTS
	};
	enum OutputIds
	{
		OUTPUT,
		NUM_OUTPUTS
	};
	enum LightIds
	{
		NUM_LIGHTS
	};

	// Four-times oversampling suppresses aliases from the nonlinear limit cycle.
	// The structure-aware split integrator uses cheap adaptive substeps only in
	// the stiff high-damping/high-frequency corner.
	using Oscillator = VdpSplitOscillator<tfdsp::X4Resampler_Order7>;
	std::array<std::unique_ptr<Oscillator>, PORT_MAX_CHANNELS> _vdpHq{};
	int _activeChannels{};

	//----------------------------------------------------------------

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	TfVDPO()
	{
		config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
		configParam(TfVDPO::FREQ, -5.0f, 5.0f, 0.0f, "Frequency offset", " oct");
		configParam(TfVDPO::DAMPING, 0.001f, 9.0f, 0.5f, "Damping");
		configParam(TfVDPO::INPUT_GAIN, 0.0, 1.0f, 1.0f, "Audio input gain", "%", 0.0f, 100.0f);
		configParam(TfVDPO::LEVEL, 0.0, 1.0f, 1.0f, "Output level", "%", 0.0f, 100.0f);
		configParam(TfVDPO::VOCT_SCALING, 0.0f, 1.0f, 1.0f, "1V/octave amount", "%", 0.0f, 100.0f);
		configParam(TfVDPO::DAMPING_ATTENUVERT, -1.0f, 1.0f, 1.0f, "Damping modulation", "%", 0.0f, 100.0f);
		configInput(VOCT_INPUT, "1V/octave pitch");
		configInput(AUDIO_INPUT, "Audio");
		configInput(DAMPING_INPUT, "Damping modulation");
		configOutput(OUTPUT, "Audio");
		for (auto& oscillator : _vdpHq)
			oscillator = std::make_unique<Oscillator>(
				tfdsp::CreateX4Resampler_Cheby7);
		//configParam(TfVDPO::HQ_MODE, -1.0f, 1.0f, -1.0f, "");
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

void TfVDPO::init(float sampleRate)
{
	//_vdp.SetSampleRate(sampleRate);
	for (auto& oscillator : _vdpHq)
		oscillator->SetSampleRate(sampleRate);
}

void TfVDPO::process(const ProcessArgs &args)
{
	const int channels = std::clamp(std::max({
		inputs[VOCT_INPUT].getChannels(),
		inputs[AUDIO_INPUT].getChannels(),
		inputs[DAMPING_INPUT].getChannels(), 1}), 1, PORT_MAX_CHANNELS);
	for (int channel = channels; channel < _activeChannels; ++channel)
		_vdpHq[channel]->Reset();
	_activeChannels = channels;
	outputs[OUTPUT].setChannels(channels);

	const double inputGain = params[INPUT_GAIN].getValue();
	const double pitchScaling = params[VOCT_SCALING].getValue();
	const double pitchOffset = params[FREQ].getValue();
	const double damping = params[DAMPING].getValue();
	const double dampingAmount = params[DAMPING_ATTENUVERT].getValue();
	const double outputLevel = params[LEVEL].getValue();
	const double log2C4AngularFrequency =
		std::log2(2.0 * tfdsp::PI * dsp::FREQ_C4);

	for (int channel = 0; channel < channels; ++channel)
	{
		auto finiteInput = [&](InputIds input)
		{
			const float value = inputs[input].getPolyVoltage(channel);
			return std::isfinite(value) ? static_cast<double>(value) : 0.0;
		};
		const double x = finiteInput(AUDIO_INPUT) * inputGain;
		const double vOct = finiteInput(VOCT_INPUT) * pitchScaling + pitchOffset;
		const double mu = damping + dampingAmount * finiteInput(DAMPING_INPUT);
		const double log2AngularFrequency = log2C4AngularFrequency + vOct;

		// TODO: menu item for low quality; leave high quality by default.
		const double y = _vdpHq[channel]->StepLogAngularFrequency(x, mu,
			log2AngularFrequency);
		const float output = static_cast<float>(y * outputLevel);
		outputs[OUTPUT].setVoltage(std::isfinite(output) ? output : 0.0f,
			channel);
	}
}
void TfVDPO::onReset(const ResetEvent& event)
{
	Module::onReset(event);
	for (auto& oscillator : _vdpHq)
		oscillator->Reset();
	_activeChannels = 0;
}
void TfVDPO::onSampleRateChange(const SampleRateChangeEvent& event)
{
	init(event.sampleRate);
}

struct TfVDPOWidget : ModuleWidget
{
	TfVDPOWidget(TfVDPO *module)
	{
		setModule(module);
		setPanel(APP->window->loadSvg(asset::plugin(pluginInstance, "res/TfVDPO.svg")));

		//Panel screws
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		//Knobs
		addParam(createParam<TfAudioKob>(Vec(14, 58), module, TfVDPO::FREQ));
		addParam(createParam<TfAudioKob>(Vec(14, 112), module, TfVDPO::DAMPING));
		addParam(createParam<TfCvKnob>(Vec(18, 175), module, TfVDPO::INPUT_GAIN));
		addParam(createParam<TfCvKnob>(Vec(76, 175), module, TfVDPO::LEVEL));

		addParam(createParam<TfTrimpot>(Vec(23, 256), module, TfVDPO::VOCT_SCALING));
		addParam(createParam<TfTrimpot>(Vec(81, 256), module, TfVDPO::DAMPING_ATTENUVERT));

		//High quality switch
		//addParam(createParam<ToggleSwitch>(Vec(50, 280), module, TfVDPO::HQ_MODE));

		//Jacks at the bottom
		addInput(createInput<PJ301MPort>(Vec(20, 280), module, TfVDPO::VOCT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(78, 280), module, TfVDPO::DAMPING_INPUT));
		addInput(createInput<PJ301MPort>(Vec(20, 324), module, TfVDPO::AUDIO_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(78, 324), module, TfVDPO::OUTPUT));
	}
};

// Specify the Module and ModuleWidget subclass, human-readable
// author name for categorization per pluginInstance, module slug (should never
// change), human-readable module name, and any number of tags
// (found in `include/tags.hpp`) separated by commas.
Model *modelTfVDPO = createModel<TfVDPO, TfVDPOWidget>("TfVDPO");
