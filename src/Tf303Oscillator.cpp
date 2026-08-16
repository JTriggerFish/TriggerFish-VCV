#include "models/Tb303Oscillator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/control.hpp"
#include "tfdsp/sampleRate.hpp"

struct Tf303Oscillator : Module
{
	enum ParamIds
	{
		OCTAVE,
		TUNE,
		SLIDE_TIME,
		SHAPE,
		WAVE,
		FM_AMOUNT,
		TIME_AMOUNT,
		SHAPE_AMOUNT,
		WAVE_AMOUNT,
		FM_MODE,
		NUM_PARAMS
	};
	enum InputIds
	{
		VOCT_INPUT,
		SLIDE_INPUT,
		TIME_INPUT,
		FM_INPUT,
		SHAPE_INPUT,
		WAVE_INPUT,
		SYNC_INPUT,
		NUM_INPUTS
	};
	enum OutputIds
	{
		CV_OUTPUT,
		AUDIO_OUTPUT,
		NUM_OUTPUTS
	};
	enum LightIds
	{
		NUM_LIGHTS
	};

	using OscillatorX2 = tfdsp::Tb303Oscillator<tfdsp::X2Resampler_Order7>;
	using OscillatorX4 = tfdsp::Tb303Oscillator<tfdsp::X4Resampler_Order7>;
	std::array<std::unique_ptr<OscillatorX2>, PORT_MAX_CHANNELS> oscillatorsX2;
	std::array<std::unique_ptr<OscillatorX4>, PORT_MAX_CHANNELS> oscillatorsX4;
	std::array<dsp::SchmittTrigger, PORT_MAX_CHANNELS> slideTriggers{};
	std::array<tfdsp::FractionalSchmittTrigger,
		PORT_MAX_CHANNELS> syncTriggers{};
	// 4x is the quality default; 2x remains available for dense polyphonic use.
	int oversampling = 1;
	int activeOversampling = 1;

	Tf303Oscillator()
	{
		config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
		configParam(OCTAVE, -3.0f, 3.0f, 0.0f, "Octave", " oct");
		getParamQuantity(OCTAVE)->snapEnabled = true;
		configParam(TUNE, -7.0f / 12.0f, 7.0f / 12.0f, 0.0f,
			"Tune", " semitones", 0.0f, 12.0f);
		configParam(SLIDE_TIME, std::log10(0.002f), std::log10(0.360f),
			std::log10(0.060f), "Slide time", " ms", 10.0f, 1000.0f);
		configParam(SHAPE, -1.0f, 1.0f, 0.0f, "Square shape", "%",
			0.0f, 100.0f);
		configParam(WAVE, 0.0f, 1.0f, 0.0f, "Saw / square morph", "%",
			0.0f, 100.0f);
		configParam(FM_AMOUNT, -1.0f, 1.0f, 0.0f, "FM amount", "%",
			0.0f, 100.0f);
		configParam(TIME_AMOUNT, -1.0f, 1.0f, 0.0f, "Slide-time CV amount",
			"%", 0.0f, 100.0f);
		configParam(SHAPE_AMOUNT, -1.0f, 1.0f, 0.0f, "Shape CV amount",
			"%", 0.0f, 100.0f);
		configParam(WAVE_AMOUNT, -1.0f, 1.0f, 0.0f, "Wave CV amount", "%",
			0.0f, 100.0f);
		configSwitch(FM_MODE, 0.0f, 1.0f, 0.0f, "FM response",
			{"Exponential", "Linear (through-zero)"});

		configInput(VOCT_INPUT, "Pitch (1V/octave)");
		configInput(SLIDE_INPUT, "Slide gate");
		configInput(TIME_INPUT, "Slide time CV");
		configInput(SYNC_INPUT, "Hard sync");
		configInput(FM_INPUT,
			"Frequency modulation (exponential or through-zero linear)");
		configInput(SHAPE_INPUT, "Square shape CV");
		configInput(WAVE_INPUT, "Saw / square morph CV");
		configOutput(CV_OUTPUT, "Post-slide pitch CV");
		configOutput(AUDIO_OUTPUT, "Audio");

		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			oscillatorsX2[channel] = std::make_unique<OscillatorX2>(
				tfdsp::CreateX2Resampler_Chebychev7);
			oscillatorsX4[channel] = std::make_unique<OscillatorX4>(
				tfdsp::CreateX4Resampler_Cheby7);
		}
		SetSampleRate(APP->engine->getSampleRate());
	}

	void SetSampleRate(float sampleRate)
	{
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			oscillatorsX2[channel]->SetSampleRate(sampleRate);
			oscillatorsX4[channel]->SetSampleRate(sampleRate);
		}
	}

	void ResetDsp()
	{
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			oscillatorsX2[channel]->Reset();
			oscillatorsX4[channel]->Reset();
			slideTriggers[channel].reset();
			syncTriggers[channel].Reset();
		}
	}

	void process(const ProcessArgs& args) override
	{
		oversampling = std::clamp(oversampling, 0, 1);
		if (activeOversampling != oversampling)
		{
			activeOversampling = oversampling;
			for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
			{
				if (activeOversampling == 0)
					oscillatorsX2[channel]->Reset();
				else
					oscillatorsX4[channel]->Reset();
			}
		}
		const int channels = std::clamp(std::max({
			inputs[VOCT_INPUT].getChannels(), inputs[SLIDE_INPUT].getChannels(),
			inputs[TIME_INPUT].getChannels(), inputs[FM_INPUT].getChannels(),
			inputs[SYNC_INPUT].getChannels(),
			inputs[SHAPE_INPUT].getChannels(),
			inputs[WAVE_INPUT].getChannels(), 1}), 1, PORT_MAX_CHANNELS);
		outputs[CV_OUTPUT].setChannels(channels);
		outputs[AUDIO_OUTPUT].setChannels(channels);

		const double octave = std::round(params[OCTAVE].getValue());
		const double tuningOffset = octave + params[TUNE].getValue();
		const double slideLog = params[SLIDE_TIME].getValue();
		const double slideRange = std::log10(0.360) - std::log10(0.002);
		const double shapeKnob = params[SHAPE].getValue();
		const double waveKnob = params[WAVE].getValue();
		const double fmAmount = params[FM_AMOUNT].getValue();
		const double timeAmount = params[TIME_AMOUNT].getValue();
		const double shapeAmount = params[SHAPE_AMOUNT].getValue();
		const double waveAmount = params[WAVE_AMOUNT].getValue();
		const bool linearFm = params[FM_MODE].getValue() > 0.5f;

		for (int channel = 0; channel < channels; ++channel)
		{
			auto finiteInput = [&](InputIds input)
			{
				const float value = inputs[input].getPolyVoltage(channel);
				return std::isfinite(value) ? static_cast<double>(value) : 0.0;
			};
			const double timeCv = finiteInput(TIME_INPUT);
			const double channelSlideLog = std::clamp(slideLog + timeAmount *
				(timeCv / 10.0) * slideRange, std::log10(0.002),
				std::log10(0.360));
			const double slideTime = std::pow(10.0, channelSlideLog);
			const double shape = std::clamp(shapeKnob + shapeAmount *
				finiteInput(SHAPE_INPUT) / 5.0, -1.0, 1.0);
			const double wave = std::clamp(waveKnob + waveAmount *
				finiteInput(WAVE_INPUT) / 10.0, 0.0, 1.0);
			const float slideVoltage = static_cast<float>(finiteInput(SLIDE_INPUT));
			slideTriggers[channel].process(slideVoltage, 0.1f, 1.0f);
			const bool slide = slideTriggers[channel].isHigh();
			const auto syncEvent = syncTriggers[channel].Process(
				finiteInput(SYNC_INPUT));
			const double syncCrossing = syncEvent.triggered ?
				syncEvent.position : -1.0;
			float renderedPitch;
			float renderedAudio;
			if (activeOversampling == 0)
			{
				const auto rendered = oscillatorsX2[channel]->Step(
					finiteInput(VOCT_INPUT), slide,
					slideTime, tuningOffset, fmAmount * finiteInput(FM_INPUT),
					linearFm, shape, wave, syncCrossing);
				renderedPitch = rendered.pitch;
				renderedAudio = rendered.mixed;
			}
			else
			{
				const auto rendered = oscillatorsX4[channel]->Step(
					finiteInput(VOCT_INPUT), slide,
					slideTime, tuningOffset, fmAmount * finiteInput(FM_INPUT),
					linearFm, shape, wave, syncCrossing);
				renderedPitch = rendered.pitch;
				renderedAudio = rendered.mixed;
			}
			outputs[CV_OUTPUT].setVoltage(std::clamp(renderedPitch,
				-12.0f, 12.0f), channel);
			outputs[AUDIO_OUTPUT].setVoltage(
				std::isfinite(renderedAudio) ? renderedAudio : 0.0f, channel);
		}
	}

	json_t* dataToJson() override
	{
		json_t* root = json_object();
		json_object_set_new(root, "oversampling", json_integer(oversampling));
		return root;
	}

	void dataFromJson(json_t* root) override
	{
		if (json_t* value = json_object_get(root, "oversampling"))
			oversampling = std::clamp(
				static_cast<int>(json_integer_value(value)), 0, 1);
	}

	void onReset(const ResetEvent& event) override
	{
		Module::onReset(event);
		oversampling = 1;
		activeOversampling = 1;
		ResetDsp();
	}

	void onSampleRateChange(const SampleRateChangeEvent& event) override
	{
		SetSampleRate(event.sampleRate);
	}
};

struct Tf303OscillatorWidget : ModuleWidget
{
	Tf303OscillatorWidget(Tf303Oscillator* module)
	{
		setModule(module);
		setPanel(APP->window->loadSvg(asset::plugin(
			pluginInstance, "res/Tf303Oscillator.svg")));

		auto* logoGraphic = new TfSvgWatermark;
		logoGraphic->setScaledSvg(APP->window->loadSvg(asset::plugin(
			pluginInstance, "res/logo.svg")), Vec(148, 80.8));
		logoGraphic->box.pos = Vec(16, 232);
		logoGraphic->opacity = 0.12f;
		addChild(logoGraphic);

		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		addParam(createParam<TfRotarySwitchKnob>(Vec(20, 49), module,
			Tf303Oscillator::OCTAVE));
		addParam(createParam<TfLargeAudioKnob>(Vec(110, 45), module,
			Tf303Oscillator::TUNE));
		addParam(createParam<TfAudioKob>(Vec(15, 120), module,
			Tf303Oscillator::SLIDE_TIME));
		addParam(createParam<TfAudioKob>(Vec(72, 120), module,
			Tf303Oscillator::SHAPE));
		addParam(createParam<TfAudioKob>(Vec(129, 120), module,
			Tf303Oscillator::WAVE));

		addParam(createParam<TfCvKnob>(Vec(10, 178), module,
			Tf303Oscillator::FM_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(54, 178), module,
			Tf303Oscillator::TIME_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(98, 178), module,
			Tf303Oscillator::SHAPE_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(142, 178), module,
			Tf303Oscillator::WAVE_AMOUNT));
		addParam(createParam<CKSS>(Vec(83, 54), module,
			Tf303Oscillator::FM_MODE));

		addInput(createInput<PJ301MPort>(Vec(12, 241), module,
			Tf303Oscillator::VOCT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(56, 241), module,
			Tf303Oscillator::SLIDE_INPUT));
		addInput(createInput<PJ301MPort>(Vec(100, 241), module,
			Tf303Oscillator::TIME_INPUT));
		addInput(createInput<PJ301MPort>(Vec(144, 241), module,
			Tf303Oscillator::SYNC_INPUT));
		addInput(createInput<PJ301MPort>(Vec(27, 284), module,
			Tf303Oscillator::FM_INPUT));
		addInput(createInput<PJ301MPort>(Vec(78, 284), module,
			Tf303Oscillator::SHAPE_INPUT));
		addInput(createInput<PJ301MPort>(Vec(129, 284), module,
			Tf303Oscillator::WAVE_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(48, 344), module,
			Tf303Oscillator::CV_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(108, 344), module,
			Tf303Oscillator::AUDIO_OUTPUT));
	}

	void appendContextMenu(Menu* menu) override
	{
		Tf303Oscillator* module = dynamic_cast<Tf303Oscillator*>(this->module);
		if (!module)
			return;
		menu->addChild(new MenuSeparator);
		menu->addChild(createIndexPtrSubmenuItem("Oversampling",
			{"2x (lower CPU)", "4x (default)"}, &module->oversampling));
	}
};

Model* modelTf303Oscillator = createModel<Tf303Oscillator,
	Tf303OscillatorWidget>("Tf303Oscillator");
