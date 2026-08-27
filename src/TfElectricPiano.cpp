#include "plugin.hpp"

#include "models/ElectricPiano.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>

#include "components.hpp"
#include "tfdsp/rail.hpp"

struct TfElectricPiano : Module
{
	enum ParamIds
	{
		VELOCITY_CURVE,
		DYNAMICS,
		BODY,
		BELL,
		COUPLING,
		HAMMER,
		TONE,
		GAP,
		DECAY,
		RELEASE,
		MECHANICS,
		DRIVE,
		NUM_PARAMS
	};

	enum InputIds
	{
		VOCT_INPUT,
		GATE_INPUT,
		VELOCITY_INPUT,
		PEDAL_INPUT,
		TONE_INPUT,
		NUM_INPUTS
	};

	enum OutputIds
	{
		DIRECT_OUTPUT,
		LEFT_OUTPUT,
		RIGHT_OUTPUT,
		NUM_OUTPUTS
	};

	enum LightIds
	{
		NUM_LIGHTS
	};

	struct TailVoice
	{
		std::unique_ptr<tfdsp::ElectricPianoVoice> voice;
		int outputChannel{};
		bool active{};
	};

	std::array<std::unique_ptr<tfdsp::ElectricPianoVoice>,
		PORT_MAX_CHANNELS> voices;
	std::array<TailVoice, PORT_MAX_CHANNELS> tails;
	tfdsp::ElectricPianoAmplifier amplifier;
	int activeChannels{};
	double sampleRate = 48000.0;
	std::uint32_t seedSequence = 1;

	TfElectricPiano()
	{
		config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
		configParam(VELOCITY_CURVE, 0.0f, 1.0f, 0.5f, "Velocity curve", "%",
			0.0f, 100.0f);
		configParam(DYNAMICS, 0.0f, 1.0f, 1.0f, "Velocity dynamics", "%",
			0.0f, 100.0f);
		configParam(BODY, 0.0f, 1.0f, 0.62f, "Tine and tone-bar body", "%",
			0.0f, 100.0f);
		configParam(BELL, 0.0f, 1.0f, 0.52f, "Inharmonic bell", "%", 0.0f,
			100.0f);
		configParam(COUPLING, 0.0f, 1.0f, 0.50f,
			"Tine/tone-bar common-base coupling", "%", 0.0f, 100.0f);
		configParam(HAMMER, 0.0f, 1.0f, 0.52f,
			"Hammer hardness (applies on the next strike)", "%", 0.0f, 100.0f);
		configParam(TONE, 0.0f, 1.0f, 0.55f, "Pickup tine alignment", "%",
			0.0f, 100.0f);
		configParam(GAP, 0.0f, 1.0f, 0.48f, "Pickup proximity", "%", 0.0f,
			100.0f);
		configParam(DECAY, 0.0f, 1.0f, 0.50f, "Natural decay", "%", 0.0f,
			100.0f);
		configParam(RELEASE, 0.0f, 1.0f, 0.24f, "Damper release", "%", 0.0f,
			100.0f);
		configParam(MECHANICS, 0.0f, 1.0f, 0.18f, "Mechanical noise", "%",
			0.0f, 100.0f);
		configParam(DRIVE, 0.0f, 1.0f, 0.32f, "Shared amplifier drive", "%",
			0.0f, 100.0f);

		configInput(VOCT_INPUT, "Pitch (1V/octave)");
		configInput(GATE_INPUT, "Key gate");
		configInput(VELOCITY_INPUT, "Strike velocity (0V to 10V)");
		configInput(PEDAL_INPUT, "Sustain pedal gate");
		configInput(TONE_INPUT, "Polyphonic pickup tine-alignment CV");
		configOutput(DIRECT_OUTPUT, "Individual pickup direct polyphonic audio");
		configOutput(LEFT_OUTPUT, "Shared amplifier left audio");
		configOutput(RIGHT_OUTPUT, "Shared amplifier right audio");
		for (auto& voice : voices)
			voice = CreateVoice();
		for (auto& tail : tails)
			tail.voice = CreateVoice();
		SetSampleRate(APP->engine->getSampleRate());
	}

	std::unique_ptr<tfdsp::ElectricPianoVoice> CreateVoice()
	{
		auto voice = std::make_unique<tfdsp::ElectricPianoVoice>();
		voice->SetNoiseSeed(0x6d2b79f5u ^ (0x9e3779b9u * seedSequence++));
		voice->SetSampleRate(sampleRate);
		return voice;
	}

	void SetSampleRate(double sampleRate)
	{
		this->sampleRate = sampleRate;
		for (auto& voice : voices)
			voice->SetSampleRate(sampleRate);
		for (auto& tail : tails)
			tail.voice->SetSampleRate(sampleRate);
		amplifier.SetSampleRate(sampleRate);
	}

	void ResetDsp()
	{
		for (auto& voice : voices)
			voice->Reset();
		for (auto& tail : tails)
		{
			tail.voice->Reset();
			tail.active = false;
		}
		amplifier.Reset();
		activeChannels = 0;
	}

	void PreserveVoiceAsTail(int channel)
	{
		if (!voices[channel]->IsAudible())
		{
			voices[channel]->Reset();
			return;
		}
		std::size_t slot = tails.size();
		double quietestActivity = std::numeric_limits<double>::infinity();
		for (std::size_t index = 0; index < tails.size(); ++index)
		{
			if (!tails[index].active)
			{
				slot = index;
				break;
			}
			if (tails[index].voice->Activity() < quietestActivity)
			{
				quietestActivity = tails[index].voice->Activity();
				slot = index;
			}
		}
		std::swap(tails[slot].voice, voices[channel]);
		tails[slot].outputChannel = channel;
		tails[slot].active = true;
		voices[channel]->Reset();
	}

	tfdsp::ElectricPianoControls Controls()
	{
		tfdsp::ElectricPianoControls controls;
		controls.velocityCurve = params[VELOCITY_CURVE].getValue();
		controls.dynamics = params[DYNAMICS].getValue();
		controls.body = params[BODY].getValue();
		controls.bell = params[BELL].getValue();
		controls.coupling = params[COUPLING].getValue();
		controls.hammer = params[HAMMER].getValue();
		controls.tone = params[TONE].getValue();
		controls.gap = params[GAP].getValue();
		controls.decay = params[DECAY].getValue();
		controls.release = params[RELEASE].getValue();
		controls.mechanics = params[MECHANICS].getValue();
		controls.drive = params[DRIVE].getValue();
		return controls;
	}

	void process(const ProcessArgs& args) override
	{
		const int channels = std::clamp(std::max({
			inputs[VOCT_INPUT].getChannels(),
			inputs[GATE_INPUT].getChannels(),
			inputs[VELOCITY_INPUT].getChannels(), 1}), 1, PORT_MAX_CHANNELS);
		for (int channel = channels; channel < activeChannels; ++channel)
			PreserveVoiceAsTail(channel);
		activeChannels = channels;

		const auto baseControls = Controls();
		const bool velocityConnected = inputs[VELOCITY_INPUT].isConnected();
		double pickupSum = 0.0;
		std::array<double, PORT_MAX_CHANNELS> directPickups{};
		for (int channel = 0; channel < channels; ++channel)
		{
			auto finitePolyVoltage = [&](InputIds input)
			{
				const float voltage = inputs[input].getPolyVoltage(channel);
				return std::isfinite(voltage) ? static_cast<double>(voltage) : 0.0;
			};
			auto controls = baseControls;
			controls.tone = std::clamp(controls.tone +
				0.1 * finitePolyVoltage(TONE_INPUT), 0.0, 1.0);
			const double velocity = velocityConnected ?
				std::clamp(finitePolyVoltage(VELOCITY_INPUT) / 10.0, 0.0, 1.0) :
				0.8;
			const bool sustain = finitePolyVoltage(PEDAL_INPUT) >= 1.0;
			const double inputPitch = finitePolyVoltage(VOCT_INPUT);
			const double gate = finitePolyVoltage(GATE_INPUT);
			if (gate >= 1.0 && !voices[channel]->GateHigh() &&
				voices[channel]->IsAudible() &&
				std::abs(inputPitch - voices[channel]->NotePitch()) > 1.0e-5)
				PreserveVoiceAsTail(channel);
			const double renderedPitch = gate < 1.0 &&
				!voices[channel]->GateHigh() && voices[channel]->IsAudible() ?
				voices[channel]->NotePitch() : inputPitch;
			const double pickup = voices[channel]->Step(
				renderedPitch, gate,
				velocity, sustain, controls);
			pickupSum += pickup;
			directPickups[channel] += pickup;
		}

		int directChannels = channels;
		for (auto& tail : tails)
		{
			if (!tail.active)
				continue;
			auto controls = baseControls;
			const float toneCv = inputs[TONE_INPUT].getPolyVoltage(
				tail.outputChannel);
			controls.tone = std::clamp(controls.tone +
				0.1 * (std::isfinite(toneCv) ? static_cast<double>(toneCv) : 0.0),
				0.0, 1.0);
			const float pedalVoltage = inputs[PEDAL_INPUT].getPolyVoltage(
				tail.outputChannel);
			const bool sustain = std::isfinite(pedalVoltage) &&
				pedalVoltage >= 1.0f;
			const double pickup = tail.voice->Step(tail.voice->NotePitch(),
				0.0, 0.0, sustain, controls);
			pickupSum += pickup;
			directPickups[tail.outputChannel] += pickup;
			directChannels = std::max(directChannels, tail.outputChannel + 1);
			if (!tail.voice->IsAudible())
			{
				tail.voice->Reset();
				tail.active = false;
			}
		}

		outputs[DIRECT_OUTPUT].setChannels(directChannels);
		for (int channel = 0; channel < directChannels; ++channel)
		{
			const double directVoltage = tfdsp::RackOutputAdapter::
				ProcessPostDecimation(20.0 * directPickups[channel]);
			outputs[DIRECT_OUTPUT].setVoltage(
				static_cast<float>(directVoltage), channel);
		}

		outputs[LEFT_OUTPUT].setChannels(1);
		outputs[RIGHT_OUTPUT].setChannels(1);
		const double amplified = amplifier.Step(5.0 * pickupSum,
			baseControls.drive);
		const float output = static_cast<float>(tfdsp::RackOutputAdapter::
			ProcessPostDecimation(5.0 * amplified));
		outputs[LEFT_OUTPUT].setVoltage(output);
		outputs[RIGHT_OUTPUT].setVoltage(output);
	}

	void onReset(const ResetEvent& event) override
	{
		Module::onReset(event);
		ResetDsp();
	}

	void onSampleRateChange(const SampleRateChangeEvent& event) override
	{
		SetSampleRate(event.sampleRate);
	}
};

struct TfElectricPianoWidget : ModuleWidget
{
	TfElectricPianoWidget(TfElectricPiano* module)
	{
		setModule(module);
		setPanel(APP->window->loadSvg(asset::plugin(pluginInstance,
			"res/TfElectricPiano.svg")));

		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(
			box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(
			box.size.x - 2 * RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		addParam(createParam<TfCvKnob>(Vec(16.0f, 48.0f), module,
			TfElectricPiano::VELOCITY_CURVE));
		addParam(createParam<TfCvKnob>(Vec(72.0f, 48.0f), module,
			TfElectricPiano::DYNAMICS));
		addParam(createParam<TfCvKnob>(Vec(128.0f, 48.0f), module,
			TfElectricPiano::BODY));
		addParam(createParam<TfCvKnob>(Vec(184.0f, 48.0f), module,
			TfElectricPiano::BELL));
		addParam(createParam<TfCvKnob>(Vec(16.0f, 105.0f), module,
			TfElectricPiano::COUPLING));
		addParam(createParam<TfCvKnob>(Vec(72.0f, 105.0f), module,
			TfElectricPiano::HAMMER));
		addParam(createParam<TfCvKnob>(Vec(128.0f, 105.0f), module,
			TfElectricPiano::TONE));
		addParam(createParam<TfCvKnob>(Vec(184.0f, 105.0f), module,
			TfElectricPiano::GAP));
		addParam(createParam<TfCvKnob>(Vec(16.0f, 162.0f), module,
			TfElectricPiano::DECAY));
		addParam(createParam<TfCvKnob>(Vec(72.0f, 162.0f), module,
			TfElectricPiano::RELEASE));
		addParam(createParam<TfCvKnob>(Vec(128.0f, 162.0f), module,
			TfElectricPiano::MECHANICS));
		addParam(createParam<TfCvKnob>(Vec(184.0f, 162.0f), module,
			TfElectricPiano::DRIVE));

		addInput(createInput<PJ301MPort>(Vec(16.0f, 259.0f), module,
			TfElectricPiano::VOCT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(62.0f, 259.0f), module,
			TfElectricPiano::GATE_INPUT));
		addInput(createInput<PJ301MPort>(Vec(108.0f, 259.0f), module,
			TfElectricPiano::VELOCITY_INPUT));
		addInput(createInput<PJ301MPort>(Vec(154.0f, 259.0f), module,
			TfElectricPiano::PEDAL_INPUT));
		addInput(createInput<PJ301MPort>(Vec(200.0f, 259.0f), module,
			TfElectricPiano::TONE_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(49.0f, 320.0f), module,
			TfElectricPiano::DIRECT_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(108.0f, 320.0f), module,
			TfElectricPiano::LEFT_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(167.0f, 320.0f), module,
			TfElectricPiano::RIGHT_OUTPUT));
	}
};

Model* modelTfElectricPiano = createModel<TfElectricPiano,
	TfElectricPianoWidget>("TfElectricPiano");
