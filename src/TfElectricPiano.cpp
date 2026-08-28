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
	// Model/Rack boundary gains are kept explicit. The amplifier converts its
	// input to physical volts internally; see ModelInputToCircuitVolts() and the
	// voltage-domain regression in tests/dsp_tests.cpp.
	static constexpr double DirectPickupToRackGain = 20.0;
	static constexpr double PickupSumToAmplifierGain = 5.0;
	static constexpr double AmplifierModelToRackGain = 5.0;

	enum ParamIds
	{
		VELOCITY_CURVE,
		DYNAMICS,
		BODY,
		BELL,
		COUPLING,
		HAMMER,
		TONE,
		PROXIMITY,
		DECAY,
		RELEASE,
		MECHANICS,
		DRIVE,
		OUTPUT_VOLUME,
		AMPLIFIER_BASS,
		AMPLIFIER_TREBLE,
		VIBRATO_SPEED,
		VIBRATO_INTENSITY,
		MODULATION_AMOUNT,
		MODULATION_MODE,
		STRIKE_POSITION,
		NUM_PARAMS
	};

	enum InputIds
	{
		VOCT_INPUT,
		GATE_INPUT,
		VELOCITY_INPUT,
		MODULATION_INPUT,
		RETRIGGER_INPUT,
		PEDAL_INPUT,
		NUM_INPUTS
	};

	enum ModulationModes
	{
		PHASE_MODULATION,
		LINEAR_THROUGH_ZERO_FM,
		EXPONENTIAL_FM
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
	int amplifierQuietSamples = 0;
	bool amplifierSleeping = false;
	int activeChannels{};
	double sampleRate = 48000.0;
	std::uint32_t seedSequence = 1;
	std::array<bool, PORT_MAX_CHANNELS> retriggerHigh{};

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
		configParam(PROXIMITY, 0.0f, 1.0f, 0.48f, "Pickup proximity", "%", 0.0f,
			100.0f);
		configParam(DECAY, 0.0f, 1.0f,
			static_cast<float>(tfdsp::ElectricPianoDefaultDecay),
			"Natural decay", "%", 0.0f, 100.0f);
		configParam(RELEASE, 0.0f, 1.0f,
			static_cast<float>(tfdsp::ElectricPianoDefaultRelease),
			"Damper release", "%", 0.0f, 100.0f);
		configParam(MECHANICS, 0.0f, 1.0f, 0.18f, "Mechanical noise", "%",
			0.0f, 100.0f);
		configParam(DRIVE, 0.0f, 1.0f, 0.32f, "Shared amplifier drive", "%",
			0.0f, 100.0f);
		configParam(OUTPUT_VOLUME, 0.0f, 1.0f, 0.50f,
			"Amplifier output volume", "%", 0.0f, 100.0f);
		configParam(AMPLIFIER_BASS, 0.0f, 1.0f, 0.50f,
			"Peterson preamplifier bass", "%", 0.0f, 100.0f);
		configParam(AMPLIFIER_TREBLE, 0.0f, 1.0f, 0.50f,
			"Peterson preamplifier treble", "%", 0.0f, 100.0f);
		configParam(VIBRATO_SPEED, 0.0f, 1.0f, 0.32f,
			"Suitcase vibrato speed", "%", 0.0f, 100.0f);
		configParam(VIBRATO_INTENSITY, 0.0f, 1.0f, 0.0f,
			"Suitcase vibrato intensity", "%", 0.0f, 100.0f);
		configParam(MODULATION_AMOUNT, -1.0f, 1.0f, 0.0f,
			"Pitch modulation depth", "%", 0.0f, 100.0f);
		// CKSSThree numbers positions from bottom (0) to top (2), hence the
		// apparently reversed enum/label order.
		configSwitch(MODULATION_MODE, PHASE_MODULATION, EXPONENTIAL_FM,
			EXPONENTIAL_FM, "Pitch modulation mode",
			{"Phase", "Linear through-zero FM", "Exponential FM"});
		configParam(STRIKE_POSITION, 0.0f, 1.0f, 0.5f,
			"Hammer strike position (applies on the next strike)", "%", 0.0f,
			100.0f);

		configInput(VOCT_INPUT, "Pitch (1V/octave)");
		configInput(GATE_INPUT, "Key gate");
		configInput(VELOCITY_INPUT, "Strike velocity (0V to 10V)");
		configInput(MODULATION_INPUT,
			"Polyphonic pitch modulation (EXP, linear through-zero FM, or phase)");
		configInput(RETRIGGER_INPUT,
			"Polyphonic note retrigger (for MIDI channel reassignment and legato)");
		configInput(PEDAL_INPUT, "Sustain pedal gate");
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
		amplifierQuietSamples = 0;
		amplifierSleeping = false;
		activeChannels = 0;
		retriggerHigh.fill(false);
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
		controls.proximity = params[PROXIMITY].getValue();
		controls.decay = params[DECAY].getValue();
		controls.release = params[RELEASE].getValue();
		controls.mechanics = params[MECHANICS].getValue();
		controls.drive = params[DRIVE].getValue();
		controls.outputVolume = params[OUTPUT_VOLUME].getValue();
		controls.amplifierBass = params[AMPLIFIER_BASS].getValue();
		controls.amplifierTreble = params[AMPLIFIER_TREBLE].getValue();
		controls.vibratoSpeed = params[VIBRATO_SPEED].getValue();
		controls.vibratoIntensity = params[VIBRATO_INTENSITY].getValue();
		controls.strikePosition = params[STRIKE_POSITION].getValue();
		return controls;
	}

	double FinitePolyVoltage(InputIds input, int channel)
	{
		const float voltage = inputs[input].getPolyVoltage(channel);
		return std::isfinite(voltage) ? static_cast<double>(voltage) : 0.0;
	}

	tfdsp::ElectricPianoModulation PitchModulation(int channel)
	{
		tfdsp::ElectricPianoModulation modulation;
		const double scaledVoltage = params[MODULATION_AMOUNT].getValue() *
			FinitePolyVoltage(MODULATION_INPUT, channel);
		const int mode = static_cast<int>(std::round(
			params[MODULATION_MODE].getValue()));
		if (mode == EXPONENTIAL_FM)
			modulation.exponentialPitch = 0.2 * scaledVoltage;
		else if (mode == LINEAR_THROUGH_ZERO_FM)
			modulation.linearFrequencyRatio = 0.2 * scaledVoltage;
		else
			modulation.phaseRadians = (3.14159265358979323846 / 5.0) *
				scaledVoltage;
		return modulation;
	}

	void process(const ProcessArgs& args) override
	{
		const int channels = std::clamp(std::max({
			inputs[VOCT_INPUT].getChannels(),
			inputs[GATE_INPUT].getChannels(),
			inputs[VELOCITY_INPUT].getChannels(), 1}), 1, PORT_MAX_CHANNELS);
		for (int channel = channels; channel < activeChannels; ++channel)
		{
			PreserveVoiceAsTail(channel);
			retriggerHigh[channel] = false;
		}
		activeChannels = channels;

		const auto baseControls = Controls();
		const bool velocityConnected = inputs[VELOCITY_INPUT].isConnected();
		double pickupSum = 0.0;
		std::array<double, PORT_MAX_CHANNELS> directPickups{};
		for (int channel = 0; channel < channels; ++channel)
		{
			const double velocity = velocityConnected ?
				std::clamp(FinitePolyVoltage(VELOCITY_INPUT, channel) / 10.0,
					0.0, 1.0) :
				0.8;
			const double inputPitch = FinitePolyVoltage(VOCT_INPUT, channel);
			const double gate = FinitePolyVoltage(GATE_INPUT, channel);
			// A disconnected Rack input reads 0 V. This deliberately retains the
			// original gate-only behavior: an upstream CC64-held MIDI gate continues
			// to ring as a held key, while raw gates can use PEDAL to report the key-up
			// and damper-up states independently.
			const bool sustain = FinitePolyVoltage(PEDAL_INPUT, channel) >= 1.0;
			const bool retrigger = FinitePolyVoltage(RETRIGGER_INPUT,
				channel) >= 1.0;
			const bool retriggerEdge = retrigger && !retriggerHigh[channel];
			retriggerHigh[channel] = retrigger;
			const bool newStrike = gate >= 1.0 &&
				(!voices[channel]->GateHigh() || retriggerEdge);
			if (newStrike &&
				voices[channel]->IsAudible() &&
				std::abs(inputPitch - voices[channel]->NotePitch()) > 1.0e-5)
				PreserveVoiceAsTail(channel);
			const double renderedPitch = gate < 1.0 &&
				!voices[channel]->GateHigh() && voices[channel]->IsAudible() ?
				voices[channel]->NotePitch() : inputPitch;
			double pickup = 0.0;
			if (gate >= 1.0 || voices[channel]->GateHigh() ||
				voices[channel]->IsAudible())
			{
				pickup = voices[channel]->Step(renderedPitch, gate,
					velocity, sustain, baseControls, PitchModulation(channel),
					retriggerEdge);
			}
			pickupSum += pickup;
			directPickups[channel] += pickup;
		}

		int directChannels = channels;
		for (auto& tail : tails)
		{
			if (!tail.active)
				continue;
			const bool sustain = FinitePolyVoltage(PEDAL_INPUT,
				tail.outputChannel) >= 1.0;
			const double pickup = tail.voice->Step(tail.voice->NotePitch(),
				0.0, 0.0, sustain, baseControls,
				PitchModulation(tail.outputChannel));
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
				ProcessPostDecimation(DirectPickupToRackGain *
					directPickups[channel]);
			outputs[DIRECT_OUTPUT].setVoltage(
				static_cast<float>(directVoltage), channel);
		}

		outputs[LEFT_OUTPUT].setChannels(1);
		outputs[RIGHT_OUTPUT].setChannels(1);
		bool anyVoiceAudible = false;
		for (int channel = 0; channel < channels; ++channel)
			anyVoiceAudible = anyVoiceAudible || voices[channel]->IsAudible();
		for (const auto& tail : tails)
			anyVoiceAudible = anyVoiceAudible ||
				(tail.active && tail.voice->IsAudible());
		if (!anyVoiceAudible && std::abs(pickupSum) < 1.0e-12)
		{
			++amplifierQuietSamples;
			const int sleepThreshold = std::max(1,
				static_cast<int>(0.25 * sampleRate));
			if (!amplifierSleeping && amplifierQuietSamples >= sleepThreshold)
			{
				amplifier.Reset();
				amplifierSleeping = true;
			}
		}
		else
		{
			amplifierQuietSamples = 0;
			amplifierSleeping = false;
		}
		const auto amplified = amplifierSleeping ?
			std::array<double, 2>{0.0, 0.0} :
			amplifier.Step(PickupSumToAmplifierGain * pickupSum, baseControls);
		const float left = static_cast<float>(tfdsp::RackOutputAdapter::
			ProcessPostDecimation(AmplifierModelToRackGain * amplified[0]));
		const float right = static_cast<float>(tfdsp::RackOutputAdapter::
			ProcessPostDecimation(AmplifierModelToRackGain * amplified[1]));
		outputs[LEFT_OUTPUT].setVoltage(left);
		outputs[RIGHT_OUTPUT].setVoltage(right);
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

		auto* logoGraphic = new TfSvgWatermark;
		logoGraphic->setScaledSvg(APP->window->loadSvg(asset::plugin(
			pluginInstance, "res/logo.svg")), Vec(210.0f, 114.6f));
		logoGraphic->box.pos = Vec(30.0f, 137.0f);
		logoGraphic->opacity = 0.12f;
		addChild(logoGraphic);

		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(
			box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(
			box.size.x - 2 * RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		// The five defining voice controls receive the largest knobs; supporting
		// physical controls, amplifier controls and utilities step down in size.
		addParam(createParam<TfAudioKob>(Vec(9.0f, 48.0f), module,
			TfElectricPiano::DYNAMICS));
		addParam(createParam<TfAudioKob>(Vec(63.0f, 48.0f), module,
			TfElectricPiano::VELOCITY_CURVE));
		addParam(createParam<TfAudioKob>(Vec(117.0f, 48.0f), module,
			TfElectricPiano::HAMMER));
		addParam(createParam<TfAudioKob>(Vec(171.0f, 48.0f), module,
			TfElectricPiano::TONE));
		addParam(createParam<TfAudioKob>(Vec(225.0f, 48.0f), module,
			TfElectricPiano::PROXIMITY));

		addParam(createParam<TfCvKnob>(Vec(4.0f, 103.0f), module,
			TfElectricPiano::BODY));
		addParam(createParam<TfCvKnob>(Vec(42.0f, 103.0f), module,
			TfElectricPiano::BELL));
		addParam(createParam<TfCvKnob>(Vec(80.0f, 103.0f), module,
			TfElectricPiano::COUPLING));
		addParam(createParam<TfCvKnob>(Vec(118.0f, 103.0f), module,
			TfElectricPiano::STRIKE_POSITION));
		addParam(createParam<TfCvKnob>(Vec(156.0f, 103.0f), module,
			TfElectricPiano::DECAY));
		addParam(createParam<TfCvKnob>(Vec(194.0f, 103.0f), module,
			TfElectricPiano::RELEASE));
		addParam(createParam<TfCvKnob>(Vec(232.0f, 103.0f), module,
			TfElectricPiano::MECHANICS));

		addParam(createParam<TfCvKnob>(Vec(8.5f, 173.0f), module,
			TfElectricPiano::DRIVE));
		addParam(createParam<TfCvKnob>(Vec(53.5f, 173.0f), module,
			TfElectricPiano::OUTPUT_VOLUME));
		addParam(createParam<TfCvKnob>(Vec(98.5f, 173.0f), module,
			TfElectricPiano::AMPLIFIER_BASS));
		addParam(createParam<TfCvKnob>(Vec(143.5f, 173.0f), module,
			TfElectricPiano::AMPLIFIER_TREBLE));
		addParam(createParam<TfCvKnob>(Vec(188.5f, 173.0f), module,
			TfElectricPiano::VIBRATO_SPEED));
		addParam(createParam<TfCvKnob>(Vec(233.5f, 173.0f), module,
			TfElectricPiano::VIBRATO_INTENSITY));

		addParam(createParam<CKSSThree>(Vec(38.5f, 241.0f), module,
			TfElectricPiano::MODULATION_MODE));
		addParam(createParam<TfCvKnob>(Vec(126.0f, 241.0f), module,
			TfElectricPiano::MODULATION_AMOUNT));
		addInput(createInput<PJ301MPort>(Vec(213.0f, 243.0f), module,
			TfElectricPiano::MODULATION_INPUT));

		addInput(createInput<PJ301MPort>(Vec(15.0f, 300.0f), module,
			TfElectricPiano::VOCT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(69.0f, 300.0f), module,
			TfElectricPiano::GATE_INPUT));
		addInput(createInput<PJ301MPort>(Vec(123.0f, 300.0f), module,
			TfElectricPiano::RETRIGGER_INPUT));
		addInput(createInput<PJ301MPort>(Vec(177.0f, 300.0f), module,
			TfElectricPiano::VELOCITY_INPUT));
		addInput(createInput<PJ301MPort>(Vec(231.0f, 300.0f), module,
			TfElectricPiano::PEDAL_INPUT));

		addOutput(createOutput<PJ301MPort>(Vec(63.0f, 342.0f), module,
			TfElectricPiano::DIRECT_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(123.0f, 342.0f), module,
			TfElectricPiano::LEFT_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(183.0f, 342.0f), module,
			TfElectricPiano::RIGHT_OUTPUT));
	}
};

Model* modelTfElectricPiano = createModel<TfElectricPiano,
	TfElectricPianoWidget>("TfElectricPiano");
