#include "models/Arp4019Vca.hpp"
#include "models/Arp4072Filter.hpp"
#include "models/ArpEnvelope.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/approx.hpp"
#include "tfdsp/sampleRate.hpp"

struct Tf4072VoiceCore : Module
{
	enum ParamIds
	{
		CUTOFF,
		RESONANCE,
		DRIVE,
		FILTER_ENV_AMOUNT,
		CUTOFF_CV_AMOUNT,
		RES_CV_AMOUNT,
		VCA_INITIAL_GAIN,
		VCA_LINEAR_AMOUNT,
		VCA_EXP_AMOUNT,
		FILTER_ATTACK,
		FILTER_DECAY,
		FILTER_SUSTAIN,
		FILTER_RELEASE,
		FILTER_ENV_MODE,
		AMP_ATTACK,
		AMP_DECAY,
		AMP_SUSTAIN,
		AMP_RELEASE,
		AMP_ENV_MODE,
		FILTER_RANGE,
		NUM_PARAMS
	};

	enum InputIds
	{
		AUDIO_INPUT,
		VCA_AUDIO_INPUT,
		VOCT_INPUT,
		CUTOFF_CV_INPUT,
		RES_CV_INPUT,
		FILTER_ENV_CV_INPUT,
		VCA_LINEAR_CV_INPUT,
		VCA_EXP_CV_INPUT,
		GATE_INPUT,
		TRIGGER_INPUT,
		NUM_INPUTS
	};

	enum OutputIds
	{
		LP_OUTPUT,
		VCA_OUTPUT,
		FILTER_ENV_OUTPUT,
		AMP_ENV_OUTPUT,
		NUM_OUTPUTS
	};

	enum LightIds
	{
		FILTER_ATTACK_LIGHT,
		FILTER_DECAY_LIGHT,
		FILTER_SUSTAIN_LIGHT,
		FILTER_RELEASE_LIGHT,
		AMP_ATTACK_LIGHT,
		AMP_DECAY_LIGHT,
		AMP_SUSTAIN_LIGHT,
		AMP_RELEASE_LIGHT,
		NUM_LIGHTS
	};

	using FilterX2 = tfdsp::Arp4072Filter<tfdsp::X2Resampler_Order7>;
	using FilterX4 = tfdsp::Arp4072Filter<tfdsp::X4Resampler_Order7>;
	using VcaX2 = tfdsp::Arp4019Vca<tfdsp::X2Resampler_Order7>;
	using VcaX4 = tfdsp::Arp4019Vca<tfdsp::X4Resampler_Order7>;

	std::array<std::unique_ptr<FilterX2>, PORT_MAX_CHANNELS> filtersX2;
	std::array<std::unique_ptr<FilterX4>, PORT_MAX_CHANNELS> filtersX4;
	std::array<std::unique_ptr<VcaX2>, PORT_MAX_CHANNELS> vcasX2;
	std::array<std::unique_ptr<VcaX4>, PORT_MAX_CHANNELS> vcasX4;
	std::array<tfdsp::ArpEnvelope, PORT_MAX_CHANNELS> filterEnvelopes{};
	std::array<tfdsp::ArpEnvelope, PORT_MAX_CHANNELS> ampEnvelopes{};
	std::array<float, 4> filterStagePeaks{};
	std::array<float, 4> ampStagePeaks{};
	dsp::ClockDivider lightDivider;

	int oversampling = 1;
	int activeOversampling = 1;
	int activeChannels = 0;

	Tf4072VoiceCore()
	{
		config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
		lightDivider.setDivision(512);
		constexpr float minimumCutoffHz = 10.0f;
		constexpr float defaultCutoffHz = 500.0f;
		constexpr float maximumCutoffHz = 20000.0f;
		const float minimumCutoffPitch = std::log2(
			minimumCutoffHz / dsp::FREQ_C4);
		const float defaultCutoffPitch = std::log2(
			defaultCutoffHz / dsp::FREQ_C4);
		const float maximumCutoffPitch = std::log2(
			maximumCutoffHz / dsp::FREQ_C4);

		configParam(CUTOFF, minimumCutoffPitch, maximumCutoffPitch,
			defaultCutoffPitch, "Initial cutoff", " Hz", 2.0f, dsp::FREQ_C4);
		configParam(RESONANCE, 0.0f, 1.0f, 0.0f, "Resonance", "%", 0.0f,
			100.0f);
		configParam(DRIVE, -60.0f, 24.0f, 0.0f, "Filter input level",
			" dB");
		configParam(FILTER_ENV_AMOUNT, -1.0f, 1.0f, 0.6f,
			"Filter envelope amount", "%", 0.0f, 100.0f);
		configParam(CUTOFF_CV_AMOUNT, -1.0f, 1.0f, 0.0f,
			"Exponential cutoff CV amount", "%", 0.0f, 100.0f);
		configParam(RES_CV_AMOUNT, -1.0f, 1.0f, 0.0f,
			"Resonance CV amount", "%", 0.0f, 100.0f);
		configParam(VCA_INITIAL_GAIN, 0.0f, 1.0f, 0.0f,
			"VCA initial gain", "%", 0.0f, 100.0f);
		configParam(VCA_LINEAR_AMOUNT, 0.0f, 1.0f, 1.0f,
			"Linear VCA CV amount", "%", 0.0f, 100.0f);
		configParam(VCA_EXP_AMOUNT, 0.0f, 1.0f, 1.0f,
			"Exponential VCA CV amount", "%", 0.0f, 100.0f);
		ConfigureEnvelope(FILTER_ATTACK, FILTER_DECAY, FILTER_SUSTAIN,
			FILTER_RELEASE, "Filter", 0.010f, 0.300f, 0.5f, 0.300f);
		ConfigureEnvelope(AMP_ATTACK, AMP_DECAY, AMP_SUSTAIN, AMP_RELEASE,
			"Amplifier", 0.005f, 0.300f, 1.0f, 0.300f);
		configSwitch(FILTER_ENV_MODE, 0.0f, 1.0f, 0.0f,
			"Filter envelope mode", {"ADSR", "AR"});
		configSwitch(AMP_ENV_MODE, 0.0f, 1.0f, 1.0f,
			"Amplifier envelope mode", {"ADSR", "AR"});
		configSwitch(FILTER_RANGE, 0.0f, 1.0f, 0.0f, "Filter range",
			{"Stock 4072", "Extended"});

		configInput(AUDIO_INPUT, "Filter audio");
		configInput(VCA_AUDIO_INPUT,
			"VCA audio override (normalled from filter)");
		configInput(VOCT_INPUT, "Filter cutoff (1V/octave)");
		configInput(CUTOFF_CV_INPUT, "Attenuverted exponential cutoff CV");
		configInput(RES_CV_INPUT, "Resonance CV");
		configInput(FILTER_ENV_CV_INPUT,
			"Filter envelope override (normalled from filter envelope)");
		configInput(VCA_LINEAR_CV_INPUT,
			"Linear VCA CV override (normalled from amplifier envelope)");
		configInput(VCA_EXP_CV_INPUT, "Exponential VCA CV");
		configInput(GATE_INPUT, "Envelope gate");
		configInput(TRIGGER_INPUT, "ADSR retrigger");
		configOutput(LP_OUTPUT, "4072 low-pass");
		configOutput(VCA_OUTPUT, "4019 VCA");
		configOutput(FILTER_ENV_OUTPUT, "Filter envelope");
		configOutput(AMP_ENV_OUTPUT, "Amplifier envelope");
		configLight(FILTER_ATTACK_LIGHT, "Filter envelope attack");
		configLight(FILTER_DECAY_LIGHT, "Filter envelope decay");
		configLight(FILTER_SUSTAIN_LIGHT, "Filter envelope sustain");
		configLight(FILTER_RELEASE_LIGHT, "Filter envelope release");
		configLight(AMP_ATTACK_LIGHT, "Amplifier envelope attack");
		configLight(AMP_DECAY_LIGHT, "Amplifier envelope decay");
		configLight(AMP_SUSTAIN_LIGHT, "Amplifier envelope sustain");
		configLight(AMP_RELEASE_LIGHT, "Amplifier envelope release");
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			filtersX2[channel] = std::make_unique<FilterX2>(
				tfdsp::CreateX2Resampler_Chebychev7);
			filtersX4[channel] = std::make_unique<FilterX4>(
				tfdsp::CreateX4Resampler_Cheby7);
			vcasX2[channel] = std::make_unique<VcaX2>(
				tfdsp::CreateX2Resampler_Chebychev7);
			vcasX4[channel] = std::make_unique<VcaX4>(
				tfdsp::CreateX4Resampler_Cheby7);
		}
		SetSampleRate(APP->engine->getSampleRate());
	}

	void ConfigureEnvelope(int attack, int decay, int sustain, int release,
		const std::string& name, float defaultAttack, float defaultDecay,
		float defaultSustain, float defaultRelease)
	{
		configParam(attack, std::log10(tfdsp::ArpEnvelope::MinimumAttackSeconds),
			std::log10(tfdsp::ArpEnvelope::MaximumAttackSeconds),
			std::log10(defaultAttack), name + " attack", " ms", 10.0f, 1000.0f);
		configParam(decay, std::log10(tfdsp::ArpEnvelope::MinimumDecaySeconds),
			std::log10(tfdsp::ArpEnvelope::MaximumDecaySeconds),
			std::log10(defaultDecay), name + " decay", " ms", 10.0f, 1000.0f);
		configParam(sustain, 0.0f, 1.0f, defaultSustain, name + " sustain",
			"%", 0.0f, 100.0f);
		configParam(release,
			std::log10(tfdsp::ArpEnvelope::MinimumReleaseSeconds),
			std::log10(tfdsp::ArpEnvelope::MaximumReleaseSeconds),
			std::log10(defaultRelease), name + " release", " ms", 10.0f,
			1000.0f);
	}

	void SetSampleRate(float sampleRate)
	{
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			filtersX2[channel]->SetSampleRate(sampleRate);
			filtersX4[channel]->SetSampleRate(sampleRate);
			vcasX2[channel]->SetSampleRate(sampleRate);
			vcasX4[channel]->SetSampleRate(sampleRate);
			filterEnvelopes[channel].SetSampleRate(sampleRate);
			ampEnvelopes[channel].SetSampleRate(sampleRate);
		}
	}

	void ResetDsp()
	{
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
			ResetChannel(channel);
		activeChannels = 0;
		filterStagePeaks.fill(0.0f);
		ampStagePeaks.fill(0.0f);
		lightDivider.reset();
		for (int light = 0; light < NUM_LIGHTS; ++light)
			lights[light].setBrightness(0.0f);
	}

	void ResetChannel(int channel)
	{
		filtersX2[channel]->Reset();
		filtersX4[channel]->Reset();
		vcasX2[channel]->Reset();
		vcasX4[channel]->Reset();
		filterEnvelopes[channel].Reset();
		ampEnvelopes[channel].Reset();
	}

	static int ActiveLightStage(const tfdsp::ArpEnvelope& envelope)
	{
		switch (envelope.GetStage())
		{
		case tfdsp::ArpEnvelope::Stage::Attack:
		case tfdsp::ArpEnvelope::Stage::Hold:
			return 0;
		case tfdsp::ArpEnvelope::Stage::Decay:
			return 1;
		case tfdsp::ArpEnvelope::Stage::Sustain:
			return 2;
		case tfdsp::ArpEnvelope::Stage::Release:
			return 3;
		case tfdsp::ArpEnvelope::Stage::Idle:
			return -1;
		}
		return -1;
	}

	static void CaptureEnvelopeLight(const tfdsp::ArpEnvelope& envelope,
		std::array<float, 4>& stagePeaks)
	{
		const int stage = ActiveLightStage(envelope);
		if (stage < 0)
			return;
		const float level = static_cast<float>(std::clamp(
			envelope.Value(), 0.0, 1.0));
		stagePeaks[stage] = std::max(stagePeaks[stage], level);
	}

	void UpdateEnvelopeLights(const ProcessArgs& args)
	{
		if (!lightDivider.process())
			return;
		const float lightTime = args.sampleTime * lightDivider.getDivision();
		for (int stage = 0; stage < 4; ++stage)
		{
			lights[FILTER_ATTACK_LIGHT + stage].setBrightnessSmooth(
				filterStagePeaks[stage], lightTime);
			lights[AMP_ATTACK_LIGHT + stage].setBrightnessSmooth(
				ampStagePeaks[stage], lightTime);
		}
		filterStagePeaks.fill(0.0f);
		ampStagePeaks.fill(0.0f);
	}

	void ResetActivePath()
	{
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			if (activeOversampling == 0)
			{
				filtersX2[channel]->Reset();
				vcasX2[channel]->Reset();
			}
			else
			{
				filtersX4[channel]->Reset();
				vcasX4[channel]->Reset();
			}
		}
	}

	void process(const ProcessArgs& args) override
	{
		oversampling = std::clamp(oversampling, 0, 1);
		if (oversampling != activeOversampling)
		{
			activeOversampling = oversampling;
			ResetActivePath();
		}

		int channels = 1;
		for (int input = 0; input < NUM_INPUTS; ++input)
			channels = std::max(channels, inputs[input].getChannels());
		for (int channel = channels; channel < activeChannels; ++channel)
			ResetChannel(channel);
		activeChannels = channels;
		for (auto output : {LP_OUTPUT, VCA_OUTPUT, FILTER_ENV_OUTPUT,
			AMP_ENV_OUTPUT})
			outputs[output].setChannels(channels);

		const double cutoffKnob = params[CUTOFF].getValue();
		const double resonanceKnob = params[RESONANCE].getValue();
		const double driveDb = params[DRIVE].getValue();
		const double driveGain = driveDb <= -59.99 ? 0.0 :
			std::pow(10.0, driveDb / 20.0);
		const double filterEnvAmount = params[FILTER_ENV_AMOUNT].getValue();
		const double cutoffCvAmount = params[CUTOFF_CV_AMOUNT].getValue();
		const double resonanceCvAmount = params[RES_CV_AMOUNT].getValue();
		const double initialGain = params[VCA_INITIAL_GAIN].getValue();
		const double linearAmount = params[VCA_LINEAR_AMOUNT].getValue();
		const double exponentialAmount = params[VCA_EXP_AMOUNT].getValue();
		const bool extendedCutoff = params[FILTER_RANGE].getValue() > 0.5f;
		const bool vcaOverride = inputs[VCA_AUDIO_INPUT].isConnected();

		const double filterAttack = std::pow(10.0,
			static_cast<double>(params[FILTER_ATTACK].getValue()));
		const double filterDecay = std::pow(10.0,
			static_cast<double>(params[FILTER_DECAY].getValue()));
		const double filterSustain = params[FILTER_SUSTAIN].getValue();
		const double filterRelease = std::pow(10.0,
			static_cast<double>(params[FILTER_RELEASE].getValue()));
		const double ampAttack = std::pow(10.0,
			static_cast<double>(params[AMP_ATTACK].getValue()));
		const double ampDecay = std::pow(10.0,
			static_cast<double>(params[AMP_DECAY].getValue()));
		const double ampSustain = params[AMP_SUSTAIN].getValue();
		const double ampRelease = std::pow(10.0,
			static_cast<double>(params[AMP_RELEASE].getValue()));

		for (int channel = 0; channel < channels; ++channel)
		{
			const double gate = inputs[GATE_INPUT].getPolyVoltage(channel);
			const double trigger = inputs[TRIGGER_INPUT].getPolyVoltage(channel);
			filterEnvelopes[channel].SetMode(
				params[FILTER_ENV_MODE].getValue() > 0.5f ?
					tfdsp::ArpEnvelope::Mode::Ar : tfdsp::ArpEnvelope::Mode::Adsr);
			ampEnvelopes[channel].SetMode(
				params[AMP_ENV_MODE].getValue() > 0.5f ?
					tfdsp::ArpEnvelope::Mode::Ar : tfdsp::ArpEnvelope::Mode::Adsr);
			const double filterEnvelope = 10.0 * filterEnvelopes[channel].Step(
				gate, trigger, filterAttack, filterDecay, filterSustain,
				filterRelease);
			const double ampEnvelope = 10.0 * ampEnvelopes[channel].Step(gate,
				trigger, ampAttack, ampDecay, ampSustain, ampRelease);
			CaptureEnvelopeLight(filterEnvelopes[channel], filterStagePeaks);
			CaptureEnvelopeLight(ampEnvelopes[channel], ampStagePeaks);
			outputs[FILTER_ENV_OUTPUT].setVoltage(filterEnvelope, channel);
			outputs[AMP_ENV_OUTPUT].setVoltage(ampEnvelope, channel);

			const double filterModulation = inputs[FILTER_ENV_CV_INPUT].isConnected() ?
				inputs[FILTER_ENV_CV_INPUT].getPolyVoltage(channel) : filterEnvelope;
			const double linearControl = linearAmount *
				(inputs[VCA_LINEAR_CV_INPUT].isConnected() ?
					inputs[VCA_LINEAR_CV_INPUT].getPolyVoltage(channel) : ampEnvelope);
			const double exponentialControl = exponentialAmount *
				inputs[VCA_EXP_CV_INPUT].getPolyVoltage(channel);
			const double pitch = cutoffKnob +
				inputs[VOCT_INPUT].getPolyVoltage(channel) + cutoffCvAmount *
				inputs[CUTOFF_CV_INPUT].getPolyVoltage(channel) +
				0.5 * filterEnvAmount * filterModulation;
			const double cutoff = dsp::FREQ_C4 * tfdsp::Exp2Taylor5(
				static_cast<float>(std::clamp(pitch, -10.0, 10.0)));
			const double resonance = resonanceKnob + resonanceCvAmount *
				inputs[RES_CV_INPUT].getPolyVoltage(channel) / 10.0;
			const double audio = inputs[AUDIO_INPUT].getPolyVoltage(channel);

			float lowPass = 0.0f;
			float vcaOutput = 0.0f;
			if (activeOversampling == 0)
			{
				if (vcaOverride)
				{
					lowPass = filtersX2[channel]->Step(audio, cutoff, resonance,
						driveGain, extendedCutoff);
					vcaOutput = vcasX2[channel]->Step(
						inputs[VCA_AUDIO_INPUT].getPolyVoltage(channel), 0.0,
						linearControl, exponentialControl, initialGain);
				}
				else
				{
					const auto rendered = filtersX2[channel]->StepWithPostProcessor(
						audio, cutoff, resonance, driveGain, extendedCutoff,
						linearControl, exponentialControl,
						[&](double filtered, double linearCv, double exponentialCv)
						{
							return vcasX2[channel]->ProcessOversampled(filtered,
								linearCv, exponentialCv, initialGain);
						});
					lowPass = rendered.lowPass;
					vcaOutput = rendered.postProcessed;
				}
			}
			else
			{
				if (vcaOverride)
				{
					lowPass = filtersX4[channel]->Step(audio, cutoff, resonance,
						driveGain, extendedCutoff);
					vcaOutput = vcasX4[channel]->Step(
						inputs[VCA_AUDIO_INPUT].getPolyVoltage(channel), 0.0,
						linearControl, exponentialControl, initialGain);
				}
				else
				{
					const auto rendered = filtersX4[channel]->StepWithPostProcessor(
						audio, cutoff, resonance, driveGain, extendedCutoff,
						linearControl, exponentialControl,
						[&](double filtered, double linearCv, double exponentialCv)
						{
							return vcasX4[channel]->ProcessOversampled(filtered,
								linearCv, exponentialCv, initialGain);
						});
					lowPass = rendered.lowPass;
					vcaOutput = rendered.postProcessed;
				}
			}
			outputs[LP_OUTPUT].setVoltage(lowPass, channel);
			outputs[VCA_OUTPUT].setVoltage(vcaOutput, channel);
		}
		UpdateEnvelopeLights(args);
	}

	void processBypass(const ProcessArgs& args) override
	{
		(void) args;
		auto route = [&](int inputId, int outputId)
		{
			Input& input = inputs[inputId];
			Output& output = outputs[outputId];
			const int channels = input.getChannels();
			for (int channel = 0; channel < channels; ++channel)
				output.setVoltage(input.getVoltage(channel), channel);
			output.setChannels(channels);
		};

		route(AUDIO_INPUT, LP_OUTPUT);
		route(inputs[VCA_AUDIO_INPUT].isConnected() ? VCA_AUDIO_INPUT :
			AUDIO_INPUT, VCA_OUTPUT);
		outputs[FILTER_ENV_OUTPUT].setChannels(0);
		outputs[AMP_ENV_OUTPUT].setChannels(0);
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

struct Tf4072VoiceCoreWidget : ModuleWidget
{
	Tf4072VoiceCoreWidget(Tf4072VoiceCore* module)
	{
		setModule(module);
		setPanel(APP->window->loadSvg(asset::plugin(
			pluginInstance, "res/Tf4072VoiceCore.svg")));

		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(
			Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(
			Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		addParam(createParam<TfLargeAudioKnob>(Vec(7, 40), module,
			Tf4072VoiceCore::CUTOFF));
		addParam(createParam<TfLargeAudioKnob>(Vec(99, 40), module,
			Tf4072VoiceCore::RESONANCE));
		addParam(createParam<TfAudioKob>(Vec(62, 112), module,
			Tf4072VoiceCore::DRIVE));
		addParam(createParam<TfCvKnob>(Vec(16, 173), module,
			Tf4072VoiceCore::FILTER_ENV_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(66, 173), module,
			Tf4072VoiceCore::CUTOFF_CV_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(116, 173), module,
			Tf4072VoiceCore::RES_CV_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(16, 231), module,
			Tf4072VoiceCore::VCA_INITIAL_GAIN));
		addParam(createParam<TfCvKnob>(Vec(66, 231), module,
			Tf4072VoiceCore::VCA_LINEAR_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(116, 231), module,
			Tf4072VoiceCore::VCA_EXP_AMOUNT));
		addParam(createParam<CKSS>(Vec(73, 65), module,
			Tf4072VoiceCore::FILTER_RANGE));

		addParam(createLightParam<TfEnvelopeSlider>(Vec(174, 52), module,
			Tf4072VoiceCore::FILTER_ATTACK,
			Tf4072VoiceCore::FILTER_ATTACK_LIGHT));
		addParam(createLightParam<TfEnvelopeSlider>(Vec(214, 52), module,
			Tf4072VoiceCore::FILTER_DECAY,
			Tf4072VoiceCore::FILTER_DECAY_LIGHT));
		addParam(createLightParam<TfEnvelopeSlider>(Vec(254, 52), module,
			Tf4072VoiceCore::FILTER_SUSTAIN,
			Tf4072VoiceCore::FILTER_SUSTAIN_LIGHT));
		addParam(createLightParam<TfEnvelopeSlider>(Vec(294, 52), module,
			Tf4072VoiceCore::FILTER_RELEASE,
			Tf4072VoiceCore::FILTER_RELEASE_LIGHT));
		addParam(createLightParam<TfEnvelopeSlider>(Vec(174, 177), module,
			Tf4072VoiceCore::AMP_ATTACK,
			Tf4072VoiceCore::AMP_ATTACK_LIGHT));
		addParam(createLightParam<TfEnvelopeSlider>(Vec(214, 177), module,
			Tf4072VoiceCore::AMP_DECAY,
			Tf4072VoiceCore::AMP_DECAY_LIGHT));
		addParam(createLightParam<TfEnvelopeSlider>(Vec(254, 177), module,
			Tf4072VoiceCore::AMP_SUSTAIN,
			Tf4072VoiceCore::AMP_SUSTAIN_LIGHT));
		addParam(createLightParam<TfEnvelopeSlider>(Vec(294, 177), module,
			Tf4072VoiceCore::AMP_RELEASE,
			Tf4072VoiceCore::AMP_RELEASE_LIGHT));
		addParam(createParam<CKSS>(Vec(329, 78), module,
			Tf4072VoiceCore::FILTER_ENV_MODE));
		addParam(createParam<CKSS>(Vec(329, 203), module,
			Tf4072VoiceCore::AMP_ENV_MODE));

		addInput(createInput<PJ301MPort>(Vec(18, 293), module,
			Tf4072VoiceCore::VOCT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(68, 293), module,
			Tf4072VoiceCore::CUTOFF_CV_INPUT));
		addInput(createInput<PJ301MPort>(Vec(118, 293), module,
			Tf4072VoiceCore::RES_CV_INPUT));
		addInput(createInput<PJ301MPort>(Vec(168, 293), module,
			Tf4072VoiceCore::GATE_INPUT));
		addInput(createInput<PJ301MPort>(Vec(218, 293), module,
			Tf4072VoiceCore::TRIGGER_INPUT));
		addInput(createInput<PJ301MPort>(Vec(268, 293), module,
			Tf4072VoiceCore::FILTER_ENV_CV_INPUT));
		addInput(createInput<PJ301MPort>(Vec(318, 293), module,
			Tf4072VoiceCore::AUDIO_INPUT));
		addInput(createInput<PJ301MPort>(Vec(18, 334), module,
			Tf4072VoiceCore::VCA_AUDIO_INPUT));
		addInput(createInput<PJ301MPort>(Vec(68, 334), module,
			Tf4072VoiceCore::VCA_LINEAR_CV_INPUT));
		addInput(createInput<PJ301MPort>(Vec(118, 334), module,
			Tf4072VoiceCore::VCA_EXP_CV_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(168, 334), module,
			Tf4072VoiceCore::LP_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(218, 334), module,
			Tf4072VoiceCore::VCA_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(268, 334), module,
			Tf4072VoiceCore::FILTER_ENV_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(318, 334), module,
			Tf4072VoiceCore::AMP_ENV_OUTPUT));
	}

	void appendContextMenu(Menu* menu) override
	{
		Tf4072VoiceCore* module = dynamic_cast<Tf4072VoiceCore*>(this->module);
		if (!module)
			return;
		menu->addChild(new MenuSeparator);
		menu->addChild(createIndexPtrSubmenuItem("Oversampling",
			{"2x (lower CPU)", "4x (default)"}, &module->oversampling));
	}
};

Model* modelTf4072VoiceCore = createModel<Tf4072VoiceCore,
	Tf4072VoiceCoreWidget>("Tf4072VoiceCore");
