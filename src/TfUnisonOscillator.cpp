#include <algorithm>
#include <array>
#include <cmath>
#include <random>

#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/approx.hpp"
#include "tfdsp/control.hpp"
#include "tfdsp/rail.hpp"
#include "tfdsp/unison.hpp"
#include "tfdsp/unison_oscillator.hpp"

namespace
{
	constexpr double MinimumDriftTimeSeconds = 0.5;
	constexpr double MaximumDriftTimeSeconds = 120.0;
	constexpr double MaximumCommonDriftCents = 10.0;
	constexpr double MaximumIndividualDriftCents = 10.0;
	constexpr double MaximumIndividualDriftHz = 1.5;
	constexpr double MaximumTrackingErrorCentsPerOctave = 20.0;
	constexpr double MaximumHumCents = 12.0;
	constexpr double MinimumPwmRateHz = 0.03;
	constexpr double MaximumPwmRateHz = 10.0;

	double DriftTimeSeconds(double speed)
	{
		const double shaped = std::pow(std::clamp(speed, 0.0, 1.0), 2.0);
		return MaximumDriftTimeSeconds * std::pow(
			MinimumDriftTimeSeconds / MaximumDriftTimeSeconds, shaped);
	}

	double PwmRateHz(double control)
	{
		return MinimumPwmRateHz * std::pow(MaximumPwmRateHz /
			MinimumPwmRateHz, std::clamp(control, 0.0, 1.0));
	}

	struct DriftSpeedQuantity : engine::ParamQuantity
	{
		float getDisplayValue() override
		{
			return static_cast<float>(DriftTimeSeconds(getValue()));
		}

		void setDisplayValue(float displayValue) override
		{
			const double seconds = std::clamp(static_cast<double>(displayValue),
				MinimumDriftTimeSeconds, MaximumDriftTimeSeconds);
			const double exponent = std::log(seconds / MaximumDriftTimeSeconds) /
				std::log(MinimumDriftTimeSeconds / MaximumDriftTimeSeconds);
			setValue(static_cast<float>(std::sqrt(
				std::clamp(exponent, 0.0, 1.0))));
		}
	};

	struct StackSpreadQuantity : engine::ParamQuantity
	{
		float getDisplayValue() override
		{
			return static_cast<float>(tfdsp::UnisonSpreadCents(getValue()));
		}

		void setDisplayValue(float displayValue) override
		{
			setValue(static_cast<float>(tfdsp::UnisonSpreadControlForCents(
				static_cast<double>(displayValue))));
		}
	};

	struct PwmRateQuantity : engine::ParamQuantity
	{
		float getDisplayValue() override
		{
			return static_cast<float>(PwmRateHz(getValue()));
		}

		void setDisplayValue(float displayValue) override
		{
			const double rate = std::clamp(static_cast<double>(displayValue),
				MinimumPwmRateHz, MaximumPwmRateHz);
			setValue(static_cast<float>(std::log(rate / MinimumPwmRateHz) /
				std::log(MaximumPwmRateHz / MinimumPwmRateHz)));
		}
	};
}

struct TfUnisonOscillator : Module
{
	enum ParamIds
	{
		OCTAVE,
		TUNE,
		VOICES,
		WAVEFORM,
		PULSE_WIDTH,
		PULSE_WIDTH_CV_AMOUNT,
		PWM_RATE,
		SPREAD,
		WIDTH,
		SUB_LEVEL,
		SUB_MODE,
		DRIFT_SPEED,
		HUM,
		COMMON_DRIFT,
		INDIVIDUAL_DRIFT,
		TRACKING,
		INDIVIDUAL_DRIFT_MODE,
		SPREAD_CV_AMOUNT,
		WIDTH_CV_AMOUNT,
		NUM_PARAMS
	};

	enum InputIds
	{
		VOCT_INPUT,
		PULSE_WIDTH_INPUT,
		SPREAD_INPUT,
		WIDTH_INPUT,
		NUM_INPUTS
	};

	enum OutputIds
	{
		MONO_OUTPUT,
		LEFT_OUTPUT,
		RIGHT_OUTPUT,
		NUM_OUTPUTS
	};

	enum LightIds
	{
		NUM_LIGHTS
	};

	std::array<std::array<tfdsp::StackedOscillatorVoice,
		tfdsp::MaximumStackedOscillatorVoices>, PORT_MAX_CHANNELS> voices{};
	std::array<tfdsp::BandlimitedFixedPulseOscillator,
		PORT_MAX_CHANNELS> centerSubs{};
	std::array<std::array<tfdsp::SmoothOrnsteinUhlenbeck,
		tfdsp::MaximumStackedOscillatorVoices>, PORT_MAX_CHANNELS>
		individualDriftProcesses{};
	tfdsp::SmoothOrnsteinUhlenbeck commonDriftProcess{};
	tfdsp::RecursiveSineOscillator humOscillator{};
	tfdsp::RecursiveSineOscillator pwmOscillator{};
	std::random_device randomSeed{};
	std::minstd_rand randomGenerator;

	std::array<double, tfdsp::MaximumStackedOscillatorVoices> voiceGains{};
	std::array<double, tfdsp::MaximumStackedOscillatorVoices> pitchPositions{};
	std::array<double, tfdsp::MaximumStackedOscillatorVoices> panPositions{};
	std::array<double, tfdsp::MaximumStackedOscillatorVoices> trackingPositions{};
	std::array<double, tfdsp::MaximumStackedOscillatorVoices> targetPitchPositions{};
	std::array<double, tfdsp::MaximumStackedOscillatorVoices> targetPanPositions{};
	std::array<double, tfdsp::MaximumStackedOscillatorVoices> targetTrackingPositions{};
	int layoutVoiceCount{7};
	double waveformMix{};
	double subModeMix{};
	double sampleRate{48000.0};
	double driftTimeSeconds{};
	double configuredPwmRateHz{};
	double transitionCoefficient{1.0};

	TfUnisonOscillator() : randomGenerator(randomSeed())
	{
		config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
		configParam(OCTAVE, -3.0f, 3.0f, 0.0f, "Octave", " oct");
		getParamQuantity(OCTAVE)->snapEnabled = true;
		configParam(TUNE, -7.0f / 12.0f, 7.0f / 12.0f, 0.0f,
			"Tune", " semitones", 0.0f, 12.0f);
		configParam(VOICES, 1.0f, 16.0f, 3.0f, "Unison voices");
		getParamQuantity(VOICES)->snapEnabled = true;
		configSwitch(WAVEFORM, 0.0f, 1.0f, 0.0f, "Waveform",
			{"Saw", "Pulse"});
		configParam(PULSE_WIDTH, 0.05f, 0.95f, 0.5f,
			"Pulse width", "%", 0.0f, 100.0f);
		configParam(PULSE_WIDTH_CV_AMOUNT, -1.0f, 1.0f, 0.0f,
			"PWM amount", "%", 0.0f, 100.0f);
		configParam<PwmRateQuantity>(PWM_RATE, 0.0f, 1.0f, 0.42f,
			"Internal PWM rate", " Hz");
		configParam<StackSpreadQuantity>(SPREAD, 0.0f, 1.0f,
			static_cast<float>(tfdsp::UnisonSpreadControlForCents(4.0)),
			"Unison spread", " cents");
		configParam(WIDTH, 0.0f, 1.0f, 0.65f,
			"Stereo width", "%", 0.0f, 100.0f);
		configParam(SUB_LEVEL, 0.0f, 1.0f, 0.0f,
			"Sub level in main outputs", "%", 0.0f, 100.0f);
		configSwitch(SUB_MODE, 0.0f, 1.0f, 0.0f, "Sub oscillator mode",
			{"Centre", "Stack"});
		configParam<DriftSpeedQuantity>(DRIFT_SPEED, 0.0f, 1.0f, 0.5f,
			"Drift time", " s");
		configParam(HUM, 0.0f, 1.0f, 0.1f,
			"Common 60 Hz hum", " cents peak", 0.0f, MaximumHumCents);
		configParam(COMMON_DRIFT, 0.0f, 1.0f, 0.1f,
			"Common drift", "%", 0.0f, 100.0f);
		configParam(INDIVIDUAL_DRIFT, 0.0f, 1.0f, 0.15f,
			"Individual drift", "%", 0.0f, 100.0f);
		configParam(TRACKING, 0.0f, 1.0f, 1.0f,
			"Tracking", "%", 0.0f, 100.0f);
		configSwitch(INDIVIDUAL_DRIFT_MODE, 0.0f, 1.0f, 1.0f,
			"Individual drift units", {"Hertz", "Cents"});
		configParam(SPREAD_CV_AMOUNT, -1.0f, 1.0f, 0.0f,
			"Spread CV amount", "%", 0.0f, 100.0f);
		configParam(WIDTH_CV_AMOUNT, -1.0f, 1.0f, 0.0f,
			"Width CV amount", "%", 0.0f, 100.0f);

		configInput(VOCT_INPUT, "Pitch (1V/octave)");
		configInput(PULSE_WIDTH_INPUT, "Pulse-width CV");
		configInput(SPREAD_INPUT, "Spread CV");
		configInput(WIDTH_INPUT, "Stereo-width CV");
		configOutput(MONO_OUTPUT, "Mono unison mix");
		configOutput(LEFT_OUTPUT, "Left unison mix");
		configOutput(RIGHT_OUTPUT, "Right unison mix");

		SetLayout(layoutVoiceCount, true);
		SetSampleRate(APP->engine->getSampleRate());
		ResetDsp();
	}

	void SetLayout(int count, bool immediate = false)
	{
		layoutVoiceCount = std::clamp(count, 1,
			tfdsp::MaximumStackedOscillatorVoices);
		targetPitchPositions =
			tfdsp::StackedOscillatorPitchPositions(layoutVoiceCount);
		targetPanPositions =
			tfdsp::StackedOscillatorPanPositions(layoutVoiceCount);
		targetTrackingPositions =
			tfdsp::StackedOscillatorTrackingPositions(layoutVoiceCount);
		if (immediate)
		{
			pitchPositions = targetPitchPositions;
			panPositions = targetPanPositions;
			trackingPositions = targetTrackingPositions;
			for (int voice = 0; voice < tfdsp::MaximumStackedOscillatorVoices;
				++voice)
				voiceGains[voice] = voice < layoutVoiceCount ? 1.0 : 0.0;
		}
	}

	void ConfigureDrift(double timeSeconds)
	{
		driftTimeSeconds = timeSeconds;
		commonDriftProcess.ConfigureStationary(sampleRate, timeSeconds, 1.0,
			100.0);
		for (auto& channel : individualDriftProcesses)
			for (auto& process : channel)
				process.ConfigureStationary(sampleRate, timeSeconds, 1.0, 100.0);
	}

	void SetSampleRate(double nextSampleRate)
	{
		sampleRate = std::max(nextSampleRate, 1.0);
		humOscillator.SetFrequency(60.0, sampleRate);
		configuredPwmRateHz = PwmRateHz(params[PWM_RATE].getValue());
		pwmOscillator.SetFrequency(configuredPwmRateHz, sampleRate);
		transitionCoefficient = -std::expm1(-1.0 / (0.010 * sampleRate));
		ConfigureDrift(DriftTimeSeconds(params[DRIFT_SPEED].getValue()));
	}

	void ResetDsp()
	{
		constexpr double GoldenConjugate = 0.6180339887498948482;
		constexpr double ChannelOffset = 0.3819660112501051518;
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
			for (int voice = 0; voice < tfdsp::MaximumStackedOscillatorVoices;
				++voice)
				voices[channel][voice].Reset(std::fmod(
					voice * GoldenConjugate + channel * ChannelOffset, 1.0));
		for (auto& oscillator : centerSubs)
			oscillator.Reset();
		commonDriftProcess.Reset();
		for (auto& channel : individualDriftProcesses)
			for (auto& process : channel)
				process.Reset();
		humOscillator.Reset();
		pwmOscillator.Reset();
		SetLayout(static_cast<int>(std::round(params[VOICES].getValue())), true);
		waveformMix = params[WAVEFORM].getValue();
		subModeMix = params[SUB_MODE].getValue();
	}

	void UpdateTransitions(int requestedVoices)
	{
		if (requestedVoices != layoutVoiceCount)
			SetLayout(requestedVoices);
		for (int voice = 0; voice < tfdsp::MaximumStackedOscillatorVoices;
			++voice)
		{
			const double gainTarget = voice < layoutVoiceCount ? 1.0 : 0.0;
			voiceGains[voice] += transitionCoefficient *
				(gainTarget - voiceGains[voice]);
			pitchPositions[voice] += transitionCoefficient *
				(targetPitchPositions[voice] - pitchPositions[voice]);
			panPositions[voice] += transitionCoefficient *
				(targetPanPositions[voice] - panPositions[voice]);
			trackingPositions[voice] += transitionCoefficient *
				(targetTrackingPositions[voice] - trackingPositions[voice]);
		}
		waveformMix = params[WAVEFORM].getValue() > 0.5f ? 1.0 : 0.0;
		subModeMix = params[SUB_MODE].getValue() > 0.5f ? 1.0 : 0.0;
	}

	void process(const ProcessArgs& args) override
	{
		const int requestedVoices = std::clamp(static_cast<int>(std::round(
			params[VOICES].getValue())), 1,
			tfdsp::MaximumStackedOscillatorVoices);
		UpdateTransitions(requestedVoices);

		const double requestedDriftTime = DriftTimeSeconds(
			params[DRIFT_SPEED].getValue());
		if (std::abs(requestedDriftTime - driftTimeSeconds) > 1.0e-9)
			ConfigureDrift(requestedDriftTime);
		const double requestedPwmRate = PwmRateHz(params[PWM_RATE].getValue());
		if (std::abs(requestedPwmRate - configuredPwmRateHz) > 1.0e-9)
		{
			configuredPwmRateHz = requestedPwmRate;
			pwmOscillator.SetFrequency(configuredPwmRateHz, sampleRate);
		}

		const int channels = std::clamp(std::max({
			inputs[VOCT_INPUT].getChannels(),
			inputs[PULSE_WIDTH_INPUT].getChannels(),
			inputs[SPREAD_INPUT].getChannels(),
			inputs[WIDTH_INPUT].getChannels(), 1}), 1, PORT_MAX_CHANNELS);
		for (auto output : {MONO_OUTPUT, LEFT_OUTPUT, RIGHT_OUTPUT})
			outputs[output].setChannels(channels);

		const bool needMono = outputs[MONO_OUTPUT].isConnected();
		const bool needLeft = outputs[LEFT_OUTPUT].isConnected();
		const bool needRight = outputs[RIGHT_OUTPUT].isConnected();
		const bool needMain = needMono || needLeft || needRight;
		const double subLevel = params[SUB_LEVEL].getValue();
		const bool needSub = needMain && subLevel > 1.0e-7;
		const bool renderCenterSub = needSub && subModeMix == 0.0;
		const bool renderStackSub = needSub && subModeMix == 1.0;

		const double commonDriftCents = MaximumCommonDriftCents *
			params[COMMON_DRIFT].getValue() *
			commonDriftProcess.Step(randomGenerator);
		const double humCents = MaximumHumCents * params[HUM].getValue() *
			humOscillator.Step();
		const double commonPitchCents = commonDriftCents + humCents;
		const double internalPwmVoltage = 5.0 * pwmOscillator.Step();
		const double trackingErrorDepth = 1.0 - params[TRACKING].getValue();
		const double individualDepth = params[INDIVIDUAL_DRIFT].getValue();
		const bool centsDrift = params[INDIVIDUAL_DRIFT_MODE].getValue() > 0.5f;
		const double pulseWidthBase = params[PULSE_WIDTH].getValue();
		const double pulseWidthAmount = params[PULSE_WIDTH_CV_AMOUNT].getValue();
		const double spreadCvAmount = params[SPREAD_CV_AMOUNT].getValue();
		const double widthCvAmount = params[WIDTH_CV_AMOUNT].getValue();
		const double tune = std::round(params[OCTAVE].getValue()) +
			params[TUNE].getValue();

		int voicesToProcess = 0;
		double sumGainSquares = 0.0;
		for (int voice = 0; voice < tfdsp::MaximumStackedOscillatorVoices;
			++voice)
		{
			if (voiceGains[voice] > 1.0e-5 || voice < layoutVoiceCount)
				voicesToProcess = voice + 1;
			sumGainSquares += voiceGains[voice] * voiceGains[voice];
		}
		const double normalization = 1.0 /
			std::sqrt(std::max(sumGainSquares, 1.0e-12));
		for (int channel = 0; channel < channels; ++channel)
		{
			auto finiteInput = [&](InputIds input)
			{
				const float value = inputs[input].getPolyVoltage(channel);
				return std::isfinite(value) ? static_cast<double>(value) : 0.0;
			};

			const double basePitch = finiteInput(VOCT_INPUT) + tune;
			const double centerPitch = basePitch + commonPitchCents / 1200.0;
			const double centerFrequency = dsp::FREQ_C4 * tfdsp::Exp2Taylor5(
				static_cast<float>(std::clamp(centerPitch, -16.0, 16.0)));
			const double pwmVoltage = inputs[PULSE_WIDTH_INPUT].isConnected() ?
				finiteInput(PULSE_WIDTH_INPUT) : internalPwmVoltage;
			const double pulseWidth = std::clamp(pulseWidthBase +
				0.45 * pulseWidthAmount * pwmVoltage / 5.0,
				0.05, 0.95);
			const double spreadControl = std::clamp(params[SPREAD].getValue() +
				spreadCvAmount * finiteInput(SPREAD_INPUT) / 5.0, 0.0, 1.0);
			const double spreadCents = tfdsp::UnisonSpreadCents(spreadControl);
			const double width = std::clamp(params[WIDTH].getValue() +
				widthCvAmount * finiteInput(WIDTH_INPUT) / 5.0, 0.0, 1.0);
			std::array<double, tfdsp::MaximumStackedOscillatorVoices>
				leftPanGains{};
			std::array<double, tfdsp::MaximumStackedOscillatorVoices>
				rightPanGains{};
			for (int voice = 0; voice < voicesToProcess; ++voice)
			{
				const double pan = std::clamp(
					width * panPositions[voice], -1.0, 1.0);
				// This equal-power law includes the sqrt(2) compensation that
				// makes a centred stereo side equal to MONO.
				leftPanGains[voice] = std::sqrt(1.0 - pan);
				rightPanGains[voice] = std::sqrt(1.0 + pan);
			}

			double monoMain = 0.0;
			double leftMain = 0.0;
			double rightMain = 0.0;
			double monoStackSub = 0.0;
			double leftStackSub = 0.0;
			double rightStackSub = 0.0;
			for (int voice = 0; voice < voicesToProcess; ++voice)
			{
				const double drift = individualDriftProcesses[channel][voice].Step(
					randomGenerator);
				const double trackingCents = MaximumTrackingErrorCentsPerOctave *
					trackingErrorDepth * trackingPositions[voice] * basePitch;
				double voicePitch = centerPitch + (spreadCents *
					pitchPositions[voice] + trackingCents) / 1200.0;
				if (centsDrift)
					voicePitch += MaximumIndividualDriftCents * individualDepth *
						drift / 1200.0;
				double frequency = dsp::FREQ_C4 * tfdsp::Exp2Taylor5(
					static_cast<float>(std::clamp(voicePitch, -16.0, 16.0)));
				if (!centsDrift)
					frequency += MaximumIndividualDriftHz * individualDepth * drift;
				frequency = std::clamp(frequency, 0.0, 0.45 * sampleRate);
				const double increment = frequency / sampleRate;
				const double gain = voiceGains[voice];
				const double leftGain = leftPanGains[voice];
				const double rightGain = rightPanGains[voice];

				if (needMain)
				{
					const double signal = voices[channel][voice].StepMain(
						increment, pulseWidth, waveformMix);
					if (needMono)
						monoMain += gain * signal;
					if (needLeft)
						leftMain += gain * leftGain * signal;
					if (needRight)
						rightMain += gain * rightGain * signal;
				}
				if (renderStackSub)
				{
					const double sub = voices[channel][voice].StepSub(increment);
					monoStackSub += gain * sub;
					leftStackSub += gain * leftGain * sub;
					rightStackSub += gain * rightGain * sub;
				}
			}

			monoMain *= normalization;
			leftMain *= normalization;
			rightMain *= normalization;
			monoStackSub *= normalization;
			leftStackSub *= normalization;
			rightStackSub *= normalization;
			const double centerSub = renderCenterSub ? centerSubs[channel].Step(
				0.5 * std::clamp(centerFrequency / sampleRate, 0.0, 0.45), 0.5) : 0.0;
			const double monoSub = centerSub + subModeMix *
				(monoStackSub - centerSub);
			const double leftSub = centerSub + subModeMix *
				(leftStackSub - centerSub);
			const double rightSub = centerSub + subModeMix *
				(rightStackSub - centerSub);

			auto writeOutput = [&](OutputIds output, double value)
			{
				outputs[output].setVoltage(static_cast<float>(
					tfdsp::RackOutputAdapter::ProcessPostDecimation(value)), channel);
			};
			writeOutput(MONO_OUTPUT, 5.0 * (monoMain + subLevel * monoSub));
			writeOutput(LEFT_OUTPUT, 5.0 * (leftMain + subLevel * leftSub));
			writeOutput(RIGHT_OUTPUT, 5.0 * (rightMain + subLevel * rightSub));
		}
	}

	void onReset(const ResetEvent& event) override
	{
		Module::onReset(event);
		ResetDsp();
		ConfigureDrift(DriftTimeSeconds(params[DRIFT_SPEED].getValue()));
	}

	void onSampleRateChange(const SampleRateChangeEvent& event) override
	{
		SetSampleRate(event.sampleRate);
	}
};

struct TfUnisonOscillatorWidget : ModuleWidget
{
	TfUnisonOscillatorWidget(TfUnisonOscillator* module)
	{
		setModule(module);
		setPanel(APP->window->loadSvg(asset::plugin(
			pluginInstance, "res/TfUnisonOscillator.svg")));

		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(
			box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(
			box.size.x - 2 * RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		addParam(createParam<TfSnapKnob>(Vec(13.83, 50), module,
			TfUnisonOscillator::OCTAVE));
		addParam(createParam<CKSS>(Vec(52, 53.85), module,
			TfUnisonOscillator::WAVEFORM));
		addParam(createParam<TfTrimpot>(Vec(81.07, 55.24), module,
			TfUnisonOscillator::TUNE));
		addParam(createParam<CKSS>(Vec(114, 53.85), module,
			TfUnisonOscillator::SUB_MODE));
		addParam(createParam<TfSnapKnob>(Vec(137.83, 50), module,
			TfUnisonOscillator::VOICES));

		addParam(createParam<TfCvKnob>(Vec(16, 110), module,
			TfUnisonOscillator::PULSE_WIDTH));
		addParam(createParam<TfCvKnob>(Vec(76, 110), module,
			TfUnisonOscillator::PULSE_WIDTH_CV_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(136, 110), module,
			TfUnisonOscillator::PWM_RATE));
		addParam(createParam<TfCvKnob>(Vec(16, 154), module,
			TfUnisonOscillator::SPREAD));
		addParam(createParam<TfCvKnob>(Vec(76, 154), module,
			TfUnisonOscillator::WIDTH));
		addParam(createParam<TfCvKnob>(Vec(136, 154), module,
			TfUnisonOscillator::SUB_LEVEL));

		addParam(createParam<TfTrimpot>(Vec(21.07, 203), module,
			TfUnisonOscillator::DRIFT_SPEED));
		addParam(createParam<TfTrimpot>(Vec(81.07, 203), module,
			TfUnisonOscillator::HUM));
		addParam(createParam<TfTrimpot>(Vec(141.07, 203), module,
			TfUnisonOscillator::COMMON_DRIFT));
		addParam(createParam<TfTrimpot>(Vec(21.07, 234), module,
			TfUnisonOscillator::INDIVIDUAL_DRIFT));
		addParam(createParam<CKSS>(Vec(83, 232.6), module,
			TfUnisonOscillator::INDIVIDUAL_DRIFT_MODE));
		addParam(createParam<TfTrimpot>(Vec(141.07, 234), module,
			TfUnisonOscillator::TRACKING));

		addParam(createParam<TfTrimpot>(Vec(21.07, 267), module,
			TfUnisonOscillator::SPREAD_CV_AMOUNT));
		addParam(createParam<TfTrimpot>(Vec(141.07, 267), module,
			TfUnisonOscillator::WIDTH_CV_AMOUNT));

		addInput(createInput<PJ301MPort>(Vec(18.15, 297), module,
			TfUnisonOscillator::VOCT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(58.15, 297), module,
			TfUnisonOscillator::PULSE_WIDTH_INPUT));
		addInput(createInput<PJ301MPort>(Vec(98.15, 297), module,
			TfUnisonOscillator::SPREAD_INPUT));
		addInput(createInput<PJ301MPort>(Vec(138.15, 297), module,
			TfUnisonOscillator::WIDTH_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(24, 333), module,
			TfUnisonOscillator::MONO_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(78, 333), module,
			TfUnisonOscillator::LEFT_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(132, 333), module,
			TfUnisonOscillator::RIGHT_OUTPUT));
	}
};

Model* modelTfUnisonOscillator = createModel<TfUnisonOscillator,
	TfUnisonOscillatorWidget>("TfUnisonOscillator");
