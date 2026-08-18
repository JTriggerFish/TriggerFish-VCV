#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <random>

#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/approx.hpp"
#include "tfdsp/control.hpp"
#include "tfdsp/rail.hpp"
#include "tfdsp/sampleRate.hpp"
#include "tfdsp/unison.hpp"
#include "tfdsp/wavefolder.hpp"

namespace
{
	constexpr double MinimumAliveTimeSeconds = 0.5;
	constexpr double MaximumAliveTimeSeconds = 120.0;

	double AliveTimeSeconds(double speed)
	{
		const double shapedSpeed = std::pow(std::clamp(speed, 0.0, 1.0), 2.0);
		return MaximumAliveTimeSeconds * std::pow(
			MinimumAliveTimeSeconds / MaximumAliveTimeSeconds, shapedSpeed);
	}

	struct AliveSpeedQuantity : engine::ParamQuantity
	{
		float getDisplayValue() override
		{
			return static_cast<float>(AliveTimeSeconds(getValue()));
		}

		void setDisplayValue(float displayValue) override
		{
			const double seconds = std::clamp(static_cast<double>(displayValue),
				MinimumAliveTimeSeconds, MaximumAliveTimeSeconds);
			const double exponent = std::log(seconds / MaximumAliveTimeSeconds) /
				std::log(MinimumAliveTimeSeconds / MaximumAliveTimeSeconds);
			setValue(static_cast<float>(std::sqrt(std::clamp(exponent, 0.0, 1.0))));
		}
	};

	struct UnisonSpreadQuantity : engine::ParamQuantity
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
}

struct TfWavefoldOscillator : Module
{
	enum ParamIds
	{
		OCTAVE,
		TUNE,
		MORPH,
		FOLD,
		SYMMETRY,
		FM_AMOUNT,
		MORPH_AMOUNT,
		FOLD_AMOUNT,
		SYMMETRY_AMOUNT,
		FM_MODE,
		CHARACTER,
		ALIVE_SPEED,
		MORPH_ALIVE,
		FOLD_ALIVE,
		SYMMETRY_ALIVE,
		UNISON_VOICES,
		UNISON_SPREAD,
		NUM_PARAMS
	};
	enum InputIds
	{
		VOCT_INPUT,
		FM_INPUT,
		MORPH_INPUT,
		FOLD_INPUT,
		SYMMETRY_INPUT,
		AUDIO_INPUT,
		NUM_INPUTS
	};
	enum OutputIds
	{
		OSCILLATOR_OUTPUT,
		FOLDED_OUTPUT,
		NUM_OUTPUTS
	};
	enum LightIds
	{
		NUM_LIGHTS
	};

	using OscillatorX2 = tfdsp::WavefoldOscillator<
		tfdsp::X2Resampler_Order7>;
	using OscillatorX4 = tfdsp::WavefoldOscillator<
		tfdsp::X4Resampler_Order7>;
	std::array<std::array<std::unique_ptr<OscillatorX2>,
		tfdsp::MaximumUnisonVoices>, PORT_MAX_CHANNELS> oscillatorsX2;
	std::array<std::array<std::unique_ptr<OscillatorX4>,
		tfdsp::MaximumUnisonVoices>, PORT_MAX_CHANNELS> oscillatorsX4;
	static constexpr int AliveProcessCount = 3;
	std::array<std::array<std::array<tfdsp::SmoothOrnsteinUhlenbeck,
		AliveProcessCount>, tfdsp::MaximumUnisonVoices>,
		PORT_MAX_CHANNELS> aliveProcesses{};
	std::random_device aliveSeed{};
	std::minstd_rand aliveRng;
	double configuredAliveTimeSeconds{};
	// 4x is the default quality mode; 2x is available as a lower-CPU fallback.
	int oversampling = 1;
	int activeOversampling = 1;
	double sampleRate = 48000.0;

	TfWavefoldOscillator() : aliveRng(aliveSeed())
	{
		config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
		configParam(OCTAVE, -3.0f, 3.0f, 0.0f, "Octave", " oct");
		getParamQuantity(OCTAVE)->snapEnabled = true;
		configParam(TUNE, -7.0f / 12.0f, 7.0f / 12.0f, 0.0f,
			"Tune", " semitones", 0.0f, 12.0f);
		configParam(MORPH, 0.0f, 1.0f, 0.5f, "Sine / triangle morph", "%",
			0.0f, 100.0f);
		configParam(FOLD, 0.0f, 1.0f, 0.4f, "Fold", "%", 0.0f, 100.0f);
		configParam(SYMMETRY, -1.0f, 1.0f, 0.0f, "Symmetry", "%",
			0.0f, 100.0f);
		configParam(FM_AMOUNT, -1.0f, 1.0f, 0.0f, "FM amount", "%",
			0.0f, 100.0f);
		configParam(MORPH_AMOUNT, -1.0f, 1.0f, 0.0f, "Wave CV amount", "%",
			0.0f, 100.0f);
		configParam(FOLD_AMOUNT, -1.0f, 1.0f, 0.0f, "Fold CV amount", "%",
			0.0f, 100.0f);
		configParam(SYMMETRY_AMOUNT, -1.0f, 1.0f, 0.0f,
			"Symmetry CV amount", "%", 0.0f, 100.0f);
		configSwitch(FM_MODE, 0.0f, 1.0f, 0.0f, "FM response",
			{"Exponential", "Linear (through-zero)"});
		configSwitch(CHARACTER, 0.0f, 2.0f, 2.0f, "Folder character",
			{"Serge", "Hinge", "Lockhart"});
		configParam<AliveSpeedQuantity>(ALIVE_SPEED, 0.0f, 1.0f, 0.5f,
			"Alive drift time", " s");
		configParam(MORPH_ALIVE, 0.0f, 1.0f, 0.5f,
			"Wave liveliness", "%", 0.0f, 100.0f);
		configParam(FOLD_ALIVE, 0.0f, 1.0f, 0.5f,
			"Fold liveliness", "%", 0.0f, 100.0f);
		configParam(SYMMETRY_ALIVE, 0.0f, 1.0f, 0.5f,
			"Symmetry liveliness", "%", 0.0f, 100.0f);
		configParam(UNISON_VOICES, 1.0f, 4.0f, 1.0f,
			"Unison voices");
		getParamQuantity(UNISON_VOICES)->snapEnabled = true;
		configParam<UnisonSpreadQuantity>(UNISON_SPREAD, 0.0f, 1.0f,
			static_cast<float>(tfdsp::UnisonSpreadControlForCents(4.0)),
			"Unison spread", " cents");

		configInput(VOCT_INPUT, "Pitch (1V/octave)");
		configInput(FM_INPUT,
			"Frequency modulation (exponential or through-zero linear)");
		configInput(MORPH_INPUT, "Wave CV (sine / triangle morph)");
		configInput(FOLD_INPUT, "Fold CV");
		configInput(SYMMETRY_INPUT, "Symmetry CV");
		configInput(AUDIO_INPUT,
			"External audio to folder (replaces internal oscillator at folder input)");
		configOutput(OSCILLATOR_OUTPUT,
			"Internal oscillator before folding");
		configOutput(FOLDED_OUTPUT,
			"Folder output (internal oscillator or external audio input)");

		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			for (int voice = 0; voice < tfdsp::MaximumUnisonVoices; ++voice)
			{
				oscillatorsX2[channel][voice] = std::make_unique<OscillatorX2>(
					tfdsp::CreateX2Resampler_Chebychev7);
				oscillatorsX4[channel][voice] = std::make_unique<OscillatorX4>(
					tfdsp::CreateX4Resampler_Cheby7);
			}
		}
		SetSampleRate(APP->engine->getSampleRate());
	}

	void ConfigureAlive(double timeSeconds)
	{
		configuredAliveTimeSeconds = timeSeconds;
		for (auto& channelProcesses : aliveProcesses)
			for (auto& voiceProcesses : channelProcesses)
				for (auto& process : voiceProcesses)
					process.ConfigureStationary(
						sampleRate, timeSeconds, 1.0, 100.0);
	}

	void ResetAlive()
	{
		for (auto& channelProcesses : aliveProcesses)
			for (auto& voiceProcesses : channelProcesses)
				for (auto& process : voiceProcesses)
					process.Reset();
	}

	void SetSampleRate(double nextSampleRate)
	{
		sampleRate = std::max(nextSampleRate, 1.0);
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			for (int voice = 0; voice < tfdsp::MaximumUnisonVoices; ++voice)
			{
				oscillatorsX2[channel][voice]->SetSampleRate(sampleRate);
				oscillatorsX4[channel][voice]->SetSampleRate(sampleRate);
			}
		}
		ConfigureAlive(AliveTimeSeconds(params[ALIVE_SPEED].getValue()));
	}

	void ResetDsp()
	{
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			for (int voice = 0; voice < tfdsp::MaximumUnisonVoices; ++voice)
			{
				oscillatorsX2[channel][voice]->Reset();
				oscillatorsX4[channel][voice]->Reset();
			}
		}
	}

	void process(const ProcessArgs& args) override
	{
		oversampling = std::clamp(oversampling, 0, 1);
		if (activeOversampling != oversampling)
		{
			activeOversampling = oversampling;
			ResetDsp();
		}

		const int channels = std::clamp(std::max({
			inputs[VOCT_INPUT].getChannels(), inputs[FM_INPUT].getChannels(),
			inputs[MORPH_INPUT].getChannels(), inputs[FOLD_INPUT].getChannels(),
			inputs[SYMMETRY_INPUT].getChannels(),
			inputs[AUDIO_INPUT].getChannels(), 1}), 1, PORT_MAX_CHANNELS);
		outputs[OSCILLATOR_OUTPUT].setChannels(channels);
		outputs[FOLDED_OUTPUT].setChannels(channels);

		const double octave = std::round(params[OCTAVE].getValue());
		const double tuningOffset = octave + params[TUNE].getValue();
		const double morphKnob = params[MORPH].getValue();
		const double foldKnob = params[FOLD].getValue();
		const double symmetryKnob = params[SYMMETRY].getValue();
		const double fmAmount = params[FM_AMOUNT].getValue();
		const double morphAmount = params[MORPH_AMOUNT].getValue();
		const double foldAmount = params[FOLD_AMOUNT].getValue();
		const double symmetryAmount = params[SYMMETRY_AMOUNT].getValue();
		const double morphAlive = params[MORPH_ALIVE].getValue();
		const double foldAlive = params[FOLD_ALIVE].getValue();
		const double symmetryAlive = params[SYMMETRY_ALIVE].getValue();
		const int requestedUnisonVoices = std::clamp(static_cast<int>(std::round(
			params[UNISON_VOICES].getValue())), 1, tfdsp::MaximumUnisonVoices);
		const double spreadCents = tfdsp::UnisonSpreadCents(
			params[UNISON_SPREAD].getValue());
		const double aliveTimeSeconds = AliveTimeSeconds(
			params[ALIVE_SPEED].getValue());
		if (std::abs(aliveTimeSeconds - configuredAliveTimeSeconds) > 1.0e-9)
			ConfigureAlive(aliveTimeSeconds);
		const bool linearFm = params[FM_MODE].getValue() > 0.5f;
		// CKSSThree numbers positions from bottom (0) to top (2).
		constexpr std::array<tfdsp::WavefolderCharacter, 3> charactersByPosition{{
			tfdsp::WavefolderCharacter::Serge,
			tfdsp::WavefolderCharacter::Hinge,
			tfdsp::WavefolderCharacter::Lockhart,
		}};
		const int characterPosition = std::clamp(static_cast<int>(std::round(
			params[CHARACTER].getValue())), 0, 2);
		const auto character = charactersByPosition[characterPosition];
		const bool externalInputConnected = inputs[AUDIO_INPUT].isConnected();

		for (int channel = 0; channel < channels; ++channel)
		{
			auto finiteInput = [&](InputIds input)
			{
				const float value = inputs[input].getPolyVoltage(channel);
				return std::isfinite(value) ? static_cast<double>(value) : 0.0;
			};

			double pitch = finiteInput(VOCT_INPUT) + tuningOffset;
			const double fmVoltage = fmAmount * finiteInput(FM_INPUT);
			double frequency;
			if (linearFm)
			{
				frequency = dsp::FREQ_C4 * tfdsp::Exp2Taylor5(
					static_cast<float>(std::clamp(pitch, -16.0, 16.0)));
				frequency += 200.0 * fmVoltage;
			}
			else
			{
				pitch += 0.2 * fmVoltage;
				frequency = dsp::FREQ_C4 * tfdsp::Exp2Taylor5(
					static_cast<float>(std::clamp(pitch, -16.0, 16.0)));
			}
			frequency = std::clamp(frequency,
				-0.45 * sampleRate, 0.45 * sampleRate);
			const double morphBase = morphKnob + morphAmount *
				finiteInput(MORPH_INPUT) / 5.0;
			const double foldBase = foldKnob + foldAmount *
				finiteInput(FOLD_INPUT) / 5.0;
			const double symmetryBase = symmetryKnob + symmetryAmount *
				finiteInput(SYMMETRY_INPUT) / 5.0;
			const double externalInput = finiteInput(AUDIO_INPUT) / 5.0;
			const auto pitchPositions =
				tfdsp::UnisonPitchPositions(requestedUnisonVoices);
			tfdsp::WavefoldOscillatorOutput rendered{};
			for (int voice = 0; voice < requestedUnisonVoices; ++voice)
			{
				const double morph = tfdsp::ApplyBoundedDrift(morphBase,
					aliveProcesses[channel][voice][0].Step(aliveRng), morphAlive);
				const double fold = tfdsp::ApplyBoundedDrift(foldBase,
					aliveProcesses[channel][voice][1].Step(aliveRng), foldAlive);
				const double symmetry = tfdsp::ApplyBoundedDrift(symmetryBase,
					aliveProcesses[channel][voice][2].Step(aliveRng), symmetryAlive,
					-1.0, 1.0);
				const double voiceFrequency = frequency * std::exp2(
					spreadCents * pitchPositions[voice] / 1200.0);
				const bool foldExternalInput = externalInputConnected && voice == 0;
				tfdsp::WavefoldOscillatorOutput voiceOutput;
				if (activeOversampling == 0)
				{
					oscillatorsX2[channel][voice]->SetCharacter(character);
					oscillatorsX2[channel][voice]->SetFolderAntialiasing(false);
					voiceOutput = oscillatorsX2[channel][voice]->StepWithInput(
						voiceFrequency, morph, fold, symmetry, externalInput,
						foldExternalInput);
				}
				else
				{
					oscillatorsX4[channel][voice]->SetCharacter(character);
					oscillatorsX4[channel][voice]->SetFolderAntialiasing(false);
					voiceOutput = oscillatorsX4[channel][voice]->StepWithInput(
						voiceFrequency, morph, fold, symmetry, externalInput,
						foldExternalInput);
				}
				rendered.oscillator += voiceOutput.oscillator;
				if (!externalInputConnected || voice == 0)
					rendered.folded += voiceOutput.folded;
			}
			const double unisonGain =
				tfdsp::UnisonOutputGain(requestedUnisonVoices);
			rendered.oscillator *= unisonGain;
			if (!externalInputConnected)
				rendered.folded *= unisonGain;
			outputs[OSCILLATOR_OUTPUT].setVoltage(static_cast<float>(
				tfdsp::RackOutputAdapter::ProcessPostDecimation(
					5.0 * rendered.oscillator)), channel);
			outputs[FOLDED_OUTPUT].setVoltage(static_cast<float>(
				tfdsp::RackOutputAdapter::ProcessPostDecimation(
					5.0 * rendered.folded)), channel);
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
		ResetAlive();
		ConfigureAlive(AliveTimeSeconds(params[ALIVE_SPEED].getValue()));
	}

	void onSampleRateChange(const SampleRateChangeEvent& event) override
	{
		SetSampleRate(event.sampleRate);
	}
};

struct TfWavefoldOscillatorWidget : ModuleWidget
{
	TfWavefoldOscillatorWidget(TfWavefoldOscillator* module)
	{
		setModule(module);
		setPanel(APP->window->loadSvg(asset::plugin(
			pluginInstance, "res/TfWavefoldOscillator.svg")));

		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(
			box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(
			box.size.x - 2 * RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		addParam(createParam<TfSnapKnob>(Vec(9.83, 48), module,
			TfWavefoldOscillator::OCTAVE));
		addParam(createParam<TfTrimpot>(Vec(99.07, 53.244), module,
			TfWavefoldOscillator::TUNE));
		addParam(createParam<CKSS>(Vec(59, 51.853), module,
			TfWavefoldOscillator::FM_MODE));
		addParam(createParam<TfSnapKnob>(Vec(135.83, 48), module,
			TfWavefoldOscillator::UNISON_VOICES));

		addParam(createParam<TfAudioKob>(Vec(12, 104), module,
			TfWavefoldOscillator::MORPH));
		addParam(createParam<TfAudioKob>(Vec(72, 104), module,
			TfWavefoldOscillator::FOLD));
		addParam(createParam<TfAudioKob>(Vec(132, 104), module,
			TfWavefoldOscillator::SYMMETRY));

		addParam(createParam<TfCvKnob>(Vec(7, 257), module,
			TfWavefoldOscillator::FM_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(53, 257), module,
			TfWavefoldOscillator::MORPH_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(99, 257), module,
			TfWavefoldOscillator::FOLD_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(145, 257), module,
			TfWavefoldOscillator::SYMMETRY_AMOUNT));

		addParam(createParam<CKSSThree>(Vec(37.5, 205), module,
			TfWavefoldOscillator::CHARACTER));
		addParam(createParam<TfCvKnob>(Vec(121, 205), module,
			TfWavefoldOscillator::UNISON_SPREAD));

		addParam(createParam<TfTrimpot>(Vec(15.07, 171), module,
			TfWavefoldOscillator::ALIVE_SPEED));
		addParam(createParam<TfTrimpot>(Vec(59.07, 171), module,
			TfWavefoldOscillator::MORPH_ALIVE));
		addParam(createParam<TfTrimpot>(Vec(103.07, 171), module,
			TfWavefoldOscillator::FOLD_ALIVE));
		addParam(createParam<TfTrimpot>(Vec(147.07, 171), module,
			TfWavefoldOscillator::SYMMETRY_ALIVE));

		addInput(createInput<PJ301MPort>(Vec(6, 299), module,
			TfWavefoldOscillator::VOCT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(42, 299), module,
			TfWavefoldOscillator::FM_INPUT));
		addInput(createInput<PJ301MPort>(Vec(78, 299), module,
			TfWavefoldOscillator::MORPH_INPUT));
		addInput(createInput<PJ301MPort>(Vec(114, 299), module,
			TfWavefoldOscillator::FOLD_INPUT));
		addInput(createInput<PJ301MPort>(Vec(150, 299), module,
			TfWavefoldOscillator::SYMMETRY_INPUT));
		addInput(createInput<PJ301MPort>(Vec(24, 334), module,
			TfWavefoldOscillator::AUDIO_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(78, 334), module,
			TfWavefoldOscillator::OSCILLATOR_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(132, 334), module,
			TfWavefoldOscillator::FOLDED_OUTPUT));
	}

	void appendContextMenu(Menu* menu) override
	{
		TfWavefoldOscillator* module =
			dynamic_cast<TfWavefoldOscillator*>(this->module);
		if (!module)
			return;
		menu->addChild(new MenuSeparator);
		menu->addChild(createIndexPtrSubmenuItem("Oversampling",
			{"2x (lower CPU)", "4x (default)"}, &module->oversampling));
	}
};

Model* modelTfWavefoldOscillator = createModel<TfWavefoldOscillator,
	TfWavefoldOscillatorWidget>("TfWavefoldOscillator");
