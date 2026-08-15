#include "models/DiodeLadderFilter.hpp"
#include "models/Tb303Voice.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

#include "plugin.hpp"
#include "components.hpp"
#include "tfdsp/approx.hpp"
#include "tfdsp/filters.hpp"
#include "tfdsp/sampleRate.hpp"

struct TfDiodeLadderFilter : Module
{
	enum ParamIds
	{
		CUTOFF,
		RESONANCE,
		HIGH_RESONANCE,
		DRIVE,
		BASS,
		CV_AMOUNT,
		FM_AMOUNT,
		RES_AMOUNT,
		ENV_AMOUNT,
		NORMAL_DECAY,
		ACCENT_DECAY,
		ACCENT_AMOUNT,
		VCA_DECAY,
		VCA_CV_AMOUNT,
		// Appended to preserve all parameter IDs used by earlier patches.
		ACCENT_SWEEP_MODE,
		NUM_PARAMS
	};
	enum InputIds
	{
		AUDIO_INPUT,
		VOCT_INPUT,
		CV_INPUT,
		FM_INPUT,
		RES_INPUT,
		GATE_INPUT,
		ACCENT_INPUT,
		VCA_CV_INPUT,
		NUM_INPUTS
	};
	enum OutputIds
	{
		LP_OUTPUT,
		VCA_OUTPUT,
		NUM_OUTPUTS
	};
	enum LightIds
	{
		NUM_LIGHTS
	};

	using FilterX2 = tfdsp::DiodeLadderFilter<tfdsp::X2Resampler_Order7>;
	using FilterX4 = tfdsp::DiodeLadderFilter<tfdsp::X4Resampler_Order7>;
	std::array<std::unique_ptr<FilterX2>, PORT_MAX_CHANNELS> filtersX2;
	std::array<std::unique_ptr<FilterX4>, PORT_MAX_CHANNELS> filtersX4;
	std::array<tfdsp::FirstOrderHighPassZdf<float>, PORT_MAX_CHANNELS> fmHighPass{};
	std::array<tfdsp::Tb303Articulation, PORT_MAX_CHANNELS> articulations{};
	std::array<tfdsp::Tb303Vca, PORT_MAX_CHANNELS> vcasX2{};
	std::array<tfdsp::Tb303Vca, PORT_MAX_CHANNELS> vcasX4{};

	// 0 = 2x, 1 = 4x. Four-times is the quality-first default; 2x roughly
	// doubles throughput and remains useful for large polyphonic patches.
	int oversampling = 1;
	int activeOversampling = 1;
	int articulationMode = 0;
	float normalizedFmHighPass{};

	TfDiodeLadderFilter()
	{
		config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
		constexpr float minimumCutoffHz = 10.0f;
		constexpr float defaultCutoffHz = 500.0f;
		constexpr float maximumCutoffHz = 20000.0f;
		constexpr float standardEnvelopePeakVolts = 10.0f;
		const float minimumCutoffPitch = std::log2(
			minimumCutoffHz / dsp::FREQ_C4);
		const float defaultCutoffPitch = std::log2(
			defaultCutoffHz / dsp::FREQ_C4);
		const float maximumCutoffPitch = std::log2(
			maximumCutoffHz / dsp::FREQ_C4);
		const float defaultCvAmount = std::log2(
			maximumCutoffHz / defaultCutoffHz) / standardEnvelopePeakVolts;

		configParam(CUTOFF, minimumCutoffPitch, maximumCutoffPitch,
			defaultCutoffPitch,
			"Cutoff", " Hz", 2.0f, dsp::FREQ_C4);
		configParam(RESONANCE, 0.0f, 1.0f, 0.0f, "Resonance", "%", 0.0f,
			100.0f);
		configSwitch(HIGH_RESONANCE, 0.0f, 1.0f, 0.0f, "Resonance range",
			{"Stock", "High"});
		configParam(DRIVE, -60.0f, 36.469f, 0.0f, "Drive", " dB");
		configParam(BASS, 0.0f, 1.0f, 0.0f, "Bass extension", "%", 0.0f,
			100.0f);
		configParam(CV_AMOUNT, -1.0f, 1.0f, defaultCvAmount,
			"Exponential cutoff CV amount", "%", 0.0f, 100.0f);
		configParam(FM_AMOUNT, -1.0f, 1.0f, 0.0f,
			"Linear audio FM amount", "%", 0.0f, 100.0f);
		configParam(RES_AMOUNT, -1.0f, 1.0f, 0.0f,
			"Resonance CV amount", "%", 0.0f, 100.0f);
		configParam(ENV_AMOUNT, 0.0f, 1.0f, 1.0f / 3.0f,
			"Internal filter envelope amount", "%", 0.0f, 100.0f);
		configParam(NORMAL_DECAY, std::log10(0.030f), std::log10(3.0f),
			std::log10(0.500f), "Normal filter envelope decay", " ms",
			10.0f, 1000.0f);
		configParam(ACCENT_DECAY, std::log10(0.030f), std::log10(3.0f),
			std::log10(0.200f), "Accented filter envelope decay", " ms",
			10.0f, 1000.0f);
		configParam(ACCENT_AMOUNT, 0.0f, 1.0f, 0.5f,
			"Accent sweep amount", "%", 0.0f, 100.0f);
		configParam(VCA_DECAY, 0.0f, 1.0f, 0.5f,
			"VCA decay / sustain", "%", 0.0f, 100.0f);
		configParam(VCA_CV_AMOUNT, 0.0f, 1.0f, 1.0f,
			"VCA CV amount", "%", 0.0f, 100.0f);
		configSwitch(ACCENT_SWEEP_MODE, 0.0f, 3.0f, 2.0f,
			"Accent sweep", {"Off", "Fast", "Normal", "Slow"});

		configInput(AUDIO_INPUT, "Audio");
		configInput(VOCT_INPUT, "Cutoff (1V/octave)");
		configInput(CV_INPUT, "Attenuverted exponential cutoff CV");
		configInput(FM_INPUT, "AC-coupled linear filter FM");
		configInput(RES_INPUT, "Resonance CV");
		configInput(GATE_INPUT, "Articulation gate");
		configInput(ACCENT_INPUT, "Accent CV");
		configInput(VCA_CV_INPUT, "VCA CV");
		configOutput(LP_OUTPUT, "Low-pass audio");
		configOutput(VCA_OUTPUT, "VCA audio");
		configBypass(AUDIO_INPUT, LP_OUTPUT);
		configBypass(AUDIO_INPUT, VCA_OUTPUT);

		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			filtersX2[channel] = std::make_unique<FilterX2>(
				tfdsp::CreateX2Resampler_Chebychev7);
			filtersX4[channel] = std::make_unique<FilterX4>(
				tfdsp::CreateX4Resampler_Cheby7);
		}
		SetSampleRate(APP->engine->getSampleRate());
	}

	void SetSampleRate(float sampleRate)
	{
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			filtersX2[channel]->SetSampleRate(sampleRate);
			filtersX4[channel]->SetSampleRate(sampleRate);
			fmHighPass[channel].Reset();
			articulations[channel].SetSampleRate(sampleRate);
			vcasX2[channel].SetSampleRate(2.0 * sampleRate);
			vcasX4[channel].SetSampleRate(4.0 * sampleRate);
		}
		normalizedFmHighPass = 5.0f / (0.5f * sampleRate);
	}

	void ResetDsp()
	{
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			filtersX2[channel]->Reset();
			filtersX4[channel]->Reset();
			fmHighPass[channel].Reset();
			articulations[channel].Reset();
			vcasX2[channel].Reset();
			vcasX4[channel].Reset();
		}
	}

	void ResetActiveFilter()
	{
		for (int channel = 0; channel < PORT_MAX_CHANNELS; ++channel)
		{
			if (activeOversampling == 0)
			{
				filtersX2[channel]->Reset();
				vcasX2[channel].Reset();
			}
			else
			{
				filtersX4[channel]->Reset();
				vcasX4[channel].Reset();
			}
		}
	}

	void process(const ProcessArgs& args) override
	{
		oversampling = std::clamp(oversampling, 0, 1);
		if (activeOversampling != oversampling)
		{
			activeOversampling = oversampling;
			// Changing quality invalidates the selected resamplers and the VCA's
			// rate-dependent C38 coupling state. Articulation remains continuous.
			ResetActiveFilter();
		}

		const int channels = std::max(inputs[AUDIO_INPUT].getChannels(), 1);
		outputs[LP_OUTPUT].setChannels(channels);
		outputs[VCA_OUTPUT].setChannels(channels);
		const float cutoffKnob = params[CUTOFF].getValue();
		const float cvAmount = params[CV_AMOUNT].getValue();
		const float fmAmount = params[FM_AMOUNT].getValue();
		const float resonanceKnob = params[RESONANCE].getValue();
		const float resonanceAmount = params[RES_AMOUNT].getValue();
		const float envelopeAmount = params[ENV_AMOUNT].getValue();
		const double normalDecay = std::pow(10.0,
			static_cast<double>(params[NORMAL_DECAY].getValue()));
		const double accentDecay = std::pow(10.0,
			static_cast<double>(params[ACCENT_DECAY].getValue()));
		const float accentAmount = params[ACCENT_AMOUNT].getValue();
		const float vcaDecay = params[VCA_DECAY].getValue();
		const float vcaCvAmount = params[VCA_CV_AMOUNT].getValue();
		const int accentSweepMode = std::clamp(static_cast<int>(std::lround(
			params[ACCENT_SWEEP_MODE].getValue())), 0, 3);
		const bool highResonance = params[HIGH_RESONANCE].getValue() > 0.5f;
		const float bass = params[BASS].getValue();
		const bool internalEnvelopePatched = inputs[GATE_INPUT].isConnected();
		const bool externalVcaPatched = inputs[VCA_CV_INPUT].isConnected();
		const float driveDb = params[DRIVE].getValue();
		const double driveGain = driveDb <= -59.99f ? 0.0 :
			std::pow(10.0, static_cast<double>(driveDb) / 20.0);

		for (int channel = 0; channel < channels; ++channel)
		{
			articulations[channel].SetMode(articulationMode == 0 ?
				tfdsp::Tb303Articulation::Mode::Stock :
				tfdsp::Tb303Articulation::Mode::DevilFish);
			articulations[channel].SetAccentSweepMode(
				static_cast<tfdsp::Tb303AccentSweep::Mode>(
					std::clamp(accentSweepMode, 0, 3)));
			const float audio = inputs[AUDIO_INPUT].getPolyVoltage(channel);
			const float voct = inputs[VOCT_INPUT].getPolyVoltage(channel);
			const float cv = inputs[CV_INPUT].getPolyVoltage(channel);
			const float fm = inputs[FM_INPUT].getPolyVoltage(channel);
			const float resonanceCv = inputs[RES_INPUT].getPolyVoltage(channel);
			const float accent = inputs[ACCENT_INPUT].getPolyVoltage(channel);
			const float gate = inputs[GATE_INPUT].getPolyVoltage(channel);

			const float finiteAudio = std::isfinite(audio) ? audio : 0.0f;
			const double resonance = std::clamp(
				static_cast<double>(resonanceKnob + resonanceAmount *
					(std::isfinite(resonanceCv) ? resonanceCv / 10.0f : 0.0f)),
				0.0, 1.0);
			const auto envelope = articulations[channel].Step(gate, accent,
				resonance, normalDecay, accentDecay, vcaDecay);
			// The Q9/R64/R65 bias makes Env Mod scale the envelope around an
			// approximately 31.37% pivot instead of simply adding a positive
			// sweep. Increasing ENV therefore opens the attack and closes the
			// tail, keeping more of the motion in the audible cutoff range.
			// The six-octave software range extends beyond ordinary stock use
			// without spending most of the knob travel against the cutoff clamp.
			constexpr double envelopePivot = 0.3137;
			const double envelopePitch = internalEnvelopePatched ?
				6.0 * envelopeAmount * (envelope.mainEnvelope - envelopePivot) +
					2.0 * accentAmount * envelope.filterAccent : 0.0;
			const float pitch = std::clamp(static_cast<float>(cutoffKnob +
				(std::isfinite(voct) ? voct : 0.0f) + cvAmount *
				(std::isfinite(cv) ? cv : 0.0f) + envelopePitch), -10.0f, 10.0f);
			double cutoff = dsp::FREQ_C4 * tfdsp::Exp2Taylor5(pitch);
			const float acFm = fmHighPass[channel](
				std::isfinite(fm) ? fm : 0.0f, normalizedFmHighPass);
			// Linear-Hz modulation. At full amount a nominal +/-5 V Rack
			// signal sweeps +/-1 kHz, keeping the negative half-cycle useful
			// at ordinary cutoff settings instead of spending it on the floor.
			cutoff += 200.0 * fmAmount * acFm;

			const double externalVca = inputs[VCA_CV_INPUT].getPolyVoltage(channel);
			const double baseVcaControl = vcaCvAmount * (externalVcaPatched ?
				(std::isfinite(externalVca) ? externalVca / 10.0 : 0.0) :
				envelope.volumeEnvelope);
			const double vcaAccentControl = accentAmount * envelope.vcaAccent;

			float output = 0.0f;
			float vcaOutput = 0.0f;
			if (activeOversampling == 0)
			{
				const auto rendered = filtersX2[channel]->StepWithPostProcessor(
					finiteAudio, cutoff, resonance, highResonance, driveGain,
					bass, baseVcaControl, [&](double audioValue, double control)
					{
						return vcasX2[channel].Step(audioValue, control,
							vcaAccentControl);
					});
				output = rendered.lowPass;
				vcaOutput = rendered.postProcessed;
			}
			else
			{
				const auto rendered = filtersX4[channel]->StepWithPostProcessor(
					finiteAudio, cutoff, resonance, highResonance, driveGain,
					bass, baseVcaControl, [&](double audioValue, double control)
					{
						return vcasX4[channel].Step(audioValue, control,
							vcaAccentControl);
					});
				output = rendered.lowPass;
				vcaOutput = rendered.postProcessed;
			}
			outputs[LP_OUTPUT].setVoltage(
				std::isfinite(output) ? std::clamp(output, -12.0f, 12.0f) : 0.0f,
				channel);
			outputs[VCA_OUTPUT].setVoltage(std::isfinite(vcaOutput) ?
				vcaOutput : 0.0f, channel);
		}
	}

	json_t* dataToJson() override
	{
		json_t* root = json_object();
		json_object_set_new(root, "oversampling", json_integer(oversampling));
		json_object_set_new(root, "articulationMode", json_integer(articulationMode));
		return root;
	}

	void dataFromJson(json_t* root) override
	{
		if (json_t* value = json_object_get(root, "oversampling"))
			oversampling = std::clamp(static_cast<int>(json_integer_value(value)), 0, 1);
		if (json_t* value = json_object_get(root, "articulationMode"))
			articulationMode = std::clamp(
				static_cast<int>(json_integer_value(value)), 0, 1);
		// Migrate development patches saved before Accent Sweep became a panel
		// parameter. New patches store it in Rack's ordinary parameter array.
		if (json_t* value = json_object_get(root, "accentSweepMode"))
			params[ACCENT_SWEEP_MODE].setValue(static_cast<float>(std::clamp(
				static_cast<int>(json_integer_value(value)), 0, 3)));
	}

	void onReset(const ResetEvent& event) override
	{
		Module::onReset(event);
		oversampling = 1;
		activeOversampling = 1;
		articulationMode = 0;
		ResetDsp();
	}

	void onSampleRateChange(const SampleRateChangeEvent& event) override
	{
		SetSampleRate(event.sampleRate);
	}
};

struct TfDiodeLadderFilterWidget : ModuleWidget
{
	TfDiodeLadderFilterWidget(TfDiodeLadderFilter* module)
	{
		setModule(module);
		setPanel(APP->window->loadSvg(asset::plugin(
			pluginInstance, "res/TfDiodeLadderFilter.svg")));

		auto* transistorGraphic = new SvgWidget;
		transistorGraphic->setSvg(APP->window->loadSvg(asset::plugin(
			pluginInstance, "res/TfDiodeConnectedTransistor.svg")));
		transistorGraphic->box.pos = Vec(1, 44);
		addChild(transistorGraphic);

		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
		addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
		addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH,
			RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

		addParam(createParam<TfLargeAudioKnob>(Vec(35, 43), module,
			TfDiodeLadderFilter::CUTOFF));
		addParam(createParam<TfLargeAudioKnob>(Vec(155, 43), module,
			TfDiodeLadderFilter::RESONANCE));
		addParam(createParam<CKSS>(Vec(113, 65), module,
			TfDiodeLadderFilter::HIGH_RESONANCE));

		addParam(createParam<TfAudioKob>(Vec(31, 123), module,
			TfDiodeLadderFilter::DRIVE));
		addParam(createParam<TfSnapKnob>(Vec(106, 123), module,
			TfDiodeLadderFilter::ACCENT_SWEEP_MODE));
		addParam(createParam<TfAudioKob>(Vec(171, 123), module,
			TfDiodeLadderFilter::BASS));

		addParam(createParam<TfCvKnob>(Vec(10, 181), module,
			TfDiodeLadderFilter::ENV_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(58, 181), module,
			TfDiodeLadderFilter::NORMAL_DECAY));
		addParam(createParam<TfCvKnob>(Vec(106, 181), module,
			TfDiodeLadderFilter::ACCENT_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(154, 181), module,
			TfDiodeLadderFilter::ACCENT_DECAY));
		addParam(createParam<TfCvKnob>(Vec(202, 181), module,
			TfDiodeLadderFilter::VCA_DECAY));

		addParam(createParam<TfCvKnob>(Vec(34, 228), module,
			TfDiodeLadderFilter::CV_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(82, 228), module,
			TfDiodeLadderFilter::FM_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(130, 228), module,
			TfDiodeLadderFilter::RES_AMOUNT));
		addParam(createParam<TfCvKnob>(Vec(178, 228), module,
			TfDiodeLadderFilter::VCA_CV_AMOUNT));

		addInput(createInput<PJ301MPort>(Vec(12, 282), module,
			TfDiodeLadderFilter::VOCT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(60, 282), module,
			TfDiodeLadderFilter::CV_INPUT));
		addInput(createInput<PJ301MPort>(Vec(108, 282), module,
			TfDiodeLadderFilter::GATE_INPUT));
		addInput(createInput<PJ301MPort>(Vec(156, 282), module,
			TfDiodeLadderFilter::ACCENT_INPUT));
		addInput(createInput<PJ301MPort>(Vec(204, 282), module,
			TfDiodeLadderFilter::VCA_CV_INPUT));
		addInput(createInput<PJ301MPort>(Vec(12, 334), module,
			TfDiodeLadderFilter::AUDIO_INPUT));
		addInput(createInput<PJ301MPort>(Vec(60, 334), module,
			TfDiodeLadderFilter::FM_INPUT));
		addInput(createInput<PJ301MPort>(Vec(108, 334), module,
			TfDiodeLadderFilter::RES_INPUT));
		addOutput(createOutput<PJ301MPort>(Vec(156, 334), module,
			TfDiodeLadderFilter::LP_OUTPUT));
		addOutput(createOutput<PJ301MPort>(Vec(204, 334), module,
			TfDiodeLadderFilter::VCA_OUTPUT));
	}

	void appendContextMenu(Menu* menu) override
	{
		TfDiodeLadderFilter* module = dynamic_cast<TfDiodeLadderFilter*>(this->module);
		if (!module)
			return;
		menu->addChild(new MenuSeparator);
		menu->addChild(createIndexPtrSubmenuItem("Oversampling",
			{"2x", "4x"}, &module->oversampling));
		menu->addChild(createIndexPtrSubmenuItem("Articulation",
			{"TB-303", "Devil Fish"}, &module->articulationMode));
	}
};

Model* modelTfDiodeLadderFilter = createModel<TfDiodeLadderFilter,
	TfDiodeLadderFilterWidget>("TfDiodeLadderFilter");
