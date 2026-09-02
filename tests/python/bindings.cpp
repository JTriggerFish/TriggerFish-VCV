#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "models/VCAcore.hpp"
#include "models/Arp4019Vca.hpp"
#include "models/Arp4072Filter.hpp"
#include "models/ArpEnvelope.hpp"
#include "models/DiodeLadderFilter.hpp"
#include "models/OtaVca.hpp"
#include "models/Tb303Voice.hpp"
#include "models/Tb303Oscillator.hpp"
#include "models/VdpOscillator.hpp"
#include "models/VdpSplitOscillator.hpp"
#include "tfdsp/control.hpp"
#include "tfdsp/late_reverb.hpp"
#include "tfdsp/noise.hpp"
#include "tfdsp/percussion/crash_cymbal.hpp"
#include "tfdsp/room_reverb.hpp"
#include "tfdsp/sampleRate.hpp"
#include "tfdsp/unison.hpp"
#include "tfdsp/unison_oscillator.hpp"
#include "tfdsp/wavefolder.hpp"

namespace py = pybind11;

namespace
{
	void RequireSameSize(const py::buffer_info& left, const py::buffer_info& right,
		const char* leftName, const char* rightName)
	{
		if (left.ndim != 1 || right.ndim != 1)
			throw std::invalid_argument("DSP inputs must be one-dimensional arrays");
		if (left.shape[0] != right.shape[0])
			throw std::invalid_argument(std::string(leftName) + " and " + rightName + " must have the same length");
	}

	py::array_t<float> RenderLateReverbWallImpulse(
		py::ssize_t sampleCount, double sampleRate, double space, double aspect,
		double decay, double damping, double diffusion, bool optimized)
	{
		if (sampleCount <= 0)
			throw std::invalid_argument("sample_count must be positive");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result({sampleCount,
			static_cast<py::ssize_t>(tfdsp::LateReverb::WallCount),
			static_cast<py::ssize_t>(tfdsp::LateReverb::WallCount)});
		auto output = result.mutable_unchecked<3>();
		tfdsp::RoomReverbControls roomControls;
		roomControls.space = static_cast<float>(space);
		roomControls.aspect = static_cast<float>(aspect);
		const auto room = tfdsp::RoomReverb::MakeRoom(roomControls);

		for (std::size_t inputWall = 0;
			 inputWall < tfdsp::LateReverb::WallCount; ++inputWall)
		{
			tfdsp::LateReverb reverb;
			reverb.SetFlavour(optimized ? tfdsp::LateReverbFlavour::Optimized
			                              : tfdsp::LateReverbFlavour::Base);
			reverb.SetSampleRate(sampleRate);
			tfdsp::LateReverbControls controls;
			controls.decay = static_cast<float>(decay);
			controls.damping = static_cast<float>(damping);
			controls.diffusion = static_cast<float>(diffusion);
			controls.modulation = 0.f;
			controls.shimmer = 0.f;
			for (std::size_t axis = 0; axis < 3; ++axis)
				controls.roomDimensionsMetres[axis] =
					static_cast<float>(room.dimensionsMetres[axis]);
			for (py::ssize_t sample = 0; sample < sampleCount; ++sample)
			{
				tfdsp::LateReverb::WallFrame input{};
				if (sample == 0)
					input[inputWall] = 1.f;
				const auto frame = reverb.ProcessWallFrame(input, controls);
				for (std::size_t outputWall = 0;
					 outputWall < tfdsp::LateReverb::WallCount; ++outputWall)
					output(sample, outputWall, inputWall) = frame[outputWall];
			}
		}
		return result;
	}

	template<typename Filter>
	py::array_t<float> RenderDiodeLadderControls(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> cutoff,
		py::array_t<double, py::array::c_style | py::array::forcecast> linearFm,
		py::array_t<double, py::array::c_style | py::array::forcecast> resonance,
		bool highResonance, double driveGain, double bass, double sampleRate)
	{
		const auto audioInfo = audio.request();
		const auto cutoffInfo = cutoff.request();
		const auto linearFmInfo = linearFm.request();
		const auto resonanceInfo = resonance.request();
		RequireSameSize(audioInfo, cutoffInfo, "audio", "cutoff");
		RequireSameSize(audioInfo, linearFmInfo, "audio", "linear_fm");
		RequireSameSize(audioInfo, resonanceInfo, "audio", "resonance");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto cutoffValues = cutoff.unchecked<1>();
		auto linearFmValues = linearFm.unchecked<1>();
		auto resonanceValues = resonance.unchecked<1>();
		Filter model([]
		{
			if constexpr (Filter::OversamplingFactor == 1)
				return tfdsp::CreateDummyResampler();
			else if constexpr (Filter::OversamplingFactor == 2)
				return tfdsp::CreateX2Resampler_Chebychev7();
			else
				return tfdsp::CreateX4Resampler_Cheby7();
		});
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
		{
			output(i) = model.StepLogCutoffModulated(audioValues(i),
				std::log2(std::max(cutoffValues(i),
					std::numeric_limits<double>::min())), linearFmValues(i),
				resonanceValues(i), highResonance, driveGain, bass);
		}
		return result;
	}

	template<typename Filter>
	py::array_t<double> RenderDiodeLadderDiagnostics(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> cutoff,
		py::array_t<double, py::array::c_style | py::array::forcecast> linearFm,
		py::array_t<double, py::array::c_style | py::array::forcecast> resonance,
		bool highResonance, double driveGain, double bass, double sampleRate)
	{
		const auto audioInfo = audio.request();
		const auto cutoffInfo = cutoff.request();
		const auto linearFmInfo = linearFm.request();
		const auto resonanceInfo = resonance.request();
		RequireSameSize(audioInfo, cutoffInfo, "audio", "cutoff");
		RequireSameSize(audioInfo, linearFmInfo, "audio", "linear_fm");
		RequireSameSize(audioInfo, resonanceInfo, "audio", "resonance");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<double> result({audioInfo.shape[0], py::ssize_t{3}});
		auto output = result.mutable_unchecked<2>();
		auto audioValues = audio.unchecked<1>();
		auto cutoffValues = cutoff.unchecked<1>();
		auto linearFmValues = linearFm.unchecked<1>();
		auto resonanceValues = resonance.unchecked<1>();
		Filter model([]
		{
			if constexpr (Filter::OversamplingFactor == 1)
				return tfdsp::CreateDummyResampler();
			else if constexpr (Filter::OversamplingFactor == 2)
				return tfdsp::CreateX2Resampler_Chebychev7();
			else
				return tfdsp::CreateX4Resampler_Cheby7();
		});
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
		{
			output(i, 0) = model.StepLogCutoffModulated(audioValues(i),
				std::log2(std::max(cutoffValues(i),
					std::numeric_limits<double>::min())), linearFmValues(i),
				resonanceValues(i), highResonance, driveGain, bass);
			output(i, 1) = static_cast<double>(model.LastIterations());
			output(i, 2) = static_cast<double>(model.SolverFailures());
		}
		return result;
	}

	enum class DetuneMethod
	{
		LegacyDouble,
		Optimized,
		ReferenceDouble
	};

	template<DetuneMethod Method>
	double DetuneLinear(double pitch, double detune, double referenceFrequency)
	{
		if (!std::isfinite(pitch) || !std::isfinite(detune) ||
			!std::isfinite(referenceFrequency) || referenceFrequency <= 0.0)
			return 0.0;
		if constexpr (Method == DetuneMethod::Optimized)
			return tfdsp::detune::linear(pitch, detune, referenceFrequency);

		const double boundedPitch = std::clamp(pitch, -100.0, 100.0);
		const double pitchRatio = std::exp2(boundedPitch);
		const double value = pitchRatio + detune / referenceFrequency;
		if (value <= 1.0e-8)
			return std::log2(1.0e-8);
		if constexpr (Method == DetuneMethod::LegacyDouble)
		{
			static const double ln2 = std::log(2.0);
			return std::log(value) / ln2;
		}

		const double relativeDetune = detune / (referenceFrequency * pitchRatio);
		return boundedPitch + std::log1p(relativeDetune) / std::log(2.0);
	}

	template<DetuneMethod Method>
	py::array_t<double> RenderDetune(
		py::array_t<double, py::array::c_style | py::array::forcecast> pitch,
		py::array_t<double, py::array::c_style | py::array::forcecast> detune,
		double referenceFrequency)
	{
		const auto pitchInfo = pitch.request();
		const auto detuneInfo = detune.request();
		RequireSameSize(pitchInfo, detuneInfo, "pitch", "detune");
		if (!(referenceFrequency > 0.0))
			throw std::invalid_argument("reference_frequency must be positive");

		py::array_t<double> result(pitchInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto pitchValues = pitch.unchecked<1>();
		auto detuneValues = detune.unchecked<1>();
		for (py::ssize_t i = 0; i < pitchInfo.shape[0]; ++i)
			output(i) = DetuneLinear<Method>(pitchValues(i), detuneValues(i), referenceFrequency);
		return result;
	}

	template<typename Oscillator>
	py::array_t<float> RenderVdpo(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> damping,
		py::array_t<double, py::array::c_style | py::array::forcecast> angularFrequency,
		double sampleRate)
	{
		const auto audioInfo = audio.request();
		const auto dampingInfo = damping.request();
		const auto frequencyInfo = angularFrequency.request();
		RequireSameSize(audioInfo, dampingInfo, "audio", "damping");
		RequireSameSize(audioInfo, frequencyInfo, "audio", "angular_frequency");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto dampingValues = damping.unchecked<1>();
		auto frequencyValues = angularFrequency.unchecked<1>();
		Oscillator model(tfdsp::CreateX4Resampler_Cheby7);
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
			output(i) = model.Step(audioValues(i), dampingValues(i), frequencyValues(i));
		return result;
	}

	template<typename Oscillator, bool UseFastExp2>
	py::array_t<float> RenderVdpoPitch(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> damping,
		py::array_t<double, py::array::c_style | py::array::forcecast> pitch,
		double sampleRate)
	{
		const auto audioInfo = audio.request();
		const auto dampingInfo = damping.request();
		const auto pitchInfo = pitch.request();
		RequireSameSize(audioInfo, dampingInfo, "audio", "damping");
		RequireSameSize(audioInfo, pitchInfo, "audio", "pitch");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto dampingValues = damping.unchecked<1>();
		auto pitchValues = pitch.unchecked<1>();
		Oscillator model(tfdsp::CreateX4Resampler_Cheby7);
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
		{
			const float pitchValue = static_cast<float>(pitchValues(i));
			const float ratio = UseFastExp2 ? tfdsp::Exp2Taylor5(pitchValue) : std::exp2(pitchValue);
			const float frequency = 261.625565f * ratio;
			output(i) = model.Step(audioValues(i), dampingValues(i), 2.0 * tfdsp::PI * frequency);
		}
		return result;
	}

	template<typename Filter, bool Order9 = false>
	py::array_t<float> RenderDiodeLadder(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		double cutoff, double resonance, bool highResonance, double driveGain,
		double bass, double sampleRate)
	{
		const auto audioInfo = audio.request();
		if (audioInfo.ndim != 1)
			throw std::invalid_argument("audio must be a one-dimensional array");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		Filter model([]
		{
			if constexpr (Filter::OversamplingFactor == 1)
			{
				return tfdsp::CreateDummyResampler();
			}
			else if constexpr (Filter::OversamplingFactor == 2)
			{
				if constexpr (Order9)
					return tfdsp::CreateX2Resampler_Chebychev9();
				else
					return tfdsp::CreateX2Resampler_Chebychev7();
			}
			else
			{
				if constexpr (Order9)
				{
					using X4Order9 = tfdsp::X4Resampler<tfdsp::X2Resampler_Order9>;
					return std::make_unique<X4Order9>(
						tfdsp::CreateX2Resampler_Chebychev9);
				}
				else
					return tfdsp::CreateX4Resampler_Cheby7();
			}
		});
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
			output(i) = model.Step(audioValues(i), cutoff, resonance,
				highResonance, driveGain, bass);
		return result;
	}

	template<typename Filter>
	py::array_t<float> RenderArp4072(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		double cutoff, double resonance, double driveGain, double sampleRate)
	{
		const auto audioInfo = audio.request();
		if (audioInfo.ndim != 1)
			throw std::invalid_argument("audio must be a one-dimensional array");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		Filter model([]
		{
			if constexpr (Filter::OversamplingFactor == 1)
				return tfdsp::CreateDummyResampler();
			else if constexpr (Filter::OversamplingFactor == 2)
				return tfdsp::CreateX2Resampler_Chebychev7();
			else
				return tfdsp::CreateX4Resampler_Cheby7();
		});
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
			output(i) = model.Step(audioValues(i), cutoff, resonance,
				driveGain);
		return result;
	}

	template<typename Filter>
	py::array_t<float> RenderArp4072Controls(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> cutoff,
		py::array_t<double, py::array::c_style | py::array::forcecast> resonance,
		double driveGain, double sampleRate)
	{
		const auto audioInfo = audio.request();
		const auto cutoffInfo = cutoff.request();
		const auto resonanceInfo = resonance.request();
		RequireSameSize(audioInfo, cutoffInfo, "audio", "cutoff");
		RequireSameSize(audioInfo, resonanceInfo, "audio", "resonance");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto cutoffValues = cutoff.unchecked<1>();
		auto resonanceValues = resonance.unchecked<1>();
		Filter model([]
		{
			if constexpr (Filter::OversamplingFactor == 1)
				return tfdsp::CreateDummyResampler();
			else if constexpr (Filter::OversamplingFactor == 2)
				return tfdsp::CreateX2Resampler_Chebychev7();
			else
				return tfdsp::CreateX4Resampler_Cheby7();
		});
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
		{
			output(i) = model.StepLogCutoff(audioValues(i), std::log2(std::max(
				cutoffValues(i), std::numeric_limits<double>::min())),
				resonanceValues(i), driveGain);
		}
		return result;
	}

	template<typename Filter>
	py::array_t<float> RenderArp4072ModulatedControls(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> cutoff,
		py::array_t<double, py::array::c_style | py::array::forcecast> linearFm,
		py::array_t<double, py::array::c_style | py::array::forcecast> resonance,
		double driveGain, double sampleRate)
	{
		const auto audioInfo = audio.request();
		const auto cutoffInfo = cutoff.request();
		const auto linearFmInfo = linearFm.request();
		const auto resonanceInfo = resonance.request();
		RequireSameSize(audioInfo, cutoffInfo, "audio", "cutoff");
		RequireSameSize(audioInfo, linearFmInfo, "audio", "linear_fm");
		RequireSameSize(audioInfo, resonanceInfo, "audio", "resonance");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto cutoffValues = cutoff.unchecked<1>();
		auto linearFmValues = linearFm.unchecked<1>();
		auto resonanceValues = resonance.unchecked<1>();
		Filter model([]
		{
			if constexpr (Filter::OversamplingFactor == 1)
				return tfdsp::CreateDummyResampler();
			else if constexpr (Filter::OversamplingFactor == 2)
				return tfdsp::CreateX2Resampler_Chebychev7();
			else
				return tfdsp::CreateX4Resampler_Cheby7();
		});
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
		{
			output(i) = model.StepModulatedLogCutoff(audioValues(i),
				std::log2(std::max(cutoffValues(i),
					std::numeric_limits<double>::min())), linearFmValues(i),
				resonanceValues(i), driveGain);
		}
		return result;
	}

	template<typename Vca>
	py::array_t<float> RenderArp4019(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> linearCv,
		py::array_t<double, py::array::c_style | py::array::forcecast> exponentialCv,
		double initialGain, double sampleRate)
	{
		const auto audioInfo = audio.request();
		const auto linearInfo = linearCv.request();
		const auto exponentialInfo = exponentialCv.request();
		RequireSameSize(audioInfo, linearInfo, "audio", "linear_cv");
		RequireSameSize(audioInfo, exponentialInfo, "audio", "exponential_cv");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto linearValues = linearCv.unchecked<1>();
		auto exponentialValues = exponentialCv.unchecked<1>();
		Vca model([]
		{
			if constexpr (Vca::OversamplingFactor == 1)
				return tfdsp::CreateDummyResampler();
			else if constexpr (Vca::OversamplingFactor == 2)
				return tfdsp::CreateX2Resampler_Chebychev7();
			else
				return tfdsp::CreateX4Resampler_Cheby7();
		});
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
			output(i) = model.Step(audioValues(i), 0.0, linearValues(i),
				exponentialValues(i), initialGain);
		return result;
	}

	template<typename Resampler>
	py::array_t<double> RenderResamplerRoundTrip(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		std::function<std::unique_ptr<Resampler>()> createResampler)
	{
		const auto audioInfo = audio.request();
		if (audioInfo.ndim != 1)
			throw std::invalid_argument("audio must be a one-dimensional array");

		py::array_t<double> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto model = createResampler();
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
			output(i) = model->Downsample(model->Upsample(audioValues(i)));
		return result;
	}

	template<typename Resampler>
	py::array_t<double> RenderUpsampled(
		py::array_t<double, py::array::c_style | py::array::forcecast> input,
		std::function<std::unique_ptr<Resampler>()> createResampler)
	{
		const auto inputInfo = input.request();
		if (inputInfo.ndim != 1)
			throw std::invalid_argument("input must be a one-dimensional array");
		py::array_t<double> result(inputInfo.shape[0] *
			Resampler::ResamplingFactor);
		auto output = result.mutable_unchecked<1>();
		auto values = input.unchecked<1>();
		auto model = createResampler();
		for (py::ssize_t i = 0; i < inputInfo.shape[0]; ++i)
		{
			const auto frame = model->Upsample(values(i));
			for (int j = 0; j < Resampler::ResamplingFactor; ++j)
				output(i * Resampler::ResamplingFactor + j) = frame(j);
		}
		return result;
	}

	template<typename Resampler>
	py::array_t<double> RenderDownsampled(
		py::array_t<double, py::array::c_style | py::array::forcecast> input,
		std::function<std::unique_ptr<Resampler>()> createResampler)
	{
		const auto inputInfo = input.request();
		if (inputInfo.ndim != 1 ||
			inputInfo.shape[0] % Resampler::ResamplingFactor != 0)
		{
			throw std::invalid_argument(
				"input must be one-dimensional and contain complete frames");
		}
		py::array_t<double> result(
			inputInfo.shape[0] / Resampler::ResamplingFactor);
		auto output = result.mutable_unchecked<1>();
		auto values = input.unchecked<1>();
		auto model = createResampler();
		for (py::ssize_t i = 0; i < result.shape(0); ++i)
		{
			Eigen::Array<double, Resampler::ResamplingFactor, 1> frame;
			for (int j = 0; j < Resampler::ResamplingFactor; ++j)
				frame(j) = values(i * Resampler::ResamplingFactor + j);
			output(i) = model->Downsample(frame);
		}
		return result;
	}

	py::array_t<double> RenderArpEnvelope(
		py::array_t<double, py::array::c_style | py::array::forcecast> gate,
		py::array_t<double, py::array::c_style | py::array::forcecast> trigger,
		double attack, double decay, double sustain, double release,
		double curve, int mode, bool autoGateTrigger, double sampleRate)
	{
		const auto gateInfo = gate.request();
		const auto triggerInfo = trigger.request();
		RequireSameSize(gateInfo, triggerInfo, "gate", "trigger");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<double> result(gateInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto gateValues = gate.unchecked<1>();
		auto triggerValues = trigger.unchecked<1>();
		tfdsp::ArpEnvelope envelope;
		envelope.SetSampleRate(sampleRate);
		switch (std::clamp(mode, 0, 2))
		{
		case 1:
			envelope.SetMode(tfdsp::ArpEnvelope::Mode::Ad);
			break;
		case 2:
			envelope.SetMode(tfdsp::ArpEnvelope::Mode::Ar);
			break;
		default:
			envelope.SetMode(tfdsp::ArpEnvelope::Mode::Adsr);
			break;
		}
		for (py::ssize_t i = 0; i < gateInfo.shape[0]; ++i)
		{
			output(i) = envelope.Step(gateValues(i), triggerValues(i), attack,
				decay, sustain, release, curve, autoGateTrigger);
		}
		return result;
	}

	template<typename Filter>
	py::array_t<float> RenderModulatedDiodeLadder(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> driveGain,
		py::array_t<double, py::array::c_style | py::array::forcecast> bass,
		double cutoff, double resonance, bool highResonance, double sampleRate)
	{
		const auto audioInfo = audio.request();
		const auto driveInfo = driveGain.request();
		const auto bassInfo = bass.request();
		RequireSameSize(audioInfo, driveInfo, "audio", "drive_gain");
		RequireSameSize(audioInfo, bassInfo, "audio", "bass");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto driveValues = driveGain.unchecked<1>();
		auto bassValues = bass.unchecked<1>();
		Filter model(tfdsp::CreateX4Resampler_Cheby7);
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
		{
			output(i) = model.Step(audioValues(i), cutoff, resonance,
				highResonance, driveValues(i), bassValues(i));
		}
		return result;
	}

	template<typename Filter>
	py::array_t<float> RenderDiodeLadderVca(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> control,
		double cutoff, double resonance, bool highResonance, double driveGain,
		double bass, double sampleRate, bool oversampledVca)
	{
		const auto audioInfo = audio.request();
		const auto controlInfo = control.request();
		RequireSameSize(audioInfo, controlInfo, "audio", "control");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto controlValues = control.unchecked<1>();
		Filter filter([]
		{
			if constexpr (Filter::OversamplingFactor == 2)
				return tfdsp::CreateX2Resampler_Chebychev7();
			else
				return tfdsp::CreateX4Resampler_Cheby7();
		});
		filter.SetSampleRate(sampleRate);
		tfdsp::Tb303Vca vca;
		vca.SetSampleRate(sampleRate * (oversampledVca ?
			Filter::OversamplingFactor : 1.0));
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
		{
			if (oversampledVca)
			{
				const auto rendered = filter.StepWithPostProcessor(
					audioValues(i), cutoff, resonance, highResonance, driveGain,
					bass, controlValues(i), [&](double audioValue, double vcaControl)
					{
						return vca.Step(audioValue, vcaControl, 0.0);
					});
				output(i) = rendered.postProcessed;
			}
			else
			{
				const double filtered = filter.Step(audioValues(i), cutoff,
					resonance, highResonance, driveGain, bass);
				output(i) = static_cast<float>(vca.Step(filtered,
					controlValues(i), 0.0));
			}
		}
		return result;
	}

	template<typename Filter>
	py::array_t<double> RenderModulatedDiodeLadderVoice(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> cutoff,
		py::array_t<double, py::array::c_style | py::array::forcecast> baseControl,
		py::array_t<double, py::array::c_style | py::array::forcecast> accentControl,
		double resonance, bool highResonance, double driveGain, double bass,
		double sampleRate)
	{
		const auto audioInfo = audio.request();
		const auto cutoffInfo = cutoff.request();
		const auto baseInfo = baseControl.request();
		const auto accentInfo = accentControl.request();
		RequireSameSize(audioInfo, cutoffInfo, "audio", "cutoff");
		RequireSameSize(audioInfo, baseInfo, "audio", "base_control");
		RequireSameSize(audioInfo, accentInfo, "audio", "accent_control");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<double> result({audioInfo.shape[0], py::ssize_t{3}});
		auto output = result.mutable_unchecked<2>();
		auto audioValues = audio.unchecked<1>();
		auto cutoffValues = cutoff.unchecked<1>();
		auto baseValues = baseControl.unchecked<1>();
		auto accentValues = accentControl.unchecked<1>();
		Filter filter([]
		{
			if constexpr (Filter::OversamplingFactor == 2)
				return tfdsp::CreateX2Resampler_Chebychev7();
			else
				return tfdsp::CreateX4Resampler_Cheby7();
		});
		filter.SetSampleRate(sampleRate);
		tfdsp::Tb303Vca vca;
		vca.SetSampleRate(sampleRate * Filter::OversamplingFactor);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
		{
			const auto rendered = filter.StepWithPostProcessor(
				audioValues(i), cutoffValues(i), resonance, highResonance,
				driveGain, bass, baseValues(i),
				[&](double audioValue, double control)
				{
					return vca.Step(audioValue, control, accentValues(i));
				});
			output(i, 0) = rendered.lowPass;
			output(i, 1) = rendered.postProcessed;
			output(i, 2) = static_cast<double>(filter.SolverFailures());
		}
		return result;
	}

	template<typename Oscillator>
	py::array_t<double> RenderTb303Oscillator(
		py::array_t<double, py::array::c_style | py::array::forcecast> pitch,
		py::array_t<double, py::array::c_style | py::array::forcecast> slide,
		py::array_t<double, py::array::c_style | py::array::forcecast> fm,
		py::array_t<double, py::array::c_style | py::array::forcecast> shape,
		py::array_t<double, py::array::c_style | py::array::forcecast> wave,
		double sampleRate, double slideTime, bool linearFm, py::object sync)
	{
		const auto pitchInfo = pitch.request();
		const auto slideInfo = slide.request();
		const auto fmInfo = fm.request();
		const auto shapeInfo = shape.request();
		const auto waveInfo = wave.request();
		RequireSameSize(pitchInfo, slideInfo, "pitch", "slide");
		RequireSameSize(pitchInfo, fmInfo, "pitch", "fm");
		RequireSameSize(pitchInfo, shapeInfo, "pitch", "shape");
		RequireSameSize(pitchInfo, waveInfo, "pitch", "wave");
		py::array_t<double, py::array::c_style | py::array::forcecast> syncArray;
		const double* syncValues = nullptr;
		if (!sync.is_none())
		{
			syncArray = py::cast<py::array_t<double,
				py::array::c_style | py::array::forcecast>>(sync);
			const auto syncInfo = syncArray.request();
			RequireSameSize(pitchInfo, syncInfo, "pitch", "sync");
			syncValues = static_cast<const double*>(syncInfo.ptr);
		}
		if (!(sampleRate > 0.0) || !(slideTime > 0.0))
			throw std::invalid_argument(
				"sample_rate and slide_time must be positive");

		py::array_t<double> result({pitchInfo.shape[0], py::ssize_t{4}});
		auto output = result.mutable_unchecked<2>();
		auto pitchValues = pitch.unchecked<1>();
		auto slideValues = slide.unchecked<1>();
		auto fmValues = fm.unchecked<1>();
		auto shapeValues = shape.unchecked<1>();
		auto waveValues = wave.unchecked<1>();
		Oscillator oscillator([]
		{
			using Resampler = typename Oscillator::Resampler;
			if constexpr (std::is_same_v<Resampler, tfdsp::DummyResampler>)
				return tfdsp::CreateDummyResampler();
			else if constexpr (std::is_same_v<Resampler,
				tfdsp::X2Resampler_Order5>)
				return tfdsp::CreateX2Resampler_Butterworth5();
			else if constexpr (std::is_same_v<Resampler,
				tfdsp::X2Resampler_Order7>)
				return tfdsp::CreateX2Resampler_Chebychev7();
			else if constexpr (std::is_same_v<Resampler,
				tfdsp::X4Resampler<tfdsp::X2Resampler_Order5>>)
				return std::make_unique<Resampler>(
					tfdsp::CreateX2Resampler_Butterworth5);
			else
				return tfdsp::CreateX4Resampler_Cheby7();
		});
		oscillator.SetSampleRate(sampleRate);
		tfdsp::FractionalSchmittTrigger syncTrigger;
		for (py::ssize_t i = 0; i < pitchInfo.shape[0]; ++i)
		{
			const auto syncEvent = syncTrigger.Process(
				syncValues ? syncValues[i] : 0.0);
			const auto value = oscillator.Step(pitchValues(i),
				slideValues(i) >= 1.0, slideTime, 0.0, fmValues(i),
				linearFm, shapeValues(i), waveValues(i),
				syncEvent.triggered ? syncEvent.position : -1.0);
			output(i, 0) = value.saw;
			output(i, 1) = value.square;
			output(i, 2) = value.mixed;
			output(i, 3) = value.pitch;
		}
		return result;
	}

	py::array_t<double> RenderTb303Q8(
		py::array_t<double, py::array::c_style | py::array::forcecast> saw,
		double frequency, double shape, double sampleRate)
	{
		const auto sawInfo = saw.request();
		if (sawInfo.ndim != 1)
			throw std::invalid_argument("saw must be a one-dimensional array");
		if (!(frequency > 0.0) || !(sampleRate > 0.0))
			throw std::invalid_argument("frequency and sample_rate must be positive");
		py::array_t<double> result(sawInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto sawValues = saw.unchecked<1>();
		tfdsp::Tb303SquareShaper shaper;
		shaper.SetSampleRate(sampleRate);
		shaper.Reset();
		for (py::ssize_t i = 0; i < sawInfo.shape[0]; ++i)
			output(i) = shaper.Step(sawValues(i), frequency, shape);
		return result;
	}

	template<typename Oscillator>
	py::array_t<double> RenderWavefoldOscillator(
		py::array_t<double, py::array::c_style | py::array::forcecast> frequency,
		py::array_t<double, py::array::c_style | py::array::forcecast> morph,
		py::array_t<double, py::array::c_style | py::array::forcecast> fold,
		py::array_t<double, py::array::c_style | py::array::forcecast> symmetry,
		double sampleRate, bool adaa, int character)
	{
		const auto frequencyInfo = frequency.request();
		const auto morphInfo = morph.request();
		const auto foldInfo = fold.request();
		const auto symmetryInfo = symmetry.request();
		RequireSameSize(frequencyInfo, morphInfo, "frequency", "morph");
		RequireSameSize(frequencyInfo, foldInfo, "frequency", "fold");
		RequireSameSize(frequencyInfo, symmetryInfo, "frequency", "symmetry");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<double> result(frequencyInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto frequencyValues = frequency.unchecked<1>();
		auto morphValues = morph.unchecked<1>();
		auto foldValues = fold.unchecked<1>();
		auto symmetryValues = symmetry.unchecked<1>();
		Oscillator oscillator([]
		{
			if constexpr (Oscillator::OversamplingFactor == 1)
				return tfdsp::CreateDummyResampler();
			else if constexpr (Oscillator::OversamplingFactor == 2)
				return tfdsp::CreateX2Resampler_Chebychev7();
			else if constexpr (Oscillator::OversamplingFactor == 4)
				return tfdsp::CreateX4Resampler_Cheby7();
			else
				return tfdsp::CreateX16Resampler_Cheby7();
		});
		oscillator.SetSampleRate(sampleRate);
		oscillator.SetFolderAntialiasing(adaa);
		if (character < 0 || character >=
			static_cast<int>(tfdsp::WavefolderCharacter::Count))
			throw std::invalid_argument("invalid wavefolder character");
		oscillator.SetCharacter(
			static_cast<tfdsp::WavefolderCharacter>(character));
		for (py::ssize_t i = 0; i < frequencyInfo.shape[0]; ++i)
			output(i) = oscillator.Step(frequencyValues(i), morphValues(i),
				foldValues(i), symmetryValues(i));
		return result;
	}

	py::array_t<double> RenderStackedOscillator(
		py::array_t<double, py::array::c_style | py::array::forcecast> frequency,
		py::array_t<double, py::array::c_style | py::array::forcecast> pulseWidth,
		int voices, double spreadCents, double pulseMix, double width,
		double sampleRate)
	{
		const auto frequencyInfo = frequency.request();
		const auto pulseWidthInfo = pulseWidth.request();
		RequireSameSize(frequencyInfo, pulseWidthInfo, "frequency", "pulse_width");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");
		const int count = std::clamp(voices, 1,
			tfdsp::MaximumStackedOscillatorVoices);
		const auto pitchPositions =
			tfdsp::StackedOscillatorPitchPositions(count);
		const auto panPositions = tfdsp::StackedOscillatorPanPositions(count);
		std::array<tfdsp::StackedOscillatorVoice,
			tfdsp::MaximumStackedOscillatorVoices> oscillators{};
		constexpr double GoldenConjugate = 0.6180339887498948482;
		for (int voice = 0; voice < count; ++voice)
			oscillators[voice].Reset(std::fmod(voice * GoldenConjugate, 1.0));
		const double normalization = 1.0 / std::sqrt(static_cast<double>(count));
		std::array<double, tfdsp::MaximumStackedOscillatorVoices> leftGains{};
		std::array<double, tfdsp::MaximumStackedOscillatorVoices> rightGains{};
		for (int voice = 0; voice < count; ++voice)
		{
			const double pan = std::clamp(width * panPositions[voice], -1.0, 1.0);
			leftGains[voice] = std::sqrt(1.0 - pan);
			rightGains[voice] = std::sqrt(1.0 + pan);
		}

		py::array_t<double> result({frequencyInfo.shape[0], py::ssize_t{4}});
		auto output = result.mutable_unchecked<2>();
		auto frequencies = frequency.unchecked<1>();
		auto pulseWidths = pulseWidth.unchecked<1>();
		for (py::ssize_t sample = 0; sample < frequencyInfo.shape[0]; ++sample)
		{
			double mono = 0.0;
			double left = 0.0;
			double right = 0.0;
			double sub = 0.0;
			for (int voice = 0; voice < count; ++voice)
			{
				const double voiceFrequency = std::clamp(frequencies(sample) *
					std::exp2(spreadCents * pitchPositions[voice] / 1200.0),
					0.0, 0.45 * sampleRate);
				const auto rendered = oscillators[voice].Step(
					voiceFrequency / sampleRate, pulseWidths(sample), pulseMix);
				mono += rendered.main;
				left += leftGains[voice] * rendered.main;
				right += rightGains[voice] * rendered.main;
				sub += rendered.sub;
			}
			output(sample, 0) = normalization * mono;
			output(sample, 1) = normalization * left;
			output(sample, 2) = normalization * right;
			output(sample, 3) = normalization * sub;
		}
		return result;
	}

	template<typename Oscillator>
	py::array_t<double> RenderWavefolderExternal(
		py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> fold,
		py::array_t<double, py::array::c_style | py::array::forcecast> symmetry,
		double sampleRate, bool adaa, int character)
	{
		const auto audioInfo = audio.request();
		const auto foldInfo = fold.request();
		const auto symmetryInfo = symmetry.request();
		RequireSameSize(audioInfo, foldInfo, "audio", "fold");
		RequireSameSize(audioInfo, symmetryInfo, "audio", "symmetry");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");
		if (character < 0 || character >=
			static_cast<int>(tfdsp::WavefolderCharacter::Count))
			throw std::invalid_argument("invalid wavefolder character");

		py::array_t<double> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto foldValues = fold.unchecked<1>();
		auto symmetryValues = symmetry.unchecked<1>();
		Oscillator oscillator([]
		{
			if constexpr (Oscillator::OversamplingFactor == 2)
				return tfdsp::CreateX2Resampler_Chebychev7();
			else if constexpr (Oscillator::OversamplingFactor == 4)
				return tfdsp::CreateX4Resampler_Cheby7();
			else
				return tfdsp::CreateX16Resampler_Cheby7();
		});
		oscillator.SetSampleRate(sampleRate);
		oscillator.SetFolderAntialiasing(adaa);
		oscillator.SetCharacter(
			static_cast<tfdsp::WavefolderCharacter>(character));
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
		{
			output(i) = oscillator.StepWithInput(261.625565, 0.0,
				foldValues(i), symmetryValues(i), audioValues(i), true).folded;
		}
		return result;
	}

	py::array_t<double> EvaluateWavefolderFunction(
		py::array_t<double, py::array::c_style | py::array::forcecast> input,
		int character, bool primitive)
	{
		const auto info = input.request();
		if (info.ndim != 1)
			throw std::invalid_argument("input must be a one-dimensional array");
		py::array_t<double> result(info.shape[0]);
		if (character < 0 || character >=
			static_cast<int>(tfdsp::WavefolderCharacter::Count))
			throw std::invalid_argument("invalid wavefolder character");
		const auto selected =
			static_cast<tfdsp::WavefolderCharacter>(character);
		auto output = result.mutable_unchecked<1>();
		auto values = input.unchecked<1>();
		for (py::ssize_t i = 0; i < info.shape[0]; ++i)
			output(i) = primitive ?
				tfdsp::Wavefolder::Primitive(values(i), selected) :
				tfdsp::Wavefolder::Transfer(values(i), selected);
		return result;
	}

	py::array_t<double> EvaluateWavefolderAdaa(
		py::array_t<double, py::array::c_style | py::array::forcecast> input,
		int character)
	{
		const auto info = input.request();
		if (info.ndim != 1)
			throw std::invalid_argument("input must be a one-dimensional array");
		if (character < 0 || character >=
			static_cast<int>(tfdsp::WavefolderCharacter::Count))
			throw std::invalid_argument("invalid wavefolder character");
		py::array_t<double> result(info.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto values = input.unchecked<1>();
		tfdsp::Wavefolder folder;
		const auto selected =
			static_cast<tfdsp::WavefolderCharacter>(character);
		for (py::ssize_t i = 0; i < info.shape[0]; ++i)
			output(i) = folder.Process(values(i), selected);
		return result;
	}
}

PYBIND11_MODULE(_triggerfish_dsp, module)
{
	using CrashFit = tfdsp::percussion::CrashCymbalFitParameters;
	py::class_<CrashFit>(module, "CrashCymbalFitParameters")
		.def(py::init<>())
		.def_readwrite("sparse_frequency_hz", &CrashFit::sparseFrequencyHz)
		.def_readwrite("sparse_decay_seconds", &CrashFit::sparseDecaySeconds)
		.def_readwrite("sparse_amplitude", &CrashFit::sparseAmplitude)
		.def_readwrite("sparse_phase_radians", &CrashFit::sparsePhaseRadians)
		.def_readwrite("sparse_tune", &CrashFit::sparseTune)
		.def_readwrite("sparse_decay_scale", &CrashFit::sparseDecayScale)
		.def_readwrite("dense_minimum_frequency_hz",
			&CrashFit::denseMinimumFrequencyHz)
		.def_readwrite("dense_maximum_frequency_hz",
			&CrashFit::denseMaximumFrequencyHz)
		.def_readwrite("dense_frequency_warp", &CrashFit::denseFrequencyWarp)
		.def_readwrite("dense_spacing_jitter", &CrashFit::denseSpacingJitter)
		.def_readwrite("dense_low_decay_seconds", &CrashFit::denseLowDecaySeconds)
		.def_readwrite("dense_high_decay_seconds", &CrashFit::denseHighDecaySeconds)
		.def_readwrite("dense_decay_curve", &CrashFit::denseDecayCurve)
		.def_readwrite("dense_decay_envelope_octaves",
			&CrashFit::denseDecayEnvelopeOctaves)
		.def_readwrite("dense_decay_spread_octaves",
			&CrashFit::denseDecaySpreadOctaves)
		.def_readwrite("dense_tilt_db_per_octave",
			&CrashFit::denseTiltDbPerOctave)
		.def_readwrite("dense_gain_envelope_db", &CrashFit::denseGainEnvelopeDb)
		.def_readwrite("dense_gain_spread_db", &CrashFit::denseGainSpreadDb)
		.def_readwrite("dense_mode_seed", &CrashFit::denseModeSeed)
		.def_readwrite("turbulence_low_gain", &CrashFit::turbulenceLowGain)
		.def_readwrite("turbulence_middle_gain", &CrashFit::turbulenceMiddleGain)
		.def_readwrite("turbulence_high_gain", &CrashFit::turbulenceHighGain)
		.def_readwrite("turbulence_low_decay_seconds",
			&CrashFit::turbulenceLowDecaySeconds)
		.def_readwrite("turbulence_middle_decay_seconds",
			&CrashFit::turbulenceMiddleDecaySeconds)
		.def_readwrite("turbulence_high_decay_seconds",
			&CrashFit::turbulenceHighDecaySeconds)
		.def_readwrite("dispersion_feedback", &CrashFit::dispersionFeedback)
		.def_readwrite("dispersion_drive", &CrashFit::dispersionDrive)
		.def_readwrite("dispersion_excursion_samples",
			&CrashFit::dispersionExcursionSamples)
		.def_readwrite("dispersion_low_decay_seconds",
			&CrashFit::dispersionLowDecaySeconds)
		.def_readwrite("dispersion_middle_decay_seconds",
			&CrashFit::dispersionMiddleDecaySeconds)
		.def_readwrite("dispersion_high_decay_seconds",
			&CrashFit::dispersionHighDecaySeconds)
		.def_readwrite("contact_duration_scale", &CrashFit::contactDurationScale)
		.def_readwrite("contact_pulse_gain", &CrashFit::contactPulseGain)
		.def_readwrite("contact_chirp_gain", &CrashFit::contactChirpGain)
		.def_readwrite("contact_chirp_frequency_scale",
			&CrashFit::contactChirpFrequencyScale)
		.def_readwrite("contact_noise_gain", &CrashFit::contactNoiseGain)
		.def_readwrite("contact_noise_duration_scale",
			&CrashFit::contactNoiseDurationScale)
		.def_readwrite("contact_noise_tilt_db", &CrashFit::contactNoiseTiltDb)
		.def_readwrite("contact_micro_gain", &CrashFit::contactMicroGain)
		.def_readwrite("contact_micro_duration_scale",
			&CrashFit::contactMicroDurationScale)
		.def_readwrite("contact_micro_density_scale",
			&CrashFit::contactMicroDensityScale)
		.def_readwrite("direct_gain", &CrashFit::directGain)
		.def_readwrite("sparse_gain", &CrashFit::sparseGain)
		.def_readwrite("dense_gain", &CrashFit::denseGain)
		.def_readwrite("sparse_bloom_gain", &CrashFit::sparseBloomGain)
		.def_readwrite("body_bypass_gain", &CrashFit::bodyBypassGain)
		.def_readwrite("output_gain", &CrashFit::outputGain)
		.def_readwrite("colour_frequency_hz", &CrashFit::colourFrequencyHz)
		.def_readwrite("colour_gain_db", &CrashFit::colourGainDb)
		.def_readwrite("high_cut_hz", &CrashFit::highCutHz)
		.def_readwrite("strength_gamma", &CrashFit::strengthGamma)
		.def_readwrite("body_strength_gamma", &CrashFit::bodyStrengthGamma)
		.def_readwrite("dense_strength_gamma", &CrashFit::denseStrengthGamma)
		.def_readwrite("dense_velocity_loss_nepers_per_second",
			&CrashFit::denseVelocityLossNepersPerSecond);
	module.def("render_crash", [](const py::ssize_t sampleCount,
		double sampleRate, float strength, float location, float hardness,
		std::uint32_t seed, const CrashFit &fit)
	{
		if (sampleCount <= 0 || !(sampleRate > 0.0))
			throw std::invalid_argument("crash render dimensions must be positive");
		tfdsp::percussion::CrashCymbal cymbal;
		const float rate = static_cast<float>(sampleRate);
		py::array_t<float> result(sampleCount);
		auto *output = result.mutable_data();
		{
			py::gil_scoped_release release;
			cymbal.Prepare(rate,
				tfdsp::percussion::DefaultCrashCymbalParameters(rate, fit));
			cymbal.Trigger({strength, location, hardness, seed});
			for (py::ssize_t sample = 0; sample < sampleCount; ++sample)
				output[sample] = cymbal.Process();
		}
		return result;
	}, py::arg("sample_count"), py::arg("sample_rate") = 48000.0,
		py::arg("strength") = .8f, py::arg("location") = 1.f,
		py::arg("hardness") = .65f, py::arg("seed") = 1u,
		py::arg("fit") = CrashFit{});
	module.def("render_crash_components", [](const py::ssize_t sampleCount,
		double sampleRate, float strength, float location, float hardness,
		std::uint32_t seed, const CrashFit &fit)
	{
		if (sampleCount <= 0 || !(sampleRate > 0.0))
			throw std::invalid_argument("crash render dimensions must be positive");
		tfdsp::percussion::CrashCymbal cymbal;
		const float rate = static_cast<float>(sampleRate);
		cymbal.Prepare(rate,
			tfdsp::percussion::DefaultCrashCymbalParameters(rate, fit));
		cymbal.Trigger({strength, location, hardness, seed});
		py::array_t<float> result({sampleCount, py::ssize_t{5}});
		auto output = result.mutable_unchecked<2>();
		for (py::ssize_t sample = 0; sample < sampleCount; ++sample) {
			const auto frame = cymbal.ProcessFrame();
			output(sample, 0) = frame.directContact;
			output(sample, 1) = frame.dispersion;
			output(sample, 2) = frame.sparseModes;
			output(sample, 3) = frame.denseResidual;
			output(sample, 4) = frame.output;
		}
		return result;
	}, py::arg("sample_count"), py::arg("sample_rate") = 48000.0,
		py::arg("strength") = .8f, py::arg("location") = 1.f,
		py::arg("hardness") = .65f, py::arg("seed") = 1u,
		py::arg("fit") = CrashFit{});
	module.def("render_crash_sequence", [](const py::ssize_t sampleCount,
		double sampleRate, py::array_t<float, py::array::c_style |
		py::array::forcecast> strengths,
		py::array_t<float, py::array::c_style |
		py::array::forcecast> locations,
		py::array_t<float, py::array::c_style |
		py::array::forcecast> hardnesses,
		py::array_t<py::ssize_t, py::array::c_style |
		py::array::forcecast> onsets,
		py::array_t<std::uint32_t, py::array::c_style |
		py::array::forcecast> seeds, const CrashFit &fit)
	{
		if (sampleCount <= 0 || !(sampleRate > 0.0))
			throw std::invalid_argument("crash render dimensions must be positive");
		const auto count = strengths.size();
		if (strengths.ndim() != 1 || locations.ndim() != 1 ||
			hardnesses.ndim() != 1 || onsets.ndim() != 1 || seeds.ndim() != 1 ||
			locations.size() != count || hardnesses.size() != count ||
			onsets.size() != count || seeds.size() != count)
			throw std::invalid_argument("crash event arrays must be equal-length vectors");
		auto strength = strengths.unchecked<1>();
		auto location = locations.unchecked<1>();
		auto hardness = hardnesses.unchecked<1>();
		auto onset = onsets.unchecked<1>();
		auto seed = seeds.unchecked<1>();
		for (py::ssize_t event = 0; event < count; ++event)
			if (onset(event) < 0 || onset(event) >= sampleCount ||
				(event > 0 && onset(event) < onset(event - 1)))
				throw std::invalid_argument(
					"crash event onsets must be ordered and inside the render");
		tfdsp::percussion::CrashCymbal cymbal;
		const float rate = static_cast<float>(sampleRate);
		cymbal.Prepare(rate,
			tfdsp::percussion::DefaultCrashCymbalParameters(rate, fit));
		py::array_t<float> result(sampleCount);
		auto output = result.mutable_unchecked<1>();
		py::ssize_t event = 0;
		for (py::ssize_t sample = 0; sample < sampleCount; ++sample) {
			while (event < count && onset(event) == sample) {
				cymbal.Trigger({strength(event), location(event), hardness(event),
					seed(event)});
				++event;
			}
			output(sample) = cymbal.Process();
		}
		return result;
	}, py::arg("sample_count"), py::arg("sample_rate"),
		py::arg("strengths"), py::arg("locations"), py::arg("hardnesses"),
		py::arg("onsets"), py::arg("seeds"), py::arg("fit") = CrashFit{});
	module.doc() = "TriggerFish DSP development bindings";
	module.def("late_reverb_wall_impulse", &RenderLateReverbWallImpulse,
		py::arg("sample_count"), py::arg("sample_rate") = 48'000.0,
		py::arg("space") = 0.5, py::arg("aspect") = 0.5,
		py::arg("decay") = 0.0, py::arg("damping") = 0.45,
		py::arg("diffusion") = 0.75, py::arg("optimized") = false);

	module.def("diode_ladder_map_cutoff", [](double requestedHz,
		double maximumHz)
	{
		return tfdsp::DiodeLadderFilter<
			tfdsp::X4Resampler_Order7>::MapCutoffControl(
			requestedHz, maximumHz);
	}, py::arg("requested_hz"), py::arg("maximum_hz") = 20000.0);

	module.def("vca_transistor", [](py::array_t<float, py::array::c_style | py::array::forcecast> audio,
		py::array_t<float, py::array::c_style | py::array::forcecast> cv,
		float sampleRate, float finalGain)
	{
		const auto audioInfo = audio.request();
		const auto cvInfo = cv.request();
		RequireSameSize(audioInfo, cvInfo, "audio", "cv");
		if (!(sampleRate > 0.0f))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto cvValues = cv.unchecked<1>();
		VCA_TransistorCore<tfdsp::X2Resampler_Order7> model(tfdsp::CreateX2Resampler_Chebychev7);
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
			output(i) = model.Step(audioValues(i), cvValues(i), finalGain);
		return result;
	}, py::arg("audio"), py::arg("cv"), py::arg("sample_rate"), py::arg("final_gain") = 1.0f);

	module.def("vca_ota_legacy", [](py::array_t<float, py::array::c_style | py::array::forcecast> audio,
		py::array_t<float, py::array::c_style | py::array::forcecast> cv,
		float sampleRate, float finalGain)
	{
		const auto audioInfo = audio.request();
		const auto cvInfo = cv.request();
		RequireSameSize(audioInfo, cvInfo, "audio", "cv");
		if (!(sampleRate > 0.0f))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto cvValues = cv.unchecked<1>();
		VCA_OTACore<tfdsp::X2Resampler_Order7> model(tfdsp::CreateX2Resampler_Chebychev7);
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
			output(i) = model.Step(audioValues(i), cvValues(i), finalGain);
		return result;
	}, py::arg("audio"), py::arg("cv"), py::arg("sample_rate"),
		py::arg("final_gain") = 1.0f);

	module.def("ota_vca_current", [](py::array_t<double, py::array::c_style | py::array::forcecast> differential,
		py::array_t<double, py::array::c_style | py::array::forcecast> control,
		double efficiency, double inputOffset, double mirrorImbalance)
	{
		const auto differentialInfo = differential.request();
		const auto controlInfo = control.request();
		RequireSameSize(differentialInfo, controlInfo, "differential", "control");
		tfdsp::OtaVcaCore::Configuration configuration;
		configuration.currentTransferEfficiency = efficiency;
		configuration.inputOffsetVolts = inputOffset;
		configuration.mirrorImbalance = mirrorImbalance;
		tfdsp::OtaVcaCore model(configuration);
		py::array_t<double> result(differentialInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto differentialValues = differential.unchecked<1>();
		auto controlValues = control.unchecked<1>();
		for (py::ssize_t i = 0; i < differentialInfo.shape[0]; ++i)
			output(i) = model.ProcessCurrent(differentialValues(i), controlValues(i));
		return result;
	}, py::arg("differential"), py::arg("control"),
		py::arg("efficiency") = 0.85, py::arg("input_offset") = 0.0,
		py::arg("mirror_imbalance") = 0.0);

	module.def("tb303_vca", [](py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> baseControl,
		py::array_t<double, py::array::c_style | py::array::forcecast> accentControl,
		double sampleRate)
	{
		const auto audioInfo = audio.request();
		const auto baseInfo = baseControl.request();
		const auto accentInfo = accentControl.request();
		RequireSameSize(audioInfo, baseInfo, "audio", "base_control");
		RequireSameSize(audioInfo, accentInfo, "audio", "accent_control");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<double> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto baseValues = baseControl.unchecked<1>();
		auto accentValues = accentControl.unchecked<1>();
		tfdsp::Tb303Vca model;
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
			output(i) = model.Step(audioValues(i), baseValues(i), accentValues(i));
		return result;
	}, py::arg("audio"), py::arg("base_control"), py::arg("accent_control"),
		py::arg("sample_rate") = 48000.0);

	module.def("tb303_articulation", [](py::array_t<double, py::array::c_style | py::array::forcecast> gate,
		py::array_t<double, py::array::c_style | py::array::forcecast> accent,
		py::array_t<double, py::array::c_style | py::array::forcecast> resonance,
		double normalDecay, double accentDecay, double vcaDecay,
		double sampleRate, bool devilFish, int accentSweep)
	{
		const auto gateInfo = gate.request();
		const auto accentInfo = accent.request();
		const auto resonanceInfo = resonance.request();
		RequireSameSize(gateInfo, accentInfo, "gate", "accent");
		RequireSameSize(gateInfo, resonanceInfo, "gate", "resonance");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<double> result({gateInfo.shape[0], py::ssize_t{4}});
		auto output = result.mutable_unchecked<2>();
		auto gateValues = gate.unchecked<1>();
		auto accentValues = accent.unchecked<1>();
		auto resonanceValues = resonance.unchecked<1>();
		tfdsp::Tb303Articulation model;
		model.SetSampleRate(sampleRate);
		model.SetMode(devilFish ? tfdsp::Tb303Articulation::Mode::DevilFish :
			tfdsp::Tb303Articulation::Mode::Stock);
		model.SetAccentSweepMode(static_cast<tfdsp::Tb303AccentSweep::Mode>(
			std::clamp(accentSweep, 0, 3)));
		for (py::ssize_t i = 0; i < gateInfo.shape[0]; ++i)
		{
			const auto value = model.Step(gateValues(i), accentValues(i),
				resonanceValues(i), normalDecay, accentDecay, vcaDecay);
			output(i, 0) = value.mainEnvelope;
			output(i, 1) = value.filterAccent;
			output(i, 2) = value.volumeEnvelope;
			output(i, 3) = value.vcaAccent;
		}
		return result;
	}, py::arg("gate"), py::arg("accent"), py::arg("resonance"),
		py::arg("normal_decay") = 0.5, py::arg("accent_decay") = 0.2,
		py::arg("vca_decay") = 0.5, py::arg("sample_rate") = 48000.0,
		py::arg("devil_fish") = false, py::arg("accent_sweep") = 2);

	module.def("vdpo_bdf", [](py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> damping,
		py::array_t<double, py::array::c_style | py::array::forcecast> angularFrequency,
		double sampleRate)
	{
		const auto audioInfo = audio.request();
		const auto dampingInfo = damping.request();
		const auto frequencyInfo = angularFrequency.request();
		RequireSameSize(audioInfo, dampingInfo, "audio", "damping");
		RequireSameSize(audioInfo, frequencyInfo, "audio", "angular_frequency");
		if (!(sampleRate > 0.0))
			throw std::invalid_argument("sample_rate must be positive");

		py::array_t<float> result(audioInfo.shape[0]);
		auto output = result.mutable_unchecked<1>();
		auto audioValues = audio.unchecked<1>();
		auto dampingValues = damping.unchecked<1>();
		auto frequencyValues = angularFrequency.unchecked<1>();
		VdpOscillator<tfdsp::X4Resampler_Order7, 3> model(tfdsp::CreateX4Resampler_Cheby7);
		model.SetSampleRate(sampleRate);
		for (py::ssize_t i = 0; i < audioInfo.shape[0]; ++i)
			output(i) = model.Step(audioValues(i), dampingValues(i), frequencyValues(i));
		return result;
	}, py::arg("audio"), py::arg("damping"), py::arg("angular_frequency"), py::arg("sample_rate"));

	using BaselineVdpo = VdpSplitOscillator<tfdsp::X4Resampler_Order7>;

	module.def("vdpo", [](py::array_t<double, py::array::c_style | py::array::forcecast> audio,
		py::array_t<double, py::array::c_style | py::array::forcecast> damping,
		py::array_t<double, py::array::c_style | py::array::forcecast> angularFrequency,
		double sampleRate)
	{
		return RenderVdpo<BaselineVdpo>(audio, damping, angularFrequency, sampleRate);
	}, py::arg("audio"), py::arg("damping"), py::arg("angular_frequency"), py::arg("sample_rate"));

	module.def("vdpo_pitch_std", &RenderVdpoPitch<BaselineVdpo, false>,
		py::arg("audio"), py::arg("damping"), py::arg("pitch"), py::arg("sample_rate"));
	module.def("vdpo_pitch_fast_exp2", &RenderVdpoPitch<BaselineVdpo, true>,
		py::arg("audio"), py::arg("damping"), py::arg("pitch"), py::arg("sample_rate"));

	module.def("detune_legacy_double", &RenderDetune<DetuneMethod::LegacyDouble>,
		py::arg("pitch"), py::arg("detune"), py::arg("reference_frequency") = 261.63);
	module.def("detune_optimized", &RenderDetune<DetuneMethod::Optimized>,
		py::arg("pitch"), py::arg("detune"), py::arg("reference_frequency") = 261.63);
	module.def("detune_reference_double", &RenderDetune<DetuneMethod::ReferenceDouble>,
		py::arg("pitch"), py::arg("detune"), py::arg("reference_frequency") = 261.63);

	using DiodeLadderX1 = tfdsp::DiodeLadderFilter<tfdsp::DummyResampler>;
	using DiodeLadderX2 = tfdsp::DiodeLadderFilter<tfdsp::X2Resampler_Order7>;
	using DiodeLadderX4 = tfdsp::DiodeLadderFilter<tfdsp::X4Resampler_Order7>;
	using DiodeLadderX2Order9 =
		tfdsp::DiodeLadderFilter<tfdsp::X2Resampler_Order9>;
	using X4ResamplerOrder9 = tfdsp::X4Resampler<tfdsp::X2Resampler_Order9>;
	using DiodeLadderX4Order9 = tfdsp::DiodeLadderFilter<X4ResamplerOrder9>;
	module.def("diode_ladder_controls_x1",
		&RenderDiodeLadderControls<DiodeLadderX1>, py::arg("audio"),
		py::arg("cutoff"), py::arg("linear_fm"), py::arg("resonance"),
		py::arg("high_resonance") = false, py::arg("drive_gain") = 1.0,
		py::arg("bass") = 0.0, py::arg("sample_rate") = 48000.0);
	module.def("diode_ladder_x1", &RenderDiodeLadder<DiodeLadderX1>,
		py::arg("audio"), py::arg("cutoff"), py::arg("resonance") = 0.0,
		py::arg("high_resonance") = false, py::arg("drive_gain") = 1.0,
		py::arg("bass") = 0.0, py::arg("sample_rate") = 48000.0);
	module.def("diode_ladder_x2", &RenderDiodeLadder<DiodeLadderX2>,
		py::arg("audio"), py::arg("cutoff"), py::arg("resonance") = 0.0,
		py::arg("high_resonance") = false, py::arg("drive_gain") = 1.0,
		py::arg("bass") = 0.0, py::arg("sample_rate") = 48000.0);
	module.def("diode_ladder_x4", &RenderDiodeLadder<DiodeLadderX4>,
		py::arg("audio"), py::arg("cutoff"), py::arg("resonance") = 0.0,
		py::arg("high_resonance") = false, py::arg("drive_gain") = 1.0,
		py::arg("bass") = 0.0, py::arg("sample_rate") = 48000.0);
	module.def("diode_ladder_modulated_x4",
		&RenderModulatedDiodeLadder<DiodeLadderX4>, py::arg("audio"),
		py::arg("drive_gain"), py::arg("bass"), py::arg("cutoff"),
		py::arg("resonance") = 0.0, py::arg("high_resonance") = false,
		py::arg("sample_rate") = 48000.0);
	module.def("diode_ladder_controls_x4",
		&RenderDiodeLadderControls<DiodeLadderX4>, py::arg("audio"),
		py::arg("cutoff"), py::arg("linear_fm"), py::arg("resonance"),
		py::arg("high_resonance") = false, py::arg("drive_gain") = 1.0,
		py::arg("bass") = 0.0, py::arg("sample_rate") = 48000.0);
	module.def("diode_ladder_diagnostics_x4",
		&RenderDiodeLadderDiagnostics<DiodeLadderX4>, py::arg("audio"),
		py::arg("cutoff"), py::arg("linear_fm"), py::arg("resonance"),
		py::arg("high_resonance") = false, py::arg("drive_gain") = 1.0,
		py::arg("bass") = 0.0, py::arg("sample_rate") = 48000.0);
	module.def("diode_ladder_x2_order9",
		&RenderDiodeLadder<DiodeLadderX2Order9, true>,
		py::arg("audio"), py::arg("cutoff"), py::arg("resonance") = 0.0,
		py::arg("high_resonance") = false, py::arg("drive_gain") = 1.0,
		py::arg("bass") = 0.0, py::arg("sample_rate") = 48000.0);
	module.def("diode_ladder_x4_order9",
		&RenderDiodeLadder<DiodeLadderX4Order9, true>,
		py::arg("audio"), py::arg("cutoff"), py::arg("resonance") = 0.0,
		py::arg("high_resonance") = false, py::arg("drive_gain") = 1.0,
		py::arg("bass") = 0.0, py::arg("sample_rate") = 48000.0);

	using Arp4072X1 = tfdsp::Arp4072Filter<tfdsp::DummyResampler>;
	using Arp4072X2 = tfdsp::Arp4072Filter<tfdsp::X2Resampler_Order7>;
	using Arp4072X4 = tfdsp::Arp4072Filter<tfdsp::X4Resampler_Order7>;
	module.def("arp4072_controls_x1", &RenderArp4072Controls<Arp4072X1>,
		py::arg("audio"), py::arg("cutoff"), py::arg("resonance"),
		py::arg("drive_gain") = 1.0, py::arg("sample_rate") = 48000.0);
	module.def("arp4072_x1", &RenderArp4072<Arp4072X1>, py::arg("audio"),
		py::arg("cutoff"), py::arg("resonance") = 0.0,
		py::arg("drive_gain") = 1.0, py::arg("sample_rate") = 48000.0);
	module.def("arp4072_x2", &RenderArp4072<Arp4072X2>, py::arg("audio"),
		py::arg("cutoff"), py::arg("resonance") = 0.0,
		py::arg("drive_gain") = 1.0, py::arg("sample_rate") = 48000.0);
	module.def("arp4072_x4", &RenderArp4072<Arp4072X4>, py::arg("audio"),
		py::arg("cutoff"), py::arg("resonance") = 0.0,
		py::arg("drive_gain") = 1.0, py::arg("sample_rate") = 48000.0);
	module.def("arp4072_controls_x4", &RenderArp4072Controls<Arp4072X4>,
		py::arg("audio"), py::arg("cutoff"), py::arg("resonance"),
		py::arg("drive_gain") = 1.0, py::arg("sample_rate") = 48000.0);
	module.def("arp4072_modulated_controls_x1",
		&RenderArp4072ModulatedControls<Arp4072X1>, py::arg("audio"),
		py::arg("cutoff"), py::arg("linear_fm"), py::arg("resonance"),
		py::arg("drive_gain") = 1.0, py::arg("sample_rate") = 48000.0);
	module.def("arp4072_modulated_controls_x4",
		&RenderArp4072ModulatedControls<Arp4072X4>, py::arg("audio"),
		py::arg("cutoff"), py::arg("linear_fm"), py::arg("resonance"),
		py::arg("drive_gain") = 1.0, py::arg("sample_rate") = 48000.0);
	module.def("arp4072_circuit_values", []
	{
		py::dict values;
		values["audio_base_scale"] = Arp4072X4::AudioBaseScale();
		values["feedback_base_scale"] = Arp4072X4::FeedbackBaseScale();
		values["audio_base_at_5v"] = Arp4072X4::AudioBaseVolts(5.0);
		values["feedback_base_at_5v"] = Arp4072X4::FeedbackBaseVolts(5.0);
		values["limiter_tail_current_amps"] =
			Arp4072X4::LimiterTailCurrentAmps();
		values["limiter_equivalent_peak_volts"] =
			Arp4072X4::LimiterEquivalentPeakVolts();
		values["nominal_limiter_equivalent_peak_volts"] =
			Arp4072X4::NominalLimiterEquivalentPeakVolts();
		values["limiter_gain_calibration"] =
			Arp4072X4::LimiterGainCalibration();
		values["limiter_differential_resistance_ohms"] =
			Arp4072X4::LimiterDifferentialResistanceOhms();
		values["stage_saturation_coefficient_per_volt"] =
			Arp4072X4::StageSaturationCoefficientPerVolt();
		values["stage_base_resistance_ohms"] =
			Arp4072X4::StageBaseResistanceOhms();
		values["output_level_shift_gain"] = Arp4072X4::OutputLevelShiftGain;
		values["small_signal_input_gain"] =
			Arp4072X4::SmallSignalInputGain();
		values["small_signal_feedback_gain"] =
			Arp4072X4::SmallSignalFeedbackGain();
		return values;
	});

	using Arp4019X1 = tfdsp::Arp4019Vca<tfdsp::DummyResampler>;
	using Arp4019X4 = tfdsp::Arp4019Vca<tfdsp::X4Resampler_Order7>;
	module.def("arp4019_x1", &RenderArp4019<Arp4019X1>, py::arg("audio"),
		py::arg("linear_cv"), py::arg("exponential_cv"),
		py::arg("initial_gain") = 0.0, py::arg("sample_rate") = 48000.0);
	module.def("arp4019_x4", &RenderArp4019<Arp4019X4>, py::arg("audio"),
		py::arg("linear_cv"), py::arg("exponential_cv"),
		py::arg("initial_gain") = 0.0, py::arg("sample_rate") = 48000.0);
	module.def("arp4019_circuit_values", []
	{
		py::dict values;
		values["audio_input_scale"] = Arp4019X4::AudioInputScale();
		values["audio_series_resistance_ohms"] =
			Arp4019X4::AudioSeriesResistanceOhms;
		values["unity_control_current_amps"] =
			Arp4019X4::UnityControlCurrentAmps();
		values["output_feedback_resistance_ohms"] =
			Arp4019X4::OutputFeedbackResistanceOhms;
		values["output_bandwidth_hz"] = Arp4019X4::OutputBandwidthHz;
		values["linear_unity_control_volts"] =
			Arp4019X4::LinearUnityControlVolts;
		values["exponential_decibels_per_volt"] =
			Arp4019X4::ExponentialDecibelsPerVolt;
		return values;
	});
	module.def("arp_envelope", &RenderArpEnvelope, py::arg("gate"),
		py::arg("trigger"), py::arg("attack"), py::arg("decay"),
		py::arg("sustain"), py::arg("release"), py::arg("curve") = 0.0,
		py::arg("mode") = 0, py::arg("auto_gate_trigger") = true,
		py::arg("sample_rate") = 48000.0);

	module.def("resampler_round_trip_x2_order7", [](py::array_t<double,
		py::array::c_style | py::array::forcecast> audio)
	{
		return RenderResamplerRoundTrip<tfdsp::X2Resampler_Order7>(audio,
			tfdsp::CreateX2Resampler_Chebychev7);
	}, py::arg("audio"));
	module.def("resampler_round_trip_x2_order9", [](py::array_t<double,
		py::array::c_style | py::array::forcecast> audio)
	{
		return RenderResamplerRoundTrip<tfdsp::X2Resampler_Order9>(audio,
			tfdsp::CreateX2Resampler_Chebychev9);
	}, py::arg("audio"));
	module.def("resampler_round_trip_x4_order7", [](py::array_t<double,
		py::array::c_style | py::array::forcecast> audio)
	{
		return RenderResamplerRoundTrip<tfdsp::X4Resampler_Order7>(audio,
			tfdsp::CreateX4Resampler_Cheby7);
	}, py::arg("audio"));
	module.def("resampler_upsample_x4_order7", [](py::array_t<double,
		py::array::c_style | py::array::forcecast> input)
	{
		return RenderUpsampled<tfdsp::X4Resampler_Order7>(input,
			tfdsp::CreateX4Resampler_Cheby7);
	}, py::arg("input"));
	module.def("resampler_downsample_x4_order7", [](py::array_t<double,
		py::array::c_style | py::array::forcecast> input)
	{
		return RenderDownsampled<tfdsp::X4Resampler_Order7>(input,
			tfdsp::CreateX4Resampler_Cheby7);
	}, py::arg("input"));
	module.def("resampler_round_trip_x4_order9", [](py::array_t<double,
		py::array::c_style | py::array::forcecast> audio)
	{
		return RenderResamplerRoundTrip<X4ResamplerOrder9>(audio, []
		{
			return std::make_unique<X4ResamplerOrder9>(
				tfdsp::CreateX2Resampler_Chebychev9);
		});
	}, py::arg("audio"));

	module.def("diode_ladder_vca_x2", &RenderDiodeLadderVca<DiodeLadderX2>,
		py::arg("audio"), py::arg("control"), py::arg("cutoff"),
		py::arg("resonance") = 0.0, py::arg("high_resonance") = false,
		py::arg("drive_gain") = 1.0, py::arg("bass") = 0.0,
		py::arg("sample_rate") = 48000.0, py::arg("oversampled_vca") = true);
	module.def("diode_ladder_vca_x4", &RenderDiodeLadderVca<DiodeLadderX4>,
		py::arg("audio"), py::arg("control"), py::arg("cutoff"),
		py::arg("resonance") = 0.0, py::arg("high_resonance") = false,
		py::arg("drive_gain") = 1.0, py::arg("bass") = 0.0,
		py::arg("sample_rate") = 48000.0, py::arg("oversampled_vca") = true);
	module.def("diode_ladder_voice_x4",
		&RenderModulatedDiodeLadderVoice<DiodeLadderX4>, py::arg("audio"),
		py::arg("cutoff"), py::arg("base_control"),
		py::arg("accent_control"), py::arg("resonance") = 0.0,
		py::arg("high_resonance") = false, py::arg("drive_gain") = 1.0,
		py::arg("bass") = 0.0, py::arg("sample_rate") = 48000.0);

	using Tb303OscillatorX1 = tfdsp::Tb303Oscillator<tfdsp::DummyResampler>;
	using Tb303OscillatorX2 =
		tfdsp::Tb303Oscillator<tfdsp::X2Resampler_Order7>;
	using Tb303OscillatorX2Order5 =
		tfdsp::Tb303Oscillator<tfdsp::X2Resampler_Order5>;
	using Tb303OscillatorX4 =
		tfdsp::Tb303Oscillator<tfdsp::X4Resampler_Order7>;
	using Tb303OscillatorX4Order5 = tfdsp::Tb303Oscillator<
		tfdsp::X4Resampler<tfdsp::X2Resampler_Order5>>;
	module.def("tb303_q8", &RenderTb303Q8, py::arg("saw"),
		py::arg("frequency"), py::arg("shape") = 0.0,
		py::arg("sample_rate") = 192000.0);
	module.def("tb303_oscillator_x1", &RenderTb303Oscillator<Tb303OscillatorX1>,
		py::arg("pitch"), py::arg("slide"), py::arg("fm"),
		py::arg("shape"), py::arg("wave"), py::arg("sample_rate") = 48000.0,
		py::arg("slide_time") = 0.060, py::arg("linear_fm") = false,
		py::arg("sync") = py::none());
	module.def("tb303_oscillator_x2", &RenderTb303Oscillator<Tb303OscillatorX2>,
		py::arg("pitch"), py::arg("slide"), py::arg("fm"),
		py::arg("shape"), py::arg("wave"), py::arg("sample_rate") = 48000.0,
		py::arg("slide_time") = 0.060, py::arg("linear_fm") = false,
		py::arg("sync") = py::none());
	module.def("tb303_oscillator_x2_order5",
		&RenderTb303Oscillator<Tb303OscillatorX2Order5>,
		py::arg("pitch"), py::arg("slide"), py::arg("fm"),
		py::arg("shape"), py::arg("wave"), py::arg("sample_rate") = 48000.0,
		py::arg("slide_time") = 0.060, py::arg("linear_fm") = false,
		py::arg("sync") = py::none());
	module.def("tb303_oscillator_x4", &RenderTb303Oscillator<Tb303OscillatorX4>,
		py::arg("pitch"), py::arg("slide"), py::arg("fm"),
		py::arg("shape"), py::arg("wave"), py::arg("sample_rate") = 48000.0,
		py::arg("slide_time") = 0.060, py::arg("linear_fm") = false,
		py::arg("sync") = py::none());
	module.def("tb303_oscillator_x4_order5",
		&RenderTb303Oscillator<Tb303OscillatorX4Order5>,
		py::arg("pitch"), py::arg("slide"), py::arg("fm"),
		py::arg("shape"), py::arg("wave"), py::arg("sample_rate") = 48000.0,
		py::arg("slide_time") = 0.060, py::arg("linear_fm") = false,
		py::arg("sync") = py::none());

	using WavefoldOscillatorX1 =
		tfdsp::WavefoldOscillator<tfdsp::DummyResampler>;
	using WavefoldOscillatorX2 =
		tfdsp::WavefoldOscillator<tfdsp::X2Resampler_Order7>;
	using WavefoldOscillatorX4 =
		tfdsp::WavefoldOscillator<tfdsp::X4Resampler_Order7>;
	using WavefoldOscillatorX16 =
		tfdsp::WavefoldOscillator<tfdsp::X16Resampler_Order7>;
	module.def("wavefolder_transfer", [](py::array_t<double,
		py::array::c_style | py::array::forcecast> input, int character)
	{
		return EvaluateWavefolderFunction(input, character, false);
	}, py::arg("input"), py::arg("character") = 0);
	module.def("wavefolder_primitive", [](py::array_t<double,
		py::array::c_style | py::array::forcecast> input, int character)
	{
		return EvaluateWavefolderFunction(input, character, true);
	}, py::arg("input"), py::arg("character") = 0);
	module.def("wavefolder_adaa", &EvaluateWavefolderAdaa,
		py::arg("input"), py::arg("character") = 0);
	module.def("unison_spread_cents", &tfdsp::UnisonSpreadCents,
		py::arg("control"));
	module.def("unison_pitch_positions", [](int voices)
	{
		const int count = std::clamp(voices, 1, tfdsp::MaximumUnisonVoices);
		const auto positions = tfdsp::UnisonPitchPositions(count);
		py::array_t<double> result(count);
		auto output = result.mutable_unchecked<1>();
		for (int voice = 0; voice < count; ++voice)
			output(voice) = positions[voice];
		return result;
	}, py::arg("voices"));
	module.def("unison_output_gain", &tfdsp::UnisonOutputGain,
		py::arg("voices"));
	module.def("stacked_oscillator", &RenderStackedOscillator,
		py::arg("frequency"), py::arg("pulse_width"), py::arg("voices") = 7,
		py::arg("spread_cents") = 4.0, py::arg("pulse_mix") = 0.0,
		py::arg("width") = 0.65, py::arg("sample_rate") = 48000.0);
	module.def("stacked_oscillator_pitch_positions", [](int voices)
	{
		const int count = std::clamp(voices, 1,
			tfdsp::MaximumStackedOscillatorVoices);
		const auto positions = tfdsp::StackedOscillatorPitchPositions(count);
		py::array_t<double> result(count);
		auto output = result.mutable_unchecked<1>();
		for (int voice = 0; voice < count; ++voice)
			output(voice) = positions[voice];
		return result;
	}, py::arg("voices"));
	module.def("stacked_oscillator_pan_positions", [](int voices)
	{
		const int count = std::clamp(voices, 1,
			tfdsp::MaximumStackedOscillatorVoices);
		const auto positions = tfdsp::StackedOscillatorPanPositions(count);
		py::array_t<double> result(count);
		auto output = result.mutable_unchecked<1>();
		for (int voice = 0; voice < count; ++voice)
			output(voice) = positions[voice];
		return result;
	}, py::arg("voices"));
	module.def("wavefold_oscillator_x1",
		&RenderWavefoldOscillator<WavefoldOscillatorX1>,
		py::arg("frequency"), py::arg("morph"), py::arg("fold"),
		py::arg("symmetry"), py::arg("sample_rate") = 48000.0,
		py::arg("adaa") = false, py::arg("character") = 0);
	module.def("wavefold_oscillator_x2",
		&RenderWavefoldOscillator<WavefoldOscillatorX2>,
		py::arg("frequency"), py::arg("morph"), py::arg("fold"),
		py::arg("symmetry"), py::arg("sample_rate") = 48000.0,
		py::arg("adaa") = false, py::arg("character") = 0);
	module.def("wavefold_oscillator_x4",
		&RenderWavefoldOscillator<WavefoldOscillatorX4>,
		py::arg("frequency"), py::arg("morph"), py::arg("fold"),
		py::arg("symmetry"), py::arg("sample_rate") = 48000.0,
		py::arg("adaa") = false, py::arg("character") = 0);
	module.def("wavefold_oscillator_x16",
		&RenderWavefoldOscillator<WavefoldOscillatorX16>,
		py::arg("frequency"), py::arg("morph"), py::arg("fold"),
		py::arg("symmetry"), py::arg("sample_rate") = 48000.0,
		py::arg("adaa") = false, py::arg("character") = 0);
	module.def("wavefolder_external_x2",
		&RenderWavefolderExternal<WavefoldOscillatorX2>,
		py::arg("audio"), py::arg("fold"), py::arg("symmetry"),
		py::arg("sample_rate") = 48000.0, py::arg("adaa") = false,
		py::arg("character") = 0);
	module.def("wavefolder_external_x4",
		&RenderWavefolderExternal<WavefoldOscillatorX4>,
		py::arg("audio"), py::arg("fold"), py::arg("symmetry"),
		py::arg("sample_rate") = 48000.0, py::arg("adaa") = false,
		py::arg("character") = 0);
	module.def("wavefolder_external_x16",
		&RenderWavefolderExternal<WavefoldOscillatorX16>,
		py::arg("audio"), py::arg("fold"), py::arg("symmetry"),
		py::arg("sample_rate") = 48000.0, py::arg("adaa") = false,
		py::arg("character") = 0);

}
