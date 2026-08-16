#include <algorithm>
#include <cstddef>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "models/VCAcore.hpp"
#include "models/Arp4019Vca.hpp"
#include "models/Arp4072Filter.hpp"
#include "models/DiodeLadderFilter.hpp"
#include "models/OtaVca.hpp"
#include "models/Tb303Voice.hpp"
#include "models/Tb303Oscillator.hpp"
#include "models/VdpOscillator.hpp"
#include "models/VdpSplitOscillator.hpp"
#include "tfdsp/control.hpp"
#include "tfdsp/noise.hpp"
#include "tfdsp/sampleRate.hpp"

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
		double cutoff, double resonance, double driveGain, bool extendedCutoff,
		double sampleRate)
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
				driveGain, extendedCutoff);
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
}

PYBIND11_MODULE(_triggerfish_dsp, module)
{
	module.doc() = "TriggerFish DSP development bindings";

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
	module.def("arp4072_x1", &RenderArp4072<Arp4072X1>, py::arg("audio"),
		py::arg("cutoff"), py::arg("resonance") = 0.0,
		py::arg("drive_gain") = 1.0, py::arg("extended_cutoff") = false,
		py::arg("sample_rate") = 48000.0);
	module.def("arp4072_x2", &RenderArp4072<Arp4072X2>, py::arg("audio"),
		py::arg("cutoff"), py::arg("resonance") = 0.0,
		py::arg("drive_gain") = 1.0, py::arg("extended_cutoff") = false,
		py::arg("sample_rate") = 48000.0);
	module.def("arp4072_x4", &RenderArp4072<Arp4072X4>, py::arg("audio"),
		py::arg("cutoff"), py::arg("resonance") = 0.0,
		py::arg("drive_gain") = 1.0, py::arg("extended_cutoff") = false,
		py::arg("sample_rate") = 48000.0);
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
		values["stage_tanh_scale_per_volt"] =
			Arp4072X4::StageTanhScalePerVolt();
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

}
