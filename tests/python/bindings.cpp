#include <algorithm>
#include <cstddef>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "models/VCAcore.hpp"
#include "models/VdpOscillator.hpp"
#include "models/VdpSplitOscillator.hpp"
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
}

PYBIND11_MODULE(_triggerfish_dsp, module)
{
	module.doc() = "TriggerFish DSP development bindings";

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

}
