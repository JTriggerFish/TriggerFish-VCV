#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "models/VCAcore.hpp"
#include "models/VdpOscillator.hpp"
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

	module.def("vdpo", [](py::array_t<double, py::array::c_style | py::array::forcecast> audio,
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
}
