#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>

#include "models/OTA1PoleIntegrator.hpp"
#include "models/Transistor1PoleIntegrator.hpp"
#include "models/VCAcore.hpp"
#include "models/VdpSplitOscillator.hpp"
#include "tfdsp/control.hpp"
#include "tfdsp/noise.hpp"
#include "tfdsp/nonlinear.hpp"
#include "tfdsp/sampleRate.hpp"

namespace
{
	int failures = 0;

	void Check(bool condition, const char* message)
	{
		if (!condition)
		{
			std::cerr << "FAIL: " << message << '\n';
			++failures;
		}
	}

	struct CountingGenerator
	{
		using result_type = std::minstd_rand::result_type;
		std::minstd_rand generator{12345};
		std::uint64_t calls{};

		static constexpr result_type min() { return std::minstd_rand::min(); }
		static constexpr result_type max() { return std::minstd_rand::max(); }
		result_type operator()()
		{
			++calls;
			return generator();
		}
	};
}

int main()
{
	using DiscreteGradient2::Tanh;
	Check(std::isfinite(Tanh<double, 1>::Value(1000.0, 999.0)),
		"stable log(cosh) discrete gradient stays finite");
	Check(std::abs(Tanh<double, 1>::Value(1000.0, -1000.0)) < 1e-12,
		"symmetric large discrete gradient is zero");

	Transistor1PoleIntegrator transistor;
	OTA1PoleIntegrator ota;
	Check(std::isfinite(transistor.Step(1.0, 0.5)), "transistor solver produces finite output");
	Check(std::isfinite(ota.Step(1.0, 0.5)), "OTA solver produces finite output");
	Check(transistor.Step(std::numeric_limits<double>::infinity(), 0.5) == 0.0,
		"transistor solver rejects non-finite input");
	Check(ota.Step(std::numeric_limits<double>::quiet_NaN(), 0.5) == 0.0,
		"OTA solver rejects non-finite input");

	tfdsp::RecursiveSineOscillator sine;
	sine.SetFrequency(60.0, 48000.0);
	double peak = 0.0;
	for (int i = 0; i < 48000; ++i)
		peak = std::max(peak, std::abs(sine.Step()));
	Check(std::abs(peak - 1.0) < 1e-10, "recursive sine oscillator holds amplitude");

	tfdsp::InterpolatedOrnsteinUhlenbeck ou;
	ou.Configure(48000.0, 60.0, 2.0, 100.0);
	CountingGenerator rng;
	double drift = 0.0;
	for (int i = 0; i < 48000; ++i)
		drift = ou.Step(rng);
	Check(std::isfinite(drift), "OU drift remains finite");
	Check(rng.calls < 1000, "OU noise generation runs at control rate");
	double maxExp2RelativeError = 0.0;
	for (int i = 0; i <= 20000; ++i)
	{
		const float exponent = -100.0f + 200.0f * i / 20000.0f;
		const double exact = std::exp2(static_cast<double>(exponent));
		maxExp2RelativeError = std::max(maxExp2RelativeError,
			std::abs(static_cast<double>(tfdsp::Exp2Taylor5(exponent)) / exact - 1.0));
	}
	Check(maxExp2RelativeError < 6.0e-6,
		"fast exp2 stays within its relative-error budget");

	Check(tfdsp::detune::linear(5.0, 0.0) == 5.0,
		"linear detune preserves pitch exactly when drift is zero");
	const double detunePitch = -4.75;
	const double detuneHz = 2.0;
	const double detuneReference = detunePitch + std::log1p(
		detuneHz / (261.63 * std::exp2(detunePitch))) / std::log(2.0);
	Check(std::abs(tfdsp::detune::linear(detunePitch, detuneHz) - detuneReference) * 1200.0 < 0.002,
		"linear detune stays within its pitch-error budget");
	Check(tfdsp::detune::linear(std::numeric_limits<double>::quiet_NaN(), 1.0) == 0.0,
		"linear detune rejects non-finite input");

	auto resampler = tfdsp::CreateX2Resampler_Chebychev7();
	resampler->Downsample(resampler->Upsample(1.0));
	resampler->Reset();
	Check(std::abs(resampler->Downsample(resampler->Upsample(0.0))) < 1e-15,
		"resampler reset clears filter history");

	VCA_TransistorCore<tfdsp::X2Resampler_Order7> vca(tfdsp::CreateX2Resampler_Chebychev7);
	vca.SetSampleRate(48000.0f);
	Check(std::isfinite(vca.Step(1.0f, 0.5f, 1.0f)), "VCA model produces finite output");
	Check(vca.Step(std::numeric_limits<float>::infinity(), 0.5f, 1.0f) == 0.0f,
		"VCA model rejects non-finite input");
	vca.Reset();

	VdpSplitOscillator<tfdsp::X4Resampler_Order7> oscillator(tfdsp::CreateX4Resampler_Cheby7);
	oscillator.SetSampleRate(48000.0);
	float oscillatorOutput = 0.0f;
	for (int i = 0; i < 100; ++i)
		oscillatorOutput = oscillator.Step(0.0, 0.5, 2.0 * tfdsp::PI * 261.625565);
	Check(std::isfinite(oscillatorOutput), "VDPO model produces finite output");
	for (int i = 0; i < 100; ++i)
		oscillatorOutput = oscillator.Step(0.0, 9.0, 2.0 * tfdsp::PI * 18000.0);
	Check(std::isfinite(oscillatorOutput), "VDPO model remains stable above the legacy frequency limit");
	Check(oscillator.Step(0.0, 0.5, std::numeric_limits<double>::quiet_NaN()) == 0.0f,
		"VDPO model rejects non-finite input");
	oscillator.Reset();

	if (failures == 0)
		std::cout << "All TriggerFish DSP tests passed\n";
	return failures == 0 ? 0 : 1;
}
