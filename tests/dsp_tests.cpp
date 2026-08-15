#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>

#include "models/OTA1PoleIntegrator.hpp"
#include "models/DiodeLadderFilter.hpp"
#include "models/Transistor1PoleIntegrator.hpp"
#include "models/Tb303Voice.hpp"
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

	tfdsp::DiodeLadderFilter<tfdsp::X2Resampler_Order7> diodeFilter(
		tfdsp::CreateX2Resampler_Chebychev7);
	diodeFilter.SetSampleRate(48000.0);
	float diodeOutput = 0.0f;
	for (int i = 0; i < 48000; ++i)
	{
		const double input = 5.0 * std::sin(2.0 * tfdsp::PI * 110.0 * i / 48000.0);
		diodeOutput = diodeFilter.Step(input, 1000.0, 0.5, false, 1.0, 0.0);
	}
	Check(std::isfinite(diodeOutput), "diode ladder produces finite output");
	Check(diodeFilter.SolverFailures() == 0,
		"diode ladder converges at stock settings");
	for (int i = 0; i < 48000; ++i)
	{
		const double input = 5.0 * std::sin(2.0 * tfdsp::PI * 997.0 * i / 48000.0);
		diodeOutput = diodeFilter.Step(input, 18000.0, 1.0, true, 66.6, 1.0);
	}
	Check(std::isfinite(diodeOutput),
		"diode ladder remains finite at maximum drive and resonance");
	Check(diodeFilter.SolverFailures() == 0,
		"diode ladder converges at maximum drive and resonance");
	Check(diodeFilter.Step(std::numeric_limits<double>::infinity(), 1000.0,
		0.5, false, 1.0, 0.0) == 0.0f,
		"diode ladder rejects non-finite input");

	// The two exposed audio paths use independent decimator state. An identity
	// post-processor must nevertheless produce the same signal as the LP path.
	tfdsp::DiodeLadderFilter<tfdsp::X2Resampler_Order7> dualOutputFilter(
		tfdsp::CreateX2Resampler_Chebychev7);
	dualOutputFilter.SetSampleRate(48000.0);
	double maximumDualOutputDifference = 0.0;
	for (int i = 0; i < 48000; ++i)
	{
		const double input = 5.0 * std::sin(2.0 * tfdsp::PI * 997.0 * i / 48000.0);
		const auto outputs = dualOutputFilter.StepWithPostProcessor(input,
			18000.0, 0.4, false, 1.0, 0.0, 1.0,
			[](double audio, double) { return audio; });
		maximumDualOutputDifference = std::max(maximumDualOutputDifference,
			std::abs(static_cast<double>(outputs.lowPass) - outputs.postProcessed));
	}
	Check(maximumDualOutputDifference < 1.0e-7,
		"diode ladder LP and post-processor decimators remain phase aligned");

	tfdsp::OtaVcaCore otaVca;
	constexpr double controlCurrent = 200.0e-6;
	const double expectedGm = 0.85 * controlCurrent / (2.0 * 0.02585);
	Check(std::abs(otaVca.SmallSignalTransconductance(controlCurrent) -
		expectedGm) < 1.0e-15,
		"OTA VCA exposes the calibrated small-signal transconductance");
	const double positiveOta = otaVca.ProcessCurrent(0.005, controlCurrent);
	const double negativeOta = otaVca.ProcessCurrent(-0.005, controlCurrent);
	Check(std::abs(positiveOta + negativeOta) < 1.0e-15,
		"OTA VCA differential-pair transfer is odd symmetric");
	Check(std::abs(otaVca.ProcessCurrent(1.0, controlCurrent)) <=
		0.85 * controlCurrent,
		"OTA VCA large-signal current is bounded by control current");
	Check(otaVca.ProcessCurrent(std::numeric_limits<double>::infinity(),
		controlCurrent) == 0.0,
		"OTA VCA rejects non-finite input");

	tfdsp::Tb303Articulation articulation;
	articulation.SetSampleRate(48000.0);
	double volumeBeforeAttack = 0.0;
	double volumeAfterAttack = 0.0;
	for (int i = 0; i < 480; ++i)
	{
		const auto output = articulation.Step(10.0, 0.0, 0.0,
			0.5, 0.2, 0.5);
		if (i == 143)
			volumeBeforeAttack = output.volumeEnvelope;
		if (i == 383)
			volumeAfterAttack = output.volumeEnvelope;
	}
	Check(volumeBeforeAttack == 0.0,
		"stock volume envelope observes its onset delay");
	Check(volumeAfterAttack > 0.99,
		"stock volume envelope reaches full level after delay and attack");

	// A tied/high gate must not retrigger the MEG. Its value should continue
	// decreasing through a second would-be note event.
	const double tiedBefore = articulation.Step(10.0, 0.0, 0.0,
		0.5, 0.2, 0.5).mainEnvelope;
	for (int i = 0; i < 480; ++i)
		articulation.Step(10.0, 0.0, 0.0, 0.5, 0.2, 0.5);
	const double tiedAfter = articulation.Step(10.0, 0.0, 0.0,
		0.5, 0.2, 0.5).mainEnvelope;
	Check(tiedAfter < tiedBefore,
		"tied gate does not retrigger the main envelope");

	tfdsp::Tb303Articulation accentedArticulation;
	accentedArticulation.SetSampleRate(48000.0);
	double accentAfterOneTimeConstant = 0.0;
	for (int i = 0; i <= 75; ++i)
		accentAfterOneTimeConstant = accentedArticulation.Step(10.0, 10.0,
			0.0, 0.5, 0.2, 0.5).vcaAccent;
	Check(std::abs(accentAfterOneTimeConstant - (1.0 - std::exp(-1.0))) < 0.02,
		"accent-to-VCA branch follows the 47k/33nF time constant");

	tfdsp::Tb303Vca tb303Vca;
	tb303Vca.SetSampleRate(48000.0);
	double vcaPeak = 0.0;
	for (int i = 0; i < 48000; ++i)
	{
		const double input = 0.1 * std::sin(2.0 * tfdsp::PI * 1000.0 * i / 48000.0);
		vcaPeak = std::max(vcaPeak, std::abs(tb303Vca.Step(input, 1.0, 0.0)));
	}
	Check(std::abs(vcaPeak / 0.1 - 1.0) < 0.03,
		"TB-303 VCA wrapper is unity-calibrated at small signal");
	Check(tb303Vca.Step(std::numeric_limits<double>::quiet_NaN(),
		1.0, 0.0) == 0.0,
		"TB-303 VCA wrapper rejects non-finite input");

	if (failures == 0)
		std::cout << "All TriggerFish DSP tests passed\n";
	return failures == 0 ? 0 : 1;
}
