#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>

#include "models/OTA1PoleIntegrator.hpp"
#include "models/Arp4019Vca.hpp"
#include "models/ArpEnvelope.hpp"
#include "models/Arp4072Filter.hpp"
#include "models/DiodeLadderFilter.hpp"
#include "models/Transistor1PoleIntegrator.hpp"
#include "models/Tb303Voice.hpp"
#include "models/Tb303Oscillator.hpp"
#include "models/VCAcore.hpp"
#include "models/VdpSplitOscillator.hpp"
#include "tfdsp/control.hpp"
#include "tfdsp/minblep.hpp"
#include "tfdsp/noise.hpp"
#include "tfdsp/nonlinear.hpp"
#include "tfdsp/oscillator.hpp"
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

	tfdsp::FractionalSchmittTrigger fractionalTrigger;
	Check(!fractionalTrigger.Process(0.0).triggered,
		"fractional trigger initializes without a false edge");
	const auto halfSampleEdge = fractionalTrigger.Process(2.0);
	Check(halfSampleEdge.triggered &&
		std::abs(halfSampleEdge.position - 0.5) < 1.0e-12,
		"fractional trigger interpolates a rising threshold crossing");
	Check(!fractionalTrigger.Process(0.5).triggered && fractionalTrigger.IsHigh(),
		"fractional trigger applies hysteresis while high");
	Check(!fractionalTrigger.Process(0.0).triggered && !fractionalTrigger.IsHigh(),
		"fractional trigger releases at its low threshold");

	using TestMinBlep = tfdsp::MinBlepGenerator<8, 32, double>;
	const auto& minBlepKernel = TestMinBlep::Kernel();
	Check(std::all_of(minBlepKernel.begin(), minBlepKernel.end(),
		[](double value) { return std::isfinite(value); }),
		"minBLEP kernel contains only finite values");
	Check(std::abs(minBlepKernel.back() - 1.0) < 1.0e-15,
		"minBLEP kernel reaches its normalized final value");
	double earlyImpulseEnergy = 0.0;
	double lateImpulseEnergy = 0.0;
	for (std::size_t i = 0; i + 1 < minBlepKernel.size(); ++i)
	{
		const double impulse = minBlepKernel[i + 1] - minBlepKernel[i];
		if (i < minBlepKernel.size() / 2)
			earlyImpulseEnergy += impulse * impulse;
		else
			lateImpulseEnergy += impulse * impulse;
	}
	Check(earlyImpulseEnergy > 100.0 * lateImpulseEnergy,
		"minBLEP reconstruction concentrates energy at minimum delay");
	TestMinBlep unitMinBlep;
	TestMinBlep doubleMinBlep;
	unitMinBlep.InsertDiscontinuity(-0.375, 1.0);
	doubleMinBlep.InsertDiscontinuity(-0.375, 2.0);
	double maximumLinearityError = 0.0;
	for (int i = 0; i < TestMinBlep::CorrectionSamples; ++i)
		maximumLinearityError = std::max(maximumLinearityError,
			std::abs(doubleMinBlep.Process() - 2.0 * unitMinBlep.Process()));
	Check(maximumLinearityError < 1.0e-14,
		"minBLEP correction is linear in discontinuity magnitude");
	Check(std::abs(unitMinBlep.Process()) < 1.0e-15,
		"minBLEP correction ends after its finite support");

	const auto startEvent = tfdsp::MapEventToOversampledFrame<4>(0.0);
	const auto middleEvent = tfdsp::MapEventToOversampledFrame<4>(0.625);
	const auto endEvent = tfdsp::MapEventToOversampledFrame<4>(1.0);
	Check(startEvent.segment == 0 && startEvent.position == 0.0 &&
		middleEvent.segment == 2 && middleEvent.position == 0.5 &&
		endEvent.segment == 3 && endEvent.position == 1.0,
		"host-frame events map to the correct oversampled segment");
	tfdsp::BandlimitedSawOscillator<> genericSaw;
	genericSaw.Step(0.2);
	genericSaw.Step(0.1, 0.25);
	Check(std::abs(genericSaw.Phase() - 0.075) < 1.0e-12,
		"generic hard sync resets phase at the fractional crossing");
	genericSaw.Reset();
	bool reverseSyncFinite = true;
	for (int i = 0; i < 1000; ++i)
		reverseSyncFinite = reverseSyncFinite && std::isfinite(genericSaw.Step(
			-0.013, i % 97 == 0 ? 0.37 : -1.0));
	Check(reverseSyncFinite,
		"generic hard sync remains finite while running backwards");

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

	using Arp4072 = tfdsp::Arp4072Filter<tfdsp::X4Resampler_Order7>;
	Check(std::abs(Arp4072::FeedbackBaseScale() /
		Arp4072::AudioBaseScale() - 6.583) < 0.01,
		"ARP 4072 resonance return drives its limiter much harder than audio");
	Check(Arp4072::AudioBaseVolts(5.0) < 0.5 * Arp4072::ThermalVoltage,
		"ARP 4072 nominal audio remains below half a thermal voltage at its limiter");
	Check(Arp4072::FeedbackBaseVolts(5.0) >
		2.5 * Arp4072::ThermalVoltage,
		"ARP 4072 full resonance return spans several thermal voltages");
	Check(std::abs(Arp4072::SmallSignalInputGain() - 1.0) < 1.0e-12,
		"ARP 4072 gain-trim calibration produces unity small-signal gain");
	Check(Arp4072::SmallSignalFeedbackGain() > 6.0 &&
		Arp4072::SmallSignalFeedbackGain() < 7.0,
		"ARP 4072 feedback loop gain follows the circuit level shifting");
	Check(Arp4072::LimiterEquivalentPeakVolts() > 3.0 &&
		Arp4072::LimiterEquivalentPeakVolts() < 3.2,
		"ARP 4072 limiter tail current sets the expected first-stage drive");
	Check(Arp4072::StageBaseResistanceOhms() > 212.2 &&
		Arp4072::StageBaseResistanceOhms() < 212.4,
		"ARP 4072 stage model includes signal-resistor loading at each base");
	Arp4072 arp4072(tfdsp::CreateX4Resampler_Cheby7);
	arp4072.SetSampleRate(48000.0);
	double arpPeak = 0.0;
	bool arpFinite = true;
	for (int i = 0; i < 48000; ++i)
	{
		const double input = 0.1 * std::sin(2.0 * 3.14159265358979323846 *
			100.0 * i / 48000.0);
		const double output = arp4072.Step(input, 5000.0, 0.0);
		arpFinite = arpFinite && std::isfinite(output);
		if (i >= 24000)
			arpPeak = std::max(arpPeak, std::abs(output));
	}
	Check(arpFinite, "ARP 4072 model remains finite");
	Check(std::abs(arpPeak - 0.1) < 0.004,
		"ARP 4072 open-filter small signal level is approximately unity");
	Check(arp4072.SolverFailures() == 0,
		"ARP 4072 nominal signal converges without solver failures");
	Check(arp4072.Step(std::numeric_limits<double>::quiet_NaN(),
		1000.0, 0.0) == 0.0f,
		"ARP 4072 rejects non-finite input");
	arp4072.SetSampleRate(48000.0);
	bool arpStressFinite = true;
	for (int i = 0; i < 48000; ++i)
	{
		const double input = (i & 16) ? 10.0 : -10.0;
		arpStressFinite = arpStressFinite && std::isfinite(arp4072.Step(input,
			20000.0, 1.0, 15.848931924611133, true));
	}
	Check(arpStressFinite,
		"ARP 4072 remains finite at maximum drive, cutoff, and resonance");
	Check(arp4072.SolverFailures() == 0,
		"ARP 4072 extreme 4x stress converges without solver failures");

	using Arp4072NoNewton = tfdsp::Arp4072Filter<tfdsp::DummyResampler, 0>;
	Arp4072NoNewton arpFallback(tfdsp::CreateDummyResampler);
	arpFallback.SetSampleRate(48000.0);
	const auto fallbackState = arpFallback.State();
	const double fallbackOutput = arpFallback.Step(5.0, 1000.0, 1.0);
	Check(fallbackOutput == 0.0 && arpFallback.State() == fallbackState &&
		arpFallback.SolverFailures() == 1,
		"ARP 4072 solver failure holds the previous valid state");

	using Arp4072X2 = tfdsp::Arp4072Filter<tfdsp::X2Resampler_Order7>;
	Arp4072X2 arpPostSafety(tfdsp::CreateX2Resampler_Chebychev7);
	arpPostSafety.SetSampleRate(48000.0);
	double arpPostSafetyPeak = 0.0;
	for (int i = 0; i < 4096; ++i)
	{
		const auto rendered = arpPostSafety.StepWithPostProcessor(0.0, 1000.0,
			0.0, 1.0, false, 0.0, 0.0,
			[](double, double, double) { return 1000.0; });
		arpPostSafetyPeak = std::max(arpPostSafetyPeak,
			std::abs(static_cast<double>(rendered.postProcessed)));
	}
	Check(arpPostSafetyPeak <= 14.500001,
		"ARP integrated post-processor output retains decimator safety");

	tfdsp::ArpEnvelope arpEnvelope;
	arpEnvelope.SetSampleRate(48000.0);
	Check(arpEnvelope.GetStage() == tfdsp::ArpEnvelope::Stage::Idle,
		"ARP envelope exposes its idle stage for panel indication");
	double envelope = arpEnvelope.Step(10.0, 0.0, 0.1, 0.2, 0.4, 0.3);
	Check(arpEnvelope.GetStage() == tfdsp::ArpEnvelope::Stage::Attack,
		"ARP envelope exposes its active attack stage");
	for (int i = 1; i < 1200; ++i)
		envelope = arpEnvelope.Step(10.0, 0.0, 0.1, 0.2, 0.4, 0.3);
	const double expectedQuarterAttack = tfdsp::ArpEnvelope::NormalizedCurve(
		0.25, tfdsp::ArpEnvelope::AttackCurve(0.0));
	Check(std::abs(envelope - expectedQuarterAttack) < 1.0e-9 &&
		envelope > 0.36 && envelope < 0.361,
		"ARP attack follows the 15 V charge target at its default curve");
	for (int i = 1200; i < 4800; ++i)
		envelope = arpEnvelope.Step(10.0, 0.0, 0.1, 0.2, 0.4, 0.3);
	Check(envelope > 0.998 && envelope <= 1.0,
		"ARP ADSR reaches its attack peak at the calibrated time");
	for (int i = 0; i < 9600; ++i)
		envelope = arpEnvelope.Step(10.0, 0.0, 0.1, 0.2, 0.4, 0.3);
	Check(std::abs(envelope - 0.4) < 0.002,
		"ARP ADSR reaches sustain through an RC decay");
	Check(arpEnvelope.GetStage() == tfdsp::ArpEnvelope::Stage::Sustain,
		"ARP envelope exposes its active sustain stage");
	envelope = arpEnvelope.Step(0.0, 0.0, 0.1, 0.2, 0.4, 0.3);
	Check(arpEnvelope.GetStage() == tfdsp::ArpEnvelope::Stage::Release,
		"ARP envelope exposes its active release stage");
	for (int i = 1; i < 14400; ++i)
		envelope = arpEnvelope.Step(0.0, 0.0, 0.1, 0.2, 0.4, 0.3);
	Check(envelope < 0.001,
		"ARP ADSR reaches zero through an RC release");
	Check(arpEnvelope.GetStage() == tfdsp::ArpEnvelope::Stage::Idle,
		"ARP envelope returns its stage indicator to idle");
	arpEnvelope.Reset();
	arpEnvelope.SetMode(tfdsp::ArpEnvelope::Mode::Ar);
	for (int i = 0; i < 4800; ++i)
		envelope = arpEnvelope.Step(10.0, 0.0, 0.1, 0.01, 0.0, 0.1);
	for (int i = 0; i < 4800; ++i)
		envelope = arpEnvelope.Step(10.0, 10.0, 0.1, 0.01, 0.0, 0.1);
	Check(envelope > 0.999,
		"ARP mode holds its peak and ignores trigger retriggering");
	Check(arpEnvelope.GetStage() == tfdsp::ArpEnvelope::Stage::Hold,
		"ARP mode exposes its held attack stage");
	arpEnvelope.SetMode(tfdsp::ArpEnvelope::Mode::Adsr);
	for (int i = 0; i < 480; ++i)
		envelope = arpEnvelope.Step(10.0, 0.0, 0.1, 0.01, 0.25, 0.1);
	Check(envelope < 0.5,
		"switching a held AR envelope to ADSR enters decay");
	arpEnvelope.SetMode(tfdsp::ArpEnvelope::Mode::Ar);
	for (int i = 0; i < 4800; ++i)
		envelope = arpEnvelope.Step(10.0, 0.0, 0.1, 0.01, 0.25, 0.1);
	Check(envelope > 0.999,
		"switching a held ADSR envelope to AR returns to its peak smoothly");
	std::array<double, 3> attackQuarterValues{};
	for (int curveIndex = 0; curveIndex < 3; ++curveIndex)
	{
		arpEnvelope.Reset();
		const double curve = static_cast<double>(curveIndex - 1);
		for (int i = 0; i < 1200; ++i)
			attackQuarterValues[curveIndex] = arpEnvelope.Step(10.0, 0.0,
				0.1, 0.2, 0.4, 0.3, curve);
	}
	Check(attackQuarterValues[0] < attackQuarterValues[1] &&
		attackQuarterValues[1] < attackQuarterValues[2],
		"ARP envelope curve spans near-linear through tighter exponentials");
	arpEnvelope.Reset();
	for (int i = 0; i < 240; ++i)
		envelope = arpEnvelope.Step(10.0, 0.0, 0.005, 0.2, 0.4, 0.3);
	Check(envelope > 0.999,
		"ARP-inspired attack remains usable at a five millisecond setting");

	using Arp4019 = tfdsp::Arp4019Vca<tfdsp::X4Resampler_Order7>;
	Check(std::abs(Arp4019::AudioInputScale() - 0.0021734835) < 1.0e-10,
		"ARP 4019 audio attenuator includes the external and internal series resistors");
	Check(std::abs(Arp4019::SmallSignalGainAtUnityControl() - 1.0) < 1.0e-12,
		"ARP 4019 unity-current calibration follows the circuit ratios");
	Check(Arp4019::UnityControlCurrentAmps() > 0.0004 &&
		Arp4019::UnityControlCurrentAmps() < 0.00045,
		"ARP 4019 unity control current is in the component-derived range");
	Check(Arp4019::OutputBandwidthHz > 28000.0 &&
		Arp4019::OutputBandwidthHz < 29000.0,
		"ARP 4019 feedback capacitor retains the original HF rolloff");
	Arp4019 arp4019(tfdsp::CreateX4Resampler_Cheby7);
	arp4019.SetSampleRate(48000.0);
	double arpVcaLinearPeak = 0.0;
	double arpVcaExponentialPeak = 0.0;
	for (int i = 0; i < 48000; ++i)
	{
		const double input = 0.1 * std::sin(2.0 *
			3.14159265358979323846 * 1000.0 * i / 48000.0);
		const double output = arp4019.Step(input, 0.0, 10.0, -10.0);
		if (i >= 24000)
			arpVcaLinearPeak = std::max(arpVcaLinearPeak, std::abs(output));
	}
	arp4019.Reset();
	for (int i = 0; i < 48000; ++i)
	{
		const double input = 0.1 * std::sin(2.0 *
			3.14159265358979323846 * 1000.0 * i / 48000.0);
		const double output = arp4019.Step(input, 0.0, -1.0, 9.0);
		if (i >= 24000)
			arpVcaExponentialPeak = std::max(arpVcaExponentialPeak,
				std::abs(output));
	}
	Check(std::abs(arpVcaLinearPeak - 0.1) < 0.004,
		"ARP 4019 linear 10 V control gives approximately unity gain");
	Check(std::abs(arpVcaExponentialPeak / arpVcaLinearPeak -
		0.316227766) < 0.01,
		"ARP 4019 exponential control follows the specified 10 dB/V law");
	Check(arp4019.Step(std::numeric_limits<double>::infinity(), 0.0,
		10.0, 0.0) == 0.0f,
		"ARP 4019 rejects non-finite input");

	Arp4072 arpVoiceFilter(tfdsp::CreateX4Resampler_Cheby7);
	Arp4019 arpVoiceVca(tfdsp::CreateX4Resampler_Cheby7);
	arpVoiceFilter.SetSampleRate(48000.0);
	arpVoiceVca.SetSampleRate(48000.0);
	double arpVoiceFilterPeak = 0.0;
	double arpVoiceVcaPeak = 0.0;
	bool arpVoiceFinite = true;
	for (int i = 0; i < 48000; ++i)
	{
		const double input = 0.01 * std::sin(2.0 *
			3.14159265358979323846 * 200.0 * i / 48000.0);
		const auto rendered = arpVoiceFilter.StepWithPostProcessor(input,
			8000.0, 0.0, 1.0, false, 10.0, -10.0,
			[&](double filtered, double linearCv, double exponentialCv)
			{
				return arpVoiceVca.ProcessOversampled(filtered, linearCv,
					exponentialCv);
			});
		arpVoiceFinite = arpVoiceFinite && std::isfinite(rendered.lowPass) &&
			std::isfinite(rendered.postProcessed);
		if (i >= 24000)
		{
			arpVoiceFilterPeak = std::max(arpVoiceFilterPeak,
				std::abs(static_cast<double>(rendered.lowPass)));
			arpVoiceVcaPeak = std::max(arpVoiceVcaPeak,
				std::abs(static_cast<double>(rendered.postProcessed)));
		}
	}
	Check(arpVoiceFinite && arpVoiceVcaPeak > 0.005,
		"ARP 4072 and 4019 share one finite oversampled signal path");
	Check(std::abs(arpVoiceVcaPeak / arpVoiceFilterPeak - 1.0) < 0.04,
		"ARP voice-core VCA preserves the filter level at unity control");

	Arp4072X2 arpVoiceFilterX2(tfdsp::CreateX2Resampler_Chebychev7);
	tfdsp::Arp4019Vca<tfdsp::X2Resampler_Order7> arpVoiceVcaX2(
		tfdsp::CreateX2Resampler_Chebychev7);
	arpVoiceFilterX2.SetSampleRate(48000.0);
	arpVoiceVcaX2.SetSampleRate(48000.0);
	double arpVoiceFilterPeakX2 = 0.0;
	double arpVoiceVcaPeakX2 = 0.0;
	bool arpVoiceFiniteX2 = true;
	for (int i = 0; i < 48000; ++i)
	{
		const double input = 0.01 * std::sin(2.0 *
			3.14159265358979323846 * 200.0 * i / 48000.0);
		const auto rendered = arpVoiceFilterX2.StepWithPostProcessor(input,
			8000.0, 0.0, 1.0, false, 10.0, -10.0,
			[&](double filtered, double linearCv, double exponentialCv)
			{
				return arpVoiceVcaX2.ProcessOversampled(filtered, linearCv,
					exponentialCv);
			});
		arpVoiceFiniteX2 = arpVoiceFiniteX2 &&
			std::isfinite(rendered.lowPass) &&
			std::isfinite(rendered.postProcessed);
		if (i >= 24000)
		{
			arpVoiceFilterPeakX2 = std::max(arpVoiceFilterPeakX2,
				std::abs(static_cast<double>(rendered.lowPass)));
			arpVoiceVcaPeakX2 = std::max(arpVoiceVcaPeakX2,
				std::abs(static_cast<double>(rendered.postProcessed)));
		}
	}
	Check(arpVoiceFiniteX2 && arpVoiceVcaPeakX2 > 0.005,
		"ARP 4072 and 4019 share one finite 2x signal path");
	Check(std::abs(arpVoiceVcaPeakX2 / arpVoiceFilterPeakX2 - 1.0) < 0.04,
		"ARP 2x voice-core VCA preserves the filter level at unity control");

	VCA_TransistorCore<tfdsp::X2Resampler_Order7> vca(tfdsp::CreateX2Resampler_Chebychev7);
	vca.SetSampleRate(48000.0f);
	Check(std::isfinite(vca.Step(1.0f, 0.5f, 1.0f)), "VCA model produces finite output");
	for (int i = 0; i < 128; ++i)
		vca.StepControls(1.0f, 0.25f, 0.5f, 20.0f, 1.0f);
	Check(std::isfinite(vca.LastControl()) && vca.LastControl() > 0.25f &&
		vca.LastControl() < 1.0f,
		"VCA reconstructs linear and exponential controls before shaping");
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
		oscillatorOutput = oscillator.StepLogAngularFrequency(0.0, 0.5,
			std::log2(2.0 * tfdsp::PI * 261.625565));
	Check(std::isfinite(oscillatorOutput),
		"VDPO reconstructs pitch before conversion to angular frequency");
	for (int i = 0; i < 100; ++i)
		oscillatorOutput = oscillator.Step(0.0, 9.0, 2.0 * tfdsp::PI * 18000.0);
	Check(std::isfinite(oscillatorOutput), "VDPO model remains stable above the legacy frequency limit");
	Check(oscillator.Step(0.0, 0.5, std::numeric_limits<double>::quiet_NaN()) == 0.0f,
		"VDPO model rejects non-finite input");
	oscillator.Reset();

	tfdsp::Tb303SquareShaper squareShaper;
	squareShaper.SetSampleRate(192000.0);
	double squarePeak = 0.0;
	for (int i = 0; i < 192000; ++i)
	{
		const double phase = std::fmod(85.0 * i / 192000.0, 1.0);
		squarePeak = std::max(squarePeak, std::abs(squareShaper.Step(
			2.0 * phase - 1.0, 85.0, 0.0)));
	}
	Check(squarePeak > 0.8 && squarePeak < 2.0,
		"TB-303 square shaper produces a bounded asymmetric waveform");
	Check(squareShaper.Step(std::numeric_limits<double>::quiet_NaN(),
		85.0, 0.0) == 0.0,
		"TB-303 square shaper rejects non-finite input");

	tfdsp::Tb303Oscillator<tfdsp::X4Resampler_Order7> tb303Oscillator(
		tfdsp::CreateX4Resampler_Cheby7);
	tb303Oscillator.SetSampleRate(48000.0);
	tfdsp::Tb303Oscillator<tfdsp::X4Resampler_Order7>::Output tb303Output;
	for (int i = 0; i < 48000; ++i)
		tb303Output = tb303Oscillator.Step(0.0, false, 0.060, 0.0,
			0.0, false, 0.0, 0.0);
	Check(std::isfinite(tb303Output.saw) &&
		std::isfinite(tb303Output.square) &&
		std::isfinite(tb303Output.mixed),
		"TB-303 oscillator produces finite oversampled outputs");
	Check(std::abs(tb303Output.mixed - tb303Output.saw) < 1.0e-6,
		"TB-303 oscillator saw endpoint matches its mixed output");
	for (int i = 0; i < 2880; ++i)
		tb303Output = tb303Oscillator.Step(1.0, true, 0.060, 0.0,
			0.0, false, 0.0, 1.0);
	Check(std::abs(tb303Output.pitch -
		(1.0 - std::exp(-0.060 / 0.022))) < 2.0e-3,
		"TB-303 stock slide follows the 22 ms pitch-CV time constant");
	const double pitchBeforeStep = tb303Output.pitch;
	tb303Output = tb303Oscillator.Step(2.0, false, 0.060, 0.0,
		0.0, false, 0.0, 1.0);
	Check(tb303Output.pitch > pitchBeforeStep && tb303Output.pitch < 2.01f,
		"TB-303 pitch steps use the oversampled reconstruction filter");
	for (int i = 0; i < 128; ++i)
		tb303Output = tb303Oscillator.Step(2.0, false, 0.060, 0.0,
			0.0, false, 0.0, 1.0);
	Check(std::abs(tb303Output.pitch - 2.0f) < 1.0e-5f,
		"TB-303 reconstructed pitch settles to the requested value");
	const auto invalidTb303Output = tb303Oscillator.Step(
		std::numeric_limits<double>::infinity(), false, 0.060, 0.0,
		0.0, false, 0.0, 0.0);
	Check(invalidTb303Output.mixed == 0.0f,
		"TB-303 oscillator rejects non-finite controls");

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
	Check(tfdsp::AnalogOutputStage::Process(5.0) == 5.0,
		"analog output stage is exactly linear at nominal Rack level");
	Check(tfdsp::AnalogOutputStage::Process(100.0) <=
		tfdsp::AnalogOutputStage::RailVolts,
		"analog output stage approaches its rail without crossing it");
	Check(tfdsp::AnalogOutputStage::ProcessSafety(100.0) <
		tfdsp::AnalogOutputStage::SafetyLimitVolts,
		"analog output safety stage approaches its limit without hard clipping");
	tb303Vca.Reset();
	double overloadedPeak = 0.0;
	for (int i = 0; i < 48000; ++i)
	{
		const double input = 100.0 * std::sin(
			2.0 * tfdsp::PI * 1000.0 * i / 48000.0);
		const double output = tb303Vca.Step(input, 1.0, 1.0);
		if (i >= 4800)
			overloadedPeak = std::max(overloadedPeak, std::abs(output));
	}
	Check(overloadedPeak > tfdsp::AnalogOutputStage::KneeVolts &&
		overloadedPeak < tfdsp::AnalogOutputStage::SafetyLimitVolts,
		"TB-303 VCA overload bends smoothly inside the safety guard");
	Check(tb303Vca.Step(std::numeric_limits<double>::quiet_NaN(),
		1.0, 0.0) == 0.0,
		"TB-303 VCA wrapper rejects non-finite input");

	if (failures == 0)
		std::cout << "All TriggerFish DSP tests passed\n";
	return failures == 0 ? 0 : 1;
}
