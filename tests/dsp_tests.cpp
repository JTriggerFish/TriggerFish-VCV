#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "models/OTA1PoleIntegrator.hpp"
#include "models/Arp4019Vca.hpp"
#include "models/ArpEnvelope.hpp"
#include "models/Arp4072Filter.hpp"
#include "models/DiodeLadderFilter.hpp"
#include "models/ElectricPiano.hpp"
#include "tfdsp/rail.hpp"
#include "models/Transistor1PoleIntegrator.hpp"
#include "models/Tb303Voice.hpp"
#include "models/Tb303Oscillator.hpp"
#include "models/VCAcore.hpp"
#include "models/VdpSplitOscillator.hpp"
#include "tfdsp/control.hpp"
#include "tfdsp/event_merge.hpp"
#include "tfdsp/minblep.hpp"
#include "tfdsp/noise.hpp"
#include "tfdsp/nonlinear.hpp"
#include "tfdsp/oscillator.hpp"
#include "tfdsp/sampleRate.hpp"
#include "tfdsp/unison.hpp"
#include "tfdsp/unison_oscillator.hpp"
#include "tfdsp/wavefolder.hpp"

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

	double HarmonicMagnitude(const std::vector<double>& signal,
		double normalizedFrequency)
	{
		constexpr double TwoPi = 6.283185307179586476925286766559;
		double real = 0.0;
		double imaginary = 0.0;
		for (std::size_t index = 0; index < signal.size(); ++index)
		{
			const double angle = TwoPi * normalizedFrequency *
				static_cast<double>(index);
			real += signal[index] * std::cos(angle);
			imaginary -= signal[index] * std::sin(angle);
		}
		return std::hypot(real, imaginary) /
			static_cast<double>(signal.size());
	}

	double WindowedHarmonicMagnitude(const std::vector<double>& signal,
		double normalizedFrequency, std::size_t begin, std::size_t end)
	{
		constexpr double TwoPi = 6.283185307179586476925286766559;
		begin = std::min(begin, signal.size());
		end = std::clamp(end, begin, signal.size());
		if (begin == end)
			return 0.0;
		double real = 0.0;
		double imaginary = 0.0;
		for (std::size_t index = begin; index < end; ++index)
		{
			const double angle = TwoPi * normalizedFrequency *
				static_cast<double>(index - begin);
			real += signal[index] * std::cos(angle);
			imaginary -= signal[index] * std::sin(angle);
		}
		return std::hypot(real, imaginary) /
			static_cast<double>(end - begin);
	}

	double LevelMatchedResidual(const std::vector<double>& first,
		const std::vector<double>& second, std::size_t begin, std::size_t end)
	{
		begin = std::min({begin, first.size(), second.size()});
		end = std::min({end, first.size(), second.size()});
		double firstEnergy = 0.0;
		double secondEnergy = 0.0;
		double cross = 0.0;
		for (std::size_t sample = begin; sample < end; ++sample)
		{
			firstEnergy += first[sample] * first[sample];
			secondEnergy += second[sample] * second[sample];
			cross += first[sample] * second[sample];
		}
		const double gain = cross / std::max(1.0e-20, secondEnergy);
		double residualEnergy = 0.0;
		for (std::size_t sample = begin; sample < end; ++sample)
		{
			const double residual = first[sample] - gain * second[sample];
			residualEnergy += residual * residual;
		}
		return std::sqrt(residualEnergy /
			std::max(1.0e-20, firstEnergy));
	}
}

int main()
{
	{
		tfdsp::EventBundle sequenced;
		sequenced.voiceCount = 4;
		tfdsp::EventBundle live;
		live.voiceCount = 13;
		for (std::size_t signal = 0; signal < tfdsp::EventSignalCount; ++signal)
		{
			for (std::size_t voice = 0; voice < sequenced.voiceCount; ++voice)
				sequenced.signals[signal][voice] =
					static_cast<float>(100 * signal + voice);
			for (std::size_t voice = 0; voice < live.voiceCount; ++voice)
				live.signals[signal][voice] =
					static_cast<float>(1000 + 100 * signal + voice);
		}
		live.signals[0][1] = std::numeric_limits<float>::quiet_NaN();
		const auto merged = tfdsp::MergeEventBundles(sequenced, live);
		Check(merged.voiceCount == tfdsp::EventVoiceLimit,
			"event merger clips two bundles at Rack's 16-channel limit");
		Check(merged.signals[3][3] == 303.f && merged.signals[3][4] == 1300.f,
			"event merger preserves signal alignment at the bundle boundary");
		Check(merged.signals[0][5] == 0.f,
			"event merger sanitizes non-finite voltages");
	}


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
	TestMinBlep::PrepareKernel();
	const auto& minBlepKernel = TestMinBlep::Kernel();
	const auto* preparedKernel = &minBlepKernel;
	TestMinBlep::PrepareKernel();
	Check(&TestMinBlep::Kernel() == preparedKernel,
		"minBLEP kernel preparation is shared and idempotent");
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

	tfdsp::BandlimitedPulseOscillator<> genericPulse;
	genericPulse.Step(0.2, 0.5);
	genericPulse.Step(0.1, 0.5, 0.25);
	Check(std::abs(genericPulse.Phase() - 0.075) < 1.0e-12,
		"generic pulse hard sync resets phase at the fractional crossing");
	genericPulse.Reset();
	constexpr int pulseSampleRate = 48000;
	constexpr double pulseFrequency = 100.0;
	constexpr double pulseDuty = 0.25;
	for (int i = 0; i < 256; ++i)
		genericPulse.Step(pulseFrequency / pulseSampleRate, pulseDuty);
	double pulseMean = 0.0;
	for (int i = 0; i < pulseSampleRate; ++i)
		pulseMean += genericPulse.Step(pulseFrequency / pulseSampleRate,
			pulseDuty);
	pulseMean /= pulseSampleRate;
	Check(std::abs(pulseMean - (2.0 * pulseDuty - 1.0)) < 2.0e-4,
		"generic pulse preserves continuous-time duty-cycle mean");

	genericPulse.Reset();
	bool reversePulseFinite = true;
	for (int i = 0; i < 10000; ++i)
	{
		const double duty = 0.5 + 0.45 * std::sin(
			2.0 * 3.14159265358979323846 * i / 997.0);
		const double value = genericPulse.Step(-0.017, duty,
			i % 113 == 0 ? 0.371 : -1.0);
		reversePulseFinite = reversePulseFinite && std::isfinite(value) &&
			genericPulse.Phase() >= 0.0 && genericPulse.Phase() < 1.0;
	}
	Check(reversePulseFinite,
		"generic pulse supports reverse PWM and fractional hard sync");

	genericSaw.Reset();
	genericPulse.Reset();
	bool pulseSawPhaseLocked = true;
	for (int i = 0; i < 20000; ++i)
	{
		const double increment = 0.03 * std::sin(
			2.0 * 3.14159265358979323846 * i / 701.0);
		const double sync = i % 137 == 0 ? 0.417 : -1.0;
		genericSaw.Step(increment, sync);
		genericPulse.Step(increment, 0.2 + 0.6 * (i % 503) / 502.0, sync);
		pulseSawPhaseLocked = pulseSawPhaseLocked &&
			std::abs(genericPulse.Phase() - genericSaw.Phase()) < 1.0e-12;
	}
	Check(pulseSawPhaseLocked,
		"generic saw and pulse stay phase-locked through FM and hard sync");

	genericPulse.Reset();
	genericPulse.Step(0.25, 0.5);
	for (int i = 0; i < 32; ++i)
		genericPulse.Step(0.0, 0.1);
	const double lowAfterMovingThreshold = genericPulse.Step(0.0, 0.1);
	for (int i = 0; i < 32; ++i)
		genericPulse.Step(0.0, 0.9);
	const double highAfterMovingThreshold = genericPulse.Step(0.0, 0.9);
	Check(lowAfterMovingThreshold < -0.999 && highAfterMovingThreshold > 0.999,
		"generic pulse detects comparator edges caused by PWM alone");

	// Compare periodic harmonic magnitudes with the continuous-time bipolar
	// pulse series. A naive sampled comparator quantizes pulse width and folds
	// its ultrasonic harmonics; the fractional minBLEP oscillator should track
	// the analytic in-band spectrum substantially more closely.
	constexpr double spectralFrequency = 1500.0;
	constexpr double spectralDuty = 0.3;
	constexpr int spectralSamples = pulseSampleRate;
	constexpr int spectralWarmup = 256;
	genericPulse.Reset();
	double naivePhase = 0.0;
	for (int i = 0; i < spectralWarmup; ++i)
	{
		genericPulse.Step(spectralFrequency / pulseSampleRate, spectralDuty);
		naivePhase += spectralFrequency / pulseSampleRate;
		naivePhase -= std::floor(naivePhase);
	}
	std::vector<double> bandlimitedPulse(spectralSamples);
	std::vector<double> naivePulse(spectralSamples);
	for (int i = 0; i < spectralSamples; ++i)
	{
		bandlimitedPulse[i] = genericPulse.Step(
			spectralFrequency / pulseSampleRate, spectralDuty);
		naivePhase += spectralFrequency / pulseSampleRate;
		naivePhase -= std::floor(naivePhase);
		naivePulse[i] = naivePhase < spectralDuty ? 1.0 : -1.0;
	}
	double bandlimitedSpectrumError = 0.0;
	double naiveSpectrumError = 0.0;
	for (int harmonic = 1; harmonic < 16; ++harmonic)
	{
		const double expected = 2.0 * std::abs(std::sin(
			3.14159265358979323846 * harmonic * spectralDuty)) /
			(3.14159265358979323846 * harmonic);
		const double normalizedFrequency = harmonic * spectralFrequency /
			pulseSampleRate;
		bandlimitedSpectrumError += std::abs(HarmonicMagnitude(
			bandlimitedPulse, normalizedFrequency) - expected);
		naiveSpectrumError += std::abs(HarmonicMagnitude(
			naivePulse, normalizedFrequency) - expected);
	}
	Check(bandlimitedSpectrumError < 0.5 * naiveSpectrumError,
		"generic pulse improves in-band spectrum over a sampled comparator");

	tfdsp::BandlimitedFixedPulseOscillator fixedPulse;
	fixedPulse.Reset();
	naivePhase = 0.0;
	for (int i = 0; i < spectralWarmup; ++i)
	{
		fixedPulse.Step(spectralFrequency / pulseSampleRate, spectralDuty);
		naivePhase += spectralFrequency / pulseSampleRate;
		naivePhase -= std::floor(naivePhase);
	}
	double fixedPulseSpectrumError = 0.0;
	bandlimitedPulse.assign(spectralSamples, 0.0);
	naivePulse.assign(spectralSamples, 0.0);
	for (int i = 0; i < spectralSamples; ++i)
	{
		bandlimitedPulse[i] = fixedPulse.Step(
			spectralFrequency / pulseSampleRate, spectralDuty);
		naivePhase += spectralFrequency / pulseSampleRate;
		naivePhase -= std::floor(naivePhase);
		naivePulse[i] = naivePhase < spectralDuty ? 1.0 : -1.0;
	}
	naiveSpectrumError = 0.0;
	for (int harmonic = 1; harmonic < 16; ++harmonic)
	{
		const double expected = 2.0 * std::abs(std::sin(
			3.14159265358979323846 * harmonic * spectralDuty)) /
			(3.14159265358979323846 * harmonic);
		const double normalizedFrequency = harmonic * spectralFrequency /
			pulseSampleRate;
		fixedPulseSpectrumError += std::abs(HarmonicMagnitude(
			bandlimitedPulse, normalizedFrequency) - expected);
		naiveSpectrumError += std::abs(HarmonicMagnitude(
			naivePulse, normalizedFrequency) - expected);
	}
	Check(fixedPulseSpectrumError < 0.65 * naiveSpectrumError,
		"fixed polyBLEP pulse improves in-band spectrum over a sampled comparator");

	genericPulse.Reset();
	Check(genericPulse.Step(std::numeric_limits<double>::quiet_NaN(), 0.5) ==
		0.0 && genericPulse.Phase() == 0.0,
		"generic pulse resets safely after non-finite phase input");
	Check(genericPulse.Step(0.01, std::numeric_limits<double>::infinity()) ==
		0.0 && genericPulse.Phase() == 0.0,
		"generic pulse resets safely after non-finite duty input");
	genericPulse.Step(0.01, -1.0);
	Check(genericPulse.DutyCycle() ==
		tfdsp::BandlimitedPulseOscillator<>::MinimumDutyCycle,
		"generic pulse safely clamps duty cycle away from degenerate endpoints");

	tfdsp::BandlimitedTriangleOscillator genericTriangle;
	constexpr double triangleIncrement = 1500.0 / 48000.0;
	for (int i = 0; i < 256; ++i)
		genericTriangle.Step(triangleIncrement);
	std::vector<double> bandlimitedTriangle(48000);
	std::vector<double> naiveTriangle(48000);
	double naiveTrianglePhase = genericTriangle.OutputPhase();
	for (int i = 0; i < 48000; ++i)
	{
		bandlimitedTriangle[i] = genericTriangle.Step(triangleIncrement);
		naiveTriangle[i] =
			tfdsp::BandlimitedTriangleOscillator::RawTriangle(naiveTrianglePhase);
		naiveTrianglePhase += triangleIncrement;
		naiveTrianglePhase -= std::floor(naiveTrianglePhase);
	}
	double bandlimitedTriangleError = 0.0;
	double naiveTriangleError = 0.0;
	for (int harmonic = 1; harmonic < 16; harmonic += 2)
	{
		const double expected = 4.0 /
			(3.14159265358979323846 * 3.14159265358979323846 *
				harmonic * harmonic);
		const double normalizedFrequency = harmonic * 1500.0 / 48000.0;
		bandlimitedTriangleError += std::abs(HarmonicMagnitude(
			bandlimitedTriangle, normalizedFrequency) - expected);
		naiveTriangleError += std::abs(HarmonicMagnitude(
			naiveTriangle, normalizedFrequency) - expected);
	}
	Check(bandlimitedTriangleError < naiveTriangleError,
		"polyBLAMP triangle improves the analytic in-band spectrum");
	genericTriangle.Reset();
	bool reverseTriangleFinite = true;
	for (int i = 0; i < 10000; ++i)
	{
		const double value = genericTriangle.Step(-0.017);
		reverseTriangleFinite = reverseTriangleFinite && std::isfinite(value) &&
			genericTriangle.Phase() >= 0.0 && genericTriangle.Phase() < 1.0;
	}
	Check(reverseTriangleFinite,
		"polyBLAMP triangle supports reverse phase motion");

	tfdsp::Wavefolder::PrepareTable();
	for (int exponent = -30; exponent <= 300; exponent += 3)
	{
		const double input = std::exp(exponent * std::log(10.0) / 10.0);
		const double w = tfdsp::wavefolder_detail::PrincipalLambertW(input);
		Check(std::abs(w * std::exp(w) / input - 1.0) < 2.0e-13,
			"Lambert W table initializer satisfies its defining equation");
	}
	Check(std::abs(tfdsp::wavefolder_detail::PrincipalLambertW(
		2.71828182845904523536) - 1.0) < 1.0e-14,
		"Lambert W initializer returns W(e) = 1");
	for (int characterIndex = 0; characterIndex <
		static_cast<int>(tfdsp::WavefolderCharacter::Count); ++characterIndex)
	{
		const auto character =
			static_cast<tfdsp::WavefolderCharacter>(characterIndex);
		double maximumOddSymmetryError = 0.0;
		double maximumPrimitiveDerivativeError = 0.0;
		for (int index = 1; index <= 1000; ++index)
		{
			const double input = 12.0 * index / 1000.0;
			maximumOddSymmetryError = std::max(maximumOddSymmetryError,
				std::abs(tfdsp::Wavefolder::Transfer(input, character) +
					tfdsp::Wavefolder::Transfer(-input, character)));
			const double h = 1.0e-5;
			const double derivative =
				(tfdsp::Wavefolder::Primitive(input + h, character) -
					tfdsp::Wavefolder::Primitive(input - h, character)) /
				(2.0 * h);
			maximumPrimitiveDerivativeError = std::max(
				maximumPrimitiveDerivativeError,
				std::abs(derivative -
					tfdsp::Wavefolder::Transfer(input, character)));
		}
		Check(maximumOddSymmetryError < 2.0e-12,
			"each wavefolder character remains odd symmetric");
		Check(maximumPrimitiveDerivativeError < 1.0e-7,
			"each wavefolder primitive differentiates to its transfer");
		tfdsp::Wavefolder testFolder;
		double maximumConstantError = 0.0;
		for (int i = 0; i < 1000; ++i)
			maximumConstantError = std::max(maximumConstantError,
				std::abs(testFolder.Process(0.37, character) -
					tfdsp::Wavefolder::Transfer(0.37, character)));
		Check(maximumConstantError < 1.0e-12,
			"wavefolder ADAA preserves a constant input for every character");
		const double h = 1.0e-5;
		const double centralDerivative =
			(tfdsp::Wavefolder::Transfer(h, character) -
				tfdsp::Wavefolder::Transfer(-h, character)) / (2.0 * h);
		Check(std::abs(centralDerivative - 1.0) < 2.0e-5,
			"wavefolder characters share unity small-signal gain");
	}
	Check(std::abs(2.0 * tfdsp::Wavefolder::Transfer(0.5) - 1.0) < 0.01,
		"zero-fold drive and makeup retain an almost-unity endpoint");
	Check(std::abs(tfdsp::Wavefolder::Transfer(1.0,
		tfdsp::WavefolderCharacter::Lockhart) - 0.88499658096) < 2.0e-6,
		"Lockhart topology retains its normalized first-fold level");
	Check(std::abs(tfdsp::Wavefolder::Transfer(1.0,
		tfdsp::WavefolderCharacter::Serge) - 0.52736535702) < 2.0e-6,
		"Serge topology retains its normalized first-fold level");

	tfdsp::WavefoldOscillator<tfdsp::X2Resampler_Order7> foldedOscillator(
		tfdsp::CreateX2Resampler_Chebychev7);
	foldedOscillator.SetSampleRate(48000.0);
	bool foldedOscillatorFinite = true;
	for (int i = 0; i < 48000; ++i)
	{
		const double value = foldedOscillator.Step(261.625565,
			0.5 + 0.5 * std::sin(2.0 * 3.14159265358979323846 * i / 997.0),
			0.5 + 0.5 * std::sin(2.0 * 3.14159265358979323846 * i / 733.0),
			0.5 * std::sin(2.0 * 3.14159265358979323846 * i / 601.0));
		foldedOscillatorFinite = foldedOscillatorFinite && std::isfinite(value);
	}
	Check(foldedOscillatorFinite,
		"oversampled wavefold oscillator remains finite under modulation");
	tfdsp::WavefoldOscillator<tfdsp::X2Resampler_Order7> normalledOscillator(
		tfdsp::CreateX2Resampler_Chebychev7);
	tfdsp::WavefoldOscillator<tfdsp::X2Resampler_Order7> detailedOscillator(
		tfdsp::CreateX2Resampler_Chebychev7);
	tfdsp::WavefoldOscillator<tfdsp::X2Resampler_Order7> foldedOnlyOscillator(
		tfdsp::CreateX2Resampler_Chebychev7);
	normalledOscillator.SetSampleRate(48000.0);
	detailedOscillator.SetSampleRate(48000.0);
	foldedOnlyOscillator.SetSampleRate(48000.0);
	double normalledDifference = 0.0;
	double omittedOutputDifference = 0.0;
	bool externalPathFinite = true;
	for (int i = 0; i < 10000; ++i)
	{
		const double normalled = normalledOscillator.Step(
			440.0, 0.4, 0.7, -0.2);
		const auto detailed = detailedOscillator.StepWithInput(
			440.0, 0.4, 0.7, -0.2, 0.0, false);
		const auto foldedOnly = foldedOnlyOscillator.StepWithInput(
			440.0, 0.4, 0.7, -0.2, 0.0, false, false, true);
		normalledDifference = std::max(normalledDifference,
			std::abs(normalled - detailed.folded));
		omittedOutputDifference = std::max(omittedOutputDifference,
			std::abs(foldedOnly.folded - detailed.folded));
		externalPathFinite = externalPathFinite &&
			std::isfinite(detailed.oscillator) &&
			std::isfinite(detailed.folded);
	}
	Check(normalledDifference < 1.0e-12,
		"normalled folder input preserves the original oscillator path");
	Check(omittedOutputDifference < 1.0e-12,
		"omitting the unused raw output does not change the folded path");
	Check(externalPathFinite,
		"wavefold oscillator exposes finite oscillator and folder outputs");
	for (int i = 0; i < 10000; ++i)
	{
		const double external = std::sin(
			2.0 * 3.14159265358979323846 * 997.0 * i / 48000.0);
		const auto rendered = detailedOscillator.StepWithInput(
			440.0, 0.4, 1.0, 0.5, external, true);
		externalPathFinite = externalPathFinite &&
			std::isfinite(rendered.oscillator) && std::isfinite(rendered.folded);
	}
	Check(externalPathFinite,
		"reconstructed external folder input remains finite");
	using X4Wavefolder =
		tfdsp::WavefoldOscillator<tfdsp::X4Resampler_Order7>;
	Check(std::abs(X4Wavefolder::FoldScaleForFrequency(100.0) - 1.0) <
		1.0e-12, "low notes retain the requested fold depth");
	Check(std::abs(X4Wavefolder::FoldScaleForFrequency(1000.0) - 0.625) <
		1.0e-12, "fold taper uses a fixed 6 kHz harmonic budget");
	Check(std::abs(X4Wavefolder::FoldScaleForFrequency(2000.0) - 0.25) <
		1.0e-12, "fold taper follows musical frequency");
	Check(std::abs(X4Wavefolder::FoldScaleForFrequency(4000.0) - 0.0625) <
		1.0e-12, "high notes retain only a shallow fold");
	Check(X4Wavefolder::FoldScaleForFrequency(6000.0) == 0.0,
		"fold taper reaches zero at its harmonic budget");

	tfdsp::InterpolatedOrnsteinUhlenbeck ou;
	ou.Configure(48000.0, 60.0, 2.0, 100.0);
	CountingGenerator rng;
	double drift = 0.0;
	for (int i = 0; i < 48000; ++i)
		drift = ou.Step(rng);
	Check(std::isfinite(drift), "OU drift remains finite");
	Check(rng.calls < 1000, "OU noise generation runs at control rate");
	const double boundedPositive = tfdsp::ApplyBoundedDrift(0.5, 1.0, 1.0);
	Check(std::abs(boundedPositive -
		(0.5 + 0.5 * std::tanh(1.0))) < 1.0e-12,
		"bounded drift scales toward the available headroom");
	Check(tfdsp::ApplyBoundedDrift(0.37, 100.0, 0.0) == 0.37,
		"zero drift depth preserves the control exactly");
	Check(tfdsp::ApplyBoundedDrift(-0.2, -100.0, 1.0, -1.0, 1.0) >= -1.0 &&
		tfdsp::ApplyBoundedDrift(-0.2, 100.0, 1.0, -1.0, 1.0) <= 1.0,
		"bounded drift approaches parameter limits without crossing them");

	tfdsp::InterpolatedOrnsteinUhlenbeck stationaryOu;
	stationaryOu.ConfigureStationary(100.0, 0.5, 2.0, 100.0);
	CountingGenerator stationaryRng;
	double sum = 0.0;
	double sumSquares = 0.0;
	constexpr int stationarySamples = 200000;
	constexpr int stationaryWarmup = 10000;
	for (int i = 0; i < stationarySamples + stationaryWarmup; ++i)
	{
		const double value = stationaryOu.Step(stationaryRng);
		if (i >= stationaryWarmup)
		{
			sum += value;
			sumSquares += value * value;
		}
	}
	const double stationaryMean = sum / stationarySamples;
	const double stationaryVariance = sumSquares / stationarySamples -
		stationaryMean * stationaryMean;
	Check(std::abs(stationaryMean) < 0.12,
		"stationary OU process remains centred");
	Check(std::abs(stationaryVariance - 4.0) < 0.35,
		"stationary OU configuration preserves requested variance");

	tfdsp::SmoothOrnsteinUhlenbeck smoothOu;
	smoothOu.ConfigureStationary(1000.0, 0.5, 1.0, 100.0);
	CountingGenerator smoothRng;
	double smoothSum = 0.0;
	double smoothSumSquares = 0.0;
	double smoothDifferenceSquares = 0.0;
	double previousSmooth = 0.0;
	constexpr int smoothSamples = 500000;
	constexpr int smoothWarmup = 20000;
	for (int i = 0; i < smoothSamples + smoothWarmup; ++i)
	{
		const double value = smoothOu.Step(smoothRng);
		if (i >= smoothWarmup)
		{
			smoothSum += value;
			smoothSumSquares += value * value;
			const double difference = value - previousSmooth;
			smoothDifferenceSquares += difference * difference;
		}
		previousSmooth = value;
	}
	const double smoothMean = smoothSum / smoothSamples;
	const double smoothVariance = smoothSumSquares / smoothSamples -
		smoothMean * smoothMean;
	const double smoothDifferenceRms = std::sqrt(
		smoothDifferenceSquares / smoothSamples);
	if (std::abs(smoothMean) >= 0.25 ||
		std::abs(smoothVariance - 1.0) >= 0.12 || smoothDifferenceRms >= 0.004)
		std::cerr << "smooth OU stats: mean=" << smoothMean <<
			", variance=" << smoothVariance <<
			", difference RMS=" << smoothDifferenceRms << '\n';
	Check(std::abs(smoothMean) < 0.25,
		"smooth OU process remains centred");
	Check(std::abs(smoothVariance - 1.0) < 0.12,
		"smooth OU process preserves requested stationary variance");
	Check(smoothDifferenceRms < 0.004,
		"smooth OU process suppresses rapid parameter movement");
	Check(tfdsp::UnisonSpreadCents(0.0) == 0.0 &&
		std::abs(tfdsp::UnisonSpreadCents(0.5) - 7.1875) < 1.0e-12 &&
		tfdsp::UnisonSpreadCents(1.0) == 50.0,
		"unison spread provides fine low-range control and full excursion");
	Check(std::abs(tfdsp::UnisonSpreadCents(
		tfdsp::UnisonSpreadControlForCents(4.0)) - 4.0) < 1.0e-8,
		"unison spread display mapping is invertible");
	for (int voices = 1; voices <= tfdsp::MaximumUnisonVoices; ++voices)
	{
		const auto positions = tfdsp::UnisonPitchPositions(voices);
		double mean = 0.0;
		for (int voice = 0; voice < voices; ++voice)
			mean += positions[voice];
		Check(std::abs(mean) < 1.0e-12,
			"unison pitch layouts preserve their tuning centre");
		Check(std::abs(tfdsp::UnisonOutputGain(voices) -
			1.0 / std::sqrt(static_cast<double>(voices))) < 1.0e-12,
			"unison output uses energy normalization");
	}
	for (int count = 1; count <= tfdsp::MaximumStackedOscillatorVoices; ++count)
	{
		const auto pitch = tfdsp::StackedOscillatorPitchPositions(count);
		const auto pan = tfdsp::StackedOscillatorPanPositions(count);
		const auto tracking = tfdsp::StackedOscillatorTrackingPositions(count);
		double pitchMean = 0.0;
		double panMean = 0.0;
		double trackingMean = 0.0;
		for (int voice = 0; voice < count; ++voice)
		{
			pitchMean += pitch[voice];
			panMean += pan[voice];
			trackingMean += tracking[voice];
			Check(std::abs(pitch[voice]) <= 1.0 + 1.0e-12 &&
				std::abs(pan[voice]) <= 1.0 + 1.0e-12 &&
				std::abs(tracking[voice]) <= 1.0 + 1.0e-12,
				"large-stack layouts stay normalized");
		}
		Check(std::abs(pitchMean) < 1.0e-12 &&
			std::abs(panMean) < 1.0e-12 &&
			std::abs(trackingMean) < 1.0e-12,
			"large-stack layouts preserve pitch, pan, and tracking centres");
	}

	tfdsp::StackedOscillatorVoice stackedVoice;
	constexpr double stackedFrequency = 440.0;
	constexpr double stackedIncrement = stackedFrequency / 48000.0;
	bool stackedFinite = true;
	for (int sample = 0; sample < 48000; ++sample)
	{
		const double pwm = 0.5 + 0.4 * std::sin(
			2.0 * 3.14159265358979323846 * sample / 16000.0);
		const auto rendered = stackedVoice.Step(stackedIncrement, pwm, 0.7);
		stackedFinite = stackedFinite && std::isfinite(rendered.main) &&
			std::isfinite(rendered.sub);
	}
	Check(stackedFinite,
		"stacked oscillator remains finite under deep pulse-width modulation");
	stackedVoice.Reset();
	double previousSub = 0.0;
	int positiveSubCrossings = 0;
	for (int sample = 0; sample < 48000; ++sample)
	{
		const double sub = stackedVoice.StepSub(stackedIncrement);
		if (previousSub <= 0.0 && sub > 0.0)
			++positiveSubCrossings;
		previousSub = sub;
	}
	Check(std::abs(positiveSubCrossings - 220) <= 1,
		"stacked oscillator sub runs exactly one octave below its parent");
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
	double maxTanPrewarpRelativeError = 0.0;
	for (int i = 1; i <= 20000; ++i)
	{
		const double argument = 0.225 * 3.14159265358979323846 * i / 20000.0;
		maxTanPrewarpRelativeError = std::max(maxTanPrewarpRelativeError,
			std::abs(tfdsp::TanPrewarp(argument) / std::tan(argument) - 1.0));
	}
	Check(maxTanPrewarpRelativeError < 2.0e-8,
		"fast filter prewarp stays within its relative-error budget");
	double maxSinTwoPiAbsoluteError = 0.0;
	for (int i = -20000; i <= 20000; ++i)
	{
		const double phase = 4.0 * i / 20000.0;
		maxSinTwoPiAbsoluteError = std::max(maxSinTwoPiAbsoluteError,
			std::abs(tfdsp::SinTwoPi(phase) - std::sin(
				6.283185307179586476925286766559 * phase)));
	}
	Check(maxSinTwoPiAbsoluteError < 7.0e-10,
		"fast control sine stays within its absolute-error budget");

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
			20000.0, 1.0, 15.848931924611133));
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
			0.0, 1.0, 0.0, 0.0,
			[](double, double, double) { return 1000.0; });
		arpPostSafetyPeak = std::max(arpPostSafetyPeak,
			std::abs(static_cast<double>(rendered.postProcessed)));
	}
	Check(arpPostSafetyPeak > 11.0 &&
		arpPostSafetyPeak < tfdsp::RackOutputAdapter::CableLimitVolts,
		"ARP integrated post-processor output stays within Rack cable headroom");
	Check(tfdsp::RackOutputAdapter::ProcessOversampled(10.0) == 10.0,
		"Rack output adapter is linear at normal full-scale level");
	Check(tfdsp::RackOutputAdapter::ProcessOversampled(10.5) == 10.5 &&
		tfdsp::RackOutputAdapter::ProcessPostDecimation(11.5) == 11.5,
		"Rack output adapter is continuous and linear through both knees");
	Check(tfdsp::RackOutputAdapter::ProcessOversampled(-100.0) ==
		-tfdsp::RackOutputAdapter::ProcessOversampled(100.0),
		"Rack output adapter is symmetric");
	Check(tfdsp::RackOutputAdapter::ProcessOversampled(1000.0) <
		tfdsp::RackOutputAdapter::OversampledLimitVolts,
		"oversampled Rack output reserves decimator headroom");
	Check(tfdsp::RackOutputAdapter::ProcessPostDecimation(1000.0) <
		tfdsp::RackOutputAdapter::CableLimitVolts,
		"post-decimation Rack output stays below the protected-rail limit");

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
	tfdsp::ArpEnvelope adEnvelope;
	adEnvelope.SetSampleRate(48000.0);
	adEnvelope.SetMode(tfdsp::ArpEnvelope::Mode::Ad);
	for (int i = 0; i < 4800; ++i)
		envelope = adEnvelope.Step(10.0, 0.0, 0.1, 0.2, 0.75, 0.5);
	Check(envelope > 0.999 &&
		adEnvelope.GetStage() == tfdsp::ArpEnvelope::Stage::Decay,
		"ARP AD reaches its peak and enters decay while Gate remains high");
	for (int i = 0; i < 9600; ++i)
		envelope = adEnvelope.Step(10.0, 0.0, 0.1, 0.2, 0.75, 0.5);
	Check(envelope < 0.001 &&
		adEnvelope.GetStage() == tfdsp::ArpEnvelope::Stage::Idle,
		"ARP AD decays to zero without waiting for Gate to fall");
	envelope = adEnvelope.Step(10.0, 10.0, 0.1, 0.2, 0.75, 0.5);
	Check(envelope > 0.0 &&
		adEnvelope.GetStage() == tfdsp::ArpEnvelope::Stage::Attack,
		"ARP AD retriggers from Trigger while Gate remains high");
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
	const auto unpatchedExponentialEnvelope =
		tfdsp::RouteArp4019ControlVoltages(10.0, 2.0, false, false, true,
			false);
	Check(unpatchedExponentialEnvelope.linear == 0.0 &&
		unpatchedExponentialEnvelope.exponential == 10.0,
		"unpatched VCA modulation retains the internal envelope law");
	const auto addLinearControls = tfdsp::RouteArp4019ControlVoltages(
		10.0, 2.0, true, true, false, false);
	Check(addLinearControls.linear == 12.0 &&
		addLinearControls.exponential == 0.0,
		"VCA add routing sums linear envelope and modulation controls");
	const auto addSplitControls = tfdsp::RouteArp4019ControlVoltages(
		10.0, 2.0, true, true, true, false);
	Check(addSplitControls.linear == 2.0 &&
		addSplitControls.exponential == 10.0,
		"VCA add routing preserves independent envelope and modulation laws");
	const auto addExponentialControls = tfdsp::RouteArp4019ControlVoltages(
		10.0, 2.0, true, true, true, true);
	Check(addExponentialControls.linear == 0.0 &&
		addExponentialControls.exponential == 12.0,
		"VCA add routing sums controls which share the exponential input");
	const auto replaceWithExternal = tfdsp::RouteArp4019ControlVoltages(
		10.0, 2.0, true, false, false, true);
	Check(replaceWithExternal.linear == 0.0 &&
		replaceWithExternal.exponential == 2.0,
		"VCA EXT routing replaces the envelope and retains modulation law");
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
			8000.0, 0.0, 1.0, 10.0, -10.0,
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
	Arp4072 arpVoiceFilterWithoutLowPass(tfdsp::CreateX4Resampler_Cheby7);
	Arp4072 arpVoiceFilterWithLowPass(tfdsp::CreateX4Resampler_Cheby7);
	Arp4019 arpVoiceVcaWithoutLowPass(tfdsp::CreateX4Resampler_Cheby7);
	Arp4019 arpVoiceVcaWithLowPass(tfdsp::CreateX4Resampler_Cheby7);
	arpVoiceFilterWithoutLowPass.SetSampleRate(48000.0);
	arpVoiceFilterWithLowPass.SetSampleRate(48000.0);
	arpVoiceVcaWithoutLowPass.SetSampleRate(48000.0);
	arpVoiceVcaWithLowPass.SetSampleRate(48000.0);
	double omittedLowPassDifference = 0.0;
	bool omittedLowPassZero = true;
	for (int i = 0; i < 4096; ++i)
	{
		const double input = 0.1 * std::sin(2.0 *
			3.14159265358979323846 * 330.0 * i / 48000.0);
		const auto withoutLowPass =
			arpVoiceFilterWithoutLowPass.StepWithPostProcessor(
				input, 2400.0, 0.35, 2.0, 0.0, 10.0,
				[&](double filtered, double linearCv, double exponentialCv)
				{
					return arpVoiceVcaWithoutLowPass.ProcessOversampled(
						filtered, linearCv, exponentialCv);
				}, false);
		const auto withLowPass = arpVoiceFilterWithLowPass.StepWithPostProcessor(
			input, 2400.0, 0.35, 2.0, 0.0, 10.0,
			[&](double filtered, double linearCv, double exponentialCv)
			{
				return arpVoiceVcaWithLowPass.ProcessOversampled(
					filtered, linearCv, exponentialCv);
			});
		omittedLowPassDifference = std::max(omittedLowPassDifference,
			std::abs(static_cast<double>(withoutLowPass.postProcessed) -
				static_cast<double>(withLowPass.postProcessed)));
		omittedLowPassZero = omittedLowPassZero &&
			withoutLowPass.lowPass == 0.0f;
	}
	Check(omittedLowPassZero, "omitted ARP low-pass output remains zero");
	Check(omittedLowPassDifference < 1.0e-12,
		"omitting the ARP low-pass decimator does not change its VCA path");

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
			8000.0, 0.0, 1.0, 10.0, -10.0,
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
	Check(tfdsp::RackOutputAdapter::ProcessPostDecimation(100.0) <
		tfdsp::RackOutputAdapter::CableLimitVolts,
		"shared Rack output stage approaches its limit without hard clipping");
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
		overloadedPeak < tfdsp::RackOutputAdapter::CableLimitVolts,
		"TB-303 VCA overload bends smoothly inside the safety guard");
	Check(tb303Vca.Step(std::numeric_limits<double>::quiet_NaN(),
		1.0, 0.0) == 0.0,
		"TB-303 VCA wrapper rejects non-finite input");

	if (failures == 0)
		std::cout << "All TriggerFish DSP tests passed\n";
	return failures == 0 ? 0 : 1;
}
