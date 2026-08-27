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
}

int main()
{
	{
		double maximumTanhError = 0.0;
		double maximumFalloffRelativeError = 0.0;
		for (int index = 0; index <= 20000; ++index)
		{
			const double argument = -10.0 + 20.0 * index / 20000.0;
			maximumTanhError = std::max(maximumTanhError, std::abs(
				tfdsp::TanhPade76(argument) - std::tanh(argument)));
			const double radius = 0.5 + 1.5 * index / 20000.0;
			const double reference = std::pow(radius, -1.3);
			maximumFalloffRelativeError = std::max(maximumFalloffRelativeError,
				std::abs(tfdsp::PowNegative1p3(radius) - reference) / reference);
		}
		Check(maximumTanhError < 1.1e-4,
			"electric piano pickup tanh approximation stays below its error bound");
		Check(maximumFalloffRelativeError < 5.7e-5,
			"electric piano pickup radial approximation stays below its error bound");

		tfdsp::ElectricPianoControls controls;
		tfdsp::ElectricPianoVoice softVoice;
		tfdsp::ElectricPianoVoice hardVoice;
		softVoice.SetSampleRate(48000.0);
		hardVoice.SetSampleRate(48000.0);
		double softEnergy = 0.0;
		double hardEnergy = 0.0;
		bool finite = true;
		for (int sample = 0; sample < 24000; ++sample)
		{
			const double gate = sample < 12000 ? 10.0 : 0.0;
			const double soft = softVoice.Step(0.0, gate, 0.18, false, controls);
			const double hard = hardVoice.Step(0.0, gate, 1.0, false, controls);
			finite = finite && std::isfinite(soft) && std::isfinite(hard);
			if (sample < 4800)
			{
				softEnergy += soft * soft;
				hardEnergy += hard * hard;
			}
		}
		Check(finite, "electric piano voice remains finite through strike and release");
		Check(softEnergy > 1.0e-8,
			"electric piano soft strike produces audible energy");
		Check(hardEnergy > 2.0 * softEnergy,
			"electric piano velocity increases note energy dynamically");

		const auto lowKey = tfdsp::MakeElectricPianoKeyProfile(-3.0);
		const auto middleKey = tfdsp::MakeElectricPianoKeyProfile(0.0);
		const auto highKey = tfdsp::MakeElectricPianoKeyProfile(3.0);
		auto displacementPerImpulse = [](const auto& key)
		{
			return 1.0 / (key.modalMassRatio * key.fundamentalHz);
		};
		Check(lowKey.modalMassRatio > middleKey.modalMassRatio &&
			middleKey.modalMassRatio > highKey.modalMassRatio &&
			displacementPerImpulse(lowKey) >
				displacementPerImpulse(middleKey) &&
			displacementPerImpulse(middleKey) >
				displacementPerImpulse(highKey) &&
			lowKey.pickupSensitivity > middleKey.pickupSensitivity &&
			highKey.pickupSensitivity > middleKey.pickupSensitivity,
			"electric piano key model separates mechanical impedance from per-key pickup calibration");
		Check(tfdsp::ElectricPianoModeBandlimitGain(0.31 * 48000.0,
			48000.0) == 1.0 &&
			tfdsp::ElectricPianoModeBandlimitGain(0.385 * 48000.0,
				48000.0) > 0.0 &&
			tfdsp::ElectricPianoModeBandlimitGain(0.385 * 48000.0,
				48000.0) < 1.0 &&
			tfdsp::ElectricPianoModeBandlimitGain(0.45 * 48000.0,
				48000.0) == 0.0,
			"electric piano high modes taper smoothly before Nyquist");

		const auto isolatedFork = tfdsp::MakeElectricPianoCoupledForkProfile(
			0.0, middleKey.keyboardPosition);
		const auto factoryFork = tfdsp::MakeElectricPianoCoupledForkProfile(
			controls.coupling, middleKey.keyboardPosition);
		const auto coupledFork = tfdsp::MakeElectricPianoCoupledForkProfile(
			1.0, middleKey.keyboardPosition);
		Check(std::abs(controls.coupling - 0.5) < 1.0e-12,
			"electric piano Coupling defaults to the calibrated control midpoint");
		auto modalCompleteness = [](const auto& profile)
		{
			return profile.inverseModalMassRatios[0] +
				profile.inverseModalMassRatios[1];
		};
		auto massOrthogonality = [](const auto& profile)
		{
			return 1.0 + profile.toneBarModalMassRatio *
				profile.toneBarDisplacementRatios[0] *
				profile.toneBarDisplacementRatios[1];
		};
		Check(std::abs(modalCompleteness(isolatedFork) - 1.0) < 1.0e-9 &&
			std::abs(modalCompleteness(factoryFork) - 1.0) < 1.0e-9 &&
			std::abs(modalCompleteness(coupledFork) - 1.0) < 1.0e-9 &&
			std::abs(massOrthogonality(factoryFork)) < 1.0e-9,
			"electric piano coupled-fork eigenvectors conserve tine impulse and are mass orthogonal");
		Check(coupledFork.frequencyRatios[0] <
			factoryFork.frequencyRatios[0] &&
			factoryFork.frequencyRatios[0] < isolatedFork.frequencyRatios[0] &&
			coupledFork.inverseModalMassRatios[0] >
				2.0 * factoryFork.inverseModalMassRatios[0] &&
			factoryFork.inverseModalMassRatios[0] >
				10.0 * isolatedFork.inverseModalMassRatios[0],
			"electric piano Coupling increases tone-bar hybridization and normal-mode separation");
		Check(isolatedFork.supportReactionLossFactors[1] >
			3.0 * factoryFork.supportReactionLossFactors[1] &&
			factoryFork.supportReactionLossFactors[1] >
			20.0 * coupledFork.supportReactionLossFactors[1],
			"electric piano tone bar cancels lossy common-support reaction at strong coupling");
		for (double sampleRate : {8000.0, 44100.0, 96000.0, 192000.0})
		{
			tfdsp::ElectricPianoVoice rateVoice;
			rateVoice.SetSampleRate(sampleRate);
			auto stressedControls = controls;
			stressedControls.hammer = 1.0;
			stressedControls.proximity = 1.0;
			double peak = 0.0;
			bool rateFinite = true;
			const int samples = static_cast<int>(0.12 * sampleRate);
			for (int sample = 0; sample < samples; ++sample)
			{
				const double output = rateVoice.Step(1.0, 10.0, 1.0, false,
					stressedControls);
				rateFinite = rateFinite && std::isfinite(output);
				peak = std::max(peak, std::abs(output));
			}
			Check(rateFinite && peak > 1.0e-5 && peak < 12.0,
				"electric piano compliant contact and pickup remain stable across sample rates");
		}

		// Continuous Rack CV can contain tiny positive residue far below MIDI's
		// minimum nonzero velocity. It must not start a long, inaudible 16x
		// collision, while a strike immediately above the calibrated hammer-speed
		// floor must still run to natural separation.
		bool subThresholdIsSilent = true;
		bool thresholdContactsComplete = true;
		bool thresholdOutputIsNonzero = true;
		double minimumThresholdDuration =
			std::numeric_limits<double>::infinity();
		double maximumThresholdDuration = 0.0;
		for (double sampleRate : {8000.0, 44100.0, 48000.0, 96000.0, 192000.0})
		{
			auto thresholdControls = controls;
			thresholdControls.hammer = 0.0;
			thresholdControls.mechanics = 0.0;
			tfdsp::ElectricPianoVoice belowThresholdVoice;
			tfdsp::ElectricPianoVoice aboveThresholdVoice;
			belowThresholdVoice.SetSampleRate(sampleRate);
			aboveThresholdVoice.SetSampleRate(sampleRate);
			double belowThresholdEnergy = 0.0;
			for (int sample = 0; sample < 128; ++sample)
			{
				const double output = belowThresholdVoice.Step(-3.0, 10.0,
					0.0010, false, thresholdControls);
				belowThresholdEnergy += output * output;
			}
			subThresholdIsSilent = subThresholdIsSilent &&
				belowThresholdEnergy < 1.0e-20 &&
				!belowThresholdVoice.ContactActive() &&
				belowThresholdVoice.ContactAge() == 0.0;

			bool contactStarted = false;
			bool contactCompleted = false;
			double aboveThresholdEnergy = 0.0;
			const int limit = static_cast<int>(std::ceil(0.100 * sampleRate));
			for (int sample = 0; sample < limit; ++sample)
			{
				const double output = aboveThresholdVoice.Step(-3.0, 10.0,
					0.0011, false, thresholdControls);
				aboveThresholdEnergy += output * output;
				contactStarted = contactStarted ||
					aboveThresholdVoice.ContactActive();
				if (contactStarted && !aboveThresholdVoice.ContactActive())
				{
					contactCompleted = true;
					break;
				}
			}
			const double duration = aboveThresholdVoice.ContactAge();
			minimumThresholdDuration = std::min(minimumThresholdDuration,
				duration);
			maximumThresholdDuration = std::max(maximumThresholdDuration,
				duration);
			thresholdContactsComplete = thresholdContactsComplete &&
				contactStarted && contactCompleted;
			thresholdOutputIsNonzero = thresholdOutputIsNonzero &&
				aboveThresholdEnergy > 1.0e-12;
		}
		Check(subThresholdIsSilent && thresholdContactsComplete &&
			thresholdOutputIsNonzero &&
			maximumThresholdDuration - minimumThresholdDuration < 0.00020,
			"electric piano rejects inaudible velocity residue without truncating a valid collision");

		// Contact must end by physical separation throughout the supported range.
		// This grid includes the soft bass combinations that previously reached
		// an arbitrary 12 ms cutoff while the hammer was still compressed.
		constexpr std::array<double, 5> ContactSampleRates{
			8000.0, 44100.0, 48000.0, 96000.0, 192000.0};
		constexpr std::array<double, 5> ContactPitches{
			-6.0, -3.0, 0.0, 3.0, 6.0};
		constexpr std::array<double, 3> ContactHardnesses{0.0, 0.52, 1.0};
		constexpr std::array<double, 5> ContactVelocities{
			1.0 / 127.0, 0.03, 0.10, 0.30, 1.0};
		bool allContactsComplete = true;
		bool contactDurationsAreRateInvariant = true;
		double longestContact = 0.0;
		for (double pitch : ContactPitches)
		{
			for (double hardness : ContactHardnesses)
			{
				for (double velocity : ContactVelocities)
				{
					double minimumDuration =
						std::numeric_limits<double>::infinity();
					double maximumDuration = 0.0;
					for (double sampleRate : ContactSampleRates)
					{
						auto contactControls = controls;
						contactControls.mechanics = 0.0;
						contactControls.hammer = hardness;
						tfdsp::ElectricPianoVoice contactVoice;
						contactVoice.SetSampleRate(sampleRate);
						bool finiteContact = true;
						bool completed = false;
						const int limit = static_cast<int>(
							std::ceil(0.050 * sampleRate));
						for (int sample = 0; sample < limit; ++sample)
						{
							const double output = contactVoice.Step(pitch, 10.0,
								velocity, false, contactControls);
							finiteContact = finiteContact && std::isfinite(output);
							if (!contactVoice.ContactActive())
							{
								completed = true;
								break;
							}
						}
						const double duration = contactVoice.ContactAge();
						minimumDuration = std::min(minimumDuration, duration);
						maximumDuration = std::max(maximumDuration, duration);
						longestContact = std::max(longestContact, duration);
						for (int sample = 0; sample < 64; ++sample)
							contactVoice.Step(pitch, 10.0, velocity, false,
								contactControls);
						allContactsComplete = allContactsComplete && finiteContact &&
							completed && !contactVoice.ContactActive();
					}
					contactDurationsAreRateInvariant =
						contactDurationsAreRateInvariant &&
						maximumDuration - minimumDuration < 0.00020;
				}
			}
		}
		if (!(allContactsComplete && contactDurationsAreRateInvariant &&
			longestContact < 0.045))
			std::cerr << "electric piano longest physical contact " <<
				longestContact << " seconds\n";
		Check(allContactsComplete && contactDurationsAreRateInvariant &&
			longestContact < 0.045,
			"electric piano hammer separates naturally and consistently across its full control range");

		auto silentControls = controls;
		silentControls.mechanics = 0.0;
		tfdsp::ElectricPianoVoice zeroVelocityVoice;
		zeroVelocityVoice.SetSampleRate(48000.0);
		double zeroVelocityEnergy = 0.0;
		for (int sample = 0; sample < 4096; ++sample)
		{
			const double output = zeroVelocityVoice.Step(0.0, 10.0, 0.0,
				false, silentControls);
			zeroVelocityEnergy += output * output;
		}
		Check(zeroVelocityEnergy < 1.0e-20 &&
			!zeroVelocityVoice.ContactActive(),
			"electric piano zero velocity produces no hammer or tine excitation");

		constexpr std::array<double, 6> TestVelocities{
			0.03, 0.08, 0.18, 0.38, 0.68, 1.0};
		std::array<double, TestVelocities.size()> velocityEnergies{};
		for (std::size_t index = 0; index < TestVelocities.size(); ++index)
		{
			tfdsp::ElectricPianoVoice voice;
			voice.SetSampleRate(48000.0);
			for (int sample = 0; sample < 4800; ++sample)
			{
				const double output = voice.Step(0.0, 10.0,
					TestVelocities[index], false, silentControls);
				velocityEnergies[index] += output * output;
			}
		}
		Check(std::is_sorted(velocityEnergies.begin(), velocityEnergies.end()) &&
			std::adjacent_find(velocityEnergies.begin(), velocityEnergies.end(),
				[](double first, double second) { return !(second > first); }) ==
				velocityEnergies.end(),
			"electric piano output energy rises monotonically across velocity");

		auto farPickupControls = silentControls;
		auto nearPickupControls = silentControls;
		farPickupControls.proximity = 0.0;
		nearPickupControls.proximity = 1.0;
		tfdsp::ElectricPianoVoice farPickupVoice;
		tfdsp::ElectricPianoVoice nearPickupVoice;
		farPickupVoice.SetSampleRate(48000.0);
		nearPickupVoice.SetSampleRate(48000.0);
		double farPickupEnergy = 0.0;
		double nearPickupEnergy = 0.0;
		for (int sample = 0; sample < 4800; ++sample)
		{
			const double farOutput = farPickupVoice.Step(0.0, 10.0, 0.3,
				false, farPickupControls);
			const double nearOutput = nearPickupVoice.Step(0.0, 10.0, 0.3,
				false, nearPickupControls);
			farPickupEnergy += farOutput * farOutput;
			nearPickupEnergy += nearOutput * nearOutput;
		}
		const double proximityEnergyRatio = nearPickupEnergy /
			std::max(1.0e-20, farPickupEnergy);
		if (!(proximityEnergyRatio > 1.3 && proximityEnergyRatio < 2.1))
			std::cerr << "electric piano pickup proximity energy ratio " <<
				proximityEnergyRatio << '\n';
		Check(proximityEnergyRatio > 1.3 && proximityEnergyRatio < 2.1,
			"electric piano pickup Proximity retains curvature without becoming a second Drive");

		tfdsp::ElectricPianoVoice restruckVoice;
		tfdsp::ElectricPianoVoice zeroRetriggerVoice;
		restruckVoice.SetSampleRate(48000.0);
		zeroRetriggerVoice.SetSampleRate(48000.0);
		for (int sample = 0; sample < 1200; ++sample)
		{
			restruckVoice.Step(0.0, 10.0, 1.0, false, silentControls);
			zeroRetriggerVoice.Step(0.0, 10.0, 1.0, false, silentControls);
		}
		restruckVoice.Step(0.0, 0.0, 0.0, false, silentControls);
		zeroRetriggerVoice.Step(0.0, 0.0, 0.0, false, silentControls);
		double retriggerDifference = 0.0;
		double zeroRetriggerEnergy = 0.0;
		for (int sample = 0; sample < 2048; ++sample)
		{
			const double restruck = restruckVoice.Step(0.0, 10.0, 0.65,
				false, silentControls);
			const double unchanged = zeroRetriggerVoice.Step(0.0, 10.0, 0.0,
				false, silentControls);
			const double difference = restruck - unchanged;
			retriggerDifference += difference * difference;
			zeroRetriggerEnergy += unchanged * unchanged;
		}
		Check(zeroRetriggerVoice.Energy() > 1.0e-6 &&
			retriggerDifference > 0.10 * zeroRetriggerEnergy,
			"electric piano restrike preserves existing tine motion and applies a new physical collision");

		auto quietControls = controls;
		quietControls.mechanics = 0.0;
		tfdsp::ElectricPianoVoice bentVoice;
		tfdsp::ElectricPianoVoice referenceVoice;
		bentVoice.SetSampleRate(48000.0);
		referenceVoice.SetSampleRate(48000.0);
		bentVoice.SetNoiseSeed(0x2468aceu);
		referenceVoice.SetNoiseSeed(0x2468aceu);
		double preBendDifference = 0.0;
		double postBendDifference = 0.0;
		for (int sample = 0; sample < 4096; ++sample)
		{
			const double bentPitch = sample < 2048 ? 0.0 : 1.0;
			const double bent = bentVoice.Step(bentPitch, 10.0, 0.8, false,
				quietControls);
			const double reference = referenceVoice.Step(0.0, 10.0, 0.8,
				false, quietControls);
			const double difference = bent - reference;
			if (sample < 2048)
				preBendDifference += difference * difference;
			else
				postBendDifference += difference * difference;
		}
		Check(preBendDifference < 1.0e-20 && postBendDifference > 1.0e-5,
			"electric piano retunes an individual held voice when its pitch changes");

		tfdsp::ElectricPianoVoice fmVoice;
		tfdsp::ElectricPianoVoice fixedPitchVoice;
		fmVoice.SetSampleRate(48000.0);
		fixedPitchVoice.SetSampleRate(48000.0);
		fmVoice.SetNoiseSeed(0x13579bdu);
		fixedPitchVoice.SetNoiseSeed(0x13579bdu);
		double fmDifference = 0.0;
		bool fmFinite = true;
		for (int sample = 0; sample < 8192; ++sample)
		{
			const double pitchModulation = 0.02 * std::sin(
				2.0 * tfdsp::PI * 5.0 * sample / 48000.0);
			const double modulated = fmVoice.Step(pitchModulation, 10.0, 0.8,
				false, quietControls);
			const double fixed = fixedPitchVoice.Step(0.0, 10.0, 0.8, false,
				quietControls);
			fmDifference += (modulated - fixed) * (modulated - fixed);
			fmFinite = fmFinite && std::isfinite(modulated);
		}
		Check(fmFinite && fmDifference > 1.0e-6,
			"electric piano preserves continuous per-sample pitch modulation");

		tfdsp::ElectricPianoVoice editedVoice;
		tfdsp::ElectricPianoVoice unchangedVoice;
		editedVoice.SetSampleRate(48000.0);
		unchangedVoice.SetSampleRate(48000.0);
		double liveControlDifference = 0.0;
		for (int sample = 0; sample < 4096; ++sample)
		{
			auto editedControls = quietControls;
			if (sample >= 2048)
			{
				editedControls.body = 0.1;
				editedControls.bell = 0.9;
				editedControls.proximity = 0.9;
				editedControls.tone = 0.2;
			}
			const double edited = editedVoice.Step(0.0, 10.0, 0.8, false,
				editedControls);
			const double unchanged = unchangedVoice.Step(0.0, 10.0, 0.8,
				false, quietControls);
			if (sample >= 2048)
			{
				const double difference = edited - unchanged;
				liveControlDifference += difference * difference;
			}
		}
		Check(liveControlDifference > 1.0e-5,
			"electric piano cached timbre coefficients follow live controls");

		auto weakCouplingControls = quietControls;
		auto strongCouplingControls = quietControls;
		weakCouplingControls.coupling = 0.0;
		strongCouplingControls.coupling = 1.0;
		tfdsp::ElectricPianoVoice weakCouplingVoice;
		tfdsp::ElectricPianoVoice strongCouplingVoice;
		weakCouplingVoice.SetSampleRate(48000.0);
		strongCouplingVoice.SetSampleRate(48000.0);
		std::vector<double> weakCouplingSignal(48000);
		std::vector<double> strongCouplingSignal(48000);
		for (std::size_t sample = 0; sample < weakCouplingSignal.size(); ++sample)
		{
			weakCouplingSignal[sample] = weakCouplingVoice.Step(0.0, 10.0,
				0.85, false, weakCouplingControls);
			strongCouplingSignal[sample] = strongCouplingVoice.Step(0.0, 10.0,
				0.85, false, strongCouplingControls);
		}
		auto couplingDifference = [&](std::size_t begin, std::size_t end)
		{
			double weakEnergy = 0.0;
			double strongEnergy = 0.0;
			double cross = 0.0;
			for (std::size_t sample = begin; sample < end; ++sample)
			{
				weakEnergy += weakCouplingSignal[sample] *
					weakCouplingSignal[sample];
				strongEnergy += strongCouplingSignal[sample] *
					strongCouplingSignal[sample];
				cross += weakCouplingSignal[sample] *
					strongCouplingSignal[sample];
			}
			const double gain = cross / std::max(1.0e-20, strongEnergy);
			double residualEnergy = 0.0;
			for (std::size_t sample = begin; sample < end; ++sample)
			{
				const double residual = weakCouplingSignal[sample] -
					gain * strongCouplingSignal[sample];
				residualEnergy += residual * residual;
			}
			return std::array<double, 2>{
				std::sqrt(residualEnergy / std::max(1.0e-20, weakEnergy)),
				10.0 * std::log10(std::max(1.0e-20,
					strongEnergy / weakEnergy))};
		};
		const auto couplingAttack = couplingDifference(0, 2400);
		const auto couplingSustain = couplingDifference(12000, 48000);
		Check(couplingAttack[0] > 0.12 && couplingSustain[0] > 0.22 &&
			std::abs(couplingAttack[1]) < 1.0 && couplingSustain[1] > 3.0,
			"electric piano Coupling audibly changes support loss and two-body evolution without acting as attack gain");

		tfdsp::ElectricPianoVoice sweptCouplingVoice;
		tfdsp::ElectricPianoVoice fixedCouplingVoice;
		sweptCouplingVoice.SetSampleRate(48000.0);
		fixedCouplingVoice.SetSampleRate(48000.0);
		double sweptCouplingDifference = 0.0;
		double fixedCouplingEnergy = 0.0;
		double preSweepStep = 0.0;
		double duringSweepStep = 0.0;
		double previousSweptOutput = 0.0;
		for (int sample = 0; sample < 24000; ++sample)
		{
			auto sweptControls = weakCouplingControls;
			if (sample >= 4096)
				sweptControls.coupling = 1.0;
			const double sweptOutput = sweptCouplingVoice.Step(0.0, 10.0, 0.85,
				false, sweptControls);
			const double fixedOutput = fixedCouplingVoice.Step(0.0, 10.0, 0.85,
				false, weakCouplingControls);
			if (sample > 2048 && sample < 4096)
				preSweepStep = std::max(preSweepStep,
					std::abs(sweptOutput - previousSweptOutput));
			if (sample >= 4096 && sample < 6144)
				duringSweepStep = std::max(duringSweepStep,
					std::abs(sweptOutput - previousSweptOutput));
			if (sample >= 8192)
			{
				const double difference = sweptOutput - fixedOutput;
				sweptCouplingDifference += difference * difference;
				fixedCouplingEnergy += fixedOutput * fixedOutput;
			}
			previousSweptOutput = sweptOutput;
		}
		if (!(sweptCouplingDifference > 0.10 * fixedCouplingEnergy &&
			duringSweepStep < 1.5 * preSweepStep))
			std::cerr << "electric piano live Coupling difference ratio " <<
				sweptCouplingDifference / std::max(1.0e-20, fixedCouplingEnergy) <<
				" step ratio " << duringSweepStep /
				std::max(1.0e-20, preSweepStep) << '\n';
		Check(sweptCouplingDifference > 0.10 * fixedCouplingEnergy &&
			duringSweepStep < 1.5 * preSweepStep,
			"electric piano live Coupling sweep preserves state while audibly changing a held note");

		std::array<tfdsp::ElectricPianoVoice, 16> chordVoices;
		for (std::size_t voice = 0; voice < chordVoices.size(); ++voice)
		{
			chordVoices[voice].SetSampleRate(48000.0);
			chordVoices[voice].SetNoiseSeed(
				0x13579bdu + static_cast<std::uint32_t>(voice));
		}
		double chordPeak = 0.0;
		bool chordFinite = true;
		std::array<double, 16> keyboardEnergies{};
		for (int sample = 0; sample < 24000; ++sample)
		{
			double chord = 0.0;
			for (std::size_t voice = 0; voice < chordVoices.size(); ++voice)
			{
				const double pitch = -3.0 + 6.0 *
					static_cast<double>(voice) / 15.0;
				const double voiceOutput = chordVoices[voice].Step(pitch, 10.0,
					0.75, false, controls);
				chord += voiceOutput;
				keyboardEnergies[voice] += voiceOutput * voiceOutput;
			}
			chordFinite = chordFinite && std::isfinite(chord);
			chordPeak = std::max(chordPeak, std::abs(chord));
		}
		Check(chordFinite && chordPeak > 1.0e-4,
			"electric piano renders a finite audible 16-voice keyboard span");
		Check(*std::min_element(keyboardEnergies.begin(),
			keyboardEnergies.end()) > 1.0e-8,
			"electric piano keeps every tested key audible across six octaves");
		const double minimumKeyboardEnergy = *std::min_element(
			keyboardEnergies.begin(), keyboardEnergies.end());
		const double maximumKeyboardEnergy = *std::max_element(
			keyboardEnergies.begin(), keyboardEnergies.end());
		if (!(maximumKeyboardEnergy < 2.0 * minimumKeyboardEnergy))
		{
			const auto minimumIterator = std::min_element(keyboardEnergies.begin(),
				keyboardEnergies.end());
			const auto maximumIterator = std::max_element(keyboardEnergies.begin(),
				keyboardEnergies.end());
			std::cerr << "electric piano keyboard energy ratio: " <<
				maximumKeyboardEnergy / minimumKeyboardEnergy << " min/max voices " <<
				std::distance(keyboardEnergies.begin(), minimumIterator) << "/" <<
				std::distance(keyboardEnergies.begin(), maximumIterator) << '\n';
		}
		Check(maximumKeyboardEnergy < 2.0 * minimumKeyboardEnergy,
			"electric piano reference pickup curve keeps equal-velocity key levels balanced");

		auto softHammerControls = silentControls;
		auto hardHammerControls = silentControls;
		softHammerControls.hammer = 0.0;
		hardHammerControls.hammer = 1.0;
		tfdsp::ElectricPianoVoice softHammerVoice;
		tfdsp::ElectricPianoVoice defaultHammerVoice;
		tfdsp::ElectricPianoVoice hardHammerVoice;
		softHammerVoice.SetSampleRate(48000.0);
		defaultHammerVoice.SetSampleRate(48000.0);
		hardHammerVoice.SetSampleRate(48000.0);
		std::vector<double> softHammerSignal(4096);
		std::vector<double> defaultHammerSignal(4096);
		std::vector<double> hardHammerSignal(4096);
		for (std::size_t sample = 0; sample < softHammerSignal.size(); ++sample)
		{
			softHammerSignal[sample] = softHammerVoice.Step(0.0, 10.0, 0.72,
				false, softHammerControls);
			defaultHammerSignal[sample] = defaultHammerVoice.Step(0.0, 10.0,
				0.72, false, silentControls);
			hardHammerSignal[sample] = hardHammerVoice.Step(0.0, 10.0, 0.72,
				false, hardHammerControls);
		}
		const double fundamentalBin = tfdsp::ElectricPianoReferenceFrequency /
			48000.0;
		const double middleKeyboardPosition = (60.0 - 28.0) / 72.0;
		const double firstAttackRatio = 6.267 *
			(1.0 + (1.0 - 2.0 * middleKeyboardPosition) * 0.035);
		const double bellBin = firstAttackRatio * fundamentalBin;
		auto attackBrightness = [&](const std::vector<double>& signal)
		{
			return WindowedHarmonicMagnitude(signal, bellBin, 0, 2048) /
				std::max(1.0e-12, WindowedHarmonicMagnitude(signal,
					fundamentalBin, 0, 2048));
		};
		const double softHammerBrightness = attackBrightness(softHammerSignal);
		const double defaultHammerBrightness = attackBrightness(
			defaultHammerSignal);
		const double hardHammerBrightness = attackBrightness(hardHammerSignal);
		Check(hardHammerBrightness > 1.35 * softHammerBrightness,
			"electric piano hammer hardness changes compliant-contact brightness");
		Check(defaultHammerBrightness > 0.025 &&
			defaultHammerBrightness < 0.080 && hardHammerBrightness < 0.14,
			"electric piano attack mode has an audible non-metallic absolute level");

		double minimumKeyboardAttack = std::numeric_limits<double>::infinity();
		double maximumKeyboardAttack = 0.0;
		for (double pitch : {-2.0, -1.0, 0.0, 1.0})
		{
			const double keyboardPosition = std::clamp(
				(60.0 + 12.0 * pitch - 28.0) / 72.0, 0.0, 1.0);
			const double fundamental = tfdsp::ElectricPianoReferenceFrequency *
				std::exp2(pitch) / 48000.0;
			const double attack = 6.267 *
				(1.0 + (1.0 - 2.0 * keyboardPosition) * 0.035) * fundamental;
			tfdsp::ElectricPianoVoice voice;
			voice.SetSampleRate(48000.0);
			std::vector<double> signal(2048);
			for (double& output : signal)
				output = voice.Step(pitch, 10.0, 0.72, false, silentControls);
			const double brightness = HarmonicMagnitude(signal, attack) /
				std::max(1.0e-12, HarmonicMagnitude(signal, fundamental));
			minimumKeyboardAttack = std::min(minimumKeyboardAttack, brightness);
			maximumKeyboardAttack = std::max(maximumKeyboardAttack, brightness);
		}
		Check(minimumKeyboardAttack > 0.018 && maximumKeyboardAttack < 0.18,
			"electric piano contact-generated first attack mode survives across the keyboard");

		auto nonlinearControls = silentControls;
		tfdsp::ElectricPianoVoice lowVelocityPickupVoice;
		tfdsp::ElectricPianoVoice highVelocityPickupVoice;
		lowVelocityPickupVoice.SetSampleRate(48000.0);
		highVelocityPickupVoice.SetSampleRate(48000.0);
		std::vector<double> lowVelocityPickupSignal(16000);
		std::vector<double> highVelocityPickupSignal(16000);
		for (std::size_t sample = 0; sample < lowVelocityPickupSignal.size(); ++sample)
		{
			lowVelocityPickupSignal[sample] = lowVelocityPickupVoice.Step(0.0,
				10.0, 0.18, false, nonlinearControls);
			highVelocityPickupSignal[sample] = highVelocityPickupVoice.Step(0.0,
				10.0, 1.0, false, nonlinearControls);
		}
		auto lateHarmonicRatio = [&](const std::vector<double>& signal,
			double harmonic)
		{
			return WindowedHarmonicMagnitude(signal, harmonic * fundamentalBin,
				6000, signal.size()) / std::max(1.0e-12,
				WindowedHarmonicMagnitude(signal, fundamentalBin, 6000,
					signal.size()));
		};
		const double lowSecond = lateHarmonicRatio(lowVelocityPickupSignal, 2.0);
		const double highSecond = lateHarmonicRatio(highVelocityPickupSignal, 2.0);
		const double lowThird = lateHarmonicRatio(lowVelocityPickupSignal, 3.0);
		const double highThird = lateHarmonicRatio(highVelocityPickupSignal, 3.0);
		if (!(highSecond > 3.0 * lowSecond && highSecond > 0.055 &&
			highSecond < 0.18 && highThird > 8.0 * lowThird &&
			highThird > 0.0074 && highThird < 0.10))
			std::cerr << "electric piano pickup ratios low/high H2 " <<
				lowSecond << "/" << highSecond << " H3 " << lowThird <<
				"/" << highThird << '\n';
		Check(highSecond > 3.0 * lowSecond && highSecond > 0.055 &&
			highSecond < 0.18 && highThird > 8.0 * lowThird &&
			highThird > 0.0074 && highThird < 0.10,
			"electric piano asymmetric-pole pickup develops progressive bark with velocity");
		double lowUpperHarmonicEnergy = 0.0;
		double highUpperHarmonicEnergy = 0.0;
		for (int harmonic = 4; harmonic <= 10; ++harmonic)
		{
			const double lowMagnitude = lateHarmonicRatio(lowVelocityPickupSignal,
				static_cast<double>(harmonic));
			const double highMagnitude = lateHarmonicRatio(highVelocityPickupSignal,
				static_cast<double>(harmonic));
			lowUpperHarmonicEnergy += lowMagnitude * lowMagnitude;
			highUpperHarmonicEnergy += highMagnitude * highMagnitude;
		}
		const double lowUpperHarmonics = std::sqrt(lowUpperHarmonicEnergy);
		const double highUpperHarmonics = std::sqrt(highUpperHarmonicEnergy);
		if (!(highUpperHarmonics > 3.5 * lowUpperHarmonics &&
			highUpperHarmonics > 0.015))
			std::cerr << "electric piano pickup H4-H10 low/high " <<
				lowUpperHarmonics << "/" << highUpperHarmonics << '\n';
		Check(highUpperHarmonics > 3.5 * lowUpperHarmonics &&
			highUpperHarmonics > 0.015,
			"electric piano asymmetric pickup generates a velocity-dependent upper harmonic ladder");
		auto earlySidebandRatio = [&](const std::vector<double>& signal)
		{
			const double fundamental = WindowedHarmonicMagnitude(signal,
				fundamentalBin, 0, 2048);
			const double lowerSideband = WindowedHarmonicMagnitude(signal,
				bellBin - fundamentalBin, 0, 2048);
			const double upperSideband = WindowedHarmonicMagnitude(signal,
				bellBin + fundamentalBin, 0, 2048);
			return std::sqrt(lowerSideband * lowerSideband +
				upperSideband * upperSideband) / std::max(1.0e-12, fundamental);
		};
		const double lowVelocitySidebands = earlySidebandRatio(
			lowVelocityPickupSignal);
		const double highVelocitySidebands = earlySidebandRatio(
			highVelocityPickupSignal);
		if (!(highVelocitySidebands > 1.5 * lowVelocitySidebands &&
			highVelocitySidebands > 0.006))
			std::cerr << "electric piano pickup sidebands low/high " <<
				lowVelocitySidebands << "/" << highVelocitySidebands << '\n';
		Check(highVelocitySidebands > 1.5 * lowVelocitySidebands &&
			highVelocitySidebands > 0.006,
			"electric piano pickup intermodulates the fundamental and attack mode");

		// Exercise the velocities musicians actually play, not only the endpoints.
		// A coupled hammer changes both contact bandwidth and pickup excursion, so
		// attack-mode and low-harmonic ratios should rise throughout the range.
		constexpr std::array<double, 6> TimbreVelocities{
			0.25, 0.40, 0.55, 0.70, 0.85, 1.0};
		std::array<double, TimbreVelocities.size()> secondHarmonics{};
		std::array<double, TimbreVelocities.size()> thirdHarmonics{};
		std::array<double, TimbreVelocities.size()> attackRatios{};
		for (std::size_t velocityIndex = 0;
			velocityIndex < TimbreVelocities.size(); ++velocityIndex)
		{
			tfdsp::ElectricPianoVoice voice;
			voice.SetSampleRate(48000.0);
			std::vector<double> signal(16000);
			for (double& output : signal)
				output = voice.Step(0.0, 10.0,
					TimbreVelocities[velocityIndex], false, nonlinearControls);
			const double fundamental = WindowedHarmonicMagnitude(signal,
				fundamentalBin, 6000, signal.size());
			secondHarmonics[velocityIndex] = WindowedHarmonicMagnitude(signal,
				2.0 * fundamentalBin, 6000, signal.size()) /
				std::max(1.0e-12, fundamental);
			thirdHarmonics[velocityIndex] = WindowedHarmonicMagnitude(signal,
				3.0 * fundamentalBin, 6000, signal.size()) /
				std::max(1.0e-12, fundamental);
			attackRatios[velocityIndex] = WindowedHarmonicMagnitude(signal,
				bellBin, 0, 2048) / std::max(1.0e-12,
					WindowedHarmonicMagnitude(signal, fundamentalBin, 0, 2048));
		}
		Check(std::is_sorted(secondHarmonics.begin(), secondHarmonics.end()) &&
			secondHarmonics.back() > 2.5 * secondHarmonics.front() &&
			std::is_sorted(thirdHarmonics.begin(), thirdHarmonics.end()) &&
			thirdHarmonics.back() > 8.0 * thirdHarmonics.front() &&
			std::is_sorted(attackRatios.begin(), attackRatios.end()) &&
			attackRatios.back() > 2.0 * attackRatios.front(),
			"electric piano timbre grows continuously across playable velocities");

		auto darkToneControls = nonlinearControls;
		auto brightToneControls = nonlinearControls;
		darkToneControls.tone = 0.0;
		brightToneControls.tone = 1.0;
		tfdsp::ElectricPianoVoice darkToneVoice;
		tfdsp::ElectricPianoVoice brightToneVoice;
		darkToneVoice.SetSampleRate(48000.0);
		brightToneVoice.SetSampleRate(48000.0);
		std::vector<double> darkToneSignal(16000);
		std::vector<double> brightToneSignal(16000);
		for (std::size_t sample = 0; sample < darkToneSignal.size(); ++sample)
		{
			darkToneSignal[sample] = darkToneVoice.Step(0.0, 10.0, 0.85,
				false, darkToneControls);
			brightToneSignal[sample] = brightToneVoice.Step(0.0, 10.0, 0.85,
				false, brightToneControls);
		}
		const double darkFirst = HarmonicMagnitude(darkToneSignal, bellBin) /
			std::max(1.0e-12, HarmonicMagnitude(darkToneSignal, fundamentalBin));
		const double brightFirst = HarmonicMagnitude(brightToneSignal, bellBin) /
			std::max(1.0e-12, HarmonicMagnitude(brightToneSignal, fundamentalBin));
		const double darkFundamental = WindowedHarmonicMagnitude(darkToneSignal,
			fundamentalBin, 6000, darkToneSignal.size());
		const double brightFundamental = WindowedHarmonicMagnitude(
			brightToneSignal, fundamentalBin, 6000, brightToneSignal.size());
		const double darkSecond = WindowedHarmonicMagnitude(darkToneSignal,
			2.0 * fundamentalBin, 6000, darkToneSignal.size()) /
			std::max(1.0e-12, darkFundamental);
		const double brightSecond = WindowedHarmonicMagnitude(brightToneSignal,
			2.0 * fundamentalBin, 6000, brightToneSignal.size()) /
			std::max(1.0e-12, brightFundamental);
		const double darkThird = WindowedHarmonicMagnitude(darkToneSignal,
			3.0 * fundamentalBin, 6000, darkToneSignal.size()) /
			std::max(1.0e-12, darkFundamental);
		const double brightThird = WindowedHarmonicMagnitude(brightToneSignal,
			3.0 * fundamentalBin, 6000, brightToneSignal.size()) /
			std::max(1.0e-12, brightFundamental);
		if (!(brightFirst > 0.80 * darkFirst && brightFirst < 1.30 * darkFirst &&
			brightSecond > 8.0 * darkSecond && brightThird > 0.020 &&
			brightThird < 1.10 * darkThird))
			std::cerr << "electric piano Tone ratios attack dark/bright " <<
				darkFirst << "/" << brightFirst << " H2 " << darkSecond <<
				"/" << brightSecond << " H3 " << darkThird << "/" <<
				brightThird << '\n';
		Check(brightFirst > 0.80 * darkFirst && brightFirst < 1.30 * darkFirst &&
			brightSecond > 8.0 * darkSecond &&
			brightThird > 0.020 && brightThird < 1.10 * darkThird,
			"electric piano Tone changes physical pickup alignment and harmonic balance");

		auto amplify = [&](const std::vector<double>& direct)
		{
			tfdsp::ElectricPianoAmplifier testAmplifier;
			testAmplifier.SetSampleRate(48000.0);
			std::vector<double> amplified(direct.size());
			for (std::size_t sample = 0; sample < direct.size(); ++sample)
				amplified[sample] = testAmplifier.Step(5.0 * direct[sample],
					controls)[0];
			return amplified;
		};
		auto levelMatchedResidual = [](const std::vector<double>& first,
			const std::vector<double>& second, std::size_t end)
		{
			end = std::min({end, first.size(), second.size()});
			double firstEnergy = 0.0;
			double secondEnergy = 0.0;
			double cross = 0.0;
			for (std::size_t sample = 0; sample < end; ++sample)
			{
				firstEnergy += first[sample] * first[sample];
				secondEnergy += second[sample] * second[sample];
				cross += first[sample] * second[sample];
			}
			const double gain = cross / std::max(1.0e-20, secondEnergy);
			double residualEnergy = 0.0;
			for (std::size_t sample = 0; sample < end; ++sample)
			{
				const double residual = first[sample] - gain * second[sample];
				residualEnergy += residual * residual;
			}
			return std::sqrt(residualEnergy /
				std::max(1.0e-20, firstEnergy));
		};
		const auto darkToneAmplified = amplify(darkToneSignal);
		const auto brightToneAmplified = amplify(brightToneSignal);
		const auto softHammerAmplified = amplify(softHammerSignal);
		const auto hardHammerAmplified = amplify(hardHammerSignal);
		const double amplifiedToneResidual = levelMatchedResidual(
			darkToneAmplified, brightToneAmplified, 2048);
		if (!(amplifiedToneResidual > 0.32))
			std::cerr << "electric piano amplified Tone residual " <<
				amplifiedToneResidual << '\n';
		Check(amplifiedToneResidual > 0.32,
			"electric piano Tone extremes remain distinct through the shared amplifier");
		Check(levelMatchedResidual(softHammerAmplified, hardHammerAmplified,
			2048) > 0.48,
			"electric piano Hammer extremes remain distinct through the shared amplifier");

		// TONE is a live pickup adjustment. Verify an isolated full-range sweep on
		// an already sounding note, independent of the note-latched HAMMER state.
		tfdsp::ElectricPianoVoice sweptToneVoice;
		tfdsp::ElectricPianoVoice fixedToneVoice;
		sweptToneVoice.SetSampleRate(48000.0);
		fixedToneVoice.SetSampleRate(48000.0);
		auto sweptControls = silentControls;
		auto fixedControls = silentControls;
		sweptControls.tone = 0.0;
		fixedControls.tone = 0.0;
		for (int sample = 0; sample < 4096; ++sample)
		{
			sweptToneVoice.Step(0.0, 10.0, 0.8, false, sweptControls);
			fixedToneVoice.Step(0.0, 10.0, 0.8, false, fixedControls);
		}
		sweptControls.tone = 1.0;
		double sweptDifferenceEnergy = 0.0;
		double fixedToneEnergy = 0.0;
		for (int sample = 0; sample < 4096; ++sample)
		{
			const double swept = sweptToneVoice.Step(0.0, 10.0, 0.8, false,
				sweptControls);
			const double fixed = fixedToneVoice.Step(0.0, 10.0, 0.8, false,
				fixedControls);
			if (sample >= 1024)
			{
				const double difference = swept - fixed;
				sweptDifferenceEnergy += difference * difference;
				fixedToneEnergy += fixed * fixed;
			}
		}
		if (!(sweptDifferenceEnergy > 0.10 * fixedToneEnergy))
			std::cerr << "electric piano active Tone sweep energy ratio " <<
				sweptDifferenceEnergy / std::max(1.0e-20, fixedToneEnergy) << '\n';
		Check(sweptDifferenceEnergy > 0.10 * fixedToneEnergy,
			"electric piano Tone sweep audibly changes an active note");

		tfdsp::ElectricPianoVoice mechanicalVoice;
		tfdsp::ElectricPianoVoice noMechanicalVoice;
		mechanicalVoice.SetSampleRate(48000.0);
		noMechanicalVoice.SetSampleRate(48000.0);
		mechanicalVoice.SetNoiseSeed(0x77112233u);
		noMechanicalVoice.SetNoiseSeed(0x77112233u);
		double dryVoiceEnergy = 0.0;
		double mechanicalDifferenceEnergy = 0.0;
		for (int sample = 0; sample < 4800; ++sample)
		{
			const double withMechanics = mechanicalVoice.Step(0.0, 10.0, 0.72,
				false, controls);
			const double withoutMechanics = noMechanicalVoice.Step(0.0, 10.0,
				0.72, false, silentControls);
			dryVoiceEnergy += withoutMechanics * withoutMechanics;
			const double difference = withMechanics - withoutMechanics;
			mechanicalDifferenceEnergy += difference * difference;
		}
		Check(mechanicalDifferenceEnergy < 0.025 * dryVoiceEnergy,
			"electric piano default mechanics remain subordinate to the pitched voice");

		// The factory decay setting is the unscaled physical calibration. Longer
		// settings must preserve modal energy, rather than acting as an output
		// envelope or changing the damper release.
		Check(std::abs(tfdsp::ElectricPianoControls{}.decay - 0.5) < 1.0e-12,
			"electric piano Decay defaults to the neutral modal-lifetime scale");
		auto shortDecayControls = silentControls;
		auto longDecayControls = silentControls;
		shortDecayControls.decay = 0.0;
		longDecayControls.decay = 1.0;
		tfdsp::ElectricPianoVoice shortDecayVoice;
		tfdsp::ElectricPianoVoice longDecayVoice;
		shortDecayVoice.SetSampleRate(48000.0);
		longDecayVoice.SetSampleRate(48000.0);
		for (int sample = 0; sample < 24000; ++sample)
		{
			shortDecayVoice.Step(0.0, 10.0, 0.8,
				false, shortDecayControls);
			longDecayVoice.Step(0.0, 10.0, 0.8,
				false, longDecayControls);
		}
		Check(longDecayVoice.Energy() > 2.0 * shortDecayVoice.Energy(),
			"electric piano Decay scales the natural loss of sounding modes");

		// Independent random carriers used to put each note's mechanical peak at
		// an arbitrary time. The resonant impact is now event-synchronous; seeds
		// alter only the small noise texture.
		std::array<double, 6> mechanicsCentroids{};
		std::array<double, 6> mechanicsPeaks{};
		for (std::size_t event = 0; event < mechanicsCentroids.size(); ++event)
		{
			auto wetControls = controls;
			auto dryControls = controls;
			wetControls.mechanics = 1.0;
			dryControls.mechanics = 0.0;
			tfdsp::ElectricPianoVoice wetVoice;
			tfdsp::ElectricPianoVoice dryVoice;
			const auto seed = 0x1234567u +
				static_cast<std::uint32_t>(7919 * event);
			wetVoice.SetNoiseSeed(seed);
			dryVoice.SetNoiseSeed(seed);
			wetVoice.SetSampleRate(48000.0);
			dryVoice.SetSampleRate(48000.0);
			double weightedEnergy = 0.0;
			double totalEnergy = 0.0;
			for (int sample = 0; sample < 240; ++sample)
			{
				const double mechanical = wetVoice.Step(0.0, 10.0, 0.8,
					false, wetControls) - dryVoice.Step(0.0, 10.0, 0.8,
					false, dryControls);
				const double eventEnergy = mechanical * mechanical;
				weightedEnergy += static_cast<double>(sample) * eventEnergy;
				totalEnergy += eventEnergy;
				mechanicsPeaks[event] = std::max(mechanicsPeaks[event],
					std::abs(mechanical));
			}
			mechanicsCentroids[event] = weightedEnergy /
				std::max(1.0e-20, totalEnergy);
		}
		const auto mechanicsCentroidRange = std::minmax_element(
			mechanicsCentroids.begin(), mechanicsCentroids.end());
		const auto mechanicsPeakRange = std::minmax_element(
			mechanicsPeaks.begin(), mechanicsPeaks.end());
		Check(*mechanicsCentroidRange.second - *mechanicsCentroidRange.first < 4.0 &&
			*mechanicsPeakRange.second < 1.10 * *mechanicsPeakRange.first,
			"electric piano polyphonic mechanics start consistently at each key event");

		tfdsp::ElectricPianoVoice releasedVoice;
		releasedVoice.SetSampleRate(48000.0);
		double earlyReleaseEnergy = 0.0;
		double lateReleaseEnergy = 0.0;
		for (int sample = 0; sample < 48000; ++sample)
		{
			const double gate = sample < 2400 ? 10.0 : 0.0;
			const double output = releasedVoice.Step(-1.0, gate, 0.8, false,
				controls);
			if (sample >= 2400 && sample < 4800)
				earlyReleaseEnergy += output * output;
			if (sample >= 43200)
				lateReleaseEnergy += output * output;
		}
		Check(lateReleaseEnergy < 1.0e-4 * earlyReleaseEnergy,
			"electric piano damper silences a released note");

		tfdsp::ElectricPianoVoice sustainedVoice;
		tfdsp::ElectricPianoVoice dampedVoice;
		sustainedVoice.SetSampleRate(48000.0);
		dampedVoice.SetSampleRate(48000.0);
		for (int sample = 0; sample < 18000; ++sample)
		{
			const double gate = sample < 2400 ? 10.0 : 0.0;
			sustainedVoice.Step(-1.0, gate, 0.8, sample >= 2400,
				silentControls);
			dampedVoice.Step(-1.0, gate, 0.8, false, silentControls);
		}
		Check(sustainedVoice.Energy() > 100.0 * dampedVoice.Energy(),
			"electric piano sustain pedal prevents damper contact after key release");
		const double heldPedalEnergy = sustainedVoice.Energy();
		for (int sample = 0; sample < 24000; ++sample)
			sustainedVoice.Step(-1.0, 0.0, 0.0, false, silentControls);
		Check(sustainedVoice.Energy() < 1.0e-3 * heldPedalEnergy,
			"electric piano pedal release applies the damper to sustained notes");

		tfdsp::ElectricPianoAmplifier amplifier;
		amplifier.SetSampleRate(48000.0);
		auto maximumAmplifierControls = controls;
		maximumAmplifierControls.drive = 1.0;
		maximumAmplifierControls.outputVolume = 1.0;
		double amplifierPeak = 0.0;
		for (int sample = 0; sample < 48000; ++sample)
		{
			const double input = 16.0 * std::sin(
				2.0 * tfdsp::PI * 220.0 * sample / 48000.0);
			const auto output = amplifier.Step(input, maximumAmplifierControls);
			amplifierPeak = std::max({amplifierPeak,
				std::abs(output[0]), std::abs(output[1])});
		}
		Check(std::isfinite(amplifierPeak) && amplifierPeak < 2.0,
			"electric piano shared amplifier overload remains bounded");

		// Direct Figure 11-8 input-pair check against the durable ngspice deck.
		// This prevents the preamp reduction from drifting into an arbitrary
		// waveshaper while combined Drive/power-stage tests still happen to pass.
		struct PreampProbe
		{
			double rms{};
			double thd{};
			double gain{};
			double averageIterations{};
			std::uint64_t failures{};
		};
		auto renderPreamp = [](double amplitude)
		{
			constexpr double rate = 96000.0;
			constexpr double frequency = 250.0;
			tfdsp::PetersonPreAmplifier preamp;
			preamp.SetSampleRate(rate);
			std::vector<double> signal(48000);
			double energy = 0.0;
			for (int sample = 0; sample < 96000; ++sample)
			{
				const double input = amplitude * std::sin(2.0 * tfdsp::PI *
					frequency * sample / rate);
				const double output = preamp.Step(input).voltage;
				if (sample >= 48000)
				{
					signal[sample - 48000] = output;
					energy += output * output;
				}
			}
			const double bin = frequency / rate;
			const double fundamental = HarmonicMagnitude(signal, bin);
			double harmonicEnergy = 0.0;
			for (int harmonic = 2; harmonic <= 7; ++harmonic)
			{
				const double magnitude = HarmonicMagnitude(signal,
					harmonic * bin);
				harmonicEnergy += magnitude * magnitude;
			}
			const double rms = std::sqrt(energy /
				static_cast<double>(signal.size()));
			return PreampProbe{rms,
				std::sqrt(harmonicEnergy) / std::max(1.0e-20, fundamental),
				rms / (amplitude / std::sqrt(2.0)),
				preamp.AverageIterations(), preamp.SolverFailures()};
		};
		const auto preampReference = renderPreamp(0.20);
		const auto preampDriven = renderPreamp(1.50);
		const auto preampExtreme = renderPreamp(10.0);
		Check(std::abs(preampReference.rms - 1.06733) < 0.025 &&
			preampReference.thd > 0.0005 && preampReference.thd < 0.0011,
			"electric piano Figure 11-8 preamp matches clean ngspice level and THD");
		Check(std::abs(preampDriven.rms - 7.68424) < 0.16 &&
			preampDriven.thd > 0.045 && preampDriven.thd < 0.060 &&
			preampDriven.gain < 0.98 * preampReference.gain,
			"electric piano Figure 11-8 preamp develops asymmetric overload and gain loss");
		Check(preampReference.failures == 0 && preampDriven.failures == 0 &&
			preampReference.averageIterations < 1.1 &&
			preampDriven.averageIterations < 1.1,
			"electric piano preamp circuit converges on its audio predictor");
		Check(preampExtreme.failures == 0 && std::isfinite(preampExtreme.rms) &&
			preampExtreme.rms < 12.0 && preampExtreme.averageIterations < 2.5,
			"electric piano preamp remains convergent beyond calibrated musical input");

		// The feedback tone amplifier must remain a circuit, not regress to the
		// former one-pole EQ followed by an unrelated transfer curve. These targets
		// come from tests/python/reference_peterson_preamp_spice.py, whose netlist
		// contains the two Darlington-connected 2N3392 devices, Q5 and the actual
		// bass/treble bridge from Figure 11-8.
		struct CompletePreampProbe
		{
			double rms{};
			double thd{};
			double averageToneIterations{};
			std::uint64_t toneFailures{};
		};
		auto renderCompletePreamp = [](double amplitude)
		{
			constexpr double rate = 96000.0;
			constexpr double frequency = 250.0;
			tfdsp::PetersonPreAmplifier inputStage;
			tfdsp::PetersonTonePreAmplifier toneStage;
			inputStage.SetSampleRate(rate);
			toneStage.SetSampleRate(rate);
			std::vector<double> signal(48000);
			double energy = 0.0;
			for (int sample = 0; sample < 96000; ++sample)
			{
				const double input = amplitude * std::sin(2.0 * tfdsp::PI *
					frequency * sample / rate);
				const double firstStage = inputStage.Step(input).voltage;
				const double output = toneStage.Step(firstStage, 0.5, 0.5).voltage;
				if (sample >= 48000)
				{
					signal[sample - 48000] = output;
					energy += output * output;
				}
			}
			const double bin = frequency / rate;
			const double fundamental = HarmonicMagnitude(signal, bin);
			double harmonicEnergy = 0.0;
			for (int harmonic = 2; harmonic <= 7; ++harmonic)
			{
				const double magnitude = HarmonicMagnitude(signal, harmonic * bin);
				harmonicEnergy += magnitude * magnitude;
			}
			return CompletePreampProbe{
				std::sqrt(energy / static_cast<double>(signal.size())),
				std::sqrt(harmonicEnergy) / std::max(1.0e-20, fundamental),
				toneStage.AverageIterations(), toneStage.SolverFailures()};
		};
		const auto completePreampSmall = renderCompletePreamp(0.001);
		const auto completePreampMedium = renderCompletePreamp(0.50);
		const auto completePreampDriven = renderCompletePreamp(1.00);
		Check(std::abs(completePreampSmall.rms - 0.00442985) < 0.00010 &&
			std::abs(completePreampMedium.rms - 2.20510) < 0.050 &&
			std::abs(completePreampDriven.rms - 4.22895) < 0.10,
			"electric piano complete Figure 11-8 preamp matches ngspice levels");
		Check(completePreampSmall.thd < 0.00015 &&
			completePreampDriven.thd > 0.045 && completePreampDriven.thd < 0.065 &&
			completePreampDriven.rms < 0.97 * 1000.0 * completePreampSmall.rms,
			"electric piano active tone-feedback transistors create measured overload");
		Check(completePreampSmall.toneFailures == 0 &&
			completePreampMedium.toneFailures == 0 &&
			completePreampDriven.toneFailures == 0 &&
			completePreampDriven.averageToneIterations < 2.2,
			"electric piano nonlinear tone-feedback solve remains convergent");

		// Validate the Figure 11-9 real-time circuit directly against the checked-in
		// ngspice reference, before Drive compensation or tone controls can obscure
		// an electrical regression. Reference targets at 250 Hz are generated by
		// tests/python/reference_peterson_power_spice.py from the same schematic.
		auto renderPowerCircuit = [](double amplitude)
		{
			constexpr double sampleRate = 192000.0;
			constexpr double frequency = 250.0;
			constexpr int warmup = 192000;
			tfdsp::PetersonPowerAmplifier circuit;
			circuit.SetSampleRate(sampleRate);
			std::vector<double> signal(192000);
			for (int sample = 0; sample < warmup +
				static_cast<int>(signal.size()); ++sample)
			{
				const auto result = circuit.Step(amplitude * std::sin(
					2.0 * tfdsp::PI * frequency * sample / sampleRate), 35.0);
				if (sample >= warmup)
					signal[sample - warmup] = result.voltage;
			}
			return std::pair{std::move(signal), circuit.SolverFailures()};
		};
		const auto [referenceLevelSignal, referenceLevelFailures] =
			renderPowerCircuit(0.05);
		const auto [referenceOverloadSignal, referenceOverloadFailures] =
			renderPowerCircuit(0.50);
		auto acRms = [](const std::vector<double>& signal)
		{
			double mean = 0.0;
			for (double sample : signal)
				mean += sample;
			mean /= static_cast<double>(signal.size());
			double energy = 0.0;
			for (double sample : signal)
				energy += (sample - mean) * (sample - mean);
			return std::sqrt(energy / static_cast<double>(signal.size()));
		};
		const double referenceLevelRms = acRms(referenceLevelSignal);
		const double referenceOverloadRms = acRms(referenceOverloadSignal);
		const double powerBin = 250.0 / 192000.0;
		const double referenceLevelFundamental = HarmonicMagnitude(
			referenceLevelSignal, powerBin);
		const double referenceLevelSecond = HarmonicMagnitude(
			referenceLevelSignal, 2.0 * powerBin) / referenceLevelFundamental;
		const double referenceLevelThird = HarmonicMagnitude(
			referenceLevelSignal, 3.0 * powerBin) / referenceLevelFundamental;
		const double referenceLevelFourth = HarmonicMagnitude(
			referenceLevelSignal, 4.0 * powerBin) / referenceLevelFundamental;
		const double overloadFundamental = HarmonicMagnitude(
			referenceOverloadSignal, powerBin);
		const double overloadSecond = HarmonicMagnitude(referenceOverloadSignal,
			2.0 * powerBin) / overloadFundamental;
		const double overloadThird = HarmonicMagnitude(referenceOverloadSignal,
			3.0 * powerBin) / overloadFundamental;
		Check(referenceLevelFailures == 0 && referenceOverloadFailures == 0,
			"Peterson Figure 11-9 Newton solve converges through rail overload");
		Check(referenceLevelRms > 1.84 && referenceLevelRms < 2.05 &&
			referenceOverloadRms > 17.3 && referenceOverloadRms < 19.3,
			"Peterson Figure 11-9 levels agree with the ngspice reference");
		Check(referenceLevelSecond > 0.0012 && referenceLevelSecond < 0.0018 &&
			referenceLevelThird > 0.00075 && referenceLevelThird < 0.00125 &&
			referenceLevelFourth > 0.00055 && referenceLevelFourth < 0.00095,
			"Peterson Figure 11-9 clean crossover harmonics agree with ngspice");
		Check(overloadSecond > 0.020 && overloadSecond < 0.060 &&
			overloadThird > 0.018 && overloadThird < 0.050 &&
			overloadSecond > overloadThird,
			"Peterson matched-PNP asymmetric topology agrees with SPICE harmonics");

		auto renderAmplifierSine = [](double drive, double amplitude)
		{
			constexpr double frequency = 250.0;
			constexpr int warmup = 24000;
			tfdsp::ElectricPianoAmplifier testAmplifier;
			testAmplifier.SetSampleRate(48000.0);
			tfdsp::ElectricPianoControls sineControls;
			sineControls.drive = drive;
			std::vector<double> signal(48000);
			for (int sample = 0; sample < warmup +
				static_cast<int>(signal.size()); ++sample)
			{
				const double input = amplitude * std::sin(
					2.0 * tfdsp::PI * frequency * sample / 48000.0);
				const double output = testAmplifier.Step(input, sineControls)[0];
				if (sample >= warmup)
					signal[sample - warmup] = output;
			}
			return signal;
		};
		auto signalRms = [](const std::vector<double>& signal)
		{
			double energy = 0.0;
			for (double sample : signal)
				energy += sample * sample;
			return std::sqrt(energy / static_cast<double>(signal.size()));
		};
		const auto cleanProbe = renderAmplifierSine(0.0, 0.05);
		const auto drivenProbe = renderAmplifierSine(1.0, 0.05);
		const double driveLevelSpan = std::abs(20.0 * std::log10(
			signalRms(drivenProbe) / signalRms(cleanProbe)));
		if (!(driveLevelSpan < 3.0))
			std::cerr << "electric piano Drive small-signal level span " <<
				driveLevelSpan << " dB\n";
		Check(driveLevelSpan < 3.0,
			"electric piano Drive cancels knob gain while preserving circuit compression");

		// Keep the circuit/Rack-domain boundary explicit. The measured model peak
		// of a hard middle-C note is about 0.0404, and the module feeds five times
		// that value to the amplifier. Drive raises the voltage into Figure 11-8,
		// while its reciprocal at the schematic's volume node keeps Figure 11-9's
		// ideal small-signal demand independent of the knob.
		double renderedHardNotePeak = 0.0;
		for (double sample : highVelocityPickupSignal)
			renderedHardNotePeak = std::max(renderedHardNotePeak,
				std::abs(sample));
		Check(renderedHardNotePeak > 0.038 && renderedHardNotePeak < 0.043,
			"electric piano hard-note pickup peak remains tied to amplifier voltage calibration");
		const double hardNoteAmplifierPeak = 5.0 * renderedHardNotePeak;
		const double hardNoteCircuitInput = tfdsp::ElectricPianoAmplifier::
			ModelInputToCircuitVolts(hardNoteAmplifierPeak, 1.0);
		const double maximumDriveGain = tfdsp::ElectricPianoAmplifier::
			DriveGainForPosition(1.0);
		const double idealHardNotePowerPeak = 57.0 * hardNoteCircuitInput /
			maximumDriveGain;
		const double stressCircuitInput = tfdsp::ElectricPianoAmplifier::
			ModelInputToCircuitVolts(0.315, 1.0);
		const double idealStressPowerPeak = 57.0 * stressCircuitInput /
			maximumDriveGain;
		Check(hardNoteCircuitInput > 0.40 && hardNoteCircuitInput < 0.43 &&
			idealHardNotePowerPeak > 1.40 && idealHardNotePowerPeak < 1.60 &&
			stressCircuitInput > 0.63 && stressCircuitInput < 0.67 &&
			idealStressPowerPeak > 2.25 && idealStressPowerPeak < 2.45,
			"electric piano Drive exercises the preamp without raising ideal power-stage demand");
		for (double drive : {0.0, 0.25, 0.50, 0.75, 1.0})
		{
			const double circuitInput = tfdsp::ElectricPianoAmplifier::
				ModelInputToCircuitVolts(0.20, drive);
			const double powerInput = circuitInput /
				tfdsp::ElectricPianoAmplifier::DriveGainForPosition(drive);
			const double compensated = tfdsp::ElectricPianoAmplifier::
				CircuitPowerOutputToModel(57.0 * powerInput);
			Check(std::abs(compensated - 1.12 * 0.20) < 1.0e-12,
				"electric piano inverse Drive gain is applied at the pre-power volume node");
		}

		// 0.315 is a calibrated near-pickup/chord stress level above the roughly
		// 0.20 peak of one hard default-Proximity note. The old arbitrary amplitude
		// of 1.0 made an unrealistically large signal define normal playing.
		const auto cleanOverload = renderAmplifierSine(0.0, 0.315);
		const auto defaultOverload = renderAmplifierSine(0.32, 0.315);
		const auto drivenOverload = renderAmplifierSine(1.0, 0.315);
		const double amplifierBin = 250.0 / 48000.0;
		auto distortion = [&](const std::vector<double>& signal)
		{
			const double fundamental = HarmonicMagnitude(signal, amplifierBin);
			double sum = 0.0;
			for (int harmonic = 2; harmonic <= 7; ++harmonic)
			{
				const double magnitude = HarmonicMagnitude(signal,
					harmonic * amplifierBin);
				sum += magnitude * magnitude;
			}
			return std::sqrt(sum) / std::max(1.0e-12, fundamental);
		};
		const double cleanDistortion = distortion(cleanOverload);
		const double defaultDistortion = distortion(defaultOverload);
		const double drivenDistortion = distortion(drivenOverload);
		const double drivenThird = HarmonicMagnitude(drivenOverload,
			3.0 * amplifierBin);
		const double drivenFifth = HarmonicMagnitude(drivenOverload,
			5.0 * amplifierBin);
		const double drivenSeventh = HarmonicMagnitude(drivenOverload,
			7.0 * amplifierBin);
		if (!(cleanDistortion < 0.007 && defaultDistortion < 0.007 &&
			drivenDistortion > 2.0 * defaultDistortion &&
			drivenDistortion > 0.025 && drivenDistortion < 0.055))
			std::cerr << "electric piano THD clean/default/driven " <<
				cleanDistortion << "/" << defaultDistortion << "/" <<
				drivenDistortion << '\n';
		Check(cleanDistortion < 0.007 && defaultDistortion < 0.007 &&
			drivenDistortion > 2.0 * defaultDistortion &&
			drivenDistortion > 0.025 && drivenDistortion < 0.055,
			"electric piano default stays clean while maximum Drive adds controlled preamp harmonics");
		if (!(drivenFifth < 3.0 * drivenThird &&
			drivenSeventh < drivenThird))
			std::cerr << "electric piano driven harmonics H3/H5/H7 " <<
				drivenThird << "/" << drivenFifth << "/" << drivenSeventh <<
				" THD clean/default/driven " << cleanDistortion << "/" <<
				defaultDistortion << "/" << drivenDistortion << '\n';
		Check(drivenFifth < 3.0 * drivenThird && drivenSeventh < drivenThird,
			"electric piano circuit overload remains concentrated in low harmonics");

		// Guard the complete Drive trajectory with the calibrated 0.315 stress.
		// Real simultaneous hard-note renders are substantially larger than this;
		// keeping the proxy clean through 80% prevents moderate chord settings from
		// reaching both shared preamp headroom limits at once, while the final 20%
		// deliberately retains an overdriven circuit range.
		const auto drive20 = renderAmplifierSine(0.20, 0.315);
		const auto drive40 = renderAmplifierSine(0.40, 0.315);
		const auto drive60 = renderAmplifierSine(0.60, 0.315);
		const auto drive80 = renderAmplifierSine(0.80, 0.315);
		const auto drive100 = renderAmplifierSine(1.00, 0.315);
		const std::array<double, 5> driveTrajectory{
			distortion(drive20), distortion(drive40), distortion(drive60),
			distortion(drive80), distortion(drive100)};
		const bool progressiveDrive = driveTrajectory[0] < 0.008 &&
			driveTrajectory[1] < 0.008 && driveTrajectory[2] < 0.008 &&
			driveTrajectory[3] < 0.008 &&
			driveTrajectory[4] > 10.0 * driveTrajectory[3];
		if (!progressiveDrive)
			std::cerr << "electric piano Drive THD trajectory " <<
				driveTrajectory[0] << "/" << driveTrajectory[1] << "/" <<
				driveTrajectory[2] << "/" << driveTrajectory[3] << "/" <<
				driveTrajectory[4] << '\n';
		Check(progressiveDrive &&
			driveTrajectory.back() > 2.0 * driveTrajectory.front() &&
			driveTrajectory.back() > 0.020 &&
			driveTrajectory.back() < 0.040,
			"electric piano Drive reserves hard shared-preamp overload for its end stop");
		const double musicalSecond = HarmonicMagnitude(drive100,
			2.0 * amplifierBin);
		const double musicalThird = HarmonicMagnitude(drive100,
			3.0 * amplifierBin);
		Check(musicalSecond > 0.50 * musicalThird,
			"electric piano driven circuit retains its asymmetric even harmonic");

		// The service drawing only specifies a nominal 16 Ohm load. The current
		// provisional electromechanical equivalent is tested explicitly so future
		// Rhodes impedance measurements can replace these targets without silently
		// changing power-stage current, feedback or supply droop.
		const auto loadDc = tfdsp::ElectricPianoAmplifier::
			ElectricalLoadImpedance(0.0);
		const auto loadResonance = tfdsp::ElectricPianoAmplifier::
			ElectricalLoadImpedance(75.0);
		const auto loadMidband = tfdsp::ElectricPianoAmplifier::
			ElectricalLoadImpedance(1000.0);
		const auto loadHigh = tfdsp::ElectricPianoAmplifier::
			ElectricalLoadImpedance(10000.0);
		Check(std::abs(loadDc.real() - 12.8) < 0.05 &&
			std::abs(loadResonance) > 4.0 * std::abs(loadDc) &&
			std::abs(loadMidband) < 1.5 * std::abs(loadDc) &&
			std::abs(loadHigh) > 2.0 * std::abs(loadDc) &&
			loadHigh.imag() > 0.0,
			"electric piano reactive load has DC resistance, motional resonance and inductive rise");

		tfdsp::ElectricPianoAmplifier cleanNoteAmplifier;
		tfdsp::ElectricPianoAmplifier drivenNoteAmplifier;
		cleanNoteAmplifier.SetSampleRate(48000.0);
		drivenNoteAmplifier.SetSampleRate(48000.0);
		auto cleanNoteControls = controls;
		auto drivenNoteControls = controls;
		cleanNoteControls.drive = 0.0;
		drivenNoteControls.drive = 1.0;
		std::vector<double> cleanNote(highVelocityPickupSignal.size());
		std::vector<double> drivenNote(highVelocityPickupSignal.size());
		for (std::size_t sample = 0; sample < highVelocityPickupSignal.size(); ++sample)
		{
			cleanNote[sample] = cleanNoteAmplifier.Step(
				5.0 * highVelocityPickupSignal[sample], cleanNoteControls)[0];
			drivenNote[sample] = drivenNoteAmplifier.Step(
				5.0 * highVelocityPickupSignal[sample], drivenNoteControls)[0];
		}
		const double drivenNoteResidual = levelMatchedResidual(cleanNote,
			drivenNote, 4800);
		if (!(drivenNoteResidual > 0.002))
			std::cerr << "electric piano full-note Drive residual " <<
				drivenNoteResidual << '\n';
		Check(drivenNoteResidual > 0.002,
			"electric piano full-range Drive changes a rendered hard note");

		// Exercise the shared electronics with an actual eight-note render. A sine
		// proxy cannot reveal the crest collapse and intermodulation caused when the
		// summed keyboard signal is incorrectly compensated only after Figure 11-9.
		constexpr std::array<double, 8> chordPitches{
			-1.0, -10.0 / 12.0, -8.0 / 12.0, -5.0 / 12.0,
			-3.0 / 12.0, 0.0, 2.0 / 12.0, 4.0 / 12.0};
		std::array<tfdsp::ElectricPianoVoice, chordPitches.size()>
			amplifierChordVoices;
		for (std::size_t voice = 0; voice < amplifierChordVoices.size(); ++voice)
		{
			amplifierChordVoices[voice].SetSampleRate(48000.0);
			amplifierChordVoices[voice].SetNoiseSeed(0x4f1bbcdcu +
				static_cast<std::uint32_t>(voice) * 0x9e3779b9u);
		}
		tfdsp::ElectricPianoAmplifier cleanChordAmplifier;
		tfdsp::ElectricPianoAmplifier drivenChordAmplifier;
		cleanChordAmplifier.SetSampleRate(48000.0);
		drivenChordAmplifier.SetSampleRate(48000.0);
		auto chordVoiceControls = controls;
		chordVoiceControls.mechanics = 0.0;
		auto cleanChordControls = controls;
		auto drivenChordControls = controls;
		cleanChordControls.drive = 0.0;
		drivenChordControls.drive = 1.0;
		double cleanChordEnergy = 0.0;
		double drivenChordEnergy = 0.0;
		double cleanChordPeak = 0.0;
		double drivenChordPeak = 0.0;
		double drivenChordRailMinimum = 35.0;
		std::size_t chordSamples = 0;
		for (int sample = 0; sample < 48000; ++sample)
		{
			double pickupSum = 0.0;
			for (std::size_t voice = 0; voice < amplifierChordVoices.size(); ++voice)
				pickupSum += amplifierChordVoices[voice].Step(chordPitches[voice], 10.0,
					1.0, false, chordVoiceControls);
			const double clean = cleanChordAmplifier.Step(5.0 * pickupSum,
				cleanChordControls)[0];
			const double driven = drivenChordAmplifier.Step(5.0 * pickupSum,
				drivenChordControls)[0];
			drivenChordRailMinimum = std::min(drivenChordRailMinimum,
				drivenChordAmplifier.SupplyRailVoltage());
			if (sample >= 256 && sample < 24000)
			{
				cleanChordPeak = std::max(cleanChordPeak, std::abs(clean));
				drivenChordPeak = std::max(drivenChordPeak, std::abs(driven));
				cleanChordEnergy += clean * clean;
				drivenChordEnergy += driven * driven;
				++chordSamples;
			}
		}
		const double cleanChordRms = std::sqrt(cleanChordEnergy /
			static_cast<double>(chordSamples));
		const double drivenChordRms = std::sqrt(drivenChordEnergy /
			static_cast<double>(chordSamples));
		const double cleanChordCrest = cleanChordPeak / cleanChordRms;
		const double drivenChordCrest = drivenChordPeak / drivenChordRms;
		if (!(drivenChordRailMinimum > 34.4 && drivenChordCrest > 2.0 &&
			drivenChordCrest > 0.60 * cleanChordCrest &&
			drivenChordRms > 0.75 * cleanChordRms))
			std::cerr << "electric piano driven chord rail/crest/rms ratio " <<
				drivenChordRailMinimum << "/" << drivenChordCrest << "/" <<
				drivenChordRms / cleanChordRms << '\n';
		Check(drivenChordRailMinimum > 34.4 && drivenChordCrest > 2.0 &&
			drivenChordCrest > 0.60 * cleanChordCrest &&
			drivenChordRms > 0.75 * cleanChordRms,
			"electric piano Drive preserves polyphonic power-stage headroom and chord crest");

		tfdsp::ElectricPianoAmplifier trebleCutAmplifier;
		tfdsp::ElectricPianoAmplifier trebleBoostAmplifier;
		trebleCutAmplifier.SetSampleRate(48000.0);
		trebleBoostAmplifier.SetSampleRate(48000.0);
		auto trebleCutNoteControls = controls;
		auto trebleBoostNoteControls = controls;
		trebleCutNoteControls.amplifierTreble = 0.0;
		trebleBoostNoteControls.amplifierTreble = 1.0;
		std::vector<double> trebleCutNote(highVelocityPickupSignal.size());
		std::vector<double> trebleBoostNote(highVelocityPickupSignal.size());
		for (std::size_t sample = 0; sample < highVelocityPickupSignal.size(); ++sample)
		{
			trebleCutNote[sample] = trebleCutAmplifier.Step(
				5.0 * highVelocityPickupSignal[sample], trebleCutNoteControls)[0];
			trebleBoostNote[sample] = trebleBoostAmplifier.Step(
				5.0 * highVelocityPickupSignal[sample], trebleBoostNoteControls)[0];
		}
		const double renderedTrebleResidual = levelMatchedResidual(trebleCutNote,
			trebleBoostNote, 4800);
		if (!(renderedTrebleResidual > 0.20))
			std::cerr << "electric piano rendered Treble residual " <<
				renderedTrebleResidual << '\n';
		Check(renderedTrebleResidual > 0.20,
			"electric piano Treble extremes are audible on a rendered hard note");

		auto amplifierRmsAt = [](double frequency,
			const tfdsp::ElectricPianoControls& amplifierControls)
		{
			tfdsp::ElectricPianoAmplifier testAmplifier;
			testAmplifier.SetSampleRate(48000.0);
			constexpr int warmup = 12000;
			constexpr int count = 48000;
			double energy = 0.0;
			for (int sample = 0; sample < warmup + count; ++sample)
			{
				const double input = 0.02 * std::sin(
					2.0 * tfdsp::PI * frequency * sample / 48000.0);
				const double output = testAmplifier.Step(input,
					amplifierControls)[0];
				if (sample >= warmup)
					energy += output * output;
			}
			return std::sqrt(energy / count);
		};
		auto bassCutControls = controls;
		auto bassBoostControls = controls;
		auto trebleCutControls = controls;
		auto trebleBoostControls = controls;
		bassCutControls.drive = bassBoostControls.drive = 0.0;
		trebleCutControls.drive = trebleBoostControls.drive = 0.0;
		bassCutControls.amplifierBass = 0.0;
		bassBoostControls.amplifierBass = 1.0;
		trebleCutControls.amplifierTreble = 0.0;
		trebleBoostControls.amplifierTreble = 1.0;
		const double bassRange = amplifierRmsAt(80.0, bassBoostControls) /
			amplifierRmsAt(80.0, bassCutControls);
		const double trebleRange = amplifierRmsAt(8000.0, trebleBoostControls) /
			amplifierRmsAt(8000.0, trebleCutControls);
		Check(bassRange > 4.0 && trebleRange > 3.5,
			"electric piano Peterson Bass and Treble controls span audible shelves");
		auto bassQuarterControls = controls;
		auto bassThreeQuarterControls = controls;
		auto trebleQuarterControls = controls;
		auto trebleThreeQuarterControls = controls;
		auto centredToneControls = controls;
		bassQuarterControls.drive = bassThreeQuarterControls.drive = 0.0;
		trebleQuarterControls.drive = trebleThreeQuarterControls.drive = 0.0;
		centredToneControls.drive = 0.0;
		bassQuarterControls.amplifierBass = 0.25;
		bassThreeQuarterControls.amplifierBass = 0.75;
		trebleQuarterControls.amplifierTreble = 0.25;
		trebleThreeQuarterControls.amplifierTreble = 0.75;
		const std::array<double, 5> bassTrajectory{
			amplifierRmsAt(80.0, bassCutControls),
			amplifierRmsAt(80.0, bassQuarterControls),
			amplifierRmsAt(80.0, centredToneControls),
			amplifierRmsAt(80.0, bassThreeQuarterControls),
			amplifierRmsAt(80.0, bassBoostControls)};
		const std::array<double, 5> trebleTrajectory{
			amplifierRmsAt(8000.0, trebleCutControls),
			amplifierRmsAt(8000.0, trebleQuarterControls),
			amplifierRmsAt(8000.0, centredToneControls),
			amplifierRmsAt(8000.0, trebleThreeQuarterControls),
			amplifierRmsAt(8000.0, trebleBoostControls)};
		auto everyQuarterAudible = [](const std::array<double, 5>& trajectory)
		{
			for (std::size_t index = 1; index < trajectory.size(); ++index)
				if (20.0 * std::log10(trajectory[index] /
					trajectory[index - 1]) < 2.5)
					return false;
			return true;
		};
		Check(everyQuarterAudible(bassTrajectory) &&
			everyQuarterAudible(trebleTrajectory),
			"electric piano tone controls remain audible throughout their travel");

		auto trebleCutWithBassCut = trebleCutControls;
		auto trebleBoostWithBassCut = trebleBoostControls;
		auto trebleCutWithBassBoost = trebleCutControls;
		auto trebleBoostWithBassBoost = trebleBoostControls;
		trebleCutWithBassCut.amplifierBass = 0.0;
		trebleBoostWithBassCut.amplifierBass = 0.0;
		trebleCutWithBassBoost.amplifierBass = 1.0;
		trebleBoostWithBassBoost.amplifierBass = 1.0;
		const double trebleRangeWithBassCut = amplifierRmsAt(700.0,
			trebleBoostWithBassCut) / amplifierRmsAt(700.0,
				trebleCutWithBassCut);
		const double trebleRangeWithBassBoost = amplifierRmsAt(700.0,
			trebleBoostWithBassBoost) / amplifierRmsAt(700.0,
				trebleCutWithBassBoost);
		const double toneInteractionDb = std::abs(20.0 * std::log10(
			trebleRangeWithBassBoost / trebleRangeWithBassCut));
		if (!(toneInteractionDb > 0.01))
			std::cerr << "electric piano tone interaction " << toneInteractionDb <<
				" dB\n";
		Check(toneInteractionDb > 0.01,
			"electric piano Peterson tone controls retain shared-network interaction");

		tfdsp::ElectricPianoAmplifier aliasingAmplifier;
		aliasingAmplifier.SetSampleRate(48000.0);
		auto aliasingControls = controls;
		aliasingControls.drive = 1.0;
		aliasingControls.amplifierTreble = 1.0;
		std::vector<double> aliasingSignal(48000);
		for (int sample = 0; sample < 72000; ++sample)
		{
			const double input = 0.315 * std::sin(
				2.0 * tfdsp::PI * 7000.0 * sample / 48000.0);
			const double output = aliasingAmplifier.Step(input,
				aliasingControls)[0];
			if (sample >= 24000)
				aliasingSignal[sample - 24000] = output;
		}
		const double aliasingFundamental = HarmonicMagnitude(aliasingSignal,
			7000.0 / 48000.0);
		// A base-rate fifth harmonic at 35 kHz would fold to 13 kHz. The entire
		// nonlinear power chain now resides between the 2x interpolator
		// and decimator, so that otherwise empty bin should remain strongly rejected.
		const double fifthHarmonicAlias = HarmonicMagnitude(aliasingSignal,
			13000.0 / 48000.0) / std::max(1.0e-12, aliasingFundamental);
		if (!(fifthHarmonicAlias < 0.015))
			std::cerr << "electric piano amplifier fifth-harmonic alias ratio " <<
				fifthHarmonicAlias << '\n';
		Check(fifthHarmonicAlias < 0.015,
			"electric piano 2x amplifier path suppresses nonlinear aliasing");

		// The earlier 7 kHz/max-Drive probe did not represent the reported upper-
		// keyboard defect: transformer-output crossover error creates a long odd
		// ladder even at low Drive. Exercise a C8-adjacent fundamental at the actual
		// default and calibrated hard-note level. Seventh and ninth harmonics above
		// host Nyquist would fold to these otherwise empty bins if the 2x path or
		// low-current germanium reduction regressed.
		tfdsp::ElectricPianoAmplifier highNoteAliasingAmplifier;
		highNoteAliasingAmplifier.SetSampleRate(48000.0);
		auto highNoteAliasingControls = controls;
		highNoteAliasingControls.drive = 0.32;
		constexpr double HighNoteProbeFrequency = 4200.0;
		std::vector<double> highNoteAliasingSignal(48000);
		for (int sample = 0; sample < 72000; ++sample)
		{
			const double input = 0.235 * std::sin(2.0 * tfdsp::PI *
				HighNoteProbeFrequency * sample / 48000.0);
			const double output = highNoteAliasingAmplifier.Step(input,
				highNoteAliasingControls)[0];
			if (sample >= 24000)
				highNoteAliasingSignal[sample - 24000] = output;
		}
		const double highNoteFundamental = HarmonicMagnitude(
			highNoteAliasingSignal, HighNoteProbeFrequency / 48000.0);
		const double seventhFold = HarmonicMagnitude(highNoteAliasingSignal,
			(48000.0 - 7.0 * HighNoteProbeFrequency) / 48000.0) /
			std::max(1.0e-12, highNoteFundamental);
		const double ninthFold = HarmonicMagnitude(highNoteAliasingSignal,
			(48000.0 - 9.0 * HighNoteProbeFrequency) / 48000.0) /
			std::max(1.0e-12, highNoteFundamental);
		if (!(seventhFold < 0.001 && ninthFold < 0.001))
			std::cerr << "electric piano default high-note alias ratios " <<
				seventhFold << "/" << ninthFold << '\n';
		Check(seventhFold < 0.001 && ninthFold < 0.001,
			"electric piano default high notes reject folded crossover harmonics");

		// A single sine can miss note-dependent intermodulation. Render the same
		// deterministic voice and four-note chord at 48/96 kHz, decimate the 96 kHz
		// result, and require the maximum-Drive electronics to add no material
		// sample-rate residual beyond that already present in the physical voice.
		struct AliasRender
		{
			std::vector<double> pickup;
			std::vector<double> output;
		};
		auto renderAliasCase = [](double rate,
			const std::vector<double>& pitches)
		{
			tfdsp::ElectricPianoControls voiceControls;
			voiceControls.mechanics = 0.0;
			tfdsp::ElectricPianoControls amplifierControls;
			amplifierControls.drive = 1.0;
			std::vector<tfdsp::ElectricPianoVoice> voices(pitches.size());
			for (std::size_t voice = 0; voice < voices.size(); ++voice)
			{
				voices[voice].SetSampleRate(rate);
				voices[voice].SetNoiseSeed(0x6d2b79f5u +
					static_cast<std::uint32_t>(voice) * 0x9e3779b9u);
			}
			tfdsp::ElectricPianoAmplifier aliasAmplifier;
			aliasAmplifier.SetSampleRate(rate);
			AliasRender render;
			render.pickup.resize(static_cast<std::size_t>(rate));
			render.output.resize(render.pickup.size());
			for (std::size_t sample = 0; sample < render.pickup.size(); ++sample)
			{
				double pickup = 0.0;
				for (std::size_t voice = 0; voice < voices.size(); ++voice)
					pickup += voices[voice].Step(pitches[voice], 10.0, 1.0,
						false, voiceControls);
				render.pickup[sample] = pickup;
				render.output[sample] = aliasAmplifier.Step(5.0 * pickup,
					amplifierControls)[0];
			}
			return render;
		};
		auto downsampleAliasCase = [](const std::vector<double>& input)
		{
			auto decimator = tfdsp::CreateX2Resampler_Chebychev7();
			std::vector<double> output(input.size() / 2);
			for (std::size_t sample = 0; sample < output.size(); ++sample)
			{
				Eigen::Array<double,
					tfdsp::X2Resampler_Order7::ResamplingFactor, 1> pair;
				pair(0) = input[2 * sample];
				pair(1) = input[2 * sample + 1];
				output[sample] = decimator->Downsample(pair);
			}
			return output;
		};
		for (const auto& pitches : std::array<std::vector<double>, 2>{{
			{0.0}, {-1.0, -8.0 / 12.0, -5.0 / 12.0, 0.0}}})
		{
			const auto at48 = renderAliasCase(48000.0, pitches);
			const auto at96 = renderAliasCase(96000.0, pitches);
			auto pickup96 = downsampleAliasCase(at96.pickup);
			auto output96 = downsampleAliasCase(at96.output);
			// The nested amplifier interpolator/decimator contributes one 48 kHz
			// sample relative to the voice-only reference.
			for (std::size_t sample = output96.size() - 1; sample > 0; --sample)
				output96[sample] = output96[sample - 1];
			const double pickupRateResidual = levelMatchedResidual(at48.pickup,
				pickup96, 4800);
			const double drivenRateResidual = levelMatchedResidual(at48.output,
				output96, 4800);
			if (!(drivenRateResidual < pickupRateResidual + 0.0015))
				std::cerr << "electric piano rendered alias residual voices/pickup/amp "
					<< pitches.size() << "/" << pickupRateResidual << "/" <<
					drivenRateResidual << '\n';
			Check(drivenRateResidual < pickupRateResidual + 0.0015,
				"electric piano maximum-Drive note/chord adds no rendered alias residual");
		}

		auto halfOutputControls = controls;
		auto fullOutputControls = controls;
		halfOutputControls.drive = fullOutputControls.drive = 0.0;
		halfOutputControls.outputVolume = 0.25;
		fullOutputControls.outputVolume = 0.50;
		const double outputVolumeRatio = amplifierRmsAt(1000.0,
			fullOutputControls) / amplifierRmsAt(1000.0, halfOutputControls);
		Check(outputVolumeRatio > 1.95 && outputVolumeRatio < 2.05,
			"electric piano Output changes level after the nonlinear amplifier");

		struct VibratoMetrics
		{
			int crossings{};
			double differenceEnergy{};
			double wetEnergy{};
			double dryEnergy{};
		};
		auto measureVibrato = [&](double speed)
		{
			tfdsp::ElectricPianoAmplifier dryAmplifier;
			tfdsp::ElectricPianoAmplifier wetAmplifier;
			dryAmplifier.SetSampleRate(48000.0);
			wetAmplifier.SetSampleRate(48000.0);
			auto dryControls = controls;
			auto wetControls = controls;
			dryControls.drive = wetControls.drive = 0.0;
			dryControls.vibratoIntensity = 0.0;
			wetControls.vibratoIntensity = 1.0;
			wetControls.vibratoSpeed = speed;
			VibratoMetrics metrics;
			int previousSign = 0;
			double stereoBalance = 0.0;
			const double balanceCoefficient = 1.0 - std::exp(
				-2.0 * tfdsp::PI * 20.0 / 48000.0);
			for (int sample = 0; sample < 120000; ++sample)
			{
				const double input = 0.05 * std::sin(
					2.0 * tfdsp::PI * 250.0 * sample / 48000.0);
				const auto dry = dryAmplifier.Step(input, dryControls);
				const auto wet = wetAmplifier.Step(input, wetControls);
				if (sample < 24000)
					continue;
				const double wetPower = wet[0] * wet[0] + wet[1] * wet[1];
				metrics.differenceEnergy += (wet[0] - wet[1]) *
					(wet[0] - wet[1]);
				metrics.wetEnergy += wetPower;
				metrics.dryEnergy += dry[0] * dry[0] + dry[1] * dry[1];
				stereoBalance += balanceCoefficient *
					(wet[0] * wet[0] - wet[1] * wet[1] - stereoBalance);
				if (wetPower > 1.0e-8 && std::abs(stereoBalance) > 1.0e-7)
				{
					const int sign = stereoBalance > 0.0 ? 1 : -1;
					if (previousSign != 0 && sign != previousSign)
						++metrics.crossings;
					previousSign = sign;
				}
			}
			return metrics;
		};
		const auto slowVibrato = measureVibrato(0.0);
		const auto fastVibrato = measureVibrato(1.0);
		const double vibratoPowerRatio = slowVibrato.wetEnergy /
			std::max(1.0e-20, slowVibrato.dryEnergy);
		Check(slowVibrato.differenceEnergy > 0.10 * slowVibrato.wetEnergy &&
			vibratoPowerRatio > 0.98 && vibratoPowerRatio < 1.02,
			"electric piano Vibrato produces an equal-power stereo pan");
		if (!(slowVibrato.crossings >= 4 &&
			fastVibrato.crossings > 3 * slowVibrato.crossings))
			std::cerr << "electric piano vibrato crossings slow/fast " <<
				slowVibrato.crossings << "/" << fastVibrato.crossings << '\n';
		Check(slowVibrato.crossings >= 4 &&
			fastVibrato.crossings > 3 * slowVibrato.crossings,
			"electric piano Vibrato Speed spans a clearly different rate");

		tfdsp::ElectricPianoAmplifier referenceAmplifier;
		tfdsp::ElectricPianoAmplifier stressedAmplifier;
		referenceAmplifier.SetSampleRate(48000.0);
		stressedAmplifier.SetSampleRate(48000.0);
		auto recoveryControls = controls;
		recoveryControls.drive = 0.75;
		double referenceEarlyEnergy = 0.0;
		double stressedEarlyEnergy = 0.0;
		double referenceLateEnergy = 0.0;
		double stressedLateEnergy = 0.0;
		double stressedRailAfterBurst = 35.0;
		double stressedRailLate = 35.0;
		double referenceRailLate = 35.0;
		for (int sample = 0; sample < 36000; ++sample)
		{
			const double phase = 2.0 * tfdsp::PI * 220.0 * sample / 48000.0;
			const double probe = 0.20 * std::sin(phase);
			// A 0.8 model-unit chord is already far beyond one hard note after the
			// module's 5x pickup sum. Do not drive the unmodelled Figure 11-8 source
			// to an impossible hundreds-of-volts output merely to test Figure 11-9.
			const double stressedInput = sample >= 2000 && sample < 6000 ?
				0.8 * std::sin(phase) : probe;
			const double reference = referenceAmplifier.Step(probe,
				recoveryControls)[0];
			const double stressed = stressedAmplifier.Step(stressedInput,
				recoveryControls)[0];
			if (sample == 6000)
				stressedRailAfterBurst = stressedAmplifier.SupplyRailVoltage();
			if (sample == 35999)
			{
				stressedRailLate = stressedAmplifier.SupplyRailVoltage();
				referenceRailLate = referenceAmplifier.SupplyRailVoltage();
			}
			if (sample >= 6100 && sample < 7600)
			{
				referenceEarlyEnergy += reference * reference;
				stressedEarlyEnergy += stressed * stressed;
			}
			if (sample >= 34000)
			{
				referenceLateEnergy += reference * reference;
				stressedLateEnergy += stressed * stressed;
			}
		}
		const double earlyRecoveryRatio = stressedEarlyEnergy /
			std::max(1.0e-20, referenceEarlyEnergy);
		const double lateRecoveryRatio = stressedLateEnergy /
			std::max(1.0e-20, referenceLateEnergy);
		const double railGapAfterBurst = std::abs(stressedRailAfterBurst -
			referenceRailLate);
		const double railGapLate = std::abs(stressedRailLate - referenceRailLate);
		// A nonlinear preamp can make a severe burst compress the subsequent
		// power-stage demand, so the shared rail may approach its steady value
		// from either side. Test recovery toward the reference, not its direction.
		if (!(railGapAfterBurst > 0.05 && railGapLate < 0.05 &&
			railGapLate < 0.25 * railGapAfterBurst))
			std::cerr << "electric piano supply rail after/late " <<
				stressedRailAfterBurst << "/" << stressedRailLate <<
				" reference " << referenceRailLate <<
				" output-energy early/late " << earlyRecoveryRatio << "/" <<
				lateRecoveryRatio << '\n';
		Check(railGapAfterBurst > 0.05 && railGapLate < 0.05 &&
			railGapLate < 0.25 * railGapAfterBurst,
			"electric piano power-stage rail headroom recovers after overload");
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
	normalledOscillator.SetSampleRate(48000.0);
	detailedOscillator.SetSampleRate(48000.0);
	double normalledDifference = 0.0;
	bool externalPathFinite = true;
	for (int i = 0; i < 10000; ++i)
	{
		const double normalled = normalledOscillator.Step(
			440.0, 0.4, 0.7, -0.2);
		const auto detailed = detailedOscillator.StepWithInput(
			440.0, 0.4, 0.7, -0.2, 0.0, false);
		normalledDifference = std::max(normalledDifference,
			std::abs(normalled - detailed.folded));
		externalPathFinite = externalPathFinite &&
			std::isfinite(detailed.oscillator) &&
			std::isfinite(detailed.folded);
	}
	Check(normalledDifference < 1.0e-12,
		"normalled folder input preserves the original oscillator path");
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
