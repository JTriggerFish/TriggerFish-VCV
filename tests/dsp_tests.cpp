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
			const double angle = TwoPi * normalizedFrequency * index;
			real += signal[index] * std::cos(angle);
			imaginary -= signal[index] * std::sin(angle);
		}
		return std::hypot(real, imaginary) / signal.size();
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
			stressedControls.gap = 1.0;
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
		farPickupControls.gap = 0.0;
		nearPickupControls.gap = 1.0;
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
		Check(nearPickupEnergy > 2.0 * farPickupEnergy,
			"electric piano pickup gap changes physical sensitivity as well as curvature");

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
				editedControls.gap = 0.9;
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
			std::cerr << "electric piano keyboard energy ratio: " <<
				maximumKeyboardEnergy / minimumKeyboardEnergy << '\n';
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
			"electric piano finite-pole pickup develops progressive bark with velocity");

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
			brightSecond > 8.0 * darkSecond && brightThird > 1.6 * darkThird))
			std::cerr << "electric piano Tone ratios attack dark/bright " <<
				darkFirst << "/" << brightFirst << " H2 " << darkSecond <<
				"/" << brightSecond << " H3 " << darkThird << "/" <<
				brightThird << '\n';
		Check(brightFirst > 0.80 * darkFirst && brightFirst < 1.30 * darkFirst &&
			brightSecond > 8.0 * darkSecond &&
			brightThird > 1.6 * darkThird,
			"electric piano Tone changes physical pickup alignment and harmonic balance");

		auto amplify = [&](const std::vector<double>& direct)
		{
			tfdsp::ElectricPianoAmplifier testAmplifier;
			testAmplifier.SetSampleRate(48000.0);
			std::vector<double> amplified(direct.size());
			for (std::size_t sample = 0; sample < direct.size(); ++sample)
				amplified[sample] = testAmplifier.Step(5.0 * direct[sample],
					controls.drive);
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
		Check(levelMatchedResidual(darkToneAmplified, brightToneAmplified,
			2048) > 0.32,
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
		double amplifierPeak = 0.0;
		for (int sample = 0; sample < 48000; ++sample)
		{
			const double input = 16.0 * std::sin(
				2.0 * tfdsp::PI * 220.0 * sample / 48000.0);
			amplifierPeak = std::max(amplifierPeak,
				std::abs(amplifier.Step(input, 1.0)));
		}
		Check(std::isfinite(amplifierPeak) && amplifierPeak < 2.0,
			"electric piano shared amplifier overload remains bounded");

		tfdsp::ElectricPianoAmplifier referenceAmplifier;
		tfdsp::ElectricPianoAmplifier stressedAmplifier;
		referenceAmplifier.SetSampleRate(48000.0);
		stressedAmplifier.SetSampleRate(48000.0);
		double referenceEarlyEnergy = 0.0;
		double stressedEarlyEnergy = 0.0;
		double referenceLateEnergy = 0.0;
		double stressedLateEnergy = 0.0;
		for (int sample = 0; sample < 18000; ++sample)
		{
			const double phase = 2.0 * tfdsp::PI * 220.0 * sample / 48000.0;
			const double probe = 0.08 * std::sin(phase);
			const double stressedInput = sample >= 2000 && sample < 6000 ?
				8.0 * std::sin(phase) : probe;
			const double reference = referenceAmplifier.Step(probe, 0.75);
			const double stressed = stressedAmplifier.Step(stressedInput, 0.75);
			if (sample >= 6100 && sample < 7600)
			{
				referenceEarlyEnergy += reference * reference;
				stressedEarlyEnergy += stressed * stressed;
			}
			if (sample >= 16000)
			{
				referenceLateEnergy += reference * reference;
				stressedLateEnergy += stressed * stressed;
			}
		}
		Check(stressedEarlyEnergy < 0.90 * referenceEarlyEnergy &&
			stressedLateEnergy > stressedEarlyEnergy,
			"electric piano shared amplifier has level-dependent recovery memory");
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
