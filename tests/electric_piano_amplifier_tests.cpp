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
		tfdsp::ElectricPianoControls controls;
		auto pickupControls = controls;
		pickupControls.mechanics = 0.0;
		tfdsp::ElectricPianoVoice highVelocityPickupVoice;
		highVelocityPickupVoice.SetSampleRate(48000.0);
		std::vector<double> highVelocityPickupSignal(16000);
		for (double& output : highVelocityPickupSignal)
			output = highVelocityPickupVoice.Step(0.0, 10.0, 1.0, false,
				pickupControls);

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

		auto completePreampSmallSignalGain = [](double frequency)
		{
			constexpr double rate = 96000.0;
			tfdsp::PetersonPreAmplifier inputStage;
			tfdsp::PetersonTonePreAmplifier toneStage;
			inputStage.SetSampleRate(rate);
			toneStage.SetSampleRate(rate);
			double inputEnergy = 0.0;
			double outputEnergy = 0.0;
			for (int sample = 0; sample < 48000; ++sample)
			{
				const double input = 0.001 * std::sin(2.0 * tfdsp::PI *
					frequency * sample / rate);
				const double output = toneStage.Step(
					inputStage.Step(input).voltage, 0.5, 0.5).voltage;
				if (sample >= 24000)
				{
					inputEnergy += input * input;
					outputEnergy += output * output;
				}
			}
			return std::sqrt(outputEnergy / inputEnergy);
		};
		const double completePreampGain1k = completePreampSmallSignalGain(1000.0);
		const double completePreampResponse5k = 20.0 * std::log10(
			completePreampSmallSignalGain(5000.0) / completePreampGain1k);
		const double completePreampResponse8k = 20.0 * std::log10(
			completePreampSmallSignalGain(8000.0) / completePreampGain1k);
		Check(std::abs(completePreampResponse5k - (-1.92)) < 0.15 &&
			std::abs(completePreampResponse8k - (-8.38)) < 0.55,
			"electric piano trapezoidal Figure 11-8 response agrees with ngspice at 5/8 kHz");

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
		// of a hard middle-C note is about 0.0438 after the tuning-spring modal
		// correction, and the module feeds five times
		// that value to the amplifier. Drive raises the voltage into Figure 11-8,
		// while its reciprocal at the schematic's volume node keeps Figure 11-9's
		// ideal small-signal demand independent of the knob.
		double renderedHardNotePeak = 0.0;
		for (double sample : highVelocityPickupSignal)
			renderedHardNotePeak = std::max(renderedHardNotePeak,
				std::abs(sample));
		Check(renderedHardNotePeak > 0.040 && renderedHardNotePeak < 0.046,
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
		// Restoring short-lived upper-mode attack energy raises the instantaneous C4
		// pickup peak slightly without changing its sustained level or the amplifier
		// conversion. Keep the voltage guard tight enough to catch gain-staging
		// drift while admitting that physically corrected transient.
		Check(hardNoteCircuitInput > 0.42 && hardNoteCircuitInput < 0.47 &&
			idealHardNotePowerPeak > 1.50 && idealHardNotePowerPeak < 1.70 &&
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
		const double drivenNoteResidual = LevelMatchedResidual(cleanNote,
			drivenNote, 0, 4800);
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
		const double renderedTrebleResidual = LevelMatchedResidual(trebleCutNote,
			trebleBoostNote, 0, 4800);
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
			const std::vector<double>& pitches, bool harmonicFm = false,
			bool forceModulationPath = false)
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
				{
					tfdsp::ElectricPianoModulation modulation;
					if (harmonicFm)
						modulation.linearFrequencyRatio = 0.18 * std::sin(
							2.0 * tfdsp::PI *
							tfdsp::ElectricPianoReferenceFrequency * sample / rate);
					else if (forceModulationPath)
						modulation.phaseRadians = 1.0e-15;
					pickup += voices[voice].Step(pitches[voice], 10.0, 1.0,
						false, voiceControls, modulation);
				}
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
			const double pickupRateResidual = LevelMatchedResidual(at48.pickup,
				pickup96, 0, 4800);
			const double drivenRateResidual = LevelMatchedResidual(at48.output,
				output96, 0, 4800);
			if (!(drivenRateResidual < pickupRateResidual + 0.0015))
				std::cerr << "electric piano rendered alias residual voices/pickup/amp "
					<< pitches.size() << "/" << pickupRateResidual << "/" <<
					drivenRateResidual << '\n';
			Check(drivenRateResidual < pickupRateResidual + 0.0015,
				"electric piano maximum-Drive note/chord adds no rendered alias residual");
		}
		const auto fmAt48 = renderAliasCase(48000.0, {0.0}, true);
		const auto fmAt96 = renderAliasCase(96000.0, {0.0}, true);
		const auto dryAt48 = renderAliasCase(48000.0, {0.0}, false, true);
		const auto dryAt96 = renderAliasCase(96000.0, {0.0}, false, true);
		auto fmPickup96 = downsampleAliasCase(fmAt96.pickup);
		auto dryPickup96 = downsampleAliasCase(dryAt96.pickup);
		const double fmRateResidual = LevelMatchedResidual(fmAt48.pickup,
			fmPickup96, 0, 4800);
		const double dryRateResidual = LevelMatchedResidual(dryAt48.pickup,
			dryPickup96, 0, 4800);
		if (!(fmRateResidual < dryRateResidual + 0.005))
			std::cerr << "electric piano dry/FM rate residual " << dryRateResidual <<
				"/" << fmRateResidual << '\n';
		Check(fmRateResidual < dryRateResidual + 0.005,
			"electric piano audio-rate FM adds no material sample-rate alias residual");

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

	if (failures == 0)
		std::cout << "All electric-piano tests passed\n";
	return failures == 0 ? 0 : 1;
}
