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
		// The production pickup uses a precomputed two-dimensional calibrated
		// field. Bound bilinear lookup error against direct evaluation of the same
		// scalar-flux gradient across ordinary and service-limit alignments.
		double maximumPickupLutRelativeError = 0.0;
		for (double vertical : {-0.8, -0.2, 0.34, 0.8, 1.4})
		{
			for (double horizontal : {0.05, 0.35})
			{
				for (double gap : {0.51, 1.0, 1.6, 3.0})
				{
					const auto lookup = tfdsp::ElectricPianoVoice::
						MagneticPickupGradient(vertical, horizontal, gap);
					const auto direct = tfdsp::ElectricPianoVoice::
						DirectMagneticPickupGradient(vertical, horizontal, gap);
					const double error = std::hypot(lookup[0] - direct[0],
						lookup[1] - direct[1]);
					maximumPickupLutRelativeError = std::max(
						maximumPickupLutRelativeError, error / std::max(1.0e-12,
							std::hypot(direct[0], direct[1])));
				}
			}
		}
		Check(maximumPickupLutRelativeError < 0.003,
			"electric piano 2D magnetic-pickup LUT follows its calibrated flux field");

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
			lowKey.tineLengthMetres > middleKey.tineLengthMetres &&
			middleKey.tineLengthMetres > highKey.tineLengthMetres &&
			lowKey.tineModalMassKg > middleKey.tineModalMassKg &&
			middleKey.tineModalMassKg > highKey.tineModalMassKg &&
			displacementPerImpulse(lowKey) >
				displacementPerImpulse(middleKey) &&
			displacementPerImpulse(middleKey) >
				displacementPerImpulse(highKey) &&
			lowKey.pickupSensitivity > 0.0 &&
			middleKey.pickupSensitivity > 0.0 &&
			highKey.pickupSensitivity > 0.0 &&
			std::abs(tfdsp::ElectricPianoPublishedMechanicalData::HammerMassKg -
				0.011) < 1.0e-12 &&
			tfdsp::ElectricPianoFactoryHammerTipDurometer(0.0) == 30.0 &&
			tfdsp::ElectricPianoFactoryHammerTipDurometer(0.5) == 50.0 &&
			tfdsp::ElectricPianoFactoryHammerTipDurometer(1.0) == 100.0,
			"electric piano key model retains published SI geometry and factory tip zones");
		Check(tfdsp::ElectricPianoModeBandlimitGain(0.39 * 48000.0,
			48000.0) == 1.0 &&
			tfdsp::ElectricPianoModeBandlimitGain(0.445 * 48000.0,
				48000.0) > 0.0 &&
			tfdsp::ElectricPianoModeBandlimitGain(0.445 * 48000.0,
				48000.0) < 1.0 &&
			tfdsp::ElectricPianoModeBandlimitGain(0.445 * 48000.0,
				48000.0) ==
				tfdsp::ElectricPianoModeBandlimitGain(0.445 * 48000.0,
					96000.0) &&
			tfdsp::ElectricPianoModeBandlimitGain(0.49 * 48000.0,
				48000.0) == 0.0,
			"electric piano high modes use the oversampled pickup bandwidth and retain the 48 kHz model at higher host rates");

		const auto isolatedFork = tfdsp::MakeElectricPianoCoupledForkProfile(
			0.0, middleKey.keyboardPosition);
		const auto factoryFork = tfdsp::MakeElectricPianoCoupledForkProfile(
			controls.coupling, middleKey.keyboardPosition);
		const auto coupledFork = tfdsp::MakeElectricPianoCoupledForkProfile(
			1.0, middleKey.keyboardPosition);
		Check(std::abs(controls.coupling - 0.5) < 1.0e-12,
			"electric piano Coupling defaults to the calibrated control midpoint");
		Check(isolatedFork.toneBarFreeFrequencyRatio > 0.995 &&
			isolatedFork.toneBarFreeFrequencyRatio < 1.0 &&
			isolatedFork.frequencyRatios[0] > 0.99 &&
			std::abs(isolatedFork.frequencyRatios[1] - 1.0) < 1.0e-12,
			"electric piano tine and tone-bar fundamental coordinates share the played pitch");
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
			2.5 * factoryFork.supportReactionLossFactors[1] &&
			factoryFork.supportReactionLossFactors[1] >
			10.0 * coupledFork.supportReactionLossFactors[1] &&
			coupledFork.supportReactionLossFactors[1] < 0.05,
			"electric piano tone bar cancels lossy common-support reaction at strong coupling");
		const double weakSubModeTransfer =
			tfdsp::ElectricPianoToneBarSubModeCouplingScale(0.0,
				middleKey.keyboardPosition);
		const double factorySubModeTransfer =
			tfdsp::ElectricPianoToneBarSubModeCouplingScale(0.5,
				middleKey.keyboardPosition);
		const double strongSubModeTransfer =
			tfdsp::ElectricPianoToneBarSubModeCouplingScale(1.0,
				middleKey.keyboardPosition);
		Check(weakSubModeTransfer < 0.80 &&
			factorySubModeTransfer == 1.0 && strongSubModeTransfer > 1.10,
			"electric piano Coupling reciprocally changes separate tone-bar transfer around an exact neutral midpoint");

		// Check the published F1/F3 measurements directly, then the reduced
		// tuning-spring point-mass extrapolation and short-tine shear correction.
		tfdsp::ElectricPianoVoice bassCalibrationVoice;
		tfdsp::ElectricPianoVoice f1CalibrationVoice;
		tfdsp::ElectricPianoVoice f3CalibrationVoice;
		tfdsp::ElectricPianoVoice a4CalibrationVoice;
		tfdsp::ElectricPianoVoice trebleCalibrationVoice;
		bassCalibrationVoice.SetSampleRate(48000.0);
		f1CalibrationVoice.SetSampleRate(48000.0);
		f3CalibrationVoice.SetSampleRate(48000.0);
		a4CalibrationVoice.SetSampleRate(48000.0);
		trebleCalibrationVoice.SetSampleRate(48000.0);
		bassCalibrationVoice.Step((28.0 - 60.0) / 12.0, 10.0, 0.7, false,
			controls);
		f1CalibrationVoice.Step((29.0 - 60.0) / 12.0, 10.0, 0.7, false,
			controls);
		f3CalibrationVoice.Step((53.0 - 60.0) / 12.0, 10.0, 0.7, false,
			controls);
		a4CalibrationVoice.Step((69.0 - 60.0) / 12.0, 10.0, 0.7, false,
			controls);
		trebleCalibrationVoice.Step(3.0, 10.0, 0.7, false, controls);
		Check(std::abs(bassCalibrationVoice.ModeFrequencyRatio(3) - 7.11) <
				1.0e-6 &&
			std::abs(bassCalibrationVoice.ModeFrequencyRatio(4) - 20.25) <
				1.0e-6 &&
			std::abs(f1CalibrationVoice.ModeFrequencyRatio(3) - 7.2) < 1.0e-6 &&
			std::abs(f1CalibrationVoice.ModeFrequencyRatio(4) - 20.6) < 1.0e-6 &&
			std::abs(f3CalibrationVoice.ModeFrequencyRatio(3) - 7.4) < 1.0e-6 &&
			std::abs(f3CalibrationVoice.ModeFrequencyRatio(4) - 20.7) < 1.0e-6 &&
			std::abs(f3CalibrationVoice.ModeFrequencyRatio(5) - 38.7) < 1.0e-6 &&
			a4CalibrationVoice.ModeFrequencyRatio(3) > 7.51 &&
			a4CalibrationVoice.ModeFrequencyRatio(3) < 7.54 &&
			trebleCalibrationVoice.ModeFrequencyRatio(3) > 7.70 &&
			trebleCalibrationVoice.ModeFrequencyRatio(3) < 7.75 &&
			trebleCalibrationVoice.ModeFrequencyRatio(9) > 140.0 &&
			trebleCalibrationVoice.ModeFrequencyRatio(9) < 146.0,
			"electric piano attack-mode ratios pass through published keyboard anchors");
		Check(std::abs(f1CalibrationVoice.ModeAmplitudeLifetimeSeconds(3) -
			8.685889638065037 / 21.1) < 1.0e-6 &&
			std::abs(f1CalibrationVoice.ModeAmplitudeLifetimeSeconds(4) -
				8.685889638065037 / 67.7) < 1.0e-6 &&
			std::abs(f3CalibrationVoice.ModeAmplitudeLifetimeSeconds(3) -
				8.685889638065037 / 294.0) < 1.0e-6 &&
			std::abs(f3CalibrationVoice.ModeAmplitudeLifetimeSeconds(4) -
				8.685889638065037 / 37.0) < 1.0e-6 &&
			std::abs(f3CalibrationVoice.ModeAmplitudeLifetimeSeconds(5) -
				8.685889638065037 / 161.0) < 1.0e-6,
			"electric piano attack-mode lifetimes match published F1/F3 decay slopes");
		// Pickup bandwidth must not also truncate the mechanical point mobility.
		// The long bass tine keeps every currently modelled attack coordinate in
		// band. By E5, mode 5 is ultrasonic at host rate but is still well resolved
		// by the 64x collision solve and must continue to load the hammer.
		tfdsp::ElectricPianoVoice e5ContactBandwidthVoice;
		e5ContactBandwidthVoice.SetSampleRate(48000.0);
		e5ContactBandwidthVoice.Step((76.0 - 60.0) / 12.0, 10.0, 0.86,
			false, controls);
		Check(bassCalibrationVoice.ModeRendered(9) &&
			!e5ContactBandwidthVoice.ModeRendered(5) &&
			e5ContactBandwidthVoice.ModeParticipatesInContact(5),
			"electric piano ultrasonic modes load the hammer without entering pickup output");
		Check(std::abs(f1CalibrationVoice.ModeFrequencyRatio(
				tfdsp::ElectricPianoToneBarSubMode) - 0.83) < 1.0e-9 &&
			std::abs(f3CalibrationVoice.ModeFrequencyRatio(
				tfdsp::ElectricPianoToneBarSubMode) - 0.58) < 1.0e-9 &&
			std::abs(f1CalibrationVoice.ModeAmplitudeLifetimeSeconds(
				tfdsp::ElectricPianoToneBarSubMode) -
				8.685889638065037 / 9.1) < 1.0e-6 &&
			std::abs(f3CalibrationVoice.ModeAmplitudeLifetimeSeconds(
				tfdsp::ElectricPianoToneBarSubMode) -
				8.685889638065037 / 138.0) < 1.0e-6,
			"electric piano separate tone-bar submode matches published F1/F3 frequency and decay anchors");
		// F3 is the last measured decay anchor. Its anomalously long fourth and
		// fifth modal lifetimes were measured on a tuning-spring-loaded structure
		// whose frequency ratios lie away from a uniform cantilever. The fitted
		// tuning mass remains in the frequency model above F3, but there is no
		// evidence that F3's anomalous modal losses remain with it. The damping
		// model therefore reaches Sonderbo's distributed beam law at A4. Freezing
		// the F3 multipliers made the 7.7--18 kHz treble modes ring for 44--99 ms.
		Check(a4CalibrationVoice.ModeAmplitudeLifetimeSeconds(3) > 0.020 &&
			a4CalibrationVoice.ModeAmplitudeLifetimeSeconds(3) < 0.030 &&
			a4CalibrationVoice.ModeAmplitudeLifetimeSeconds(4) > 0.007 &&
			a4CalibrationVoice.ModeAmplitudeLifetimeSeconds(4) < 0.011 &&
			a4CalibrationVoice.ModeAmplitudeLifetimeSeconds(5) > 0.0035 &&
			a4CalibrationVoice.ModeAmplitudeLifetimeSeconds(5) < 0.006 &&
			a4CalibrationVoice.ModeAmplitudeLifetimeSeconds(3) >
				a4CalibrationVoice.ModeAmplitudeLifetimeSeconds(4) &&
			a4CalibrationVoice.ModeAmplitudeLifetimeSeconds(4) >
				a4CalibrationVoice.ModeAmplitudeLifetimeSeconds(5),
			"electric piano treble attack-mode damping returns to the sourced beam law at A4");
		// The measured damping and the spring-loaded eigenfunctions jointly set the
		// reduced modal residue. Lock their finite-contact projection/lifetime
		// products at both measurement anchors so frequency, shape, mass and damping
		// cannot silently become four unrelated approximations again.
		auto modalContactEnergy = [](const auto& voice, std::size_t mode)
		{
			const double projection = voice.ContactModeProjection(mode);
			return projection * projection *
				voice.ModeAmplitudeLifetimeSeconds(mode);
		};
		Check(modalContactEnergy(f1CalibrationVoice, 3) > 0.120 &&
			modalContactEnergy(f1CalibrationVoice, 3) < 0.135 &&
			modalContactEnergy(f1CalibrationVoice, 4) > 0.024 &&
			modalContactEnergy(f1CalibrationVoice, 4) < 0.030 &&
			modalContactEnergy(f3CalibrationVoice, 3) > 0.0085 &&
			modalContactEnergy(f3CalibrationVoice, 3) < 0.0100 &&
			modalContactEnergy(f3CalibrationVoice, 4) > 0.0060 &&
			modalContactEnergy(f3CalibrationVoice, 4) < 0.0070 &&
			modalContactEnergy(f3CalibrationVoice, 5) > 0.0027 &&
			modalContactEnergy(f3CalibrationVoice, 5) < 0.0034,
			"electric piano spring-loaded modes retain calibrated contact residues");

		// Table 8.6 reports the separate submode at -28 dB relative to the played
		// fundamental for F1 and -4.5 dB for F3. Complete the physical collision
		// before reading analytic states; the same participation enters hammer and
		// pickup projections, so this locks a reciprocal residue rather than an
		// output-only gain.
		tfdsp::ElectricPianoVoice f1ResidueVoice;
		tfdsp::ElectricPianoVoice f3ResidueVoice;
		tfdsp::ElectricPianoVoice upperSubModeVoice;
		f1ResidueVoice.SetSampleRate(48000.0);
		f3ResidueVoice.SetSampleRate(48000.0);
		upperSubModeVoice.SetSampleRate(48000.0);
		auto completeCollision = [&](auto& voice, double pitch)
		{
			voice.Step(pitch, 10.0, 0.86, false, controls);
			for (int sample = 0; sample < 4096 && voice.ContactActive(); ++sample)
				voice.Step(pitch, 10.0, 0.86, false, controls);
		};
		completeCollision(f1ResidueVoice, (29.0 - 60.0) / 12.0);
		completeCollision(f3ResidueVoice, (53.0 - 60.0) / 12.0);
		completeCollision(upperSubModeVoice, (76.0 - 60.0) / 12.0);
		auto pickupDisplacementRatioDb = [](const auto& voice,
			std::size_t mode)
		{
			return 20.0 * std::log10(std::max(1.0e-20,
				voice.ModePickupDisplacementAmplitude(mode)) /
				std::max(1.0e-20,
					voice.ModePickupDisplacementAmplitude(1)));
		};
		Check(std::abs(pickupDisplacementRatioDb(f1ResidueVoice,
				tfdsp::ElectricPianoToneBarSubMode) + 28.0) < 0.15 &&
			std::abs(pickupDisplacementRatioDb(f3ResidueVoice,
				tfdsp::ElectricPianoToneBarSubMode) + 4.5) < 0.15,
			"electric piano reciprocal tone-bar submode matches measured F1/F3 magnitude");
		Check(upperSubModeVoice.ModeFrequencyRatio(
				tfdsp::ElectricPianoToneBarSubMode) > 0.42 &&
			upperSubModeVoice.ModeFrequencyRatio(
				tfdsp::ElectricPianoToneBarSubMode) < 0.45 &&
			pickupDisplacementRatioDb(upperSubModeVoice,
				tfdsp::ElectricPianoToneBarSubMode) < -34.0 &&
			upperSubModeVoice.ModeAmplitudeLifetimeSeconds(
				tfdsp::ElectricPianoToneBarSubMode) < 0.055,
			"electric piano provisional upper tone-bar submode cannot restore the persistent false sideband");

		// Contact geometry is independent of playing velocity. The nonlinear force
		// duration supplies velocity-dependent bandwidth; an earlier trial deriving
		// width from indentation assumed an unsupported spherical tip profile.
		tfdsp::ElectricPianoVoice softProjectionVoice;
		tfdsp::ElectricPianoVoice hardProjectionVoice;
		softProjectionVoice.SetSampleRate(48000.0);
		hardProjectionVoice.SetSampleRate(48000.0);
		softProjectionVoice.Step(0.0, 10.0, 0.2, false, controls);
		hardProjectionVoice.Step(0.0, 10.0, 1.0, false, controls);
		bool velocityIndependentProjection = true;
		bool finiteContactProjectsHigherModes = false;
		for (std::size_t mode = tfdsp::ElectricPianoAttackModeBegin;
			mode < tfdsp::ElectricPianoAttackModeEnd; ++mode)
		{
			const double softProjection =
				softProjectionVoice.ContactModeProjection(mode);
			const double hardProjection =
				hardProjectionVoice.ContactModeProjection(mode);
			velocityIndependentProjection = velocityIndependentProjection &&
				softProjection == hardProjection;
			finiteContactProjectsHigherModes = finiteContactProjectsHigherModes ||
				std::abs(softProjection) > 1.0e-3;
		}
		Check(velocityIndependentProjection && finiteContactProjectsHigherModes &&
			softProjectionVoice.ContactWidthMetres() >= 0.006 &&
			softProjectionVoice.ContactWidthMetres() <= 0.012 &&
			softProjectionVoice.ContactWidthMetres() ==
				hardProjectionVoice.ContactWidthMetres(),
			"electric piano finite contact footprint is geometric and velocity-independent");
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
		// minimum nonzero velocity. It must not start a long, inaudible oversampled
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
		double maximumContactRateSpan = 0.0;
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
							std::ceil(0.080 * sampleRate));
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
						// At 8 kHz the deliberately bandlimited beam contains fewer
						// modes and therefore has a genuinely different contact-point
						// impedance. Compare integration convergence only across normal
						// audio rates; the 8 kHz case is still required to finish safely.
						if (sampleRate >= 44100.0)
						{
							minimumDuration = std::min(minimumDuration, duration);
							maximumDuration = std::max(maximumDuration, duration);
						}
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
					maximumContactRateSpan = std::max(maximumContactRateSpan,
						maximumDuration - minimumDuration);
				}
			}
		}
		if (!(allContactsComplete && contactDurationsAreRateInvariant &&
			longestContact < 0.060))
			std::cerr << "electric piano longest physical contact " <<
				longestContact << " seconds, rate span " <<
				maximumContactRateSpan << " seconds\n";
		Check(allContactsComplete && contactDurationsAreRateInvariant &&
			longestContact < 0.060,
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
		double pickupCrossProduct = 0.0;
		for (int sample = 0; sample < 4800; ++sample)
		{
			const double farOutput = farPickupVoice.Step(0.0, 10.0, 0.3,
				false, farPickupControls);
			const double nearOutput = nearPickupVoice.Step(0.0, 10.0, 0.3,
				false, nearPickupControls);
			farPickupEnergy += farOutput * farOutput;
			nearPickupEnergy += nearOutput * nearOutput;
			pickupCrossProduct += farOutput * nearOutput;
		}
		const double proximityEnergyRatio = nearPickupEnergy /
			std::max(1.0e-20, farPickupEnergy);
		const double proximityResidual = std::sqrt(std::max(0.0,
			nearPickupEnergy - pickupCrossProduct * pickupCrossProduct /
				std::max(1.0e-20, farPickupEnergy)) /
			std::max(1.0e-20, nearPickupEnergy));
		if (!(proximityEnergyRatio > 0.70 && proximityEnergyRatio < 1.55 &&
			proximityResidual > 0.025))
			std::cerr << "electric piano pickup proximity energy ratio " <<
				proximityEnergyRatio << " residual " << proximityResidual << '\n';
		Check(proximityEnergyRatio > 0.70 && proximityEnergyRatio < 1.55 &&
			proximityResidual > 0.025,
			"electric piano pickup Proximity changes curvature without becoming a second Drive");

		// The old additive per-key gap offset made full Proximity cut bass notes but
		// boost treble notes by roughly 6 dB. The trajectory normalization must keep
		// the close setting within a narrow level window in every register; its
		// useful effect is magnetic curvature, checked above, rather than loudness.
		for (double pitch : {-8.0 / 3.0, -2.0, -1.0, 0.0, 1.0, 2.0,
			10.0 / 3.0})
		{
			tfdsp::ElectricPianoVoice defaultGapVoice;
			tfdsp::ElectricPianoVoice closeGapVoice;
			defaultGapVoice.SetSampleRate(48000.0);
			closeGapVoice.SetSampleRate(48000.0);
			auto defaultGapControls = silentControls;
			auto closeGapControls = silentControls;
			closeGapControls.proximity = 1.0;
			double defaultGapEnergy = 0.0;
			double closeGapEnergy = 0.0;
			for (int sample = 0; sample < 4800; ++sample)
			{
				const double ordinary = defaultGapVoice.Step(pitch, 10.0, 0.86,
					false, defaultGapControls);
				const double close = closeGapVoice.Step(pitch, 10.0, 0.86,
					false, closeGapControls);
				defaultGapEnergy += ordinary * ordinary;
				closeGapEnergy += close * close;
			}
			const double closeLevelDb = 10.0 * std::log10(closeGapEnergy /
				std::max(1.0e-20, defaultGapEnergy));
			Check(closeLevelDb > -1.3 && closeLevelDb < 1.5,
				"electric piano Proximity level remains calibrated across the keyboard");
		}

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

		// Rack MIDI-to-CV keeps a stolen polyphonic channel's gate high and
		// reports the replacement note on its separate retrigger output. The model
		// must therefore accept an explicit strike event while gate remains high.
		tfdsp::ElectricPianoVoice heldGateRetriggerVoice;
		tfdsp::ElectricPianoVoice heldGateReferenceVoice;
		heldGateRetriggerVoice.SetSampleRate(48000.0);
		heldGateReferenceVoice.SetSampleRate(48000.0);
		heldGateRetriggerVoice.SetNoiseSeed(0x13579bdfu);
		heldGateReferenceVoice.SetNoiseSeed(0x13579bdfu);
		for (int sample = 0; sample < 1200; ++sample)
		{
			heldGateRetriggerVoice.Step(0.0, 10.0, 0.8, false,
				silentControls);
			heldGateReferenceVoice.Step(0.0, 10.0, 0.8, false,
				silentControls);
		}
		double heldGateRetriggerDifference = 0.0;
		for (int sample = 0; sample < 2048; ++sample)
		{
			const bool retrigger = sample == 0;
			const double struck = heldGateRetriggerVoice.Step(0.0, 10.0, 0.65,
				false, silentControls, {}, retrigger);
			const double reference = heldGateReferenceVoice.Step(0.0, 10.0, 0.65,
				false, silentControls);
			heldGateRetriggerDifference += (struck - reference) *
				(struck - reference);
		}
		Check(heldGateRetriggerDifference > 1.0e-6,
			"electric piano explicit retrigger strikes while the note gate remains high");

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

		// The default modulation object must take the same fast path exactly;
		// disconnected modulation must not move or colour normal piano output.
		tfdsp::ElectricPianoVoice implicitBypassVoice;
		tfdsp::ElectricPianoVoice explicitBypassVoice;
		implicitBypassVoice.SetSampleRate(48000.0);
		explicitBypassVoice.SetSampleRate(48000.0);
		implicitBypassVoice.SetNoiseSeed(0x31415926u);
		explicitBypassVoice.SetNoiseSeed(0x31415926u);
		bool modulationBypassExact = true;
		for (int sample = 0; sample < 4096; ++sample)
		{
			const double implicit = implicitBypassVoice.Step(0.0, 10.0, 0.8,
				false, quietControls);
			const double explicitZero = explicitBypassVoice.Step(0.0, 10.0, 0.8,
				false, quietControls, tfdsp::ElectricPianoModulation{});
			modulationBypassExact = modulationBypassExact &&
				implicit == explicitZero;
		}
		Check(modulationBypassExact,
			"electric piano disconnected modulation is sample-exact bypass");

		tfdsp::ElectricPianoVoice exponentialFmVoice;
		tfdsp::ElectricPianoVoice linearFmVoice;
		tfdsp::ElectricPianoVoice phaseModVoice;
		tfdsp::ElectricPianoVoice modulationReferenceVoice;
		for (auto* voice : {&exponentialFmVoice, &linearFmVoice,
			&phaseModVoice, &modulationReferenceVoice})
		{
			voice->SetSampleRate(48000.0);
			voice->SetNoiseSeed(0x27182818u);
		}
		double exponentialDifference = 0.0;
		double linearDifference = 0.0;
		double phaseDifference = 0.0;
		bool creativeModulationFinite = true;
		for (int sample = 0; sample < 8192; ++sample)
		{
			const double lfo = std::sin(2.0 * tfdsp::PI * 137.0 * sample /
				48000.0);
			tfdsp::ElectricPianoModulation exponential;
			exponential.exponentialPitch = 0.18 * lfo;
			tfdsp::ElectricPianoModulation linear;
			linear.linearFrequencyRatio = 0.42 * lfo;
			tfdsp::ElectricPianoModulation phase;
			phase.phaseRadians = 0.75 * lfo;
			const double reference = modulationReferenceVoice.Step(0.0, 10.0,
				0.8, false, quietControls);
			const double exponentialOutput = exponentialFmVoice.Step(0.0, 10.0,
				0.8, false, quietControls, exponential);
			const double linearOutput = linearFmVoice.Step(0.0, 10.0, 0.8,
				false, quietControls, linear);
			const double phaseOutput = phaseModVoice.Step(0.0, 10.0, 0.8,
				false, quietControls, phase);
			exponentialDifference += (exponentialOutput - reference) *
				(exponentialOutput - reference);
			linearDifference += (linearOutput - reference) *
				(linearOutput - reference);
			phaseDifference += (phaseOutput - reference) *
				(phaseOutput - reference);
			creativeModulationFinite = creativeModulationFinite &&
				std::isfinite(exponentialOutput) && std::isfinite(linearOutput) &&
				std::isfinite(phaseOutput);
		}
		Check(creativeModulationFinite && exponentialDifference > 1.0e-5 &&
			linearDifference > 1.0e-5 && phaseDifference > 1.0e-5,
			"electric piano EXP FM, linear TZ FM, and phase modulation are finite and audible");

		// At -100% linear deviation the three coupled body coordinates are
		// stationary. The deliberately protected inharmonic attack modes retain a
		// small residual, while going farther reverses the body through zero.
		tfdsp::ElectricPianoVoice zeroFrequencyVoice;
		tfdsp::ElectricPianoVoice reverseFrequencyVoice;
		zeroFrequencyVoice.SetSampleRate(48000.0);
		reverseFrequencyVoice.SetSampleRate(48000.0);
		zeroFrequencyVoice.SetNoiseSeed(0xabcdef01u);
		reverseFrequencyVoice.SetNoiseSeed(0xabcdef01u);
		for (int sample = 0; sample < 48000; ++sample)
		{
			zeroFrequencyVoice.Step(0.0, 10.0, 0.8, false, quietControls);
			reverseFrequencyVoice.Step(0.0, 10.0, 0.8, false, quietControls);
		}
		double zeroFrequencyEnergy = 0.0;
		double reverseFrequencyEnergy = 0.0;
		for (int sample = 0; sample < 4096; ++sample)
		{
			tfdsp::ElectricPianoModulation zeroFrequency;
			zeroFrequency.linearFrequencyRatio = -1.0;
			tfdsp::ElectricPianoModulation reverseFrequency;
			reverseFrequency.linearFrequencyRatio = -1.5;
			const double zeroOutput = zeroFrequencyVoice.Step(0.0, 10.0, 0.8,
				false, quietControls, zeroFrequency);
			const double reverseOutput = reverseFrequencyVoice.Step(0.0, 10.0,
				0.8, false, quietControls, reverseFrequency);
			if (sample >= 3072)
			{
				zeroFrequencyEnergy += zeroOutput * zeroOutput;
				reverseFrequencyEnergy += reverseOutput * reverseOutput;
			}
		}
		Check(zeroFrequencyEnergy < 0.01 * reverseFrequencyEnergy &&
			reverseFrequencyEnergy > 1.0e-10,
			"electric piano linear FM crosses continuously through its zero-frequency body state");

		// Constant linear FM repeatedly wraps every modal phase accumulator. A
		// former shared wrapped phase was multiplied by non-integer modal ratios,
		// creating isolated output steps at every wrap. Require a normal bounded
		// derivative crest under both repeated wraps and low-rate PM.
		auto modulationDerivativeCrest = [&](bool phaseModulation)
		{
			tfdsp::ElectricPianoVoice voice;
			voice.SetSampleRate(48000.0);
			voice.SetNoiseSeed(0x31415926u);
			double previous = 0.0;
			double maximumDifference = 0.0;
			double differenceEnergy = 0.0;
			int measuredSamples = 0;
			for (int sample = 0; sample < 48000; ++sample)
			{
				tfdsp::ElectricPianoModulation mod;
				if (phaseModulation)
					mod.phaseRadians = 0.10 * std::sin(2.0 * tfdsp::PI *
						5.0 * sample / 48000.0);
				else
					mod.linearFrequencyRatio = 0.37;
				const double output = voice.Step(-1.0, 10.0, 0.82, false,
					quietControls, mod);
				if (sample >= 4096)
				{
					const double difference = output - previous;
					maximumDifference = std::max(maximumDifference,
						std::abs(difference));
					differenceEnergy += difference * difference;
					++measuredSamples;
				}
				previous = output;
			}
			return maximumDifference / std::max(1.0e-12,
				std::sqrt(differenceEnergy / measuredSamples));
		};
		const double wrappedFmCrest = modulationDerivativeCrest(false);
		const double slowPmCrest = modulationDerivativeCrest(true);
		if (!(wrappedFmCrest < 20.0 && slowPmCrest < 20.0))
			std::cerr << "electric piano modulation derivative crests " <<
				wrappedFmCrest << "/" << slowPmCrest << '\n';
		Check(wrappedFmCrest < 20.0 && slowPmCrest < 20.0,
			"electric piano FM phase wraps and low-rate PM remain continuous");

		auto nearStrikeControls = quietControls;
		auto farStrikeControls = quietControls;
		nearStrikeControls.strikePosition = 0.0;
		farStrikeControls.strikePosition = 1.0;
		tfdsp::ElectricPianoVoice nearStrikeVoice;
		tfdsp::ElectricPianoVoice middleStrikeVoice;
		tfdsp::ElectricPianoVoice farStrikeVoice;
		nearStrikeVoice.SetSampleRate(48000.0);
		middleStrikeVoice.SetSampleRate(48000.0);
		farStrikeVoice.SetSampleRate(48000.0);
		nearStrikeVoice.SetNoiseSeed(0x10203040u);
		middleStrikeVoice.SetNoiseSeed(0x10203040u);
		farStrikeVoice.SetNoiseSeed(0x10203040u);
		double nearStrikeEnergy = 0.0;
		double middleStrikeEnergy = 0.0;
		double farStrikeEnergy = 0.0;
		double nearMiddleProduct = 0.0;
		double farMiddleProduct = 0.0;
		for (int sample = 0; sample < 4096; ++sample)
		{
			const double nearOutput = nearStrikeVoice.Step(0.0, 10.0, 0.85,
				false, nearStrikeControls);
			const double middleOutput = middleStrikeVoice.Step(0.0, 10.0, 0.85,
				false, quietControls);
			const double farOutput = farStrikeVoice.Step(0.0, 10.0, 0.85,
				false, farStrikeControls);
			nearStrikeEnergy += nearOutput * nearOutput;
			middleStrikeEnergy += middleOutput * middleOutput;
			farStrikeEnergy += farOutput * farOutput;
			nearMiddleProduct += nearOutput * middleOutput;
			farMiddleProduct += farOutput * middleOutput;
		}
		auto strikeLevelMatchedResidual = [&](double energy, double product)
		{
			const double residual = std::max(0.0, energy - product * product /
				std::max(1.0e-20, middleStrikeEnergy));
			return std::sqrt(residual / std::max(1.0e-20, energy));
		};
		const double nearStrikeResidual = strikeLevelMatchedResidual(
			nearStrikeEnergy, nearMiddleProduct);
		const double farStrikeResidual = strikeLevelMatchedResidual(
			farStrikeEnergy, farMiddleProduct);
		if (!(nearStrikeResidual > 0.04 && farStrikeResidual > 0.09))
			std::cerr << "electric piano strike residuals min/max versus default " <<
				nearStrikeResidual << "/" << farStrikeResidual << '\n';
		Check(nearStrikeResidual > 0.04 && farStrikeResidual > 0.09,
			"electric piano Strike has a perceptually useful range around its calibrated midpoint");

		bool strikeLawContinuous = true;
		double previousFactoryPosition = 0.0;
		for (int midi = 28; midi <= 100; ++midi)
		{
			auto leftControls = quietControls;
			auto centreControls = quietControls;
			auto rightControls = quietControls;
			leftControls.strikePosition = 0.49;
			centreControls.strikePosition = 0.50;
			rightControls.strikePosition = 0.51;
			tfdsp::ElectricPianoVoice leftVoice;
			tfdsp::ElectricPianoVoice centreVoice;
			tfdsp::ElectricPianoVoice rightVoice;
			leftVoice.SetSampleRate(48000.0);
			centreVoice.SetSampleRate(48000.0);
			rightVoice.SetSampleRate(48000.0);
			const double pitch = (static_cast<double>(midi) - 60.0) / 12.0;
			leftVoice.Step(pitch, 10.0, 0.8, false, leftControls);
			centreVoice.Step(pitch, 10.0, 0.8, false, centreControls);
			rightVoice.Step(pitch, 10.0, 0.8, false, rightControls);
			const double leftDelta = centreVoice.StrikePosition() -
				leftVoice.StrikePosition();
			const double rightDelta = rightVoice.StrikePosition() -
				centreVoice.StrikePosition();
			const auto key = tfdsp::MakeElectricPianoKeyProfile(pitch);
			const double physicalDelta = 0.5 * (leftDelta + rightDelta) *
				key.tineLengthMetres;
			strikeLawContinuous = strikeLawContinuous && leftDelta > 0.0 &&
				rightDelta > 0.0 &&
				std::abs(leftDelta - rightDelta) < 1.0e-12 &&
				physicalDelta >= 0.000019 && physicalDelta <= 0.000121;
			if (midi > 28)
				strikeLawContinuous = strikeLawContinuous &&
					std::abs(centreVoice.StrikePosition() -
						previousFactoryPosition) < 0.025;
			previousFactoryPosition = centreVoice.StrikePosition();
		}
		Check(strikeLawContinuous,
			"electric piano Strike is symmetric and continuous across all 73 keys");

		tfdsp::ElectricPianoVoice heldStrikeEditVoice;
		tfdsp::ElectricPianoVoice heldStrikeReferenceVoice;
		heldStrikeEditVoice.SetSampleRate(48000.0);
		heldStrikeReferenceVoice.SetSampleRate(48000.0);
		heldStrikeEditVoice.SetNoiseSeed(0x50607080u);
		heldStrikeReferenceVoice.SetNoiseSeed(0x50607080u);
		bool heldStrikeEditExact = true;
		for (int sample = 0; sample < 4096; ++sample)
		{
			auto editedStrikeControls = quietControls;
			if (sample > 0)
				editedStrikeControls.strikePosition = 1.0;
			const double edited = heldStrikeEditVoice.Step(0.0, 10.0, 0.85,
				false, editedStrikeControls);
			const double reference = heldStrikeReferenceVoice.Step(0.0, 10.0,
				0.85, false, quietControls);
			heldStrikeEditExact = heldStrikeEditExact && edited == reference;
		}
		Check(heldStrikeEditExact,
			"electric piano strike position is sampled only on key-down");

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

		// Bell is a residue balance, not an output-level control. Its range grows
		// toward the upper keyboard where progressively fewer mechanical attack
		// modes remain below pickup bandwidth. Guard a clearly different, level-
		// matched onset at C4, E5 and E6 without requiring a large RMS change.
		constexpr std::array<int, 3> BellKeys{60, 76, 88};
		constexpr std::array<double, 3> MinimumBellResiduals{0.25, 0.22, 0.32};
		for (std::size_t key = 0; key < BellKeys.size(); ++key)
		{
			auto lowBellControls = quietControls;
			auto highBellControls = quietControls;
			lowBellControls.bell = 0.0;
			highBellControls.bell = 1.0;
			tfdsp::ElectricPianoVoice lowBellVoice;
			tfdsp::ElectricPianoVoice highBellVoice;
			lowBellVoice.SetSampleRate(48000.0);
			highBellVoice.SetSampleRate(48000.0);
			std::vector<double> lowBellSignal(2048);
			std::vector<double> highBellSignal(2048);
			const double pitch = (static_cast<double>(BellKeys[key]) - 60.0) /
				12.0;
			for (std::size_t sample = 0; sample < lowBellSignal.size(); ++sample)
			{
				lowBellSignal[sample] = lowBellVoice.Step(pitch, 10.0, 0.90,
					false, lowBellControls);
				highBellSignal[sample] = highBellVoice.Step(pitch, 10.0, 0.90,
					false, highBellControls);
			}
			const double residual = LevelMatchedResidual(lowBellSignal,
				highBellSignal, 0, 960);
			double lowBellEnergy = 0.0;
			double highBellEnergy = 0.0;
			for (std::size_t sample = 0; sample < 960; ++sample)
			{
				lowBellEnergy += lowBellSignal[sample] * lowBellSignal[sample];
				highBellEnergy += highBellSignal[sample] * highBellSignal[sample];
			}
			const double bellLevelChangeDb = 10.0 * std::log10(
				std::max(1.0e-20, highBellEnergy) /
				std::max(1.0e-20, lowBellEnergy));
			if (!(residual > MinimumBellResiduals[key] &&
				std::abs(bellLevelChangeDb) < 1.8))
				std::cerr << "electric piano Bell onset residual key " <<
					BellKeys[key] << " " << residual << " level " <<
					bellLevelChangeDb << " dB\n";
			Check(residual > MinimumBellResiduals[key] &&
				std::abs(bellLevelChangeDb) < 1.8,
				"electric piano Bell has an audible level-matched onset range without becoming another VCA");
		}
		double maximumBellInteractionLevelDb = 0.0;
		bool bellInteractionsFinite = true;
		for (int midi : {76, 88})
			for (double tone : {0.0, 1.0})
				for (double proximity : {0.0, 1.0})
				{
					auto lowBellControls = quietControls;
					auto highBellControls = quietControls;
					lowBellControls.bell = 0.0;
					highBellControls.bell = 1.0;
					lowBellControls.tone = highBellControls.tone = tone;
					lowBellControls.proximity = highBellControls.proximity =
						proximity;
					tfdsp::ElectricPianoVoice lowBellVoice;
					tfdsp::ElectricPianoVoice highBellVoice;
					lowBellVoice.SetSampleRate(48000.0);
					highBellVoice.SetSampleRate(48000.0);
					double lowEnergy = 0.0;
					double highEnergy = 0.0;
					const double pitch = (static_cast<double>(midi) - 60.0) / 12.0;
					for (int sample = 0; sample < 960; ++sample)
					{
						const double low = lowBellVoice.Step(pitch, 10.0, 0.90,
							false, lowBellControls);
						const double high = highBellVoice.Step(pitch, 10.0, 0.90,
							false, highBellControls);
						bellInteractionsFinite = bellInteractionsFinite &&
							std::isfinite(low) && std::isfinite(high);
						lowEnergy += low * low;
						highEnergy += high * high;
					}
					maximumBellInteractionLevelDb = std::max(
						maximumBellInteractionLevelDb, std::abs(10.0 * std::log10(
							std::max(1.0e-20, highEnergy) /
							std::max(1.0e-20, lowEnergy))));
				}
		if (!(bellInteractionsFinite && maximumBellInteractionLevelDb < 3.2))
			std::cerr << "electric piano Bell interaction maximum level " <<
				maximumBellInteractionLevelDb << " dB\n";
		Check(bellInteractionsFinite && maximumBellInteractionLevelDb < 3.2,
			"electric piano Bell remains bounded at extreme Tone and Proximity settings");

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
		if (!(couplingAttack[0] > 0.12 && couplingSustain[0] > 0.22 &&
			couplingAttack[1] > -6.0 && couplingAttack[1] < 1.0 &&
			couplingSustain[1] > couplingAttack[1] + 3.0))
			std::cerr << "electric piano coupling attack/sustain residual/dB " <<
				couplingAttack[0] << "/" << couplingAttack[1] << " " <<
				couplingSustain[0] << "/" << couplingSustain[1] << '\n';
		Check(couplingAttack[0] > 0.12 && couplingSustain[0] > 0.22 &&
			couplingAttack[1] > -6.0 && couplingAttack[1] < 1.0 &&
			couplingSustain[1] > couplingAttack[1] + 3.0,
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
		if (!(sweptCouplingDifference > 0.05 * fixedCouplingEnergy &&
			duringSweepStep < 1.5 * preSweepStep))
			std::cerr << "electric piano live Coupling difference ratio " <<
				sweptCouplingDifference / std::max(1.0e-20, fixedCouplingEnergy) <<
				" step ratio " << duringSweepStep /
				std::max(1.0e-20, preSweepStep) << '\n';
		Check(sweptCouplingDifference > 0.05 * fixedCouplingEnergy &&
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
		if (!(maximumKeyboardEnergy < 3.0 * minimumKeyboardEnergy))
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
		Check(maximumKeyboardEnergy < 3.0 * minimumKeyboardEnergy,
			"electric piano reference pickup curve stays balanced while retaining longer bass sustain");

		// Real pickups are adjusted one key at a time. Verify the calibration at the
		// same four-key checkpoints used by the model, but average soft, medium and
		// hard strikes geometrically so the test cannot flatten the physical
		// register-dependent velocity response. The intended -0.42 dB/octave voltage
		// tilt spans about 2.5 dB over all 73 keys and remains smooth between screws.
		std::array<double, 19> calibratedKeyLevels{};
		for (std::size_t key = 0; key < calibratedKeyLevels.size(); ++key)
		{
			const double pitch = (-32.0 + 4.0 * static_cast<double>(key)) / 12.0;
			double levelProduct = 1.0;
			for (double velocity : {0.30, 0.60, 0.90})
			{
				tfdsp::ElectricPianoVoice voice;
				voice.SetSampleRate(48000.0);
				double energy = 0.0;
				for (int sample = 0; sample < 4800; ++sample)
				{
					const double output = voice.Step(pitch, 10.0, velocity,
						false, silentControls);
					energy += output * output;
				}
				levelProduct *= std::sqrt(energy / 4800.0);
			}
			calibratedKeyLevels[key] = std::cbrt(levelProduct);
		}
		const auto calibratedLevelLimits = std::minmax_element(
			calibratedKeyLevels.begin(), calibratedKeyLevels.end());
		double maximumCheckpointStepDb = 0.0;
		for (std::size_t key = 1; key < calibratedKeyLevels.size(); ++key)
			maximumCheckpointStepDb = std::max(maximumCheckpointStepDb,
				std::abs(20.0 * std::log10(calibratedKeyLevels[key] /
					calibratedKeyLevels[key - 1])));
		Check(*calibratedLevelLimits.second <
			1.38 * *calibratedLevelLimits.first && maximumCheckpointStepDb < 0.40,
			"electric piano per-key pickup voicing is smooth across the 73-key harp and multiple velocities");

		// Track lower-register attack level, decay continuity, low harmonics and the
		// velocity-dependent H4--H10 ladder at four separated keys. These bounds are
		// intentionally relaxed during the fork-topology correction so the previous
		// render cannot force incorrect physics; tighten them again after the next
		// listening approval rather than silently accepting later drift.
		for (int midi : {29, 40, 48, 53})
		{
			const double pitch = (static_cast<double>(midi) - 60.0) / 12.0;
			const double fundamental = tfdsp::ElectricPianoReferenceFrequency *
				std::exp2(pitch) / 48000.0;
			std::array<std::vector<double>, 2> signals{
				std::vector<double>(12288), std::vector<double>(12288)};
			for (std::size_t velocity = 0; velocity < signals.size(); ++velocity)
			{
				tfdsp::ElectricPianoVoice voice;
				voice.SetSampleRate(48000.0);
				for (double& output : signals[velocity])
					output = voice.Step(pitch, 10.0,
						velocity == 0 ? 0.30 : 0.86, false, silentControls);
			}
			auto rms = [](const std::vector<double>& signal, std::size_t begin,
				std::size_t end)
			{
				double energy = 0.0;
				for (std::size_t sample = begin; sample < end; ++sample)
					energy += signal[sample] * signal[sample];
				return std::sqrt(energy / static_cast<double>(end - begin));
			};
			auto harmonicRatio = [&](const std::vector<double>& signal, int harmonic)
			{
				return WindowedHarmonicMagnitude(signal,
					harmonic * fundamental, 6000, signal.size()) /
					std::max(1.0e-12, WindowedHarmonicMagnitude(signal,
						fundamental, 6000, signal.size()));
			};
			auto upperHarmonics = [&](const std::vector<double>& signal)
			{
				double energy = 0.0;
				for (int harmonic = 4; harmonic <= 10; ++harmonic)
				{
					const double ratio = harmonicRatio(signal, harmonic);
					energy += ratio * ratio;
				}
				return std::sqrt(energy);
			};
			const double hardRms = rms(signals[1], 0, 4096);
			const double sustainRatio = rms(signals[1], 4096, 12288) /
				std::max(1.0e-12, hardRms);
			const double hardSecond = harmonicRatio(signals[1], 2);
			const double hardThird = harmonicRatio(signals[1], 3);
			const double softUpper = upperHarmonics(signals[0]);
			const double hardUpper = upperHarmonics(signals[1]);
			if (!(hardRms > 0.0095 && hardRms < 0.0195 &&
				sustainRatio > 0.80 && sustainRatio < 1.12 &&
				hardSecond > 0.10 && hardSecond < 0.40 &&
				hardThird > 0.30 && hardThird < 1.10 && hardUpper > 0.18 &&
				hardUpper > 1.8 * softUpper))
				std::cerr << "electric piano low-register fingerprint key " << midi <<
					" rms/sustain/H2/H3/upper soft-hard " << hardRms << '/' <<
					sustainRatio << '/' << hardSecond << '/' << hardThird << '/' <<
					softUpper << '-' << hardUpper << '\n';
			Check(hardRms > 0.0095 && hardRms < 0.0195 &&
				sustainRatio > 0.80 && sustainRatio < 1.12 &&
				hardSecond > 0.10 && hardSecond < 0.40 &&
				hardThird > 0.30 && hardThird < 1.10 && hardUpper > 0.18 &&
				hardUpper > 1.8 * softUpper,
				"electric piano lower-register level, decay and velocity bark retain their calibrated character");
		}

		// The upper register must not expose one otherwise valid inharmonic attack
		// line against an empty spectrum. The closer per-key neutral pickup setting
		// should generate a coherent integer ladder before the signal reaches the
		// amplifier, while the short mechanical line remains subordinate.
		constexpr std::array<int, 3> TrebleKeys{84, 96, 100};
		constexpr std::array<double, 3> MinimumUpperLadders{0.065, 0.009, 0.0035};
		for (std::size_t keyIndex = 0; keyIndex < TrebleKeys.size(); ++keyIndex)
		{
			const double pitch = (static_cast<double>(TrebleKeys[keyIndex]) -
				60.0) / 12.0;
			const double fundamentalFrequency =
				tfdsp::ElectricPianoReferenceFrequency * std::exp2(pitch);
			tfdsp::ElectricPianoVoice voice;
			voice.SetSampleRate(48000.0);
			std::vector<double> signal(4096);
			for (double& output : signal)
				output = voice.Step(pitch, 10.0, 0.86, false, silentControls);
			const double fundamental = WindowedHarmonicMagnitude(signal,
				fundamentalFrequency / 48000.0, 0, 2048);
			double harmonicEnergy = 0.0;
			double upperHarmonicEnergy = 0.0;
			for (int harmonic = 2;
				harmonic * fundamentalFrequency < 22000.0; ++harmonic)
			{
				const double magnitude = WindowedHarmonicMagnitude(signal,
					harmonic * fundamentalFrequency / 48000.0, 0, 2048);
				harmonicEnergy += magnitude * magnitude;
				if (harmonic >= 4)
					upperHarmonicEnergy += magnitude * magnitude;
			}
			const double harmonicRatio = std::sqrt(harmonicEnergy) /
				std::max(1.0e-12, fundamental);
			const double upperRatio = std::sqrt(upperHarmonicEnergy) /
				std::max(1.0e-12, fundamental);
			const double attackModeRatio = WindowedHarmonicMagnitude(signal,
				voice.ModeFrequencyRatio(3) * fundamentalFrequency / 48000.0,
				0, 2048) / std::max(1.0e-12, fundamental);
			if (!(harmonicRatio > 0.15 &&
				upperRatio > MinimumUpperLadders[keyIndex] &&
				attackModeRatio < 0.10 * harmonicRatio))
				std::cerr << "electric piano treble harmonic bed key " <<
					TrebleKeys[keyIndex] << " all/upper/mode " << harmonicRatio <<
					'/' << upperRatio << '/' << attackModeRatio << '\n';
			Check(harmonicRatio > 0.15 &&
				upperRatio > MinimumUpperLadders[keyIndex] &&
				attackModeRatio < 0.10 * harmonicRatio,
				"electric piano treble pickup harmonic bed masks isolated inharmonic attack lines");
		}

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
		const double firstAttackRatio =
			defaultHammerVoice.ModeFrequencyRatio(3);
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
		auto upperHammerBrightness = [&](const std::vector<double>& signal)
		{
			double energy = 0.0;
			for (int harmonic = 20; harmonic <= 60; ++harmonic)
			{
				const double magnitude = WindowedHarmonicMagnitude(signal,
					harmonic * fundamentalBin, 0, 2048);
				energy += magnitude * magnitude;
			}
			return std::sqrt(energy) / std::max(1.0e-12,
				WindowedHarmonicMagnitude(signal, fundamentalBin, 0, 2048));
		};
		const double softHammerUpper = upperHammerBrightness(softHammerSignal);
		const double defaultHammerUpper = upperHammerBrightness(defaultHammerSignal);
		const double hardHammerUpper = upperHammerBrightness(hardHammerSignal);
		if (!(defaultHammerBrightness < 0.06 &&
			hardHammerBrightness < 0.06))
			std::cerr << "electric piano Hammer attack ratios soft/default/hard " <<
				softHammerBrightness << "/" << defaultHammerBrightness << "/" <<
				hardHammerBrightness << '\n';
		if (!(defaultHammerUpper > 1.30 * softHammerUpper &&
			hardHammerUpper > 1.55 * softHammerUpper))
			std::cerr << "electric piano Hammer upper ratios soft/default/hard " <<
				softHammerUpper << "/" << defaultHammerUpper << "/" <<
				hardHammerUpper << '\n';
		Check(defaultHammerUpper > 1.30 * softHammerUpper &&
			hardHammerUpper > 1.55 * softHammerUpper,
			"electric piano hammer hardness changes compliant-contact brightness");
		Check(defaultHammerBrightness < 0.06 && hardHammerBrightness < 0.06 &&
			hardHammerUpper < 0.02,
			"electric piano first attack mode remains non-metallic at every Hammer setting");

		double minimumKeyboardAttack = std::numeric_limits<double>::infinity();
		double maximumKeyboardAttack = 0.0;
		for (double pitch : {-2.0, -1.0, 0.0, 1.0})
		{
			const double fundamental = tfdsp::ElectricPianoReferenceFrequency *
				std::exp2(pitch) / 48000.0;
			tfdsp::ElectricPianoVoice voice;
			voice.SetSampleRate(48000.0);
			std::vector<double> signal(2048);
			for (double& output : signal)
				output = voice.Step(pitch, 10.0, 0.72, false, silentControls);
			const double attack = voice.ModeFrequencyRatio(3) * fundamental;
			const double brightness = HarmonicMagnitude(signal, attack) /
				std::max(1.0e-12, HarmonicMagnitude(signal, fundamental));
			minimumKeyboardAttack = std::min(minimumKeyboardAttack, brightness);
			maximumKeyboardAttack = std::max(maximumKeyboardAttack, brightness);
		}
		if (!(maximumKeyboardAttack > 0.01 && maximumKeyboardAttack < 0.35))
			std::cerr << "electric piano keyboard attack range " <<
				minimumKeyboardAttack << "/" << maximumKeyboardAttack << '\n';
		Check(maximumKeyboardAttack > 0.01 && maximumKeyboardAttack < 0.35,
			"electric piano keyed strike line remains present without dominating");

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
		if (!(highSecond > 1.6 * lowSecond && highSecond > 0.16 &&
			highSecond < 0.28 && highThird > 6.0 * lowThird &&
			highThird > 0.30 && highThird < 0.65))
			std::cerr << "electric piano pickup ratios low/high H2 " <<
				lowSecond << "/" << highSecond << " H3 " << lowThird <<
				"/" << highThird << '\n';
		Check(highSecond > 1.6 * lowSecond && highSecond > 0.16 &&
			highSecond < 0.28 && highThird > 6.0 * lowThird &&
			highThird > 0.30 && highThird < 0.65,
			"electric piano asymmetric pickup develops progressive low-harmonic bark with velocity");
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
		if (!(std::is_sorted(secondHarmonics.begin(), secondHarmonics.end()) &&
			secondHarmonics.back() > 1.35 * secondHarmonics.front() &&
			std::is_sorted(thirdHarmonics.begin(), thirdHarmonics.end()) &&
			thirdHarmonics.back() > 4.0 * thirdHarmonics.front() &&
			*std::max_element(attackRatios.begin(), attackRatios.end()) > 0.003 &&
			*std::max_element(attackRatios.begin(), attackRatios.end()) < 0.06))
		{
			std::cerr << "electric piano velocity H2/H3/attack";
			for (std::size_t index = 0; index < TimbreVelocities.size(); ++index)
				std::cerr << " [" << secondHarmonics[index] << ',' <<
					thirdHarmonics[index] << ',' << attackRatios[index] << ']';
			std::cerr << '\n';
		}
		Check(std::is_sorted(secondHarmonics.begin(), secondHarmonics.end()) &&
			secondHarmonics.back() > 1.35 * secondHarmonics.front() &&
			std::is_sorted(thirdHarmonics.begin(), thirdHarmonics.end()) &&
			thirdHarmonics.back() > 4.0 * thirdHarmonics.front() &&
			*std::max_element(attackRatios.begin(), attackRatios.end()) > 0.003 &&
			*std::max_element(attackRatios.begin(), attackRatios.end()) < 0.06,
			"electric piano pickup harmonic timbre grows across playable velocities");

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
		const double toneSecondLogSpan = std::abs(std::log(
			brightSecond / std::max(1.0e-12, darkSecond)));
		if (!(toneSecondLogSpan > 0.35 && brightThird > 0.020 &&
			darkThird > 0.020))
			std::cerr << "electric piano Tone ratios attack dark/bright " <<
				darkFirst << "/" << brightFirst << " H2 " << darkSecond <<
				"/" << brightSecond << " H3 " << darkThird << "/" <<
				brightThird << '\n';
		Check(toneSecondLogSpan > 0.35 && brightThird > 0.020 && darkThird > 0.020,
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
		if (!(amplifiedToneResidual > 0.20))
			std::cerr << "electric piano amplified Tone residual " <<
				amplifiedToneResidual << '\n';
		Check(amplifiedToneResidual > 0.20,
			"electric piano Tone extremes remain distinct through the shared amplifier");
		const double amplifiedHammerResidual = levelMatchedResidual(
			softHammerAmplified, hardHammerAmplified, 2048);
		if (!(amplifiedHammerResidual > 0.25))
			std::cerr << "electric piano amplified Hammer residual " <<
				amplifiedHammerResidual << '\n';
		Check(amplifiedHammerResidual > 0.25,
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
		if (!(sweptDifferenceEnergy > 0.06 * fixedToneEnergy))
			std::cerr << "electric piano active Tone sweep energy ratio " <<
				sweptDifferenceEnergy / std::max(1.0e-20, fixedToneEnergy) << '\n';
		Check(sweptDifferenceEnergy > 0.06 * fixedToneEnergy,
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
		std::array<double, 3> defaultPlayedModeLifetimes{};
		for (std::size_t key = 0; key < defaultPlayedModeLifetimes.size(); ++key)
		{
			const std::array<double, 3> pitches{-8.0 / 3.0, 0.0, 10.0 / 3.0};
			tfdsp::ElectricPianoVoice lifetimeVoice;
			lifetimeVoice.SetSampleRate(48000.0);
			lifetimeVoice.Step(pitches[key], 10.0, 0.8, false, silentControls);
			defaultPlayedModeLifetimes[key] =
				lifetimeVoice.ModeAmplitudeLifetimeSeconds(1);
		}
		Check(defaultPlayedModeLifetimes[0] > 2.5 &&
			defaultPlayedModeLifetimes[0] < 4.0 &&
			defaultPlayedModeLifetimes[1] > 1.0 &&
			defaultPlayedModeLifetimes[1] < 1.3 &&
			defaultPlayedModeLifetimes[2] > 0.75 &&
			defaultPlayedModeLifetimes[2] < 1.0 &&
			std::is_sorted(defaultPlayedModeLifetimes.rbegin(),
				defaultPlayedModeLifetimes.rend()),
			"electric piano neutral Decay retains long bass sustain and progressively shorter upper-key lifetimes");

		const double defaultDamperRelease =
			tfdsp::ElectricPianoDamperReleaseSeconds(
				tfdsp::ElectricPianoControls{}.release);
		Check(std::abs(tfdsp::ElectricPianoControls{}.release -
			tfdsp::ElectricPianoDefaultRelease) < 1.0e-12 &&
			defaultDamperRelease > 0.008 && defaultDamperRelease < 0.012 &&
			std::abs(tfdsp::ElectricPianoDamperReleaseSeconds(0.0) - 0.005) <
				1.0e-12 &&
			std::abs(tfdsp::ElectricPianoDamperReleaseSeconds(1.0) - 1.2) <
				1.0e-12,
			"electric piano Release defaults near the published damped-spring relaxation while preserving its creative range");
		tfdsp::ElectricPianoVoice defaultDamperVoice;
		defaultDamperVoice.SetSampleRate(48000.0);
		defaultDamperVoice.Step(0.0, 10.0, 0.8, false, silentControls);
		defaultDamperVoice.Step(0.0, 0.0, 0.8, false, silentControls);
		Check(std::abs(defaultDamperVoice.ModeAmplitudeLifetimeSeconds(1) -
			defaultDamperRelease) < 1.0e-9,
			"electric piano default key-up applies the calibrated damper lifetime to the played mode");
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

		// Dormant Rack voices are not stepped. Controls edited after their previous
		// note becomes silent must therefore be re-primed on the next key-down,
		// rather than spending the new attack smoothing from stale settings. Decay
		// is the most sensitive probe because its deliberate smoothing time is long.
		auto dormantControls = silentControls;
		dormantControls.decay = 0.0;
		dormantControls.release = 0.0;
		auto editedDormantControls = dormantControls;
		editedDormantControls.decay = 1.0;
		tfdsp::ElectricPianoVoice reusedSilentVoice;
		reusedSilentVoice.SetSampleRate(48000.0);
		reusedSilentVoice.Step(0.0, 10.0, 0.8, false, dormantControls);
		for (int sample = 0; sample < 48000 && reusedSilentVoice.IsAudible();
			++sample)
			reusedSilentVoice.Step(0.0, 0.0, 0.0, false, dormantControls);
		const bool reusedVoiceReachedSilence = !reusedSilentVoice.IsAudible();
		reusedSilentVoice.Step(0.0, 10.0, 0.8, false,
			editedDormantControls);
		tfdsp::ElectricPianoVoice freshEditedVoice;
		freshEditedVoice.SetSampleRate(48000.0);
		freshEditedVoice.Step(0.0, 10.0, 0.8, false, editedDormantControls);
		const double reusedEditedLifetime =
			reusedSilentVoice.ModeAmplitudeLifetimeSeconds(0);
		const double freshEditedLifetime =
			freshEditedVoice.ModeAmplitudeLifetimeSeconds(0);
		Check(reusedVoiceReachedSilence &&
			std::abs(reusedEditedLifetime - freshEditedLifetime) < 1.0e-12,
			"electric piano silent voices snap edited controls before the next attack");

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
		tfdsp::ElectricPianoVoice upstreamHeldGateVoice;
		tfdsp::ElectricPianoVoice dampedVoice;
		sustainedVoice.SetSampleRate(48000.0);
		upstreamHeldGateVoice.SetSampleRate(48000.0);
		dampedVoice.SetSampleRate(48000.0);
		for (int sample = 0; sample < 18000; ++sample)
		{
			const double gate = sample < 2400 ? 10.0 : 0.0;
			sustainedVoice.Step(-1.0, gate, 0.8, sample >= 2400,
				silentControls);
			upstreamHeldGateVoice.Step(-1.0, 10.0, 0.8, false,
				silentControls);
			dampedVoice.Step(-1.0, gate, 0.8, false, silentControls);
		}
		Check(sustainedVoice.Energy() > 100.0 * dampedVoice.Energy(),
			"electric piano sustain pedal prevents damper contact after key release");
		Check(std::abs(upstreamHeldGateVoice.Energy() - sustainedVoice.Energy()) <
			1.0e-10 * std::max(1.0, sustainedVoice.Energy()),
			"electric piano held upstream gate preserves pedal-disconnected sustain fallback");
		const double heldPedalEnergy = sustainedVoice.Energy();
		for (int sample = 0; sample < 24000; ++sample)
			sustainedVoice.Step(-1.0, 0.0, 0.0, false, silentControls);
		Check(sustainedVoice.Energy() < 1.0e-3 * heldPedalEnergy,
			"electric piano pedal release applies the damper to sustained notes");
	}

	if (failures == 0)
		std::cout << "All electric-piano tests passed\n";
	return failures == 0 ? 0 : 1;
}
