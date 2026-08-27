#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <memory>

#include "models/PetersonPowerAmplifier.hpp"
#include "models/PetersonPreAmplifier.hpp"
#include "models/PetersonTonePreAmplifier.hpp"
#include "tfdsp/approx.hpp"
#include "tfdsp/sampleRate.hpp"

namespace tfdsp
{

struct ElectricPianoControls
{
	double velocityCurve = 0.5;
	double dynamics = 1.0;
	double body = 0.62;
	double bell = 0.52;
	double coupling = 0.50;
	double hammer = 0.52;
	double tone = 0.55;
	double proximity = 0.48;
	double decay = 0.50;
	double release = 0.24;
	double mechanics = 0.18;
	double drive = 0.32;
	double outputVolume = 0.50;
	double amplifierBass = 0.50;
	double amplifierTreble = 0.50;
	double vibratoSpeed = 0.32;
	double vibratoIntensity = 0.0;
};

inline constexpr double ElectricPianoReferenceFrequency =
	261.6255653005986;

struct ElectricPianoCoupledForkProfile
{
	std::array<double, 2> frequencyRatios{};
	std::array<double, 2> toneBarDisplacementRatios{};
	std::array<double, 2> inverseModalMassRatios{};
	std::array<double, 2> supportReactionLossFactors{};
	double toneBarModalMassRatio{};
};

inline ElectricPianoCoupledForkProfile MakeElectricPianoCoupledForkProfile(
	double coupling, double keyboardPosition)
{
	coupling = std::clamp(std::isfinite(coupling) ? coupling : 0.0,
		0.0, 1.0);
	keyboardPosition = std::clamp(std::isfinite(keyboardPosition) ?
		keyboardPosition : 0.5, 0.0, 1.0);

	// The patent specifies a tone bar whose total mass is at least ten times
	// the tine's. A cantilever's distributed mass does not move uniformly, so
	// its first-mode effective mass is substantially smaller; four tine modal
	// masses is the reduced-coordinate value used here.
	constexpr double ToneBarModalMassRatio = 4.0;
	const double toneBarFrequencyRatio = 0.9970 +
		0.0010 * keyboardPosition;
	// A logarithmic stiffness span moves from an almost isolated tine to a
	// strongly hybridized tuning fork. The gentle quadratic remap places the
	// calibrated physical fork at the control midpoint while retaining both
	// sound-design extremes (normalized log positions 0.0, 0.58, and 1.0).
	const double normalizedLogStiffness = coupling *
		(1.32 - 0.32 * coupling);
	const double couplingStiffnessRatio = 0.00020 *
		std::pow(200.0, normalizedLogStiffness);

	// Solve the symmetric mass-normalized generalized eigenproblem
	// K phi = omega^2 M phi. Stiffnesses are normalized to m_tine*omega_tine^2.
	const double diagonalTine = 1.0 + couplingStiffnessRatio;
	const double diagonalBar = toneBarFrequencyRatio * toneBarFrequencyRatio +
		couplingStiffnessRatio / ToneBarModalMassRatio;
	const double offDiagonal = -couplingStiffnessRatio /
		std::sqrt(ToneBarModalMassRatio);
	const double discriminant = std::sqrt(
		(diagonalTine - diagonalBar) * (diagonalTine - diagonalBar) +
		4.0 * offDiagonal * offDiagonal);
	const std::array<double, 2> eigenvalues{
		0.5 * (diagonalTine + diagonalBar - discriminant),
		0.5 * (diagonalTine + diagonalBar + discriminant)};
	const double frequencyCalibration = 1.0 /
		std::sqrt(std::max(1.0e-12, eigenvalues[1]));

	ElectricPianoCoupledForkProfile profile;
	profile.toneBarModalMassRatio = ToneBarModalMassRatio;
	for (std::size_t mode = 0; mode < 2; ++mode)
	{
		profile.frequencyRatios[mode] = frequencyCalibration *
			std::sqrt(std::max(1.0e-12, eigenvalues[mode]));
		// Normalize each eigenvector to unit tine-tip displacement. This makes
		// hammer and pickup projections explicit; its generalized modal mass is
		// m_tine * (1 + massRatio * barRatio^2).
		const double barRatio = (diagonalTine - eigenvalues[mode]) /
			couplingStiffnessRatio;
		profile.toneBarDisplacementRatios[mode] = barRatio;
		profile.inverseModalMassRatios[mode] = 1.0 /
			(1.0 + ToneBarModalMassRatio * barRatio * barRatio);
		// The mounting block loses energy in proportion to its reaction force.
		// Tine and tone-bar bending forces cancel in the desired asymmetric-fork
		// mode; an isolated tine instead drives the lossy support almost fully.
		// Dividing squared reaction by modal mass makes this invariant to the
		// unit-tine eigenvector normalization above.
		const double supportReaction = 1.0 + ToneBarModalMassRatio *
			toneBarFrequencyRatio * toneBarFrequencyRatio * barRatio;
		profile.supportReactionLossFactors[mode] = supportReaction *
			supportReaction * profile.inverseModalMassRatios[mode];
	}
	return profile;
}

struct ElectricPianoKeyProfile
{
	double fundamentalHz{};
	double modalMassRatio{};
	double pickupSensitivity{};
	double keyboardPosition{};
};

inline ElectricPianoKeyProfile MakeElectricPianoKeyProfile(double pitchVolts)
{
	const double boundedPitch = std::clamp(
		std::isfinite(pitchVolts) ? pitchVolts : 0.0, -6.0, 6.0);
	const double frequencyRatio = std::exp2(boundedPitch);
	// First calibration model: a uniform circular cantilever has length
	// proportional to f^-1/2. With constant diameter its effective modal mass
	// follows the same law. A hammer impulse therefore produces displacement
	// proportional to 1 / (mass * frequency), without any bark-specific curve.
	// The pickup-sensitivity term represents the per-key pickup adjustment used
	// to equalize the complete hammer/tine collision, rather than only the bare
	// cantilever impedance. Real pickups are individually positioned; the
	// shallow frequency tilt and end-range correction are an initial
	// whole-keyboard level calibration.
	const double modalMassRatio = std::pow(frequencyRatio, -0.5);
	const double midiNote = 60.0 + 12.0 * boundedPitch;
	const double keyboardPosition = std::clamp(
		(midiNote - 28.0) / 72.0, 0.0, 1.0);
	return {
		ElectricPianoReferenceFrequency * frequencyRatio,
		modalMassRatio,
		std::pow(frequencyRatio, -0.09) *
			(1.0 + std::pow(2.0 * keyboardPosition - 1.0, 2.0)),
		keyboardPosition};
}

inline double ElectricPianoModeBandlimitGain(double frequency,
	double sampleRate)
{
	if (!std::isfinite(frequency) || !std::isfinite(sampleRate) ||
		!(sampleRate > 0.0))
		return 0.0;
	const double normalizedFrequency = std::abs(frequency) / sampleRate;
	const double taper = std::clamp((normalizedFrequency - 0.32) /
		(0.45 - 0.32), 0.0, 1.0);
	const double smoothTaper = taper * taper * (3.0 - 2.0 * taper);
	return 1.0 - smoothTaper;
}

// Provenance, equations, evidence levels and known calibration gaps for every
// section below are recorded in docs/TfElectricPiano-modelling-notes.md.
class ElectricPianoVoice
{
public:
	ElectricPianoVoice()
		: verticalPositionInterpolator_(CreateX4Resampler_Cheby7()),
		  verticalVelocityInterpolator_(CreateX4Resampler_Cheby7()),
		  horizontalPositionInterpolator_(CreateX4Resampler_Cheby7()),
		  horizontalVelocityInterpolator_(CreateX4Resampler_Cheby7()),
		  pickupDecimator_(CreateX4Resampler_Cheby7())
	{
	}

	void SetNoiseSeed(std::uint32_t seed)
	{
		noiseState_ = seed != 0 ? seed : 0x6d2b79f5u;
		modeCoefficientUpdateCountdown_ = static_cast<int>(noiseState_ %
			static_cast<std::uint32_t>(modeCoefficientUpdatePeriod_));
	}

	void SetSampleRate(double sampleRate)
	{
		sampleRate_ = std::clamp(sampleRate, 8000.0, 768000.0);
		pickupLowPass_ = 0.0;
		hammerNoiseDecay_ = std::exp(-1.0 / (0.0045 * sampleRate_));
		keyReleaseNoiseDecay_ = std::exp(-1.0 / (0.012 * sampleRate_));
		damperNoiseDecay_ = std::exp(-1.0 / (0.018 * sampleRate_));
		mechanicsLowPassCoefficient_ = 1.0 - std::exp(
			-TwoPi * 420.0 / sampleRate_);
		mechanicsHighPassCoefficient_ = 1.0 - std::exp(
			-TwoPi * 3100.0 / sampleRate_);
		controlSmoothingCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.006 * sampleRate_));
		// A decay change alters the radii of every resonant mode. Treating it
		// like a fast tone control creates audible amplitude modulation on held
		// notes, so it gets its own deliberately slow, still-live transition.
		decaySmoothingCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.180 * sampleRate_));
		modeCoefficientUpdatePeriod_ = std::max(1,
			static_cast<int>(std::round(sampleRate_ / 1000.0)));
		modeCoefficientUpdateCountdown_ = static_cast<int>(noiseState_ %
			static_cast<std::uint32_t>(modeCoefficientUpdatePeriod_));
		for (auto& resonator : mechanicalResonators_)
			resonator.state = {};
		coefficientsDirty_ = true;
		cachedTone_ = -1.0;
		verticalPositionInterpolator_->Reset();
		verticalVelocityInterpolator_->Reset();
		horizontalPositionInterpolator_->Reset();
		horizontalVelocityInterpolator_->Reset();
		pickupDecimator_->Reset();
	}

	void Reset()
	{
		for (auto& mode : modes_)
			mode = {};
		lastGate_ = false;
		lastSustain_ = false;
		keyHeld_ = false;
		latchedVelocity_ = 0.0;
		keyPosition_ = 0.5;
		pickupLowPass_ = 0.0;
		hammerNoise_ = 0.0;
		keyReleaseNoise_ = 0.0;
		damperNoise_ = 0.0;
		mechanicsLowPass_ = 0.0;
		mechanicsHighPass_ = 0.0;
		for (auto& resonator : mechanicalResonators_)
			resonator = {};
		hammerPosition_ = 0.0;
		hammerVelocity_ = 0.0;
		hammerIncomingVelocity_ = 1.0;
		hammerInverseMass_ = 1.0;
		contactStiffness_ = 0.0;
		contactExponent_ = 2.1;
		contactLoss_ = 0.0;
		contactAge_ = 0.0;
		contactActive_ = false;
		contactEngaged_ = false;
		contactModeShape_.fill(0.0);
		modeInverseMass_.fill(1.0);
		coupledToneBarRatio_.fill(0.0);
		tineModalMass_ = 1.0;
		coupledForkInitialized_ = false;
		energy_ = 0.0;
		controlsInitialized_ = false;
		coefficientsDirty_ = true;
		timbreDirty_ = true;
		verticalPositionInterpolator_->Reset();
		verticalVelocityInterpolator_->Reset();
		horizontalPositionInterpolator_->Reset();
		horizontalVelocityInterpolator_->Reset();
		pickupDecimator_->Reset();
	}

	double Step(double pitchVolts, double gateVolts, double velocity,
		bool sustain, const ElectricPianoControls& controls)
	{
		if (!std::isfinite(pitchVolts))
			pitchVolts = 0.0;
		if (!std::isfinite(gateVolts))
			gateVolts = 0.0;
		velocity = std::clamp(std::isfinite(velocity) ? velocity : 0.8,
			0.0, 1.0);

		const bool gate = gateVolts >= 1.0;
		if (gate && !lastGate_)
			Strike(pitchVolts, velocity, controls);
		else if (!gate && lastGate_)
		{
			keyHeld_ = false;
			keyReleaseNoise_ = 0.035 + 0.075 * latchedVelocity_;
			TriggerKeyReleaseMechanics();
			if (!sustain)
			{
				damperNoise_ = 0.030 + 0.090 * latchedVelocity_;
				TriggerDamperMechanics();
			}
		}
		if (lastSustain_ && !sustain && !keyHeld_)
		{
			damperNoise_ = std::max(damperNoise_,
				0.030 + 0.090 * latchedVelocity_);
			TriggerDamperMechanics();
		}
		lastGate_ = gate;
		lastSustain_ = sustain;
		if (gate)
			keyHeld_ = true;

		const auto activeControls = SmoothActiveControls(controls);
		const double boundedPitch = std::clamp(pitchVolts, -6.0, 6.0);
		bool pitchChanged = false;
		if (boundedPitch != cachedPitch_ || coefficientsDirty_)
		{
			currentFundamental_ = ElectricPianoReferenceFrequency *
				std::exp2(boundedPitch);
			cachedPitch_ = boundedPitch;
			pitchChanged = true;
		}
		const bool coefficientUpdateTick =
			modeCoefficientUpdateCountdown_ <= 0;
		if (coefficientUpdateTick)
			modeCoefficientUpdateCountdown_ = modeCoefficientUpdatePeriod_ - 1;
		else
			--modeCoefficientUpdateCountdown_;
		const bool damped = !keyHeld_ && !sustain;
		RefreshModeCoefficients(currentFundamental_, damped, activeControls,
			coefficientUpdateTick || pitchChanged);
		RefreshTimbreCoefficients(activeControls);

		if (contactActive_)
			AdvanceCoupledHammerAndModes();
		else
			AdvanceFreeModes();

		double verticalPosition = 0.0;
		double verticalVelocity = 0.0;
		double horizontalPosition = 0.0;
		double horizontalVelocity = 0.0;
		energy_ = 0.0;
		for (std::size_t index = 0; index < modes_.size(); ++index)
		{
			if (!modeActive_[index])
				continue;
			const double verticalWeight = modeOutputWeight_[index] *
				modeBandlimitGain_[index];
			const double horizontalWeight = modeHorizontalWeight_[index] *
				modeBandlimitGain_[index];
			verticalPosition += verticalWeight * modes_[index].real;
			verticalVelocity -= verticalWeight * frequencyVelocityScale_ *
				modeRatio_[index] * modes_[index].imaginary;
			horizontalPosition += horizontalWeight * modes_[index].real;
			horizontalVelocity -= horizontalWeight * frequencyVelocityScale_ *
				modeRatio_[index] * modes_[index].imaginary;
			energy_ += modes_[index].real * modes_[index].real +
				modes_[index].imaginary * modes_[index].imaginary;
		}

		const double white = WhiteNoise();
		mechanicsLowPass_ += mechanicsLowPassCoefficient_ *
			(white - mechanicsLowPass_);
		mechanicsHighPass_ += mechanicsHighPassCoefficient_ *
			(white - mechanicsHighPass_);
		const double middleNoise = mechanicsHighPass_ - mechanicsLowPass_;
		const double hammerColour = (1.0 - hammerNoiseColour_) *
			mechanicsLowPass_ + hammerNoiseColour_ * middleNoise;
		// Key mechanisms are impacts, not continuously running noise sources.
		// These paired, inharmonic resonators restart at a defined phase for
		// every event, keeping chords tight. A little filtered noise prevents
		// the transients from sounding like tiny tuned oscillators.
		const double hammerCarrier = hammerNoise_ > 1.0e-8 ?
			0.26 * AdvanceMechanicalResonator(mechanicalResonators_[0]) +
			0.13 * AdvanceMechanicalResonator(mechanicalResonators_[1]) +
			0.06 * hammerColour : 0.0;
		const double releaseCarrier = keyReleaseNoise_ > 1.0e-8 ?
			0.22 * AdvanceMechanicalResonator(mechanicalResonators_[2]) +
			0.11 * AdvanceMechanicalResonator(mechanicalResonators_[3]) +
			0.05 * mechanicsLowPass_ : 0.0;
		const double damperCarrier = damperNoise_ > 1.0e-8 ?
			0.20 * AdvanceMechanicalResonator(mechanicalResonators_[4]) +
			0.13 * AdvanceMechanicalResonator(mechanicalResonators_[5]) +
			0.07 * middleNoise : 0.0;
		const double mechanicalSignal = mechanicsLevel_ *
			(0.74 * hammerNoise_ * hammerCarrier +
				0.20 * keyReleaseNoise_ * releaseCarrier +
				0.28 * damperNoise_ * damperCarrier);

		const auto verticalPositions = verticalPositionInterpolator_->Upsample(
			TineDisplacementScale * verticalPosition);
		const auto verticalVelocities = verticalVelocityInterpolator_->Upsample(
			verticalVelocity);
		const auto horizontalPositions = horizontalPositionInterpolator_->Upsample(
			TineDisplacementScale * horizontalPosition);
		const auto horizontalVelocities = horizontalVelocityInterpolator_->Upsample(
			horizontalVelocity);
		Eigen::Array<double, PickupOversamplingFactor, 1> pickupValues;
		for (int index = 0; index < PickupOversamplingFactor; ++index)
		{
			const auto gradient = MagneticFluxGradient(
				pickupVerticalOffset_ + verticalPositions(index),
				pickupHorizontalOffset_ + horizontalPositions(index), pickupGap_);
			const double emf = gradient[0] * verticalVelocities(index) +
				gradient[1] * horizontalVelocities(index);
			const double pickup = pickupVoltageScale_ * emf;
			pickupValues(index) = std::clamp(pickup, -12.0, 12.0) +
				mechanicalSignal;
		}
		const double pickup = pickupDecimator_->Downsample(pickupValues);
		pickupLowPass_ += pickupLowPassCoefficient_ *
			(pickup - pickupLowPass_);
		hammerNoise_ *= hammerNoiseDecay_;
		keyReleaseNoise_ *= keyReleaseNoiseDecay_;
		damperNoise_ *= damperNoiseDecay_;

		if (energy_ < 1.0e-14 && hammerNoise_ < 1.0e-8 &&
			keyReleaseNoise_ < 1.0e-8 && damperNoise_ < 1.0e-8)
		{
			for (auto& mode : modes_)
				mode = {};
			pickupLowPass_ *= 0.9;
		}

		const double result = pickupLowPass_;
		return std::isfinite(result) ? result : 0.0;
	}

	double Energy() const { return energy_; }
	double Activity() const
	{
		return energy_ + hammerNoise_ * hammerNoise_ +
			keyReleaseNoise_ * keyReleaseNoise_ +
			damperNoise_ * damperNoise_ +
			pickupLowPass_ * pickupLowPass_;
	}
	bool IsAudible() const { return Activity() > 1.0e-12; }
	bool GateHigh() const { return lastGate_; }
	double NotePitch() const { return notePitch_; }
	bool ContactActive() const { return contactActive_; }
	double ContactAge() const { return contactAge_; }

private:
	struct Mode
	{
		double real{};
		double imaginary{};
	};

	struct MechanicalResonator
	{
		Mode state{};
		double realCoefficient = 1.0;
		double imaginaryCoefficient{};
	};

	static double Clamp01(double value)
	{
		return std::clamp(std::isfinite(value) ? value : 0.0, 0.0, 1.0);
	}

	void TriggerMechanicalResonator(MechanicalResonator& resonator,
		double frequency, double amplitude = 1.0)
	{
		const double boundedFrequency = std::clamp(frequency, 20.0,
			0.43 * sampleRate_);
		const double angle = TwoPi * boundedFrequency / sampleRate_;
		resonator.realCoefficient = std::cos(angle);
		resonator.imaginaryCoefficient = std::sin(angle);
		// A cosine-phase start puts the physical impact at the event sample;
		// it cannot wait for an arbitrary random carrier crossing.
		resonator.state.real = amplitude;
		resonator.state.imaginary = 0.0;
	}

	static double AdvanceMechanicalResonator(
		MechanicalResonator& resonator)
	{
		const double output = resonator.state.real;
		const double oldReal = resonator.state.real;
		const double oldImaginary = resonator.state.imaginary;
		resonator.state.real = resonator.realCoefficient * oldReal -
			resonator.imaginaryCoefficient * oldImaginary;
		resonator.state.imaginary = resonator.imaginaryCoefficient * oldReal +
			resonator.realCoefficient * oldImaginary;
		return output;
	}

	void TriggerKeyReleaseMechanics()
	{
		const double baseFrequency = 310.0 + 450.0 * keyPosition_;
		TriggerMechanicalResonator(mechanicalResonators_[2], baseFrequency);
		TriggerMechanicalResonator(mechanicalResonators_[3],
			1.86 * baseFrequency);
	}

	void TriggerDamperMechanics()
	{
		const double baseFrequency = 950.0 + 850.0 * keyPosition_;
		TriggerMechanicalResonator(mechanicalResonators_[4], baseFrequency);
		TriggerMechanicalResonator(mechanicalResonators_[5],
			1.52 * baseFrequency);
	}

	ElectricPianoControls SmoothActiveControls(
		const ElectricPianoControls& controls)
	{
		if (!controlsInitialized_)
		{
			smoothedControls_ = controls;
			controlsInitialized_ = true;
		}
		auto smooth = [&](double& current, double target, double coefficient)
		{
			target = Clamp01(target);
			current += coefficient * (target - current);
		};
		smooth(smoothedControls_.body, controls.body,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.bell, controls.bell,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.coupling, controls.coupling,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.tone, controls.tone,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.proximity, controls.proximity,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.decay, controls.decay,
			decaySmoothingCoefficient_);
		smooth(smoothedControls_.release, controls.release,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.mechanics, controls.mechanics,
			controlSmoothingCoefficient_);
		return smoothedControls_;
	}

	void Strike(double pitchVolts, double velocity,
		const ElectricPianoControls& controls)
	{
		const auto key = MakeElectricPianoKeyProfile(pitchVolts);
		notePitch_ = std::clamp(pitchVolts, -6.0, 6.0);
		keyPosition_ = key.keyboardPosition;
		keyPickupSensitivity_ = key.pickupSensitivity;
		coefficientsDirty_ = true;
		timbreDirty_ = true;
		const double velocityCurve = Clamp01(controls.velocityCurve);
		const double dynamics = Clamp01(controls.dynamics);
		const double hammer = Clamp01(controls.hammer);
		const double gamma = std::exp2(1.0 - 2.0 * velocityCurve);
		const double curvedVelocity = velocity > 0.0 ?
			std::pow(velocity, gamma) : 0.0;
		latchedVelocity_ = curvedVelocity;
		const double compressedVelocity = curvedVelocity /
			(0.28 + 0.72 * curvedVelocity);
		const double dynamicAmplitude = (1.0 - dynamics) *
			compressedVelocity + dynamics * curvedVelocity;

		// The hammer is an independent finite mass. MIDI velocity sets its
		// incoming speed; no post-contact amplitude or brightness envelope is
		// prescribed. The key-dependent term represents the graduated neoprene
		// tips used across the real keyboard, while HAMMER moves around that
		// baseline by changing contact material properties only.
		const double contactHardness = std::clamp(0.12 + 0.52 * hammer +
			0.36 * keyPosition_, 0.0, 1.0);
		const double hammerMass = 0.42 * (1.08 - 0.18 * keyPosition_);
		hammerInverseMass_ = 1.0 / hammerMass;
		hammerPosition_ = 0.0;
		const double requestedHammerVelocity = 1150.0 * dynamicAmplitude *
			(0.86 + 0.14 * std::sqrt(std::max(0.0, dynamicAmplitude)));
		// Below this speed the calibrated Hunt-Crossley loss law already clamps
		// its incoming-velocity denominator. Such strikes are many orders below
		// the audible output range but can otherwise spend seconds in an almost
		// static 16x collision when driven by tiny positive CV residue. Classify
		// them as zero-energy strikes before contact begins; never truncate an
		// active collision based on elapsed time.
		const bool audibleHammerStrike =
			requestedHammerVelocity >= MinimumHammerContactVelocity;
		hammerVelocity_ = audibleHammerStrike ? requestedHammerVelocity : 0.0;
		hammerIncomingVelocity_ = std::max(MinimumHammerContactVelocity,
			requestedHammerVelocity);
		// The original instrument progresses from soft neoprene in the bass to
		// wrapped, extra-hard tips in the treble. That change is substantially
		// larger than the front-panel adjustment around any one key.
		// Keep the default physical calibration fixed while giving the panel
		// control a useful six-octave stiffness span. The upper limit prevents
		// the already wrapped treble tips from becoming numerically impulsive.
		contactStiffness_ = std::min(1.6e9, 4.0e5 * std::exp2(
			6.0 * hammer + 9.0 * keyPosition_ - 2.08));
		contactExponent_ = 1.95 + 0.32 * contactHardness;
		const double restitution = 0.24 + 0.30 * contactHardness;
		contactLoss_ = 1.5 * (1.0 - restitution * restitution) /
			hammerIncomingVelocity_;
		contactAge_ = 0.0;
		contactActive_ = audibleHammerStrike;
		contactEngaged_ = false;
		contactModeShape_[0] = 0.82;
		contactModeShape_[1] = 0.82;
		contactModeShape_[2] = -0.055;
		tineModalMass_ = key.modalMassRatio;
		modeInverseMass_[2] = 1.0 / (key.modalMassRatio * 1.65);
		for (std::size_t index = 3; index < modes_.size(); ++index)
		{
			contactModeShape_[index] = CantileverModeShapeAtStrike(index,
				keyPosition_);
			modeInverseMass_[index] = 1.0 /
				(key.modalMassRatio *
					AttackModeModalMassMultipliers[index - 3]);
		}

		// Mechanics remains a separate sound-design layer. Tip hardness has only
		// a restrained influence on its colour so it cannot duplicate the modal
		// attack as a second, disconnected click.
		hammerNoiseColour_ = 0.24 + 0.22 * contactHardness;
		const double hammerMechanicalFrequency = 680.0 +
			1250.0 * keyPosition_ + 620.0 * contactHardness;
		TriggerMechanicalResonator(mechanicalResonators_[0],
			hammerMechanicalFrequency);
		TriggerMechanicalResonator(mechanicalResonators_[1],
			2.17 * hammerMechanicalFrequency);
		hammerNoise_ = dynamicAmplitude * (0.085 + 0.075 * curvedVelocity) *
			(0.82 + 0.18 * keyPosition_);
		keyReleaseNoise_ = 0.0;
		damperNoise_ = 0.0;
		keyHeld_ = true;
		energy_ = 1.0;
	}

	void RefreshTimbreCoefficients(const ElectricPianoControls& controls)
	{
		const double body = Clamp01(controls.body);
		const double bell = Clamp01(controls.bell);
		const double proximity = Clamp01(controls.proximity);
		const double tone = Clamp01(controls.tone);
		const double mechanics = Clamp01(controls.mechanics);
		if (!timbreDirty_ && body == cachedBodyWeight_ &&
			bell == cachedBellWeight_ && proximity == cachedProximity_ &&
			tone == cachedTone_ &&
			mechanics == cachedMechanics_)
			return;

		// Both coupled normal coordinates are normalized to the same unit tine
		// displacement. The pickup must therefore observe them with identical
		// gain; changing their weights independently would break the underlying
		// two-body system and turn Coupling back into a hidden level control.
		modeOutputWeight_[0] = 0.58 + 0.52 * body;
		modeOutputWeight_[1] = modeOutputWeight_[0];
		modeOutputWeight_[2] = 0.06 + 0.12 * body;
		for (std::size_t index = 3; index < modeOutputWeight_.size(); ++index)
			modeOutputWeight_[index] = 0.10 + 0.38 * bell;

		// Imperfect transverse coupling gives the tine tip a shallow elliptical
		// orbit. The near-fundamental third mode carries most of the horizontal
		// motion; attack modes alternate orientation as observed on a real fork.
		modeHorizontalWeight_[0] = 0.010;
		modeHorizontalWeight_[1] = -0.012;
		modeHorizontalWeight_[2] = 0.42 + 0.20 * body;
		for (std::size_t index = 3; index < modeHorizontalWeight_.size(); ++index)
			modeHorizontalWeight_[index] = (index % 2 == 0 ? -1.0 : 1.0) *
			(0.045 + 0.10 * bell);

		// The service adjustment called "timbre" is the tine's vertical alignment
		// to the pickup pole. Moving toward the pole centre suppresses the linear
		// fundamental relative to curvature-generated harmonics. PROXIMITY controls
		// the independent front-to-back distance and therefore magnetic curvature.
		pickupGap_ = 1.32 - 0.78 * proximity;
		// At minimum Tone the tine rests close to the primary edge's maximum
		// gradient. Raising Tone moves it onto one flank, trading odd-harmonic
		// symmetry for the stronger even/sideband mixture used when voicing a
		// Rhodes for bite.
		pickupVerticalOffset_ = 0.34 + 0.22 * tone * tone +
			0.020 * keyPosition_;
		pickupHorizontalOffset_ = 0.10 + 0.035 * keyPosition_;
		const auto alignmentGradient = MagneticFluxGradient(
			pickupVerticalOffset_, pickupHorizontalOffset_, ReferencePickupGap);
		const double alignmentMagnitude = std::sqrt(
			alignmentGradient[0] * alignmentGradient[0] +
			alignmentGradient[1] * alignmentGradient[1]);
		const auto proximityGradient = MagneticFluxGradient(
			pickupVerticalOffset_, pickupHorizontalOffset_, pickupGap_);
		const double proximityMagnitude = std::sqrt(
			proximityGradient[0] * proximityGradient[0] +
			proximityGradient[1] * proximityGradient[1]);
		// A real service adjustment changes both sensitivity and curvature, but the
		// uncorrected pole field spans well over 16 dB and turns Proximity into a
		// second amplifier Drive control. Retaining only a small part of that raw dB
		// span keeps the adjustment audible while leaving the excursion-dependent
		// magnetic curvature intact.
		const double proximityGain = std::pow(alignmentMagnitude /
			std::max(1.0e-6, proximityMagnitude), 0.85);
		pickupAlignmentGain_ = std::clamp(ReferenceGradientMagnitude() /
			std::max(1.0e-6, alignmentMagnitude) * proximityGain, 0.40, 5.0);
		pickupVoltageScale_ = 0.135 * keyPickupSensitivity_ *
			pickupAlignmentGain_ * inverseReferenceGradient_;
		const double cutoff = std::min(0.42 * sampleRate_, 16500.0);
		pickupLowPassCoefficient_ = 1.0 - std::exp(
			-TwoPi * cutoff / sampleRate_);
		// Preserve a wide sound-design range while making the calibrated default
		// a quiet mechanical contribution rather than a parallel percussion layer.
		mechanicsLevel_ = 0.060 * mechanics * (0.30 + 0.70 * mechanics);
		cachedBodyWeight_ = body;
		cachedBellWeight_ = bell;
		cachedProximity_ = proximity;
		cachedTone_ = tone;
		cachedMechanics_ = mechanics;
		timbreDirty_ = false;
	}

	void RefreshModeCoefficients(double fundamental, bool damped,
		const ElectricPianoControls& controls, bool controlTick)
	{
		const double decay = Clamp01(controls.decay);
		const double bell = Clamp01(controls.bell);
		const double coupling = Clamp01(controls.coupling);
		const double release = Clamp01(controls.release);
		const bool stateChange = coefficientsDirty_ ||
			fundamental != cachedFundamental_ || damped != cachedDamped_;
		if (!stateChange && !controlTick)
			return;
		const bool couplingChanged = coefficientsDirty_ ||
			coupling != cachedCoupling_;
		const bool frequencyChanged = coefficientsDirty_ ||
			fundamental != cachedFundamental_ || couplingChanged;
		const bool decayChanged = coefficientsDirty_ ||
			damped != cachedDamped_ || decay != cachedDecay_ ||
			bell != cachedBell_ || couplingChanged || release != cachedRelease_;
		if (!frequencyChanged && !decayChanged)
			return;

		const bool transformCoupledState = coupledForkInitialized_ &&
			frequencyChanged &&
			modeAngularFrequency_[0] > 0.0 && modeAngularFrequency_[1] > 0.0;
		double oldTinePosition = 0.0;
		double oldToneBarPosition = 0.0;
		double oldTineVelocity = 0.0;
		double oldToneBarVelocity = 0.0;
		if (transformCoupledState)
		{
			for (std::size_t mode = 0; mode < 2; ++mode)
			{
				const double modalVelocity = -modeAngularFrequency_[mode] *
					modes_[mode].imaginary;
				oldTinePosition += modes_[mode].real;
				oldToneBarPosition += coupledToneBarRatio_[mode] *
					modes_[mode].real;
				oldTineVelocity += modalVelocity;
				oldToneBarVelocity += coupledToneBarRatio_[mode] * modalVelocity;
			}
		}

		const double keyboardDecayScale = 1.18 - 0.50 * keyPosition_;
		if (couplingChanged)
		{
			// Coupled-fork eigenvectors and support reaction depend on Coupling and
			// key geometry, but not instantaneous pitch. Cache them so future FM only
			// updates oscillator angles and the state-preserving velocity transform.
			const auto coupledFork = MakeElectricPianoCoupledForkProfile(coupling,
				keyPosition_);
			for (std::size_t index = 0; index < 2; ++index)
			{
				modeRatio_[index] = coupledFork.frequencyRatios[index];
				coupledToneBarRatio_[index] =
					coupledFork.toneBarDisplacementRatios[index];
				coupledSupportLossFactor_[index] =
					coupledFork.supportReactionLossFactors[index];
				modeInverseMass_[index] =
					coupledFork.inverseModalMassRatios[index] /
					std::max(1.0e-9, tineModalMass_);
				contactModeShape_[index] = 0.82;
			}
			modeRatio_[2] = 1.0018 + 0.0010 * (keyPosition_ - 0.5);
			for (std::size_t index = 3; index < modes_.size(); ++index)
				modeRatio_[index] = KeyboardModeRatio(index, keyPosition_);
		}

		if (decayChanged)
		{
			const double releaseSeconds = 0.012 * std::pow(100.0, release);
			const double decayScale = std::pow(4.0, 2.0 * (decay - 0.5));
			const double bellDecay = 0.82 + 0.38 * bell;
			// Material loss is small. Most fundamental loss occurs when an
			// unbalanced fork drives its resilient mounting block and harp. The
			// support coordinate has been analytically reduced into the reaction
			// factors above; anti-phase tine/bar motion cancels that force and raises
			// Q, while an isolated tine rings briefly. The keyboard term represents
			// increasing mounting loss toward the short, stiff treble assemblies.
			const double intrinsicForkLifetime = IntrinsicForkDecaySeconds *
				keyboardDecayScale;
			const double supportLossRate = 0.45 + 3.55 * keyPosition_;
			for (std::size_t index = 0; index < modes_.size(); ++index)
			{
				double lifetime;
				if (index < 2)
				{
					const double lossRate = 1.0 / intrinsicForkLifetime +
						supportLossRate * coupledSupportLossFactor_[index];
					lifetime = decayScale / std::max(1.0e-9, lossRate);
				}
				else
					lifetime = HigherModeBaseDecaySeconds[index - 2] *
						decayScale * keyboardDecayScale;
				if (index >= 3)
					lifetime *= bellDecay;
				if (damped)
					lifetime = std::min(lifetime, releaseSeconds *
						(index >= 3 ? 0.55 : 1.0));
				modeRadius_[index] = std::exp(-1.0 /
					(std::max(0.002, lifetime) * sampleRate_));
				modeSubstepRadius_[index] = std::exp(-1.0 /
					(std::max(0.002, lifetime) * sampleRate_ *
						static_cast<double>(ContactOversamplingFactor)));
			}
		}

		frequencyVelocityScale_ = fundamental /
			ElectricPianoReferenceFrequency;
		for (std::size_t index = 0; index < modes_.size(); ++index)
		{
			const double frequency = fundamental * modeRatio_[index];
			const double normalizedFrequency = frequency / sampleRate_;
			modeBandlimitGain_[index] = ElectricPianoModeBandlimitGain(
				frequency, sampleRate_);
			modeActive_[index] = normalizedFrequency < 0.49;
			if (!modeActive_[index])
			{
				modes_[index] = {};
				modeRealCoefficient_[index] = 0.0;
				modeImaginaryCoefficient_[index] = 0.0;
				modeSubstepRealCoefficient_[index] = 0.0;
				modeSubstepImaginaryCoefficient_[index] = 0.0;
				modeAngularFrequency_[index] = 0.0;
				continue;
			}
			const double angle = TwoPi * frequency / sampleRate_;
			modeRealCoefficient_[index] = modeRadius_[index] * std::cos(angle);
			modeImaginaryCoefficient_[index] = modeRadius_[index] * std::sin(angle);
			const double substepAngle = angle /
				static_cast<double>(ContactOversamplingFactor);
			modeSubstepRealCoefficient_[index] = modeSubstepRadius_[index] *
				std::cos(substepAngle);
			modeSubstepImaginaryCoefficient_[index] = modeSubstepRadius_[index] *
				std::sin(substepAngle);
			modeAngularFrequency_[index] = TwoPi * frequency;
		}

		if (transformCoupledState && modeActive_[0] && modeActive_[1])
		{
			const double ratioDifference = coupledToneBarRatio_[0] -
				coupledToneBarRatio_[1];
			if (std::abs(ratioDifference) > 1.0e-9)
			{
				modes_[0].real = (oldToneBarPosition -
					coupledToneBarRatio_[1] * oldTinePosition) /
					ratioDifference;
				modes_[1].real = oldTinePosition - modes_[0].real;
				const double modeZeroVelocity = (oldToneBarVelocity -
					coupledToneBarRatio_[1] * oldTineVelocity) /
					ratioDifference;
				const double modeOneVelocity = oldTineVelocity - modeZeroVelocity;
				modes_[0].imaginary = -modeZeroVelocity /
					modeAngularFrequency_[0];
				modes_[1].imaginary = -modeOneVelocity /
					modeAngularFrequency_[1];
			}
		}
		coupledForkInitialized_ = true;
		cachedFundamental_ = fundamental;
		cachedDecay_ = decay;
		cachedBell_ = bell;
		cachedCoupling_ = coupling;
		cachedRelease_ = release;
		cachedDamped_ = damped;
		coefficientsDirty_ = false;
	}

	void AdvanceFreeModes()
	{
		for (std::size_t index = 0; index < modes_.size(); ++index)
		{
			if (!modeActive_[index])
				continue;
			const double oldReal = modes_[index].real;
			const double oldImaginary = modes_[index].imaginary;
			modes_[index].real = modeRealCoefficient_[index] * oldReal -
				modeImaginaryCoefficient_[index] * oldImaginary;
			modes_[index].imaginary = modeImaginaryCoefficient_[index] *
				oldReal + modeRealCoefficient_[index] * oldImaginary;
		}
	}

	void AdvanceCoupledHammerAndModes()
	{
		// A free hammer has no preload that could sustain contact, so a valid
		// collision must end through physical separation. Do not impose an elapsed-
		// time cutoff: soft bass strikes legitimately remain compressed longest.
		// Only non-finite state is treated as a numerical failure below.
		const double timeStep = 1.0 /
			(sampleRate_ * static_cast<double>(ContactOversamplingFactor));
		auto contactPointState = [&]()
		{
			double position = 0.0;
			double velocity = 0.0;
			for (std::size_t index = 0; index < modes_.size(); ++index)
			{
				if (!modeActive_[index])
					continue;
				position -= contactModeShape_[index] * modes_[index].real;
				velocity += contactModeShape_[index] *
					modeAngularFrequency_[index] * modes_[index].imaginary;
			}
			return std::array<double, 2>{position, velocity};
		};
		for (int substep = 0; substep < ContactOversamplingFactor; ++substep)
		{
			if (contactActive_)
			{
				const auto tine = contactPointState();
				const double compression = hammerPosition_ - tine[0];
				const double relativeVelocity = hammerVelocity_ - tine[1];
				if (!std::isfinite(compression) ||
					!std::isfinite(relativeVelocity))
				{
					for (auto& mode : modes_)
						mode = {};
					hammerPosition_ = 0.0;
					hammerVelocity_ = 0.0;
					contactActive_ = false;
					contactEngaged_ = false;
				}
				else
				{
					double force = 0.0;
					if (compression > 0.0)
					{
						contactEngaged_ = true;
						const double hysteresis = std::clamp(
							1.0 + contactLoss_ * relativeVelocity, 0.04, 2.2);
						force = contactStiffness_ *
							std::pow(compression, contactExponent_) * hysteresis;
					}
					if (!std::isfinite(force))
					{
						for (auto& mode : modes_)
							mode = {};
						hammerPosition_ = 0.0;
						hammerVelocity_ = 0.0;
						contactActive_ = false;
						contactEngaged_ = false;
					}
					else
					{
						const double impulse = force * timeStep;
						hammerVelocity_ -= impulse * hammerInverseMass_;
						for (std::size_t index = 0; index < modes_.size(); ++index)
						{
							if (!modeActive_[index] ||
								!(modeAngularFrequency_[index] > 0.0))
								continue;
							modes_[index].imaginary += contactModeShape_[index] *
								modeInverseMass_[index] * impulse /
								modeAngularFrequency_[index];
						}
						hammerPosition_ += timeStep * hammerVelocity_;
					}
				}
			}

			for (std::size_t index = 0; index < modes_.size(); ++index)
			{
				if (!modeActive_[index])
					continue;
				const double oldReal = modes_[index].real;
				const double oldImaginary = modes_[index].imaginary;
				modes_[index].real = modeSubstepRealCoefficient_[index] *
					oldReal - modeSubstepImaginaryCoefficient_[index] *
					oldImaginary;
				modes_[index].imaginary = modeSubstepImaginaryCoefficient_[index] *
					oldReal + modeSubstepRealCoefficient_[index] * oldImaginary;
			}

			if (contactActive_)
			{
				contactAge_ += timeStep;
				const auto tine = contactPointState();
				const double compression = hammerPosition_ - tine[0];
				const double relativeVelocity = hammerVelocity_ - tine[1];
				if (!std::isfinite(hammerPosition_) ||
					!std::isfinite(hammerVelocity_) ||
					!std::isfinite(compression) ||
					!std::isfinite(relativeVelocity))
				{
					for (auto& mode : modes_)
						mode = {};
					hammerPosition_ = 0.0;
					hammerVelocity_ = 0.0;
					contactActive_ = false;
					contactEngaged_ = false;
				}
				else if (contactEngaged_ && compression <= 0.0 &&
					relativeVelocity < 0.0)
				{
					contactActive_ = false;
					contactEngaged_ = false;
				}
			}
		}
	}

	static double KeyboardModeRatio(std::size_t index,
		double keyboardPosition)
	{
		if (index < 3)
			return 1.0;
		const double position = std::clamp(keyboardPosition, 0.0, 1.0);
		const double scale = 1.0 + (1.0 - 2.0 * position) *
			ModeRatioKeyboardSpread[index];
		return ModeRatios[index] * scale;
	}

	static double CantileverModeShapeAtStrike(std::size_t index,
		double keyboardPosition)
	{
		if (index < 3)
			return 1.0;
		const double beta = CantileverWaveNumbers[index - 2];
		const double strikePosition = 0.31 + 0.07 *
			std::clamp(keyboardPosition, 0.0, 1.0);
		const double sigma = (std::cosh(beta) + std::cos(beta)) /
			(std::sinh(beta) + std::sin(beta));
		auto shape = [&](double position)
		{
			const double argument = beta * position;
			return std::cosh(argument) - std::cos(argument) - sigma *
				(std::sinh(argument) - std::sin(argument));
		};
		return std::clamp(shape(strikePosition) / shape(1.0), -1.5, 1.5);
	}

	static std::array<double, 2> MagneticFluxGradient(double vertical,
		double horizontal, double gap)
	{
		// The Peterson-era pole is not a rounded, symmetric point source. Its two
		// offset protuberances concentrate flux into a steep transition whose
		// gradient remains unidirectional across the normal tine excursion. Model
		// the linked flux as two softened magnetic edges: the narrow primary edge
		// supplies the harmonic curvature, while the broader return edge gives the
		// published field models their strong spatial asymmetry. Differentiating
		// this scalar flux map preserves Faraday's-law consistency in both planes.
		constexpr double HorizontalFieldScale = 0.62;
		constexpr double GapRegularization = 0.020;
		const double radialDistance = std::sqrt(gap * gap +
			HorizontalFieldScale * horizontal * horizontal + GapRegularization);
		const double radialDerivative = HorizontalFieldScale * horizontal /
			std::max(1.0e-9, radialDistance);
		// Both pole edges share the same radial falloff. This fractional power is
		// one of the dominant steady-state costs because the pickup is evaluated
		// four times per sample for every voice; calculate it once per field point.
		const double radialFalloff = PowNegative1p3(radialDistance);

		double verticalGradient = 0.0;
		double radialGradient = 0.0;
		auto accumulateEdge = [&](double edgePosition, double weight,
			double edgeRadius)
		{
			constexpr double GapBroadening = 0.180;
			const double width = edgeRadius + GapBroadening * radialDistance;
			const double inverseWidth = 1.0 / std::max(1.0e-9, width);
			const double displacement = vertical - edgePosition;
			const double argument = displacement * inverseWidth;
			const double transition = TanhPade76(argument);
			const double transitionDerivative = 1.0 -
				transition * transition;
			const double amplitude = weight * radialFalloff;
			verticalGradient += amplitude * transitionDerivative * inverseWidth;
			const double amplitudeDerivative = -FieldFalloff * amplitude /
				std::max(1.0e-9, radialDistance);
			const double argumentDerivative = -displacement * GapBroadening *
				inverseWidth * inverseWidth;
			radialGradient += amplitudeDerivative * transition + amplitude *
				transitionDerivative * argumentDerivative;
		};

		accumulateEdge(0.34, 1.0, 0.030);
		accumulateEdge(-0.46, 0.20, 0.095);
		return {verticalGradient, radialGradient * radialDerivative};
	}

	static double ReferenceGradientMagnitude()
	{
		const auto gradient = MagneticFluxGradient(ReferencePickupVertical,
			ReferencePickupHorizontal, ReferencePickupGap);
		return std::sqrt(gradient[0] * gradient[0] +
			gradient[1] * gradient[1]);
	}

	double WhiteNoise()
	{
		noiseState_ ^= noiseState_ << 13;
		noiseState_ ^= noiseState_ >> 17;
		noiseState_ ^= noiseState_ << 5;
		return 2.0 * (static_cast<double>(noiseState_) /
			static_cast<double>(UINT32_MAX)) - 1.0;
	}

	static constexpr double TwoPi = 6.283185307179586476925286766559;
	static constexpr double FieldFalloff = 1.30;
	static constexpr double TineDisplacementScale = 0.56;
	static constexpr int ContactOversamplingFactor = 16;
	static constexpr double MinimumHammerContactVelocity = 1.0;
	static constexpr int PickupOversamplingFactor =
		X4Resampler_Order7::ResamplingFactor;
	static constexpr double ReferencePickupGap = 1.32 - 0.78 * 0.48;
	static constexpr double ReferencePickupVertical = 0.34 + 0.22 * 0.55 * 0.55 +
		0.020 * ((60.0 - 28.0) / 72.0);
	static constexpr double ReferencePickupHorizontal = 0.10 + 0.035 *
		((60.0 - 28.0) / 72.0);
	static constexpr std::array<double, 8> ModeRatios{
		1.0, 1.0, 1.0, 6.267, 17.55, 34.39, 56.84, 83.0};
	static constexpr std::array<double, 8> ModeRatioKeyboardSpread{
		0.0, 0.0, 0.0, 0.035, 0.050, 0.065, 0.080, 0.095};
	static constexpr std::array<double, 6> CantileverWaveNumbers{
		1.875104, 4.694091, 7.854757, 10.995541, 14.137168, 17.278760};
	static constexpr double IntrinsicForkDecaySeconds = 6.8;
	static constexpr std::array<double, 6> HigherModeBaseDecaySeconds{
		1.9, 0.20, 0.120, 0.070, 0.040, 0.024};
	static constexpr std::array<double, 5> AttackModeModalMassMultipliers{
		0.92, 0.96, 1.02, 1.10, 1.22};

	std::array<Mode, 8> modes_{};
	std::array<double, 8> modeRealCoefficient_{};
	std::array<double, 8> modeImaginaryCoefficient_{};
	std::array<double, 8> modeRadius_{};
	std::array<double, 8> modeSubstepRealCoefficient_{};
	std::array<double, 8> modeSubstepImaginaryCoefficient_{};
	std::array<double, 8> modeSubstepRadius_{};
	std::array<double, 8> modeAngularFrequency_{};
	std::array<double, 8> modeOutputWeight_{};
	std::array<double, 8> modeHorizontalWeight_{};
	std::array<double, 8> modeRatio_{};
	std::array<double, 8> modeBandlimitGain_{};
	std::array<bool, 8> modeActive_{};
	std::array<double, 8> contactModeShape_{};
	std::array<double, 8> modeInverseMass_{};
	std::array<double, 2> coupledToneBarRatio_{};
	std::array<double, 2> coupledSupportLossFactor_{};
	double sampleRate_ = 48000.0;
	double latchedVelocity_{};
	double keyPosition_ = 0.5;
	double tineModalMass_ = 1.0;
	double notePitch_{};
	double frequencyVelocityScale_ = 1.0;
	double keyPickupSensitivity_ = 1.0;
	double pickupLowPass_{};
	double pickupLowPassCoefficient_{};
	double pickupGap_ = 1.0;
	double pickupVerticalOffset_ = ReferencePickupVertical;
	double pickupHorizontalOffset_ = ReferencePickupHorizontal;
	double inverseReferenceGradient_ = 1.0 / ReferenceGradientMagnitude();
	double pickupAlignmentGain_ = 1.0;
	double pickupVoltageScale_{};
	double mechanicsLevel_{};
	double hammerNoiseColour_ = 0.5;
	double hammerNoiseDecay_ = 0.9979;
	double keyReleaseNoiseDecay_ = 0.9986;
	double damperNoiseDecay_ = 0.9992;
	double hammerNoise_{};
	double keyReleaseNoise_{};
	double damperNoise_{};
	double mechanicsLowPass_{};
	double mechanicsHighPass_{};
	double mechanicsLowPassCoefficient_ = 0.054;
	double mechanicsHighPassCoefficient_ = 0.334;
	std::array<MechanicalResonator, 6> mechanicalResonators_{};
	double hammerPosition_{};
	double hammerVelocity_{};
	double hammerIncomingVelocity_ = 1.0;
	double hammerInverseMass_ = 1.0;
	double contactStiffness_{};
	double contactExponent_ = 2.1;
	double contactLoss_{};
	double contactAge_{};
	double controlSmoothingCoefficient_ = 0.0035;
	double decaySmoothingCoefficient_ = 0.00012;
	int modeCoefficientUpdatePeriod_ = 48;
	int modeCoefficientUpdateCountdown_{};
	ElectricPianoControls smoothedControls_{};
	double energy_{};
	double cachedFundamental_{};
	double currentFundamental_{};
	double cachedPitch_ = -1000.0;
	double cachedDecay_{};
	double cachedBell_{};
	double cachedCoupling_ = -1.0;
	double cachedRelease_{};
	double cachedBodyWeight_ = -1.0;
	double cachedBellWeight_ = -1.0;
	double cachedProximity_ = -1.0;
	double cachedTone_ = -1.0;
	double cachedMechanics_ = -1.0;
	std::uint32_t noiseState_ = 0x6d2b79f5u;
	bool cachedDamped_{};
	bool coefficientsDirty_ = true;
	bool timbreDirty_ = true;
	bool controlsInitialized_{};
	bool coupledForkInitialized_{};
	bool contactActive_{};
	bool contactEngaged_{};
	bool lastGate_{};
	bool lastSustain_{};
	bool keyHeld_{};
	std::unique_ptr<X4Resampler_Order7> verticalPositionInterpolator_;
	std::unique_ptr<X4Resampler_Order7> verticalVelocityInterpolator_;
	std::unique_ptr<X4Resampler_Order7> horizontalPositionInterpolator_;
	std::unique_ptr<X4Resampler_Order7> horizontalVelocityInterpolator_;
	std::unique_ptr<X4Resampler_Order7> pickupDecimator_;
};

// Figure 11-8's Q1/Q2 input pair, its Q3/Q4 Darlington active-tone feedback
// amplifier and Q5 follower, plus Figure 11-9, are nonlinear circuit solves.
// Only the post-volume vibrato-feed buffer and lamp/LDR routing are reduced;
// no fitted transfer curve stands in for preamp or power-stage overload. Keep
// measured-data gaps synchronized with the modelling notes and SPICE references.
class ElectricPianoAmplifier
{
public:
	ElectricPianoAmplifier()
		: inputInterpolator_(CreateX2Resampler_Chebychev7()),
		  leftOutputDecimator_(CreateX2Resampler_Chebychev7()),
		  rightOutputDecimator_(CreateX2Resampler_Chebychev7())
	{
	}

	void SetSampleRate(double sampleRate)
	{
		sampleRate_ = std::clamp(sampleRate, 8000.0, 768000.0);
		const double oversampledRate = sampleRate_ *
			X2Resampler_Order7::ResamplingFactor;
		// Figure 11-10 specifies 3000 uF reservoirs. Effective winding,
		// rectifier and 0.5 Ohm emitter-path resistance determine the fast
		// droop; recharge is slower than discharge between mains peaks.
		supplyDischargeCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.0030 * sampleRate_));
		supplyRechargeCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.055 * sampleRate_));
		for (auto& channel : powerChannels_)
			channel.circuit.SetSampleRate(oversampledRate);
		preamp_.SetSampleRate(oversampledRate);
		tonePreamp_.SetSampleRate(oversampledRate);
		driveSmoothingCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.008 * sampleRate_));
		controlSmoothingCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.012 * sampleRate_));
		vibratoLampCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.010 * sampleRate_));
		Reset();
	}

	void Reset()
	{
		for (auto& channel : powerChannels_)
		{
			channel.circuit.Reset(NominalSupplyRail);
			channel.supplyCurrent = 0.0;
		}
		preamp_.Reset();
		tonePreamp_.Reset();
		supplyRail_ = NominalSupplyRail;
		vibratoPhase_ = 0.0;
		vibratoLamp_ = 0.0;
		cachedDrive_ = -1.0;
		controlsInitialized_ = false;
		rightPowerDormant_ = true;
		inputInterpolator_->Reset();
		leftOutputDecimator_->Reset();
		rightOutputDecimator_->Reset();
	}

	// Small-signal impedance of the electrical speaker equivalent used by the
	// power-stage solve. This is deliberately public so tests and future fitting
	// tools can validate the load independently of the nonlinear amplifier.
	// It is not an acoustic cabinet response.
	static std::complex<double> ElectricalLoadImpedance(double frequency)
	{
		return PetersonPowerAmplifier::ElectricalLoadImpedance(frequency);
	}

	static double DriveGainForPosition(double drive)
	{
		drive = std::clamp(std::isfinite(drive) ? drive : 0.32, 0.0, 1.0);
		// This is an explicit extension, not a historical Peterson control. The
		// eased-quadratic dB law keeps the default and middle of the travel in the
		// clean headroom region, then approaches compression progressively. Its
		// reciprocal is
		// applied at Figure 11-8's real interstage volume node, before Figure 11-9:
		// Drive can exercise the preamp without turning the power amp into the
		// compensator or making polyphonic peaks hit its rails. No additional
		// clipper is introduced.
		// With the complete tone-feedback circuit present, the former smoothstep
		// law already delivered 20 dB by 70% travel. Simultaneous hard notes then
		// crossed both preamp headroom boundaries together. A 24 dB eased span
		// retains a deliberately overdriven end stop while keeping moderate chord
		// settings out of that abrupt shared overload region. Blending 40% of a
		// quadratic with 60% smoothstep avoids moving all audible change to the last
		// few degrees of travel.
		const double driveShape = drive * drive * (2.2 - 1.2 * drive);
		const double driveDb = 24.0 * driveShape;
		return std::pow(10.0, driveDb / 20.0);
	}

	static double BassPotPositionForControl(double control)
	{
		control = std::clamp(std::isfinite(control) ? control : 0.5, 0.0, 1.0);
		// The historical 100 k linear pot is highly bunched electrically: most
		// of both cut and boost occurs in the last few degrees at either end.
		// Preserve its centre and end stops, but use a reverse-quarter-power UI
		// law so equal panel travel produces a much more even dB trajectory.
		const double displacement = 2.0 * std::abs(control - 0.5);
		const double mappedDisplacement = std::sqrt(std::sqrt(displacement));
		return 0.5 + std::copysign(0.5 * mappedDisplacement,
			control - 0.5);
	}

	static double TreblePotPositionForControl(double control)
	{
		control = std::clamp(std::isfinite(control) ? control : 0.5, 0.0, 1.0);
		const double displacement = 2.0 * std::abs(control - 0.5);
		const double mappedDisplacement = std::pow(displacement, 0.70);
		return 0.5 + std::copysign(0.5 * mappedDisplacement,
			control - 0.5);
	}

	// Explicit circuit/Rack-domain boundary used by calibration tests. The
	// caller supplies the amplifier's model-domain input (the module currently
	// feeds five times the summed pickup signal).
	static double ModelInputToCircuitVolts(double modelInput, double drive)
	{
		return InputVoltsPerModelUnit * DriveGainForPosition(drive) * modelInput;
	}

	static double CircuitPowerOutputToModel(double circuitVolts)
	{
		return NominalGain * circuitVolts /
			(PowerClosedLoopGain * InputVoltsPerModelUnit);
	}

	double SupplyRailVoltage() const { return supplyRail_; }

	std::array<double, 2> Step(double input,
		const ElectricPianoControls& controls)
	{
		if (!std::isfinite(input))
			input = 0.0;
		auto control = [](double value, double fallback)
		{
			return std::clamp(std::isfinite(value) ? value : fallback, 0.0, 1.0);
		};
		const double drive = control(controls.drive, 0.32);
		const double outputVolume = control(controls.outputVolume, 0.50);
		const double bass = control(controls.amplifierBass, 0.50);
		const double treble = control(controls.amplifierTreble, 0.50);
		const double vibratoSpeed = control(controls.vibratoSpeed, 0.32);
		const double vibratoIntensity = control(controls.vibratoIntensity, 0.0);
		if (!controlsInitialized_)
		{
			smoothedDrive_ = drive;
			smoothedOutputVolume_ = outputVolume;
			smoothedBass_ = bass;
			smoothedTreble_ = treble;
			smoothedVibratoSpeed_ = vibratoSpeed;
			smoothedVibratoIntensity_ = vibratoIntensity;
			controlsInitialized_ = true;
		}
		smoothedDrive_ += driveSmoothingCoefficient_ *
			(drive - smoothedDrive_);
		smoothedOutputVolume_ += controlSmoothingCoefficient_ *
			(outputVolume - smoothedOutputVolume_);
		smoothedBass_ += controlSmoothingCoefficient_ *
			(bass - smoothedBass_);
		smoothedTreble_ += controlSmoothingCoefficient_ *
			(treble - smoothedTreble_);
		smoothedVibratoSpeed_ += controlSmoothingCoefficient_ *
			(vibratoSpeed - smoothedVibratoSpeed_);
		smoothedVibratoIntensity_ += controlSmoothingCoefficient_ *
			(vibratoIntensity - smoothedVibratoIntensity_);
		const bool stereoPowerActive = smoothedVibratoIntensity_ > 1.0e-7;
		if (stereoPowerActive && rightPowerDormant_)
		{
			// Both channels have received identical drive while vibrato was off.
			// Materialize the second circuit state only when stereo routing starts.
			powerChannels_[1] = powerChannels_[0];
			rightPowerDormant_ = false;
		}
		else if (!stereoPowerActive)
		{
			rightPowerDormant_ = true;
		}
		if (smoothedDrive_ != cachedDrive_)
		{
			// The historical volume control did not create an independent clipper.
			// This extended Drive control raises the signal presented to the preamp;
			// its small-signal gain is compensated at the real volume node before the
			// power modules. The eased law is deliberately shallow around the clean
			// default and steeper through the deliberate upper-range transition.
			driveGain_ = DriveGainForPosition(smoothedDrive_);
			makeupGain_ = NominalGain;
			cachedDrive_ = smoothedDrive_;
		}

		// Interpolate before Figure 11-8. The complete preamp -> tone -> Figure
		// 11-9 path runs at 2x; only controls and the slow supply remain at host
		// rate. The pickup/voice model has its own independent 4x path.
		const auto inputs = inputInterpolator_->Upsample(input);
		Eigen::Array<double, X2Resampler_Order7::ResamplingFactor, 1> leftOutputs;
		Eigen::Array<double, X2Resampler_Order7::ResamplingFactor, 1> rightOutputs;
		for (int index = 0; index < X2Resampler_Order7::ResamplingFactor; ++index)
		{
			// The piano/pickup signal is expressed in model units. Convert it to
			// the physical voltage presented to the +25 V Peterson input transistor.
			// Figure 11-8 places that stage before its tone/feedback and volume
			// sections. Extended Drive is compensated at the pre-power volume node.
			const double physicalInput = InputVoltsPerModelUnit * inputs(index);
			// Drive raises the voltage presented to the Figure 11-8 Q1/Q2 circuit.
			// Its reciprocal is applied at Figure 11-8's volume node below, so the
			// control changes preamp headroom rather than becoming a volume control.
			// Q1/Q2 feed the real active feedback network. Do not normalize between
			// those stages: their physical interstage voltage is what determines when
			// the Darlington feedback pair and Q5 leave their small-signal region.
			const double inputStageOutput = preamp_.Step(
				driveGain_ * physicalInput).voltage;
			const double toneOutput = tonePreamp_.Step(inputStageOutput,
				BassPotPositionForControl(smoothedBass_),
				TreblePotPositionForControl(smoothedTreble_)).voltage;
			// Figure 11-8 places its 100 kOhm volume control before the channel
			// output buffers and Figure 11-9. Use that physical node for Drive's
			// reciprocal gain. Putting the compensation after Figure 11-9 preserved
			// the output level but still drove chord peaks into the power rails.
			const double preampOutput = toneOutput /
				(driveGain_ * NominalFullPreampGain);

			// In Figure 11-8 the optical network routes the preamp signal to two
			// outputs, each of which feeds its own Figure 11-9 power module. Keep
			// those power stages independent so their crossover and overload history
			// follow the channel currents rather than a post-distortion pan law.
			const double pan = smoothedVibratoIntensity_ * vibratoLamp_;
			const double leftDrive = std::sqrt(std::max(0.0, 1.0 + pan));
			const double rightDrive = std::sqrt(std::max(0.0, 1.0 - pan));
			const double leftVoltage = ProcessPowerModule(preampOutput * leftDrive,
				powerChannels_[0], supplyRail_);
			double rightVoltage = leftVoltage;
			if (stereoPowerActive)
			{
				rightVoltage = ProcessPowerModule(preampOutput * rightDrive,
					powerChannels_[1], supplyRail_);
			}
			const double modelScale = makeupGain_ /
				(PowerClosedLoopGain * InputVoltsPerModelUnit);
			leftOutputs(index) = modelScale * leftVoltage;
			rightOutputs(index) = modelScale * rightVoltage;
		}
		const double left = leftOutputDecimator_->Downsample(leftOutputs);
		const double right = rightOutputDecimator_->Downsample(rightOutputs);
		const double totalCurrent = stereoPowerActive ?
			powerChannels_[0].supplyCurrent + powerChannels_[1].supplyCurrent :
			2.0 * powerChannels_[0].supplyCurrent;
		const double loadedRail = std::clamp(NominalSupplyRail -
			SupplyResistance * totalCurrent, MinimumSupplyRail,
			NominalSupplyRail);
		const double supplyCoefficient = loadedRail < supplyRail_ ?
			supplyDischargeCoefficient_ : supplyRechargeCoefficient_;
		supplyRail_ += supplyCoefficient * (loadedRail - supplyRail_);

		const double vibratoFrequency = 1.5 * std::pow(8.0,
			smoothedVibratoSpeed_);
		vibratoPhase_ += vibratoFrequency / sampleRate_;
		if (vibratoPhase_ >= 1.0)
			vibratoPhase_ -= 1.0;
		const double sine = std::sin(6.2831853071795864769 * vibratoPhase_);
		const double lampTarget = std::tanh(1.6 * sine) / std::tanh(1.6);
		vibratoLamp_ += vibratoLampCoefficient_ *
			(lampTarget - vibratoLamp_);
		const double outputGain = 2.0 * smoothedOutputVolume_;
		return {
			outputGain * (std::isfinite(left) ? left : 0.0),
			outputGain * (std::isfinite(right) ? right : 0.0)};
	}

private:
	struct PowerChannel
	{
		PetersonPowerAmplifier circuit;
		double supplyCurrent{};
	};

	double ProcessPowerModule(double input, PowerChannel& channel,
		double supplyRail)
	{
		const auto result = channel.circuit.Step(input, supplyRail);
		channel.supplyCurrent = result.positiveRailCurrent +
			result.negativeRailCurrent;
		return result.voltage;
	}

	// With the module's 5x pickup-sum feed, 0.13 V/model-unit makes a calibrated
	// hard single note present about 0.42 Vpk to the preamp at maximum Drive.
	// After nominal preamp-gain normalization, Figure 11-9 retains the existing
	// voltage calibration while real preamp compression reduces its overload
	// demand. The default remains in both circuits' feedback-linear range.
	static constexpr double InputVoltsPerModelUnit = 0.13;
	static constexpr double PowerClosedLoopGain = 57.0; // 1 + 5.6k / 100
	// Small-signal gain of Q1/Q2 plus the centred Figure 11-8 nonlinear tone
	// feedback stage at 250 Hz. This is a unit calibration only: there is no
	// normalization between the two physical stages, so interstage overload and
	// gain compression remain in the circuit path.
	static constexpr double NominalFullPreampGain = 6.3230;
	static constexpr double NominalSupplyRail = 35.0;
	static constexpr double MinimumSupplyRail = 27.0;
	static constexpr double SupplyResistance = 0.72;
	static constexpr double NominalGain = 1.12;

	double sampleRate_ = 48000.0;
	PetersonPreAmplifier preamp_{};
	PetersonTonePreAmplifier tonePreamp_{};
	std::array<PowerChannel, 2> powerChannels_{};
	double supplyRail_ = NominalSupplyRail;
	double vibratoPhase_{};
	double vibratoLamp_{};
	double cachedDrive_ = -1.0;
	double smoothedDrive_{};
	double driveSmoothingCoefficient_ = 0.0035;
	double controlSmoothingCoefficient_ = 0.002;
	double vibratoLampCoefficient_ = 0.002;
	double supplyDischargeCoefficient_ = 0.0069;
	double supplyRechargeCoefficient_ = 0.00038;
	double driveGain_ = 1.0;
	double makeupGain_ = NominalGain;
	double smoothedOutputVolume_ = 0.50;
	double smoothedBass_ = 0.50;
	double smoothedTreble_ = 0.50;
	double smoothedVibratoSpeed_ = 0.32;
	double smoothedVibratoIntensity_{};
	bool controlsInitialized_{};
	bool rightPowerDormant_{true};
	std::unique_ptr<X2Resampler_Order7> inputInterpolator_;
	std::unique_ptr<X2Resampler_Order7> leftOutputDecimator_;
	std::unique_ptr<X2Resampler_Order7> rightOutputDecimator_;
};

} // namespace tfdsp
