#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>

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
	double gap = 0.48;
	double decay = 0.50;
	double release = 0.24;
	double mechanics = 0.18;
	double drive = 0.32;
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
		std::pow(frequencyRatio, -0.08) *
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
			const double pickup = 0.135 * keyPickupSensitivity_ *
				pickupAlignmentGain_ * inverseReferenceGradient_ * emf;
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
		smooth(smoothedControls_.gap, controls.gap,
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
		const double gap = Clamp01(controls.gap);
		const double tone = Clamp01(controls.tone);
		const double mechanics = Clamp01(controls.mechanics);
		if (!timbreDirty_ && body == cachedBodyWeight_ &&
			bell == cachedBellWeight_ && gap == cachedGap_ && tone == cachedTone_ &&
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
		// fundamental relative to curvature-generated harmonics. GAP controls the
		// independent front-to-back distance and therefore level and dynamics.
		pickupGap_ = 1.32 - 0.78 * gap;
		pickupVerticalOffset_ = 0.988 - 0.90 * tone +
			0.055 * keyPosition_;
		pickupHorizontalOffset_ = 0.10 + 0.035 * keyPosition_;
		const auto alignmentGradient = MagneticFluxGradient(
			pickupVerticalOffset_, pickupHorizontalOffset_, ReferencePickupGap);
		const double alignmentMagnitude = std::sqrt(
			alignmentGradient[0] * alignmentGradient[0] +
			alignmentGradient[1] * alignmentGradient[1]);
		pickupAlignmentGain_ = std::clamp(ReferenceGradientMagnitude() /
			std::max(1.0e-6, alignmentMagnitude), 0.50, 4.0);
		const double cutoff = std::min(0.42 * sampleRate_, 16500.0);
		pickupLowPassCoefficient_ = 1.0 - std::exp(
			-TwoPi * cutoff / sampleRate_);
		// Preserve a wide sound-design range while making the calibrated default
		// a quiet mechanical contribution rather than a parallel percussion layer.
		mechanicsLevel_ = 0.060 * mechanics * (0.30 + 0.70 * mechanics);
		cachedBodyWeight_ = body;
		cachedBellWeight_ = bell;
		cachedGap_ = gap;
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
		if (!coefficientsDirty_ && fundamental == cachedFundamental_ &&
			damped == cachedDamped_ && decay == cachedDecay_ &&
			bell == cachedBell_ && coupling == cachedCoupling_ &&
			release == cachedRelease_)
			return;

		const double releaseSeconds = 0.012 * std::pow(100.0, release);
		const double decayScale = std::pow(4.0, 2.0 * (decay - 0.5));
		const double bellDecay = 0.82 + 0.38 * bell;
		frequencyVelocityScale_ = fundamental /
			ElectricPianoReferenceFrequency;

		const bool transformCoupledState = coupledForkInitialized_ &&
			(fundamental != cachedFundamental_ || coupling != cachedCoupling_) &&
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

		const auto coupledFork = MakeElectricPianoCoupledForkProfile(coupling,
			keyPosition_);
		const double keyboardDecayScale = 1.18 - 0.50 * keyPosition_;
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
			if (index < 2)
			{
				modeRatio_[index] = coupledFork.frequencyRatios[index];
				coupledToneBarRatio_[index] =
					coupledFork.toneBarDisplacementRatios[index];
				modeInverseMass_[index] =
					coupledFork.inverseModalMassRatios[index] /
					std::max(1.0e-9, tineModalMass_);
				contactModeShape_[index] = 0.82;
			}
			else if (index == 2)
				modeRatio_[index] = 1.0018 +
					0.0010 * (keyPosition_ - 0.5);
			else
				modeRatio_[index] = KeyboardModeRatio(index, keyPosition_);
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
			double lifetime;
			if (index < 2)
			{
				const double lossRate = 1.0 / intrinsicForkLifetime +
					supportLossRate *
					coupledFork.supportReactionLossFactors[index];
				lifetime = decayScale / std::max(1.0e-9, lossRate);
			}
			else
				lifetime = HigherModeBaseDecaySeconds[index - 2] * decayScale *
					keyboardDecayScale;
			if (index >= 3)
				lifetime *= bellDecay;
			if (damped)
				lifetime = std::min(lifetime, releaseSeconds *
					(index >= 3 ? 0.55 : 1.0));
			const double radius = std::exp(-1.0 /
				(std::max(0.002, lifetime) * sampleRate_));
			const double angle = TwoPi * frequency / sampleRate_;
			modeRealCoefficient_[index] = radius * std::cos(angle);
			modeImaginaryCoefficient_[index] = radius * std::sin(angle);
			const double substepRadius = std::exp(-1.0 /
				(std::max(0.002, lifetime) * sampleRate_ *
					static_cast<double>(ContactOversamplingFactor)));
			const double substepAngle = angle /
				static_cast<double>(ContactOversamplingFactor);
			modeSubstepRealCoefficient_[index] = substepRadius *
				std::cos(substepAngle);
			modeSubstepImaginaryCoefficient_[index] = substepRadius *
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
		// Five-point quadrature approximates the finite rounded pole face. This
		// preserves the inexpensive spatial transfer-function approach while
		// avoiding the unrealistically sharp curvature of a magnetic point source.
		constexpr double HorizontalFieldScale = 0.72;
		constexpr std::array<double, 5> PolePositions{
			-0.34, -0.17, 0.0, 0.17, 0.34};
		constexpr std::array<double, 5> PoleWeights{
			1.0 / 12.0, 3.0 / 12.0, 4.0 / 12.0, 3.0 / 12.0, 1.0 / 12.0};
		double verticalGradient = 0.0;
		double horizontalGradient = 0.0;
		for (std::size_t index = 0; index < PolePositions.size(); ++index)
		{
			const double verticalDistance = vertical - PolePositions[index];
			const double squaredDistance = gap * gap +
				verticalDistance * verticalDistance +
				HorizontalFieldScale * horizontal * horizontal + 0.018;
			const double inverseCubedDistance = 1.0 /
				(squaredDistance * std::sqrt(squaredDistance));
			verticalGradient += PoleWeights[index] * verticalDistance *
				inverseCubedDistance;
			horizontalGradient += PoleWeights[index] * HorizontalFieldScale *
				horizontal * inverseCubedDistance;
		}
		return {verticalGradient, horizontalGradient};
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
	static constexpr double TineDisplacementScale = 0.56;
	static constexpr int ContactOversamplingFactor = 16;
	static constexpr double MinimumHammerContactVelocity = 1.0;
	static constexpr int PickupOversamplingFactor =
		X4Resampler_Order7::ResamplingFactor;
	static constexpr double ReferencePickupGap = 1.32 - 0.78 * 0.48;
	static constexpr double ReferencePickupVertical = 0.988 - 0.90 * 0.55 +
		0.055 * ((60.0 - 28.0) / 72.0);
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
	std::array<double, 8> modeSubstepRealCoefficient_{};
	std::array<double, 8> modeSubstepImaginaryCoefficient_{};
	std::array<double, 8> modeAngularFrequency_{};
	std::array<double, 8> modeOutputWeight_{};
	std::array<double, 8> modeHorizontalWeight_{};
	std::array<double, 8> modeRatio_{};
	std::array<double, 8> modeBandlimitGain_{};
	std::array<bool, 8> modeActive_{};
	std::array<double, 8> contactModeShape_{};
	std::array<double, 8> modeInverseMass_{};
	std::array<double, 2> coupledToneBarRatio_{};
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
	double cachedGap_ = -1.0;
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

class ElectricPianoAmplifier
{
public:
	ElectricPianoAmplifier()
		: inputInterpolator_(CreateX4Resampler_Cheby7()),
		  outputDecimator_(CreateX4Resampler_Cheby7())
	{
	}

	void SetSampleRate(double sampleRate)
	{
		sampleRate_ = std::clamp(sampleRate, 8000.0, 768000.0);
		const double rc = 1.0 / (6.2831853071795864769 * 5.0);
		highPassCoefficient_ = rc / (rc + 1.0 / sampleRate_);
		inputBassCoefficient_ = OnePoleCoefficient(180.0, sampleRate_);
		inputBodyCoefficient_ = OnePoleCoefficient(1250.0, sampleRate_);
		cabinetLowCoefficient_ = OnePoleCoefficient(190.0, sampleRate_);
		cabinetBodyCoefficient_ = OnePoleCoefficient(1150.0, sampleRate_);
		cabinetHighCoefficient_ = OnePoleCoefficient(10500.0, sampleRate_);
		detectorAttackCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.0035 * sampleRate_));
		detectorReleaseCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.140 * sampleRate_));
		stageSmoothingCoefficient_ = OnePoleCoefficient(10500.0,
			sampleRate_ * X4Resampler_Order7::ResamplingFactor);
		driveSmoothingCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.006 * sampleRate_));
		Reset();
	}

	void Reset()
	{
		previousInput_ = 0.0;
		previousHighPass_ = 0.0;
		inputBass_ = 0.0;
		inputBody_ = 0.0;
		cabinetLow_ = 0.0;
		cabinetBody_ = 0.0;
		cabinetHigh_ = 0.0;
		detectorEnvelope_ = 0.0;
		stageMemory_ = 0.0;
		cachedDrive_ = -1.0;
		driveInitialized_ = false;
		inputInterpolator_->Reset();
		outputDecimator_->Reset();
	}

	double Step(double input, double drive)
	{
		if (!std::isfinite(input))
			return 0.0;
		drive = std::clamp(std::isfinite(drive) ? drive : 0.0, 0.0, 1.0);
		if (!driveInitialized_)
		{
			smoothedDrive_ = drive;
			driveInitialized_ = true;
		}
		smoothedDrive_ += driveSmoothingCoefficient_ *
			(drive - smoothedDrive_);
		if (smoothedDrive_ != cachedDrive_)
		{
			gain_ = std::pow(10.0, (-6.0 + 22.0 * smoothedDrive_) / 20.0);
			bias_ = 0.015 + 0.045 * smoothedDrive_;
			biasShape_ = AsymmetricSoftClip(bias_);
			feedback_ = 0.10 + 0.12 * smoothedDrive_;
			supplySensitivity_ = 0.08 + 0.30 * smoothedDrive_;
			cachedDrive_ = smoothedDrive_;
		}

		// Approximate the harp/load and preamplifier frequency response before
		// overload. A small low-mid lift and restrained top end keep the amplifier
		// from being used as a broadband distortion pedal.
		inputBass_ += inputBassCoefficient_ * (input - inputBass_);
		inputBody_ += inputBodyCoefficient_ * (input - inputBody_);
		const double highBand = input - inputBody_;
		const double voicedInput = input + 0.10 * inputBass_ +
			0.14 * (inputBody_ - inputBass_) - 0.04 * highBand;

		const double level = std::abs(voicedInput);
		const double detectorCoefficient = level > detectorEnvelope_ ?
			detectorAttackCoefficient_ : detectorReleaseCoefficient_;
		detectorEnvelope_ += detectorCoefficient *
			(level - detectorEnvelope_);
		const double supply = 1.0 /
			(1.0 + supplySensitivity_ * detectorEnvelope_);

		const auto inputs = inputInterpolator_->Upsample(voicedInput);
		Eigen::Array<double, X4Resampler_Order7::ResamplingFactor, 1> outputs;
		for (int index = 0; index < X4Resampler_Order7::ResamplingFactor; ++index)
		{
			const double stageInput = gain_ * supply * inputs(index) + bias_ -
				feedback_ * stageMemory_;
			const double nonlinear = AsymmetricSoftClip(stageInput) - biasShape_;
			stageMemory_ += stageSmoothingCoefficient_ *
				(nonlinear - stageMemory_);
			outputs(index) = 0.78 * nonlinear + 0.22 * stageMemory_;
		}
		const double shaped = outputDecimator_->Downsample(outputs);
		cabinetLow_ += cabinetLowCoefficient_ * (shaped - cabinetLow_);
		cabinetBody_ += cabinetBodyCoefficient_ * (shaped - cabinetBody_);
		cabinetHigh_ += cabinetHighCoefficient_ * (shaped - cabinetHigh_);
		const double cabinet = cabinetHigh_ + 0.08 * cabinetLow_ +
			0.14 * (cabinetBody_ - cabinetLow_);
		const double highPassed = highPassCoefficient_ *
			(previousHighPass_ + cabinet - previousInput_);
		previousInput_ = cabinet;
		previousHighPass_ = highPassed;
		return std::isfinite(highPassed) ? highPassed : 0.0;
	}

private:
	static double OnePoleCoefficient(double frequency, double sampleRate)
	{
		return 1.0 - std::exp(-6.2831853071795864769 * frequency /
			sampleRate);
	}

	static double AsymmetricSoftClip(double value)
	{
		if (value >= 0.0)
			return value / (1.0 + 0.55 * value);
		return value / (1.0 - 0.72 * value);
	}

	double sampleRate_ = 48000.0;
	double highPassCoefficient_ = 0.9993;
	double previousInput_{};
	double previousHighPass_{};
	double inputBass_{};
	double inputBody_{};
	double cabinetLow_{};
	double cabinetBody_{};
	double cabinetHigh_{};
	double detectorEnvelope_{};
	double stageMemory_{};
	double cachedDrive_ = -1.0;
	double smoothedDrive_{};
	double driveSmoothingCoefficient_ = 0.0035;
	double inputBassCoefficient_ = 0.023;
	double inputBodyCoefficient_ = 0.15;
	double cabinetLowCoefficient_ = 0.025;
	double cabinetBodyCoefficient_ = 0.14;
	double cabinetHighCoefficient_ = 0.63;
	double detectorAttackCoefficient_ = 0.006;
	double detectorReleaseCoefficient_ = 0.00015;
	double stageSmoothingCoefficient_ = 0.29;
	double gain_ = 1.0;
	double bias_{};
	double biasShape_{};
	double feedback_ = 0.1;
	double supplySensitivity_ = 0.08;
	bool driveInitialized_{};
	std::unique_ptr<X4Resampler_Order7> inputInterpolator_;
	std::unique_ptr<X4Resampler_Order7> outputDecimator_;
};

} // namespace tfdsp
