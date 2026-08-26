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
	double touch = 0.5;
	double dynamics = 0.75;
	double body = 0.62;
	double bell = 0.52;
	double coupling = 0.58;
	double hammer = 0.52;
	double tone = 0.55;
	double gap = 0.48;
	double decay = 0.58;
	double release = 0.24;
	double mechanics = 0.18;
	double drive = 0.32;
};

inline constexpr double ElectricPianoReferenceFrequency =
	261.6255653005986;

struct ElectricPianoKeyProfile
{
	double fundamentalHz{};
	double modalMassRatio{};
	double displacementPerImpulse{};
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
	// to equalize small-signal velocity; it cancels modal-mass level variation
	// but deliberately leaves the larger low-key displacement intact.
	const double modalMassRatio = std::pow(frequencyRatio, -0.5);
	const double midiNote = 60.0 + 12.0 * boundedPitch;
	return {
		ElectricPianoReferenceFrequency * frequencyRatio,
		modalMassRatio,
		1.0 / (modalMassRatio * frequencyRatio),
		modalMassRatio,
		std::clamp((midiNote - 28.0) / 72.0, 0.0, 1.0)};
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
		: positionInterpolator_(CreateX4Resampler_Cheby7()),
		  velocityInterpolator_(CreateX4Resampler_Cheby7()),
		  pickupDecimator_(CreateX4Resampler_Cheby7())
	{
	}

	void SetNoiseSeed(std::uint32_t seed)
	{
		noiseState_ = seed != 0 ? seed : 0x6d2b79f5u;
	}

	void SetSampleRate(double sampleRate)
	{
		sampleRate_ = std::clamp(sampleRate, 8000.0, 768000.0);
		pickupLowPass_ = 0.0;
		hammerNoiseDecay_ = std::exp(-1.0 / (0.010 * sampleRate_));
		keyReleaseNoiseDecay_ = std::exp(-1.0 / (0.015 * sampleRate_));
		damperNoiseDecay_ = std::exp(-1.0 / (0.025 * sampleRate_));
		controlSmoothingCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.006 * sampleRate_));
		coefficientsDirty_ = true;
		cachedTone_ = -1.0;
		positionInterpolator_->Reset();
		velocityInterpolator_->Reset();
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
		noiseFilter_ = 0.0;
		energy_ = 0.0;
		controlsInitialized_ = false;
		coefficientsDirty_ = true;
		timbreDirty_ = true;
		positionInterpolator_->Reset();
		velocityInterpolator_->Reset();
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
			keyReleaseNoise_ = 0.08 + 0.18 * latchedVelocity_;
			if (!sustain)
				damperNoise_ = 0.06 + 0.22 * latchedVelocity_;
		}
		if (lastSustain_ && !sustain && !keyHeld_)
			damperNoise_ = std::max(damperNoise_,
				0.06 + 0.22 * latchedVelocity_);
		lastGate_ = gate;
		lastSustain_ = sustain;
		if (gate)
			keyHeld_ = true;

		const auto activeControls = SmoothActiveControls(controls);
		const double boundedPitch = std::clamp(pitchVolts, -6.0, 6.0);
		if (boundedPitch != cachedPitch_ || coefficientsDirty_)
		{
			currentFundamental_ = ElectricPianoReferenceFrequency *
				std::exp2(boundedPitch);
			cachedPitch_ = boundedPitch;
		}
		const bool damped = !keyHeld_ && !sustain;
		RefreshModeCoefficients(currentFundamental_, damped, activeControls);
		RefreshTimbreCoefficients(activeControls);

		double position = 0.0;
		double normalizedVelocity = 0.0;
		energy_ = 0.0;
		for (std::size_t index = 0; index < modes_.size(); ++index)
		{
			if (!modeActive_[index])
				continue;
			const double oldReal = modes_[index].real;
			const double oldImaginary = modes_[index].imaginary;
			modes_[index].real = modeRealCoefficient_[index] * oldReal -
				modeImaginaryCoefficient_[index] * oldImaginary;
			modes_[index].imaginary = modeImaginaryCoefficient_[index] * oldReal +
				modeRealCoefficient_[index] * oldImaginary;

			const double outputWeight = modeOutputWeight_[index] *
				modeBandlimitGain_[index];
			position += outputWeight * modes_[index].real;
			normalizedVelocity -= outputWeight * frequencyVelocityScale_ *
				modeRatio_[index] * modes_[index].imaginary;
			energy_ += modes_[index].real * modes_[index].real +
				modes_[index].imaginary * modes_[index].imaginary;
		}

		const auto pickupPositions = positionInterpolator_->Upsample(
			TineDisplacementScale * position);
		const auto pickupVelocities = velocityInterpolator_->Upsample(
			normalizedVelocity);
		Eigen::Array<double, PickupOversamplingFactor, 1> pickupValues;
		for (int index = 0; index < PickupOversamplingFactor; ++index)
		{
			const double movingGradient = MagneticGradient(
				pickupOffset_ + pickupPositions(index), pickupGap_);
			const double pickupSensitivity = movingGradient *
				inverseReferenceGradient_;
			double pickup = 0.115 * keyPickupSensitivity_ *
				pickupVelocities(index) *
				std::clamp(pickupSensitivity, -8.0, 8.0);
			pickupValues(index) = std::tanh(pickupDrive_ * pickup) /
				pickupDrive_;
		}
		const double pickup = pickupDecimator_->Downsample(pickupValues);
		pickupLowPass_ += pickupLowPassCoefficient_ *
			(pickup - pickupLowPass_);

		const double white = WhiteNoise();
		noiseFilter_ += 0.34 * (white - noiseFilter_);
		const double mechanicalSignal = mechanicsLevel_ *
			(0.34 * hammerNoise_ * noiseFilter_ +
				0.12 * keyReleaseNoise_ * noiseFilter_ +
				0.18 * damperNoise_ * (white - 0.6 * noiseFilter_));
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

		const double result = pickupLowPass_ + mechanicalSignal;
		return std::isfinite(result) ? result : 0.0;
	}

	double Energy() const { return energy_; }
	double Activity() const
	{
		return energy_ + hammerNoise_ * hammerNoise_ +
			keyReleaseNoise_ * keyReleaseNoise_ +
			damperNoise_ * damperNoise_ + pickupLowPass_ * pickupLowPass_;
	}
	bool IsAudible() const { return Activity() > 1.0e-12; }
	bool GateHigh() const { return lastGate_; }
	double NotePitch() const { return notePitch_; }

private:
	struct Mode
	{
		double real{};
		double imaginary{};
	};

	static double Clamp01(double value)
	{
		return std::clamp(std::isfinite(value) ? value : 0.0, 0.0, 1.0);
	}

	ElectricPianoControls SmoothActiveControls(
		const ElectricPianoControls& controls)
	{
		if (!controlsInitialized_)
		{
			smoothedControls_ = controls;
			controlsInitialized_ = true;
		}
		auto smooth = [&](double& current, double target)
		{
			target = Clamp01(target);
			current += controlSmoothingCoefficient_ * (target - current);
		};
		smooth(smoothedControls_.body, controls.body);
		smooth(smoothedControls_.bell, controls.bell);
		smooth(smoothedControls_.coupling, controls.coupling);
		smooth(smoothedControls_.tone, controls.tone);
		smooth(smoothedControls_.gap, controls.gap);
		smooth(smoothedControls_.decay, controls.decay);
		smooth(smoothedControls_.release, controls.release);
		smooth(smoothedControls_.mechanics, controls.mechanics);
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
		const double touch = Clamp01(controls.touch);
		const double dynamics = Clamp01(controls.dynamics);
		const double body = Clamp01(controls.body);
		const double bell = Clamp01(controls.bell);
		const double coupling = Clamp01(controls.coupling);
		const double hammer = Clamp01(controls.hammer);
		const double gamma = std::exp2(1.0 - 2.0 * touch);
		const double curvedVelocity = velocity > 0.0 ?
			std::pow(velocity, gamma) : 0.0;
		latchedVelocity_ = curvedVelocity;
		const double compressedVelocity = curvedVelocity /
			(0.28 + 0.72 * curvedVelocity);
		const double dynamicAmplitude = (1.0 - dynamics) *
			compressedVelocity + dynamics * curvedVelocity;
		const double hammerHardness = 0.38 + 0.82 * hammer +
			0.48 * keyPosition_ + 0.70 * curvedVelocity;
		const double bodyImpulse = dynamicAmplitude *
			(0.72 + 0.42 * body);
		const double bellImpulse = dynamicAmplitude * bell *
			(0.055 + 0.78 * curvedVelocity * curvedVelocity) * hammerHardness;

		modes_[0].imaginary += bodyImpulse * key.displacementPerImpulse;
		modes_[1].imaginary += bodyImpulse *
			(0.06 + 0.54 * coupling) * key.displacementPerImpulse /
			std::max(0.1, 1.0 - (0.00015 + 0.00135 * coupling));
		modes_[2].imaginary += bodyImpulse *
			(0.015 + 0.10 * coupling) * (0.55 + 0.45 * curvedVelocity) *
			key.displacementPerImpulse /
			(1.0 + 0.00040 + 0.00410 * coupling);
		const double rolloffStep = std::exp(-(0.56 + 0.16 * (1.0 - hammer)));
		double rolloff = 1.0;
		for (std::size_t index = 3; index < modes_.size(); ++index)
		{
			const double strikeNotch = std::abs(std::sin(
				ModeRatios[index] * (0.115 - 0.052 * keyPosition_)));
			modes_[index].imaginary += bellImpulse * rolloff *
				(0.32 + 0.68 * strikeNotch) /
				ModeRatios[index] * key.displacementPerImpulse;
			rolloff *= rolloffStep;
		}

		hammerNoise_ = dynamicAmplitude * (0.28 + 0.72 * curvedVelocity) *
			(0.72 + 0.35 * keyPosition_);
		keyReleaseNoise_ = 0.0;
		damperNoise_ = 0.0;
		keyHeld_ = true;
		energy_ = 1.0;
	}

	void RefreshTimbreCoefficients(const ElectricPianoControls& controls)
	{
		const double body = Clamp01(controls.body);
		const double bell = Clamp01(controls.bell);
		const double coupling = Clamp01(controls.coupling);
		const double gap = Clamp01(controls.gap);
		const double tone = Clamp01(controls.tone);
		const double mechanics = Clamp01(controls.mechanics);
		if (!timbreDirty_ && body == cachedBodyWeight_ &&
			bell == cachedBellWeight_ && coupling == cachedCouplingWeight_ &&
			gap == cachedGap_ && tone == cachedTone_ &&
			mechanics == mechanicsLevel_)
			return;

		modeOutputWeight_[0] = 0.64 + 0.48 * body;
		modeOutputWeight_[1] = (0.20 + 0.56 * body) *
			(0.20 + 0.80 * coupling);
		modeOutputWeight_[2] = 0.68 * (0.15 + 0.85 * coupling);
		for (std::size_t index = 3; index < modeOutputWeight_.size(); ++index)
			modeOutputWeight_[index] = 0.44 + 0.92 * bell;

		pickupGap_ = 1.40 - 0.88 * gap;
		pickupOffset_ = 0.22 + 0.86 * tone;
		const double cutoff = std::min(0.42 * sampleRate_,
			6200.0 + 6800.0 * tone);
		pickupLowPassCoefficient_ = 1.0 - std::exp(
			-TwoPi * cutoff / sampleRate_);
		mechanicsLevel_ = mechanics;
		cachedBodyWeight_ = body;
		cachedBellWeight_ = bell;
		cachedCouplingWeight_ = coupling;
		cachedGap_ = gap;
		cachedTone_ = tone;
		timbreDirty_ = false;
	}

	void RefreshModeCoefficients(double fundamental, bool damped,
		const ElectricPianoControls& controls)
	{
		const double decay = Clamp01(controls.decay);
		const double bell = Clamp01(controls.bell);
		const double coupling = Clamp01(controls.coupling);
		const double release = Clamp01(controls.release);
		if (!coefficientsDirty_ && fundamental == cachedFundamental_ &&
			damped == cachedDamped_ && decay == cachedDecay_ &&
			bell == cachedBell_ && coupling == cachedCoupling_ &&
			release == cachedRelease_)
			return;

		const double releaseSeconds = 0.012 * std::pow(100.0, release);
		const double decayScale = std::pow(4.0, 2.0 * (decay - 0.5));
		const double bellDecay = 0.72 + 0.78 * bell;
		frequencyVelocityScale_ = fundamental /
			ElectricPianoReferenceFrequency;
		for (std::size_t index = 0; index < modes_.size(); ++index)
		{
			modeRatio_[index] = ModeRatios[index];
			if (index == 1)
				modeRatio_[index] = 1.0 - (0.00015 + 0.00135 * coupling);
			else if (index == 2)
				modeRatio_[index] = 1.0 + (0.00040 + 0.00410 * coupling);
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
				continue;
			}
			double lifetime = BaseDecaySeconds[index] * decayScale;
			if (index >= 3)
				lifetime *= bellDecay;
			lifetime *= 1.18 - 0.50 * keyPosition_;
			if (damped)
				lifetime = std::min(lifetime, releaseSeconds *
					(index >= 3 ? 0.55 : 1.0));
			const double radius = std::exp(-1.0 /
				(std::max(0.002, lifetime) * sampleRate_));
			const double angle = TwoPi * frequency / sampleRate_;
			modeRealCoefficient_[index] = radius * std::cos(angle);
			modeImaginaryCoefficient_[index] = radius * std::sin(angle);
		}
		cachedFundamental_ = fundamental;
		cachedDecay_ = decay;
		cachedBell_ = bell;
		cachedCoupling_ = coupling;
		cachedRelease_ = release;
		cachedDamped_ = damped;
		coefficientsDirty_ = false;
	}

	static double MagneticGradient(double position, double gap)
	{
		const double squaredDistance = gap * gap + position * position;
		return position / (squaredDistance * std::sqrt(squaredDistance));
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
	static constexpr double TineDisplacementScale = 0.34;
	static constexpr int PickupOversamplingFactor =
		X4Resampler_Order7::ResamplingFactor;
	static constexpr double ReferencePickupGap = 1.40 - 0.88 * 0.48;
	static constexpr double ReferencePickupOffset = 0.22 + 0.86 * 0.55;
	static constexpr std::array<double, 8> ModeRatios{
		1.0, 0.99915, 1.0027, 6.267, 17.55, 34.39, 56.84, 83.0};
	static constexpr std::array<double, 8> BaseDecaySeconds{
		4.6, 6.4, 2.3, 0.30, 0.18, 0.105, 0.062, 0.038};

	std::array<Mode, 8> modes_{};
	std::array<double, 8> modeRealCoefficient_{};
	std::array<double, 8> modeImaginaryCoefficient_{};
	std::array<double, 8> modeOutputWeight_{};
	std::array<double, 8> modeRatio_{};
	std::array<double, 8> modeBandlimitGain_{};
	std::array<bool, 8> modeActive_{};
	double sampleRate_ = 48000.0;
	double latchedVelocity_{};
	double keyPosition_ = 0.5;
	double notePitch_{};
	double frequencyVelocityScale_ = 1.0;
	double keyPickupSensitivity_ = 1.0;
	double pickupLowPass_{};
	double pickupLowPassCoefficient_{};
	double pickupDrive_ = 0.85;
	double pickupGap_ = 1.0;
	double pickupOffset_ = 0.7;
	double inverseReferenceGradient_ = 1.0 / std::abs(MagneticGradient(
		ReferencePickupOffset, ReferencePickupGap));
	double mechanicsLevel_{};
	double hammerNoiseDecay_ = 0.9979;
	double keyReleaseNoiseDecay_ = 0.9986;
	double damperNoiseDecay_ = 0.9992;
	double hammerNoise_{};
	double keyReleaseNoise_{};
	double damperNoise_{};
	double noiseFilter_{};
	double controlSmoothingCoefficient_ = 0.0035;
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
	double cachedCouplingWeight_ = -1.0;
	double cachedGap_ = -1.0;
	double cachedTone_ = -1.0;
	std::uint32_t noiseState_ = 0x6d2b79f5u;
	bool cachedDamped_{};
	bool coefficientsDirty_ = true;
	bool timbreDirty_ = true;
	bool controlsInitialized_{};
	bool lastGate_{};
	bool lastSustain_{};
	bool keyHeld_{};
	std::unique_ptr<X4Resampler_Order7> positionInterpolator_;
	std::unique_ptr<X4Resampler_Order7> velocityInterpolator_;
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
		driveSmoothingCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.006 * sampleRate_));
		Reset();
	}

	void Reset()
	{
		previousInput_ = 0.0;
		previousHighPass_ = 0.0;
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
			gain_ = std::pow(10.0, (-4.0 + 22.0 * smoothedDrive_) / 20.0);
			bias_ = 0.035 + 0.085 * smoothedDrive_;
			biasTanh_ = std::tanh(bias_);
			normalization_ = 0.72 * (1.0 - biasTanh_ * biasTanh_);
			cachedDrive_ = smoothedDrive_;
		}
		const auto inputs = inputInterpolator_->Upsample(input);
		Eigen::Array<double, X4Resampler_Order7::ResamplingFactor, 1> outputs;
		for (int index = 0; index < X4Resampler_Order7::ResamplingFactor; ++index)
		{
			const double biased = std::tanh(
				0.72 * gain_ * inputs(index) + bias_) - biasTanh_;
			outputs(index) = biased / std::max(0.1, normalization_);
		}
		const double shaped = outputDecimator_->Downsample(outputs);
		const double highPassed = highPassCoefficient_ *
			(previousHighPass_ + shaped - previousInput_);
		previousInput_ = shaped;
		previousHighPass_ = highPassed;
		return std::isfinite(highPassed) ? highPassed : 0.0;
	}

private:
	double sampleRate_ = 48000.0;
	double highPassCoefficient_ = 0.9993;
	double previousInput_{};
	double previousHighPass_{};
	double cachedDrive_ = -1.0;
	double smoothedDrive_{};
	double driveSmoothingCoefficient_ = 0.0035;
	double gain_ = 1.0;
	double bias_{};
	double biasTanh_{};
	double normalization_ = 0.72;
	bool driveInitialized_{};
	std::unique_ptr<X4Resampler_Order7> inputInterpolator_;
	std::unique_ptr<X4Resampler_Order7> outputDecimator_;
};

} // namespace tfdsp
