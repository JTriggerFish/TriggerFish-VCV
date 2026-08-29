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

inline constexpr double ElectricPianoDefaultDecay = 0.50;
inline constexpr double ElectricPianoDefaultRelease = 0.12;

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
	double decay = ElectricPianoDefaultDecay;
	double release = ElectricPianoDefaultRelease;
	double mechanics = 0.18;
	double drive = 0.32;
	double outputVolume = 0.50;
	double amplifierBass = 0.50;
	double amplifierTreble = 0.50;
	double vibratoSpeed = 0.32;
	double vibratoIntensity = 0.0;
	// There is no historical front-panel equivalent. At 0.5 the graduated
	// factory strike point is unchanged; the control is sampled on key-down.
	double strikePosition = 0.5;
};

// Creative pitch modulation is kept separate from the physical controls.
// Exponential FM, linear FM and phase modulation rotate analytic modal
// coordinates at oversampled pickup readout. V/oct remains the physical
// resonator tuning, so through-zero motion never gives the hammer/contact solve
// a non-physical negative tine frequency. See the modelling notes for scope.
struct ElectricPianoModulation
{
	double exponentialPitch{}; // octaves
	// Additive deviation as a fraction of the positive modal frequency. A value
	// of -1 reaches zero and values below -1 reverse phase direction.
	double linearFrequencyRatio{};
	double phaseRadians{};
};

inline constexpr double ElectricPianoReferenceFrequency =
	261.6255653005986;
inline constexpr std::size_t ElectricPianoAttackModeBegin = 3;
inline constexpr std::size_t ElectricPianoAttackModeEnd = 10;
inline constexpr std::size_t ElectricPianoToneBarSubMode = 10;
inline constexpr std::size_t ElectricPianoModeCount = 11;

// Published quantities are deliberately kept in SI units and separate from
// the small, named calibration trims below.  This prevents a later listening
// pass from silently turning a measured Hunt--Crossley parameter into another
// arbitrary model-space constant.  Values are from Sonderbo (2024), tables
// 3.1/3.2 and section 3.2.1.
struct ElectricPianoPublishedMechanicalData
{
	static constexpr double HammerMassKg = 0.011;
	static constexpr double MaximumHammerVelocityMetresPerSecond = 4.0;
	static constexpr double ContactExponent = 2.8;
	static constexpr double ContactStiffnessNewtonPerMetrePower = 1.5e11;
	static constexpr double ContactDampingWeight = 9.0e10;
	static constexpr double TineRadiusMetres = 0.0005;
	static constexpr double TineDensityKgPerCubicMetre = 7850.0;
	static constexpr double TineYoungsModulusPascal = 2.0e11;
	static constexpr double LongestTineMetres = 0.1564;
	static constexpr double ShortestTineMetres = 0.0226;
	static constexpr double FrequencyIndependentLossPerSecond = 0.0001;
	static constexpr double FrequencyDependentLossSquareMetresPerSecond = 0.005;
	// Sonderbo (2024), table 3.3. These are the damped-spring connection
	// parameters used in the published real-time Rhodes model, not measurements
	// of a particular restored instrument.
	static constexpr double DamperLinearSpringNewtonPerMetre = 100.0;
	static constexpr double DamperNonlinearSpringNewtonPerCubicMetre = 100000.0;
	static constexpr double DamperViscousKgPerSecond = 0.5;
};

inline double ElectricPianoDamperReleaseSeconds(double release)
{
	release = std::clamp(std::isfinite(release) ? release : 0.0, 0.0, 1.0);
	// In the strongly overdamped published connection, its slow relaxation is
	// approximately R/K1 = 5 ms. Use that as the fast, properly adjusted damper
	// endpoint. The logarithmic panel travel extends to the existing 1.2 s
	// deliberate sound-design limit; the factory default retains a small amount
	// of tone-bar transfer rather than imposing the unattainable ideal endpoint.
	constexpr double PublishedDamperRelaxationSeconds =
		ElectricPianoPublishedMechanicalData::DamperViscousKgPerSecond /
		ElectricPianoPublishedMechanicalData::DamperLinearSpringNewtonPerMetre;
	constexpr double MaximumReleaseSeconds = 1.2;
	return PublishedDamperRelaxationSeconds * std::pow(
		MaximumReleaseSeconds / PublishedDamperRelaxationSeconds, release);
}

// These are intentionally the only unsourced voicing degrees of freedom in
// the mechanical calibration.  Keep them close to unity and adjust them only
// from controlled renders/listening tests; the UI controls remain sound-design
// controls rather than compensation for hidden model errors.
struct ElectricPianoMechanicalTrim
{
	static constexpr double HammerVelocity = 1.0;
	// MIDI velocity is not a measurement of hammer speed.  This is the one
	// explicit action-law trim pending optical hammer-velocity measurements.
	static constexpr double HammerVelocityCurveExponent = 0.85;
	static constexpr double ContactStiffness = 1.0;
	static constexpr double TineModalMass = 1.0;
	// A measured pickup flux map is not available. These two independent trims
	// set the excursion-to-field scale and final voltage respectively, so future
	// listening can change magnetic curvature without disturbing gain staging.
	static constexpr double PickupExcursion = 0.05;
	static constexpr double PickupOutput = 0.00033;
	static constexpr double AttackModeOutput = 1.0;
};

struct ElectricPianoCoupledForkProfile
{
	std::array<double, 2> frequencyRatios{};
	std::array<double, 2> toneBarDisplacementRatios{};
	std::array<double, 2> inverseModalMassRatios{};
	std::array<double, 2> supportReactionLossFactors{};
	double toneBarModalMassRatio{};
	double toneBarFreeFrequencyRatio{};
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
	// Service documentation and SLDV both show a strong played fundamental in
	// each prong. The separate sub-fundamental belongs to another tone-bar mode;
	// detuning this fundamental coordinate down to the submode was a structural
	// error that generated a persistent false sideband family in the pickup.
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
	profile.toneBarFreeFrequencyRatio = toneBarFrequencyRatio;
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

inline double ElectricPianoToneBarSubModeCouplingScale(double coupling,
	double keyboardPosition)
{
	// The separate tone-bar submode is transmitted through the same common base
	// as the played fundamental pair. The old reduced coordinate kept its tine
	// participation fixed while Coupling changed the fundamental eigenvectors,
	// which made the control almost exclusively a late-decay adjustment. Use the
	// played tone-bar component as the observable common-base transfer and retain
	// a small intrinsic path for the real fork's distributed mounting. The square
	// root keeps the deliberately wide panel endpoints controlled; this scale is
	// applied reciprocally at hammer projection and pickup observation, so modal
	// residue changes by its square. The calibrated midpoint is exactly unity.
	constexpr double IntrinsicDistributedTransfer = 0.10;
	const auto profile = MakeElectricPianoCoupledForkProfile(coupling,
		keyboardPosition);
	const auto midpoint = MakeElectricPianoCoupledForkProfile(0.50,
		keyboardPosition);
	const double transfer = IntrinsicDistributedTransfer +
		std::abs(profile.toneBarDisplacementRatios[1]);
	const double midpointTransfer = IntrinsicDistributedTransfer +
		std::abs(midpoint.toneBarDisplacementRatios[1]);
	return std::sqrt(transfer / std::max(1.0e-12, midpointTransfer));
}

struct ElectricPianoKeyProfile
{
	double fundamentalHz{};
	double modalMassRatio{};
	double tineLengthMetres{};
	double tineModalMassKg{};
	double pickupSensitivity{};
	double keyboardPosition{};
};

inline double ElectricPianoFactoryHammerTipDurometer(double keyboardPosition)
{
	keyboardPosition = std::clamp(std::isfinite(keyboardPosition) ?
		keyboardPosition : 0.5, 0.0, 1.0);
	const int keyNumber = 1 + static_cast<int>(std::lround(
		72.0 * keyboardPosition));
	if (keyNumber <= 30)
		return 30.0;
	if (keyNumber <= 40)
		return 50.0;
	if (keyNumber <= 50)
		return 70.0;
	if (keyNumber <= 64)
		return 90.0;
	// The service manual calls keys 65--88 wrapped "extra hard" rather than
	// assigning another durometer.  100 is only a convenient normalized tag;
	// the branch itself, not that number, is the sourced classification.
	return 100.0;
}

inline ElectricPianoKeyProfile MakeElectricPianoKeyProfile(double pitchVolts)
{
	const double boundedPitch = std::clamp(
		std::isfinite(pitchVolts) ? pitchVolts : 0.0, -6.0, 6.0);
	const double frequencyRatio = std::exp2(boundedPitch);
	// The published 73-key cutting-chart endpoints are interpolated
	// logarithmically.  A tip-normalized uniform cantilever has an effective
	// modal mass close to one quarter of its distributed mass.  The tuning
	// spring perturbs individual modes, but this gives the contact solve an SI
	// mass scale without inventing a register-dependent loudness law.
	// The pickup-sensitivity term represents the per-key pickup adjustment used
	// to equalize the complete hammer/tine collision, rather than only the bare
	// cantilever impedance. Real pickups are individually positioned; the
	// shallow frequency tilt and end-range correction are an initial
	// whole-keyboard level calibration.
	const double midiNote = 60.0 + 12.0 * boundedPitch;
	const double keyboardPosition = std::clamp(
		(midiNote - 28.0) / 72.0, 0.0, 1.0);
	const double tineLength =
		ElectricPianoPublishedMechanicalData::LongestTineMetres * std::pow(
			ElectricPianoPublishedMechanicalData::ShortestTineMetres /
			ElectricPianoPublishedMechanicalData::LongestTineMetres,
			keyboardPosition);
	constexpr double Pi = 3.1415926535897932384626433832795;
	const double tineArea = Pi *
		ElectricPianoPublishedMechanicalData::TineRadiusMetres *
		ElectricPianoPublishedMechanicalData::TineRadiusMetres;
	const double tineModalMass = 0.25 *
		ElectricPianoPublishedMechanicalData::TineDensityKgPerCubicMetre *
		tineArea * tineLength * ElectricPianoMechanicalTrim::TineModalMass;
	const double referenceKeyboardPosition = (60.0 - 28.0) / 72.0;
	const double referenceLength =
		ElectricPianoPublishedMechanicalData::LongestTineMetres * std::pow(
			ElectricPianoPublishedMechanicalData::ShortestTineMetres /
			ElectricPianoPublishedMechanicalData::LongestTineMetres,
			referenceKeyboardPosition);
	const double modalMassRatio = tineLength / referenceLength;
	// Every real pickup is moved individually during voicing.  The service
	// manual supplies the allowable gaps but no per-key voltage curve, so these
	// nineteen four-key checkpoints remain an intentionally exposed listening
	// fit. They were calibrated from the geometric-mean RMS of soft, medium and
	// hard strikes, not one velocity, with a restrained -0.42 dB/octave voltage
	// tilt. That preserves the contact model's growing treble dynamic range while
	// removing the former 16 dB upper-register hump and abrupt top-key collapse.
	constexpr std::array<double, 19> PickupLevelTrims{
		0.8374, 0.8237, 0.8009, 0.7278, 0.6768, 0.6367, 0.6092,
		0.5759, 0.5479, 0.5117, 0.4838, 0.4601, 0.4519, 0.4187,
		0.4148, 0.4409, 0.4540, 0.5184, 0.6504};
	const double pickupTrimPosition = 18.0 * keyboardPosition;
	const std::size_t pickupTrimLower = std::min<std::size_t>(17,
		static_cast<std::size_t>(pickupTrimPosition));
	const double pickupTrimFraction = pickupTrimPosition -
		static_cast<double>(pickupTrimLower);
	const double pickupLevelTrim = PickupLevelTrims[pickupTrimLower] +
		pickupTrimFraction * (PickupLevelTrims[pickupTrimLower + 1] -
			PickupLevelTrims[pickupTrimLower]);
	return {
		ElectricPianoReferenceFrequency * frequencyRatio,
		modalMassRatio,
		tineLength,
		tineModalMass,
		pickupLevelTrim * std::pow(frequencyRatio, -0.19) *
			(1.0 + std::pow(2.0 * keyboardPosition - 1.0, 2.0)),
		keyboardPosition};
}

inline double ElectricPianoModeBandlimitGain(double frequency,
	double sampleRate)
{
	if (!std::isfinite(frequency) || !std::isfinite(sampleRate) ||
		!(sampleRate > 0.0))
		return 0.0;
	// Do not expose a different mechanical model merely because the host runs
	// above 48 kHz: ultrasonic modes can intermodulate through the pickup and
	// fold into different audible results after a later rate conversion.  The
	// pickup is evaluated at 4x and passes through the seventh-order 4x
	// decimator, so audible attack modes need not be discarded at the former
	// 15.36 kHz threshold before they reach that anti-alias path. Retain them
	// unchanged through 19.2 kHz and taper smoothly to zero below the common
	// 48 kHz Nyquist boundary. The fixed ceiling gives every host rate at or above
	// 48 kHz the same mechanical model; lower rates necessarily use their own
	// Nyquist limit.
	const double modalBandwidthSampleRate = std::min(sampleRate, 48000.0);
	const double normalizedFrequency = std::abs(frequency) /
		modalBandwidthSampleRate;
	const double taper = std::clamp((normalizedFrequency - 0.40) /
		(0.49 - 0.40), 0.0, 1.0);
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
		  exponentialFmInterpolator_(CreateX4Resampler_Cheby7()),
		  linearFmInterpolator_(CreateX4Resampler_Cheby7()),
		  phaseModulationInterpolator_(CreateX4Resampler_Cheby7()),
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
		pickupOutput_ = 0.0;
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
		exponentialFmInterpolator_->Reset();
		linearFmInterpolator_->Reset();
		phaseModulationInterpolator_->Reset();
		verticalPositionInterpolator_->Reset();
		verticalVelocityInterpolator_->Reset();
		horizontalPositionInterpolator_->Reset();
		horizontalVelocityInterpolator_->Reset();
		exponentialFmInterpolatorPrimed_ = false;
		linearFmInterpolatorPrimed_ = false;
		phaseModulationInterpolatorPrimed_ = false;
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
		pickupExcursionScale_ = 1.0;
		pickupOutput_ = 0.0;
		hammerNoise_ = 0.0;
		keyReleaseNoise_ = 0.0;
		damperNoise_ = 0.0;
		mechanicsLowPass_ = 0.0;
		mechanicsHighPass_ = 0.0;
		for (auto& resonator : mechanicalResonators_)
			resonator = {};
		hammerPosition_ = 0.0;
		hammerVelocity_ = 0.0;
		hammerIncomingVelocity_ = MinimumHammerContactVelocity;
		hammerInverseMass_ = 1.0 /
			ElectricPianoPublishedMechanicalData::HammerMassKg;
		contactStiffness_ = 0.0;
		contactExponent_ =
			ElectricPianoPublishedMechanicalData::ContactExponent;
		contactLoss_ = 0.0;
		contactAge_ = 0.0;
		contactActive_ = false;
		contactEngaged_ = false;
		contactModeShape_.fill(0.0);
		modeActive_.fill(false);
		contactModeActive_.fill(false);
		modeInverseMass_.fill(1.0);
		coupledToneBarRatio_.fill(0.0);
		tineModalMass_ = 1.0e-4;
		coupledForkInitialized_ = false;
		energy_ = 0.0;
		modulationPhasorReal_.fill(1.0);
		modulationPhasorImaginary_.fill(0.0);
		modulationPhasorNormalizationCountdown_ = 4096;
		previousOversampledPhaseModulation_ = 0.0;
		previousHostExponentialPitch_ = 0.0;
		previousHostLinearFrequencyRatio_ = 0.0;
		previousHostPhaseModulation_ = 0.0;
		modulationPathActive_ = false;
		modulationControlsPrimed_ = false;
		exponentialFmPathActive_ = false;
		linearFmPathActive_ = false;
		phaseModulationPathActive_ = false;
		exponentialFmInterpolatorPrimed_ = false;
		linearFmInterpolatorPrimed_ = false;
		phaseModulationInterpolatorPrimed_ = false;
		modulationCrossfade_ = 0.0;
		for (auto& frame : pickupModeFrames_)
			frame.fill({});
		controlsInitialized_ = false;
		coefficientsDirty_ = true;
		timbreDirty_ = true;
		exponentialFmInterpolator_->Reset();
		linearFmInterpolator_->Reset();
		phaseModulationInterpolator_->Reset();
		verticalPositionInterpolator_->Reset();
		verticalVelocityInterpolator_->Reset();
		horizontalPositionInterpolator_->Reset();
		horizontalVelocityInterpolator_->Reset();
		pickupDecimator_->Reset();
	}

	double Step(double pitchVolts, double gateVolts, double velocity,
		bool sustain, const ElectricPianoControls& controls,
		const ElectricPianoModulation& modulation = {}, bool retrigger = false)
	{
		if (!std::isfinite(pitchVolts))
			pitchVolts = 0.0;
		if (!std::isfinite(gateVolts))
			gateVolts = 0.0;
		velocity = std::clamp(std::isfinite(velocity) ? velocity : 0.8,
			0.0, 1.0);
		const double exponentialPitch = std::clamp(
			std::isfinite(modulation.exponentialPitch) ?
				modulation.exponentialPitch : 0.0, -4.0, 4.0);
		const double linearFrequencyRatio = std::clamp(
			std::isfinite(modulation.linearFrequencyRatio) ?
				modulation.linearFrequencyRatio : 0.0, -4.0, 4.0);
		const double phaseModulation = std::clamp(
			std::isfinite(modulation.phaseRadians) ?
				modulation.phaseRadians : 0.0, -2.0 * TwoPi, 2.0 * TwoPi);

		const bool gate = gateVolts >= 1.0;
		bool newSilentStrike = false;
		if (gate && (!lastGate_ || retrigger))
		{
			const bool wasAudible = IsAudible();
			if (!wasAudible)
			{
				newSilentStrike = true;
				// A dormant voice is not stepped, so its smoothed controls otherwise
				// remain frozen at the preceding note's values. Re-prime them from the
				// current panel state before this fresh attack; in particular, do not
				// spend the next note's transient gliding old pickup/coupling settings
				// or the slow Decay smoother toward values chosen during silence.
				controlsInitialized_ = false;
				// A silent key has no phase history to preserve. The oversampled
				// modulation controls are primed below, so an already-present DC PM
				// voltage becomes initial phase rather than a one-sample impulse.
				modulationPhasorReal_.fill(1.0);
				modulationPhasorImaginary_.fill(0.0);
				modulationPhasorNormalizationCountdown_ = 4096;
				modulationPathActive_ = false;
				modulationControlsPrimed_ = false;
				exponentialFmPathActive_ = false;
				linearFmPathActive_ = false;
				phaseModulationPathActive_ = false;
				exponentialFmInterpolatorPrimed_ = false;
				linearFmInterpolatorPrimed_ = false;
				phaseModulationInterpolatorPrimed_ = false;
				modulationCrossfade_ = 0.0;
				verticalPositionInterpolator_->Reset();
				verticalVelocityInterpolator_->Reset();
				horizontalPositionInterpolator_->Reset();
				horizontalVelocityInterpolator_->Reset();
				previousHostExponentialPitch_ = exponentialPitch;
				previousHostLinearFrequencyRatio_ = linearFrequencyRatio;
				previousHostPhaseModulation_ = phaseModulation;
			}
			Strike(pitchVolts, velocity, controls);
		}
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
		RefreshTimbreCoefficients(activeControls,
			coefficientUpdateTick || pitchChanged);

		if (exponentialPitch != 0.0)
			exponentialFmPathActive_ = true;
		if (linearFrequencyRatio != 0.0)
			linearFmPathActive_ = true;
		if (phaseModulation != 0.0)
			phaseModulationPathActive_ = true;
		modulationPathActive_ = exponentialFmPathActive_ ||
			linearFmPathActive_ || phaseModulationPathActive_;
		if (newSilentStrike && modulationPathActive_)
			modulationCrossfade_ = 1.0;

		if (contactActive_)
			AdvanceCoupledHammerAndModes();
		else
			AdvanceFreeModes();

		Eigen::Array<double, PickupOversamplingFactor, 1> exponentialPitchFrames;
		Eigen::Array<double, PickupOversamplingFactor, 1> linearRatioFrames;
		Eigen::Array<double, PickupOversamplingFactor, 1> phaseModulationFrames;
		if (modulationPathActive_)
		{
			exponentialPitchFrames.setZero();
			linearRatioFrames.setZero();
			phaseModulationFrames.setZero();
			if (!modulationControlsPrimed_)
			{
				previousOversampledPhaseModulation_ =
					previousHostPhaseModulation_;
				for (std::size_t mode = 0; mode < modes_.size(); ++mode)
				{
					const double initialPhase = ModulationModeWeights[mode] *
						previousOversampledPhaseModulation_;
					modulationPhasorReal_[mode] = std::cos(initialPhase);
					modulationPhasorImaginary_[mode] = std::sin(initialPhase);
				}
				modulationControlsPrimed_ = true;
			}
			if (exponentialFmPathActive_)
			{
				if (!exponentialFmInterpolatorPrimed_)
				{
					exponentialFmInterpolator_->PrimeUpsample(
						previousHostExponentialPitch_);
					exponentialFmInterpolatorPrimed_ = true;
				}
				exponentialPitchFrames = exponentialFmInterpolator_->Upsample(
					exponentialPitch);
			}
			if (linearFmPathActive_)
			{
				if (!linearFmInterpolatorPrimed_)
				{
					linearFmInterpolator_->PrimeUpsample(
						previousHostLinearFrequencyRatio_);
					linearFmInterpolatorPrimed_ = true;
				}
				linearRatioFrames = linearFmInterpolator_->Upsample(
					linearFrequencyRatio);
			}
			if (phaseModulationPathActive_)
			{
				if (!phaseModulationInterpolatorPrimed_)
				{
					phaseModulationInterpolator_->PrimeUpsample(
						previousHostPhaseModulation_);
					phaseModulationInterpolatorPrimed_ = true;
				}
				phaseModulationFrames = phaseModulationInterpolator_->Upsample(
					phaseModulation);
			}
		}
		previousHostExponentialPitch_ = exponentialPitch;
		previousHostLinearFrequencyRatio_ = linearFrequencyRatio;
		previousHostPhaseModulation_ = phaseModulation;

		energy_ = 0.0;
		for (std::size_t index = 0; index < modes_.size(); ++index)
		{
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

		// Preserve the established, cheaper aggregate interpolation path whenever
		// FM/PM has never been active on this note. If modulation arrives on a held
		// note, keep it alive briefly and crossfade to the per-mode 4x readout; this
		// avoids both a timing discontinuity and permanent baseline CPU overhead.
		const bool legacyPickupNeeded = !modulationPathActive_ ||
			modulationCrossfade_ < 1.0;
		Eigen::Array<double, PickupOversamplingFactor, 1> legacyVerticalPositions;
		Eigen::Array<double, PickupOversamplingFactor, 1> legacyVerticalVelocities;
		Eigen::Array<double, PickupOversamplingFactor, 1> legacyHorizontalPositions;
		Eigen::Array<double, PickupOversamplingFactor, 1> legacyHorizontalVelocities;
		if (legacyPickupNeeded)
		{
			double verticalPosition = 0.0;
			double verticalVelocity = 0.0;
			double horizontalPosition = 0.0;
			double horizontalVelocity = 0.0;
			for (std::size_t mode = 0; mode < modes_.size(); ++mode)
			{
				if (!modeActive_[mode])
					continue;
				const double verticalWeight = modeOutputWeight_[mode] *
					modeBandlimitGain_[mode];
				const double horizontalWeight = modeHorizontalWeight_[mode] *
					modeBandlimitGain_[mode];
				const double velocityScale = TineDisplacementScale *
					pickupExcursionScale_ *
					modeAngularFrequency_[mode];
				verticalPosition += verticalWeight * modes_[mode].real;
				verticalVelocity -= verticalWeight * velocityScale *
					modes_[mode].imaginary;
				horizontalPosition += horizontalWeight * modes_[mode].real;
				horizontalVelocity -= horizontalWeight * velocityScale *
					modes_[mode].imaginary;
			}
			legacyVerticalPositions = verticalPositionInterpolator_->Upsample(
				TineDisplacementScale * pickupExcursionScale_ * verticalPosition);
			legacyVerticalVelocities = verticalVelocityInterpolator_->Upsample(
				verticalVelocity);
			legacyHorizontalPositions = horizontalPositionInterpolator_->Upsample(
				TineDisplacementScale * pickupExcursionScale_ * horizontalPosition);
			legacyHorizontalVelocities = horizontalVelocityInterpolator_->Upsample(
				horizontalVelocity);
		}
		Eigen::Array<double, PickupOversamplingFactor, 1> pickupValues;
		const double pickupRate = sampleRate_ *
			static_cast<double>(PickupOversamplingFactor);
		for (int frame = 0; frame < PickupOversamplingFactor; ++frame)
		{
			if (!modulationPathActive_)
			{
				const auto legacyGradient = MagneticFluxGradient(
					pickupVerticalOffset_ + legacyVerticalPositions(frame),
					pickupHorizontalOffset_ + legacyHorizontalPositions(frame),
					pickupGap_);
				const double legacyEmf = legacyGradient[0] *
					legacyVerticalVelocities(frame) + legacyGradient[1] *
					legacyHorizontalVelocities(frame);
				pickupValues(frame) = std::clamp(
					pickupVoltageScale_ * legacyEmf, -12.0, 12.0) +
					mechanicalSignal;
				continue;
			}
			double verticalPosition = 0.0;
			double verticalVelocity = 0.0;
			double horizontalPosition = 0.0;
			double horizontalVelocity = 0.0;
			const double frameExponentialPitch = std::clamp(
				exponentialPitchFrames(frame), -4.0, 4.0);
			const double frameLinearRatio = std::clamp(linearRatioFrames(frame),
				-4.0, 4.0);
			const double framePhaseModulation = std::clamp(
				phaseModulationFrames(frame), -2.0 * TwoPi, 2.0 * TwoPi);
			const double bodyExponentialRatio =
				frameExponentialPitch == 0.0 ? 1.0 :
				Exp2Taylor5(static_cast<float>(frameExponentialPitch));
			const double phaseModulationFrequency = modulationPathActive_ ?
				(framePhaseModulation - previousOversampledPhaseModulation_) *
					pickupRate / TwoPi : 0.0;
			for (std::size_t mode = 0;
				mode < ElectricPianoToneBarSubMode; ++mode)
			{
				if (!modeActive_[mode])
					continue;
				const double physicalFrequency = currentFundamental_ *
					modeRatio_[mode];
				double renderedFrequency = physicalFrequency;
				double renderedReal = pickupModeFrames_[frame][mode].real;
				double renderedImaginary =
					pickupModeFrames_[frame][mode].imaginary;
				double renderedBandlimit = modeBandlimitGain_[mode];
				if (modulationPathActive_ && ModulationModeWeights[mode] != 0.0)
				{
					const double modulationWeight = ModulationModeWeights[mode];
					const double weightedExponentialPitch = std::clamp(
						modulationWeight * frameExponentialPitch, -4.0, 4.0);
					const double exponentialRatio = modulationWeight == 1.0 ?
						bodyExponentialRatio :
						(weightedExponentialPitch == 0.0 ? 1.0 :
							Exp2Taylor5(static_cast<float>(weightedExponentialPitch)));
					const double frequencyOffset = physicalFrequency *
						(exponentialRatio - 1.0 +
							modulationWeight * frameLinearRatio);
					const double phaseIncrement = TwoPi * frequencyOffset / pickupRate +
						modulationWeight * (framePhaseModulation -
							previousOversampledPhaseModulation_);
					const auto rotation = ModulationPhaseRotation(phaseIncrement);
					const double oldPhasorReal = modulationPhasorReal_[mode];
					const double oldPhasorImaginary =
						modulationPhasorImaginary_[mode];
					modulationPhasorReal_[mode] = rotation[0] * oldPhasorReal -
						rotation[1] * oldPhasorImaginary;
					modulationPhasorImaginary_[mode] = rotation[1] * oldPhasorReal +
						rotation[0] * oldPhasorImaginary;
					const double phaseCosine = modulationPhasorReal_[mode];
					const double phaseSine = modulationPhasorImaginary_[mode];
					const double physicalReal = renderedReal;
					const double physicalImaginary = renderedImaginary;
					renderedReal = phaseCosine * physicalReal -
						phaseSine * physicalImaginary;
					renderedImaginary = phaseSine * physicalReal +
						phaseCosine * physicalImaginary;
					renderedFrequency = physicalFrequency + frequencyOffset +
						modulationWeight * phaseModulationFrequency;
					renderedBandlimit = std::min(renderedBandlimit,
						ElectricPianoModeBandlimitGain(
							std::abs(renderedFrequency), sampleRate_));
				}
				const double verticalWeight = modeOutputWeight_[mode] *
					renderedBandlimit;
				const double horizontalWeight = modeHorizontalWeight_[mode] *
					renderedBandlimit;
				const double velocityScale = TineDisplacementScale *
					pickupExcursionScale_ * TwoPi * renderedFrequency;
				verticalPosition += verticalWeight * renderedReal;
				verticalVelocity -= verticalWeight * velocityScale *
					renderedImaginary;
				horizontalPosition += horizontalWeight * renderedReal;
				horizontalVelocity -= horizontalWeight * velocityScale *
					renderedImaginary;
			}
			// The separately measured tone-bar submode remains physical under the
			// creative FM/PM inputs. Keep it out of the hot ten-coordinate modulation
			// loop: besides avoiding meaningless phasor work, this preserves the
			// compiler's established vectorization for polyphonic audio-rate FM.
			constexpr std::size_t subMode = ElectricPianoToneBarSubMode;
			if (modeActive_[subMode])
			{
				const double renderedBandlimit = modeBandlimitGain_[subMode];
				const double verticalWeight = modeOutputWeight_[subMode] *
					renderedBandlimit;
				const double horizontalWeight = modeHorizontalWeight_[subMode] *
					renderedBandlimit;
				const double renderedFrequency = currentFundamental_ *
					modeRatio_[subMode];
				const double velocityScale = TineDisplacementScale *
					pickupExcursionScale_ * TwoPi * renderedFrequency;
				const double renderedReal =
					pickupModeFrames_[frame][subMode].real;
				const double renderedImaginary =
					pickupModeFrames_[frame][subMode].imaginary;
				verticalPosition += verticalWeight * renderedReal;
				verticalVelocity -= verticalWeight * velocityScale *
					renderedImaginary;
				horizontalPosition += horizontalWeight * renderedReal;
				horizontalVelocity -= horizontalWeight * velocityScale *
					renderedImaginary;
			}
			if (modulationPathActive_)
				previousOversampledPhaseModulation_ = framePhaseModulation;
			double modulatedPickup = 0.0;
			if (modulationPathActive_)
			{
				const auto gradient = MagneticFluxGradient(
					pickupVerticalOffset_ + TineDisplacementScale *
						pickupExcursionScale_ * verticalPosition,
					pickupHorizontalOffset_ + TineDisplacementScale *
						pickupExcursionScale_ * horizontalPosition,
					pickupGap_);
				const double emf = gradient[0] * verticalVelocity +
					gradient[1] * horizontalVelocity;
				modulatedPickup = pickupVoltageScale_ * emf;
			}
			double selectedPickup = modulatedPickup;
			if (legacyPickupNeeded)
			{
				const auto legacyGradient = MagneticFluxGradient(
					pickupVerticalOffset_ + legacyVerticalPositions(frame),
					pickupHorizontalOffset_ + legacyHorizontalPositions(frame),
					pickupGap_);
				const double legacyEmf = legacyGradient[0] *
					legacyVerticalVelocities(frame) + legacyGradient[1] *
					legacyHorizontalVelocities(frame);
				const double legacyPickup = pickupVoltageScale_ * legacyEmf;
				selectedPickup = modulationPathActive_ ?
					legacyPickup + modulationCrossfade_ *
						(modulatedPickup - legacyPickup) : legacyPickup;
			}
			pickupValues(frame) = std::clamp(selectedPickup, -12.0, 12.0) +
				mechanicalSignal;
		}
		if (modulationPathActive_ && modulationCrossfade_ < 1.0)
			modulationCrossfade_ = std::min(1.0, modulationCrossfade_ +
				1.0 / std::max(1.0, 0.006 * sampleRate_));
		if (modulationPathActive_ &&
			--modulationPhasorNormalizationCountdown_ <= 0)
		{
			for (std::size_t mode = 0; mode < modes_.size(); ++mode)
			{
				const double magnitude = std::hypot(modulationPhasorReal_[mode],
					modulationPhasorImaginary_[mode]);
				if (magnitude > 0.0)
				{
					modulationPhasorReal_[mode] /= magnitude;
					modulationPhasorImaginary_[mode] /= magnitude;
				}
			}
			modulationPhasorNormalizationCountdown_ = 4096;
		}
		// The seventh-order 4x decimator is the pickup reconstruction filter.
		// A former host-rate 16.5 kHz one-pole was not part of the pickup or
		// Peterson circuit and added 0.6--1.9 dB of avoidable upper-band loss.
		// Keep the decimated signal directly; alias rejection is verified at the
		// actual 4x boundary by the register-wide spectral tests.
		pickupOutput_ = pickupDecimator_->Downsample(pickupValues);
		hammerNoise_ *= hammerNoiseDecay_;
		keyReleaseNoise_ *= keyReleaseNoiseDecay_;
		damperNoise_ *= damperNoiseDecay_;

		if (energy_ < 1.0e-14 && hammerNoise_ < 1.0e-8 &&
			keyReleaseNoise_ < 1.0e-8 && damperNoise_ < 1.0e-8)
		{
			for (auto& mode : modes_)
				mode = {};
		}

		const double result = pickupOutput_;
		return std::isfinite(result) ? result : 0.0;
	}

	double Energy() const { return energy_; }
	double Activity() const
	{
		return energy_ + hammerNoise_ * hammerNoise_ +
			keyReleaseNoise_ * keyReleaseNoise_ +
			damperNoise_ * damperNoise_ +
			pickupOutput_ * pickupOutput_;
	}
	bool IsAudible() const { return Activity() > 1.0e-12; }
	bool GateHigh() const { return lastGate_; }
	double NotePitch() const { return notePitch_; }
	bool ContactActive() const { return contactActive_; }
	double ContactAge() const { return contactAge_; }
	double ModeFrequencyRatio(std::size_t index) const
	{
		return index < modeRatio_.size() ? modeRatio_[index] : 0.0;
	}
	double ModeAmplitudeLifetimeSeconds(std::size_t index) const
	{
		if (index >= modeRadius_.size() || modeRadius_[index] <= 0.0 ||
			modeRadius_[index] >= 1.0)
			return 0.0;
		return -1.0 / (sampleRate_ * std::log(modeRadius_[index]));
	}
	double ModeDisplacementAmplitude(std::size_t index) const
	{
		return index < modes_.size() ?
			std::hypot(modes_[index].real, modes_[index].imaginary) : 0.0;
	}
	double ModePickupDisplacementAmplitude(std::size_t index) const
	{
		return index < modes_.size() ? ModeDisplacementAmplitude(index) *
			std::abs(modeOutputWeight_[index]) : 0.0;
	}
	bool ModeRendered(std::size_t index) const
	{
		return index < modeActive_.size() && modeActive_[index];
	}
	bool ModeParticipatesInContact(std::size_t index) const
	{
		return index < contactModeActive_.size() && contactModeActive_[index];
	}
	double ContactModeProjection(std::size_t index) const
	{
		return index < contactModeShape_.size() ? contactModeShape_[index] : 0.0;
	}
	double StrikePosition() const { return currentStrikePosition_; }
	double ContactWidthMetres() const { return contactWidthMetres_; }
	static std::array<double, 2> MagneticPickupGradient(double vertical,
		double horizontal, double gap)
	{
		return MagneticFluxGradient(vertical, horizontal, gap);
	}
	static std::array<double, 2> DirectMagneticPickupGradient(
		double vertical, double horizontal, double gap)
	{
		constexpr double HorizontalFieldScale = 0.62;
		const double radial = std::sqrt(gap * gap +
			HorizontalFieldScale * horizontal * horizontal + 0.020);
		const auto field = EvaluateCalibratedPoleGradient(vertical, radial);
		return {field.vertical, field.radial * HorizontalFieldScale * horizontal /
			std::max(1.0e-9, radial)};
	}

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
			// Only these continuously active controls are consumed below. Sanitize
			// them while priming so a non-finite API caller cannot poison the
			// smoother's persistent state before Clamp01() sees the target.
			smoothedControls_.body = Clamp01(controls.body);
			smoothedControls_.bell = Clamp01(controls.bell);
			smoothedControls_.coupling = Clamp01(controls.coupling);
			smoothedControls_.tone = Clamp01(controls.tone);
			smoothedControls_.proximity = Clamp01(controls.proximity);
			smoothedControls_.decay = Clamp01(controls.decay);
			smoothedControls_.release = Clamp01(controls.release);
			smoothedControls_.mechanics = Clamp01(controls.mechanics);
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
		const double bassExcursionPosition = std::clamp(keyPosition_ /
			((60.0 - 28.0) / 72.0), 0.0, 1.0);
		const double bassExcursionSmooth = bassExcursionPosition *
			bassExcursionPosition * (3.0 - 2.0 * bassExcursionPosition);
		// The published gap range tops out at 1/8 inch, yet the reduced uniform-
		// beam coordinates over-predict the excursion of the longest, tapered
		// tines. Keep this visible as a bass-only calibration awaiting measured
		// tine-tip trajectories; it reaches unity at C4 and never boosts motion.
		pickupExcursionScale_ = 0.55 + 0.45 * bassExcursionSmooth;
		coefficientsDirty_ = true;
		timbreDirty_ = true;
		const double velocityCurve = Clamp01(controls.velocityCurve);
		const double dynamics = Clamp01(controls.dynamics);
		const double hammer = Clamp01(controls.hammer);
		const double strikePosition = Clamp01(controls.strikePosition);
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
		const double factoryDurometer =
			ElectricPianoFactoryHammerTipDurometer(keyPosition_);
		const double factoryHardness = (factoryDurometer - 30.0) / 70.0;
		const double contactHardness = std::clamp(factoryHardness +
			0.70 * (hammer - 0.52), 0.0, 1.0);
		hammerInverseMass_ = 1.0 /
			ElectricPianoPublishedMechanicalData::HammerMassKg;
		hammerPosition_ = 0.0;
		const double requestedHammerVelocity =
			ElectricPianoPublishedMechanicalData::
				MaximumHammerVelocityMetresPerSecond *
			ElectricPianoMechanicalTrim::HammerVelocity *
			std::pow(dynamicAmplitude,
				ElectricPianoMechanicalTrim::HammerVelocityCurveExponent);
		// Strikes below this physical speed are many orders below the audible
		// output range but can otherwise spend seconds in an almost static
		// oversampled collision when driven by tiny positive CV residue. Classify
		// them as zero-energy before contact begins; never truncate an active
		// collision based on elapsed time.
		const bool audibleHammerStrike =
			requestedHammerVelocity >= MinimumHammerContactVelocity;
		hammerVelocity_ = audibleHammerStrike ? requestedHammerVelocity : 0.0;
		hammerIncomingVelocity_ = std::max(MinimumHammerContactVelocity,
			requestedHammerVelocity);
		// The service manual documents five graduated tip grades (30, 50, 70,
		// 90 durometer, then wrapped extra-hard) but not their force curves.  The
		// published k is therefore the middle-C/default anchor.  The modest zoned
		// slope and the symmetric panel span are explicit listening trims around
		// that anchor, not substitutes for the sourced SI stiffness.
		const double graduatedTipOctaves = 1.5 * std::log2(
			factoryDurometer / 50.0);
		const double panelHardnessOctaves = 5.0 * (hammer - 0.52);
		contactStiffness_ =
			ElectricPianoPublishedMechanicalData::
				ContactStiffnessNewtonPerMetrePower *
			ElectricPianoMechanicalTrim::ContactStiffness *
			std::exp2(graduatedTipOctaves + panelHardnessOctaves);
		contactExponent_ =
			ElectricPianoPublishedMechanicalData::ContactExponent;
		// Sonderbo reports lambda separately from k. Scale lambda with the same
		// tip-grade factor so lambda/k (seconds/metre) remains the measured 0.6
		// while Hammer changes compliance rather than inventing restitution.
		contactLoss_ =
			ElectricPianoPublishedMechanicalData::ContactDampingWeight /
			ElectricPianoPublishedMechanicalData::
				ContactStiffnessNewtonPerMetrePower;
		contactAge_ = 0.0;
		contactActive_ = audibleHammerStrike;
		contactEngaged_ = false;
		currentStrikePosition_ = StrikePositionFromControl(keyPosition_,
			strikePosition);
		// Falaize's distributed-force reference uses a 15 mm hammer width. The
		// effective loaded strip is smaller than the complete tip, so retain the
		// calibrated 6--12 mm hard-to-soft interval. Deriving this longitudinal
		// width from Hunt--Crossley indentation assumed a spherical 4 mm tip that
		// is neither present in the source nor representative of the block-shaped
		// Rhodes neoprene tip; it narrowed the footprint to sub-5.5 mm values and
		// changed the previously validated bass modal balance.
		contactWidthMetres_ = 0.006 + 0.006 * (1.0 - contactHardness);
		const double fundamentalProjection = FiniteContactModeProjection(0,
			keyPosition_, currentStrikePosition_, contactWidthMetres_);
		contactModeShape_[0] = fundamentalProjection;
		contactModeShape_[1] = fundamentalProjection;
		contactModeShape_[2] = -0.055 * fundamentalProjection;
		tineModalMass_ = key.tineModalMassKg;
		modeInverseMass_[2] = 1.0 / (key.tineModalMassKg * 1.65);
		for (std::size_t index = ElectricPianoAttackModeBegin;
			index < ElectricPianoAttackModeEnd; ++index)
		{
			// Generalized force projection is a property of contact geometry and the
			// resonator mode, never of MIDI velocity.  Integrate each mode over a
			// finite neoprene patch; contact duration then supplies velocity-dependent
			// bandwidth through the Hunt-Crossley force itself.
			contactModeShape_[index] = FiniteContactModeProjection(index,
				keyPosition_, currentStrikePosition_, contactWidthMetres_) *
				AttackModeEnergyCalibration(index, keyPosition_);
			modeInverseMass_[index] = 1.0 /
				(key.tineModalMassKg *
					AttackModeModalMassMultiplier(index, keyPosition_));
		}
		// The tone bar has a separately observed sub-fundamental mode in addition
		// to its strong played fundamental. In normal-coordinate form the tine
		// participation projects both the hammer force and pickup observation;
		// applying it reciprocally avoids an arbitrary output-only oscillator.
		const double subModeParticipation =
			ToneBarSubModeTineParticipation(keyPosition_) *
			ElectricPianoToneBarSubModeCouplingScale(controls.coupling,
				keyPosition_);
		contactModeShape_[ElectricPianoToneBarSubMode] =
			subModeParticipation * fundamentalProjection;
		modeInverseMass_[ElectricPianoToneBarSubMode] =
			1.0 / key.tineModalMassKg;

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

	void RefreshTimbreCoefficients(const ElectricPianoControls& controls,
		bool controlTick)
	{
		// The nonlinear trajectory normalizer below is deliberately control-rate.
		// All inputs have already been smoothed at audio rate, and 1 kHz coefficient
		// refresh is far above knob/CV bandwidth; evaluating dozens of field points
		// for every active voice on every host sample would make this calibration a
		// significant steady CPU regression while a control is moving.
		if (!timbreDirty_ && !controlTick)
			return;
		const double body = Clamp01(controls.body);
		const double bell = Clamp01(controls.bell);
		const double coupling = Clamp01(controls.coupling);
		const double proximity = Clamp01(controls.proximity);
		const double tone = Clamp01(controls.tone);
		const double mechanics = Clamp01(controls.mechanics);
		if (!timbreDirty_ && body == cachedBodyWeight_ &&
			bell == cachedBellWeight_ &&
			coupling == cachedCouplingWeight_ &&
			proximity == cachedProximity_ &&
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
		// Bell is a deliberate balance control around the calibrated physical
		// residue, not a correction for the default model. The upper keyboard has
		// progressively fewer attack coordinates below Nyquist, so normalize panel
		// sensitivity from a 24 dB bass/middle span toward 36 dB at C6 rather than
		// leaving the same knob effectively dead there. The midpoint remains unity.
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		constexpr double MiddleCPosition = (60.0 - 28.0) / 72.0;
		constexpr double C6Position = (84.0 - 28.0) / 72.0;
		const double upperBellLinear = std::clamp(
			(keyPosition_ - F3Position) / (C6Position - F3Position), 0.0, 1.0);
		const double upperBellSmooth = upperBellLinear * upperBellLinear *
			(3.0 - 2.0 * upperBellLinear);
		const double bellOctaveSpan = 4.0 + 2.0 * upperBellSmooth;
		const double bellGain = std::exp2(bellOctaveSpan * (bell - 0.52));
		for (std::size_t index = ElectricPianoAttackModeBegin;
			index < ElectricPianoAttackModeEnd; ++index)
			modeOutputWeight_[index] =
				ElectricPianoMechanicalTrim::AttackModeOutput *
				bellGain;

		// Imperfect transverse coupling gives the tine tip a shallow elliptical
		// orbit. The near-fundamental third mode carries most of the horizontal
		// motion; attack modes alternate orientation as observed on a real fork.
		modeHorizontalWeight_[0] = 0.010;
		modeHorizontalWeight_[1] = -0.012;
		modeHorizontalWeight_[2] = 0.42 + 0.20 * body;
		modeOutputWeight_[ElectricPianoToneBarSubMode] =
			modeOutputWeight_[1] * ToneBarSubModeTineParticipation(keyPosition_) *
			ElectricPianoToneBarSubModeCouplingScale(coupling,
				keyPosition_);
		for (std::size_t index = ElectricPianoAttackModeBegin;
			index < ElectricPianoAttackModeEnd; ++index)
			modeHorizontalWeight_[index] = (index % 2 == 0 ? -1.0 : 1.0) *
			0.097 * bellGain;
		modeHorizontalWeight_[ElectricPianoToneBarSubMode] = 0.0;

		// The service adjustment called "timbre" is the tine's vertical alignment
		// to the pickup pole. Moving toward the pole centre suppresses the linear
		// fundamental relative to curvature-generated harmonics. PROXIMITY controls
		// the independent front-to-back distance and therefore magnetic curvature.
		// Each real key is voiced by moving its pickup. Trajectory normalization
		// below makes this a dynamic-timbre calibration, not a register-dependent
		// gain contour.
		// The service range is expressed directly in millimetres: ordinary
		// voicing is 1/16--1/8 inch, while post-1972 middle/upper pickups can be
		// brought as close as 0.020 inch. The panel midpoint is the keyed neutral
		// setup below; each half is logarithmic so the close end retains useful
		// resolution.
		constexpr double MinimumServiceGapMillimetres = 0.020 * 25.4;
		constexpr double DefaultServiceGapMillimetres = 1.6;
		constexpr double MaximumServiceGapMillimetres = 0.125 * 25.4;
		// A serviced harp is individually voiced. Long bass tines need the ordinary
		// 1/16-inch clearance, while the service manual explicitly permits 0.020 inch
		// in the middle/upper range to retain dynamic response from their much smaller
		// travel. Preserve the successful bass and bark calibration exactly through
		// middle C, then graduate the neutral per-key screw setting toward 0.52 mm at
		// the top, still just above the documented 0.020-inch service limit. The
		// former 0.60 mm endpoint left the highest keys with little pickup-curvature
		// sparkle even after the mechanical topology was corrected. This remains a
		// physical gap change, not a register EQ or output crossfade.
		constexpr double TrebleNeutralGapMillimetres = 0.52;
		const double trebleGapPosition = std::clamp(
			(keyPosition_ - MiddleCPosition) / (1.0 - MiddleCPosition), 0.0, 1.0);
		const double trebleGapSmooth = trebleGapPosition * trebleGapPosition *
			(3.0 - 2.0 * trebleGapPosition);
		const double neutralKeyGap = DefaultServiceGapMillimetres * std::pow(
			TrebleNeutralGapMillimetres / DefaultServiceGapMillimetres,
			trebleGapSmooth);
		// Only the closest endpoint is additionally zoned: bass keys retain the
		// ordinary 1/16-inch service limit.
		const double closeGapTransition = std::clamp(keyPosition_ /
			((52.0 - 28.0) / 72.0), 0.0, 1.0);
		const double closeGapSmooth = closeGapTransition * closeGapTransition *
			(3.0 - 2.0 * closeGapTransition);
		const double minimumKeyGap = 0.0625 * 25.4 + closeGapSmooth *
			(MinimumServiceGapMillimetres - 0.0625 * 25.4);
		pickupGap_ = proximity < 0.48 ?
			MaximumServiceGapMillimetres * std::pow(
				neutralKeyGap / MaximumServiceGapMillimetres,
				proximity / 0.48) :
			neutralKeyGap * std::pow(
				minimumKeyGap / neutralKeyGap,
				(proximity - 0.48) / 0.52);
		// At minimum Tone the tine rests close to the primary edge's maximum
		// gradient. Raising Tone moves it onto one flank, trading odd-harmonic
		// symmetry for the stronger even/sideband mixture used when voicing a
		// Rhodes for bite.
		// In the treble, most inharmonic attack coordinates lie above pickup
		// bandwidth and raw modal Bell gain cannot create a useful control. Real
		// Rhodes voicing also changes bell character by moving the pickup vertically
		// relative to the tine. Add that second physical route only above middle C,
		// where it smoothly takes over from the mechanical residues. The signed
		// offset is exactly zero at Bell's calibrated 0.52 default. Trajectory and
		// alignment normalization below remove the associated broadband level shift,
		// leaving curvature-generated harmonic balance rather than another VCA.
		const double pickupBellLinear = std::clamp(
			(keyPosition_ - MiddleCPosition) /
				(C6Position - MiddleCPosition), 0.0, 1.0);
		const double pickupBellBlend = pickupBellLinear * pickupBellLinear *
			(3.0 - 2.0 * pickupBellLinear);
		constexpr double PickupBellAlignmentSpanMillimetres = 0.09;
		const double neutralPickupVerticalOffset = 0.34 + 0.22 * tone * tone +
			0.020 * keyPosition_;
		const double minimumPickupVerticalOffset = 0.34 + 0.020 * keyPosition_;
		const double maximumPickupVerticalOffset = 0.56 + 0.020 * keyPosition_;
		pickupVerticalOffset_ = std::clamp(neutralPickupVerticalOffset +
			PickupBellAlignmentSpanMillimetres * pickupBellBlend * (bell - 0.52),
			minimumPickupVerticalOffset, maximumPickupVerticalOffset);
		pickupHorizontalOffset_ = 0.10 + 0.035 * keyPosition_;
		const auto alignmentGradient = MagneticFluxGradient(
			neutralPickupVerticalOffset, pickupHorizontalOffset_,
			ReferencePickupGap);
		const double alignmentMagnitude = std::sqrt(
			alignmentGradient[0] * alignmentGradient[0] +
			alignmentGradient[1] * alignmentGradient[1]);
		// Small-signal normalization alone cannot calibrate a nonlinear pickup: it
		// made the same panel movement cut bass keys while boosting upper keys by
		// roughly 6 dB. Normalize the RMS EMF of a representative keyed trajectory
		// instead. The nominal excursion follows the existing, documented bass
		// excursion trim; sampling three amplitudes makes the correction robust to
		// velocity without imposing an envelope on the actual signal. Curvature and
		// harmonic balance are retained because this is one scalar per pickup setup,
		// never an instantaneous waveshaper inversion.
		const double nominalExcursionMillimetres =
			(0.08 + 0.32 * std::pow(1.0 - keyPosition_, 1.5)) *
			pickupExcursionScale_;
		const double referenceTrajectory = MagneticTrajectoryRms(
			neutralPickupVerticalOffset, pickupHorizontalOffset_,
			ReferencePickupGap, nominalExcursionMillimetres);
		const double currentTrajectory = MagneticTrajectoryRms(
			pickupVerticalOffset_, pickupHorizontalOffset_, pickupGap_,
			nominalExcursionMillimetres);
		// Compare Bell's shifted observation against the same key's neutral pickup
		// placement, rather than allowing the shifted placement to establish a new
		// gain reference. This jointly normalizes Bell and Proximity across the
		// representative trajectories and prevents extreme Tone/Proximity settings
		// from turning their interaction into several decibels of hidden gain.
		const double trajectoryGain = referenceTrajectory /
			std::max(1.0e-9, currentTrajectory);
		// Preserve a restrained service-like sensitivity change: the close end may
		// become at most about 1 dB louder after trajectory normalization, instead
		// of acting as another output-level control.
		// The three-amplitude normalizer samples a generic ellipse; real keyed
		// trajectories retain a small register-dependent residual at close gaps.
		// This measured correction raises the under-normalized middle register
		// without reviving the former upper-key level jump.
		const double intentionalSensitivity = std::pow(
			neutralKeyGap / pickupGap_, 0.08);
		pickupAlignmentGain_ = std::clamp(ReferenceGradientMagnitude() /
			std::max(1.0e-6, alignmentMagnitude) * trajectoryGain *
			intentionalSensitivity, 0.04, 7.0);
		pickupVoltageScale_ = 0.130 *
			ElectricPianoMechanicalTrim::PickupOutput * keyPickupSensitivity_ *
			pickupAlignmentGain_ * inverseReferenceGradient_;
		// Preserve a wide sound-design range while making the calibrated default
		// a quiet mechanical contribution rather than a parallel percussion layer.
		mechanicsLevel_ = 0.060 * mechanics * (0.30 + 0.70 * mechanics);
		cachedBodyWeight_ = body;
		cachedBellWeight_ = bell;
		cachedCouplingWeight_ = coupling;
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
			}
			modeRatio_[2] = 1.0018 + 0.0010 * (keyPosition_ - 0.5);
			for (std::size_t index = ElectricPianoAttackModeBegin;
				index < ElectricPianoAttackModeEnd; ++index)
				modeRatio_[index] = KeyboardModeRatio(index, keyPosition_);
			modeRatio_[ElectricPianoToneBarSubMode] =
				ToneBarSubModeFrequencyRatio(keyPosition_);
		}

		if (decayChanged)
		{
			const double releaseSeconds =
				ElectricPianoDamperReleaseSeconds(release);
			const double decayScale = std::pow(4.0, 2.0 * (decay - 0.5));
			// Published modal lifetimes are the default calibration. Bell may extend
			// or shorten the attack modes deliberately, but its factory 0.52 position
			// must not silently add the former 1.76% hidden correction.
			const double bellDecay = 1.0 + 0.38 * (bell - 0.52);
			// Material loss is small. Most fundamental loss occurs when an
			// unbalanced fork drives its resilient mounting block and harp. The
			// support coordinate has been analytically reduced into the reaction
			// factors above; anti-phase tine/bar motion cancels that force and raises
			// Q, while an isolated tine rings briefly. The keyboard term represents
			// increasing mounting loss toward the short, stiff treble assemblies.
			const double intrinsicForkLifetime = IntrinsicForkDecaySeconds *
				keyboardDecayScale;
			const double supportLossRate = 0.55 + 4.45 * keyPosition_;
			for (std::size_t index = 0; index < modes_.size(); ++index)
			{
				double lifetime;
				if (index < 2)
				{
					const double lossRate = 1.0 / intrinsicForkLifetime +
						supportLossRate * coupledSupportLossFactor_[index];
					lifetime = decayScale / std::max(1.0e-9, lossRate);
				}
				else if (index == ElectricPianoToneBarSubMode)
					lifetime = ToneBarSubModeDecaySeconds(keyPosition_) *
						decayScale;
				else
					lifetime = (index == 2 ? TransverseModeDecaySeconds *
						keyboardDecayScale :
						AttackModeDecaySeconds(index, keyPosition_)) *
						decayScale;
				if (index >= ElectricPianoAttackModeBegin &&
					index < ElectricPianoAttackModeEnd)
					lifetime *= bellDecay;
				if (damped)
					lifetime = std::min(lifetime, releaseSeconds *
						(index >= ElectricPianoAttackModeBegin &&
							index < ElectricPianoAttackModeEnd ? 0.55 : 1.0));
				modeRadius_[index] = std::exp(-1.0 /
					(std::max(0.002, lifetime) * sampleRate_));
				modePickupRadius_[index] = std::exp(-1.0 /
					(std::max(0.002, lifetime) * sampleRate_ *
						static_cast<double>(PickupOversamplingFactor)));
				modeSubstepRadius_[index] = std::exp(-1.0 /
					(std::max(0.002, lifetime) * sampleRate_ *
						static_cast<double>(ContactOversamplingFactor)));
			}
		}

		for (std::size_t index = 0; index < modes_.size(); ++index)
		{
			const double frequency = fundamental * modeRatio_[index];
			const double normalizedFrequency = frequency / sampleRate_;
			modeBandlimitGain_[index] = ElectricPianoModeBandlimitGain(
				frequency, sampleRate_);
			modeActive_[index] = normalizedFrequency < 0.49 &&
				modeBandlimitGain_[index] > 0.0;
			// A mode need not be audible to load the hammer. The contact solve runs
			// at 64x host rate, so retain the already modelled ultrasonic modes in
			// the point mobility during collision while keeping them out of pickup
			// reconstruction. The old single activity flag removed most higher
			// modes from upper-key collisions and concentrated their energy in the
			// last audible coordinate. This split changes no bass mode that was
			// already below the pickup bandwidth.
			if (contactActive_)
				contactModeActive_[index] = frequency < 0.49 * sampleRate_ *
					static_cast<double>(ContactOversamplingFactor);
			// Contact-only coefficients are needed only during the short collision.
			// Avoid paying their trigonometric update cost on every pitch-modulated
			// free-decay sample after the hammer has separated.
			const bool collisionModeNeeded = contactActive_ &&
				contactModeActive_[index];
			if (!modeActive_[index] && !collisionModeNeeded)
			{
				modes_[index] = {};
				modeRealCoefficient_[index] = 0.0;
				modeImaginaryCoefficient_[index] = 0.0;
				modePickupRealCoefficient_[index] = 0.0;
				modePickupImaginaryCoefficient_[index] = 0.0;
				modeSubstepRealCoefficient_[index] = 0.0;
				modeSubstepImaginaryCoefficient_[index] = 0.0;
				modeAngularFrequency_[index] = 0.0;
				continue;
			}
			const double angle = TwoPi * frequency / sampleRate_;
			if (modeActive_[index])
			{
				modeRealCoefficient_[index] = modeRadius_[index] * std::cos(angle);
				modeImaginaryCoefficient_[index] = modeRadius_[index] * std::sin(angle);
				const double pickupAngle = angle /
					static_cast<double>(PickupOversamplingFactor);
				modePickupRealCoefficient_[index] = modePickupRadius_[index] *
					std::cos(pickupAngle);
				modePickupImaginaryCoefficient_[index] = modePickupRadius_[index] *
					std::sin(pickupAngle);
			}
			else
			{
				modeRealCoefficient_[index] = 0.0;
				modeImaginaryCoefficient_[index] = 0.0;
				modePickupRealCoefficient_[index] = 0.0;
				modePickupImaginaryCoefficient_[index] = 0.0;
				if (!contactActive_)
					modes_[index] = {};
			}
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
		if (!modulationPathActive_)
		{
			for (std::size_t index = 0;
				index < ElectricPianoToneBarSubMode; ++index)
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
			constexpr std::size_t subMode = ElectricPianoToneBarSubMode;
			if (modeActive_[subMode])
			{
				const double oldReal = modes_[subMode].real;
				const double oldImaginary = modes_[subMode].imaginary;
				modes_[subMode].real = modeRealCoefficient_[subMode] * oldReal -
					modeImaginaryCoefficient_[subMode] * oldImaginary;
				modes_[subMode].imaginary = modeImaginaryCoefficient_[subMode] *
					oldReal + modeRealCoefficient_[subMode] * oldImaginary;
			}
			return;
		}
		// The pickup nonlinearity and creative FM/PM are evaluated at 4x. Advance
		// the physical resonators at that same cadence so their modal coordinates
		// do not have to be reconstructed from a host-rate aggregate signal.
		for (int frame = 0; frame < PickupOversamplingFactor; ++frame)
		{
			for (std::size_t index = 0;
				index < ElectricPianoToneBarSubMode; ++index)
			{
				if (!modeActive_[index])
				{
					pickupModeFrames_[frame][index] = {};
					continue;
				}
				const double oldReal = modes_[index].real;
				const double oldImaginary = modes_[index].imaginary;
				modes_[index].real = modePickupRealCoefficient_[index] * oldReal -
					modePickupImaginaryCoefficient_[index] * oldImaginary;
				modes_[index].imaginary = modePickupImaginaryCoefficient_[index] *
					oldReal + modePickupRealCoefficient_[index] * oldImaginary;
				pickupModeFrames_[frame][index] = modes_[index];
			}
			constexpr std::size_t subMode = ElectricPianoToneBarSubMode;
			if (!modeActive_[subMode])
				pickupModeFrames_[frame][subMode] = {};
			else
			{
				const double oldReal = modes_[subMode].real;
				const double oldImaginary = modes_[subMode].imaginary;
				modes_[subMode].real = modePickupRealCoefficient_[subMode] * oldReal -
					modePickupImaginaryCoefficient_[subMode] * oldImaginary;
				modes_[subMode].imaginary = modePickupImaginaryCoefficient_[subMode] *
					oldReal + modePickupRealCoefficient_[subMode] * oldImaginary;
				pickupModeFrames_[frame][subMode] = modes_[subMode];
			}
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
				if (!contactModeActive_[index])
					continue;
				position -= contactModeShape_[index] * modes_[index].real;
				velocity += contactModeShape_[index] *
					modeAngularFrequency_[index] * modes_[index].imaginary;
			}
			return std::array<double, 2>{position, velocity};
		};
		if (modulationPathActive_)
			for (auto& frame : pickupModeFrames_)
				frame.fill({});
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
						// Hunt--Crossley damping may reach zero on a sufficiently fast
						// unloading stroke, but a contact cannot pull the hammer back
						// toward the tine.  Do not cap the measured loading branch.
						const double hysteresis =
							1.0 + contactLoss_ * relativeVelocity;
						if (hysteresis <= 0.0 && relativeVelocity < 0.0)
						{
							// The dissipative fit has crossed the point at which its
							// algebraic force would become tensile. A real neoprene tip
							// cannot pull the receding hammer, so this is separation,
							// not a zero-force state held at positive penetration.
							hammerPosition_ = tine[0];
							contactActive_ = false;
							contactEngaged_ = false;
						}
						else
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
							if (!contactModeActive_[index] ||
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
				if (!contactModeActive_[index])
					continue;
				const double oldReal = modes_[index].real;
				const double oldImaginary = modes_[index].imaginary;
				modes_[index].real = modeSubstepRealCoefficient_[index] *
					oldReal - modeSubstepImaginaryCoefficient_[index] *
					oldImaginary;
				modes_[index].imaginary = modeSubstepImaginaryCoefficient_[index] *
					oldReal + modeSubstepRealCoefficient_[index] * oldImaginary;
			}
			if (modulationPathActive_ && (substep + 1) %
				(ContactOversamplingFactor / PickupOversamplingFactor) == 0)
			{
				const int frame = (substep + 1) /
					(ContactOversamplingFactor / PickupOversamplingFactor) - 1;
				for (std::size_t index = 0; index < modes_.size(); ++index)
					pickupModeFrames_[frame][index] = modeActive_[index] ?
						modes_[index] : Mode{};
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
		// Contact-only coordinates cannot be advanced at host rate once the
		// hammer separates. Their remaining ultrasonic energy is outside the
		// pickup model and is intentionally discarded here, after it has affected
		// the complete physical collision rather than before it.
		if (!contactActive_)
			for (std::size_t index = 0; index < modes_.size(); ++index)
				if (!modeActive_[index])
					modes_[index] = {};
	}

	static double KeyboardModeRatio(std::size_t index,
		double keyboardPosition)
	{
		if (index < 3)
			return 1.0;
		const std::size_t waveNumberIndex = index - 2;
		const double uniformRatio = std::pow(
			CantileverWaveNumbers[waveNumberIndex] /
			CantileverWaveNumbers[0], 2.0);
		auto applyShearCorrection = [&](double ratio)
		{
			// The measured F1/F3 ratios already contain the real tine's shear and
			// rotary-inertia effects. Above the final measurement, correct the
			// Euler--Bernoulli extrapolation with the first-order Timoshenko/Rayleigh
			// factor cited by Pfeifle and Sonderbo. It is normalized by the
			// fundamental so pitch is unchanged; only the increasingly sensitive
			// upper-mode ratios of short tines move. Blend in only beyond F3 so the
			// measured anchors and the successful bass remain exact.
			constexpr double PoissonRatio = 0.29;
			constexpr double CircularShearCoefficient =
				6.0 * (1.0 + PoissonRatio) / (7.0 + 6.0 * PoissonRatio);
			const double youngsToShearRatio = 2.0 * (1.0 + PoissonRatio);
			const double tineLength =
				ElectricPianoPublishedMechanicalData::LongestTineMetres * std::pow(
					ElectricPianoPublishedMechanicalData::ShortestTineMetres /
						ElectricPianoPublishedMechanicalData::LongestTineMetres,
					std::clamp(keyboardPosition, 0.0, 1.0));
			const double radiusToLength =
				ElectricPianoPublishedMechanicalData::TineRadiusMetres /
				std::max(1.0e-9, tineLength);
			const double correctionCoefficient = 0.25 * radiusToLength *
				radiusToLength * (1.0 + youngsToShearRatio /
					CircularShearCoefficient);
			const double correction = std::sqrt(
				(1.0 + correctionCoefficient * CantileverWaveNumbers[0] *
					CantileverWaveNumbers[0]) /
				(1.0 + correctionCoefficient *
					CantileverWaveNumbers[waveNumberIndex] *
					CantileverWaveNumbers[waveNumberIndex]));
			constexpr double F3Position = (53.0 - 28.0) / 72.0;
			constexpr double A4Position = (69.0 - 28.0) / 72.0;
			const double blendLinear = std::clamp((keyboardPosition - F3Position) /
				(A4Position - F3Position), 0.0, 1.0);
			const double blend = blendLinear * blendLinear *
				(3.0 - 2.0 * blendLinear);
			return ratio * std::exp(blend * std::log(correction));
		};
		// Gabrielli/Cantarini report the first measured overtones of F1 at
		// 7.2/20.6 and those of F3 at 7.4/20.7/38.7. A Rayleigh--Ritz cantilever with
		// one sliding point mass fits both anchors with tuning-mass fractions
		// 0.146/0.177 and positions 0.859L/0.840L. Extrapolating those physical
		// parameters (mass fraction versus tine length and spring position versus
		// key) yields the compact ratio polynomials below. This replaces the former
		// return to an unloaded A4 beam even though every real key retains a tuning
		// spring. The fit is explicitly provisional above F3 until the article's
		// full-key data or new SLDV measurements are available. The very lowest key
		// retains the calibrated 7.11/20.25 ratios.
		constexpr double F1Position = (29.0 - 28.0) / 72.0;
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		auto smoothInterpolate = [](double position, double lowerPosition,
			double upperPosition, double lowerValue, double upperValue)
		{
			const double linear = std::clamp((position - lowerPosition) /
				(upperPosition - lowerPosition), 0.0, 1.0);
			const double smooth = linear * linear * (3.0 - 2.0 * linear);
			return lowerValue + smooth * (upperValue - lowerValue);
		};
		if (waveNumberIndex >= 1 && waveNumberIndex <= 3)
		{
			const std::size_t measuredIndex = waveNumberIndex - 1;
			constexpr std::array<double, 3> LowestKeyRatios{
				7.11, 20.25, 34.386};
			constexpr std::array<double, 3> F1Ratios{7.2, 20.6, 34.386};
			constexpr std::array<double, 3> F3Ratios{7.4, 20.7, 38.7};
			if (keyboardPosition < F1Position)
				return applyShearCorrection(smoothInterpolate(keyboardPosition, 0.0,
					F1Position, LowestKeyRatios[measuredIndex],
					F1Ratios[measuredIndex]));
			if (keyboardPosition < F3Position)
				return applyShearCorrection(smoothInterpolate(keyboardPosition,
					F1Position, F3Position, F1Ratios[measuredIndex],
					F3Ratios[measuredIndex]));
			const double pointMassPosition = std::clamp(
				(keyboardPosition - F3Position) / (1.0 - F3Position), 0.0, 1.0);
			const double squared = pointMassPosition * pointMassPosition;
			constexpr std::array<double, 3> LinearPointMassCoefficients{
				0.3887, 0.2667, -1.1750};
			constexpr std::array<double, 3> QuadraticPointMassCoefficients{
				0.0022, -0.3675, 0.6190};
			return applyShearCorrection(F3Ratios[measuredIndex] +
				LinearPointMassCoefficients[measuredIndex] * pointMassPosition +
				QuadraticPointMassCoefficients[measuredIndex] * squared);
		}
		return applyShearCorrection(uniformRatio);
	}

	static std::array<double, 9> SpringLoadedModeReduction(std::size_t index,
		double keyboardPosition)
	{
		// Offline Rayleigh--Ritz reduction of a uniform cantilever plus the sliding
		// tuning mass fitted at F1/F3. Components 0--7 expand a tip-normalized loaded
		// eigenfunction in the first eight tip-normalized cantilever shapes; component
		// 8 is generalized modal mass divided by the quarter-beam reference mass.
		// Fifth-order keyed polynomials reproduce the 73-key eigensolve with maximum
		// coefficient error below 6.5e-6 and mass-ratio error below 2.3e-5. Keeping
		// this solve out of the audio path also preserves cheap continuous tuning.
		constexpr std::array<std::array<std::array<double, 6>, 9>, 3>
			ReductionPolynomials{{
			{{
				{{-0.0956129, 0.0456194, 0.00777974, -0.000116554,
					-9.46565e-05, -3.09906e-06}},
				{{1.10324, -0.0455643, -0.0104353, -9.22499e-05,
					0.000148787, 1.19656e-05}},
				{{-0.00398375, -0.00128955, 0.00183961, 0.000294257,
					-1.9185e-05, -8.84302e-06}},
				{{-0.00220983, 0.000352183, 0.000614798, 2.54259e-05,
					-1.9119e-05, -2.65702e-06}},
				{{-0.000963916, 0.000374115, 0.000208078, -2.96839e-05,
					-1.20937e-05, -3.05088e-07}},
				{{-0.000370971, 0.000262402, 4.75427e-05, -3.62169e-05,
					-5.54762e-06, 7.79864e-07}},
				{{-0.000102739, 0.000160785, -1.67049e-05, -2.82371e-05,
					-6.24785e-07, 1.14076e-06}},
				{{8.68671e-06, 8.48909e-05, -3.78117e-05, -1.67413e-05,
					2.43924e-06, 1.01795e-06}},
				{{1.24536, -0.130942, -0.0154338, 0.00239892,
					0.00042329, -1.50138e-05}}
			}},
			{{
				{{0.0856358, 0.0718697, 0.00239745, -0.000800856,
					0.000103003, 2.45192e-05}},
				{{0.0264046, 0.00725223, -0.0122908, -0.00115264,
					9.91097e-05, -1.35701e-05}},
				{{0.856621, -0.10777, 0.0105401, 0.00509314,
					-0.000285909, -0.000179625}},
				{{0.0204186, 0.0217557, 0.00242873, -0.00193329,
					-0.00012666, 0.000105396}},
				{{0.00748901, 0.00635952, -0.000622121, -0.00081385,
					2.26418e-05, 4.0268e-05}},
				{{0.00274673, 0.00146424, -0.00106307, -0.000356234,
					6.33963e-05, 1.89727e-05}},
				{{0.000746877, -0.000241582, -0.00086532, -9.23743e-05,
					6.8684e-05, 6.35504e-06}},
				{{-6.27144e-05, -0.000689937, -0.000524997, 5.61063e-05,
					5.57334e-05, -2.31667e-06}},
				{{0.758331, -0.147027, 0.0427399, 0.00513243,
					-0.00150902, -5.12784e-05}}
			}},
			{{
				{{0.179255, 0.069498, -0.000834428, -9.1712e-05,
					-0.000210303, -0.000145463}},
				{{0.0500679, -0.00904389, -0.0123546, -0.000331001,
					5.78134e-07, 7.40167e-07}},
				{{-0.0772817, -0.0930271, -0.0253526, 0.000664447,
					0.000802763, 0.000230014}},
				{{0.750814, 0.0193002, 0.0563247, -0.000129465,
					-0.00115623, 0.000269954}},
				{{0.0702597, 0.022177, -0.00752142, -0.000515268,
					0.000267363, -0.000236803}},
				{{0.021748, 0.000528111, -0.00541273, -0.000221897,
					0.000127871, -7.32499e-05}},
				{{0.00559402, -0.00461574, -0.00332659, 0.000191674,
					0.000100441, -2.94211e-05}},
				{{-0.000457813, -0.00481666, -0.00152236, 0.000433223,
					6.75157e-05, -1.57717e-05}},
				{{0.680285, 0.112127, 0.101452, 0.00504365,
					0.00165225, -7.54704e-05}}
			}}
		}};
		std::array<double, 9> reduction{};
		if (index < 3 || index > 5)
		{
			reduction[8] = 1.0;
			return reduction;
		}
		const double coordinate = 2.0 * std::clamp(keyboardPosition, 0.0, 1.0) - 1.0;
		const auto& polynomials = ReductionPolynomials[index - 3];
		for (std::size_t component = 0; component < reduction.size(); ++component)
		{
			double value = polynomials[component].back();
			for (std::size_t term = polynomials[component].size() - 1;
				term-- > 0;)
				value = value * coordinate + polynomials[component][term];
			reduction[component] = value;
		}
		reduction[8] = std::clamp(reduction[8], 0.5, 1.5);
		return reduction;
	}

	static double AttackModeModalMassMultiplier(std::size_t index,
		double keyboardPosition)
	{
		// Modes 3--5 use the spring-loaded Rayleigh--Ritz eigenvectors fitted to
		// the measured F1/F3 ratios. Their generalized masses therefore differ from
		// the quarter-beam mass used by KeyProfile. Higher, unconstrained modes keep
		// the uniform-beam normalization until measurements justify loading them.
		if (index < 3 || index > 5)
			return 1.0;
		return SpringLoadedModeReduction(index, keyboardPosition)[8];
	}

	static double BeamAttackModeDecaySeconds(std::size_t index,
		double keyboardPosition)
	{
		const double tineLength =
			ElectricPianoPublishedMechanicalData::LongestTineMetres * std::pow(
				ElectricPianoPublishedMechanicalData::ShortestTineMetres /
				ElectricPianoPublishedMechanicalData::LongestTineMetres,
				std::clamp(keyboardPosition, 0.0, 1.0));
		const std::size_t waveNumberIndex = index - 2;
		const double spatialWaveNumber =
			CantileverWaveNumbers[waveNumberIndex] / tineLength;
		const double lossRate =
			ElectricPianoPublishedMechanicalData::FrequencyIndependentLossPerSecond +
			ElectricPianoPublishedMechanicalData::
				FrequencyDependentLossSquareMetresPerSecond *
			spatialWaveNumber * spatialWaveNumber;
		return 1.0 / std::max(1.0e-9, lossRate);
	}

	static double AttackModeDecaySeconds(std::size_t index,
		double keyboardPosition)
	{
		const double beamDecay = BeamAttackModeDecaySeconds(index,
			keyboardPosition);
		const std::size_t waveNumberIndex = index - 2;
		if (waveNumberIndex < 1 || waveNumberIndex > 3)
			return beamDecay;

		// Table 8.6 of Cantarini/Gabrielli gives amplitude-decay slopes in
		// dB/s. Convert with tau = 20 log10(e)/|slope| and use the resulting
		// measured lifetime as a multiplier on Sonderbo's distributed-loss law.
		// Retaining the beam law between anchors preserves register scaling; the
		// multiplier admits the observed non-monotonic mode-dependent losses that a
		// single sigma0 + sigma1*k^2 curve cannot reproduce. The F3 multiplier must
		// not remain frozen above the last measured damping anchor: doing that made
		// modes four and five roughly eleven and five times too long throughout the
		// treble. In the absence of treble damping measurements, relax the measured
		// anomaly independently by A4 to the sourced distributed-loss law. Do not tie
		// this provisional loss coordinate to the point-mass frequency extrapolation:
		// a tuning spring can shift frequency without preserving F3's anomalous loss.
		constexpr double DbPerNeper = 8.685889638065037;
		constexpr double F1Position = (29.0 - 28.0) / 72.0;
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		constexpr std::array<double, 3> F1DecaySeconds{
			DbPerNeper / 21.1, DbPerNeper / 67.7, 0.0};
		constexpr std::array<double, 3> F3DecaySeconds{
			DbPerNeper / 294.0, DbPerNeper / 37.0, DbPerNeper / 161.0};
		const std::size_t measuredIndex = waveNumberIndex - 1;
		const double f1Multiplier = F1DecaySeconds[measuredIndex] > 0.0 ?
			F1DecaySeconds[measuredIndex] /
				BeamAttackModeDecaySeconds(index, F1Position) : 1.0;
		const double f3Multiplier = F3DecaySeconds[measuredIndex] /
			BeamAttackModeDecaySeconds(index, F3Position);
		if (keyboardPosition > F3Position)
		{
			constexpr double A4Position = (69.0 - 28.0) / 72.0;
			const double relaxationLinear = std::clamp(
				(keyboardPosition - F3Position) / (A4Position - F3Position),
				0.0, 1.0);
			const double relaxation = relaxationLinear * relaxationLinear *
				(3.0 - 2.0 * relaxationLinear);
			const double measuredLossPerturbation = 1.0 - relaxation;
			return beamDecay * std::exp(measuredLossPerturbation *
				std::log(std::max(1.0e-12, f3Multiplier)));
		}
		const double linear = std::clamp((keyboardPosition - F1Position) /
			(F3Position - F1Position), 0.0, 1.0);
		const double smooth = linear * linear * (3.0 - 2.0 * linear);
		const double logMultiplier = std::log(std::max(1.0e-12, f1Multiplier)) +
			smooth * (std::log(std::max(1.0e-12, f3Multiplier)) -
				std::log(std::max(1.0e-12, f1Multiplier)));
		return beamDecay * std::exp(logMultiplier);
	}

	static double AttackModeEnergyCalibration(std::size_t index,
		double keyboardPosition)
	{
		// The measured F1/F3 decay constants describe modal residues whose
		// magnitudes were fitted at the same time. Substituting only those decay
		// constants changes a mode's total observed energy in direct proportion to
		// its lifetime. Preserve the existing reduced model's broadband impulse
		// calibration by scaling the contact projection by sqrt(tau_beam/tau_fit).
		// This is an explicit residue normalization; it does not alter the loaded
		// eigenfunction, generalized mass, frequency, or measured decay constant.
		const double fittedLifetime = AttackModeDecaySeconds(index,
			keyboardPosition);
		const double beamLifetime = BeamAttackModeDecaySeconds(index,
			keyboardPosition);
		return std::sqrt(std::clamp(beamLifetime /
			std::max(1.0e-9, fittedLifetime), 0.05, 4.0));
	}

	static double ToneBarSubModeFrequencyRatio(double keyboardPosition)
	{
		// Gabrielli/Cantarini SLDV measurements locate the tone-bar submode at
		// 0.83*f0 for F1 and 0.58*f0 for F3. The openly available measurements
		// stop there. Above F3, aligned direct-harp recordings show the strongest
		// early submode settling near 0.4--0.5*f0 rather than rising toward f0.
		constexpr double F1Position = (29.0 - 28.0) / 72.0;
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		if (keyboardPosition <= F1Position)
			return 0.83;
		if (keyboardPosition <= F3Position)
		{
			const double linear = (keyboardPosition - F1Position) /
				(F3Position - F1Position);
			const double smooth = linear * linear * (3.0 - 2.0 * linear);
			return std::exp(std::log(0.83) + smooth *
				(std::log(0.58) - std::log(0.83)));
		}
		const double upperPosition = (keyboardPosition - F3Position) /
			(1.0 - F3Position);
		return 0.42 + 0.16 * std::exp(-5.0 * upperPosition);
	}

	static double ToneBarSubModeDecaySeconds(double keyboardPosition)
	{
		// Table 8.6 gives amplitude slopes of -9.1 dB/s at F1 and
		// -138 dB/s at F3. Convert with tau=20*log10(e)/|slope|. The
		// upper continuation is constrained by direct-harp spectrograms, where
		// the line falls by roughly 25--40 dB between 100 and 250 ms.
		constexpr double DbPerNeper = 8.685889638065037;
		constexpr double F1Lifetime = DbPerNeper / 9.1;
		constexpr double F3Lifetime = DbPerNeper / 138.0;
		constexpr double F1Position = (29.0 - 28.0) / 72.0;
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		if (keyboardPosition <= F1Position)
			return F1Lifetime;
		if (keyboardPosition <= F3Position)
		{
			const double linear = (keyboardPosition - F1Position) /
				(F3Position - F1Position);
			const double smooth = linear * linear * (3.0 - 2.0 * linear);
			return std::exp(std::log(F1Lifetime) + smooth *
				(std::log(F3Lifetime) - std::log(F1Lifetime)));
		}
		const double upperPosition = (keyboardPosition - F3Position) /
			(1.0 - F3Position);
		return F3Lifetime * std::pow(0.035 / F3Lifetime, upperPosition);
	}

	static double ToneBarSubModeTineParticipation(double keyboardPosition)
	{
		// A normal mode is observed and forced through the same tine component.
		// The two SLDV anchors imply a weak F1 residue (-28 dB relative to the
		// played mode) and a much stronger F3 residue (-4.5 dB). Above F3 the
		// direct-harp corpus shows a rapid reduction; by the upper register the
		// submode and its pickup sidebands are 40--60 dB below the fundamental.
		// These participation values are an initial reciprocal modal fit; unlike
		// the former observation gain they enter twice, at force and pickup.
		constexpr double F1Position = (29.0 - 28.0) / 72.0;
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		constexpr std::array<double, 8> KeyPositions{
			F1Position, F3Position, (55.0 - 28.0) / 72.0,
			(59.0 - 28.0) / 72.0, (65.0 - 28.0) / 72.0,
			(76.0 - 28.0) / 72.0, (88.0 - 28.0) / 72.0, 1.0};
		constexpr std::array<double, 8> Participations{
			0.177, 0.56, 0.40, 0.25, 0.15, 0.060, 0.045, 0.040};
		if (keyboardPosition <= KeyPositions.front())
			return Participations.front();
		for (std::size_t upper = 1; upper < KeyPositions.size(); ++upper)
		{
			if (keyboardPosition > KeyPositions[upper])
				continue;
			const double linear = (keyboardPosition - KeyPositions[upper - 1]) /
				(KeyPositions[upper] - KeyPositions[upper - 1]);
			const double smooth = linear * linear * (3.0 - 2.0 * linear);
			return std::exp(std::log(Participations[upper - 1]) + smooth *
				(std::log(Participations[upper]) -
					std::log(Participations[upper - 1])));
		}
		return Participations.back();
	}

	static double StrikePositionFromControl(double keyboardPosition,
		double strikeControl)
	{
		// The service procedure does not describe a straight line: it asks the
		// technician to set C4, F3 and C3 independently, then accept the intervening
		// keys when they fall close to maximum power without a thunk. A linear
		// 0.40L--0.20L approximation put C4 and much of the upper register close to
		// higher-mode nodes. A full-keyboard render sweep found the same optimum the
		// manual describes: a fairly flat bass line, a pronounced bend around the
		// C3--C4 checkpoints, and a gentler approach toward the treble clamp.
		//
		// Positions remain explicit listening trims because no published factory
		// jig dimensions identify the contact point on every tine. Keeping the six
		// checkpoints here makes a later measured-harp calibration local and avoids
		// hiding the correction in per-mode gains.
		constexpr std::array<double, 6> KeyboardCheckpoints{
			0.0, (48.0 - 28.0) / 72.0, (53.0 - 28.0) / 72.0,
			(60.0 - 28.0) / 72.0, (84.0 - 28.0) / 72.0, 1.0};
		constexpr std::array<double, 6> FactoryStrikePositions{
			0.38, 0.29, 0.22, 0.205, 0.16, 0.14};
		const double position = std::clamp(keyboardPosition, 0.0, 1.0);
		std::size_t lower = 0;
		while (lower + 2 < KeyboardCheckpoints.size() &&
			position > KeyboardCheckpoints[lower + 1])
			++lower;
		const double interval = KeyboardCheckpoints[lower + 1] -
			KeyboardCheckpoints[lower];
		const double linearFraction = std::clamp(
			(position - KeyboardCheckpoints[lower]) / interval, 0.0, 1.0);
		const double smoothFraction = linearFraction * linearFraction *
			(3.0 - 2.0 * linearFraction);
		const double factoryPosition = FactoryStrikePositions[lower] +
			smoothFraction * (FactoryStrikePositions[lower + 1] -
				FactoryStrikePositions[lower]);
		// The panel now represents a signed physical displacement of the harp's
		// striking line about the factory point. A full-range movement tapers from
		// 6 mm in the bass to 1 mm in the treble, consistent with the much shorter
		// upper tines. Unlike interpolation toward 0.04L/0.96L, this is smooth at
		// centre and cannot fling adjacent high-register settings across several
		// modal nodes. Centre remains bit-for-bit the calibrated factory position.
		const double tineLength =
			ElectricPianoPublishedMechanicalData::LongestTineMetres * std::pow(
				ElectricPianoPublishedMechanicalData::ShortestTineMetres /
					ElectricPianoPublishedMechanicalData::LongestTineMetres,
				position);
		const double maximumOffsetMetres = 0.006 + position * (0.001 - 0.006);
		const double signedControl = 2.0 *
			(std::clamp(strikeControl, 0.0, 1.0) - 0.5);
		return std::clamp(factoryPosition + signedControl *
			maximumOffsetMetres / tineLength, 0.04, 0.96);
	}

	static double UniformCantileverModeShape(std::size_t waveNumberIndex,
		double position)
	{
		const double beta = CantileverWaveNumbers[waveNumberIndex];
		const double sigma = (std::cosh(beta) + std::cos(beta)) /
			(std::sinh(beta) + std::sin(beta));
		auto shape = [&](double normalizedPosition)
		{
			const double argument = beta * normalizedPosition;
			return std::cosh(argument) - std::cos(argument) - sigma *
				(std::sinh(argument) - std::sin(argument));
		};
		return std::clamp(shape(std::clamp(position, 0.0, 1.0)) / shape(1.0),
			-1.5, 1.5);
	}

	static double CantileverModeShape(std::size_t index, double position)
	{
		const std::size_t waveNumberIndex = index < 3 ? 0 : index - 2;
		return UniformCantileverModeShape(waveNumberIndex, position);
	}

	static double FiniteContactModeProjection(std::size_t index,
		double keyboardPosition, double centre, double contactWidthMetres)
	{
		// Integrate the calibrated physical contact strip over the modal shape. The
		// symmetric binomial weights approximate the peaked pressure of a compliant
		// contact; this runs once per strike, not in the contact substeps.
		const double tineLength =
			ElectricPianoPublishedMechanicalData::LongestTineMetres * std::pow(
				ElectricPianoPublishedMechanicalData::ShortestTineMetres /
					ElectricPianoPublishedMechanicalData::LongestTineMetres,
				std::clamp(keyboardPosition, 0.0, 1.0));
		const double halfWidth = std::clamp(0.5 * contactWidthMetres /
			std::max(1.0e-6, tineLength), 0.001, 0.14);
		constexpr std::array<double, 5> Offsets{-1.0, -0.5, 0.0, 0.5, 1.0};
		constexpr std::array<double, 5> Weights{1.0, 4.0, 6.0, 4.0, 1.0};
		const bool springLoadedMode = index >= 3 && index <= 5;
		const auto loadedReduction = springLoadedMode ?
			SpringLoadedModeReduction(index, keyboardPosition) :
			std::array<double, 9>{};
		auto modeShape = [&](double position)
		{
			if (!springLoadedMode)
				return CantileverModeShape(index, position);
			double shape = 0.0;
			for (std::size_t basis = 0; basis < 8; ++basis)
				shape += loadedReduction[basis] *
					UniformCantileverModeShape(basis, position);
			return std::clamp(shape, -1.5, 1.5);
		};
		double projection = 0.0;
		for (std::size_t sample = 0; sample < Offsets.size(); ++sample)
			projection += Weights[sample] * modeShape(
				centre + halfWidth * Offsets[sample]);
		return projection / 16.0;
	}

	struct MagneticFieldSample
	{
		double vertical{};
		double radial{};
	};

	static MagneticFieldSample EvaluateCalibratedPoleGradient(double vertical,
		double radialDistance)
	{
		// Pfeifle derives a three-part integral over an idealised pole
		// cross-section, but publishes neither production dimensions nor a flux scan.
		// Falaize's sourced circular-pole reduction was also evaluated here; as a
		// keyboard-wide replacement for this asymmetric pole it over-produced bass H2
		// by 5--14x and reduced rather than restored the treble harmonic ladder. Keep
		// the measured-listening-calibrated asymmetric edge reduction until a real 2D
		// pickup scan exists. It remains one scalar flux construction, differentiated
		// consistently for arbitrary two-dimensional tine motion.
		const double radialFalloff = PowNegative1p3(radialDistance);
		MagneticFieldSample result;
		auto accumulateEdge = [&](double edgePosition, double weight,
			double edgeRadius)
		{
			constexpr double GapBroadening = 0.180;
			const double width = edgeRadius + GapBroadening * radialDistance;
			const double inverseWidth = 1.0 / std::max(1.0e-9, width);
			const double displacement = vertical - edgePosition;
			const double argument = displacement * inverseWidth;
			const double transition = TanhPade76(argument);
			const double transitionDerivative = 1.0 - transition * transition;
			const double amplitude = weight * radialFalloff;
			result.vertical += amplitude * transitionDerivative * inverseWidth;
			const double amplitudeDerivative = -FieldFalloff * amplitude /
				std::max(1.0e-9, radialDistance);
			const double argumentDerivative = -displacement * GapBroadening *
				inverseWidth * inverseWidth;
			result.radial += amplitudeDerivative * transition + amplitude *
				transitionDerivative * argumentDerivative;
		};
		accumulateEdge(0.34, 1.0, 0.030);
		accumulateEdge(-0.46, 0.20, 0.095);
		return result;
	}

	struct MagneticFieldTable
	{
		// The primary edge is only 0.03 mm wide at zero clearance. A dense
		// process-wide table keeps interpolation below the direct-field error bound
		// without changing the four-lookups-per-sample runtime cost.
		static constexpr std::size_t VerticalSamples = 241;
		static constexpr std::size_t RadialSamples = 135;
		static constexpr double VerticalMinimum = -1.2;
		static constexpr double VerticalMaximum = 1.8;
		static constexpr double RadialMinimum = 0.45;
		static constexpr double RadialMaximum = 3.8;
		std::array<MagneticFieldSample,
			VerticalSamples * RadialSamples> samples{};

		MagneticFieldTable()
		{
			for (std::size_t verticalIndex = 0;
				verticalIndex < VerticalSamples; ++verticalIndex)
			{
				const double vertical = VerticalMinimum +
					(VerticalMaximum - VerticalMinimum) *
						static_cast<double>(verticalIndex) /
						static_cast<double>(VerticalSamples - 1);
				for (std::size_t radialIndex = 0;
					radialIndex < RadialSamples; ++radialIndex)
				{
					const double radial = RadialMinimum +
						(RadialMaximum - RadialMinimum) *
							static_cast<double>(radialIndex) /
							static_cast<double>(RadialSamples - 1);
					samples[verticalIndex * RadialSamples + radialIndex] =
						EvaluateCalibratedPoleGradient(vertical, radial);
				}
			}
		}

		MagneticFieldSample Interpolate(double vertical, double radial) const
		{
			const double verticalCoordinate = std::clamp(
				(vertical - VerticalMinimum) * (VerticalSamples - 1) /
					(VerticalMaximum - VerticalMinimum), 0.0,
				static_cast<double>(VerticalSamples - 1));
			const double radialCoordinate = std::clamp(
				(radial - RadialMinimum) * (RadialSamples - 1) /
					(RadialMaximum - RadialMinimum), 0.0,
				static_cast<double>(RadialSamples - 1));
			const std::size_t vertical0 = std::min(
				static_cast<std::size_t>(verticalCoordinate), VerticalSamples - 2);
			const std::size_t radial0 = std::min(
				static_cast<std::size_t>(radialCoordinate), RadialSamples - 2);
			const double verticalFraction = verticalCoordinate -
				static_cast<double>(vertical0);
			const double radialFraction = radialCoordinate -
				static_cast<double>(radial0);
			auto at = [&](std::size_t v, std::size_t r) -> const MagneticFieldSample&
			{
				return samples[v * RadialSamples + r];
			};
			auto blend = [&](double MagneticFieldSample::*component)
			{
				const double lower = at(vertical0, radial0).*component +
					radialFraction * (at(vertical0, radial0 + 1).*component -
						at(vertical0, radial0).*component);
				const double upper = at(vertical0 + 1, radial0).*component +
					radialFraction * (at(vertical0 + 1, radial0 + 1).*component -
						at(vertical0 + 1, radial0).*component);
				return lower + verticalFraction * (upper - lower);
			};
			return {blend(&MagneticFieldSample::vertical),
				blend(&MagneticFieldSample::radial)};
		}
	};

	static std::array<double, 2> MagneticFluxGradient(double vertical,
		double horizontal, double gap)
	{
		// The table is two-dimensional in alignment and radial clearance. Horizontal
		// tine motion changes that clearance, so the second component follows by the
		// chain rule. Construction occurs once per process; every voice then pays
		// only a bilinear lookup at the 4x pickup rate.
		constexpr double HorizontalFieldScale = 0.62;
		const double radialDistance = std::sqrt(gap * gap +
			HorizontalFieldScale * horizontal * horizontal + 0.020);
		static const MagneticFieldTable Table;
		const auto field = Table.Interpolate(vertical, radialDistance);
		return {field.vertical, field.radial * HorizontalFieldScale * horizontal /
			std::max(1.0e-9, radialDistance)};
	}

	static double MagneticTrajectoryRms(double verticalOffset,
		double horizontalOffset, double gap, double nominalExcursion)
	{
		// The normalization trajectory is dimensionless in time, so its angular
		// frequency cancels from the ratio. A small quadrature component represents
		// the measured two-polarization orbit without claiming a per-key ellipse.
		// Average low, nominal and high excursions geometrically: one arbitrary
		// MIDI velocity then cannot dominate the pickup calibration.
		constexpr std::array<double, 3> ExcursionScales{0.45, 1.0, 1.8};
		constexpr std::array<double, 12> PhaseCosines{
			0.9659258263, 0.7071067812, 0.2588190451, -0.2588190451,
			-0.7071067812, -0.9659258263, -0.9659258263, -0.7071067812,
			-0.2588190451, 0.2588190451, 0.7071067812, 0.9659258263};
		constexpr std::array<double, 12> PhaseSines{
			0.2588190451, 0.7071067812, 0.9659258263, 0.9659258263,
			0.7071067812, 0.2588190451, -0.2588190451, -0.7071067812,
			-0.9659258263, -0.9659258263, -0.7071067812, -0.2588190451};
		double logarithmicRms = 0.0;
		for (double excursionScale : ExcursionScales)
		{
			const double verticalAmplitude = nominalExcursion * excursionScale;
			const double horizontalAmplitude = 0.12 * verticalAmplitude;
			double energy = 0.0;
			for (std::size_t sample = 0; sample < PhaseCosines.size(); ++sample)
			{
				const double cosine = PhaseCosines[sample];
				const double sine = PhaseSines[sample];
				const auto gradient = MagneticFluxGradient(
					verticalOffset + verticalAmplitude * cosine,
					horizontalOffset + horizontalAmplitude * sine, gap);
				const double emf = gradient[0] * (-verticalAmplitude * sine) +
					gradient[1] * (horizontalAmplitude * cosine);
				energy += emf * emf;
			}
			const double rms = std::sqrt(energy /
				static_cast<double>(PhaseCosines.size()));
			// The explicit branch both handles a hypothetical NaN and makes the
			// strictly-positive logarithm precondition visible to static analysis.
			const double positiveRms = rms > 1.0e-12 ? rms : 1.0e-12;
			logarithmicRms += std::log(positiveRms);
		}
		return std::exp(logarithmicRms /
			static_cast<double>(ExcursionScales.size()));
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

	static std::array<double, 2> ModulationPhaseRotation(double angle)
	{
		// Normal musical FM moves only a small fraction of a cycle per 4x frame.
		// This seventh/sixth-order local rotation is effectively exact in that
		// range and avoids two transcendental calls per mode. Extreme modulation
		// retains the full library implementation rather than range-reducing a
		// polynomial outside its error bound.
		if (std::abs(angle) < 0.35)
		{
			const double squared = angle * angle;
			const double sine = angle * (1.0 + squared * (-1.0 / 6.0 +
				squared * (1.0 / 120.0 - squared / 5040.0)));
			const double cosine = 1.0 + squared * (-1.0 / 2.0 + squared *
				(1.0 / 24.0 - squared / 720.0));
			return {cosine, sine};
		}
		return {std::cos(angle), std::sin(angle)};
	}

	static constexpr double TwoPi = 6.283185307179586476925286766559;
	static constexpr double FieldFalloff = 1.30;
	// Modal coordinates and velocities are now metres and metres/second. The
	// pickup field below is expressed in millimetres. PickupExcursion is the
	// separately documented fit needed because no measured flux map is available.
	static constexpr double TineDisplacementScale = 1000.0 *
		ElectricPianoMechanicalTrim::PickupExcursion;
	// The published SI stiffness is roughly five orders above the former
	// model-space fit. Sixty-four symplectic substeps keep collision spectra
	// sample-rate invariant; this path runs only for the few milliseconds of
	// actual contact and does not affect sustained-voice CPU.
	static constexpr int ContactOversamplingFactor = 64;
	static constexpr double MinimumHammerContactVelocity = 0.012;
	static constexpr int PickupOversamplingFactor =
		X4Resampler_Order7::ResamplingFactor;
	static constexpr double ReferencePickupGap = 1.6;
	static constexpr double ReferencePickupVertical = 0.34 + 0.22 * 0.55 * 0.55 +
		0.020 * ((60.0 - 28.0) / 72.0);
	static constexpr double ReferencePickupHorizontal = 0.10 + 0.035 *
		((60.0 - 28.0) / 72.0);
	// The first three coordinates form the audible coupled tine/tone-bar body.
	// Broadly FM-modulating the intentionally inharmonic, short-lived attack
	// modes produces a sideband cloud rather than the coherent overtones players
	// expect from FM. Keep those impact coordinates physical and let the nonlinear
	// magnetic pickup derive the upper modulation spectrum from the body motion.
	static constexpr std::array<double, ElectricPianoModeCount>
		ModulationModeWeights{
			1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	static constexpr std::array<double, 8> CantileverWaveNumbers{
		1.875104, 4.694091, 7.854757, 10.995541, 14.137168, 17.278760,
		20.420352, 23.561945};
	static constexpr double IntrinsicForkDecaySeconds = 6.8;
	static constexpr double TransverseModeDecaySeconds = 1.9;

	std::array<Mode, ElectricPianoModeCount> modes_{};
	std::array<double, ElectricPianoModeCount> modeRealCoefficient_{};
	std::array<double, ElectricPianoModeCount> modeImaginaryCoefficient_{};
	std::array<double, ElectricPianoModeCount> modeRadius_{};
	std::array<double, ElectricPianoModeCount> modePickupRealCoefficient_{};
	std::array<double, ElectricPianoModeCount> modePickupImaginaryCoefficient_{};
	std::array<double, ElectricPianoModeCount> modePickupRadius_{};
	std::array<double, ElectricPianoModeCount> modeSubstepRealCoefficient_{};
	std::array<double, ElectricPianoModeCount> modeSubstepImaginaryCoefficient_{};
	std::array<double, ElectricPianoModeCount> modeSubstepRadius_{};
	std::array<double, ElectricPianoModeCount> modeAngularFrequency_{};
	std::array<double, ElectricPianoModeCount> modeOutputWeight_{};
	std::array<double, ElectricPianoModeCount> modeHorizontalWeight_{};
	std::array<double, ElectricPianoModeCount> modeRatio_{};
	std::array<double, ElectricPianoModeCount> modeBandlimitGain_{};
	std::array<bool, ElectricPianoModeCount> modeActive_{};
	std::array<bool, ElectricPianoModeCount> contactModeActive_{};
	std::array<double, ElectricPianoModeCount> contactModeShape_{};
	std::array<double, ElectricPianoModeCount> modeInverseMass_{};
	std::array<double, 2> coupledToneBarRatio_{};
	std::array<double, 2> coupledSupportLossFactor_{};
	double sampleRate_ = 48000.0;
	double latchedVelocity_{};
	double keyPosition_ = 0.5;
	double currentStrikePosition_ = 0.22;
	double contactWidthMetres_ = 0.009;
	double tineModalMass_ = 1.0e-4;
	double notePitch_{};
	double keyPickupSensitivity_ = 1.0;
	double pickupExcursionScale_ = 1.0;
	double pickupOutput_{};
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
	double hammerIncomingVelocity_ = 0.004;
	double hammerInverseMass_ =
		1.0 / ElectricPianoPublishedMechanicalData::HammerMassKg;
	double contactStiffness_{};
	double contactExponent_ =
		ElectricPianoPublishedMechanicalData::ContactExponent;
	double contactLoss_{};
	double contactAge_{};
	double controlSmoothingCoefficient_ = 0.0035;
	double decaySmoothingCoefficient_ = 0.00012;
	int modeCoefficientUpdatePeriod_ = 48;
	int modeCoefficientUpdateCountdown_{};
	ElectricPianoControls smoothedControls_{};
	double energy_{};
	std::array<std::array<Mode, ElectricPianoModeCount>, PickupOversamplingFactor>
		pickupModeFrames_{};
	std::array<double, ElectricPianoModeCount> modulationPhasorReal_{};
	std::array<double, ElectricPianoModeCount> modulationPhasorImaginary_{};
	double previousOversampledPhaseModulation_{};
	double previousHostExponentialPitch_{};
	double previousHostLinearFrequencyRatio_{};
	double previousHostPhaseModulation_{};
	double modulationCrossfade_{};
	int modulationPhasorNormalizationCountdown_ = 4096;
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
	bool modulationPathActive_{};
	bool modulationControlsPrimed_{};
	bool exponentialFmPathActive_{};
	bool linearFmPathActive_{};
	bool phaseModulationPathActive_{};
	bool exponentialFmInterpolatorPrimed_{};
	bool linearFmInterpolatorPrimed_{};
	bool phaseModulationInterpolatorPrimed_{};
	std::unique_ptr<X4Resampler_Order7> verticalPositionInterpolator_;
	std::unique_ptr<X4Resampler_Order7> verticalVelocityInterpolator_;
	std::unique_ptr<X4Resampler_Order7> horizontalPositionInterpolator_;
	std::unique_ptr<X4Resampler_Order7> horizontalVelocityInterpolator_;
	std::unique_ptr<X4Resampler_Order7> exponentialFmInterpolator_;
	std::unique_ptr<X4Resampler_Order7> linearFmInterpolator_;
	std::unique_ptr<X4Resampler_Order7> phaseModulationInterpolator_;
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
		cachedBass_ = -1.0;
		cachedTreble_ = -1.0;
		cachedVibratoSpeed_ = -1.0;
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
		if (smoothedVibratoSpeed_ != cachedVibratoSpeed_)
		{
			vibratoFrequency_ = 1.5 * std::pow(8.0, smoothedVibratoSpeed_);
			cachedVibratoSpeed_ = smoothedVibratoSpeed_;
		}
		if (smoothedBass_ != cachedBass_)
		{
			bassPotPosition_ = BassPotPositionForControl(smoothedBass_);
			cachedBass_ = smoothedBass_;
		}
		if (smoothedTreble_ != cachedTreble_)
		{
			treblePotPosition_ = TreblePotPositionForControl(smoothedTreble_);
			cachedTreble_ = smoothedTreble_;
		}
		// Interpolate before Figure 11-8. The complete preamp -> tone -> Figure
		// 11-9 path runs at 2x; only controls and the slow supply remain at host
		// rate. The pickup/voice model has its own independent 4x path.
		const auto inputs = inputInterpolator_->Upsample(input);
		Eigen::Array<double, X2Resampler_Order7::ResamplingFactor, 1> leftOutputs;
		Eigen::Array<double, X2Resampler_Order7::ResamplingFactor, 1> rightOutputs;
		const double compensatedPreampDenominator =
			driveGain_ * NominalFullPreampGain;
		const double pan = smoothedVibratoIntensity_ * vibratoLamp_;
		const double leftDrive = stereoPowerActive ?
			std::sqrt(std::max(0.0, 1.0 + pan)) : 1.0;
		const double rightDrive = stereoPowerActive ?
			std::sqrt(std::max(0.0, 1.0 - pan)) : 1.0;
		const double modelScale = makeupGain_ /
			(PowerClosedLoopGain * InputVoltsPerModelUnit);
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
				bassPotPosition_, treblePotPosition_).voltage;
			// Figure 11-8 places its 100 kOhm volume control before the channel
			// output buffers and Figure 11-9. Use that physical node for Drive's
			// reciprocal gain. Putting the compensation after Figure 11-9 preserved
			// the output level but still drove chord peaks into the power rails.
			const double preampOutput = toneOutput /
				compensatedPreampDenominator;

			// In Figure 11-8 the optical network routes the preamp signal to two
			// outputs, each of which feeds its own Figure 11-9 power module. Keep
			// those power stages independent so their crossover and overload history
			// follow the channel currents rather than a post-distortion pan law.
			const double leftVoltage = ProcessPowerModule(preampOutput * leftDrive,
				powerChannels_[0], supplyRail_);
			double rightVoltage = leftVoltage;
			if (stereoPowerActive)
			{
				rightVoltage = ProcessPowerModule(preampOutput * rightDrive,
					powerChannels_[1], supplyRail_);
			}
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

		vibratoPhase_ += vibratoFrequency_ / sampleRate_;
		if (vibratoPhase_ >= 1.0)
			vibratoPhase_ -= 1.0;
		// Keep the physical lamp/LDR state running while Intensity is down. This
		// costs one slow nonlinear evaluation but makes enabling Vibrato resume the
		// continuously evolving oscillator instead of restarting its lamp envelope.
		const double sine = tfdsp::SinTwoPi(vibratoPhase_);
		const double lampTarget = tfdsp::TanhPade76(1.6 * sine) /
			tfdsp::TanhPade76(1.6);
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
	double cachedBass_ = -1.0;
	double cachedTreble_ = -1.0;
	double cachedVibratoSpeed_ = -1.0;
	double smoothedDrive_{};
	double driveSmoothingCoefficient_ = 0.0035;
	double controlSmoothingCoefficient_ = 0.002;
	double vibratoLampCoefficient_ = 0.002;
	double supplyDischargeCoefficient_ = 0.0069;
	double supplyRechargeCoefficient_ = 0.00038;
	double driveGain_ = 1.0;
	double makeupGain_ = NominalGain;
	double bassPotPosition_ = 0.5;
	double treblePotPosition_ = 0.5;
	double vibratoFrequency_ = 3.0;
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
