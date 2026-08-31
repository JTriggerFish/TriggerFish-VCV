#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <memory>
#include <stdexcept>

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


} // namespace tfdsp
