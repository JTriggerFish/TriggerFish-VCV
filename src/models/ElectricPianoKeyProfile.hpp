#pragma once

#include "ElectricPianoCommon.hpp"

namespace tfdsp {
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


} // namespace tfdsp
