#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace tfdsp
{
	constexpr int MaximumUnisonVoices = 4;
	constexpr int MaximumStackedOscillatorVoices = 16;
	constexpr double MaximumUnisonSpreadCents = 50.0;

	/** Nonlinear musical spread with fine control near unison and a one-semitone
	 * maximum excursion for the outer voices.
	 */
	inline double UnisonSpreadCents(double control)
	{
		const double position = std::clamp(control, 0.0, 1.0);
		return MaximumUnisonSpreadCents *
			(0.05 * position + 0.95 * position * position * position);
	}

	inline double UnisonSpreadControlForCents(double cents)
	{
		const double target = std::clamp(
			cents, 0.0, MaximumUnisonSpreadCents);
		double lower = 0.0;
		double upper = 1.0;
		for (int iteration = 0; iteration < 32; ++iteration)
		{
			const double middle = 0.5 * (lower + upper);
			if (UnisonSpreadCents(middle) < target)
				lower = middle;
			else
				upper = middle;
		}
		return 0.5 * (lower + upper);
	}

	/** Pitch positions have zero mean. The irregular four-voice spacing avoids
	 * concentrating all pairwise beating at one rate.
	 */
	inline std::array<double, MaximumUnisonVoices> UnisonPitchPositions(int voices)
	{
		switch (std::clamp(voices, 1, MaximumUnisonVoices))
		{
		case 2:
			return {{-1.0, 1.0, 0.0, 0.0}};
		case 3:
			return {{-1.0, 0.0, 1.0, 0.0}};
		case 4:
			return {{-1.0, -0.18, 0.20, 0.98}};
		default:
			return {{0.0, 0.0, 0.0, 0.0}};
		}
	}

	inline double UnisonOutputGain(int voices)
	{
		return 1.0 / std::sqrt(static_cast<double>(
			std::clamp(voices, 1, MaximumUnisonVoices)));
	}

	/** Symmetric detune positions for a larger oscillator stack.
	 *
	 * The power curve keeps several oscillators near the tuning centre while the
	 * outer pair reaches the full Spread setting. Every layout has exactly zero
	 * mean, including even voice counts.
	 */
	inline std::array<double, MaximumStackedOscillatorVoices>
	StackedOscillatorPitchPositions(int voices)
	{
		std::array<double, MaximumStackedOscillatorVoices> positions{};
		const int count = std::clamp(voices, 1,
			MaximumStackedOscillatorVoices);
		if (count == 1)
			return positions;
		for (int voice = 0; voice < count; ++voice)
		{
			const double linear = -1.0 + 2.0 * voice / (count - 1.0);
			positions[voice] = std::copysign(
				std::pow(std::abs(linear), 1.6), linear);
		}
		return positions;
	}

	/** Quasi-random stereo positions independent of the pitch ordering. */
	inline std::array<double, MaximumStackedOscillatorVoices>
	StackedOscillatorPanPositions(int voices)
	{
		std::array<double, MaximumStackedOscillatorVoices> positions{};
		const int count = std::clamp(voices, 1,
			MaximumStackedOscillatorVoices);
		if (count == 1)
			return positions;

		constexpr double GoldenConjugate = 0.6180339887498948482;
		double mean = 0.0;
		for (int voice = 0; voice < count; ++voice)
		{
			const double wrapped = std::fmod(
				(voice + 0.5) * GoldenConjugate, 1.0);
			positions[voice] = 2.0 * wrapped - 1.0;
			mean += positions[voice];
		}
		mean /= count;
		double maximum = 0.0;
		for (int voice = 0; voice < count; ++voice)
		{
			positions[voice] -= mean;
			maximum = std::max(maximum, std::abs(positions[voice]));
		}
		if (maximum > 0.0)
			for (int voice = 0; voice < count; ++voice)
				positions[voice] /= maximum;
		return positions;
	}

	/** Fixed zero-mean oscillator calibration errors used by Tracking. */
	inline std::array<double, MaximumStackedOscillatorVoices>
	StackedOscillatorTrackingPositions(int voices)
	{
		std::array<double, MaximumStackedOscillatorVoices> positions{};
		const int count = std::clamp(voices, 1,
			MaximumStackedOscillatorVoices);
		if (count == 1)
			return positions;

		constexpr double IrrationalStep = 0.4142135623730950488;
		double mean = 0.0;
		for (int voice = 0; voice < count; ++voice)
		{
			const double wrapped = std::fmod(
				(voice + 0.25) * IrrationalStep, 1.0);
			positions[voice] = 2.0 * wrapped - 1.0;
			mean += positions[voice];
		}
		mean /= count;
		double maximum = 0.0;
		for (int voice = 0; voice < count; ++voice)
		{
			positions[voice] -= mean;
			maximum = std::max(maximum, std::abs(positions[voice]));
		}
		if (maximum > 0.0)
			for (int voice = 0; voice < count; ++voice)
				positions[voice] /= maximum;
		return positions;
	}
}
