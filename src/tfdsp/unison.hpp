#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace tfdsp
{
	constexpr int MaximumUnisonVoices = 4;
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
}
