#pragma once

#include <algorithm>
#include <cmath>

#include "tfdsp/minblep.hpp"

namespace tfdsp
{

struct OversampledEvent
{
	int segment{-1};
	double position{-1.0};

	explicit operator bool() const { return segment >= 0; }
};

template<int OversamplingFactor>
OversampledEvent MapEventToOversampledFrame(double framePosition)
{
	static_assert(OversamplingFactor >= 1,
		"An oversampled frame needs at least one segment");
	if (!(framePosition >= 0.0 && framePosition <= 1.0))
		return {};
	const double scaled = framePosition * OversamplingFactor;
	const int segment = std::min(static_cast<int>(scaled),
		OversamplingFactor - 1);
	return {segment, std::clamp(scaled - segment, 0.0, 1.0)};
}

/** Through-zero saw oscillator with fractional hard sync.
 *
 * Ordinary phase wraps use a short polyBLEP. A hard reset can have any step
 * height, so it and any wrap in the same sample are reconstructed with a
 * minimum-phase band-limited step at their exact sub-sample positions.
 */
template<int MinBlepZeroCrossings = 8, int MinBlepTableOversampling = 32>
class BandlimitedSawOscillator
{
public:
	void Reset()
	{
		_phase = 0.0;
		_discontinuityBlep.Reset();
	}

	double Step(double phaseIncrement, double syncPosition = -1.0)
	{
		if (!std::isfinite(phaseIncrement))
			return 0.0;
		const double antiAliasIncrement = AdvanceWithSync(
			phaseIncrement, syncPosition);
		const double saw = 2.0 * _phase - 1.0 -
			SignedPolyBlep(_phase, antiAliasIncrement) +
			_discontinuityBlep.Process();
		if (std::isfinite(saw))
			return saw;
		Reset();
		return 0.0;
	}

	double Phase() const { return _phase; }

private:
	MinBlepGenerator<MinBlepZeroCrossings,
		MinBlepTableOversampling, double> _discontinuityBlep;
	double _phase{};

	static double WrapPhase(double phase)
	{
		return phase - std::floor(phase);
	}

	static double PolyBlep(double phase, double increment)
	{
		if (!(increment > 0.0) || increment >= 1.0)
			return 0.0;
		if (phase < increment)
		{
			const double x = phase / increment;
			return x + x - x * x - 1.0;
		}
		if (phase > 1.0 - increment)
		{
			const double x = (phase - 1.0) / increment;
			return x * x + x + x + 1.0;
		}
		return 0.0;
	}

	static double SignedPolyBlep(double phase, double increment)
	{
		if (increment > 0.0)
			return PolyBlep(phase, increment);
		if (increment < 0.0)
			return -PolyBlep(1.0 - phase, -increment);
		return 0.0;
	}

	void AdvanceBandlimited(double increment, double startPosition,
		double endPosition)
	{
		const double startPhase = _phase;
		const double endPhase = startPhase + increment;
		if (increment > 0.0 && endPhase >= 1.0)
		{
			const double fraction = (1.0 - startPhase) / increment;
			const double eventPosition = startPosition + fraction *
				(endPosition - startPosition);
			_discontinuityBlep.InsertDiscontinuity(eventPosition - 1.0, -2.0);
		}
		else if (increment < 0.0 && endPhase < 0.0)
		{
			const double fraction = -startPhase / increment;
			const double eventPosition = startPosition + fraction *
				(endPosition - startPosition);
			_discontinuityBlep.InsertDiscontinuity(eventPosition - 1.0, 2.0);
		}
		_phase = WrapPhase(endPhase);
	}

	double AdvanceWithSync(double increment, double syncPosition)
	{
		if (!(syncPosition >= 0.0 && syncPosition <= 1.0))
		{
			_phase = WrapPhase(_phase + increment);
			return increment;
		}

		AdvanceBandlimited(increment * syncPosition, 0.0, syncPosition);
		const double resetPhase = increment < 0.0 ?
			std::nextafter(1.0, 0.0) : 0.0;
		const double discontinuity = 2.0 * (resetPhase - _phase);
		_discontinuityBlep.InsertDiscontinuity(
			syncPosition - 1.0, discontinuity);
		const double remainingIncrement = increment * (1.0 - syncPosition);
		_phase = resetPhase;
		AdvanceBandlimited(remainingIncrement, syncPosition, 1.0);
		return 0.0;
	}
};

} // namespace tfdsp
