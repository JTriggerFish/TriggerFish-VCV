#pragma once

#include <algorithm>
#include <cmath>

#include "tfdsp/minblep.hpp"

namespace tfdsp
{

namespace oscillator_detail
{
	// Compact integrated-polyBLEP residuals. Applied to a discontinuity in
	// first derivative, these form a two-sample polynomial BLAMP correction.
	// This polynomial form is also used by Mutable Instruments Plaits.
	inline double NextPolyBlampSample(double elapsed)
	{
		const double half = 0.5 * std::clamp(elapsed, 0.0, 1.0);
		const double square = half * half;
		return 0.1875 - half + 1.5 * square - square * square;
	}

	inline double ThisPolyBlampSample(double elapsed)
	{
		return NextPolyBlampSample(1.0 - elapsed);
	}
}

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

/** Bipolar pulse oscillator with fractional PWM and hard sync.
 *
 * The raw pulse is +1 while phase is below duty cycle and -1 otherwise.
 * Phase wraps, moving-duty comparator crossings, and hard-sync resets are
 * reconstructed with minimum-phase band-limited steps at their exact
 * sub-sample positions. Duty cycle is interpolated linearly over each call,
 * matching a comparator threshold driven by a reconstructed control signal.
 */
template<int MinBlepZeroCrossings = 8, int MinBlepTableOversampling = 32>
class BandlimitedPulseOscillator
{
public:
	static constexpr double MinimumDutyCycle = 1.0e-6;
	static constexpr double MaximumDutyCycle = 1.0 - MinimumDutyCycle;

	void Reset()
	{
		_phase = 0.0;
		_dutyCycle = 0.5;
		_dutyInitialized = false;
		_discontinuityBlep.Reset();
	}

	double Step(double phaseIncrement, double dutyCycle = 0.5,
		double syncPosition = -1.0)
	{
		if (!std::isfinite(phaseIncrement) || !std::isfinite(dutyCycle))
		{
			Reset();
			return 0.0;
		}

		const double nextDuty = ClampDutyCycle(dutyCycle);
		if (!_dutyInitialized)
		{
			_dutyCycle = nextDuty;
			_dutyInitialized = true;
		}

		if (syncPosition >= 0.0 && syncPosition <= 1.0)
			AdvanceWithSync(phaseIncrement, nextDuty, syncPosition);
		else
			AdvanceContinuous(phaseIncrement, _dutyCycle, nextDuty, 0.0, 1.0);

		_dutyCycle = nextDuty;
		const double pulse = RawPulse(_phase, _dutyCycle) +
			_discontinuityBlep.Process();
		if (std::isfinite(pulse))
			return pulse;
		Reset();
		return 0.0;
	}

	double Phase() const { return _phase; }
	double DutyCycle() const { return _dutyCycle; }

private:
	MinBlepGenerator<MinBlepZeroCrossings,
		MinBlepTableOversampling, double> _discontinuityBlep;
	double _phase{};
	double _dutyCycle{0.5};
	bool _dutyInitialized{};

	static double ClampDutyCycle(double dutyCycle)
	{
		return std::clamp(dutyCycle, MinimumDutyCycle, MaximumDutyCycle);
	}

	static double WrapPhase(double phase)
	{
		return phase - std::floor(phase);
	}

	static double RawPulse(double phase, double dutyCycle)
	{
		return phase < dutyCycle ? 1.0 : -1.0;
	}

	void InsertStep(double framePosition, double magnitude)
	{
		if (magnitude != 0.0)
			_discontinuityBlep.InsertDiscontinuity(framePosition - 1.0,
				magnitude);
	}

	template<typename Callback>
	static void ForEachIntegerCrossing(double start, double end,
		Callback&& callback)
	{
		if (end > start)
		{
			const auto first = static_cast<long long>(std::floor(start)) + 1;
			const auto last = static_cast<long long>(std::floor(end));
			for (long long integer = first; integer <= last; ++integer)
			{
				const double fraction = (static_cast<double>(integer) - start) /
					(end - start);
				if (fraction > 0.0 && fraction <= 1.0)
					callback(fraction);
			}
		}
		else if (end < start)
		{
			const auto first = static_cast<long long>(std::ceil(start)) - 1;
			const auto last = static_cast<long long>(std::ceil(end));
			for (long long integer = first; integer >= last; --integer)
			{
				const double fraction = (static_cast<double>(integer) - start) /
					(end - start);
				if (fraction > 0.0 && fraction <= 1.0)
					callback(fraction);
			}
		}
	}

	void AdvanceContinuous(double phaseIncrement, double startDuty,
		double endDuty, double startPosition, double endPosition)
	{
		if (!(endPosition > startPosition))
			return;

		const double unwrappedStart = _phase;
		const double unwrappedEnd = unwrappedStart + phaseIncrement;
		const double comparatorStart = unwrappedStart - startDuty;
		const double comparatorEnd = unwrappedEnd - endDuty;

		// Each integer crossing of phase-duty is one comparator transition.
		// Increasing phase-duty changes +1 to -1; decreasing does the reverse.
		const double comparatorStep = comparatorEnd > comparatorStart ?
			-2.0 : 2.0;
		ForEachIntegerCrossing(comparatorStart, comparatorEnd,
			[&](double fraction)
			{
				InsertStep(startPosition + fraction *
					(endPosition - startPosition), comparatorStep);
			});

		// A forward phase wrap changes the pulse from low to high. A reverse
		// wrap performs the opposite transition.
		const double wrapStep = phaseIncrement > 0.0 ? 2.0 : -2.0;
		ForEachIntegerCrossing(unwrappedStart, unwrappedEnd,
			[&](double fraction)
			{
				InsertStep(startPosition + fraction *
					(endPosition - startPosition), wrapStep);
			});

		_phase = WrapPhase(unwrappedEnd);
	}

	void AdvanceWithSync(double phaseIncrement, double nextDuty,
		double syncPosition)
	{
		const double syncDuty = _dutyCycle + syncPosition *
			(nextDuty - _dutyCycle);
		AdvanceContinuous(phaseIncrement * syncPosition, _dutyCycle,
			syncDuty, 0.0, syncPosition);

		const double pulseBefore = RawPulse(_phase, syncDuty);
		const double resetPhase = phaseIncrement < 0.0 ?
			std::nextafter(1.0, 0.0) : 0.0;
		const double pulseAfter = RawPulse(resetPhase, syncDuty);
		InsertStep(syncPosition, pulseAfter - pulseBefore);
		_phase = resetPhase;

		AdvanceContinuous(phaseIncrement * (1.0 - syncPosition), syncDuty,
			nextDuty, syncPosition, 1.0);
	}
};

/** Bipolar triangle oscillator corrected at its fractional corner positions.
 *
 * The waveform is continuous, but its first derivative changes abruptly at
 * phase 0 and 0.5. A compact polynomial BLAMP corrects those slope changes.
 * The correction is causal with one sample of latency; OutputPhase() exposes
 * the correspondingly delayed phase so another waveform can remain aligned.
 */
class BandlimitedTriangleOscillator
{
public:
	void Reset(double phase = 0.0)
	{
		_phase = WrapPhase(std::isfinite(phase) ? phase : 0.0);
		_outputPhase = _phase;
		_nextSample = RawTriangle(_phase);
	}

	double Step(double phaseIncrement)
	{
		if (!std::isfinite(phaseIncrement))
		{
			Reset();
			return 0.0;
		}

		double output = _nextSample;
		_nextSample = 0.0;
		_outputPhase = _phase;

		const double start = _phase;
		const double end = start + phaseIncrement;
		ForEachCorner(start, end, [&](double corner, bool integerCorner)
		{
			const double fraction = (corner - start) / phaseIncrement;
			const double elapsed = 1.0 - fraction;
			double derivativeJump;
			if (phaseIncrement > 0.0)
				derivativeJump = integerCorner ?
					8.0 * phaseIncrement : -8.0 * phaseIncrement;
			else
				derivativeJump = integerCorner ?
					-8.0 * phaseIncrement : 8.0 * phaseIncrement;

			output += oscillator_detail::ThisPolyBlampSample(elapsed) *
				derivativeJump;
			_nextSample += oscillator_detail::NextPolyBlampSample(elapsed) *
				derivativeJump;
		});

		_phase = WrapPhase(end);
		_nextSample += RawTriangle(_phase);
		if (std::isfinite(output))
			return output;
		Reset();
		return 0.0;
	}

	double Phase() const { return _phase; }
	double OutputPhase() const { return _outputPhase; }

	static double RawTriangle(double phase)
	{
		const double wrapped = WrapPhase(phase);
		return 1.0 - 4.0 * std::abs(wrapped - 0.5);
	}

private:
	double _phase{};
	double _outputPhase{};
	double _nextSample{-1.0};

	static double WrapPhase(double phase)
	{
		return phase - std::floor(phase);
	}

	template<typename Callback>
	static void ForEachCorner(double start, double end, Callback&& callback)
	{
		if (end > start)
		{
			long long cornerIndex = static_cast<long long>(
				std::floor(2.0 * start)) + 1;
			const long long last = static_cast<long long>(std::floor(2.0 * end));
			for (; cornerIndex <= last; ++cornerIndex)
			{
				const double corner = 0.5 * static_cast<double>(cornerIndex);
				if (corner > start && corner <= end)
					callback(corner, cornerIndex % 2 == 0);
			}
		}
		else if (end < start)
		{
			long long cornerIndex = static_cast<long long>(
				std::ceil(2.0 * start)) - 1;
			const long long last = static_cast<long long>(std::ceil(2.0 * end));
			for (; cornerIndex >= last; --cornerIndex)
			{
				const double corner = 0.5 * static_cast<double>(cornerIndex);
				if (corner < start && corner >= end)
					callback(corner, cornerIndex % 2 == 0);
			}
		}
	}
};

} // namespace tfdsp
