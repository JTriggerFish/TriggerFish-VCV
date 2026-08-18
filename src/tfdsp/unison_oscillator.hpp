#pragma once

#include <algorithm>
#include <cmath>

#include "tfdsp/oscillator.hpp"

namespace tfdsp
{
	struct StackedOscillatorSample
	{
		double main{};
		double sub{};
	};

	/** One independently band-limited voice in an oscillator stack.
	 *
	 * Saw and pulse selections remain phase-aligned, and callers may also use a
	 * continuous blend. The sub oscillator runs at exactly half the parent frequency;
	 * this is the digital equivalent of a clock divider when neither oscillator
	 * is externally synchronized or frequency-modulated.
	 */
	class StackedOscillatorVoice
	{
	public:
		void Reset(double phase = 0.0)
		{
			const double finitePhase = std::isfinite(phase) ? phase : 0.0;
			const double wrapped = finitePhase - std::floor(finitePhase);
			_saw.Reset(wrapped);
			_pulse.Reset(wrapped);
			_sub.Reset(0.5 * wrapped);
			_mainMode = MainMode::Saw;
		}

		double StepMain(double phaseIncrement, double pulseWidth, double pulseMix)
		{
			const double boundedIncrement = std::clamp(phaseIncrement, 0.0, 0.45);
			const double mix = std::clamp(pulseMix, 0.0, 1.0);
			constexpr double SettledThreshold = 1.0e-4;
			if (mix <= SettledThreshold)
			{
				if (_mainMode == MainMode::Pulse)
					_saw.Reset(_pulse.Phase());
				_mainMode = MainMode::Saw;
				return _saw.Step(boundedIncrement);
			}
			if (mix >= 1.0 - SettledThreshold)
			{
				if (_mainMode == MainMode::Saw)
					_pulse.Reset(_saw.Phase());
				_mainMode = MainMode::Pulse;
				return _pulse.Step(boundedIncrement,
					std::clamp(pulseWidth, 0.05, 0.95));
			}
			if (_mainMode == MainMode::Saw)
				_pulse.Reset(_saw.Phase());
			else if (_mainMode == MainMode::Pulse)
				_saw.Reset(_pulse.Phase());
			_mainMode = MainMode::Blend;
			const double saw = _saw.Step(boundedIncrement);
			const double pulse = _pulse.Step(boundedIncrement,
				std::clamp(pulseWidth, 0.05, 0.95));
			const double result = saw + mix * (pulse - saw);
			if (std::isfinite(result))
				return result;
			Reset();
			return 0.0;
		}

		double StepSub(double phaseIncrement)
		{
			const double result = _sub.Step(0.5 * std::clamp(
				phaseIncrement, 0.0, 0.45), 0.5);
			if (std::isfinite(result))
				return result;
			Reset();
			return 0.0;
		}

		StackedOscillatorSample Step(double phaseIncrement, double pulseWidth,
			double pulseMix, bool renderSub = true)
		{
			return {
				StepMain(phaseIncrement, pulseWidth, pulseMix),
				renderSub ? StepSub(phaseIncrement) : 0.0,
			};
		}

		double Phase() const
		{
			return _mainMode == MainMode::Pulse ? _pulse.Phase() : _saw.Phase();
		}

	private:
		enum class MainMode
		{
			Saw,
			Blend,
			Pulse,
		};

		BandlimitedSawOscillator<> _saw{};
		BandlimitedPulseOscillator<> _pulse{};
		BandlimitedFixedPulseOscillator _sub{};
		MainMode _mainMode{MainMode::Saw};
	};
}
