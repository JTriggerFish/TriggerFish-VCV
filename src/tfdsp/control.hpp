#pragma once

#include <algorithm>
#include <cmath>
#include <random>

#include "filters.hpp"

namespace tfdsp
{
	class FractionalSchmittTrigger
	{
	public:
		struct Event
		{
			bool triggered{};
			double position{1.0};
		};

		void Reset()
		{
			_initialized = false;
			_high = false;
			_previous = 0.0;
		}

		Event Process(double input, double lowThreshold = 0.1,
			double highThreshold = 1.0)
		{
			if (!std::isfinite(input))
				input = 0.0;
			if (!std::isfinite(lowThreshold) || !std::isfinite(highThreshold) ||
				lowThreshold >= highThreshold)
			{
				lowThreshold = 0.1;
				highThreshold = 1.0;
			}

			if (!_initialized)
			{
				_initialized = true;
				_high = input >= highThreshold;
				_previous = input;
				return {};
			}

			Event event;
			if (!_high && input >= highThreshold)
			{
				_high = true;
				event.triggered = true;
				const double difference = input - _previous;
				event.position = difference > 1.0e-15 ?
					std::clamp((highThreshold - _previous) / difference, 0.0, 1.0) :
					1.0;
			}
			else if (_high && input <= lowThreshold)
			{
				_high = false;
			}
			_previous = input;
			return event;
		}

		bool IsHigh() const { return _high; }

	private:
		double _previous{};
		bool _initialized{};
		bool _high{};
	};

	/** A fixed-frequency sine oscillator using a complex rotation.
	 *
	 * Compared with calling sin() for every sample, this trades a tiny amount of
	 * accumulated floating-point error for a few multiplies and additions. Double
	 * precision plus periodic normalization keeps that error negligible.
	 */
	class RecursiveSineOscillator
	{
		double _sin{};
		double _cos{1.0};
		double _sinStep{};
		double _cosStep{1.0};
		unsigned int _samplesUntilNormalization{4096};

	public:
		void SetFrequency(double frequency, double sampleRate)
		{
			const double phaseStep = 2.0 * PI * frequency / sampleRate;
			_sinStep = std::sin(phaseStep);
			_cosStep = std::cos(phaseStep);
		}

		void Reset()
		{
			_sin = 0.0;
			_cos = 1.0;
			_samplesUntilNormalization = 4096;
		}

		double Step()
		{
			const double nextSin = _sin * _cosStep + _cos * _sinStep;
			const double nextCos = _cos * _cosStep - _sin * _sinStep;
			_sin = nextSin;
			_cos = nextCos;

			if (--_samplesUntilNormalization == 0)
			{
				const double inverseMagnitude = 1.0 / std::hypot(_sin, _cos);
				_sin *= inverseMagnitude;
				_cos *= inverseMagnitude;
				_samplesUntilNormalization = 4096;
			}
			return _sin;
		}
	};

	/** Apply a zero-centred drift value to a bounded control without clipping.
	 *
	 * A unit Gaussian drift reaches about 76% of the available headroom at one
	 * standard deviation and 96% at two. Depth scales that excursion from zero
	 * to the full available range on either side of the control's current value.
	 */
	inline double ApplyBoundedDrift(double value, double drift, double depth,
		double minimum = 0.0, double maximum = 1.0)
	{
		if (!std::isfinite(minimum) || !std::isfinite(maximum) ||
			!(maximum > minimum))
			return std::isfinite(value) ? value : 0.0;
		if (!std::isfinite(value) || !std::isfinite(drift) ||
			!std::isfinite(depth))
			return std::clamp(std::isfinite(value) ? value : minimum,
				minimum, maximum);

		const double normalized = std::clamp(
			(value - minimum) / (maximum - minimum), 0.0, 1.0);
		const double boundedDepth = std::clamp(depth, 0.0, 1.0);
		if (boundedDepth == 0.0)
			return std::clamp(value, minimum, maximum);
		const double shapedDrift = std::tanh(drift);
		const double headroom = shapedDrift >= 0.0 ?
			1.0 - normalized : normalized;
		const double modulated = normalized + boundedDepth *
			headroom * shapedDrift;
		return minimum + (maximum - minimum) * modulated;
	}

	/** Exact sampled Ornstein-Uhlenbeck process with interpolated control-rate updates. */
	class InterpolatedOrnsteinUhlenbeck
	{
		double _start{};
		double _target{};
		double _phase{};
		double _phaseIncrement{};
		double _decay{};
		double _innovationStdDev{};
		std::normal_distribution<double> _normal{0.0, 1.0};

	public:
		void Configure(double sampleRate, double timeConstant, double diffusion,
			double controlRate = 100.0)
		{
			const double boundedControlRate = std::clamp(controlRate, 1.0, sampleRate);
			const double controlInterval = 1.0 / boundedControlRate;
			_phaseIncrement = boundedControlRate / sampleRate;
			_decay = std::exp(-controlInterval / timeConstant);
			_innovationStdDev = diffusion * std::sqrt(
				0.5 * timeConstant * (1.0 - _decay * _decay));
		}

		/** Configure by stationary standard deviation instead of diffusion.
		 * This keeps excursion independent of the mean-reversion time constant.
		 */
		void ConfigureStationary(double sampleRate, double timeConstant,
			double stationaryStdDev = 1.0, double controlRate = 100.0)
		{
			const double boundedTimeConstant = std::isfinite(timeConstant) ?
				std::max(timeConstant, 1.0e-6) : 1.0;
			const double boundedStdDev = std::isfinite(stationaryStdDev) ?
				std::max(stationaryStdDev, 0.0) : 0.0;
			Configure(sampleRate, boundedTimeConstant,
				boundedStdDev * std::sqrt(2.0 / boundedTimeConstant), controlRate);
		}

		void Reset()
		{
			_start = 0.0;
			_target = 0.0;
			_phase = 0.0;
			_normal.reset();
		}

		template <typename Generator>
		double Step(Generator& generator)
		{
			_phase += _phaseIncrement;
			if (_phase >= 1.0)
			{
				_phase -= 1.0;
				_start = _target;
				_target = _decay * _target + _innovationStdDev * _normal(generator);
			}
			return _start + _phase * (_target - _start);
		}
	};

	/** Smooth stationary OU motion for slowly moving musical controls.
	 *
	 * A one-pole smoother with the same time constant follows the OU process.
	 * The resulting second-order motion has a continuous slope instead of the
	 * fine Brownian roughness of a first-order OU path. The sqrt(2) gain restores
	 * the stationary variance lost by cascading two equal-time-constant poles.
	 */
	class SmoothOrnsteinUhlenbeck
	{
		InterpolatedOrnsteinUhlenbeck _ou{};
		double _smoothed{};
		double _smoothingCoefficient{1.0};

	public:
		void ConfigureStationary(double sampleRate, double timeConstant,
			double stationaryStdDev = 1.0, double controlRate = 100.0)
		{
			const double boundedSampleRate = std::isfinite(sampleRate) ?
				std::max(sampleRate, 1.0) : 1.0;
			const double boundedTimeConstant = std::isfinite(timeConstant) ?
				std::max(timeConstant, 1.0e-6) : 1.0;
			_ou.ConfigureStationary(boundedSampleRate, boundedTimeConstant,
				stationaryStdDev, controlRate);
			_smoothingCoefficient = -std::expm1(
				-1.0 / (boundedSampleRate * boundedTimeConstant));
		}

		void Reset()
		{
			_ou.Reset();
			_smoothed = 0.0;
		}

		template <typename Generator>
		double Step(Generator& generator)
		{
			const double raw = _ou.Step(generator);
			_smoothed += _smoothingCoefficient * (raw - _smoothed);
			return std::sqrt(2.0) * _smoothed;
		}
	};
}
