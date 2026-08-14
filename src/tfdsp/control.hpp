#pragma once

#include <algorithm>
#include <cmath>
#include <random>

#include "filters.hpp"

namespace tfdsp
{
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
}
