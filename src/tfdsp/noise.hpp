#pragma once
#include <algorithm>
#include <array>
#include <random>
#include "filters.hpp"


/**
 * References:
 *  -Pinking filter for pink noise : https://ccrma.stanford.edu/~jos/sasp/Example_Synthesis_1_F_Noise.html
 */

namespace tfdsp 
{
	class WhiteNoiseSource
	{
	private:
		std::random_device _seed{};
		std::minstd_rand _rng;
		std::normal_distribution<float> _gaussian{ 0.0, 1.0 };
	public:
		WhiteNoiseSource() : _rng(_seed())
		{
		}
		float Step()
		{
			return _gaussian(_rng);
		}
		void Reset()
		{
			_gaussian.reset();
		}
	};

	class PinkNoiseSource
	{
	private:
		WhiteNoiseSource _white{};
		FirstOrderLowPassZdf<float> _filter;
		std::array<float, 4> _x{};
		std::array<float, 4> _y{};
		std::array<float, 4> _a{{ 1.0f, -2.494956002f, 2.017265875f, -0.522189400f }};
		std::array<float, 4> _b{{ 0.049922035, -0.095993537, 0.050612699, -0.004408786 }};
	public:
		PinkNoiseSource() {}

		float Filter3dbPerOctave(float x)
		{
			_x[0] = x;
			auto y = _b[0] * _x[0];
			for(std::size_t i=1; i < 4; ++i)
				y += _b[i] * _x[i] - _a[i] * _y[i];

			_y[0] = y;
			for(std::size_t i=3; i > 0; --i)
			{
				_x[i] = _x[i - 1];
				_y[i] = _y[i - 1];
			}
			return y;
		}

		float Step()
		{
			auto x = _white.Step();
			return Filter3dbPerOctave(x);
		}
		void Reset()
		{
			_white.Reset();
			_filter.Reset();
			_x.fill(0.0f);
			_y.fill(0.0f);
		}
	};
	class detune
	{
		public:
		/*
		Return a detuned v/oct value such that after being converted to a frequency
		by f(vOct) = f0 * pow(2, vOct), f( output ) = f(vOct) + det
		*/
		static double linear(double vOct, double det, double f0 = 261.63)
		{
			if (!std::isfinite(vOct) || !std::isfinite(det) ||
				!std::isfinite(f0) || f0 <= 0.0)
				return 0.0;

			static const double ln2 = std::log(2.0);
			const double boundedVOct = std::clamp(vOct, -100.0, 100.0);
			double v = det / f0 + std::exp2(boundedVOct);
			v = std::max<double>(1.0e-8, v);
			const double result = std::log(v) / ln2;
			return std::isfinite(result) ? result : 0.0;
		}
	};

}
