#pragma once
#include <cmath>

/**
 * Mid point discrete gradients for integration of ODEs with nonlinear functions
 * References: Reducing
 * References:
 *     [1] : REDUCING THE ALIASING OF NONLINEAR WAVESHAPING USING CONTINUOUS-TIME CONVOLUTION, Julian Parker et al, DAFX 16
 *     [2] : Geometric intergration using discrete gradients, R. McLahlan, G.R.W Quispel, Nicolas Robidoux, 1998
 */
namespace DiscreteGradient2
{
	template<typename Float>
	struct _TanhEpsilon
	{
		static constexpr Float Value = 1.0e-12;
	};
	template<>
	struct _TanhEpsilon<float>
	{
		static constexpr float Value = 1.0e-6f;
	};
	template<>
	struct _TanhEpsilon<double>
	{
		static constexpr double Value = 1.0e-12;
	};

	template<typename Float>
	class Tanh
	{
	private:
		static constexpr Float eps = _TanhEpsilon<Float>::Value;
		static Float ValueLarge(const Float x, const Float xPrev)
		{
			return std::log(std::cosh(x)/std::cosh(xPrev)) / (x - xPrev);
		}

	public:

		static Float Value(const Float x, const Float xPrev)
		{
			return std::abs(x - xPrev) <= eps ?
				std::tanh(0.5*(x + xPrev)) :
				ValueLarge(x, xPrev);
		}
		static Float Derivative(const Float x, const Float xPrev)
		{
			if (std::abs(x - xPrev) <= eps)
			{
				auto t = std::tanh(0.5*(x + xPrev));
				return 0.5*(1.0 - t * t);
			}
			return (std::tanh(x) - ValueLarge(x, xPrev)) / (x - xPrev);
		}
	};

}