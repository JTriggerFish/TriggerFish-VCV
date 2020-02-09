#pragma once
#include "../Eigen/Dense"
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

	template<typename Float, int blockSize>
	class Tanh
	{
	private:
		static constexpr Float eps = _TanhEpsilon<Float>::Value;
		static Float ValueLarge(const Float x, const Float xPrev)
		{
			return std::log(std::cosh(x) / std::cosh(xPrev)) / (x - xPrev);
		}

	public:
		static Float Value(const Float x, const Float xPrev)
		{
			return std::abs(x - xPrev) <= eps ?
				std::tanh(0.5*(x + xPrev)) :
				ValueLarge(x, xPrev);
		}
		static Eigen::Array<Float, blockSize, 1> Value(const Eigen::Array<Float, blockSize, 1>&x, const Eigen::Array<Float, blockSize,1>& xPrev)
		{
			const Eigen::Array<Float, blockSize, 1> epsV = Eigen::Array<Float, blockSize, 1>::Constant(eps);
			//Ternary operator as expressed in eigen:
			return (Eigen::abs(x - xPrev) <= epsV).select(
				Eigen::tanh(0.5*(x + xPrev)),
				Eigen::log(Eigen::cosh(x) / Eigen::cosh(xPrev)) / (x - xPrev));
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
	template<typename Float, int blockSize>
	class TanhBlock
	{
	private:
		Eigen::Array<Float, blockSize, 1>  _x1;
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		TanhBlock()
		{
			_x1 = Eigen::Array<Float, blockSize, 1>::Zero();
		}
		Eigen::Array<Float, blockSize, 1> Process(const Eigen::Array<Float, blockSize, 1>& x)
		{
			//_x1 is the z^-1 or lag(1) operator on the input
			_x1.template segment<blockSize - 1>(1) = x.template segment<blockSize - 1>(0);
			Eigen::Array<Float, blockSize, 1> y = Tanh<Float, blockSize>::Value(x, _x1);
			//Need to keep the last element for when called on the next block
			_x1(0) = x(blockSize - 1);
			return y;
		}
	};

}