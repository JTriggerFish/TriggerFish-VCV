#pragma once

#include "../Eigen/Dense"
#include <cmath>
#include <array>
#include "../dsp/nonlinear.hpp"

/**
 * References:
 *     [1] : REDUCING THE ALIASING OF NONLINEAR WAVESHAPING USING CONTINUOUS-TIME CONVOLUTION, Julian Parker et al, DAFX 16
 *     [2] : Geometric intergration using discrete gradients, R. McLahlan, G.R.W Quispel, Nicolas Robidoux, 1998
 */

using namespace DiscreteGradient2;

/**
 * \brief Model of 1 pole Transistor based integrator / low pass filter
 */
class Transistor1PoleIntegrator
{
private:

	//Previous value of internal state
	double _y1{};
	//Previous value of input
	double _x1{};

	//solve for y[n] numerically : f(y) = g* Grad2[tanh](y, y1) - g* phiX + y - y1 = 0
	//phiX = Grad2[tanh](x, x1)
	static double func(const double y, const double y1, const double phiX, const double g)
	{
		return g * Tanh<double, 1>::Value(y, y1) - g*phiX + y - y1;
	}

public:
	Transistor1PoleIntegrator() = default;
	~Transistor1PoleIntegrator() = default;
	Transistor1PoleIntegrator(const Transistor1PoleIntegrator&) = delete;
	Transistor1PoleIntegrator& operator=(const Transistor1PoleIntegrator&) = delete;

	static double constexpr DefaultRolloff{ 20000.0 };

private:
	/**
	 * \brief Compute the initial 2 guesses of the secant method for each of the models
	 * \tparam N number of independent models
	 * \param phiX processed input for each model
	 * \param g normalised cutoff gain for each model, must be prewarped: g = w^~_c * T  = 2* tan( wc *T / 2 )
	 * \param y1 previous state for each model
	 * \param y interleaved state guesses ( 2 guesses ) for each model
	 * \param f interleaved function values for the 2 guesses of each model
	 */
	template<int N>
	static void InitGuesses(
		const Eigen::Array<double, N, 1 >& phiX,
		const Eigen::Array<double, N, 1>& g,
		const Eigen::Array<double, N, 1>& y1,
		Eigen::Ref<Eigen::Array<double, 2*N, 1>> y,
		Eigen::Ref<Eigen::Array<double, 2*N, 1>> f)
	{
		for (int i = 0; i < N; ++i)
		{
			//initial guess linearised around u ~ 0
			double guess = (g(i) * phiX(i) + (1.0 - 0.5*g(i)) * y1(i)) / (1 + 0.5*g(i));
			y(2 * i) = guess;
			//y(2 * i + 1) = guess + 1.0e-6;
			y(2 * i + 1) = y1(i);
		}
		for(int i=0; i < N; ++i)
		{
			if (std::abs(y(2 * i + 1) - y(2 * i)) < 1.0e-8)
				y(2 * i + 1) += 1.0e-8;
		}
		for(int i=0; i < 2*N; ++i)
		{
			f(i) = func(y(i), y1(i / 2), phiX(i / 2), g(i / 2));
		}
	}

	/**
	 * \brief Compute 1 step of the secant method for each of the models
	 * \tparam N number of independent models
	 * \param phiX processed input for each model
	 * \param g normalised cutoff gain for each model, must be prewarped: g = w^~_c * T  = 2* tan( wc *T / 2 )
	 * \param y1 previous state for each model
	 * \param y interleaved state guesses ( 2 guesses ) for each model
	 * \param f interleaved function values for the 2 guesses of each model
	 */
	template<int N>
	static void SecantIteration(
		const Eigen::Array<double, N, 1 >& phiX,
		const Eigen::Array<double, N, 1>& g,
		const Eigen::Array<double, N, 1>& y1,
		Eigen::Ref<Eigen::Array<double, 2*N, 1>> y,
		Eigen::Ref<Eigen::Array<double, 2*N, 1>> f)
	{
		for (int i = 0; i < N; ++i)
		{
			const auto y_cur = y(2 * i);
			y(2 * i) = y_cur - f(2 * i) * (y_cur - y(2 * i + 1)) / (f(2 * i) - f(2 * i + 1));
			f(2 * i + 1) = f(2 * i);
			f(2 * i) = func(y(2 * i), y1(i), phiX(i), g(i));
			y(2 * i + 1) = y_cur;
		}
	}

	double SolveSecantAndUpdateState( const double x, const double phiX, const double g,
		Eigen::Ref<Eigen::Array<double, 2, 1>> y,
		Eigen::Ref<Eigen::Array<double, 2, 1>> f)
	{
		Eigen::Array<double, 1, 1> v_phiX; v_phiX << phiX;
		Eigen::Array<double, 1, 1> v_g   ; v_g	   << g;
		Eigen::Array<double, 1, 1> v_y1  ; v_y1   << _y1;
		//Secant method
		while (true)
		{
			//Note:we know input x should be in a fairly narrow range so a fixed constant stopping criterion seems
			//reasonable here, this should correspond to about 120dB accuracy
			if (std::fabs(f(0)) < 1.0e-6 || std::fabs(f(0) - f(1)) < 1.0e-12)
				break;
			SecantIteration<1>(v_phiX, v_g, v_y1, y, f);
		}

		_y1 = y(0);
		_x1 = x;

		return y(0);
	}


public:
	/**
	 * \brief Process one sample of the discretized version of the system :
	 * dy/dt = w_c * tanh(x-y)
	 * \param x input, expected between -10 and 10
	 * \param g normalised cutoff gain, must be prewarped: g = w^~_c * T  = 2* tan( wc *T / 2 ) 
	 * \return filtered value
	 */
	double Step(const double x, const double g)
	{
		//Solve the ODE discretized with a second order gradient method:
		//y[n] - y[n-1] = g * [ Grad2[tanh](x[n], x[n-1]) - Grad2[tanh](y[n], y[n-1]) ]
		// y[n-1] = _y1
		// x[n-1] = _x1
		//solve for y[n] numerically : f(y) = g* Grad2[tanh](y, y1) - g* Grad2[tanh](x, x1) + y - y1 = 0

		//[0] : current guess
		//[1] : previous guess
		Eigen::Array<double, 2, 1> u;
		Eigen::Array<double, 2, 1> f;

		const auto phiX = Tanh<double,1>::Value(x, _x1);

		//Note: Newton method is much slower than the secant method here
		Eigen::Array<double, 1, 1> v_phiX; v_phiX << phiX;
		Eigen::Array<double, 1, 1> v_g   ; v_g	   << g;
		Eigen::Array<double, 1, 1> v_y1  ; v_y1   << _y1;

		InitGuesses<1>(v_phiX, v_g, v_y1, u, f);
		SecantIteration<1>(v_phiX, v_g, v_y1, u, f);

		return SolveSecantAndUpdateState(x, phiX, g, u, f);
	}
	/**
	 * \brief Process one sample of the discretized version of the system, but for 2 models and inputs :
	 * dy/dt = w_c * tanh(x-y)
	 * \param models the two models
	 * \param x inputs / outputs, expected between -10 and 10
	 * \param g normalised cutoff gains, must be prewarped: g = w^~_c * T  = 2* tan( wc *T / 2 ) 
	 * \return filtered values in place of input
	 */
	static void StepDual (std::array<Transistor1PoleIntegrator, 2>& models,  Eigen::Ref<Eigen::Array<double,2, 1>> x, const Eigen::Array<double,2, 1>& g)
	{
		Eigen::Array<double, 4, 1> uAudioCv;
		Eigen::Array<double, 4, 1> fAudioCv;
		Eigen::Array<double, 2, 1> y1;
		y1 << models[0]._y1, models[1]._y1;
		Eigen::Array<double, 2, 1> x1;
		x1 << models[0]._x1, models[1]._x1;
		Eigen::Array<double, 2, 1> phiX;

		for (int j = 0; j < 2; ++j)
			phiX(j) = Tanh<double, 1>::Value(x(j), x1(j));

		//Interleave the 2 guesses for the secant method for the audio and cv to
		//make the compiler vectorize:
		InitGuesses<2>(phiX, g, y1, uAudioCv, fAudioCv);

		//Similarly, Do one step of secant method in parallel
		SecantIteration<2>(phiX, g, y1, uAudioCv, fAudioCv);

		//Do rest of the secant method, compiler probably can't unroll this
		//due to the branches in each inner loop 
		for (int j = 0; j < 2; ++j)
		{
			Eigen::Array<double, 2, 1> u;
			Eigen::Array<double, 2, 1> f;
			u(0) = uAudioCv(2 * j);
			u(1) = uAudioCv(2 * j + 1);
			f(0) = fAudioCv(2 * j);
			f(1) = fAudioCv(2 * j + 1);
			x(j) = models[j].SolveSecantAndUpdateState(x(j), phiX(j), g(j), u, f);
		}
	}
};
