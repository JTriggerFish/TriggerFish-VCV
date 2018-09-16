#pragma once

#include <cmath>
#include "../dsp/nonlinear.hpp"
#include <array>

/**
 * References:
 *     [1] : REDUCING THE ALIASING OF NONLINEAR WAVESHAPING USING CONTINUOUS-TIME CONVOLUTION, Julian Parker et al, DAFX 16
 *     [2] : Geometric intergration using discrete gradients, R. McLahlan, G.R.W Quispel, Nicolas Robidoux, 1998
 */

using namespace DiscreteGradient2;

/**
 * \brief Model of 1 pole OTA based integrator / low pass filter
 */
class OTA1PoleIntegrator
{
private:
	//Previous value of internal state
	double _u1{};
	//Previous value of input
	double _x1{};

	//solve for u[n] numerically : f(u) = g* Grad2[tanh](u, u1) + u - u1 - x + x1 = 0
	static double func(const double u, const double u1, const double x, const double x1, const double g)
	{
		return g * Tanh<double>::Value(u, u1) + u - u1 - x + x1;
	}

public:
	OTA1PoleIntegrator() = default;
	~OTA1PoleIntegrator() = default;
	OTA1PoleIntegrator(const OTA1PoleIntegrator&) = delete;
	OTA1PoleIntegrator& operator=(const OTA1PoleIntegrator&) = delete;

	static double constexpr DefaultRolloff{ 10000.0 };

private:
	/**
	 * \brief Compute the initial 2 guesses of the secant method for each of the models
	 * \tparam N number of independent models
	 * \param x input for each model, expected between -10 and 10
	 * \param g normalised cutoff gain for each model, must be prewarped: g = w^~_c * T  = 2* tan( wc *T / 2 )
	 * \param u1 previous state for each model
	 * \param x1 previous input for each model
	 * \param u interleaved state guesses ( 2 guesses ) for each model
	 * \param f interleaved function values for the 2 guesses of each model
	 */
	template<std::size_t N>
	static void InitGuesses(const std::array<double, N>& x,
		const std::array<double, N>& g,
		const std::array<double, N>& u1,
		const std::array<double, N>& x1,
		std::array<double, 2*N>& u,
		std::array<double, 2*N>& f)
	{
		for (std::size_t i = 0; i < N; ++i)
		{
			//initial guess linearised around u ~ 0
			auto guess = (x[i] - x1[i] + u1[i] * (1.0 - 0.5*g[i])) / (1 + 0.5*g[i]);
			u[2 * i] = guess;
			//u[2 * i + 1] = guess + 1.0e-6;
			u[2 * i + 1] = u1[i];
		}
		for(std::size_t i=0; i < N; ++i)
		{
			if (std::abs(u[2 * i + 1] - u[2 * i]) < 1.0e-8)
				u[2 * i + 1] += 1.0e-8;
		}
		for(std::size_t i=0; i < 2*N; ++i)
		{
			f[i] = func(u[i], u1[i / 2], x[i / 2], x1[i / 2], g[i/2]);
		}
	}

	/**
	 * \brief Compute 1 step of the secant method for each of the models
	 * \tparam N number of independent models
	 * \param x input for each model, expected between -10 and 10
	 * \param g normalised cutoff gain for each model, must be prewarped: g = w^~_c * T  = 2* tan( wc *T / 2 )
	 * \param u1 previous state for each model
	 * \param x1 previous input for each model
	 * \param u interleaved state guesses ( 2 guesses ) for each model
	 * \param f interleaved function values ( 2 ) for the guesses of each model
	 */
	template<std::size_t N>
	static void SecantIteration(const std::array<double, N>& x,
		const std::array<double, N>& g,
		const std::array<double, N>& u1,
		const std::array<double, N>& x1,
		std::array<double, 2*N>& u,
		std::array<double, 2*N>& f)
	{
		for (std::size_t i = 0; i < N; ++i)
		{
			const auto u_cur = u[2 * i];
			u[2 * i] = u_cur - f[2 * i] * (u_cur - u[2 * i + 1]) / (f[2 * i] - f[2 * i + 1]);
			f[2 * i + 1] = f[2 * i];
			f[2 * i] = func(u[2 * i], u1[i], x[i], x1[i], g[i]);
			u[2 * i + 1] = u_cur;
		}
	}

	double SolveSecantAndUpdateState(const double x, const double g, std::array<double, 2>& u, std::array<double, 2>& f)
	{
		//Secant method
		while (true)
		{
			//Note:we know input x should be in a fairly narrow range so a fixed constant stopping criterion seems
			//reasonable here, this should correspond to about 120dB accuracy
			if (std::fabs(f[0]) < 1.0e-6 || std::fabs(f[0] - f[1]) < 1.0e-12)
				break;
			SecantIteration<1>({{ x }}, {{ g }}, {{ _u1 }}, {{ _x1 }}, u, f);
		}

		_u1 = u[0];
		_x1 = x;

		return x - u[0];
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
		//u[n] - u[n-1] = x[n]-x[n-1] - g * Grad2[tanh](u[n], u[n-1])
		//u[n] = x[n] - y[n]
		// u[n-1] = _u1
		// x[n-1] = _x1
		//solve for u[n] numerically : f(u) = g* Grad2[tanh](u, u1) + u - u1 - x + x1 = 0

		//[0] : current guess
		//[1] : previous guess
		std::array<double, 2> u;
		std::array<double, 2> f;

		//Note: Newton method is much slower
		InitGuesses<1>({{ x }}, {{ g }}, {{ _u1 }}, {{ _x1 }}, u, f);
		SecantIteration<1>({{ x }}, {{ g }}, {{ _u1 }}, {{ _x1 }}, u, f);

		return SolveSecantAndUpdateState(x, g, u, f);
	}
	/**
	 * \brief Process one sample of the discretized version of the system, but for 2 models and inputs :
	 * dy/dt = w_c * tanh(x-y)
	 * \param models the two models
	 * \param x inputs / outputs, expected between -10 and 10
	 * \param g normalised cutoff gains, must be prewarped: g = w^~_c * T  = 2* tan( wc *T / 2 ) 
	 * \return filtered values in place of input
	 */
	static void StepDual (std::array<OTA1PoleIntegrator,2>& models,  std::array<double,2>& x, const std::array<double,2>& g)
	{
		std::array<double, 4> uAudioCv;
		std::array<double, 4> fAudioCv;
		std::array<double, 2> u1{ {models[0]._u1, models[1]._u1} };
		std::array<double, 2> x1{ {models[0]._x1, models[1]._x1} };

		//Interleave the 2 guesses for the secant method for the audio and cv to
		//make the compiler vectorize:
		InitGuesses(x, g, u1, x1, uAudioCv, fAudioCv);

		//Similarly, Do one step of secant method in parallel
		SecantIteration(x, g, u1, x1, uAudioCv, fAudioCv);

		//Do rest of the secant method, compiler probably can't unroll this
		//due to the branches in each inner loop 
		for (std::size_t j = 0; j < 2; ++j)
		{
			std::array<double, 2> u;
			std::array<double, 2> f;
			u[0] = uAudioCv[2 * j];
			u[1] = uAudioCv[2 * j + 1];
			f[0] = fAudioCv[2 * j];
			f[1] = fAudioCv[2 * j + 1];
			x[j] = models[j].SolveSecantAndUpdateState(x[j], g[j], u, f);
		}
	}
};
