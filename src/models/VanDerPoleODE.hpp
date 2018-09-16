#pragma once
#include "../Eigen/Dense"
#include "../dsp/filters.hpp"
using namespace Eigen;

namespace ode
{
	/**
	 * \brief The classic Van der Pole 2d ODE
	 */
	template<typename Scalar>
	struct VanDerPoleODE
	{
		//Damping parameter
		Scalar _mu{ Scalar(0.1) };
		//Radian frequency ( 2* pi * f)
		Scalar _w{ Scalar(200 * 2 * PI) };

		explicit VanDerPoleODE(Scalar mu, Scalar w)
			: _mu(mu), _w(w)
		{
		}
		VanDerPoleODE() {};

		/**
		 * \brief 
		 * \param y current state
		 * \param x forcing input
		 * \return Derivative
		 */
        Matrix<Scalar, 2, 1> DyDt(const Matrix<Scalar, 2, 1> &y, Scalar x) const
        {
			Matrix<Scalar, 2, 1> deriv;
			deriv << y(1),
				_mu * (1.0f - y(0)*y(0)) * y(1) * _w + _w * _w * (x - y(0));
			return deriv;
		}
        Matrix<Scalar, 2, 2> Jacobian(const Matrix<Scalar, 2, 1> &y, Scalar x) const
        {
			Matrix<Scalar, 2, 2> J;
			J << 0.0f, 1.0f,
				-2.0f*_mu*y(0)*y(1)*_w - _w * _w, _mu * (1.0f - y(0)*y(0)) *_w;
			return J;
		}
	};
}