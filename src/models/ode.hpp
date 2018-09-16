#pragma once
#include "../Eigen/Dense"

using namespace Eigen;

namespace ode
{
	/**
	 * \brief Representation of current and past states for a linear multistep ODE methods
	 * like BDF.
	 * Note that the number of past states is fixed at 8 even though the max useable order of BDF is 6 which would require storing 7 elements,
	 * this is done for alignment reasons and fixed to 8 for simplicity.
	 * TODO : review if this should be reduced for low order methods performance
	 * \tparam Scalar float or double
	 * \tparam Dim state space dimension, should be <= 4
	 */
	template<typename Scalar, int Dim>
	using StateMatrix = Eigen::Matrix<Scalar, Dim, 8>;

	template<typename Scalar, int Dim>
	using StateVector = Eigen::Matrix<Scalar, Dim, 1>;

	template<typename Function, typename Scalar, int Dim, int Order >
	class SubStep
	{
	public:
		static void StepSolve(Ref<StateMatrix<Scalar, Dim>> state, int& numPastSteps, const Function& f, Scalar T, Scalar x);
	};
	template<typename Scalar>
	struct BDFNewtonConstants
	{
		static constexpr Scalar incrementEps{1.0e-12};
	};
	template<typename Function, typename Scalar, int Dim>
	void BDFNewtonSolve(Ref<StateMatrix<Scalar, Dim>> state, const Function& f, Scalar x, const StateVector<Scalar,Dim>& midTerm, Scalar rightCoeff, Scalar T, Scalar tolerance = 1.0e-7,
		int maxIterations = 10);

	/**
	 * Backward differentiation linear multistep ODE solver.
	 * Important note : in the current version this should only be used for small Dimensions
	 * i.e Dim <= 4, a dynamic allocation version should be written for larger Dimensions.
	 */
	template<typename Function, typename Scalar, int Dim, int Order >
	class BDF
	{
	private:
		Scalar _T;
		StateMatrix<Scalar, Dim> _state;
		int _numStepsKnown{};

	public:
	  BDF(const BDF &) = delete;
	  const BDF &operator=(const BDF &) = delete;
	  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	  BDF()
	  {
		  _state = StateMatrix<Scalar, Dim>::Zero();
		}
		void SetSampleRate(Scalar sampleRate)
		{
			_T = Scalar(1.0) / sampleRate;
		}
		void SetInitConditions(const StateVector<Scalar,Dim>& init)
		{
			_state = StateMatrix<Scalar,Dim>::Zero();
			_state.col(1) = init;
			//Note setting the first column has no effect unless we take the current state before solving a step
			_state.col(0) = init;
			_numStepsKnown = 1;
		}
		void Step(const Function& f, Scalar x)
		{
			SubStep<Function, Scalar, Dim, Order>::StepSolve(_state, _numStepsKnown, f, _T, x);
		}
		StateVector<Scalar,Dim> CurrentState() const
		{
			return _state.col(0);
		}
		Ref<StateMatrix<Scalar,Dim>> FullState()
		{
			return _state;
		}
	};

	//Order 1 BDF ie implicit Euler method
	template<typename Function, typename Scalar, int Dim>
	class SubStep<Function,Scalar,Dim,1>
	{
	public:
		static void StepSolve(Ref<StateMatrix<Scalar, Dim>> state, int& numPastSteps, const Function& f, Scalar T, Scalar x)
		{
			//Solve  y[n] - y[n-1] = T * ( f(y[n] , x[n]) )
			BDFNewtonSolve<Function,Scalar,Dim>(state, f,x,
				-state.col(1),
				Scalar(1),
				T);
 
			//Note : always keep one more column in case the method is used as an initialization for a 
			//higher order method
			state.col(2) = state.col(1);
			state.col(1) = state.col(0);
		}
	};
	template<typename Function, typename Scalar, int Dim>
	class SubStep<Function,Scalar,Dim,2>
	{
	public:
		static void StepSolve(Ref<StateMatrix<Scalar, Dim>> state, int& numPastSteps, const Function& f, Scalar T, Scalar x)
		{
			if(numPastSteps < 2)
			{
				SubStep<Function, Scalar, Dim, 1>::StepSolve(state, numPastSteps, f, T, x);
				numPastSteps = 2;
			}
			//Solve  y[n] - 4/3*y[n-1] + 1/3*y[n-2] = 2/3 * T * ( f(y[n] , x[n]) )

			BDFNewtonSolve<Function, Scalar, Dim>(state, f, x,
				-Scalar(4) / 3  * state.col(1)
				+ Scalar(1) / 3 * state.col(2),
				Scalar(2) / 3,
				T);

			//Note : always keep one more column in case the method is used as an initialization for a 
			//higher order method
			state.col(3) = state.col(2);
			state.col(2) = state.col(1);
			state.col(1) = state.col(0);
		}
	};
	template<typename Function, typename Scalar, int Dim>
	class SubStep<Function,Scalar,Dim,3>
	{
	public:
		static void StepSolve(Ref<StateMatrix<Scalar, Dim>> state, int& numPastSteps, const Function& f, Scalar T, Scalar x)
		{
			if(numPastSteps < 3)
			{
				SubStep<Function, Scalar, Dim, 2>::StepSolve(state, numPastSteps, f, T, x);
				numPastSteps = 3;
			}

			BDFNewtonSolve<Function, Scalar, Dim>(state, f, x,
				-Scalar(18) / 11 * state.col(1)
				+ Scalar(9) / 11 * state.col(2)
				- Scalar(2) / 11 * state.col(3),
				Scalar(6) / 11,
				T);

			//Note : always keep one more column in case the method is used as an initialization for a 
			//higher order method
			state.col(4) = state.col(3);
			state.col(3) = state.col(2);
			state.col(2) = state.col(1);
			state.col(1) = state.col(0);
		}
	};
	template<typename Function, typename Scalar, int Dim>
	class SubStep<Function,Scalar,Dim,4>
	{
	public:
		static void StepSolve(Ref<StateMatrix<Scalar, Dim>> state, int& numPastSteps, const Function& f, Scalar T, Scalar x)
		{
			if(numPastSteps < 4)
			{
				SubStep<Function, Scalar, Dim, 3>::StepSolve(state, numPastSteps, f, T, x);
				numPastSteps = 4;
			}

			BDFNewtonSolve<Function, Scalar, Dim>(state, f, x,
				-Scalar(48) / 25 * state.col(1)
				+ Scalar(36) / 25 * state.col(2)
				- Scalar(16) / 25 * state.col(3)
				+ Scalar(3) / 25 * state.col(4),
				Scalar(12) / 25,
				T);

			//Note : always keep one more column in case the method is used as an initialization for a 
			//higher order method
			state.col(5) = state.col(4);
			state.col(4) = state.col(3);
			state.col(3) = state.col(2);
			state.col(2) = state.col(1);
			state.col(1) = state.col(0);
		}
	};
	template<typename Function, typename Scalar, int Dim>
	class SubStep<Function,Scalar,Dim,5>
	{
	public:
		static void StepSolve(Ref<StateMatrix<Scalar, Dim>> state, int& numPastSteps, const Function& f, Scalar T, Scalar x)
		{
			if(numPastSteps < 5)
			{
				SubStep<Function, Scalar, Dim, 4>::StepSolve(state, numPastSteps, f, T, x);
				numPastSteps = 5;
			}

			BDFNewtonSolve<Function, Scalar, Dim>(state, f, x,
				-Scalar(300) / 137 * state.col(1)
				+ Scalar(300) / 137 * state.col(2)
				- Scalar(200) / 137 * state.col(3)
				+ Scalar(75) / 137 * state.col(4)
				- Scalar(12) / 137 * state.col(5),
				Scalar(60) / 137,
				T);

			//Note : always keep one more column in case the method is used as an initialization for a 
			//higher order method
			state.col(6) = state.col(5);
			state.col(5) = state.col(4);
			state.col(4) = state.col(3);
			state.col(3) = state.col(2);
			state.col(2) = state.col(1);
			state.col(1) = state.col(0);
		}
	};
	//Note this is the maximum stable order for BDF
	template<typename Function, typename Scalar, int Dim>
	class SubStep<Function,Scalar,Dim,6>
	{
	public:
		static void StepSolve(Ref<StateMatrix<Scalar, Dim>> state, int& numPastSteps, const Function& f, Scalar T, Scalar x)
		{
			if(numPastSteps < 6)
			{
				SubStep<Function, Scalar, Dim, 5>::StepSolve(state, numPastSteps, f, T, x);
				numPastSteps = 6;
			}

			BDFNewtonSolve<Function, Scalar, Dim>(state, f, x,
				-Scalar(360) / 147 * state.col(1)
				+ Scalar(450) / 147 * state.col(2)
				- Scalar(400) / 147 * state.col(3)
				+ Scalar(225) / 147 * state.col(4)
				- Scalar(72) / 147 * state.col(5)
				+ Scalar(10) / 147 * state.col(6),
				Scalar(60) / 147,
				T);

			state.col(7) = state.col(6);
			state.col(6) = state.col(5);
			state.col(5) = state.col(4);
			state.col(4) = state.col(3);
			state.col(3) = state.col(2);
			state.col(2) = state.col(1);
			state.col(1) = state.col(0);
		}
	};
	template<>
	struct BDFNewtonConstants<float>
	{
		static constexpr float incrementEps{ 1.0e-10f };
	};
	template<>
	struct BDFNewtonConstants<double>
	{
		static constexpr double incrementEps{ 1.0e-12 };
	};

	template<typename Function, typename Scalar, int Dim>
	void BDFNewtonSolve(Ref<StateMatrix<Scalar, Dim>> state,
		const Function& f,
		Scalar x,
		const StateVector<Scalar,Dim>& midTerm,
		Scalar rightCoeff,
		Scalar T,
		Scalar tolerance,
		int maxIterations)
	{
		//Solve Phi(Y) = Y + midTerm - T* rightCoeff * f(Y) = 0
		//where Y = the first column of state
		//and midterm is a function of the other columns

		//1 step euler for the initial guess
		StateVector<Scalar, Dim> y = state.col(1) + T * f.DyDt(state.col(1), x);
		StateVector<Scalar, Dim> phi = y + midTerm - T * rightCoeff * f.DyDt(y, x);

		int i = 0;

		while( i++ < maxIterations && phi.norm() > tolerance)
		{
			Matrix<Scalar, Dim, Dim> J = Matrix<Scalar, Dim, Dim>::Identity() - T * rightCoeff * f.Jacobian(y, x);
			//Solve J * (Y[n+1] - Y[n]) = - Phi(Y[n])
			StateVector<Scalar, Dim> increment = J.colPivHouseholderQr().solve(-phi);
			if (increment.norm() <= BDFNewtonConstants<Scalar>::incrementEps)
				break;
			//if(increment(0) != increment(0))
			//{
			//	int breakHere = 0;
			//}
			y = y + increment;
			phi = y + midTerm - T * rightCoeff * f.DyDt(y, x);
		}
		//Store solution
		state.col(0) = y;
	}
}