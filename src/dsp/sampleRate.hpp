#pragma once
#include "../Eigen/Dense"
#include <memory>
#include <functional>
#include "util.hpp"

/***
 * References : Digital signal processing schemes for efficient interpolation and decimation, R.A Valenzuela, A.G Constantinides, IEE 1983
 */

namespace dsp
{

	template<typename Impl, int Factor>
	class Resampler : public enable_down_cast<Impl>
	{
	private:
		using enable_down_cast<Impl>::Self;
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		static constexpr int ResamplingFactor{ Factor };
		inline Eigen::Array<double, Factor, 1> Upsample(const double x)
		{
			return Self()->_Upsample(x);
		}
		inline double Downsample(const Eigen::Array<double, Factor, 1> &x2)
		{
			return Self()->_Downsample(x2);
		}
	};

	class DummyResampler : public Resampler<DummyResampler, 1>
	{
	public:
		DummyResampler() {}
		friend class Resampler<DummyResampler, 1>;
	protected:
		inline Eigen::Array<double, 1, 1> _Upsample(const double x)
		{
			Eigen::Array<double, 1, 1> xA;
			xA << x;
			return xA;
		}
		inline double _Downsample(const Eigen::Array<double, 1, 1>& xA)
		{
			return xA(0);
		}
	};

	/// N direct stages, M delayed stages
	template<int N, int M>
	class PolyphaseIIR_X2Resampler : public Resampler<PolyphaseIIR_X2Resampler<N, M>, 2>
	{
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
			PolyphaseIIR_X2Resampler(const Eigen::Array<double, N, 1>& coeffsDirect,
				const Eigen::Array<double, M, 1>& coeffsDelayed) : _coeffsDirect(coeffsDirect), _coeffsDelayed(coeffsDelayed)
		{
			_sInDirect = Eigen::Array<double, N, 1>::Zero();
			_sInDelayed = Eigen::Array<double, M, 1>::Zero();
			_sOutDirect = Eigen::Array<double, N, 1>::Zero();
			_sOutDelayed = Eigen::Array<double, M, 1>::Zero();
		}


	private:
		friend class Resampler<PolyphaseIIR_X2Resampler<N, M>, 2>;

		Eigen::Array<double, N, 1> _sInDirect;
		Eigen::Array<double, M, 1> _sInDelayed;
		Eigen::Array<double, N, 1> _sOutDirect;
		Eigen::Array<double, M, 1> _sOutDelayed;

		Eigen::Array<double, N, 1> _coeffsDirect;
		Eigen::Array<double, M, 1> _coeffsDelayed;

		double _delay{};

	protected:
		Eigen::Array<double, 2, 1> _Upsample(const double x)
		{
			Eigen::Array<double, 2, 1> x2;

			auto v = x;
			//Direct path
			for (int i = 0; i < N; ++i)
			{
				auto y = _coeffsDirect(i) * v + _sInDirect(i);
				_sInDirect(i) = v - _coeffsDirect(i) * y;
				v = y;
			}
			x2(0) = v;
			v = x;
			//Delayed path
			for (int i = 0; i < M; ++i)
			{
				auto y = _coeffsDelayed(i) * v + _sInDelayed(i);
				_sInDelayed(i) = v - _coeffsDelayed(i) * y;
				v = y;
			}
			x2(1) = v;

			return x2;
		}
		double _Downsample(const Eigen::Array<double, 2, 1> &x2)
		{
			double x = 0.;
			auto v = x2(0);
			//Direct path
			for (int i = 0; i < N; ++i)
			{
				auto y = _coeffsDirect(i) * v + _sOutDirect(i);
				_sOutDirect(i) = v - _coeffsDirect(i) * y;
				v = y;
			}
			x += 0.5 * v;
			v = x2(1);
			//Delayed path
			for (int i = 0; i < M; ++i)
			{
				auto y = _coeffsDelayed(i) * v + _sOutDelayed(i);
				_sOutDelayed(i) = v - _coeffsDelayed(i) * y;
				v = y;
			}
			x += 0.5 * _delay;
			_delay = v;

			return x;
		}
	};
	template<typename X2Type>
	class X4Resampler : public Resampler<X4Resampler<X2Type>, 4>
	{
		std::unique_ptr<X2Type> _stage1;
		std::unique_ptr<X2Type> _stage2;

	public:
		X4Resampler(std::function<X2Type*()> resamplerCreator) : _stage1{ resamplerCreator() }, _stage2(resamplerCreator())
		{
		}
	private:
		friend class Resampler<X4Resampler<X2Type>, 4>;
		Eigen::Array<double, 4, 1> _Upsample(const double x)
		{
			Eigen::Array<double, 4 ,1> x4;
			auto x1 = _stage1->Upsample(x);
			for (int i = 0; i < 2; ++i)
			{
				auto x2 = _stage2->Upsample(x1(i));
				x4(2 * i) = x2(0);
				x4(2 * i + 1) = x2(1);
			}
			return x4;
		}
		double _Downsample(const Eigen::Array<double, 4, 1> &x4)
		{
			Eigen::Array<double, 2, 1> x2;
			x2(0) = x4(0);
			x2(1) = x4(1);
			auto s1 = _stage2->Downsample(x2);
			x2(0) = x4(2);
			x2(1) = x4(3);
			auto s2 = _stage2->Downsample(x2);
			x2(0) = s1;
			x2(1) = s2;

			return _stage1->Downsample(x2);
		}
	};



	using X2Resampler_Order5 = PolyphaseIIR_X2Resampler<1, 1>;
	using X2Resampler_Order7 = PolyphaseIIR_X2Resampler<2, 1>;
	using X2Resampler_Order9 = PolyphaseIIR_X2Resampler<2, 2>;

	using X4Resampler_Order7 = X4Resampler<X2Resampler_Order7>;

	X2Resampler_Order5* CreateX2Resampler_Butterworth5();
	X2Resampler_Order7* CreateX2Resampler_Chebychev7();
	X2Resampler_Order9* CreateX2Resampler_Chebychev9();
	DummyResampler* CreateDummyResampler();
	X4Resampler_Order7* CreateX4Resampler_Cheby7();

}