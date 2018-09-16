#pragma once
#include <array>
#include <memory>
#include <functional>
#include "filters.hpp"
#include "util.hpp"

/***
 * References : Digital signal processing schemes for efficient interpolation and decimation, R.A Valenzuela, A.G Constantinides, IEE 1983
 */

namespace dsp {


template<typename Impl, unsigned int Factor>
class Resampler : public enable_down_cast<Impl>
{
private:
    using enable_down_cast< Impl >::Self;
public:
    static constexpr unsigned int ResamplingFactor{Factor};
    inline std::array<double, Factor> Upsample(const double x)
    {
        return Self()->_Upsample(x);
    }
    inline double Downsample(const std::array<double, Factor> &x2)
    {
        return Self()->_Downsample(x2);
    }
};

class DummyResampler : public Resampler<DummyResampler, 1>
{
public:
    DummyResampler() {}
    friend class Resampler<DummyResampler,1>;
protected:
    inline std::array<double,1> _Upsample(const double x)
    {
        std::array<double,1> xA{{x}};
        return xA;
    }
    inline double _Downsample(const std::array<double,1>& xA)
    {
        return xA[0];
    }
};

/// N direct stages, M delayed stages
template<unsigned int N, unsigned int M>
class PolyphaseIIR_X2Resampler : public Resampler<PolyphaseIIR_X2Resampler<N,M>,2>
{
public:
    PolyphaseIIR_X2Resampler(const std::array<double,N>& coeffsDirect,
    const std::array<double,M>& coeffsDelayed) : _coeffsDirect(coeffsDirect), _coeffsDelayed(coeffsDelayed) {}


private:
    friend class Resampler<PolyphaseIIR_X2Resampler<N,M>, 2>;

    std::array<double, N> _sInDirect{};
    std::array<double, M> _sInDelayed{};
    std::array<double, N> _sOutDirect{};
    std::array<double, M> _sOutDelayed{};
    std::array<double, N> _coeffsDirect{};
    std::array<double, M> _coeffsDelayed{};

	double _delay{};

protected:
    std::array<double, 2> _Upsample(const double x)
    {
        std::array<double, 2> x2{};

        auto v = x;
        //Direct path
        for (unsigned int i = 0; i < N; ++i)
        {
            auto y = _coeffsDirect[i] * v + _sInDirect[i];
            _sInDirect[i] = v - _coeffsDirect[i] * y;
            v = y;
        }
        x2[0] = v;
        v = x;
        //Delayed path
        for (unsigned int i = 0; i < M; ++i)
        {
            auto y = _coeffsDelayed[i] * v + _sInDelayed[i];
            _sInDelayed[i] = v - _coeffsDelayed[i] * y;
            v = y;
        }
        x2[1] = v;

        return x2;
    }
    double _Downsample(const std::array<double, 2> &x2)
    {
        double x = 0.;
        auto v = x2[0];
        //Direct path
        for (unsigned int i = 0; i < N; ++i)
        {
            auto y = _coeffsDirect[i] * v + _sOutDirect[i];
            _sOutDirect[i] = v - _coeffsDirect[i] * y;
            v = y;
        }
        x += 0.5 * v;
        v = x2[1];
        //Delayed path
        for (unsigned int i = 0; i < M; ++i)
        {
            auto y = _coeffsDelayed[i] * v + _sOutDelayed[i];
            _sOutDelayed[i] = v - _coeffsDelayed[i] * y;
            v = y;
        }
        x += 0.5 * _delay;
		_delay = v;

        return x;
    }
};
template<typename X2Type>
class X4Resampler : public Resampler<X4Resampler<X2Type>,4>
{
	std::unique_ptr<X2Type> _stage1;
	std::unique_ptr<X2Type> _stage2;
	
public:
    X4Resampler(std::function<X2Type*()> resamplerCreator) : _stage1{resamplerCreator()}, _stage2(resamplerCreator())
    {
    }
private:
	friend class Resampler<X4Resampler<X2Type>, 4>;
    std::array<double, 4> _Upsample(const double x)
    {
		std::array<double, 4> x4;
		auto x1 = _stage1->Upsample(x);
		for(unsigned int i=0; i < 2; ++i)
		{
			auto x2 = _stage2->Upsample(x1[i]);
			x4[2 * i] = x2[0];
			x4[2 * i + 1] = x2[1];
		}
		return x4;
    }
    double _Downsample(const std::array<double, 4> &x4)
    {
		std::array<double, 2> x2;
		x2[0] = x4[0];
		x2[1] = x4[1];
		auto s1 = _stage2->Downsample(x2);
		x2[0] = x4[2];
		x2[1] = x4[3];
		auto s2 = _stage2->Downsample(x2);
		x2[0] = s1;
		x2[1] = s2;

		return _stage1->Downsample(x2);
    }
};



using X2Resampler_Order5 = PolyphaseIIR_X2Resampler<1,1>;
using X2Resampler_Order7 = PolyphaseIIR_X2Resampler<2,1>;
using X2Resampler_Order9 = PolyphaseIIR_X2Resampler<2,2>;

using X4Resampler_Order7 = X4Resampler<X2Resampler_Order7>;

 X2Resampler_Order5* CreateX2Resampler_Butterworth5();
 X2Resampler_Order7* CreateX2Resampler_Chebychev7();
 X2Resampler_Order9* CreateX2Resampler_Chebychev9();
 DummyResampler* CreateDummyResampler();
 X4Resampler_Order7* CreateX4Resampler_Cheby7();

}