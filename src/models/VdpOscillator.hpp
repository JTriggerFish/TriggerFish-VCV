#pragma once
#include "VanDerPoleODE.hpp"
#include "ode.hpp"
#include <algorithm>
#include "../tfdsp/sampleRate.hpp"

using namespace ode;

/**
 * \brief Van der Pole Oscillator
 * Note using double as float is too unstable for the way the ode is parameterised.
 * \tparam IntegrationOrder less than 5 is not recommended ( which leaves 5 or 6 as choices )
 */
template<typename Oversampler, int IntegrationOrder>
class VdpOscillator
{
private:
	VanDerPoleODE<double> _vdp{};
	BDF<VanDerPoleODE<double>, double, 2, IntegrationOrder> _integrator{};
	double _sampleRate{};
	StateVector<double,2> _initConditions;

	//More or less C8 on a piano, no point going higher and things start to blow up anyway
	static constexpr double maxAngularFreq = double(2 * PI * 4200);
	//Clamp the output to avoid blow ups
	static constexpr double maxOutput = 12.0;

	static constexpr double _initY0 = double(0.0);
	static constexpr double _initY1 = double(1.0);

	std::unique_ptr<Oversampler> _resamplerX;
	std::unique_ptr<Oversampler> _resamplerMu;
	std::unique_ptr<Oversampler> _resamplerW;

	static constexpr int ResamplingFactor{ Oversampler::ResamplingFactor };

	//tfdsp::PinkNoiseSource _noise{};

	float ModelStep(double x, double mu, double w)
	{
		_vdp._mu = std::max<double>(1.0e-8, mu);
		_vdp._w = std::max<double>(1.0e-4, std::min<double>(maxAngularFreq, w));
		_integrator.Step(_vdp, x);// +_noise.Step());

		//Clamp the state to avoid exploding if the ODE becomes unstable, typically when both mu and w are high.
		auto state = _integrator.FullState();
		for (int i = 0; i < 2; ++i)
		{
			state(0, i) = std::max<double>(std::min<double>(state(0, i), maxOutput), -maxOutput);
			auto derivMax = 2 * maxOutput * _sampleRate;
			state(1, i) = std::max<double>(std::min<double>(state(1, i), derivMax), -derivMax);
		}

		return _integrator.CurrentState()(0);
	}
public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	explicit VdpOscillator(std::function<Oversampler*()> resamplerCreator) : _resamplerX(resamplerCreator()),
		_resamplerMu(resamplerCreator()),
		_resamplerW(resamplerCreator())
	{
		_initConditions << _initY0, _initY1;
		_integrator.SetInitConditions(_initConditions);
	}

	void SetSampleRate(const double f0)
	{
		_sampleRate = f0 * ResamplingFactor;
		StateVector<double, 2> currentState = _integrator.CurrentState();
		_initConditions << currentState(0), currentState(1);
		_integrator.SetInitConditions(_initConditions);
		_integrator.SetSampleRate(_sampleRate);
	}

	float Step(double x, double mu, double w)
	{
		if (_sampleRate <= 0.)
		{
			throw std::runtime_error("Sample rate invalid or not initialized");
		}
		auto xA = _resamplerX->Upsample(x);
		auto muA = _resamplerMu->Upsample(mu);
		auto wA = _resamplerW->Upsample(w);

		Eigen::Matrix<double, ResamplingFactor, 1> output;

		for(int i=0; i < ResamplingFactor; ++i)
		{
			output(i) =  ModelStep(xA(i), muA(i), wA(i));
		}
		return _resamplerX->Downsample(output);

	}
};

//Fix linking errors for old versions of gcc
template<typename Oversampler, int IntegrationOrder>
constexpr double VdpOscillator<Oversampler, IntegrationOrder>::maxAngularFreq;
template<typename Oversampler, int IntegrationOrder>
constexpr double VdpOscillator<Oversampler, IntegrationOrder>::maxOutput;
template<typename Oversampler, int IntegrationOrder>
constexpr double VdpOscillator<Oversampler, IntegrationOrder>::_initY0;
template<typename Oversampler, int IntegrationOrder>
constexpr double VdpOscillator<Oversampler, IntegrationOrder>::_initY1;
template<typename Oversampler, int IntegrationOrder>
constexpr int VdpOscillator<Oversampler, IntegrationOrder>::ResamplingFactor;
