#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>

#include "tfdsp/rail.hpp"
#include "../tfdsp/filters.hpp"
#include "../tfdsp/approx.hpp"
#include "../tfdsp/sampleRate.hpp"

/**
 * Van der Pol oscillator using a symmetric, conditionally-linear split.
 *
 * Position drift is exact. With position held constant, the velocity equation
 * is linear and is advanced with an analytic implicit-midpoint update. The
 * resulting step is second-order, time-symmetric, and requires no Newton or
 * matrix solve. Frequency prewarping gives exact pitch when damping is zero.
 */
template<typename Oversampler>
class VdpSplitOscillator
{
private:
	static constexpr int ResamplingFactor{Oversampler::ResamplingFactor};
	static constexpr double maxOutput{12.0};
	static constexpr double maxDampingPhasePerStep{0.4};
	static constexpr int maxSubsteps{24};

	double _position{};
	double _velocity{1.0};
	double _sampleRate{};
	double _maxAngularFrequency{};
	std::unique_ptr<Oversampler> _resamplerX;
	std::unique_ptr<Oversampler> _resamplerMu;
	std::unique_ptr<Oversampler> _resamplerW;

	bool VelocityStep(double input, double damping, double interval, double& normalizedVelocity)
	{
		const double rate = damping * (1.0 - _position * _position);
		const double halfRateInterval = 0.5 * rate * interval;
		const double denominator = 1.0 - halfRateInterval;
		if (!std::isfinite(denominator) || std::abs(denominator) < 1.0e-12)
			return false;
		normalizedVelocity = ((1.0 + halfRateInterval) * normalizedVelocity
			+ interval * (input - _position)) / denominator;
		return std::isfinite(normalizedVelocity);
	}

	double ModelStep(double input, double damping, double angularFrequency)
	{
		const double requestedW = std::clamp(angularFrequency, 1.0e-4, _maxAngularFrequency);
		const double mu = std::clamp(damping, 1.0e-8, 9.0);
		const double requestedPhase = requestedW / _sampleRate;
		const int substeps = std::clamp(
			static_cast<int>(std::ceil(mu * requestedPhase / maxDampingPhasePerStep)),
			1, maxSubsteps);
		const double targetPhaseStep = requestedPhase / substeps;
		const double phaseStep = 2.0 * std::sin(0.5 * targetPhaseStep);
		const double effectiveW = phaseStep * _sampleRate * substeps;
		double normalizedVelocity = _velocity / effectiveW;

		for (int i = 0; i < substeps; ++i)
		{
			if (!VelocityStep(input, mu, 0.5 * phaseStep, normalizedVelocity))
				return ResetAfterFailure();
			_position += phaseStep * normalizedVelocity;
			if (!std::isfinite(_position))
				return ResetAfterFailure();
			_position = std::clamp(_position, -maxOutput, maxOutput);
			if (!VelocityStep(input, mu, 0.5 * phaseStep, normalizedVelocity))
				return ResetAfterFailure();
		}

		_velocity = std::clamp(
			effectiveW * normalizedVelocity,
			-2.0 * maxOutput * _sampleRate,
			2.0 * maxOutput * _sampleRate);
		return _position;
	}

	double ResetAfterFailure()
	{
		_position = 0.0;
		_velocity = 1.0;
		return 0.0;
	}

public:
	explicit VdpSplitOscillator(std::function<std::unique_ptr<Oversampler>()> resamplerCreator)
		: _resamplerX(resamplerCreator()),
		  _resamplerMu(resamplerCreator()),
		  _resamplerW(resamplerCreator())
	{
	}

	void SetSampleRate(double sampleRate)
	{
		_sampleRate = sampleRate * ResamplingFactor;
		_maxAngularFrequency = tfdsp::PI * 0.9 * sampleRate;
	}

	void Reset()
	{
		ResetAfterFailure();
		_resamplerX->Reset();
		_resamplerMu->Reset();
		_resamplerW->Reset();
	}

	float Step(double input, double damping, double angularFrequency)
	{
		if (!(_sampleRate > 0.0))
			throw std::runtime_error("Sample rate invalid or not initialized");
		if (!std::isfinite(input) || !std::isfinite(damping) || !std::isfinite(angularFrequency))
		{
			Reset();
			return 0.0f;
		}
		const auto inputValues = _resamplerX->Upsample(input);
		const auto dampingValues = _resamplerMu->Upsample(damping);
		const auto frequencyValues = _resamplerW->Upsample(angularFrequency);
		Eigen::Array<double, ResamplingFactor, 1> output;
		for (int i = 0; i < ResamplingFactor; ++i)
			output(i) = tfdsp::RackOutputAdapter::ProcessOversampled(
				ModelStep(inputValues(i), dampingValues(i), frequencyValues(i)));
		const float result = _resamplerX->Downsample(output);
		if (!std::isfinite(result))
		{
			Reset();
			return 0.0f;
		}
		return static_cast<float>(
			tfdsp::RackOutputAdapter::ProcessPostDecimation(result));
	}

	float StepLogAngularFrequency(double input, double damping,
		double log2AngularFrequency)
	{
		if (!(_sampleRate > 0.0))
			throw std::runtime_error("Sample rate invalid or not initialized");
		if (!std::isfinite(input) || !std::isfinite(damping) ||
			!std::isfinite(log2AngularFrequency))
		{
			Reset();
			return 0.0f;
		}
		const auto inputValues = _resamplerX->Upsample(input);
		const auto dampingValues = _resamplerMu->Upsample(damping);
		const auto frequencyLogValues =
			_resamplerW->Upsample(log2AngularFrequency);
		Eigen::Array<double, ResamplingFactor, 1> output;
		for (int i = 0; i < ResamplingFactor; ++i)
		{
			const double angularFrequency = tfdsp::Exp2Taylor5(
				static_cast<float>(std::clamp(frequencyLogValues(i), -100.0,
					100.0)));
			output(i) = tfdsp::RackOutputAdapter::ProcessOversampled(
				ModelStep(inputValues(i), dampingValues(i), angularFrequency));
		}
		const float result = _resamplerX->Downsample(output);
		if (!std::isfinite(result))
		{
			Reset();
			return 0.0f;
		}
		return static_cast<float>(
			tfdsp::RackOutputAdapter::ProcessPostDecimation(result));
	}
};
