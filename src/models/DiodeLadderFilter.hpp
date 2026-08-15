#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <utility>

#include "AnalogOutputStage.hpp"
#include "tfdsp/filters.hpp"
#include "tfdsp/sampleRate.hpp"

namespace tfdsp
{

// Bilinear/TPT realization of the analog first-order section
//
//              s + zero
//     H(s) =  ---------- .
//              s + pole
//
// It exposes the current output as an affine function of the current input.
// This lets the complete AC-coupled resonance path participate in the
// zero-delay nonlinear ladder solve instead of acquiring a sample delay.
class AnalogRatioSection
{
public:
	struct Affine
	{
		double gain{};
		double offset{};
	};

	void Configure(double poleRadians, double zeroRadians, double sampleRate)
	{
		const double safePole = std::max(poleRadians, 1.0e-12);
		const double g = std::tan(safePole / (2.0 * sampleRate));
		_integratorGain = g / (1.0 + g);
		_stateGain = 1.0 / (1.0 + g);
		_mix = zeroRadians / safePole - 1.0;
	}

	Affine Preview() const
	{
		return {1.0 + _mix * _integratorGain,
			_mix * _stateGain * _state};
	}

	double Process(double input)
	{
		const double lowPass = _integratorGain * input + _stateGain * _state;
		const double output = input + _mix * lowPass;
		_state = 2.0 * lowPass - _state;
		return output;
	}

	void Reset()
	{
		_state = 0.0;
	}

private:
	double _integratorGain{};
	double _stateGain{1.0};
	double _mix{};
	double _state{};
};

template<std::size_t Size>
class AnalogRatioCascade
{
public:
	using PoleZero = std::pair<double, double>;

	void Configure(const std::array<PoleZero, Size>& poleZeros,
		double sampleRate, double gain)
	{
		_gain = gain;
		for (std::size_t i = 0; i < Size; ++i)
			_sections[i].Configure(poleZeros[i].first,
				poleZeros[i].second, sampleRate);
	}

	typename AnalogRatioSection::Affine Preview() const
	{
		double gain = 1.0;
		double offset = 0.0;
		for (const auto& section : _sections)
		{
			const auto affine = section.Preview();
			gain = affine.gain * gain;
			offset = affine.gain * offset + affine.offset;
		}
		return {_gain * gain, _gain * offset};
	}

	double Process(double input)
	{
		for (auto& section : _sections)
			input = section.Process(input);
		return _gain * input;
	}

	void Reset()
	{
		for (auto& section : _sections)
			section.Reset();
	}

private:
	std::array<AnalogRatioSection, Size> _sections{};
	double _gain{1.0};
};

// Circuit-structured nonlinear model of the four-capacitor TB-303 diode
// ladder. The first capacitor is 18 nF rather than the other stages' 33 nF.
// Stinchcombe's idealized 1:2 ratio produces the small-signal polynomial
//
// s^4 + 2^(11/4)s^3 + 10*sqrt(2)s^2 + 2^(13/4)s + 1.
//
// The surrounding input and resonance coupling sections implement the six
// poles and six zeroes of the published complete linear model. The ladder
// itself uses nonlinear diode currents and an implicit midpoint update with an
// analytic Jacobian. This implementation uses the actual 33/18 stage-current
// ratio; its coefficients differ from the idealized polynomial by less than
// 0.2%. Feedback-filter feedthrough is included in that solve.
template<typename ResamplerType>
class DiodeLadderFilter
{
public:
	static constexpr int OversamplingFactor = ResamplerType::ResamplingFactor;

	explicit DiodeLadderFilter(
		std::function<std::unique_ptr<ResamplerType>()> resamplerCreator)
		: _resampler(resamplerCreator()),
		  _postResampler(resamplerCreator())
	{
	}

	struct ProcessedOutputs
	{
		float lowPass{};
		float postProcessed{};
	};

	void SetSampleRate(double sampleRate)
	{
		if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
			return;
		_hostSampleRate = sampleRate;
		_sampleRate = sampleRate * OversamplingFactor;

		// The constants are analog angular frequencies (rad/s) from the
		// complete TB-303 transfer function. One common (s + 7.41) factor in
		// the forward path is cancelled before realization.
		_forward.Configure({{
			{578.1, 0.0},
			{97.5, 0.0},
			{20.0, 109.9},
			{4.45, 34.0},
		}}, _sampleRate, 1.06);
		// Place the lowest forward-path coupling section after the nonlinear
		// ladder. This leaves the published small-signal transfer unchanged,
		// while correctly rejecting DC generated inside the asymmetric driven
		// ladder before the Devil Fish bass shelves can amplify it.
		_outputCoupling.Configure(38.5, 0.0, _sampleRate);

		_feedback.Configure({{
			{578.1, 0.0},
			{97.5, 0.0},
			{38.5, 0.0},
			{20.0, 0.0},
			{7.41, 46.5},
			{4.45, 4.40},
		}}, _sampleRate, 18.7);

		_bassSmoothing = 1.0 - std::exp(-1.0 / (0.010 * _sampleRate));
		_driveSmoothing = 1.0 - std::exp(-1.0 / (0.010 * _sampleRate));
		ConfigureBassCorrection();
		Reset();
	}

	void Reset()
	{
		_state = {1.0e-12, 0.0, 0.0, 0.0};
		_forward.Reset();
		_feedback.Reset();
		_outputCoupling.Reset();
		for (auto& section : _bassCorrection)
			section.Reset();
		_resampler->Reset();
		_postResampler->Reset();
		_smoothedBass = 0.0;
		_smoothedDrive = 0.0;
		_configuredBass = -1.0;
		_lastIterations = 0;
		_solverFailures = 0;
	}

	float Step(double inputVolts, double cutoffHz, double resonance,
		bool highResonance, double driveGain, double bass)
	{
		if (!std::isfinite(inputVolts) || !std::isfinite(cutoffHz) ||
			!std::isfinite(resonance) || !std::isfinite(driveGain) ||
			!std::isfinite(bass) || !(_sampleRate > 0.0))
		{
			Reset();
			return 0.0f;
		}

		const double maximumCutoff = std::min(20000.0, 0.45 * _hostSampleRate);
		cutoffHz = MapCutoffControl(cutoffHz, maximumCutoff);
		resonance = std::clamp(resonance, 0.0, 1.0);
		driveGain = std::clamp(driveGain, 0.0, 66.6);
		bass = std::clamp(bass, 0.0, 1.0);

		// The service schematic AC-couples the approximately 5.5 Vpp stock saw
		// through R62 = 220k into the input node held by R70 = 2.2k. In the
		// pair's tanh(v / 2VT) coordinates this is about 1.05 normalized Vpp.
		// Map a nominal 10 Vpp Rack oscillator to that circuit drive. The
		// Devil Fish range then extends to 66.6 times the stock level.
		const auto upsampled = _resampler->Upsample(inputVolts * StockInputScale);
		Eigen::Array<double, OversamplingFactor, 1> output;
		const double resonanceMakeup = 1.0 + resonance *
			(highResonance ? HighResonanceMakeup : StockResonanceMakeup);
		const double outputScale = RackOutputScale * resonanceMakeup;
		for (int i = 0; i < OversamplingFactor; ++i)
			output(i) = AnalogOutputStage::Process(outputScale *
				ProcessOversampled(upsampled(i), cutoffHz, resonance,
					highResonance, driveGain, bass));

		// Open303 and tbvcf both use resonance-dependent output makeup. The
		// calibration is based on AC signal level after the output coupling
		// section, retaining some authentic thinning without counting nonlinear
		// DC offset as useful output level.
		const double result = _resampler->Downsample(output);
		if (!std::isfinite(result))
		{
			Reset();
			return 0.0f;
		}
		return static_cast<float>(AnalogOutputStage::ProcessSafety(result));
	}

	// Process a second nonlinear stage before decimation. The second resampler
	// provides independent output-filter state and interpolates its control
	// signal, while the audio input is upsampled only once.
	template<typename PostProcessor>
	ProcessedOutputs StepWithPostProcessor(double inputVolts, double cutoffHz,
		double resonance, bool highResonance, double driveGain, double bass,
		double postControl, PostProcessor&& postProcessor)
	{
		if (!std::isfinite(inputVolts) || !std::isfinite(cutoffHz) ||
			!std::isfinite(resonance) || !std::isfinite(driveGain) ||
			!std::isfinite(bass) || !std::isfinite(postControl) ||
			!(_sampleRate > 0.0))
		{
			Reset();
			return {};
		}

		const double maximumCutoff = std::min(20000.0, 0.45 * _hostSampleRate);
		cutoffHz = MapCutoffControl(cutoffHz, maximumCutoff);
		resonance = std::clamp(resonance, 0.0, 1.0);
		driveGain = std::clamp(driveGain, 0.0, 66.6);
		bass = std::clamp(bass, 0.0, 1.0);

		const auto upsampled = _resampler->Upsample(inputVolts * StockInputScale);
		const auto upsampledControl = _postResampler->Upsample(postControl);
		Eigen::Array<double, OversamplingFactor, 1> lowPass;
		Eigen::Array<double, OversamplingFactor, 1> postProcessed;
		const double resonanceMakeup = 1.0 + resonance *
			(highResonance ? HighResonanceMakeup : StockResonanceMakeup);
		const double lowPassOutputScale = RackOutputScale * resonanceMakeup;
		const double vcaInputScale = RackOutputScale;
		for (int i = 0; i < OversamplingFactor; ++i)
		{
			const double filtered = ProcessOversampled(upsampled(i), cutoffHz,
				resonance, highResonance, driveGain, bass);
			lowPass(i) = AnalogOutputStage::Process(
				lowPassOutputScale * filtered);
			// Resonance makeup is a Rack output calibration rather than part of
			// the circuit. Apply it after the nonlinear VCA so it cannot change
			// the BA662 drive and then pass it through the modeled output rail.
			postProcessed(i) = AnalogOutputStage::Process(resonanceMakeup *
				postProcessor(vcaInputScale * filtered, upsampledControl(i)));
		}

		const double lowPassResult = _resampler->Downsample(lowPass);
		const double postResult = _postResampler->Downsample(postProcessed);
		if (!std::isfinite(lowPassResult) || !std::isfinite(postResult))
		{
			Reset();
			return {};
		}
		return {
			static_cast<float>(
				AnalogOutputStage::ProcessSafety(lowPassResult)),
			static_cast<float>(AnalogOutputStage::ProcessSafety(postResult)),
		};
	}

	int LastIterations() const { return _lastIterations; }
	std::size_t SolverFailures() const { return _solverFailures; }

	// Cutoff is proportional to a unidirectional transistor control current.
	// Linear FM can request a negative current; the hardware pinches off through
	// the device curve instead of stopping at an arbitrary frequency boundary.
	// A second soft knee protects the discrete model near its usable ceiling.
	static double MapCutoffControl(double requestedHz, double maximumHz)
	{
		const double positive = Softplus(requestedHz, CutoffPinchKneeHz);
		return maximumHz - Softplus(maximumHz - positive,
			CutoffCeilingKneeHz);
	}

private:
	static constexpr double FirstStageScale = 33.0 / 18.0;
	static constexpr double CapacitorScale = 0.8593887047640296;
	static constexpr double StockInputScale = 0.10532968190436065;
	// Convert the ladder's normalized voltage back to the nominal Rack scale.
	// This is the reciprocal of the stock input mapping, so an open, linear
	// ladder has approximately unity signal gain before the surrounding network.
	static constexpr double RackOutputScale = 9.494;
	static constexpr double StockResonanceScale = 0.78;
	static constexpr double HighResonanceMultiplier = 2.0;
	static constexpr double StockResonanceMakeup = 2.0;
	static constexpr double HighResonanceMakeup = 3.0;
	static constexpr double BassPole = 2.0 * PI * 24.66;
	static constexpr double CutoffPinchKneeHz = 1.0;
	static constexpr double CutoffCeilingKneeHz = 10.0;

	std::unique_ptr<ResamplerType> _resampler;
	std::unique_ptr<ResamplerType> _postResampler;
	AnalogRatioCascade<4> _forward;
	AnalogRatioCascade<6> _feedback;
	AnalogRatioSection _outputCoupling;
	std::array<AnalogRatioSection, 2> _bassCorrection{};
	std::array<double, 4> _state{};
	double _hostSampleRate{};
	double _sampleRate{};
	double _smoothedBass{};
	double _smoothedDrive{};
	double _configuredBass{-1.0};
	double _bassSmoothing{};
	double _driveSmoothing{};
	int _lastIterations{};
	std::size_t _solverFailures{};

	static double Softplus(double value, double knee)
	{
		const double normalized = value / knee;
		if (normalized > 40.0)
			return value;
		if (normalized < -40.0)
			return knee * std::exp(normalized);
		return knee * std::log1p(std::exp(normalized));
	}

	void ConfigureBassCorrection()
	{
		// C20/C21 are increased by a factor of ten in the bass modification.
		// Two (s + stockPole)/(s + variedPole) shelves correct the complete
		// stock network, giving approximately the documented change at 32 Hz
		// from -5 dB to -1 dB while retaining DC blocking.
		const double variedPole = BassPole * std::pow(0.1, _smoothedBass);
		for (auto& section : _bassCorrection)
			section.Configure(variedPole, BassPole, _sampleRate);
		_configuredBass = _smoothedBass;
	}

	double ProcessOversampled(double input, double cutoffHz, double resonance,
		bool highResonance, double driveGain, double bass)
	{
		_smoothedDrive += _driveSmoothing * (driveGain - _smoothedDrive);
		_smoothedBass += _bassSmoothing * (bass - _smoothedBass);
		if (std::abs(_smoothedBass - bass) < 1.0e-8)
			_smoothedBass = bass;
		if (_configuredBass < 0.0 ||
			std::abs(_smoothedBass - _configuredBass) > 1.0e-7)
			ConfigureBassCorrection();

		const double forward = _forward.Process(input * _smoothedDrive);
		const double feedbackAmount = resonance * StockResonanceScale *
			(highResonance ? HighResonanceMultiplier : 1.0);
		const auto feedbackAffine = _feedback.Preview();

		// Prewarp at the oversampled rate. Multiplication by the sample period
		// is folded in, leaving the dimensionless midpoint step coefficient.
		const double step = CapacitorScale * 2.0 *
			std::tan(PI * cutoffHz / _sampleRate);

		const auto previous = _state;
		auto next = previous;
		bool converged = false;
		_lastIterations = 0;

		for (int iteration = 0; iteration < 8; ++iteration)
		{
			++_lastIterations;
			std::array<double, 4> midpoint{};
			for (int i = 0; i < 4; ++i)
				midpoint[i] = 0.5 * (previous[i] + next[i]);

			const double feedback = feedbackAffine.gain * next[3] +
				feedbackAffine.offset;
			const double inputJunction = forward - feedbackAmount * feedback;
			const double junction0 = std::tanh(inputJunction);
			const double junction1 = std::tanh(midpoint[0] - midpoint[1]);
			const double junction2 = std::tanh(midpoint[1] - midpoint[2]);
			const double junction3 = std::tanh(midpoint[2] - midpoint[3]);
			const double junction4 = std::tanh(midpoint[3]);

			const std::array<double, 4> derivative{{
				FirstStageScale * (junction0 - junction1),
				junction1 - junction2,
				junction2 - junction3,
				junction3 - junction4,
			}};
			std::array<double, 4> residual{};
			double maximumResidual = 0.0;
			for (int i = 0; i < 4; ++i)
			{
				residual[i] = next[i] - previous[i] - step * derivative[i];
				maximumResidual = std::max(maximumResidual, std::abs(residual[i]));
			}
			if (maximumResidual < 1.0e-11)
			{
				converged = true;
				break;
			}

			const double slope0 = 1.0 - junction0 * junction0;
			const double slope1 = 1.0 - junction1 * junction1;
			const double slope2 = 1.0 - junction2 * junction2;
			const double slope3 = 1.0 - junction3 * junction3;
			const double slope4 = 1.0 - junction4 * junction4;

			double jacobian[4][4]{
				{1.0 + 0.5 * step * FirstStageScale * slope1,
					-0.5 * step * FirstStageScale * slope1, 0.0,
					step * FirstStageScale * slope0 * feedbackAmount * feedbackAffine.gain},
				{-0.5 * step * slope1,
					1.0 + 0.5 * step * (slope1 + slope2),
					-0.5 * step * slope2, 0.0},
				{0.0, -0.5 * step * slope2,
					1.0 + 0.5 * step * (slope2 + slope3),
					-0.5 * step * slope3},
				{0.0, 0.0, -0.5 * step * slope3,
					1.0 + 0.5 * step * (slope3 + slope4)},
			};
			for (double& value : residual)
				value = -value;
			if (!Solve4x4(jacobian, residual))
				break;

			double maximumDelta = 0.0;
			for (double value : residual)
				maximumDelta = std::max(maximumDelta, std::abs(value));
			const double damping = maximumDelta > 1.0 ? 1.0 / maximumDelta : 1.0;
			for (int i = 0; i < 4; ++i)
				next[i] += damping * residual[i];
		}

		if (!converged)
			++_solverFailures;
		for (double value : next)
		{
			if (!std::isfinite(value) || std::abs(value) > 100.0)
			{
				Reset();
				return 0.0;
			}
		}

		_state = next;
		_feedback.Process(next[3]);
		double output = _outputCoupling.Process(next[3]);
		for (auto& section : _bassCorrection)
			output = section.Process(output);
		return output;
	}

	static bool Solve4x4(double matrix[4][4], std::array<double, 4>& rhs)
	{
		for (int column = 0; column < 4; ++column)
		{
			int pivot = column;
			for (int row = column + 1; row < 4; ++row)
			{
				if (std::abs(matrix[row][column]) > std::abs(matrix[pivot][column]))
					pivot = row;
			}
			if (std::abs(matrix[pivot][column]) < 1.0e-14)
				return false;
			if (pivot != column)
			{
				for (int i = column; i < 4; ++i)
					std::swap(matrix[column][i], matrix[pivot][i]);
				std::swap(rhs[column], rhs[pivot]);
			}

			for (int row = column + 1; row < 4; ++row)
			{
				const double factor = matrix[row][column] / matrix[column][column];
				for (int i = column + 1; i < 4; ++i)
					matrix[row][i] -= factor * matrix[column][i];
				rhs[row] -= factor * rhs[column];
			}
		}

		for (int row = 3; row >= 0; --row)
		{
			double value = rhs[row];
			for (int column = row + 1; column < 4; ++column)
				value -= matrix[row][column] * rhs[column];
			rhs[row] = value / matrix[row][row];
		}
		return true;
	}
};

} // namespace tfdsp
