#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>

#include "tfdsp/sampleRate.hpp"

namespace tfdsp
{

// Circuit-scaled model of the ARP 4072 four-pole low-pass filter.
//
// The four filter stages are LM3900 current integrators driven by matched PNP
// pairs. Local resistive feedback keeps those pairs close to their small-signal
// region. A separate PNP pair in front of stage one receives the audio input
// and the resonance return. That pair deliberately operates over a much wider
// base-to-base range and is the circuit's smooth resonance limiter.
//
// States are the AC components of the four LM3900 output voltages. Their large
// negative quiescent voltage is removed analytically; the even number of
// inverting stages and the final level-shifter cancel it in the hardware too.
template<typename ResamplerType, int MaximumNewtonIterations = 8>
class Arp4072Filter
{
public:
	static constexpr int OversamplingFactor = ResamplerType::ResamplingFactor;
	static_assert(MaximumNewtonIterations >= 0,
		"Newton iteration limit must be non-negative");

	explicit Arp4072Filter(
		std::function<std::unique_ptr<ResamplerType>()> resamplerCreator)
		: _resampler(resamplerCreator()),
		  _postOutputResampler(resamplerCreator()),
		  _postLinearCvResampler(resamplerCreator()),
		  _postExponentialCvResampler(resamplerCreator())
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
		Reset();
	}

	void Reset()
	{
		_state = {};
		_resampler->Reset();
		_postOutputResampler->Reset();
		_postLinearCvResampler->Reset();
		_postExponentialCvResampler->Reset();
		_lastIterations = 0;
		_solverFailures = 0;
	}

	float Step(double inputRackVolts, double cutoffHz, double resonance,
		double driveGain = 1.0, bool extendedCutoff = false)
	{
		if (!std::isfinite(inputRackVolts) || !std::isfinite(cutoffHz) ||
			!std::isfinite(resonance) || !std::isfinite(driveGain) ||
			!(_sampleRate > 0.0))
		{
			Reset();
			return 0.0f;
		}

		resonance = std::clamp(resonance, 0.0, 1.0);
		driveGain = std::clamp(driveGain, 0.0, MaximumDriveGain);
		const double circuitCeiling = extendedCutoff ?
			ExtendedCutoffCeilingHz : StockCutoffCeilingHz;
		const double numericalCeiling = 0.45 * _hostSampleRate;
		cutoffHz = SoftLimitCutoff(cutoffHz,
			std::min(circuitCeiling, numericalCeiling));

		const auto upsampled = _resampler->Upsample(inputRackVolts * driveGain);
		Eigen::Array<double, OversamplingFactor, 1> output;
		for (int i = 0; i < OversamplingFactor; ++i)
		{
			const double physicalOutput = OutputLevelShiftGain *
				ProcessOversampled(upsampled(i), cutoffHz, resonance);
			output(i) = SoftOutputCompliance(physicalOutput);
		}

		const double result = _resampler->Downsample(output);
		if (!std::isfinite(result))
		{
			Reset();
			return 0.0f;
		}
		return static_cast<float>(SoftDecimatorSafety(result));
	}

	// Run a second nonlinear stage at the same internal rate before either
	// signal is decimated. Raw linear and exponential control voltages are
	// interpolated independently so audio-rate VCA modulation is also evaluated
	// inside the oversampled path.
	template<typename PostProcessor>
	ProcessedOutputs StepWithPostProcessor(double inputRackVolts,
		double cutoffHz, double resonance, double driveGain,
		bool extendedCutoff, double linearControlVolts,
		double exponentialControlVolts, PostProcessor&& postProcessor)
	{
		if (!std::isfinite(inputRackVolts) || !std::isfinite(cutoffHz) ||
			!std::isfinite(resonance) || !std::isfinite(driveGain) ||
			!std::isfinite(linearControlVolts) ||
			!std::isfinite(exponentialControlVolts) || !(_sampleRate > 0.0))
		{
			Reset();
			return {};
		}

		resonance = std::clamp(resonance, 0.0, 1.0);
		driveGain = std::clamp(driveGain, 0.0, MaximumDriveGain);
		const double circuitCeiling = extendedCutoff ?
			ExtendedCutoffCeilingHz : StockCutoffCeilingHz;
		const double numericalCeiling = 0.45 * _hostSampleRate;
		cutoffHz = SoftLimitCutoff(cutoffHz,
			std::min(circuitCeiling, numericalCeiling));

		const auto audio = _resampler->Upsample(inputRackVolts * driveGain);
		const auto linearCv = _postLinearCvResampler->Upsample(
			linearControlVolts);
		const auto exponentialCv = _postExponentialCvResampler->Upsample(
			exponentialControlVolts);
		Eigen::Array<double, OversamplingFactor, 1> lowPass;
		Eigen::Array<double, OversamplingFactor, 1> postProcessed;
		for (int i = 0; i < OversamplingFactor; ++i)
		{
			const double physicalOutput = OutputLevelShiftGain *
				ProcessOversampled(audio(i), cutoffHz, resonance);
			// The final ARP level shifter and its supply compliance precede the
			// normalled connection to the VCA. Both output paths therefore see
			// exactly the same filter-node level and overload behavior.
			const double limitedOutput = SoftOutputCompliance(physicalOutput);
			lowPass(i) = limitedOutput;
			postProcessed(i) = postProcessor(limitedOutput, linearCv(i),
				exponentialCv(i));
		}

		const double lowPassResult = _resampler->Downsample(lowPass);
		const double postResult = _postOutputResampler->Downsample(postProcessed);
		if (!std::isfinite(lowPassResult) || !std::isfinite(postResult))
		{
			Reset();
			return {};
		}
		return {
			static_cast<float>(SoftDecimatorSafety(lowPassResult)),
			static_cast<float>(SoftDecimatorSafety(postResult)),
		};
	}

	int LastIterations() const { return _lastIterations; }
	std::size_t SolverFailures() const { return _solverFailures; }
	const std::array<double, 4>& State() const { return _state; }

	// Values below are kept public so regression tests can distinguish circuit
	// measurements from later Rack-facing calibration choices.
	static constexpr double ThermalVoltage = 25.85e-3;
	static constexpr double StageInputResistanceOhms = 12100.0;
	static constexpr double StageShuntResistanceOhms = 220.0;
	static constexpr double AudioInputResistanceOhms = 100000.0;
	static constexpr double AudioInputShuntOhms = 220.0;
	static constexpr double FeedbackResistanceOhms = 150000.0;
	static constexpr double FeedbackInputShuntOhms = 2200.0;
	static constexpr double LimiterTailResistanceOhms = 56000.0;
	static constexpr double PositiveSupplyVolts = 15.0;
	static constexpr double LimiterEmitterDropVolts = 0.65;
	static constexpr double OutputLevelShiftGain = 100.0 / 13.0;
	static constexpr double StockCutoffCeilingHz = 12000.0;
	static constexpr double ExtendedCutoffCeilingHz = 20000.0;

	static constexpr double AudioBaseScale()
	{
		return AudioInputShuntOhms /
			(AudioInputResistanceOhms + AudioInputShuntOhms);
	}

	static constexpr double FeedbackBaseScale()
	{
		return FeedbackInputShuntOhms /
			(FeedbackResistanceOhms + FeedbackInputShuntOhms);
	}

	static constexpr double AudioBaseVolts(double externalVolts)
	{
		return AudioBaseScale() * externalVolts;
	}

	static constexpr double FeedbackBaseVolts(double outputVolts)
	{
		return FeedbackBaseScale() * outputVolts;
	}

	static constexpr double LimiterTailCurrentAmps()
	{
		return (PositiveSupplyVolts - LimiterEmitterDropVolts) /
			LimiterTailResistanceOhms;
	}

	static constexpr double LimiterEquivalentPeakVolts()
	{
		return LimiterTailCurrentAmps() * StageInputResistanceOhms;
	}

	static constexpr double StageTanhScalePerVolt()
	{
		return StageShuntResistanceOhms /
			(StageInputResistanceOhms * 2.0 * ThermalVoltage);
	}

	static constexpr double SmallSignalInputGain()
	{
		return OutputLevelShiftGain * LimiterEquivalentPeakVolts() *
			AudioBaseScale() / (2.0 * ThermalVoltage);
	}

	static constexpr double SmallSignalFeedbackGain()
	{
		return OutputLevelShiftGain * LimiterEquivalentPeakVolts() *
			FeedbackBaseScale() / (2.0 * ThermalVoltage);
	}

private:
	static constexpr double Pi = 3.14159265358979323846;
	static constexpr double MaximumDriveGain = 15.848931924611133; // +24 dB
	static constexpr double CutoffFloorKneeHz = 1.0;
	static constexpr double CutoffCeilingKneeHz = 20.0;
	static constexpr double OutputKneeVolts = 10.0;
	static constexpr double OutputRailVolts = 13.5;

	std::unique_ptr<ResamplerType> _resampler;
	std::unique_ptr<ResamplerType> _postOutputResampler;
	std::unique_ptr<ResamplerType> _postLinearCvResampler;
	std::unique_ptr<ResamplerType> _postExponentialCvResampler;
	std::array<double, 4> _state{};
	double _hostSampleRate{};
	double _sampleRate{};
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

	static double SoftLimitCutoff(double requestedHz, double ceilingHz)
	{
		const double positive = Softplus(requestedHz, CutoffFloorKneeHz);
		return ceilingHz - Softplus(ceilingHz - positive,
			CutoffCeilingKneeHz);
	}

	static double SoftOutputCompliance(double voltage)
	{
		const double magnitude = std::abs(voltage);
		if (magnitude <= OutputKneeVolts)
			return voltage;
		const double headroom = OutputRailVolts - OutputKneeVolts;
		return std::copysign(OutputKneeVolts + headroom * std::tanh(
			(magnitude - OutputKneeVolts) / headroom), voltage);
	}

	static double SoftDecimatorSafety(double voltage)
	{
		const double magnitude = std::abs(voltage);
		if (magnitude <= OutputRailVolts)
			return voltage;
		const double excess = magnitude - OutputRailVolts;
		return std::copysign(OutputRailVolts + excess /
			std::sqrt(1.0 + excess * excess), voltage);
	}

	double ProcessOversampled(double input, double cutoffHz, double resonance)
	{
		// The midpoint coefficient is prewarped so the small-signal one-pole
		// sections reach their requested analog cutoff at the oversampled rate.
		const double gamma = 2.0 * std::tan(Pi * cutoffHz / _sampleRate);
		const double stageScale = StageTanhScalePerVolt();
		const double stageStep = gamma / stageScale;
		const double limiterPeak = LimiterEquivalentPeakVolts();
		const double feedbackScale = resonance * FeedbackBaseScale() *
			OutputLevelShiftGain;

		const auto previous = _state;
		auto next = previous;
		bool converged = false;
		_lastIterations = 0;

		for (int iteration = 0; iteration < MaximumNewtonIterations; ++iteration)
		{
			++_lastIterations;
			std::array<double, 4> midpoint{};
			for (int i = 0; i < 4; ++i)
				midpoint[i] = 0.5 * (previous[i] + next[i]);

			const double limiterVoltage = input * AudioBaseScale() -
				feedbackScale * midpoint[3];
			const double limiterTanh = std::tanh(
				limiterVoltage / (2.0 * ThermalVoltage));
			const double firstInput = limiterPeak * limiterTanh;

			std::array<double, 4> stageTanh{};
			stageTanh[0] = std::tanh(stageScale *
				(firstInput + midpoint[0]));
			for (int i = 1; i < 4; ++i)
				stageTanh[i] = std::tanh(stageScale *
					(midpoint[i - 1] + midpoint[i]));

			std::array<double, 4> residual{};
			double maximumResidual = 0.0;
			for (int i = 0; i < 4; ++i)
			{
				residual[i] = next[i] - previous[i] +
					stageStep * stageTanh[i];
				maximumResidual = std::max(maximumResidual,
					std::abs(residual[i]));
			}
			if (maximumResidual < 1.0e-11)
			{
				converged = true;
				break;
			}

			const double limiterSlope = 1.0 - limiterTanh * limiterTanh;
			const double firstInputDerivative = -0.5 * limiterPeak *
				limiterSlope * feedbackScale /
				(2.0 * ThermalVoltage);
			std::array<double, 4> stageSlope{};
			for (int i = 0; i < 4; ++i)
				stageSlope[i] = 1.0 - stageTanh[i] * stageTanh[i];

			double jacobian[4][4]{
				{1.0 + 0.5 * gamma * stageSlope[0], 0.0, 0.0,
					gamma * stageSlope[0] * firstInputDerivative},
				{0.5 * gamma * stageSlope[1],
					1.0 + 0.5 * gamma * stageSlope[1], 0.0, 0.0},
				{0.0, 0.5 * gamma * stageSlope[2],
					1.0 + 0.5 * gamma * stageSlope[2], 0.0},
				{0.0, 0.0, 0.5 * gamma * stageSlope[3],
					1.0 + 0.5 * gamma * stageSlope[3]},
			};
			for (double& value : residual)
				value = -value;
			if (!Solve4x4(jacobian, residual))
				break;

			double maximumDelta = 0.0;
			for (double value : residual)
				maximumDelta = std::max(maximumDelta, std::abs(value));
			const double damping = maximumDelta > 2.0 ?
				2.0 / maximumDelta : 1.0;
			for (int i = 0; i < 4; ++i)
				next[i] += damping * residual[i];
		}

		if (!converged)
		{
			// Never commit an iterate whose residual has not met the solver
			// tolerance. Holding the previous state for one internal sample keeps
			// the result bounded and avoids an arbitrary state jump. The following
			// sample can resume from a known-valid state.
			++_solverFailures;
			return previous[3];
		}
		for (double value : next)
		{
			if (!std::isfinite(value) || std::abs(value) > 100.0)
			{
				Reset();
				return 0.0;
			}
		}
		_state = next;
		return next[3];
	}

	static bool Solve4x4(double matrix[4][4], std::array<double, 4>& rhs)
	{
		for (int column = 0; column < 4; ++column)
		{
			int pivot = column;
			for (int row = column + 1; row < 4; ++row)
			{
				if (std::abs(matrix[row][column]) >
					std::abs(matrix[pivot][column]))
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
				const double factor = matrix[row][column] /
					matrix[column][column];
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
