#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <utility>

#include "tfdsp/sampleRate.hpp"
#include "tfdsp/approx.hpp"
#include "tfdsp/rail.hpp"

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
		: _audioInputResampler(resamplerCreator()),
		  _lowPassOutputResampler(resamplerCreator()),
		  _cutoffPitchResampler(resamplerCreator()),
		  _linearFmResampler(resamplerCreator()),
		  _resonanceResampler(resamplerCreator()),
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
		_audioInputResampler->Reset();
		_lowPassOutputResampler->Reset();
		_cutoffPitchResampler->Reset();
		_linearFmResampler->Reset();
		_resonanceResampler->Reset();
		_postOutputResampler->Reset();
		_postLinearCvResampler->Reset();
		_postExponentialCvResampler->Reset();
		_lastIterations = 0;
		_solverFailures = 0;
		_renderedLowPassLastSample = true;
	}

	float Step(double inputRackVolts, double cutoffHz, double resonance,
		double driveGain = 1.0)
	{
		if (!std::isfinite(cutoffHz))
		{
			Reset();
			return 0.0f;
		}
		return StepLogCutoff(inputRackVolts, std::log2(std::max(cutoffHz,
			std::numeric_limits<double>::min())), resonance, driveGain);
	}

	float StepLogCutoff(double inputRackVolts, double log2CutoffHz,
		double resonance, double driveGain = 1.0)
	{
		return StepModulatedLogCutoff(inputRackVolts, log2CutoffHz, 0.0,
			resonance, driveGain);
	}

	float StepModulatedLogCutoff(double inputRackVolts, double log2CutoffHz,
		double linearFmHz, double resonance, double driveGain = 1.0)
	{
		if (!std::isfinite(inputRackVolts) || !std::isfinite(log2CutoffHz) ||
			!std::isfinite(linearFmHz) ||
			!std::isfinite(resonance) || !std::isfinite(driveGain) ||
			!(_sampleRate > 0.0))
		{
			Reset();
			return 0.0f;
		}

		driveGain = std::clamp(driveGain, 0.0, MaximumDriveGain);
		const double numericalCeiling = 0.45 * _hostSampleRate;
		const auto controls = UpsampleControls(log2CutoffHz, linearFmHz, resonance,
			std::min(CutoffCeilingHz, numericalCeiling));

		const auto upsampled = _audioInputResampler->Upsample(
			inputRackVolts * driveGain);
		Eigen::Array<double, OversamplingFactor, 1> output;
		for (int i = 0; i < OversamplingFactor; ++i)
		{
			const double physicalOutput = OutputLevelShiftGain *
				ProcessOversampled(upsampled(i), controls.cutoffHz(i),
					controls.resonance(i));
			output(i) = RackOutputAdapter::ProcessOversampled(
				SoftOutputCompliance(physicalOutput));
		}

		const double result = DownsampleLowPass(output, true);
		if (!std::isfinite(result))
		{
			Reset();
			return 0.0f;
		}
		return static_cast<float>(
			RackOutputAdapter::ProcessPostDecimation(result));
	}

	// Run a second nonlinear stage at the same internal rate before either
	// signal is decimated. Raw linear and exponential control voltages are
	// interpolated independently so audio-rate VCA modulation is also evaluated
	// inside the oversampled path.
	template<typename PostProcessor>
	ProcessedOutputs StepWithPostProcessor(double inputRackVolts,
		double cutoffHz, double resonance, double driveGain,
		double linearControlVolts,
		double exponentialControlVolts, PostProcessor&& postProcessor,
		bool renderLowPass = true)
	{
		if (!std::isfinite(cutoffHz))
		{
			Reset();
			return {};
		}
		return StepWithPostProcessorLogCutoff(inputRackVolts,
			std::log2(std::max(cutoffHz,
				std::numeric_limits<double>::min())), resonance, driveGain,
			linearControlVolts, exponentialControlVolts,
			std::forward<PostProcessor>(postProcessor), renderLowPass);
	}

	template<typename PostProcessor>
	ProcessedOutputs StepWithPostProcessorLogCutoff(double inputRackVolts,
		double log2CutoffHz, double resonance, double driveGain,
		double linearControlVolts,
		double exponentialControlVolts, PostProcessor&& postProcessor,
		bool renderLowPass = true)
	{
		return StepWithPostProcessorModulatedLogCutoff(inputRackVolts,
			log2CutoffHz, 0.0, resonance, driveGain, linearControlVolts,
			exponentialControlVolts, std::forward<PostProcessor>(postProcessor),
			renderLowPass);
	}

	template<typename PostProcessor>
	ProcessedOutputs StepWithPostProcessorModulatedLogCutoff(
		double inputRackVolts, double log2CutoffHz, double linearFmHz,
		double resonance, double driveGain, double linearControlVolts,
		double exponentialControlVolts, PostProcessor&& postProcessor,
		bool renderLowPass = true)
	{
		if (!std::isfinite(inputRackVolts) || !std::isfinite(log2CutoffHz) ||
			!std::isfinite(linearFmHz) ||
			!std::isfinite(resonance) || !std::isfinite(driveGain) ||
			!std::isfinite(linearControlVolts) ||
			!std::isfinite(exponentialControlVolts) || !(_sampleRate > 0.0))
		{
			Reset();
			return {};
		}

		driveGain = std::clamp(driveGain, 0.0, MaximumDriveGain);
		const double numericalCeiling = 0.45 * _hostSampleRate;
		const auto controls = UpsampleControls(log2CutoffHz, linearFmHz, resonance,
			std::min(CutoffCeilingHz, numericalCeiling));

		const auto audio = _audioInputResampler->Upsample(
			inputRackVolts * driveGain);
		const auto linearCv = _postLinearCvResampler->Upsample(
			linearControlVolts);
		const auto exponentialCv = _postExponentialCvResampler->Upsample(
			exponentialControlVolts);
		Eigen::Array<double, OversamplingFactor, 1> lowPass;
		Eigen::Array<double, OversamplingFactor, 1> postProcessed;
		for (int i = 0; i < OversamplingFactor; ++i)
		{
			const double physicalOutput = OutputLevelShiftGain *
				ProcessOversampled(audio(i), controls.cutoffHz(i),
					controls.resonance(i));
			// The final ARP level shifter and its supply compliance precede the
			// normalled connection to the VCA. Both output paths therefore see
			// exactly the same filter-node level and overload behavior.
			const double limitedOutput = SoftOutputCompliance(physicalOutput);
			if (renderLowPass)
				lowPass(i) = RackOutputAdapter::ProcessOversampled(limitedOutput);
			postProcessed(i) = RackOutputAdapter::ProcessOversampled(
				postProcessor(limitedOutput, linearCv(i), exponentialCv(i)));
		}

		const double lowPassResult = DownsampleLowPass(lowPass, renderLowPass);
		const double postResult = _postOutputResampler->Downsample(postProcessed);
		if (!std::isfinite(lowPassResult) || !std::isfinite(postResult))
		{
			Reset();
			return {};
		}
		return {
			static_cast<float>(
				RackOutputAdapter::ProcessPostDecimation(lowPassResult)),
			static_cast<float>(
				RackOutputAdapter::ProcessPostDecimation(postResult)),
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
	static constexpr double CutoffCeilingHz = 20000.0;

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
		return NominalLimiterEquivalentPeakVolts() * LimiterGainCalibration();
	}

	static constexpr double LimiterDifferentialResistanceOhms()
	{
		// The reference-side limiter base sees 220 ohm; the driven side sees
		// the loaded stage resistance. Differential current develops voltage
		// across their mean resistance.
		return 0.5 * (StageShuntResistanceOhms + StageBaseResistanceOhms());
	}

	static constexpr double NominalLimiterEquivalentPeakVolts()
	{
		return LimiterTailCurrentAmps() * StageInputResistanceOhms *
			LimiterDifferentialResistanceOhms() / StageBaseResistanceOhms();
	}

	static constexpr double LimiterGainCalibration()
	{
		// The module's pin-10 gain trim and the console's R163 adjustment set
		// the open filter to unity. Keep that service calibration separate from
		// the nominal component-ratio estimate above.
		const double unityPeak = 2.0 * ThermalVoltage /
			(OutputLevelShiftGain * AudioBaseScale());
		return unityPeak / NominalLimiterEquivalentPeakVolts();
	}

	static constexpr double StageSaturationCoefficientPerVolt()
	{
		return StageBaseResistanceOhms() /
			(StageInputResistanceOhms * 2.0 * ThermalVoltage);
	}

	static constexpr double StageBaseResistanceOhms()
	{
		// Each 220 ohm base shunt is loaded by the two 12.1 kohm signal
		// resistors connected to the local differential-pair node.
		return 1.0 / (1.0 / StageShuntResistanceOhms +
			2.0 / StageInputResistanceOhms);
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

	std::unique_ptr<ResamplerType> _audioInputResampler;
	std::unique_ptr<ResamplerType> _lowPassOutputResampler;
	std::unique_ptr<ResamplerType> _cutoffPitchResampler;
	std::unique_ptr<ResamplerType> _linearFmResampler;
	std::unique_ptr<ResamplerType> _resonanceResampler;
	std::unique_ptr<ResamplerType> _postOutputResampler;
	std::unique_ptr<ResamplerType> _postLinearCvResampler;
	std::unique_ptr<ResamplerType> _postExponentialCvResampler;
	std::array<double, 4> _state{};
	double _hostSampleRate{};
	double _sampleRate{};
	int _lastIterations{};
	std::size_t _solverFailures{};
	bool _renderedLowPassLastSample{true};

	struct OversampledControls
	{
		Eigen::Array<double, OversamplingFactor, 1> cutoffHz;
		Eigen::Array<double, OversamplingFactor, 1> resonance;
	};

	double DownsampleLowPass(
		const Eigen::Array<double, OversamplingFactor, 1>& input,
		bool render)
	{
		if (!render)
		{
			_renderedLowPassLastSample = false;
			return 0.0;
		}
		if (!_renderedLowPassLastSample)
			_lowPassOutputResampler->Reset();
		_renderedLowPassLastSample = true;
		return _lowPassOutputResampler->Downsample(input);
	}

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

	OversampledControls UpsampleControls(double log2CutoffHz,
		double linearFmHz, double resonance, double ceilingHz)
	{
		// Reconstruct cutoff in its exponential control domain. Mapping to hertz
		// after interpolation keeps audio-rate 1 V/octave modulation band-limited
		// before it changes the nonlinear solver coefficients.
		const auto cutoffPitch = _cutoffPitchResampler->Upsample(log2CutoffHz);
		const auto linearFm = _linearFmResampler->Upsample(linearFmHz);
		auto resonanceValues = _resonanceResampler->Upsample(resonance);
		OversampledControls controls;
		for (int i = 0; i < OversamplingFactor; ++i)
		{
			const double reconstructedHz = tfdsp::Exp2Taylor5(
				static_cast<float>(std::clamp(cutoffPitch(i), -100.0, 100.0)));
			controls.cutoffHz(i) = SoftLimitCutoff(reconstructedHz + linearFm(i),
				ceilingHz);
			controls.resonance(i) = std::clamp(resonanceValues(i), 0.0, 1.0);
		}
		return controls;
	}

	static double SoftOutputCompliance(double voltage)
	{
		const double magnitude = std::abs(voltage);
		if (magnitude <= OutputKneeVolts)
			return voltage;
		const double headroom = OutputRailVolts - OutputKneeVolts;
		return std::copysign(OutputKneeVolts + headroom * tfdsp::TanhPade76(
			(magnitude - OutputKneeVolts) / headroom), voltage);
	}

	double ProcessOversampled(double input, double cutoffHz, double resonance)
	{
		// The midpoint coefficient is prewarped so the small-signal one-pole
		// sections reach their requested analog cutoff at the oversampled rate.
		const double gamma = 2.0 * tfdsp::TanPrewarp(
			Pi * cutoffHz / _sampleRate);
		const double stageSaturationCoefficient =
			StageSaturationCoefficientPerVolt();
		const double stageStep = gamma / stageSaturationCoefficient;
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
			const double limiterTanh = tfdsp::TanhPade76(
				limiterVoltage / (2.0 * ThermalVoltage));
			const double firstInput = limiterPeak * limiterTanh;

			std::array<double, 4> stageTanh{};
			stageTanh[0] = tfdsp::TanhPade76(stageSaturationCoefficient *
				(firstInput + midpoint[0]));
			for (int i = 1; i < 4; ++i)
				stageTanh[i] = tfdsp::TanhPade76(stageSaturationCoefficient *
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

			for (double& value : residual)
				value = -value;
			if (!SolveCyclicBidiagonal4(
				1.0 + 0.5 * gamma * stageSlope[0],
				gamma * stageSlope[0] * firstInputDerivative,
				0.5 * gamma * stageSlope[1],
				1.0 + 0.5 * gamma * stageSlope[1],
				0.5 * gamma * stageSlope[2],
				1.0 + 0.5 * gamma * stageSlope[2],
				0.5 * gamma * stageSlope[3],
				1.0 + 0.5 * gamma * stageSlope[3], residual))
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

	/** Solve the fixed filter Jacobian without treating its twelve structural
	 * zeros as arbitrary dense entries. Rows 1..3 form a lower bidiagonal chain
	 * and the only wraparound term is row 0, column 3. Expressing each chained
	 * unknown as an affine function of x0 keeps the Newton system unchanged.
	 */
	static bool SolveCyclicBidiagonal4(double diagonal0, double wrap03,
		double lower10, double diagonal1, double lower21, double diagonal2,
		double lower32, double diagonal3, std::array<double, 4>& rhs)
	{
		constexpr double MinimumPivot = 1.0e-14;
		if (std::abs(diagonal1) < MinimumPivot ||
			std::abs(diagonal2) < MinimumPivot ||
			std::abs(diagonal3) < MinimumPivot)
			return false;

		const double intercept1 = rhs[1] / diagonal1;
		const double slope1 = -lower10 / diagonal1;
		const double intercept2 =
			(rhs[2] - lower21 * intercept1) / diagonal2;
		const double slope2 = -lower21 * slope1 / diagonal2;
		const double intercept3 =
			(rhs[3] - lower32 * intercept2) / diagonal3;
		const double slope3 = -lower32 * slope2 / diagonal3;
		const double reducedPivot = diagonal0 + wrap03 * slope3;
		if (std::abs(reducedPivot) < MinimumPivot)
			return false;

		rhs[0] = (rhs[0] - wrap03 * intercept3) / reducedPivot;
		rhs[1] = intercept1 + slope1 * rhs[0];
		rhs[2] = intercept2 + slope2 * rhs[0];
		rhs[3] = intercept3 + slope3 * rhs[0];
		return true;
	}
};

} // namespace tfdsp
