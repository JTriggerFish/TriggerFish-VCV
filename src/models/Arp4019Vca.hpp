#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>

#include "tfdsp/rail.hpp"
#include "tfdsp/sampleRate.hpp"

namespace tfdsp
{

struct Arp4019ControlVoltages
{
	double linear{};
	double exponential{};
};

inline Arp4019ControlVoltages RouteArp4019ControlVoltages(
	double envelopeVolts, double modulationVolts, bool modulationConnected,
	bool addEnvelope, bool envelopeIsExponential,
	bool modulationIsExponential)
{
	Arp4019ControlVoltages routed;
	if (!modulationConnected || addEnvelope)
	{
		if (envelopeIsExponential)
			routed.exponential += envelopeVolts;
		else
			routed.linear += envelopeVolts;
	}
	if (modulationConnected)
	{
		if (modulationIsExponential)
			routed.exponential += modulationVolts;
		else
			routed.linear += modulationVolts;
	}
	return routed;
}

// Circuit-scaled reduction of the ARP 4019 VCA.
//
// The audio core is a PNP differential pair driven through the input
// attenuator and loaded by the LM301 transimpedance stage. The pair's tail
// current is the sum of the linear control converter and the exponential
// converter. This retains the mechanisms which matter at audio rate: smooth
// differential-pair compression, the 10 dB/V exponential law, finite control
// current and the original output-stage bandwidth.
template<typename ResamplerType>
class Arp4019Vca
{
public:
	static constexpr int OversamplingFactor = ResamplerType::ResamplingFactor;

	explicit Arp4019Vca(
		std::function<std::unique_ptr<ResamplerType>()> resamplerCreator)
		: _audioResampler(resamplerCreator()),
		  _linearCvResampler(resamplerCreator()),
		  _exponentialCvResampler(resamplerCreator())
	{
	}

	void SetSampleRate(double sampleRate)
	{
		if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
			return;
		_hostSampleRate = sampleRate;
		_sampleRate = sampleRate * OversamplingFactor;
		_outputCoefficient = -std::expm1(-2.0 * Pi *
			OutputBandwidthHz / _sampleRate);
		Reset();
	}

	void Reset()
	{
		_audioResampler->Reset();
		_linearCvResampler->Reset();
		_exponentialCvResampler->Reset();
		_outputLowPass = 0.0;
		_lastControlGain = 0.0;
		_lastControlCurrent = 0.0;
	}

	float Step(double noninvertingInputVolts, double invertingInputVolts,
		double linearControlVolts, double exponentialControlVolts,
		double initialGain = 0.0)
	{
		if (!std::isfinite(noninvertingInputVolts) ||
			!std::isfinite(invertingInputVolts) ||
			!std::isfinite(linearControlVolts) ||
			!std::isfinite(exponentialControlVolts) ||
			!std::isfinite(initialGain) || !(_sampleRate > 0.0))
		{
			Reset();
			return 0.0f;
		}

		const auto audio = _audioResampler->Upsample(
			noninvertingInputVolts - invertingInputVolts);
		const auto linearCv = _linearCvResampler->Upsample(linearControlVolts);
		const auto exponentialCv = _exponentialCvResampler->Upsample(
			exponentialControlVolts);
		Eigen::Array<double, OversamplingFactor, 1> output;
		for (int i = 0; i < OversamplingFactor; ++i)
			output(i) = RackOutputAdapter::ProcessOversampled(
				ProcessOversampled(audio(i), linearCv(i), exponentialCv(i),
					initialGain));

		const double result = _audioResampler->Downsample(output);
		if (!std::isfinite(result))
		{
			Reset();
			return 0.0f;
		}
		return static_cast<float>(
			RackOutputAdapter::ProcessPostDecimation(result));
	}

	// Process one sample already running at the module's internal rate. This is
	// used when the VCA follows another oversampled nonlinear block, avoiding a
	// redundant downsample/upsample boundary.
	double ProcessOversampled(double audioDifferenceVolts,
		double linearControlVolts, double exponentialControlVolts,
		double initialGain = 0.0)
	{
		if (!std::isfinite(audioDifferenceVolts) ||
			!std::isfinite(linearControlVolts) ||
			!std::isfinite(exponentialControlVolts) ||
			!std::isfinite(initialGain))
		{
			_outputLowPass = 0.0;
			_lastControlGain = 0.0;
			_lastControlCurrent = 0.0;
			return 0.0;
		}

		const double linearGain = SoftPositive(initialGain +
			linearControlVolts / LinearUnityControlVolts);
		const double exponentialGain = ExponentialGain(
			exponentialControlVolts);
		_lastControlGain = LimitControlGain(linearGain + exponentialGain);
		_lastControlCurrent = UnityControlCurrentAmps() * _lastControlGain;

		const double differentialBaseVolts = AudioInputScale() *
			audioDifferenceVolts;
		const double outputCurrent = _lastControlCurrent * std::tanh(
			differentialBaseVolts / (2.0 * ThermalVoltage));
		const double target = OutputFeedbackResistanceOhms * outputCurrent;
		_outputLowPass += _outputCoefficient * (target - _outputLowPass);
		return SoftOutputCompliance(_outputLowPass);
	}

	double LastControlGain() const { return _lastControlGain; }
	double LastControlCurrent() const { return _lastControlCurrent; }

	static constexpr double ThermalVoltage = 25.85e-3;
	static constexpr double Pi = 3.14159265358979323846;
	static constexpr double AudioInputResistanceOhms = 100000.0;
	static constexpr double AudioSeriesResistanceOhms = 1000.0;
	static constexpr double AudioInputShuntOhms = 220.0;
	static constexpr double OutputFeedbackResistanceOhms = 56000.0;
	static constexpr double OutputFeedbackCapacitanceFarads = 100.0e-12;
	static constexpr double LinearUnityControlVolts = 10.0;
	static constexpr double ExponentialDecibelsPerVolt = 10.0;

	static constexpr double AudioInputScale()
	{
		return AudioInputShuntOhms /
			(AudioInputResistanceOhms + AudioSeriesResistanceOhms +
				AudioInputShuntOhms);
	}

	static constexpr double UnityControlCurrentAmps()
	{
		return 2.0 * ThermalVoltage /
			(OutputFeedbackResistanceOhms * AudioInputScale());
	}

	static constexpr double SmallSignalGainAtUnityControl()
	{
		return OutputFeedbackResistanceOhms * UnityControlCurrentAmps() *
			AudioInputScale() / (2.0 * ThermalVoltage);
	}

	static constexpr double OutputBandwidthHz = 1.0 /
		(2.0 * Pi * OutputFeedbackResistanceOhms *
			OutputFeedbackCapacitanceFarads);

private:
	static constexpr double LogTen = 2.30258509299404568402;
	static constexpr double LinearCutoffKnee = 1.0e-6;
	static constexpr double MaximumControlGain = 16.0;
	static constexpr double OutputKneeVolts = 10.0;
	static constexpr double OutputRailVolts = 13.5;

	std::unique_ptr<ResamplerType> _audioResampler;
	std::unique_ptr<ResamplerType> _linearCvResampler;
	std::unique_ptr<ResamplerType> _exponentialCvResampler;
	double _hostSampleRate{};
	double _sampleRate{};
	double _outputCoefficient{};
	double _outputLowPass{};
	double _lastControlGain{};
	double _lastControlCurrent{};

	static double SoftPositive(double value)
	{
		const double normalized = value / LinearCutoffKnee;
		if (normalized > 40.0)
			return value;
		if (normalized < -40.0)
			return LinearCutoffKnee * std::exp(normalized);
		return LinearCutoffKnee * std::log1p(std::exp(normalized));
	}

	static double ExponentialGain(double controlVolts)
	{
		// The owner's manual specifies 10 dB/V and a 100 dB range. The
		// resulting calibration places unity at 10 V and the residual gain at
		// 0 V at 1e-5. Extreme external CV is limited later by the finite tail
		// current rather than clipped at this input.
		const double gainExponent = LogTen * 0.5 *
			(controlVolts - LinearUnityControlVolts);
		if (gainExponent < -50.0)
			return std::exp(-50.0);
		if (gainExponent > 50.0)
			return std::exp(50.0);
		return std::exp(gainExponent);
	}

	static double LimitControlGain(double gain)
	{
		if (!(gain > 0.0))
			return 0.0;
		const double ratio = gain / MaximumControlGain;
		return gain / std::pow(1.0 + ratio * ratio * ratio * ratio, 0.25);
	}

	static double SoftOutputCompliance(double voltage)
	{
		if (!std::isfinite(voltage))
			return 0.0;
		const double magnitude = std::abs(voltage);
		if (magnitude <= OutputKneeVolts)
			return voltage;
		const double headroom = OutputRailVolts - OutputKneeVolts;
		return std::copysign(OutputKneeVolts + headroom * std::tanh(
			(magnitude - OutputKneeVolts) / headroom), voltage);
	}

};

} // namespace tfdsp
