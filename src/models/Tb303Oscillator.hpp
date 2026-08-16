#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>

#include "tfdsp/approx.hpp"
#include "tfdsp/oscillator.hpp"
#include "tfdsp/sampleRate.hpp"

namespace tfdsp
{

/** Reduced model of the TB-303 Q8 saw-to-square wave shaper.
 *
 * The original common-emitter PNP stage is deliberately represented by its
 * dominant behavior rather than a per-sample transistor-network solve. The
 * emitter-bias state uses the 2.2 ms R104/C11 time constant. The collector
 * state uses the approximately 100 us R36/C10 time constant, with a faster
 * discharge that reproduces the visibly asymmetric edges. Its stock switching
 * point follows published hardware measurements: roughly 71% positive duty at
 * 10 Hz, 50% near 85 Hz, and 45% near 1 kHz.
 */
class Tb303SquareShaper
{
	double _sampleRate{192000.0};
	double _emitterBias{0.5};
	double _collector{};

	static double StockThreshold(double frequency)
	{
		const double octaves = std::log2(std::max(frequency, 0.01) / 85.0);
		const double threshold = octaves < 0.0 ?
			-0.1359 * octaves : -0.0281 * octaves;
		return std::clamp(threshold, -0.15, 0.45);
	}

public:
	void SetSampleRate(double sampleRate)
	{
		_sampleRate = std::max(sampleRate, 1.0);
	}

	void Reset()
	{
		_emitterBias = 0.5;
		_collector = 0.0;
	}

	double Step(double saw, double frequency, double shape)
	{
		if (!std::isfinite(saw) || !std::isfinite(frequency) ||
			!std::isfinite(shape))
			return 0.0;

		// SHAPE shifts the Q8 switching bias around the measured stock point.
		const double threshold = std::clamp(StockThreshold(frequency) +
			0.55 * std::clamp(shape, -1.0, 1.0), -0.82, 0.82);

		// Q8's emitter network makes the switching point depend slightly on its
		// recent conduction history. This produces the characteristic top droop
		// and makes fast pitch changes differ from a static tanh waveshaper.
		const double effectiveThreshold = threshold +
			0.12 * (0.5 - _emitterBias);
		const double conduction = 0.5 * (1.0 + std::tanh(
			(effectiveThreshold - saw) / 0.055));
		const double emitterAlpha = -std::expm1(
			-1.0 / (_sampleRate * 2.2e-3));
		_emitterBias += emitterAlpha * (conduction - _emitterBias);

		// C10 is connected around the collector/output node. Different charge
		// and discharge impedances give the hardware waveform one rounded edge
		// and one appreciably sharper edge.
		const double collectorTarget = 2.0 * conduction - 1.0;
		const double physicalCollectorTau = collectorTarget > _collector ?
			100.0e-6 : 22.0e-6;
		// Above the original oscillator's practical range, limit the edge time
		// to a fraction of one cycle. This retains a usable VCV audio oscillator
		// instead of allowing the vintage coupling network to erase the square.
		const double collectorTau = std::min(physicalCollectorTau,
			0.08 / std::max(frequency, 0.01));
		const double collectorAlpha = -std::expm1(
			-1.0 / (_sampleRate * collectorTau));
		_collector += collectorAlpha * (collectorTarget - _collector);

		// Remove the pitch-dependent DC component while preserving pulse-width
		// asymmetry. The Rack output gain is applied by Tb303Oscillator.
		return _collector - threshold;
	}
};

template<typename ResamplerType>
class Tb303Oscillator
{
public:
	using Resampler = ResamplerType;
	static constexpr int OversamplingFactor = ResamplerType::ResamplingFactor;

	struct Output
	{
		float saw{};
		float square{};
		float mixed{};
		float pitch{};
	};

private:
	std::unique_ptr<ResamplerType> _pitchInterpolator;
	std::unique_ptr<ResamplerType> _slideTimeInterpolator;
	std::unique_ptr<ResamplerType> _fmInterpolator;
	std::unique_ptr<ResamplerType> _shapeInterpolator;
	std::unique_ptr<ResamplerType> _waveInterpolator;
	std::unique_ptr<ResamplerType> _sawDecimator;
	std::unique_ptr<ResamplerType> _squareDecimator;
	std::unique_ptr<ResamplerType> _mixedDecimator;
	Tb303SquareShaper _squareShaper;
	tfdsp::BandlimitedSawOscillator<> _sawOscillator;
	double _sampleRate{48000.0};
	double _pitch{};
	bool _pitchInitialized{};

	static double LimitFrequency(double frequency, double sampleRate)
	{
		// Preserve exact tuning through the normal audio range, then approach the
		// discrete-time ceiling with a unity-slope soft knee. Linear FM remains
		// signed, so crossing zero reverses phase rather than pinning the VCO.
		const double knee = 0.40 * sampleRate;
		const double limit = 0.45 * sampleRate;
		const double magnitude = std::abs(frequency);
		if (magnitude <= knee)
			return frequency;
		const double curved = knee + (limit - knee) * std::tanh(
			(magnitude - knee) / (limit - knee));
		return std::copysign(curved, frequency);
	}

public:
	explicit Tb303Oscillator(
		std::function<std::unique_ptr<ResamplerType>()> createResampler) :
		_pitchInterpolator(createResampler()),
		_slideTimeInterpolator(createResampler()),
		_fmInterpolator(createResampler()),
		_shapeInterpolator(createResampler()),
		_waveInterpolator(createResampler()),
		_sawDecimator(createResampler()),
		_squareDecimator(createResampler()),
		_mixedDecimator(createResampler())
	{
		SetSampleRate(_sampleRate);
	}

	void SetSampleRate(double sampleRate)
	{
		_sampleRate = std::max(sampleRate, 1.0);
		_squareShaper.SetSampleRate(_sampleRate * OversamplingFactor);
	}

	void Reset()
	{
		_pitchInterpolator->Reset();
		_slideTimeInterpolator->Reset();
		_fmInterpolator->Reset();
		_shapeInterpolator->Reset();
		_waveInterpolator->Reset();
		_sawDecimator->Reset();
		_squareDecimator->Reset();
		_mixedDecimator->Reset();
		_squareShaper.Reset();
		_sawOscillator.Reset();
		_pitch = 0.0;
		_pitchInitialized = false;
	}

	Output Step(double targetPitch, bool slide, double slideTime,
		double tuningOffset, double fmVoltage, bool linearFm, double shape,
		double wave, double syncCrossing = -1.0)
	{
		if (!std::isfinite(targetPitch) || !std::isfinite(slideTime) ||
			!std::isfinite(tuningOffset) || !std::isfinite(fmVoltage) ||
			!std::isfinite(shape) || !std::isfinite(wave))
			return {};

		const double slideTimeLog = std::log(std::max(slideTime,
			std::numeric_limits<double>::min()));
		if (!_pitchInitialized)
		{
			_pitchInterpolator->PrimeUpsample(targetPitch);
			_slideTimeInterpolator->PrimeUpsample(slideTimeLog);
			_fmInterpolator->PrimeUpsample(fmVoltage);
			_shapeInterpolator->PrimeUpsample(shape);
			_waveInterpolator->PrimeUpsample(wave);
			_pitch = targetPitch;
			_pitchInitialized = true;
		}
		const auto targetPitchValues = _pitchInterpolator->Upsample(targetPitch);
		const auto slideTimeLogValues = _slideTimeInterpolator->Upsample(
			slideTimeLog);
		const auto fm = _fmInterpolator->Upsample(fmVoltage);
		const auto shapeValues = _shapeInterpolator->Upsample(shape);
		const auto waveValues = _waveInterpolator->Upsample(wave);
		Eigen::Array<double, OversamplingFactor, 1> sawValues;
		Eigen::Array<double, OversamplingFactor, 1> squareValues;
		Eigen::Array<double, OversamplingFactor, 1> mixedValues;

		const double internalRate = _sampleRate * OversamplingFactor;
		const auto syncEvent = tfdsp::MapEventToOversampledFrame<
			OversamplingFactor>(syncCrossing);
		for (int index = 0; index < OversamplingFactor; ++index)
		{
			if (slide)
			{
				const double tau = std::clamp(std::exp(slideTimeLogValues(index)) *
					(22.0 / 60.0), 0.0001, 2.0);
				const double slideAlpha = -std::expm1(
					-1.0 / (internalRate * tau));
				_pitch += slideAlpha * (targetPitchValues(index) - _pitch);
			}
			else
				_pitch = targetPitchValues(index);

			double pitch = _pitch + tuningOffset;
			double frequency;
			if (linearFm)
			{
				frequency = 261.625565 * tfdsp::Exp2Taylor5(
					static_cast<float>(std::clamp(pitch, -16.0, 16.0)));
				frequency += 200.0 * fm(index);
			}
			else
			{
				pitch += 0.2 * fm(index);
				frequency = 261.625565 * tfdsp::Exp2Taylor5(
					static_cast<float>(std::clamp(pitch, -16.0, 16.0)));
			}
			frequency = LimitFrequency(frequency, _sampleRate);
			const double phaseIncrement = frequency / internalRate;
			const double saw = _sawOscillator.Step(phaseIncrement,
				index == syncEvent.segment ? syncEvent.position : -1.0);
			const double square = _squareShaper.Step(saw, std::abs(frequency),
				std::clamp(shapeValues(index), -1.0, 1.0));
			const double blend = std::clamp(waveValues(index), 0.0, 1.0);
			const double sawRack = 5.0 * saw;
			// Original hardware square is about 4 Vpp versus 5.5 Vpp saw.
			// Q8's collector waveform is inverted relative to the saw. That
			// polarity is immaterial on the hardware selector, but retaining it in
			// a continuous morph creates a deep cancellation notch near the middle.
			const double squareRack = -(20.0 / 5.5) * square;
			sawValues(index) = sawRack;
			squareValues(index) = squareRack;
			mixedValues(index) = (1.0 - blend) * sawRack +
				blend * squareRack;
		}

		Output output;
		output.saw = static_cast<float>(_sawDecimator->Downsample(sawValues));
		output.square = static_cast<float>(
			_squareDecimator->Downsample(squareValues));
		output.mixed = static_cast<float>(
			_mixedDecimator->Downsample(mixedValues));
		output.pitch = static_cast<float>(_pitch);
		if (!std::isfinite(output.saw) || !std::isfinite(output.square) ||
			!std::isfinite(output.mixed) || !std::isfinite(output.pitch))
			return {};
		return output;
	}
};

} // namespace tfdsp
