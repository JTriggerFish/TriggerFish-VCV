#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>

#include "tfdsp/rail.hpp"
#include "tfdsp/approx.hpp"
#include "tfdsp/oscillator.hpp"
#include "tfdsp/sampleRate.hpp"

namespace tfdsp
{

/** Circuit-scaled model of the TB-303 Q8 saw-to-square wave shaper.
 *
 * R34/R35/C10 drive the base, R45/C11 bias the emitter, and R36 loads the
 * collector. The two capacitors use trapezoidal companion models. A compact
 * two-junction Ebers-Moll solve captures cutoff, forward-active operation, and
 * saturation of the original 2SA733P.
 */
class Tb303SquareShaper
{
	double _sampleRate{192000.0};
	double _c10Voltage{};
	double _c10Current{};
	double _emitterVoltage{8.4};
	double _c11Current{};
	double _forwardJunctionExponent{10.12667110305036};
	double _reverseJunctionExponent{-10.0};

	static constexpr double SupplyVoltage = 12.0;
	static constexpr double BiasVoltage = 5.333;
	static constexpr double CollectorReferenceVoltage = 6.8;
	static constexpr double SawBiasVoltage = 8.5;
	static constexpr double SawHalfRangeVoltage = 2.75;
	static constexpr double ShapeRangeVoltage = 3.0;
	static constexpr double R34 = 10.0e3;
	static constexpr double R35 = 100.0e3;
	static constexpr double R36 = 10.0e3;
	static constexpr double R45 = 22.0e3;
	static constexpr double C10 = 10.0e-9;
	static constexpr double C11 = 1.0e-6;
	static constexpr double SaturationCurrent = 4.0e-14;
	static constexpr double ForwardBeta = 300.0;
	static constexpr double ReverseBeta = 10.0;
	static constexpr double ForwardAlpha = ForwardBeta / (ForwardBeta + 1.0);
	static constexpr double ReverseAlpha = ReverseBeta / (ReverseBeta + 1.0);
	static constexpr double ReverseSaturationCurrent =
		ForwardAlpha / ReverseAlpha * SaturationCurrent;
	static constexpr double ThermalVoltage = 8.617333262e-5 * (273.15 + 27.0);
	static constexpr double StockMaximumFrequency = 1000.0;
	static constexpr double ExtendedSquareFullFrequency = 2000.0;
	static constexpr int MaximumNewtonIterations = 8;
	static constexpr double Log2E = 1.4426950408889634;

	struct JunctionState
	{
		double forwardCurrent{};
		double reverseCurrent{};
		double emitterCurrent{};
		double collectorCurrent{};
		double baseCurrent{};
		double baseVoltage{};
		double emitterVoltage{};
		double collectorVoltage{};
		double forwardDerivative{};
		double reverseDerivative{};
	};

	static JunctionState EvaluateJunctions(double forwardExponent,
		double reverseExponent, double baseOpenVoltage, double baseResistance,
		double emitterOpenVoltage, double emitterResistance)
	{
		const double forwardExponential = tfdsp::Exp2Taylor5(
			static_cast<float>(std::clamp(forwardExponent, -30.0, 40.0) * Log2E));
		const double reverseExponential = tfdsp::Exp2Taylor5(
			static_cast<float>(std::clamp(reverseExponent, -30.0, 40.0) * Log2E));
		JunctionState state;
		state.forwardCurrent = SaturationCurrent *
			(forwardExponential - 1.0);
		state.reverseCurrent = ReverseSaturationCurrent *
			(reverseExponential - 1.0);
		state.forwardDerivative = SaturationCurrent * forwardExponential;
		state.reverseDerivative = ReverseSaturationCurrent * reverseExponential;
		state.emitterCurrent = state.forwardCurrent -
			ReverseAlpha * state.reverseCurrent;
		state.collectorCurrent = ForwardAlpha * state.forwardCurrent -
			state.reverseCurrent;
		state.baseCurrent = (1.0 - ForwardAlpha) * state.forwardCurrent +
			(1.0 - ReverseAlpha) * state.reverseCurrent;
		state.baseVoltage = baseOpenVoltage +
			baseResistance * state.baseCurrent;
		state.emitterVoltage = emitterOpenVoltage -
			emitterResistance * state.emitterCurrent;
		state.collectorVoltage = BiasVoltage +
			R36 * state.collectorCurrent;
		return state;
	}

public:
	void SetSampleRate(double sampleRate)
	{
		_sampleRate = std::max(sampleRate, 1.0);
	}

	void Reset()
	{
		_c10Voltage = 0.0;
		_c10Current = 0.0;
		_emitterVoltage = 8.4;
		_c11Current = 0.0;
		_forwardJunctionExponent =
			std::log(1.0 + 1.0e-9 / SaturationCurrent);
		_reverseJunctionExponent = -10.0;
	}

	double Step(double saw, double frequency, double shape)
	{
		if (!std::isfinite(saw) || !std::isfinite(frequency) ||
			!std::isfinite(shape))
			return 0.0;

		const double samplePeriod = 1.0 / _sampleRate;
		const double shapeControl = std::clamp(shape, -1.0, 1.0);
		const double sawVoltage = SawBiasVoltage -
			SawHalfRangeVoltage * saw - ShapeRangeVoltage *
			shapeControl;

		// Trapezoidal companion for the series C10/R34 path.
		const double c10Resistance = samplePeriod / (2.0 * C10);
		const double c10History = _c10Voltage +
			c10Resistance * _c10Current;
		const double acBranchResistance = R34 + c10Resistance;
		const double baseResistance = 1.0 /
			(1.0 / R35 + 1.0 / acBranchResistance);
		const double baseOpenVoltage = baseResistance *
			(sawVoltage / R35 +
				(sawVoltage - c10History) / acBranchResistance);

		// Trapezoidal companion for the C11 emitter-bias network.
		const double c11Conductance = 2.0 * C11 / samplePeriod;
		const double c11History = -_c11Current -
			c11Conductance * _emitterVoltage;
		const double emitterResistance = 1.0 /
			(1.0 / R45 + c11Conductance);
		const double emitterOpenVoltage = emitterResistance *
			(SupplyVoltage / R45 - c11History);

		double forwardExponent = _forwardJunctionExponent;
		double reverseExponent = _reverseJunctionExponent;
		JunctionState junction;
		for (int iteration = 0; iteration < MaximumNewtonIterations; ++iteration)
		{
			junction = EvaluateJunctions(forwardExponent, reverseExponent,
				baseOpenVoltage, baseResistance, emitterOpenVoltage,
				emitterResistance);
			const double residualForward = forwardExponent -
				(junction.emitterVoltage - junction.baseVoltage) /
					ThermalVoltage;
			const double residualReverse = reverseExponent -
				(junction.collectorVoltage - junction.baseVoltage) /
					ThermalVoltage;

			const double emitterForward = junction.forwardDerivative;
			const double emitterReverse = -ReverseAlpha *
				junction.reverseDerivative;
			const double collectorForward = ForwardAlpha *
				junction.forwardDerivative;
			const double collectorReverse = -junction.reverseDerivative;
			const double baseForward = (1.0 - ForwardAlpha) *
				junction.forwardDerivative;
			const double baseReverse = (1.0 - ReverseAlpha) *
				junction.reverseDerivative;

			const double j00 = 1.0 +
				(emitterResistance * emitterForward +
					baseResistance * baseForward) / ThermalVoltage;
			const double j01 =
				(emitterResistance * emitterReverse +
					baseResistance * baseReverse) / ThermalVoltage;
			const double j10 = -
				(R36 * collectorForward - baseResistance * baseForward) /
					ThermalVoltage;
			const double j11 = 1.0 -
				(R36 * collectorReverse - baseResistance * baseReverse) /
					ThermalVoltage;
			const double determinant = j00 * j11 - j01 * j10;
			if (!std::isfinite(determinant) || std::abs(determinant) < 1.0e-20)
				break;
			double deltaForward =
				(residualForward * j11 - residualReverse * j01) / determinant;
			double deltaReverse =
				(j00 * residualReverse - j10 * residualForward) / determinant;
			const double damping = std::max({1.0,
				std::abs(deltaForward) / 2.0, std::abs(deltaReverse) / 2.0});
			deltaForward /= damping;
			deltaReverse /= damping;
			forwardExponent -= deltaForward;
			reverseExponent -= deltaReverse;
			if (std::max(std::abs(deltaForward), std::abs(deltaReverse)) < 1.0e-4)
				break;
		}

		junction = EvaluateJunctions(forwardExponent, reverseExponent,
			baseOpenVoltage, baseResistance, emitterOpenVoltage,
			emitterResistance);
		if (!std::isfinite(junction.collectorVoltage) ||
			!std::isfinite(junction.emitterVoltage) ||
			!std::isfinite(junction.baseVoltage))
		{
			Reset();
			return 0.0;
		}
		_forwardJunctionExponent = forwardExponent;
		_reverseJunctionExponent = reverseExponent;

		_c10Current = (sawVoltage - c10History - junction.baseVoltage) /
			acBranchResistance;
		_c10Voltage = c10History + c10Resistance * _c10Current;
		_c11Current = (SupplyVoltage - junction.emitterVoltage) / R45 -
			junction.emitterCurrent;
		_emitterVoltage = junction.emitterVoltage;

		// Preserve duty-dependent DC movement with a fixed circuit-to-Rack
		// reference. The hardware C17/R62 coupling network is modeled at the input
		// of the 303 Voice Core, where it appears in the original signal path.
		const double circuitOutput =
			0.5 * (junction.collectorVoltage - CollectorReferenceVoltage);

		// Q8's two RC networks progressively narrow and attenuate the pulse above
		// the original oscillator's useful range. The module keeps that circuit
		// response through 1 kHz, then crosses smoothly to an oversampled
		// Schmitt-like extension so that its additional high octaves remain useful.
		const double extensionPosition = std::clamp(
			std::log2(std::max(frequency, StockMaximumFrequency) /
				StockMaximumFrequency) /
				std::log2(ExtendedSquareFullFrequency /
					StockMaximumFrequency),
			0.0, 1.0);
		const double extensionBlend = extensionPosition * extensionPosition *
			(3.0 - 2.0 * extensionPosition);
		// Keep the stock switching point and the existing centre sensitivity,
		// while curving the end stops to +/-0.8 of the normalized saw range.
		// This gives the extended comparator roughly 10% to 90% duty instead
		// of allowing the negative Shape end to move beyond the saw maximum.
		const double extensionThreshold = 0.36 - 0.8 * shapeControl -
			0.36 * shapeControl * shapeControl;
		const double extensionOutput =
			std::tanh((saw - extensionThreshold) / 0.055);
		return circuitOutput +
			extensionBlend * (extensionOutput - circuitOutput);
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
			// The physical saw is mapped with descending polarity before Q8, so its
			// collector output is already phase-aligned with the Rack-facing saw.
			const double squareRack = (20.0 / 5.5) * square;
			sawValues(index) = RackOutputAdapter::ProcessOversampled(sawRack);
			squareValues(index) = RackOutputAdapter::ProcessOversampled(squareRack);
			mixedValues(index) = RackOutputAdapter::ProcessOversampled(
				(1.0 - blend) * sawRack + blend * squareRack);
		}

		Output output;
		output.saw = static_cast<float>(RackOutputAdapter::ProcessPostDecimation(
			_sawDecimator->Downsample(sawValues)));
		output.square = static_cast<float>(
			RackOutputAdapter::ProcessPostDecimation(
				_squareDecimator->Downsample(squareValues)));
		output.mixed = static_cast<float>(
			RackOutputAdapter::ProcessPostDecimation(
				_mixedDecimator->Downsample(mixedValues)));
		output.pitch = static_cast<float>(_pitch);
		if (!std::isfinite(output.saw) || !std::isfinite(output.square) ||
			!std::isfinite(output.mixed) || !std::isfinite(output.pitch))
			return {};
		return output;
	}
};

} // namespace tfdsp
