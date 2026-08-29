#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <utility>

#include "tfdsp/oscillator.hpp"
#include "tfdsp/sampleRate.hpp"

namespace tfdsp
{

enum class WavefolderCharacter
{
	Hinge,
	Lockhart,
	Serge,
	Count
};

namespace wavefolder_detail
{
	/** Principal branch W0(x) for finite non-negative x.
	 *
	 * Wavefolder tables are generated once during initialization, so a robust
	 * Halley iteration is preferable to a more complicated audio-rate
	 * approximation. Small arguments use the local series to retain precision.
	 */
	inline double PrincipalLambertW(double input)
	{
		if (!(input > 0.0))
			return 0.0;
		if (input < 1.0e-8)
		{
			const double square = input * input;
			return input - square + 1.5 * square * input;
		}

		double value;
		if (input < 3.0)
			value = std::log1p(input);
		else
		{
			const double first = std::log(input);
			const double second = std::log(first);
			value = first - second + second / first;
		}
		for (int iteration = 0; iteration < 8; ++iteration)
		{
			const double exponential = std::exp(value);
			const double residual = value * exponential - input;
			const double denominator = exponential * (value + 1.0) -
				(value + 2.0) * residual / (2.0 * value + 2.0);
			const double correction = residual / denominator;
			value -= correction;
			if (std::abs(correction) <= 2.0e-15 *
				std::max(1.0, std::abs(value)))
				break;
		}
		return value;
	}
}

/** Selectable, odd-symmetric wavefolder with tabulated first-order ADAA.
 *
 * Hinge is a generic four-stage cascade of smooth square-root hinges. Lockhart
 * and Serge implement the four- and six-stage circuit models in Esqueda et al.,
 * "Virtual Analog Models of the Lockhart and Serge Wavefolders" (2017),
 * https://doi.org/10.3390/app7121328. Each fixed cascade and its derivative are
 * sampled once, then represented as a cubic Hermite curve whose integral is
 * evaluated analytically. Fold and symmetry remain external gain and offset
 * controls, so ordinary one-dimensional ADAA applies to every character.
 */
class Wavefolder
{
public:
	static constexpr int HingeStageCount = 4;
	static constexpr int TableIntervals = 4096;
	static constexpr double TableMinimum = -16.0;
	static constexpr double TableMaximum = 16.0;

	Wavefolder()
	{
		PrepareTable();
	}

	static void PrepareTable()
	{
		(void) Table();
	}

	void Reset()
	{
		_initialized = false;
		_previousInput = 0.0;
		_previousPrimitive = 0.0;
	}

	double Process(double input,
		WavefolderCharacter character = WavefolderCharacter::Hinge)
	{
		if (!std::isfinite(input))
		{
			Reset();
			return 0.0;
		}

		const Evaluation current = Evaluate(input, character);
		if (!_initialized || character != _character)
		{
			_initialized = true;
			_character = character;
			_previousInput = input;
			_previousPrimitive = current.primitive;
			return current.value;
		}

		const double difference = input - _previousInput;
		const double scale = std::max({1.0, std::abs(input),
			std::abs(_previousInput)});
		double output;
		if (std::abs(difference) > 1.0e-7 * scale)
			output = (current.primitive - _previousPrimitive) / difference;
		else
			output = Evaluate(0.5 * (input + _previousInput), character).value;

		_previousInput = input;
		_previousPrimitive = current.primitive;
		if (std::isfinite(output))
			return output;
		Reset();
		return 0.0;
	}

	static double Transfer(double input,
		WavefolderCharacter character = WavefolderCharacter::Hinge)
	{
		return Evaluate(input, character).value;
	}

	static double Primitive(double input,
		WavefolderCharacter character = WavefolderCharacter::Hinge)
	{
		return Evaluate(input, character).primitive;
	}

private:
	struct Stage
	{
		double threshold;
		double softness;
	};

	struct Node
	{
		double value{};
		double derivative{};
		double primitive{};
	};

	struct Evaluation
	{
		double value{};
		double primitive{};
	};

	using LookupTable = std::array<Node, TableIntervals + 1>;

	bool _initialized{};
	WavefolderCharacter _character{WavefolderCharacter::Hinge};
	double _previousInput{};
	double _previousPrimitive{};

	static constexpr std::array<Stage, HingeStageCount> Stages{{
		{1.000, 0.055},
		{0.965, 0.065},
		{1.035, 0.075},
		{0.985, 0.085},
	}};

	static double SmoothHinge(double value, double softness)
	{
		const double radius = std::hypot(value, softness);
		if (value >= 0.0)
			return 0.5 * (value + radius);
		return 0.5 * softness * softness / (radius - value);
	}

	static double SmoothHingeDerivative(double value, double softness)
	{
		return 0.5 * (1.0 + value / std::hypot(value, softness));
	}

	static std::pair<double, double> HingeTransferAndDerivative(double input)
	{
		double value = input;
		double derivative = 1.0;
		for (const Stage& stage : Stages)
		{
			const double positive = value - stage.threshold;
			const double negative = -value - stage.threshold;
			const double stageDerivative = 1.0 -
				2.0 * SmoothHingeDerivative(positive, stage.softness) -
				2.0 * SmoothHingeDerivative(negative, stage.softness);
			value = value - 2.0 * SmoothHinge(positive, stage.softness) +
				2.0 * SmoothHinge(negative, stage.softness);
			derivative *= stageDerivative;
		}

		// Restore unity small-signal gain after the small losses caused by the
		// soft tails of the four folding stages.
		static const double normalization = []
		{
			double centralDerivative = 1.0;
			double centralValue = 0.0;
			for (const Stage& stage : Stages)
			{
				const double positive = centralValue - stage.threshold;
				const double negative = -centralValue - stage.threshold;
				centralDerivative *= 1.0 -
					2.0 * SmoothHingeDerivative(positive, stage.softness) -
					2.0 * SmoothHingeDerivative(negative, stage.softness);
				centralValue = centralValue -
					2.0 * SmoothHinge(positive, stage.softness) +
					2.0 * SmoothHinge(negative, stage.softness);
			}
			return 1.0 / centralDerivative;
		}();
		return {normalization * value, normalization * derivative};
	}

	static std::pair<double, double> LockhartStage(double input)
	{
		constexpr double ThermalVoltage = 0.025864;
		constexpr double Resistance = 15000.0;
		constexpr double LoadResistance = 7500.0;
		constexpr double SaturationCurrent = 1.0e-17;
		constexpr double Alpha = 2.0 * LoadResistance / Resistance;
		constexpr double Beta = (2.0 * LoadResistance + Resistance) /
			(ThermalVoltage * Resistance);
		constexpr double Delta = LoadResistance * SaturationCurrent /
			ThermalVoltage;
		const double polarity = input >= 0.0 ? 1.0 : -1.0;
		const double w = wavefolder_detail::PrincipalLambertW(
			Delta * std::exp(polarity * Beta * input));
		const double value = Alpha * input -
			polarity * ThermalVoltage * w;
		const double derivative = Alpha -
			ThermalVoltage * Beta * w / (1.0 + w);
		return {value, derivative};
	}

	static std::pair<double, double> SergeStage(double input)
	{
		constexpr double ThermalVoltage = 0.025864;
		constexpr double EmissionCoefficient = 1.752;
		constexpr double Resistance = 33000.0;
		constexpr double SaturationCurrent = 2.52e-9;
		constexpr double JunctionVoltage =
			EmissionCoefficient * ThermalVoltage;
		constexpr double SeriesDrop = Resistance * SaturationCurrent;
		constexpr double Delta = Resistance * SaturationCurrent /
			JunctionVoltage;
		const double polarity = input >= 0.0 ? 1.0 : -1.0;
		const double w = wavefolder_detail::PrincipalLambertW(
			Delta * std::exp((polarity * input + SeriesDrop) /
				JunctionVoltage));
		// The paper simplifies the diode law before its equation (39). Retaining
		// the 83.16 uV series term from the circuit equation makes the transfer
		// continuous and odd through zero.
		const double value = input - 2.0 * polarity *
			(JunctionVoltage * w - SeriesDrop);
		const double derivative = (1.0 - w) / (1.0 + w);
		return {value, derivative};
	}

	template<int StagesInCascade, typename StageFunction>
	static std::pair<double, double> Cascade(double input,
		StageFunction&& stageFunction)
	{
		double value = input;
		double derivative = 1.0;
		for (int stage = 0; stage < StagesInCascade; ++stage)
		{
			const auto [next, stageDerivative] = stageFunction(value);
			value = next;
			derivative *= stageDerivative;
		}
		return {value, derivative};
	}

	static std::pair<double, double> LockhartTransferAndDerivative(double input)
	{
		// Equation (28), arranged as the paper's four-stage topology with its
		// 1/3 input and 3x output gains. Global scales align the first fold to
		// normalized input 1 and restore unity small-signal gain.
		constexpr double InputScale = 1.069465474121151;
		constexpr double OutputScale = 0.9350465482253767;
		const auto [value, derivative] = Cascade<4>(
			input * InputScale / 3.0, LockhartStage);
		return {3.0 * OutputScale * value,
			OutputScale * InputScale * derivative};
	}

	static std::pair<double, double> SergeTransferAndDerivative(double input)
	{
		// The circuit transfer preceding equation (39), cascaded through the six
		// stages of the Serge middle wave multiplier. The paper's 4x output gain
		// is included before global input and output calibration.
		constexpr double InputScale = 0.33073419047751157;
		constexpr double OutputScale = 0.7727253519539319;
		const auto [value, derivative] = Cascade<6>(
			input * InputScale, SergeStage);
		return {4.0 * OutputScale * value,
			4.0 * OutputScale * InputScale * derivative};
	}

	static std::pair<double, double> DirectTransferAndDerivative(double input,
		WavefolderCharacter character)
	{
		switch (character)
		{
		case WavefolderCharacter::Lockhart:
			return LockhartTransferAndDerivative(input);
		case WavefolderCharacter::Serge:
			return SergeTransferAndDerivative(input);
		default:
			return HingeTransferAndDerivative(input);
		}
	}

	static LookupTable BuildTable(WavefolderCharacter character)
	{
		LookupTable table{};
		constexpr double step =
			(TableMaximum - TableMinimum) / TableIntervals;
		for (int index = 0; index <= TableIntervals; ++index)
		{
			const double input = TableMinimum + step * index;
			const auto [value, derivative] =
				DirectTransferAndDerivative(input, character);
			table[index].value = value;
			table[index].derivative = derivative;
		}

		table[0].primitive = 0.0;
		for (int index = 0; index < TableIntervals; ++index)
		{
			const Node& left = table[index];
			const Node& right = table[index + 1];
			const double area = step * (0.5 * (left.value + right.value) +
				step * (left.derivative - right.derivative) / 12.0);
			table[index + 1].primitive = left.primitive + area;
		}
		return table;
	}

	static const LookupTable& Table(WavefolderCharacter character =
		WavefolderCharacter::Hinge)
	{
		static const std::array<LookupTable,
			static_cast<std::size_t>(WavefolderCharacter::Count)> tables{{
			BuildTable(WavefolderCharacter::Hinge),
			BuildTable(WavefolderCharacter::Lockhart),
			BuildTable(WavefolderCharacter::Serge),
		}};
		const std::size_t index = std::min(static_cast<std::size_t>(character),
			tables.size() - 1);
		return tables[index];
	}

	static Evaluation Evaluate(double input, WavefolderCharacter character)
	{
		const LookupTable& table = Table(character);
		constexpr double step =
			(TableMaximum - TableMinimum) / TableIntervals;

		// The cascade is already in its final nearly-linear outer branch at the
		// table boundaries. Continue with the boundary tangent rather than
		// imposing a digital clamp if an unexpected input exceeds the table.
		if (input <= TableMinimum)
		{
			const Node& node = table.front();
			const double distance = input - TableMinimum;
			return {node.value + node.derivative * distance,
				node.primitive + node.value * distance +
					0.5 * node.derivative * distance * distance};
		}
		if (input >= TableMaximum)
		{
			const Node& node = table.back();
			const double distance = input - TableMaximum;
			return {node.value + node.derivative * distance,
				node.primitive + node.value * distance +
					0.5 * node.derivative * distance * distance};
		}

		const double position = (input - TableMinimum) / step;
		const int index = std::clamp(static_cast<int>(std::floor(position)),
			0, TableIntervals - 1);
		const double t = position - index;
		const double t2 = t * t;
		const double t3 = t2 * t;
		const double t4 = t2 * t2;
		const Node& left = table[index];
		const Node& right = table[index + 1];

		const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
		const double h10 = t3 - 2.0 * t2 + t;
		const double h01 = -2.0 * t3 + 3.0 * t2;
		const double h11 = t3 - t2;
		const double value = h00 * left.value + h10 * step * left.derivative +
			h01 * right.value + h11 * step * right.derivative;

		const double i00 = 0.5 * t4 - t3 + t;
		const double i10 = 0.25 * t4 - (2.0 / 3.0) * t3 + 0.5 * t2;
		const double i01 = -0.5 * t4 + t3;
		const double i11 = 0.25 * t4 - (1.0 / 3.0) * t3;
		const double primitive = left.primitive + step *
			(i00 * left.value + i10 * step * left.derivative +
				i01 * right.value + i11 * step * right.derivative);
		return {value, primitive};
	}
};

struct WavefoldOscillatorOutput
{
	double oscillator{};
	double folded{};
};

/** Triangle/sine morph oscillator feeding the selectable wavefolder. */
template<typename ResamplerType>
class WavefoldOscillator
{
public:
	static constexpr int OversamplingFactor = ResamplerType::ResamplingFactor;
	static constexpr double FoldHarmonicBudgetHz = 6000.0;
	using Output = WavefoldOscillatorOutput;

	explicit WavefoldOscillator(
		std::function<std::unique_ptr<ResamplerType>()> createResampler) :
		_frequencyInterpolator(createResampler()),
		_morphInterpolator(createResampler()),
		_foldInterpolator(createResampler()),
		_symmetryInterpolator(createResampler()),
		_externalInputInterpolator(createResampler()),
		_oscillatorDecimator(createResampler()),
		_foldedDecimator(createResampler())
	{
		SetSampleRate(_sampleRate);
	}

	void SetSampleRate(double sampleRate)
	{
		_sampleRate = std::max(sampleRate, 1.0);
	}

	void Reset()
	{
		_frequencyInterpolator->Reset();
		_morphInterpolator->Reset();
		_foldInterpolator->Reset();
		_symmetryInterpolator->Reset();
		_externalInputInterpolator->Reset();
		_oscillatorDecimator->Reset();
		_foldedDecimator->Reset();
		_triangle.Reset();
		_folder.Reset();
		_previousFolderSource = 0.0;
		_folderSourceInitialized = false;
		_controlsInitialized = false;
		_externalInputInitialized = false;
		_renderedOscillatorLastSample = true;
		_renderedFoldedLastSample = true;
	}

	void SetFolderAntialiasing(bool enabled)
	{
		if (enabled != _useAdaa)
		{
			_useAdaa = enabled;
			_folder.Reset();
			_folderSourceInitialized = false;
		}
	}

	void SetCharacter(WavefolderCharacter character)
	{
		if (character != _character)
		{
			_character = character;
			_folder.Reset();
		}
	}

	WavefolderCharacter Character() const { return _character; }

	/** Reduce fold when its approximate harmonic density would exceed a fixed
	 * musical-frequency budget. The mapping is independent of host sample rate,
	 * so changing the Rack engine rate does not change the patch voicing.
	 */
	static double FoldScaleForFrequency(double frequencyHz)
	{
		const double frequency = std::abs(frequencyHz);
		if (!(frequency > 0.0) ||
			frequency <= FoldHarmonicBudgetHz / 9.0)
			return 1.0;
		if (frequency >= FoldHarmonicBudgetHz)
			return 0.0;
		return (FoldHarmonicBudgetHz / frequency - 1.0) / 8.0;
	}

	double Step(double frequencyHz, double morph, double fold,
		double symmetry)
	{
		return StepWithInput(frequencyHz, morph, fold, symmetry, 0.0,
			false).folded;
	}

	/** Render the oscillator and folded paths together.
	 *
	 * When useExternalInput is false, the internal oversampled oscillator feeds
	 * the folder directly. An external input is reconstructed through the same
	 * interpolation filter as the CV controls before entering the nonlinear
	 * path. Inputs and outputs use a normalized peak level of one.
	 */
	Output StepWithInput(double frequencyHz, double morph, double fold,
		double symmetry, double externalInput, bool useExternalInput,
		bool renderOscillator = true, bool renderFolded = true)
	{
		if (!std::isfinite(frequencyHz) || !std::isfinite(morph) ||
			!std::isfinite(fold) || !std::isfinite(symmetry) ||
			!std::isfinite(externalInput))
		{
			Reset();
			return {};
		}

		if (!_controlsInitialized)
		{
			_frequencyInterpolator->PrimeUpsample(frequencyHz);
			_morphInterpolator->PrimeUpsample(morph);
			_foldInterpolator->PrimeUpsample(fold);
			_symmetryInterpolator->PrimeUpsample(symmetry);
			_controlsInitialized = true;
		}
		if (useExternalInput && !_externalInputInitialized)
		{
			_externalInputInterpolator->PrimeUpsample(externalInput);
			_externalInputInitialized = true;
		}
		else if (!useExternalInput)
			_externalInputInitialized = false;

		const auto frequencies = _frequencyInterpolator->Upsample(frequencyHz);
		const auto morphs = _morphInterpolator->Upsample(morph);
		const auto folds = _foldInterpolator->Upsample(fold);
		const auto symmetries = _symmetryInterpolator->Upsample(symmetry);
		Eigen::Array<double, OversamplingFactor, 1> externalInputs;
		if (useExternalInput)
			externalInputs = _externalInputInterpolator->Upsample(externalInput);
		Eigen::Array<double, OversamplingFactor, 1> oscillatorOutput;
		Eigen::Array<double, OversamplingFactor, 1> foldedOutput;
		const double internalRate = _sampleRate * OversamplingFactor;
		for (int index = 0; index < OversamplingFactor; ++index)
		{
			const double increment = std::clamp(frequencies(index) / internalRate,
				-0.45, 0.45);
			const double triangle = _triangle.Step(increment);
			const double sine = -std::cos(
				6.283185307179586476925286766559 * _triangle.OutputPhase());
			const double shape = std::clamp(morphs(index), 0.0, 1.0);
			const double source = sine + shape * (triangle - sine);
			const double foldTaper = useExternalInput ? 1.0 :
				FoldScaleForFrequency(frequencies(index));
			const double foldAmount =
				std::clamp(folds(index), 0.0, 1.0) * foldTaper;
			// Keep the zero-fold endpoint almost linear, then traverse all four
			// folding stages over the rest of the control. Complementary makeup
			// keeps the unfolded and fully folded endpoints near the same peak level.
			const double drive = 0.5 + 8.5 * foldAmount;
			const double makeup = 2.0 / (1.0 + foldAmount);
			const double folderSource = useExternalInput ?
				externalInputs(index) : source;
			const double alignedDrySource = _useAdaa &&
				_folderSourceInitialized ?
				0.5 * (folderSource + _previousFolderSource) : folderSource;
			_previousFolderSource = folderSource;
			_folderSourceInitialized = true;
			const double foldedInput = drive * folderSource +
				std::clamp(symmetries(index), -1.0, 1.0);
			if (renderOscillator)
				oscillatorOutput(index) = source;
			const double wet = makeup * (_useAdaa ?
				_folder.Process(foldedInput, _character) :
				Wavefolder::Transfer(foldedInput, _character));
			// Fold at zero is an exact bypass. The wet contribution grows with the
			// same control that drives the cascade, avoiding residual character at
			// the minimum setting. ADAA's dry term uses its matching half-sample
			// average so intermediate blends remain phase coherent.
			if (renderFolded)
				foldedOutput(index) = foldAmount <= 0.0 ? folderSource :
					alignedDrySource + foldAmount * (wet - alignedDrySource);
		}
		if (renderOscillator && !_renderedOscillatorLastSample)
			_oscillatorDecimator->Reset();
		if (renderFolded && !_renderedFoldedLastSample)
			_foldedDecimator->Reset();
		Output result{
			renderOscillator ?
				_oscillatorDecimator->Downsample(oscillatorOutput) : 0.0,
			renderFolded ? _foldedDecimator->Downsample(foldedOutput) : 0.0,
		};
		_renderedOscillatorLastSample = renderOscillator;
		_renderedFoldedLastSample = renderFolded;
		if (std::isfinite(result.oscillator) && std::isfinite(result.folded))
			return result;
		Reset();
		return {};
	}

private:
	std::unique_ptr<ResamplerType> _frequencyInterpolator;
	std::unique_ptr<ResamplerType> _morphInterpolator;
	std::unique_ptr<ResamplerType> _foldInterpolator;
	std::unique_ptr<ResamplerType> _symmetryInterpolator;
	std::unique_ptr<ResamplerType> _externalInputInterpolator;
	std::unique_ptr<ResamplerType> _oscillatorDecimator;
	std::unique_ptr<ResamplerType> _foldedDecimator;
	BandlimitedTriangleOscillator _triangle;
	Wavefolder _folder;
	double _sampleRate{48000.0};
	bool _controlsInitialized{};
	bool _externalInputInitialized{};
	bool _renderedOscillatorLastSample{true};
	bool _renderedFoldedLastSample{true};
	bool _useAdaa{};
	double _previousFolderSource{};
	bool _folderSourceInitialized{};
	WavefolderCharacter _character{WavefolderCharacter::Hinge};
};

} // namespace tfdsp
