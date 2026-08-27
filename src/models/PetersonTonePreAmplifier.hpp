#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "tfdsp/approx.hpp"

namespace tfdsp
{

// Nonlinear reduction of the active bass/treble feedback amplifier in the
// Peterson 80 W suitcase preamp (service-manual Figure 11-8).
//
// This is not a post-EQ waveshaper.  The two selected 2N3392 devices drawn as
// a Darlington common-emitter stage and the following selected-2N3392 emitter
// follower are stamped into the same Newton solve as the actual 100 k pots,
// 390 k / 1 M feedback path and frequency-dependent feedback capacitors.  Tone
// position consequently changes both the linear response and the way the
// transistor stage approaches its limited 25 V headroom.
//
// Figure 11-8 places a low-impedance emitter follower immediately before this
// block.  Treating that boundary as a Thevenin voltage drive avoids needlessly
// coupling the preceding two-transistor solve to thirteen passive tone nodes;
// the two 68 k / 6.8 k input arms and every component inside the feedback loop
// remain explicit here. Above the centred Treble position, the 68 k arm has a
// documented extended-range parallel equivalent; centre and cut retain the
// production value. The final vibrato-feed buffer is outside this feedback loop
// and is not used to manufacture overload.
class PetersonTonePreAmplifier
{
public:
	struct Result
	{
		double voltage{};
		int iterations{};
		bool converged{};
	};

	void SetSampleRate(double sampleRate)
	{
		sampleRate_ = std::clamp(sampleRate, 32000.0, 3072000.0);
		outputDcCoefficient_ = 1.0 - std::exp(-TwoPi * OutputCouplingCorner /
			sampleRate_);
		Reset();
	}

	void Reset()
	{
		// The values are only Newton seeds. SolveOperatingPoint() establishes the
		// actual Figure 11-8 bias for centred controls.
		unknowns_ = {0.25, 0.55, 0.85, 0.25, 0.55, 0.30, 0.85,
			1.25, 0.62, 4.45, 3.15, 2.30, 1.62};
		SolveOperatingPoint(0.5, 0.5);
		ResetCapacitors();
		outputDc_ = unknowns_[EmitterQ5];
		controlsInitialized_ = false;
		solverFailures_ = 0;
		maximumIterationsUsed_ = 0;
		totalIterations_ = 0;
		processedSamples_ = 0;
	}

	Result Step(double inputVoltage, double bass, double treble)
	{
		inputVoltage = std::clamp(std::isfinite(inputVoltage) ? inputVoltage :
			0.0, -10.0, 10.0);
		bass = std::clamp(std::isfinite(bass) ? bass : 0.5, 0.0, 1.0);
		treble = std::clamp(std::isfinite(treble) ? treble : 0.5, 0.0, 1.0);
		if (!controlsInitialized_)
		{
			SolveOperatingPoint(bass, treble);
			ResetCapacitors();
			outputDc_ = unknowns_[EmitterQ5];
			controlsInitialized_ = true;
		}
		const double timeStep = 1.0 / sampleRate_;
		std::array<double, CapacitorCount> histories{};
		for (std::size_t index = 0; index < CapacitorCount; ++index)
			histories[index] = capacitors_[index].History(Capacitances[index],
				timeStep);

		bool converged = false;
		int iterationsPerformed = 0;
		for (int iteration = 0; iteration < MaximumNewtonIterations; ++iteration)
		{
			++iterationsPerformed;
			const auto residual = Evaluate(unknowns_, inputVoltage, bass, treble,
				histories, timeStep, false);
			double maximumResidual = 0.0;
			for (const auto& value : residual)
				maximumResidual = std::max(maximumResidual, std::abs(value.value));
			if (maximumResidual < (iteration == 0 ? NewtonTolerance :
				RealtimeRefinementTolerance))
			{
				converged = true;
				break;
			}

			double matrix[UnknownCount][UnknownCount]{};
			std::array<double, UnknownCount> correction{};
			for (std::size_t row = 0; row < UnknownCount; ++row)
			{
				correction[row] = -residual[row].value;
				for (std::size_t column = 0; column < UnknownCount; ++column)
					matrix[row][column] = residual[row].derivative[column];
			}
			if (!SolveLinear(matrix, correction))
				break;

			double damping = 1.0;
			for (double value : correction)
			{
				if (std::abs(value) > MaximumNewtonStep)
					damping = std::min(damping,
						MaximumNewtonStep / std::abs(value));
			}
			const std::array<double, 6> junctionCorrections{
				correction[FeedbackBase] - correction[BaseQ3],
				correction[FeedbackBase] - correction[CollectorOutput],
				correction[BaseQ3],
				correction[BaseQ3] - correction[CollectorOutput],
				correction[BaseQ5] - correction[EmitterQ5],
				correction[BaseQ5]};
			for (double value : junctionCorrections)
			{
				if (std::abs(value) > MaximumJunctionStep)
					damping = std::min(damping,
						MaximumJunctionStep / std::abs(value));
			}

			if (iteration >= 2)
			{
				for (int lineSearch = 0; lineSearch < 5; ++lineSearch)
				{
					auto candidate = unknowns_;
					for (std::size_t index = 0; index < UnknownCount; ++index)
						candidate[index] += damping * correction[index];
					const auto candidateResidual = Evaluate(candidate, inputVoltage,
						bass, treble, histories, timeStep, false);
					double candidateMaximum = 0.0;
					for (const auto& value : candidateResidual)
						candidateMaximum = std::max(candidateMaximum,
							std::abs(value.value));
					if (candidateMaximum < maximumResidual || lineSearch == 4)
					{
						unknowns_ = candidate;
						break;
					}
					damping *= 0.5;
				}
			}
			else
			{
				for (std::size_t index = 0; index < UnknownCount; ++index)
					unknowns_[index] += damping * correction[index];
			}
			// At audio rate the previous solution is the predictor and the linear
			// passive network is exact. An undamped first correction means every BJT
			// junction stayed inside the locally validated Ebers-Moll interval; the
			// checked ngspice level/THD regressions bound the resulting one-step error.
			// Junction-limited or line-searched samples continue to a residual check.
			if (iteration == 0 && damping >= 1.0)
			{
				converged = true;
				break;
			}
		}

		if (!converged)
		{
			const auto residual = Evaluate(unknowns_, inputVoltage, bass, treble,
				histories, timeStep, false);
			double maximumResidual = 0.0;
			for (const auto& value : residual)
				maximumResidual = std::max(maximumResidual, std::abs(value.value));
			converged = maximumResidual < RealtimeFailureTolerance;
		}
		if (!converged)
			++solverFailures_;
		maximumIterationsUsed_ = std::max(maximumIterationsUsed_,
			iterationsPerformed);
		totalIterations_ += static_cast<std::uint64_t>(iterationsPerformed);
		++processedSamples_;

		CommitCapacitors();
		outputDc_ += outputDcCoefficient_ * (unknowns_[EmitterQ5] - outputDc_);
		return {unknowns_[EmitterQ5] - outputDc_, iterationsPerformed,
			converged};
	}

	std::uint64_t SolverFailures() const { return solverFailures_; }
	int MaximumIterationsUsed() const { return maximumIterationsUsed_; }
	double AverageIterations() const
	{
		return processedSamples_ == 0 ? 0.0 :
			static_cast<double>(totalIterations_) /
			static_cast<double>(processedSamples_);
	}

private:
	template <std::size_t Size>
	struct Dual
	{
		double value{};
		std::array<double, Size> derivative{};

		Dual() = default;
		Dual(double scalar) : value(scalar) {}

		static Dual Variable(double scalar, std::size_t index)
		{
			Dual result(scalar);
			result.derivative[index] = 1.0;
			return result;
		}

		Dual& operator+=(const Dual& other)
		{
			value += other.value;
			for (std::size_t index = 0; index < Size; ++index)
				derivative[index] += other.derivative[index];
			return *this;
		}

		Dual& operator-=(const Dual& other)
		{
			value -= other.value;
			for (std::size_t index = 0; index < Size; ++index)
				derivative[index] -= other.derivative[index];
			return *this;
		}
	};

	template <std::size_t Size>
	friend Dual<Size> operator+(Dual<Size> left, const Dual<Size>& right)
	{
		left += right;
		return left;
	}

	template <std::size_t Size>
	friend Dual<Size> operator-(Dual<Size> left, const Dual<Size>& right)
	{
		left -= right;
		return left;
	}

	template <std::size_t Size>
	friend Dual<Size> operator-(const Dual<Size>& value)
	{
		Dual<Size> result(-value.value);
		for (std::size_t index = 0; index < Size; ++index)
			result.derivative[index] = -value.derivative[index];
		return result;
	}

	template <std::size_t Size>
	friend Dual<Size> operator*(const Dual<Size>& value, double scalar)
	{
		Dual<Size> result(value.value * scalar);
		for (std::size_t index = 0; index < Size; ++index)
			result.derivative[index] = value.derivative[index] * scalar;
		return result;
	}

	template <std::size_t Size>
	friend Dual<Size> operator*(double scalar, const Dual<Size>& value)
	{
		return value * scalar;
	}

	template <std::size_t Size>
	friend Dual<Size> operator/(const Dual<Size>& value, double scalar)
	{
		return value * (1.0 / scalar);
	}

	struct CapacitorCompanion
	{
		double previousVoltage{};

		void Reset(double voltage) { previousVoltage = voltage; }
		double History(double capacitance, double timeStep) const
		{
			return -capacitance * previousVoltage / timeStep;
		}
		void Commit(double voltage) { previousVoltage = voltage; }
	};

	struct DeviceCurrents
	{
		Dual<13> collector;
		Dual<13> base;
		Dual<13> emitter;
	};

	enum Unknown : std::size_t
	{
		TrebleTop,
		TrebleWiper,
		TrebleBottom,
		TrebleJunction,
		BassWiper,
		BassTop,
		BassBottom,
		FeedbackBase,
		BaseQ3,
		CollectorOutput,
		CompensationJunction,
		BaseQ5,
		EmitterQ5,
		UnknownCount
	};

	enum Capacitor : std::size_t
	{
		TrebleSeriesCapacitor,
		TrebleFeedbackCapacitor,
		BassTopCapacitor,
		BassBottomCapacitor,
		FollowerCompensationCapacitor,
		FollowerBaseCapacitor,
		CapacitorCount
	};

	using CircuitDual = Dual<UnknownCount>;
	using Residual = std::array<CircuitDual, UnknownCount>;

	static CircuitDual LimitedExp(const CircuitDual& argument)
	{
		if (argument.value <= -20.0)
			return CircuitDual(0.0);
		const double limited = std::min(argument.value, 40.0);
		const double exponential = static_cast<double>(Exp2Taylor5(
			static_cast<float>(limited * Log2E)));
		CircuitDual result(exponential);
		if (argument.value < 40.0)
		{
			for (std::size_t index = 0; index < UnknownCount; ++index)
				result.derivative[index] = exponential *
					argument.derivative[index];
		}
		return result;
	}

	static DeviceCurrents Npn(const CircuitDual& collector,
		const CircuitDual& base, const CircuitDual& emitter)
	{
		const auto forwardExp = LimitedExp((base - emitter) / ThermalVoltage);
		const auto reverseExp = LimitedExp((base - collector) / ThermalVoltage);
		const auto forward = (forwardExp - CircuitDual(1.0)) *
			(SaturationCurrent / ForwardAlpha);
		const auto reverse = (reverseExp - CircuitDual(1.0)) *
			(SaturationCurrent / ReverseAlpha);
		DeviceCurrents result;
		result.collector = ForwardAlpha * forward - reverse;
		result.emitter = -forward + ReverseAlpha * reverse;
		result.base = -(result.collector + result.emitter);
		return result;
	}

	static void Stamp(Residual& residual, int positiveNode, int negativeNode,
		const CircuitDual& current)
	{
		if (positiveNode >= 0)
			residual[static_cast<std::size_t>(positiveNode)] += current;
		if (negativeNode >= 0)
			residual[static_cast<std::size_t>(negativeNode)] -= current;
	}

	static void StampResistor(Residual& residual,
		const std::array<CircuitDual, UnknownCount>& node, int positiveNode,
		int negativeNode, double resistance, double positiveFixed = 0.0,
		double negativeFixed = 0.0)
	{
		const CircuitDual positive = positiveNode >= 0 ?
			node[static_cast<std::size_t>(positiveNode)] :
			CircuitDual(positiveFixed);
		const CircuitDual negative = negativeNode >= 0 ?
			node[static_cast<std::size_t>(negativeNode)] :
			CircuitDual(negativeFixed);
		Stamp(residual, positiveNode, negativeNode,
			(positive - negative) / resistance);
	}

	static void StampCapacitor(Residual& residual,
		const std::array<CircuitDual, UnknownCount>& node, int positiveNode,
		int negativeNode, double capacitance, double history, double timeStep)
	{
		const CircuitDual positive = positiveNode >= 0 ?
			node[static_cast<std::size_t>(positiveNode)] : CircuitDual(0.0);
		const CircuitDual negative = negativeNode >= 0 ?
			node[static_cast<std::size_t>(negativeNode)] : CircuitDual(0.0);
		Stamp(residual, positiveNode, negativeNode,
			(positive - negative) * (capacitance / timeStep) +
				CircuitDual(history));
	}

	Residual Evaluate(const std::array<double, UnknownCount>& values,
		double inputVoltage, double bass, double treble,
		const std::array<double, CapacitorCount>& histories, double timeStep,
		bool dcMode) const
	{
		std::array<CircuitDual, UnknownCount> node{};
		for (std::size_t index = 0; index < UnknownCount; ++index)
			node[index] = CircuitDual::Variable(values[index], index);
		Residual residual{};

		const double trebleTopResistance = PotMinimumResistance +
			(1.0 - treble) * TonePotResistance;
		const double trebleBottomResistance = PotMinimumResistance +
			treble * TonePotResistance;
		const double bassTopResistance = PotMinimumResistance +
			(1.0 - bass) * TonePotResistance;
		const double bassBottomResistance = PotMinimumResistance +
			bass * TonePotResistance;

		// Figure 11-8's input arms and the two actual three-terminal controls.
		// The original 68 k treble input arm yields only about 4 dB of boost from
		// the centred setting. Above centre the module progressively parallels
		// that arm down to 32 k, extending boost while retaining the real active-
		// feedback topology, original centre calibration and original cut side.
		// This is an explicit sound-design range extension, not a claim about the
		// production Peterson component value.
		const double trebleBoostExtension = std::max(0.0,
			2.0 * (treble - 0.5));
		const double trebleInputResistance = TrebleInputResistance * std::pow(
			TrebleBoostInputResistanceRatio, trebleBoostExtension);
		StampResistor(residual, node, TrebleTop, -1, trebleInputResistance,
			0.0, inputVoltage);
		StampResistor(residual, node, BassTop, -1, BassInputResistance,
			0.0, inputVoltage);
		StampResistor(residual, node, TrebleTop, TrebleWiper,
			trebleTopResistance);
		StampResistor(residual, node, TrebleWiper, TrebleBottom,
			trebleBottomResistance);
		StampResistor(residual, node, TrebleBottom, CollectorOutput,
			ToneEndResistance);
		StampResistor(residual, node, BassTop, BassWiper, bassTopResistance);
		StampResistor(residual, node, BassWiper, BassBottom,
			bassBottomResistance);
		StampResistor(residual, node, BassBottom, CollectorOutput,
			ToneEndResistance);
		StampResistor(residual, node, TrebleJunction, BassWiper,
			ToneBridgeResistance);

		// The bias/feedback divider drives Q4, whose emitter drives Q3. Their
		// collectors share the 12 k output node exactly as drawn.
		StampResistor(residual, node, FeedbackBase, -1,
			FeedbackGroundResistance);
		StampResistor(residual, node, FeedbackBase, CollectorOutput,
			FeedbackResistance);
		StampResistor(residual, node, CollectorOutput, -1,
			CollectorResistance, 0.0, SupplyVoltage);

		// Q5 follows the Darlington collector through the compensated 6.8 k /
		// 68 k network. The 6.8 k emitter resistor is the dominant load seen by
		// the historical volume/output network.
		StampResistor(residual, node, CollectorOutput, CompensationJunction,
			FollowerSeriesResistance);
		StampResistor(residual, node, CompensationJunction, BaseQ5,
			FollowerBaseResistance);
		StampResistor(residual, node, EmitterQ5, -1,
			FollowerEmitterResistance);
		StampResistor(residual, node, EmitterQ5, -1,
			HistoricalVolumeResistance);

		if (!dcMode)
		{
			StampCapacitor(residual, node, TrebleWiper, TrebleJunction,
				Capacitances[TrebleSeriesCapacitor],
				histories[TrebleSeriesCapacitor], timeStep);
			StampCapacitor(residual, node, TrebleJunction, FeedbackBase,
				Capacitances[TrebleFeedbackCapacitor],
				histories[TrebleFeedbackCapacitor], timeStep);
			StampCapacitor(residual, node, BassWiper, BassTop,
				Capacitances[BassTopCapacitor],
				histories[BassTopCapacitor], timeStep);
			StampCapacitor(residual, node, BassWiper, BassBottom,
				Capacitances[BassBottomCapacitor],
				histories[BassBottomCapacitor], timeStep);
			StampCapacitor(residual, node, CompensationJunction, EmitterQ5,
				Capacitances[FollowerCompensationCapacitor],
				histories[FollowerCompensationCapacitor], timeStep);
			StampCapacitor(residual, node, BaseQ5, -1,
				Capacitances[FollowerBaseCapacitor],
				histories[FollowerBaseCapacitor], timeStep);
		}

		const auto q4 = Npn(node[CollectorOutput], node[FeedbackBase],
			node[BaseQ3]);
		const auto q3 = Npn(node[CollectorOutput], node[BaseQ3],
			CircuitDual(0.0));
		const auto q5 = Npn(CircuitDual(SupplyVoltage), node[BaseQ5],
			node[EmitterQ5]);
		residual[CollectorOutput] += q4.collector;
		residual[FeedbackBase] += q4.base;
		residual[BaseQ3] += q4.emitter;
		residual[CollectorOutput] += q3.collector;
		residual[BaseQ3] += q3.base;
		residual[BaseQ5] += q5.base;
		residual[EmitterQ5] += q5.emitter;
		return residual;
	}

	static bool SolveLinear(double matrix[UnknownCount][UnknownCount],
		std::array<double, UnknownCount>& rightHandSide)
	{
		for (std::size_t column = 0; column < UnknownCount; ++column)
		{
			std::size_t pivot = column;
			for (std::size_t row = column + 1; row < UnknownCount; ++row)
			{
				if (std::abs(matrix[row][column]) >
					std::abs(matrix[pivot][column]))
					pivot = row;
			}
			if (std::abs(matrix[pivot][column]) < 1.0e-18)
				return false;
			if (pivot != column)
			{
				for (std::size_t index = column; index < UnknownCount; ++index)
					std::swap(matrix[column][index], matrix[pivot][index]);
				std::swap(rightHandSide[column], rightHandSide[pivot]);
			}
			for (std::size_t row = column + 1; row < UnknownCount; ++row)
			{
				const double factor = matrix[row][column] / matrix[column][column];
				for (std::size_t index = column + 1; index < UnknownCount; ++index)
					matrix[row][index] -= factor * matrix[column][index];
				rightHandSide[row] -= factor * rightHandSide[column];
			}
		}
		for (std::size_t reverse = UnknownCount; reverse-- > 0;)
		{
			double value = rightHandSide[reverse];
			for (std::size_t column = reverse + 1; column < UnknownCount; ++column)
				value -= matrix[reverse][column] * rightHandSide[column];
			rightHandSide[reverse] = value / matrix[reverse][reverse];
		}
		return true;
	}

	void SolveOperatingPoint(double bass, double treble)
	{
		const std::array<double, CapacitorCount> histories{};
		for (int iteration = 0; iteration < 100; ++iteration)
		{
			const auto residual = Evaluate(unknowns_, 0.0, bass, treble,
				histories, 1.0, true);
			double maximumResidual = 0.0;
			for (const auto& value : residual)
				maximumResidual = std::max(maximumResidual, std::abs(value.value));
			if (maximumResidual < NewtonTolerance)
				break;
			double matrix[UnknownCount][UnknownCount]{};
			std::array<double, UnknownCount> correction{};
			for (std::size_t row = 0; row < UnknownCount; ++row)
			{
				correction[row] = -residual[row].value;
				for (std::size_t column = 0; column < UnknownCount; ++column)
					matrix[row][column] = residual[row].derivative[column];
			}
			if (!SolveLinear(matrix, correction))
				break;
			double damping = 1.0;
			for (double value : correction)
			{
				if (std::abs(value) > MaximumJunctionStep)
					damping = std::min(damping,
						MaximumJunctionStep / std::abs(value));
			}
			for (std::size_t index = 0; index < UnknownCount; ++index)
				unknowns_[index] += damping * correction[index];
		}
	}

	void ResetCapacitors()
	{
		capacitors_[TrebleSeriesCapacitor].Reset(
			unknowns_[TrebleWiper] - unknowns_[TrebleJunction]);
		capacitors_[TrebleFeedbackCapacitor].Reset(
			unknowns_[TrebleJunction] - unknowns_[FeedbackBase]);
		capacitors_[BassTopCapacitor].Reset(
			unknowns_[BassWiper] - unknowns_[BassTop]);
		capacitors_[BassBottomCapacitor].Reset(
			unknowns_[BassWiper] - unknowns_[BassBottom]);
		capacitors_[FollowerCompensationCapacitor].Reset(
			unknowns_[CompensationJunction] - unknowns_[EmitterQ5]);
		capacitors_[FollowerBaseCapacitor].Reset(unknowns_[BaseQ5]);
	}

	void CommitCapacitors()
	{
		capacitors_[TrebleSeriesCapacitor].Commit(
			unknowns_[TrebleWiper] - unknowns_[TrebleJunction]);
		capacitors_[TrebleFeedbackCapacitor].Commit(
			unknowns_[TrebleJunction] - unknowns_[FeedbackBase]);
		capacitors_[BassTopCapacitor].Commit(
			unknowns_[BassWiper] - unknowns_[BassTop]);
		capacitors_[BassBottomCapacitor].Commit(
			unknowns_[BassWiper] - unknowns_[BassBottom]);
		capacitors_[FollowerCompensationCapacitor].Commit(
			unknowns_[CompensationJunction] - unknowns_[EmitterQ5]);
		capacitors_[FollowerBaseCapacitor].Commit(unknowns_[BaseQ5]);
	}

	static constexpr double SupplyVoltage = 25.0;
	static constexpr double TonePotResistance = 100000.0;
	static constexpr double PotMinimumResistance = 1.0;
	static constexpr double TrebleInputResistance = 68000.0;
	static constexpr double TrebleBoostInputResistanceRatio = 0.47;
	static constexpr double BassInputResistance = 6800.0;
	static constexpr double ToneEndResistance = 6800.0;
	static constexpr double ToneBridgeResistance = 22000.0;
	static constexpr double FeedbackGroundResistance = 390000.0;
	static constexpr double FeedbackResistance = 1000000.0;
	static constexpr double CollectorResistance = 12000.0;
	static constexpr double FollowerSeriesResistance = 6800.0;
	static constexpr double FollowerBaseResistance = 68000.0;
	static constexpr double FollowerEmitterResistance = 6800.0;
	static constexpr double HistoricalVolumeResistance = 100000.0;
	static constexpr std::array<double, CapacitorCount> Capacitances{
		4.7e-9, 47.0e-9, 0.1e-6, 0.1e-6, 4.7e-9, 470.0e-12};
	static constexpr double SaturationCurrent = 2.0e-14;
	static constexpr double ThermalVoltage = 0.02585;
	static constexpr double ForwardAlpha = 220.0 / 221.0;
	static constexpr double ReverseAlpha = 3.0 / 4.0;
	static constexpr double Log2E = 1.4426950408889634074;
	static constexpr double TwoPi = 6.2831853071795864769;
	static constexpr double OutputCouplingCorner = 2.0;
	static constexpr int MaximumNewtonIterations = 8;
	static constexpr double NewtonTolerance = 1.0e-8;
	static constexpr double RealtimeRefinementTolerance = 2.0e-6;
	static constexpr double RealtimeFailureTolerance = 2.0e-4;
	static constexpr double MaximumNewtonStep = 8.0;
	static constexpr double MaximumJunctionStep = 0.30;

	double sampleRate_ = 96000.0;
	std::array<double, UnknownCount> unknowns_{};
	std::array<CapacitorCompanion, CapacitorCount> capacitors_{};
	double outputDc_{};
	double outputDcCoefficient_ = 0.00013;
	std::uint64_t solverFailures_{};
	int maximumIterationsUsed_{};
	std::uint64_t totalIterations_{};
	std::uint64_t processedSamples_{};
	bool controlsInitialized_{};
};

} // namespace tfdsp
