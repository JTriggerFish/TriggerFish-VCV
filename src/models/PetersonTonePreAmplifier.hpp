#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

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
		UpdateRealtimeControls(bass, treble);
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
			const auto residual = EvaluateRealtime(unknowns_, inputVoltage,
				histories, timeStep);
			double maximumResidual = 0.0;
			for (const auto& value : residual)
				maximumResidual = std::max(maximumResidual, std::abs(value));
			if (maximumResidual < (iteration == 0 ? NewtonTolerance :
				RealtimeRefinementTolerance))
			{
				converged = true;
				break;
			}

			double matrix[UnknownCount][UnknownCount]{};
			std::array<double, UnknownCount> correction{};
			BuildRealtimeJacobian(matrix, timeStep);
			for (std::size_t row = 0; row < UnknownCount; ++row)
				correction[row] = -residual[row];
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
					const auto candidateResidual = EvaluateRealtime(candidate,
						inputVoltage, histories, timeStep);
					double candidateMaximum = 0.0;
					for (const auto& value : candidateResidual)
						candidateMaximum = std::max(candidateMaximum,
							std::abs(value));
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
			const auto residual = EvaluateRealtime(unknowns_, inputVoltage,
				histories, timeStep);
			double maximumResidual = 0.0;
			for (const auto& value : residual)
				maximumResidual = std::max(maximumResidual, std::abs(value));
			converged = maximumResidual < RealtimeFailureTolerance;
		}
		if (!converged)
			++solverFailures_;
		maximumIterationsUsed_ = std::max(maximumIterationsUsed_,
			iterationsPerformed);
		totalIterations_ += static_cast<std::uint64_t>(iterationsPerformed);
		++processedSamples_;

		CommitCapacitors(timeStep);
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
		double previousCurrent{};

		void Reset(double voltage)
		{
			previousVoltage = voltage;
			previousCurrent = 0.0;
		}
		double History(double capacitance, double timeStep) const
		{
			return -previousCurrent -
				2.0 * capacitance * previousVoltage / timeStep;
		}
		void Commit(double voltage, double capacitance, double timeStep)
		{
			previousCurrent = 2.0 * capacitance * voltage / timeStep +
				History(capacitance, timeStep);
			previousVoltage = voltage;
		}
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

	struct DeviceLinearization
	{
		std::array<double, 3> current{}; // collector, base, emitter
		double derivative[3][3]{}; // terminal current by c/b/e voltage
	};

	struct ExponentialLinearization
	{
		double value{};
		double slope{};
	};

	static ExponentialLinearization LimitedExpLinearized(double argument)
	{
		if (argument <= -20.0)
			return {};
		const double limited = std::min(argument, 40.0);
		const double value = static_cast<double>(Exp2Taylor5(
			static_cast<float>(limited * Log2E)));
		return {value, argument < 40.0 ? value : 0.0};
	}

	static DeviceLinearization NpnLinearization(double collector, double base,
		double emitter)
	{
		const auto forwardExp = LimitedExpLinearized(
			(base - emitter) / ThermalVoltage);
		const auto reverseExp = LimitedExpLinearized(
			(base - collector) / ThermalVoltage);
		const double forward = SaturationCurrent / ForwardAlpha *
			(forwardExp.value - 1.0);
		const double reverse = SaturationCurrent / ReverseAlpha *
			(reverseExp.value - 1.0);
		const double gf = SaturationCurrent /
			(ForwardAlpha * ThermalVoltage) * forwardExp.slope;
		const double gr = SaturationCurrent /
			(ReverseAlpha * ThermalVoltage) * reverseExp.slope;
		DeviceLinearization result;
		result.current[0] = ForwardAlpha * forward - reverse;
		result.current[2] = -forward + ReverseAlpha * reverse;
		result.current[1] = -(result.current[0] + result.current[2]);
		result.derivative[0][0] = gr;
		result.derivative[0][1] = ForwardAlpha * gf - gr;
		result.derivative[0][2] = -ForwardAlpha * gf;
		result.derivative[2][0] = -ReverseAlpha * gr;
		result.derivative[2][1] = -gf + ReverseAlpha * gr;
		result.derivative[2][2] = gf;
		for (int column = 0; column < 3; ++column)
			result.derivative[1][column] = -result.derivative[0][column] -
				result.derivative[2][column];
		return result;
	}

	using RealtimeResidual = std::array<double, UnknownCount>;

	void UpdateRealtimeControls(double bass, double treble)
	{
		if (bass != cachedRealtimeBass_)
		{
			cachedRealtimeBass_ = bass;
			realtimeBassTopResistance_ = PotMinimumResistance +
				(1.0 - bass) * TonePotResistance;
			realtimeBassBottomResistance_ = PotMinimumResistance +
				bass * TonePotResistance;
		}
		if (treble != cachedRealtimeTreble_)
		{
			cachedRealtimeTreble_ = treble;
			realtimeTrebleTopResistance_ = PotMinimumResistance +
				(1.0 - treble) * TonePotResistance;
			realtimeTrebleBottomResistance_ = PotMinimumResistance +
				treble * TonePotResistance;
			const double extension = std::max(0.0, 2.0 * (treble - 0.5));
			realtimeTrebleInputResistance_ = TrebleInputResistance * std::pow(
				TrebleBoostInputResistanceRatio, extension);
		}
	}

	RealtimeResidual EvaluateRealtime(
		const std::array<double, UnknownCount>& node, double inputVoltage,
		const std::array<double, CapacitorCount>& histories,
		double timeStep) const
	{
		RealtimeResidual residual{};
		const auto stamp = [&](int positive, int negative, double current)
		{
			if (positive >= 0)
				residual[static_cast<std::size_t>(positive)] += current;
			if (negative >= 0)
				residual[static_cast<std::size_t>(negative)] -= current;
		};
		const auto voltage = [&](int index, double fixed)
		{
			return index >= 0 ? node[static_cast<std::size_t>(index)] : fixed;
		};
		const auto resistor = [&](int positive, int negative, double resistance,
			double positiveFixed = 0.0, double negativeFixed = 0.0)
		{
			stamp(positive, negative,
				(voltage(positive, positiveFixed) -
					voltage(negative, negativeFixed)) / resistance);
		};
		const auto capacitor = [&](int positive, int negative,
			double capacitance, double history)
		{
			stamp(positive, negative,
				(voltage(positive, 0.0) - voltage(negative, 0.0)) *
					(2.0 * capacitance / timeStep) + history);
		};

		resistor(TrebleTop, -1, realtimeTrebleInputResistance_, 0.0,
			inputVoltage);
		resistor(BassTop, -1, BassInputResistance, 0.0, inputVoltage);
		resistor(TrebleTop, TrebleWiper, realtimeTrebleTopResistance_);
		resistor(TrebleWiper, TrebleBottom,
			realtimeTrebleBottomResistance_);
		resistor(TrebleBottom, CollectorOutput, ToneEndResistance);
		resistor(BassTop, BassWiper, realtimeBassTopResistance_);
		resistor(BassWiper, BassBottom, realtimeBassBottomResistance_);
		resistor(BassBottom, CollectorOutput, ToneEndResistance);
		resistor(TrebleJunction, BassWiper, ToneBridgeResistance);
		resistor(FeedbackBase, -1, FeedbackGroundResistance);
		resistor(FeedbackBase, CollectorOutput, FeedbackResistance);
		resistor(CollectorOutput, -1, CollectorResistance, 0.0,
			SupplyVoltage);
		resistor(CollectorOutput, CompensationJunction,
			FollowerSeriesResistance);
		resistor(CompensationJunction, BaseQ5, FollowerBaseResistance);
		resistor(EmitterQ5, -1, FollowerEmitterResistance);
		resistor(EmitterQ5, -1, HistoricalVolumeResistance);

		capacitor(TrebleWiper, TrebleJunction,
			Capacitances[TrebleSeriesCapacitor],
			histories[TrebleSeriesCapacitor]);
		capacitor(TrebleJunction, FeedbackBase,
			Capacitances[TrebleFeedbackCapacitor],
			histories[TrebleFeedbackCapacitor]);
		capacitor(BassWiper, BassTop, Capacitances[BassTopCapacitor],
			histories[BassTopCapacitor]);
		capacitor(BassWiper, BassBottom, Capacitances[BassBottomCapacitor],
			histories[BassBottomCapacitor]);
		capacitor(CompensationJunction, EmitterQ5,
			Capacitances[FollowerCompensationCapacitor],
			histories[FollowerCompensationCapacitor]);
		capacitor(BaseQ5, -1, Capacitances[FollowerBaseCapacitor],
			histories[FollowerBaseCapacitor]);

		const auto q4 = NpnLinearization(node[CollectorOutput],
			node[FeedbackBase], node[BaseQ3]);
		const auto q3 = NpnLinearization(node[CollectorOutput], node[BaseQ3],
			0.0);
		const auto q5 = NpnLinearization(SupplyVoltage, node[BaseQ5],
			node[EmitterQ5]);
		residual[CollectorOutput] += q4.current[0] + q3.current[0];
		residual[FeedbackBase] += q4.current[1];
		residual[BaseQ3] += q4.current[2] + q3.current[1];
		residual[BaseQ5] += q5.current[1];
		residual[EmitterQ5] += q5.current[2];
		return residual;
	}

	void BuildRealtimeJacobian(double matrix[UnknownCount][UnknownCount],
		double timeStep) const
	{
		const auto stampConductance = [&](int positive, int negative,
			double conductance)
		{
			if (positive >= 0)
			{
				matrix[positive][positive] += conductance;
				if (negative >= 0)
					matrix[positive][negative] -= conductance;
			}
			if (negative >= 0)
			{
				matrix[negative][negative] += conductance;
				if (positive >= 0)
					matrix[negative][positive] -= conductance;
			}
		};
		const auto resistor = [&](int positive, int negative, double resistance)
		{
			stampConductance(positive, negative, 1.0 / resistance);
		};
		const auto capacitor = [&](int positive, int negative,
			double capacitance)
		{
			stampConductance(positive, negative,
				2.0 * capacitance / timeStep);
		};
		const auto stampDevice = [&](const DeviceLinearization& device,
			const std::array<int, 3>& nodes)
		{
			for (int terminal = 0; terminal < 3; ++terminal)
				if (nodes[terminal] >= 0)
					for (int voltage = 0; voltage < 3; ++voltage)
						if (nodes[voltage] >= 0)
							matrix[nodes[terminal]][nodes[voltage]] +=
								device.derivative[terminal][voltage];
		};

		resistor(TrebleTop, -1, realtimeTrebleInputResistance_);
		resistor(BassTop, -1, BassInputResistance);
		resistor(TrebleTop, TrebleWiper, realtimeTrebleTopResistance_);
		resistor(TrebleWiper, TrebleBottom,
			realtimeTrebleBottomResistance_);
		resistor(TrebleBottom, CollectorOutput, ToneEndResistance);
		resistor(BassTop, BassWiper, realtimeBassTopResistance_);
		resistor(BassWiper, BassBottom, realtimeBassBottomResistance_);
		resistor(BassBottom, CollectorOutput, ToneEndResistance);
		resistor(TrebleJunction, BassWiper, ToneBridgeResistance);
		resistor(FeedbackBase, -1, FeedbackGroundResistance);
		resistor(FeedbackBase, CollectorOutput, FeedbackResistance);
		resistor(CollectorOutput, -1, CollectorResistance);
		resistor(CollectorOutput, CompensationJunction,
			FollowerSeriesResistance);
		resistor(CompensationJunction, BaseQ5, FollowerBaseResistance);
		resistor(EmitterQ5, -1, FollowerEmitterResistance);
		resistor(EmitterQ5, -1, HistoricalVolumeResistance);
		capacitor(TrebleWiper, TrebleJunction,
			Capacitances[TrebleSeriesCapacitor]);
		capacitor(TrebleJunction, FeedbackBase,
			Capacitances[TrebleFeedbackCapacitor]);
		capacitor(BassWiper, BassTop, Capacitances[BassTopCapacitor]);
		capacitor(BassWiper, BassBottom, Capacitances[BassBottomCapacitor]);
		capacitor(CompensationJunction, EmitterQ5,
			Capacitances[FollowerCompensationCapacitor]);
		capacitor(BaseQ5, -1, Capacitances[FollowerBaseCapacitor]);

		stampDevice(NpnLinearization(unknowns_[CollectorOutput],
			unknowns_[FeedbackBase], unknowns_[BaseQ3]),
			{CollectorOutput, FeedbackBase, BaseQ3});
		stampDevice(NpnLinearization(unknowns_[CollectorOutput],
			unknowns_[BaseQ3], 0.0), {CollectorOutput, BaseQ3, -1});
		stampDevice(NpnLinearization(SupplyVoltage, unknowns_[BaseQ5],
			unknowns_[EmitterQ5]), {-1, BaseQ5, EmitterQ5});
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
			(positive - negative) * (2.0 * capacitance / timeStep) +
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

	void CommitCapacitors(double timeStep)
	{
		capacitors_[TrebleSeriesCapacitor].Commit(
			unknowns_[TrebleWiper] - unknowns_[TrebleJunction],
			Capacitances[TrebleSeriesCapacitor], timeStep);
		capacitors_[TrebleFeedbackCapacitor].Commit(
			unknowns_[TrebleJunction] - unknowns_[FeedbackBase],
			Capacitances[TrebleFeedbackCapacitor], timeStep);
		capacitors_[BassTopCapacitor].Commit(
			unknowns_[BassWiper] - unknowns_[BassTop],
			Capacitances[BassTopCapacitor], timeStep);
		capacitors_[BassBottomCapacitor].Commit(
			unknowns_[BassWiper] - unknowns_[BassBottom],
			Capacitances[BassBottomCapacitor], timeStep);
		capacitors_[FollowerCompensationCapacitor].Commit(
			unknowns_[CompensationJunction] - unknowns_[EmitterQ5],
			Capacitances[FollowerCompensationCapacitor], timeStep);
		capacitors_[FollowerBaseCapacitor].Commit(unknowns_[BaseQ5],
			Capacitances[FollowerBaseCapacitor], timeStep);
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
	double cachedRealtimeBass_{std::numeric_limits<double>::quiet_NaN()};
	double cachedRealtimeTreble_{std::numeric_limits<double>::quiet_NaN()};
	double realtimeTrebleTopResistance_{50'001.0};
	double realtimeTrebleBottomResistance_{50'001.0};
	double realtimeBassTopResistance_{50'001.0};
	double realtimeBassBottomResistance_{50'001.0};
	double realtimeTrebleInputResistance_{TrebleInputResistance};
};

} // namespace tfdsp
