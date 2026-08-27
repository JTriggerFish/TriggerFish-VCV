#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "tfdsp/approx.hpp"

namespace tfdsp
{

// Nonlinear audio-rate reduction of the input section in Peterson Figure 11-8.
// The two selected 2N3392 devices are solved as the circuit drawn: Q1 is a
// 390k/33k-biased common-emitter stage with 12k collector and 1.5k emitter
// resistors; Q2 is its direct-coupled emitter follower with a 4.7k emitter
// resistor. The 0.22 uF input, 330 pF base/collector and 5 uF output capacitors
// are backward-Euler companions in the same Newton solve.
//
// Figure 11-8's later active tone-feedback section is solved separately by
// PetersonTonePreAmplifier at the low-impedance Q2 emitter-follower boundary.
// It is not represented by a post-EQ waveshaper. This class remains deliberately
// scoped to Q1/Q2 so each solve stays small enough for the audio thread.
class PetersonPreAmplifier
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
		Reset();
	}

	void Reset()
	{
		unknowns_ = {1.95, 1.25, 14.0, 13.35, 0.0};
		SolveOperatingPoint();
		inputCapacitor_.Reset(unknowns_[BaseQ1]);
		millerCapacitor_.Reset(unknowns_[BaseQ1] - unknowns_[CollectorQ1BaseQ2]);
		outputCapacitor_.Reset(unknowns_[EmitterQ2] - unknowns_[Output]);
		solverFailures_ = 0;
		maximumIterationsUsed_ = 0;
		totalIterations_ = 0;
		processedSamples_ = 0;
	}

	Result Step(double inputVoltage)
	{
		inputVoltage = std::clamp(std::isfinite(inputVoltage) ? inputVoltage :
			0.0, -10.0, 10.0);
		// Ten volts is over 7.5 times the calibrated maximum hard-note input and
		// above every tested musical/polyphonic case. It is only a boundary for
		// impossible Rack/CV excursions, where the Ebers-Moll proxy lacks the real
		// 2N3392's reverse-junction avalanche and must not be extrapolated.
		const double timeStep = 1.0 / sampleRate_;
		const Histories histories{
			inputCapacitor_.History(InputCapacitance, timeStep),
			millerCapacitor_.History(MillerCapacitance, timeStep),
			outputCapacitor_.History(OutputCapacitance, timeStep)};

		bool converged = false;
		int iterationsPerformed = 0;
		for (int iteration = 0; iteration < MaximumNewtonIterations; ++iteration)
		{
			++iterationsPerformed;
			const auto residual = Evaluate(unknowns_, inputVoltage, histories,
				timeStep, false);
			double maximumResidual = 0.0;
			for (double value : residual)
				maximumResidual = std::max(maximumResidual, std::abs(value));
			if (maximumResidual < (iteration == 0 ? NewtonTolerance :
				RealtimeRefinementTolerance))
			{
				converged = true;
				break;
			}

			double matrix[UnknownCount][UnknownCount]{};
			std::array<double, UnknownCount> correction{};
			BuildJacobian(matrix, timeStep, false);
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
			const std::array<double, 4> junctionCorrections{
				correction[BaseQ1] - correction[EmitterQ1],
				correction[BaseQ1] - correction[CollectorQ1BaseQ2],
				correction[CollectorQ1BaseQ2] - correction[EmitterQ2],
				correction[CollectorQ1BaseQ2]};
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
					for (std::size_t i = 0; i < UnknownCount; ++i)
						candidate[i] += damping * correction[i];
					const auto candidateResidual = Evaluate(candidate,
						inputVoltage, histories, timeStep, false);
					double candidateMaximum = 0.0;
					for (double value : candidateResidual)
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
				for (std::size_t i = 0; i < UnknownCount; ++i)
					unknowns_[i] += damping * correction[i];
			}
			if (iteration == 0 && damping >= 1.0)
			{
				converged = true;
				break;
			}
		}

		if (!converged)
		{
			const auto residual = Evaluate(unknowns_, inputVoltage, histories,
				timeStep, false);
			double maximumResidual = 0.0;
			for (double value : residual)
				maximumResidual = std::max(maximumResidual, std::abs(value));
			converged = maximumResidual < RealtimeFailureTolerance;
		}
		if (!converged)
			++solverFailures_;
		maximumIterationsUsed_ = std::max(maximumIterationsUsed_,
			iterationsPerformed);
		totalIterations_ += static_cast<std::uint64_t>(iterationsPerformed);
		++processedSamples_;

		inputCapacitor_.Commit(unknowns_[BaseQ1] - inputVoltage,
			InputCapacitance, timeStep);
		millerCapacitor_.Commit(unknowns_[BaseQ1] -
			unknowns_[CollectorQ1BaseQ2], MillerCapacitance, timeStep);
		outputCapacitor_.Commit(unknowns_[EmitterQ2] - unknowns_[Output],
			OutputCapacitance, timeStep);
		return {unknowns_[Output], iterationsPerformed, converged};
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
	struct CapacitorCompanion
	{
		double previousVoltage{};

		void Reset(double voltage) { previousVoltage = voltage; }
		double History(double capacitance, double timeStep) const
		{
			return -capacitance * previousVoltage / timeStep;
		}
		void Commit(double voltage, double, double)
		{
			previousVoltage = voltage;
		}
	};

	struct DeviceParameters
	{
		double saturationCurrent;
		double thermalVoltage;
		double forwardAlpha;
		double reverseAlpha;
	};

	struct DeviceLinearization
	{
		std::array<double, 3> current{}; // collector, base, emitter
		double derivative[3][3]{}; // terminal current by c/b/e voltage
	};

	struct Histories
	{
		double input;
		double miller;
		double output;
	};

	enum Unknown : std::size_t
	{
		BaseQ1,
		EmitterQ1,
		CollectorQ1BaseQ2,
		EmitterQ2,
		Output,
		UnknownCount
	};

	static double LimitedExp(double argument)
	{
		if (argument <= -20.0)
			return 0.0;
		return static_cast<double>(Exp2Taylor5(static_cast<float>(
			std::min(argument, 40.0) * Log2E)));
	}

	static DeviceLinearization NpnLinearization(double collector, double base,
		double emitter)
	{
		const double forwardExp = LimitedExp((base - emitter) /
			Transistor.thermalVoltage);
		const double reverseExp = LimitedExp((base - collector) /
			Transistor.thermalVoltage);
		const double forward = Transistor.saturationCurrent /
			Transistor.forwardAlpha * (forwardExp - 1.0);
		const double reverse = Transistor.saturationCurrent /
			Transistor.reverseAlpha * (reverseExp - 1.0);
		const double gf = Transistor.saturationCurrent /
			(Transistor.forwardAlpha * Transistor.thermalVoltage) * forwardExp;
		const double gr = Transistor.saturationCurrent /
			(Transistor.reverseAlpha * Transistor.thermalVoltage) * reverseExp;
		DeviceLinearization result;
		result.current[0] = Transistor.forwardAlpha * forward - reverse;
		result.current[2] = -forward + Transistor.reverseAlpha * reverse;
		result.current[1] = -(result.current[0] + result.current[2]);
		result.derivative[0][0] = gr;
		result.derivative[0][1] = Transistor.forwardAlpha * gf - gr;
		result.derivative[0][2] = -Transistor.forwardAlpha * gf;
		result.derivative[2][0] = -Transistor.reverseAlpha * gr;
		result.derivative[2][1] = -gf + Transistor.reverseAlpha * gr;
		result.derivative[2][2] = gf;
		for (int column = 0; column < 3; ++column)
			result.derivative[1][column] = -result.derivative[0][column] -
				result.derivative[2][column];
		return result;
	}

	std::array<double, UnknownCount> Evaluate(
		const std::array<double, UnknownCount>& x, double inputVoltage,
		const Histories& histories, double timeStep, bool dcMode) const
	{
		const auto q1 = NpnLinearization(x[CollectorQ1BaseQ2], x[BaseQ1],
			x[EmitterQ1]);
		const auto q2 = NpnLinearization(SupplyVoltage, x[CollectorQ1BaseQ2],
			x[EmitterQ2]);
		const double inputCurrent = dcMode ? 0.0 :
			InputCapacitance / timeStep * (x[BaseQ1] - inputVoltage) +
			histories.input;
		const double millerCurrent = dcMode ? 0.0 :
			MillerCapacitance / timeStep *
				(x[BaseQ1] - x[CollectorQ1BaseQ2]) + histories.miller;
		const double outputCurrent = dcMode ? 0.0 :
			OutputCapacitance / timeStep * (x[EmitterQ2] - x[Output]) +
			histories.output;
		return {
			(x[BaseQ1] - SupplyVoltage) / BaseBiasResistance +
				x[BaseQ1] / BaseGroundResistance + inputCurrent +
				millerCurrent + q1.current[1],
			x[EmitterQ1] / Q1EmitterResistance + q1.current[2],
			(x[CollectorQ1BaseQ2] - SupplyVoltage) / Q1CollectorResistance -
				millerCurrent + q1.current[0] + q2.current[1],
			x[EmitterQ2] / Q2EmitterResistance + outputCurrent + q2.current[2],
			x[Output] / ToneInputResistance - outputCurrent};
	}

	void BuildJacobian(double matrix[UnknownCount][UnknownCount],
		double timeStep, bool dcMode) const
	{
		const auto q1 = NpnLinearization(unknowns_[CollectorQ1BaseQ2],
			unknowns_[BaseQ1], unknowns_[EmitterQ1]);
		const auto q2 = NpnLinearization(SupplyVoltage,
			unknowns_[CollectorQ1BaseQ2], unknowns_[EmitterQ2]);
		const double inputConductance = dcMode ? 0.0 :
			InputCapacitance / timeStep;
		const double millerConductance = dcMode ? 0.0 :
			MillerCapacitance / timeStep;
		const double outputConductance = dcMode ? 0.0 :
			OutputCapacitance / timeStep;

		matrix[BaseQ1][BaseQ1] = 1.0 / BaseBiasResistance +
			1.0 / BaseGroundResistance + inputConductance +
			millerConductance + q1.derivative[1][1];
		matrix[BaseQ1][EmitterQ1] = q1.derivative[1][2];
		matrix[BaseQ1][CollectorQ1BaseQ2] = -millerConductance +
			q1.derivative[1][0];

		matrix[EmitterQ1][BaseQ1] = q1.derivative[2][1];
		matrix[EmitterQ1][EmitterQ1] = 1.0 / Q1EmitterResistance +
			q1.derivative[2][2];
		matrix[EmitterQ1][CollectorQ1BaseQ2] = q1.derivative[2][0];

		matrix[CollectorQ1BaseQ2][BaseQ1] = -millerConductance +
			q1.derivative[0][1];
		matrix[CollectorQ1BaseQ2][EmitterQ1] = q1.derivative[0][2];
		matrix[CollectorQ1BaseQ2][CollectorQ1BaseQ2] =
			1.0 / Q1CollectorResistance + millerConductance +
			q1.derivative[0][0] + q2.derivative[1][1];
		matrix[CollectorQ1BaseQ2][EmitterQ2] = q2.derivative[1][2];

		matrix[EmitterQ2][CollectorQ1BaseQ2] = q2.derivative[2][1];
		matrix[EmitterQ2][EmitterQ2] = 1.0 / Q2EmitterResistance +
			outputConductance + q2.derivative[2][2];
		matrix[EmitterQ2][Output] = -outputConductance;

		matrix[Output][EmitterQ2] = -outputConductance;
		matrix[Output][Output] = 1.0 / ToneInputResistance +
			outputConductance;
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
				for (std::size_t i = column; i < UnknownCount; ++i)
					std::swap(matrix[column][i], matrix[pivot][i]);
				std::swap(rightHandSide[column], rightHandSide[pivot]);
			}
			for (std::size_t row = column + 1; row < UnknownCount; ++row)
			{
				const double factor = matrix[row][column] /
					matrix[column][column];
				for (std::size_t i = column + 1; i < UnknownCount; ++i)
					matrix[row][i] -= factor * matrix[column][i];
				rightHandSide[row] -= factor * rightHandSide[column];
			}
		}
		for (std::size_t reverse = UnknownCount; reverse-- > 0;)
		{
			double value = rightHandSide[reverse];
			for (std::size_t column = reverse + 1; column < UnknownCount;
				++column)
				value -= matrix[reverse][column] * rightHandSide[column];
			rightHandSide[reverse] = value / matrix[reverse][reverse];
		}
		return true;
	}

	void SolveOperatingPoint()
	{
		const Histories histories{};
		for (int iteration = 0; iteration < 80; ++iteration)
		{
			const auto residual = Evaluate(unknowns_, 0.0, histories, 1.0, true);
			double maximumResidual = 0.0;
			for (double value : residual)
				maximumResidual = std::max(maximumResidual, std::abs(value));
			if (maximumResidual < NewtonTolerance)
				break;
			double matrix[UnknownCount][UnknownCount]{};
			std::array<double, UnknownCount> correction{};
			BuildJacobian(matrix, 1.0, true);
			for (std::size_t row = 0; row < UnknownCount; ++row)
				correction[row] = -residual[row];
			if (!SolveLinear(matrix, correction))
				break;
			double damping = 1.0;
			for (double value : correction)
			{
				if (std::abs(value) > MaximumJunctionStep)
					damping = std::min(damping,
						MaximumJunctionStep / std::abs(value));
			}
			for (std::size_t i = 0; i < UnknownCount; ++i)
				unknowns_[i] += damping * correction[i];
		}
		unknowns_[Output] = 0.0;
	}

	static constexpr double SupplyVoltage = 25.0;
	static constexpr double BaseBiasResistance = 390000.0;
	static constexpr double BaseGroundResistance = 33000.0;
	static constexpr double Q1CollectorResistance = 12000.0;
	static constexpr double Q1EmitterResistance = 1500.0;
	static constexpr double Q2EmitterResistance = 4700.0;
	static constexpr double ToneInputResistance = 64000.0;
	static constexpr double InputCapacitance = 0.22e-6;
	static constexpr double MillerCapacitance = 330.0e-12;
	static constexpr double OutputCapacitance = 5.0e-6;
	static constexpr double Log2E = 1.4426950408889634074;
	static constexpr int MaximumNewtonIterations = 7;
	static constexpr double NewtonTolerance = 1.0e-8;
	static constexpr double RealtimeRefinementTolerance = 2.0e-6;
	static constexpr double RealtimeFailureTolerance = 2.0e-4;
	static constexpr double MaximumNewtonStep = 8.0;
	static constexpr double MaximumJunctionStep = 0.30;
	static constexpr DeviceParameters Transistor{
		2.0e-14, 0.02585, 220.0 / 221.0, 3.0 / 4.0};

	double sampleRate_ = 96000.0;
	std::array<double, UnknownCount> unknowns_{};
	CapacitorCompanion inputCapacitor_{};
	CapacitorCompanion millerCapacitor_{};
	CapacitorCompanion outputCapacitor_{};
	std::uint64_t solverFailures_{};
	int maximumIterationsUsed_{};
	std::uint64_t totalIterations_{};
	std::uint64_t processedSamples_{};
};

} // namespace tfdsp
