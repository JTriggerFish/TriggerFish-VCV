#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include <Eigen/Dense>

#include "tfdsp/approx.hpp"

namespace tfdsp
{

// Transistor-level real-time model of one Figure 11-9 Peterson 80 W power
// module. The topology, rather than a fitted transfer curve, produces overload:
// Q1 is the globally-fed-back error stage, Q2 is the class-A transformer
// driver, and the transformer drives an asymmetric pair of matched PNP
// germanium devices. The stiff transformer/junction loop deliberately retains
// backward-Euler capacitor companions: undamped trapezoidal companions excite
// a non-physical sample-to-sample mode in this nonlinear solve. The preceding
// Peterson voltage/tone stages use trapezoidal companions where their measured
// audio-band response benefits and no such stiff loop exists.
class PetersonPowerAmplifier
{
public:
	struct Result
	{
		double voltage{};
		double positiveRailCurrent{};
		double negativeRailCurrent{};
		int iterations{};
		bool converged{};
	};

	void SetSampleRate(double sampleRate)
	{
		sampleRate_ = std::clamp(sampleRate, 32000.0, 3072000.0);
		ConfigureLoad();
		Reset(35.0);
	}

	void Reset(double railVoltage = 35.0)
	{
		railVoltage = std::clamp(railVoltage, 20.0, 40.0);
		unknowns_ = {0.68, 0.05, 1.25, 0.58, 0.0,
			railVoltage, 0.0, railVoltage - 0.01,
			-railVoltage + 0.01, 0.0};
		SolveOperatingPoint(railVoltage);
		inputCapacitor_.Reset(unknowns_[BaseQ1]);
		feedbackCapacitor_.Reset(unknowns_[EmitterQ1] - unknowns_[Output]);
		driverMillerCapacitor_.Reset(unknowns_[CollectorQ1BaseQ2] -
			(PositiveRegulatedRail - unknowns_[PrimaryVoltage]));
		upperMillerCapacitor_.Reset(unknowns_[UpperBaseNode] -
			unknowns_[Output]);
		lowerMillerCapacitor_.Reset(unknowns_[LowerBaseNode] -
			unknowns_[LowerCollector]);
		// The class-A DC primary current is the operating offset of the gapped
		// driver transformer. Integrate audio flux around that operating point.
		// Core saturation is intentionally not synthesized: the service data do
		// not identify it, and the transformer remains in its designed linear
		// region while the transistor/output circuit produces overload.
		primaryFlux_ = 0.0;
		loadState_.setZero();
		previousLoadVoltage_ = unknowns_[Output];
		solverFailures_ = 0;
		maximumIterationsUsed_ = 0;
		maximumRejectedResidual_ = 0.0;
		totalIterations_ = 0;
		processedSamples_ = 0;
	}

	Result Step(double inputVoltage, double railVoltage)
	{
		inputVoltage = std::isfinite(inputVoltage) ? inputVoltage : 0.0;
		railVoltage = std::clamp(std::isfinite(railVoltage) ? railVoltage : 35.0,
			20.0, 40.0);
		const double timeStep = 1.0 / sampleRate_;
		const CompanionHistories histories{
			inputCapacitor_.History(InputCapacitance, timeStep),
			feedbackCapacitor_.History(FeedbackCapacitance, timeStep),
			driverMillerCapacitor_.History(DriverMillerCapacitance, timeStep),
			upperMillerCapacitor_.History(OutputMillerCapacitance, timeStep),
			lowerMillerCapacitor_.History(OutputMillerCapacitance, timeStep)};
		const Eigen::Vector3d loadHistory = loadTransition_ * loadState_ +
			loadInput_ * previousLoadVoltage_;
		const double loadCurrentSlope = loadInput_(0) + 1.0 / BleederResistance;

		bool converged = false;
		int iterations = 0;
		int iterationsPerformed = 0;
		for (; iterations < MaximumNewtonIterations; ++iterations)
		{
			++iterationsPerformed;
			const auto residual = EvaluateCircuit(unknowns_, inputVoltage,
				railVoltage, histories, loadHistory(0), loadCurrentSlope,
				timeStep, false);
			double maximumResidual = 0.0;
			for (std::size_t row = 0; row < UnknownCount; ++row)
			{
				maximumResidual = std::max(maximumResidual,
					std::abs(residual[row]));
			}
			if (maximumResidual < NewtonTolerance)
			{
				converged = true;
				break;
			}
			double matrix[UnknownCount][UnknownCount]{};
			std::array<double, UnknownCount> correction{};
			BuildJacobian(matrix, railVoltage, loadCurrentSlope, timeStep,
				false);
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
			const std::array<double, 8> junctionCorrections{
				correction[BaseQ1] - correction[EmitterQ1],
				correction[BaseQ1] - correction[CollectorQ1BaseQ2],
				correction[CollectorQ1BaseQ2] - correction[EmitterQ2],
				correction[CollectorQ1BaseQ2] + correction[PrimaryVoltage],
				correction[UpperBaseNode] - correction[UpperEmitter],
				correction[UpperBaseNode] - correction[Output],
				correction[LowerBaseNode] - correction[Output],
				correction[LowerBaseNode] - correction[LowerCollector]};
			for (double junctionCorrection : junctionCorrections)
			{
				if (std::abs(junctionCorrection) > MaximumJunctionStep)
					damping = std::min(damping, MaximumJunctionStep /
						std::abs(junctionCorrection));
			}
			// The ordinary audio-rate path converges before this is needed. Near
			// saturation the germanium exponentials can make an undamped Newton
			// direction overshoot, so after two unsuccessful updates require the
			// candidate to reduce the actual KCL residual.
			if (iterations >= 2)
			{
				for (int lineSearch = 0; lineSearch < 5; ++lineSearch)
				{
					auto candidate = unknowns_;
					for (std::size_t i = 0; i < UnknownCount; ++i)
						candidate[i] += damping * correction[i];
					const auto candidateResidual = EvaluateCircuit(candidate,
						inputVoltage, railVoltage, histories, loadHistory(0),
						loadCurrentSlope, timeStep, false);
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
			// Always re-evaluate KCL after a Newton correction. Treating one
			// undamped correction as converged left a small absolute residual that
			// was insignificant near overload but became dominant crossover-like
			// buzz on quiet and high-register notes.
		}

		if (!converged)
		{
			const auto residual = EvaluateCircuit(unknowns_,
				inputVoltage, railVoltage, histories, loadHistory(0),
				loadCurrentSlope, timeStep, false);
			double maximumResidual = 0.0;
			for (const auto& value : residual)
				maximumResidual = std::max(maximumResidual, std::abs(value));
			converged = maximumResidual < RealtimeFailureTolerance;
			maximumRejectedResidual_ = std::max(maximumRejectedResidual_,
				maximumResidual);
		}
		if (!converged)
			++solverFailures_;
		maximumIterationsUsed_ = std::max(maximumIterationsUsed_,
			iterationsPerformed);
		totalIterations_ += static_cast<std::uint64_t>(iterationsPerformed);
		++processedSamples_;

		const double collectorQ2 = PositiveRegulatedRail -
			unknowns_[PrimaryVoltage];
		inputCapacitor_.Commit(unknowns_[BaseQ1] - inputVoltage,
			InputCapacitance, timeStep);
		feedbackCapacitor_.Commit(unknowns_[EmitterQ1] - unknowns_[Output],
			FeedbackCapacitance, timeStep);
		driverMillerCapacitor_.Commit(unknowns_[CollectorQ1BaseQ2] - collectorQ2,
			DriverMillerCapacitance, timeStep);
		upperMillerCapacitor_.Commit(unknowns_[UpperBaseNode] -
			unknowns_[Output],
			OutputMillerCapacitance, timeStep);
		lowerMillerCapacitor_.Commit(unknowns_[LowerBaseNode] -
			unknowns_[LowerCollector],
			OutputMillerCapacitance, timeStep);
		primaryFlux_ += timeStep * unknowns_[PrimaryVoltage];
		loadState_ = loadHistory + loadInput_ * unknowns_[Output];
		previousLoadVoltage_ = unknowns_[Output];

		const double positiveCurrent = std::max(0.0,
			(railVoltage - unknowns_[UpperEmitter]) / OutputEmitterResistance);
		const double negativeCurrent = std::max(0.0,
			(unknowns_[LowerCollector] + railVoltage) / OutputEmitterResistance);
		return {std::isfinite(unknowns_[Output]) ? unknowns_[Output] : 0.0,
			positiveCurrent, negativeCurrent, iterationsPerformed, converged};
	}

	std::size_t SolverFailures() const { return solverFailures_; }
	int MaximumIterationsUsed() const { return maximumIterationsUsed_; }
	double OutputVoltage() const { return unknowns_[Output]; }
	double MaximumRejectedResidual() const { return maximumRejectedResidual_; }
	double AverageIterations() const
	{
		return processedSamples_ == 0 ? 0.0 :
			static_cast<double>(totalIterations_) /
			static_cast<double>(processedSamples_);
	}

	static std::complex<double> ElectricalLoadImpedance(double frequency)
	{
		if (!(frequency > 0.0) || !std::isfinite(frequency))
			return {LoadDcResistance, 0.0};
		const double angularFrequency = TwoPi * frequency;
		const std::complex<double> mechanicalImpedance{
			LoadMechanicalResistance,
			angularFrequency * LoadMovingMass -
				1.0 / (angularFrequency * LoadCompliance)};
		return std::complex<double>{LoadDcResistance,
			angularFrequency * LoadVoiceCoilInductance} +
			LoadForceFactor * LoadForceFactor / mechanicalImpedance;
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
	};

	template <std::size_t N>
	friend Dual<N> operator+(const Dual<N>& a, const Dual<N>& b)
	{
		Dual<N> result(a.value + b.value);
		for (std::size_t i = 0; i < N; ++i)
			result.derivative[i] = a.derivative[i] + b.derivative[i];
		return result;
	}
	template <std::size_t N>
	friend Dual<N> operator-(const Dual<N>& a, const Dual<N>& b)
	{
		Dual<N> result(a.value - b.value);
		for (std::size_t i = 0; i < N; ++i)
			result.derivative[i] = a.derivative[i] - b.derivative[i];
		return result;
	}
	template <std::size_t N>
	friend Dual<N> operator-(const Dual<N>& value)
	{
		Dual<N> result(-value.value);
		for (std::size_t i = 0; i < N; ++i)
			result.derivative[i] = -value.derivative[i];
		return result;
	}
	template <std::size_t N>
	friend Dual<N> operator*(const Dual<N>& a, const Dual<N>& b)
	{
		Dual<N> result(a.value * b.value);
		for (std::size_t i = 0; i < N; ++i)
			result.derivative[i] = a.derivative[i] * b.value +
				a.value * b.derivative[i];
		return result;
	}
	template <std::size_t N>
	friend Dual<N> operator/(const Dual<N>& a, const Dual<N>& b)
	{
		Dual<N> result(a.value / b.value);
		const double inverseSquare = 1.0 / (b.value * b.value);
		for (std::size_t i = 0; i < N; ++i)
			result.derivative[i] = (a.derivative[i] * b.value -
				a.value * b.derivative[i]) * inverseSquare;
		return result;
	}
	template <std::size_t N>
	friend Dual<N> operator+(const Dual<N>& a, double b)
	{
		return a + Dual<N>(b);
	}
	template <std::size_t N>
	friend Dual<N> operator+(double a, const Dual<N>& b)
	{
		return Dual<N>(a) + b;
	}
	template <std::size_t N>
	friend Dual<N> operator-(const Dual<N>& a, double b)
	{
		return a - Dual<N>(b);
	}
	template <std::size_t N>
	friend Dual<N> operator-(double a, const Dual<N>& b)
	{
		return Dual<N>(a) - b;
	}
	template <std::size_t N>
	friend Dual<N> operator*(const Dual<N>& a, double b)
	{
		return a * Dual<N>(b);
	}
	template <std::size_t N>
	friend Dual<N> operator*(double a, const Dual<N>& b)
	{
		return Dual<N>(a) * b;
	}
	template <std::size_t N>
	friend Dual<N> operator/(const Dual<N>& a, double b)
	{
		return a / Dual<N>(b);
	}
	template <std::size_t N>
	friend Dual<N> operator/(double a, const Dual<N>& b)
	{
		return Dual<N>(a) / b;
	}

	template <std::size_t N>
	static Dual<N> LimitedExp(const Dual<N>& argument)
	{
		// Below -20 the exponential contributes less than 2.1e-9 of Is.
		// Returning zero is below the KCL tolerance even for the leaky germanium
		// proxy and avoids costly exp() calls for reverse-biased junctions.
		if (argument.value <= -20.0)
			return Dual<N>(0.0);
		const double limited = std::min(argument.value, 40.0);
		const double exponential = std::exp(limited);
		Dual<N> result(exponential);
		if (argument.value < 40.0)
		{
			for (std::size_t i = 0; i < N; ++i)
				result.derivative[i] = exponential * argument.derivative[i];
		}
		return result;
	}
	static double LimitedExp(double argument)
	{
		return argument <= -20.0 ? 0.0 :
			static_cast<double>(Exp2Taylor5(static_cast<float>(
				std::min(argument, 40.0) * Log2E)));
	}

	template <typename T>
	struct TerminalCurrents
	{
		T collector;
		T base;
		T emitter;
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

	struct RealtimeDeviceLinearizations
	{
		DeviceLinearization q1;
		DeviceLinearization q2;
		DeviceLinearization upper;
		DeviceLinearization lower;
	};

	static DeviceLinearization NpnLinearization(double collector, double base,
		double emitter, const DeviceParameters& parameters)
	{
		const double forwardArgument = (base - emitter) /
			parameters.thermalVoltage;
		const double reverseArgument = (base - collector) /
			parameters.thermalVoltage;
		const double forwardExp = LimitedExp(forwardArgument);
		const double reverseExp = LimitedExp(reverseArgument);
		const double forward = parameters.saturationCurrent /
			parameters.forwardAlpha * (forwardExp - 1.0);
		const double reverse = parameters.saturationCurrent /
			parameters.reverseAlpha * (reverseExp - 1.0);
		const double gf = parameters.saturationCurrent /
			(parameters.forwardAlpha * parameters.thermalVoltage) * forwardExp;
		const double gr = parameters.saturationCurrent /
			(parameters.reverseAlpha * parameters.thermalVoltage) * reverseExp;
		DeviceLinearization result;
		result.current[0] = parameters.forwardAlpha * forward - reverse;
		result.current[2] = -forward + parameters.reverseAlpha * reverse;
		result.current[1] = -(result.current[0] + result.current[2]);
		result.derivative[0][0] = gr;
		result.derivative[0][1] = parameters.forwardAlpha * gf - gr;
		result.derivative[0][2] = -parameters.forwardAlpha * gf;
		result.derivative[2][0] = -parameters.reverseAlpha * gr;
		result.derivative[2][1] = -gf + parameters.reverseAlpha * gr;
		result.derivative[2][2] = gf;
		for (int column = 0; column < 3; ++column)
			result.derivative[1][column] = -result.derivative[0][column] -
				result.derivative[2][column];
		return result;
	}

	static DeviceLinearization PnpLinearization(double collector, double base,
		double emitter, const DeviceParameters& parameters)
	{
		auto result = NpnLinearization(-collector, -base, -emitter,
			parameters);
		for (double& current : result.current)
			current = -current;
		// Negating both terminal voltage and current leaves dI/dV unchanged.
		return result;
	}

	template <typename T>
	static TerminalCurrents<T> Npn(const T& collector, const T& base,
		const T& emitter, const DeviceParameters& parameters)
	{
		const T forward = (parameters.saturationCurrent /
			parameters.forwardAlpha) *
			(LimitedExp((base - emitter) / parameters.thermalVoltage) - 1.0);
		const T reverse = (parameters.saturationCurrent /
			parameters.reverseAlpha) *
			(LimitedExp((base - collector) / parameters.thermalVoltage) - 1.0);
		const T collectorCurrent = parameters.forwardAlpha * forward - reverse;
		const T emitterCurrent = -forward + parameters.reverseAlpha * reverse;
		return {collectorCurrent, -(collectorCurrent + emitterCurrent),
			emitterCurrent};
	}

	template <typename T>
	static TerminalCurrents<T> Pnp(const T& collector, const T& base,
		const T& emitter, const DeviceParameters& parameters)
	{
		const auto mirrored = Npn(-collector, -base, -emitter, parameters);
		return {-mirrored.collector, -mirrored.base, -mirrored.emitter};
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
			return -capacitance * previousVoltage / timeStep;
		}
		void Commit(double voltage, double capacitance, double timeStep)
		{
			previousCurrent = capacitance * voltage / timeStep +
				History(capacitance, timeStep);
			previousVoltage = voltage;
		}
	};

	struct CompanionHistories
	{
		double input;
		double feedback;
		double driverMiller;
		double upperMiller;
		double lowerMiller;
	};

	enum Unknown : std::size_t
	{
		BaseQ1,
		EmitterQ1,
		CollectorQ1BaseQ2,
		EmitterQ2,
		PrimaryVoltage,
		UpperBaseNode,
		LowerBaseNode,
		UpperEmitter,
		LowerCollector,
		Output,
		UnknownCount
	};

	template <typename T>
	static T CapacitorCurrent(const T& voltage, double capacitance,
		double timeStep, double history)
	{
		return (capacitance / timeStep) * voltage + history;
	}

	template <typename T>
	static T UpperWindingSource(double railVoltage, const T& primaryVoltage)
	{
		return railVoltage + primaryVoltage / TransformerTurnsRatio;
	}
	template <typename T>
	static T LowerWindingSource(const T& outputVoltage,
		const T& primaryVoltage)
	{
		// The lower DTG110B is also PNP. Its secondary is referenced to the
		// audio output and wound with the opposite phase, making this device a
		// common-collector sink while the upper PNP is common-emitter.
		return outputVoltage - primaryVoltage / TransformerTurnsRatio;
	}

	template <typename T>
	static T MagnetizingCurrent(const T& flux)
	{
		// Preserve the transformer's dynamic magnetising branch and its reflected
		// loading without inventing an unmeasured saturation knee. A nonlinear
		// core belongs here only after identification from an original T1.
		return flux / MagnetizingInductance;
	}
	template <typename T>
	std::array<T, UnknownCount> EvaluateCircuit(
		const std::array<T, UnknownCount>& x, double inputVoltage,
		double railVoltage, const CompanionHistories& histories,
		double loadHistory, double loadCurrentSlope, double timeStep,
		bool dcMode) const
	{
		const T collectorQ2 = PositiveRegulatedRail - x[PrimaryVoltage];
		const T upperBase = x[UpperBaseNode];
		const T lowerBase = x[LowerBaseNode];
		const T upperWindingSource = UpperWindingSource(railVoltage,
			x[PrimaryVoltage]);
		const T lowerWindingSource = LowerWindingSource(x[Output],
			x[PrimaryVoltage]);
		TerminalCurrents<T> q1;
		TerminalCurrents<T> q2;
		TerminalCurrents<T> upper;
		TerminalCurrents<T> lower;
		if constexpr (std::is_same_v<T, double>)
		{
			realtimeDevices_.q1 = NpnLinearization(x[CollectorQ1BaseQ2],
				x[BaseQ1], x[EmitterQ1], SmallSignalNpn);
			realtimeDevices_.q2 = NpnLinearization(collectorQ2,
				x[CollectorQ1BaseQ2], x[EmitterQ2], SmallSignalNpn);
			realtimeDevices_.upper = PnpLinearization(x[Output], upperBase,
				x[UpperEmitter], GermaniumPowerDevice);
			realtimeDevices_.lower = PnpLinearization(x[LowerCollector],
				lowerBase, x[Output], GermaniumPowerDevice);
			auto currents = [](const DeviceLinearization& device)
			{
				return TerminalCurrents<double>{device.current[0],
					device.current[1], device.current[2]};
			};
			q1 = currents(realtimeDevices_.q1);
			q2 = currents(realtimeDevices_.q2);
			upper = currents(realtimeDevices_.upper);
			lower = currents(realtimeDevices_.lower);
		}
		else
		{
			q1 = Npn(x[CollectorQ1BaseQ2], x[BaseQ1], x[EmitterQ1],
				SmallSignalNpn);
			q2 = Npn(collectorQ2, x[CollectorQ1BaseQ2], x[EmitterQ2],
				SmallSignalNpn);
			upper = Pnp(x[Output], upperBase, x[UpperEmitter],
				GermaniumPowerDevice);
			lower = Pnp(x[LowerCollector], lowerBase, x[Output],
				GermaniumPowerDevice);
		}

		const T inputCapCurrent = dcMode ? T(0.0) : CapacitorCurrent(
			x[BaseQ1] - inputVoltage, InputCapacitance, timeStep,
			histories.input);
		const T feedbackCapCurrent = dcMode ? T(0.0) : CapacitorCurrent(
			x[EmitterQ1] - x[Output], FeedbackCapacitance, timeStep,
			histories.feedback);
		const T driverMillerCurrent = dcMode ? T(0.0) : CapacitorCurrent(
			x[CollectorQ1BaseQ2] - collectorQ2, DriverMillerCapacitance,
			timeStep, histories.driverMiller);
		const T upperMillerCurrent = dcMode ? T(0.0) : CapacitorCurrent(
			upperBase - x[Output], OutputMillerCapacitance, timeStep,
			histories.upperMiller);
		const T lowerMillerCurrent = dcMode ? T(0.0) : CapacitorCurrent(
			lowerBase - x[LowerCollector], OutputMillerCapacitance, timeStep,
			histories.lowerMiller);

		const T upperWindingCurrent =
			(upperBase - upperWindingSource) / SecondaryWindingResistance;
		const T lowerWindingCurrent =
			(lowerBase - lowerWindingSource) / SecondaryWindingResistance;
		const T flux = dcMode ? T(primaryFlux_) : T(primaryFlux_) +
			timeStep * x[PrimaryVoltage];
		const T primaryCurrent = dcMode ? T(operatingMagnetizingCurrent_) :
			operatingMagnetizingCurrent_ + MagnetizingCurrent(flux) -
				(upperWindingCurrent - lowerWindingCurrent) /
				TransformerTurnsRatio;

		std::array<T, UnknownCount> residual;
		residual[BaseQ1] = q1.base +
			(x[BaseQ1] - x[CollectorQ1BaseQ2]) / InputBiasResistance +
			inputCapCurrent;
		residual[EmitterQ1] = q1.emitter + x[EmitterQ1] /
			FeedbackGroundResistance +
			(x[EmitterQ1] - x[Output]) / FeedbackResistance +
			feedbackCapCurrent;
		residual[CollectorQ1BaseQ2] = q1.collector + q2.base +
			(x[CollectorQ1BaseQ2] - PositiveRegulatedRail) /
				Q1CollectorResistance + driverMillerCurrent;
		residual[EmitterQ2] = q2.emitter + x[EmitterQ2] /
			Q2EmitterResistance +
			(x[EmitterQ2] - collectorQ2) / Q2CollectorEmitterResistance;
		residual[PrimaryVoltage] = dcMode ? x[PrimaryVoltage] :
			q2.collector +
			(collectorQ2 - x[EmitterQ2]) / Q2CollectorEmitterResistance +
			(collectorQ2 - PositiveRegulatedRail) / CoreLossResistance -
			driverMillerCurrent - primaryCurrent;
		residual[UpperBaseNode] = upper.base +
			(upperBase - x[UpperEmitter]) / BaseEmitterShuntResistance +
			upperMillerCurrent + upperWindingCurrent;
		residual[LowerBaseNode] = lower.base +
			(lowerBase - x[Output]) / BaseEmitterShuntResistance +
			lowerMillerCurrent + lowerWindingCurrent;
		residual[UpperEmitter] = upper.emitter +
			(x[UpperEmitter] - railVoltage) / OutputEmitterResistance +
			(x[UpperEmitter] - upperBase) / BaseEmitterShuntResistance;
		residual[LowerCollector] = lower.collector +
			(x[LowerCollector] + railVoltage) / OutputEmitterResistance -
			lowerMillerCurrent;
		const T loadCurrent = dcMode ? x[Output] /
			LoadDcResistance : loadHistory + loadCurrentSlope * x[Output];
		residual[Output] = upper.collector + lower.emitter + loadCurrent +
			(x[Output] - x[EmitterQ1]) / FeedbackResistance -
			feedbackCapCurrent - upperMillerCurrent +
			(x[Output] - lowerBase) / BaseEmitterShuntResistance -
			lowerWindingCurrent;
		return residual;
	}

	void BuildJacobian(double matrix[UnknownCount][UnknownCount],
		double railVoltage, double loadCurrentSlope, double timeStep,
		bool dcMode) const
	{
		(void)railVoltage;
		auto add = [&](std::size_t row, std::size_t column, double value)
		{
			matrix[row][column] += value;
		};
		auto stampDevice = [&](const DeviceLinearization& device,
			const std::array<std::size_t, 3>& rows,
			const std::array<std::size_t, 3>& columns,
			const std::array<double, 3>& columnScale)
		{
			for (int terminal = 0; terminal < 3; ++terminal)
				for (int voltage = 0; voltage < 3; ++voltage)
					add(rows[terminal], columns[voltage],
						device.derivative[terminal][voltage] *
						columnScale[voltage]);
		};

		const auto& q1 = realtimeDevices_.q1;
		stampDevice(q1,
			{CollectorQ1BaseQ2, BaseQ1, EmitterQ1},
			{CollectorQ1BaseQ2, BaseQ1, EmitterQ1}, {1.0, 1.0, 1.0});
		const auto& q2 = realtimeDevices_.q2;
		stampDevice(q2,
			{PrimaryVoltage, CollectorQ1BaseQ2, EmitterQ2},
			{PrimaryVoltage, CollectorQ1BaseQ2, EmitterQ2}, {-1.0, 1.0, 1.0});
		const auto& upper = realtimeDevices_.upper;
		stampDevice(upper,
			{Output, UpperBaseNode, UpperEmitter},
			{Output, UpperBaseNode, UpperEmitter}, {1.0, 1.0, 1.0});
		const auto& lower = realtimeDevices_.lower;
		stampDevice(lower,
			{LowerCollector, LowerBaseNode, Output},
			{LowerCollector, LowerBaseNode, Output}, {1.0, 1.0, 1.0});

		const double inputConductance = dcMode ? 0.0 :
			InputCapacitance / timeStep;
		const double feedbackConductance = dcMode ? 0.0 :
			FeedbackCapacitance / timeStep;
		const double driverMillerConductance = dcMode ? 0.0 :
			DriverMillerCapacitance / timeStep;
		const double outputMillerConductance = dcMode ? 0.0 :
			OutputMillerCapacitance / timeStep;
		const double inputBiasConductance = 1.0 / InputBiasResistance;
		const double q1CollectorConductance = 1.0 / Q1CollectorResistance;
		const double feedbackGroundConductance = 1.0 /
			FeedbackGroundResistance;
		const double feedbackConductanceR = 1.0 / FeedbackResistance;
		const double q2EmitterConductance = 1.0 / Q2EmitterResistance;
		const double q2ShuntConductance = 1.0 /
			Q2CollectorEmitterResistance;
		const double coreLossConductance = 1.0 / CoreLossResistance;
		const double windingConductance = 1.0 / SecondaryWindingResistance;
		const double baseShuntConductance = 1.0 /
			BaseEmitterShuntResistance;
		const double emitterPathConductance = 1.0 / OutputEmitterResistance;
		const double inverseTurns = 1.0 / TransformerTurnsRatio;

		add(BaseQ1, BaseQ1, inputBiasConductance + inputConductance);
		add(BaseQ1, CollectorQ1BaseQ2, -inputBiasConductance);
		add(EmitterQ1, EmitterQ1, feedbackGroundConductance +
			feedbackConductanceR + feedbackConductance);
		add(EmitterQ1, Output, -feedbackConductanceR - feedbackConductance);
		add(CollectorQ1BaseQ2, CollectorQ1BaseQ2,
			q1CollectorConductance + driverMillerConductance);
		add(CollectorQ1BaseQ2, PrimaryVoltage, driverMillerConductance);
		add(EmitterQ2, EmitterQ2, q2EmitterConductance + q2ShuntConductance);
		add(EmitterQ2, PrimaryVoltage, q2ShuntConductance);

		if (dcMode)
		{
			for (std::size_t column = 0; column < UnknownCount; ++column)
				matrix[PrimaryVoltage][column] = 0.0;
			matrix[PrimaryVoltage][PrimaryVoltage] = 1.0;
		}
		else
		{
			add(PrimaryVoltage, EmitterQ2, -q2ShuntConductance);
			add(PrimaryVoltage, PrimaryVoltage, -q2ShuntConductance -
				coreLossConductance - driverMillerConductance);
			add(PrimaryVoltage, CollectorQ1BaseQ2,
				-driverMillerConductance);
			add(PrimaryVoltage, UpperBaseNode,
				windingConductance * inverseTurns);
			add(PrimaryVoltage, LowerBaseNode,
				-windingConductance * inverseTurns);
			add(PrimaryVoltage, Output,
				windingConductance * inverseTurns);
			add(PrimaryVoltage, PrimaryVoltage,
				-2.0 * windingConductance * inverseTurns * inverseTurns);
			add(PrimaryVoltage, PrimaryVoltage,
				-timeStep / MagnetizingInductance);
		}

		add(UpperBaseNode, UpperBaseNode, baseShuntConductance +
			outputMillerConductance + windingConductance);
		add(UpperBaseNode, UpperEmitter, -baseShuntConductance);
		add(UpperBaseNode, Output, -outputMillerConductance);
		add(UpperBaseNode, PrimaryVoltage,
			-windingConductance * inverseTurns);
		add(LowerBaseNode, LowerBaseNode, baseShuntConductance +
			outputMillerConductance + windingConductance);
		add(LowerBaseNode, Output, -baseShuntConductance - windingConductance);
		add(LowerBaseNode, LowerCollector, -outputMillerConductance);
		add(LowerBaseNode, PrimaryVoltage,
			windingConductance * inverseTurns);
		add(UpperEmitter, UpperEmitter,
			emitterPathConductance + baseShuntConductance);
		add(UpperEmitter, UpperBaseNode, -baseShuntConductance);
		add(LowerCollector, LowerCollector,
			emitterPathConductance + outputMillerConductance);
		add(LowerCollector, LowerBaseNode, -outputMillerConductance);
		add(Output, Output, (dcMode ? 1.0 / LoadDcResistance :
			loadCurrentSlope) + feedbackConductanceR + feedbackConductance +
			outputMillerConductance + baseShuntConductance + windingConductance);
		add(Output, EmitterQ1, -feedbackConductanceR - feedbackConductance);
		add(Output, UpperBaseNode, -outputMillerConductance);
		add(Output, LowerBaseNode,
			-baseShuntConductance - windingConductance);
		add(Output, PrimaryVoltage, -windingConductance * inverseTurns);
	}

	void SolveOperatingPoint(double railVoltage)
	{
		const CompanionHistories histories{};
		for (int iteration = 0; iteration < 40; ++iteration)
		{
			std::array<Dual<UnknownCount>, UnknownCount> variables;
			for (std::size_t i = 0; i < UnknownCount; ++i)
				variables[i] = Dual<UnknownCount>::Variable(unknowns_[i], i);
			const auto residual = EvaluateCircuit(variables, 0.0, railVoltage,
				histories, 0.0, 1.0 / LoadDcResistance, 1.0 / sampleRate_, true);
			double matrix[UnknownCount][UnknownCount]{};
			std::array<double, UnknownCount> correction{};
			double maximumResidual = 0.0;
			for (std::size_t row = 0; row < UnknownCount; ++row)
			{
				correction[row] = -residual[row].value;
				maximumResidual = std::max(maximumResidual,
					std::abs(residual[row].value));
				for (std::size_t column = 0; column < UnknownCount; ++column)
					matrix[row][column] = residual[row].derivative[column];
			}
			if (maximumResidual < NewtonTolerance)
				break;
			if (!SolveLinear(matrix, correction))
				break;
			double damping = 1.0;
			for (double value : correction)
				if (std::abs(value) > 1.0)
					damping = std::min(damping, 1.0 / std::abs(value));
			for (std::size_t i = 0; i < UnknownCount; ++i)
				unknowns_[i] += damping * correction[i];
		}

		const double collectorQ2 = PositiveRegulatedRail -
			unknowns_[PrimaryVoltage];
		const auto q2 = Npn(collectorQ2, unknowns_[CollectorQ1BaseQ2],
			unknowns_[EmitterQ2], SmallSignalNpn);
		operatingMagnetizingCurrent_ = q2.collector +
			(collectorQ2 - unknowns_[EmitterQ2]) /
				Q2CollectorEmitterResistance +
			(collectorQ2 - PositiveRegulatedRail) / CoreLossResistance;
	}

	static bool SolveLinear(double matrix[UnknownCount][UnknownCount],
		std::array<double, UnknownCount>& right)
	{
		for (std::size_t column = 0; column < UnknownCount; ++column)
		{
			std::size_t pivot = column;
			for (std::size_t row = column + 1; row < UnknownCount; ++row)
				if (std::abs(matrix[row][column]) >
					std::abs(matrix[pivot][column]))
					pivot = row;
			if (std::abs(matrix[pivot][column]) < 1.0e-18)
				return false;
			if (pivot != column)
			{
				for (std::size_t entry = column; entry < UnknownCount; ++entry)
					std::swap(matrix[column][entry], matrix[pivot][entry]);
				std::swap(right[column], right[pivot]);
			}
			for (std::size_t row = column + 1; row < UnknownCount; ++row)
			{
				const double factor = matrix[row][column] /
					matrix[column][column];
				for (std::size_t entry = column + 1; entry < UnknownCount; ++entry)
					matrix[row][entry] -= factor * matrix[column][entry];
				right[row] -= factor * right[column];
			}
		}
		for (std::size_t reverse = UnknownCount; reverse-- > 0;)
		{
			double value = right[reverse];
			for (std::size_t column = reverse + 1; column < UnknownCount; ++column)
				value -= matrix[reverse][column] * right[column];
			right[reverse] = value / matrix[reverse][reverse];
		}
		return true;
	}

	void ConfigureLoad()
	{
		Eigen::Matrix3d continuous;
		continuous <<
			-LoadDcResistance / LoadVoiceCoilInductance,
			-LoadForceFactor / LoadVoiceCoilInductance, 0.0,
			LoadForceFactor / LoadMovingMass,
			-LoadMechanicalResistance / LoadMovingMass,
			-1.0 / (LoadMovingMass * LoadCompliance),
			0.0, 1.0, 0.0;
		Eigen::Vector3d input;
		input << 1.0 / LoadVoiceCoilInductance, 0.0, 0.0;
		const double halfStep = 0.5 / sampleRate_;
		const Eigen::Matrix3d inverse =
			(Eigen::Matrix3d::Identity() - halfStep * continuous).inverse();
		loadTransition_ = inverse *
			(Eigen::Matrix3d::Identity() + halfStep * continuous);
		loadInput_ = inverse * (halfStep * input);
	}

	static constexpr double TwoPi = 6.2831853071795864769;
	static constexpr double Log2E = 1.4426950408889634074;
	// Two updates cover normal operation; up to eight are available only on the
	// steep germanium/core knee. This removes rejected overload samples without
	// paying the extra solves in the clean region.
	static constexpr int MaximumNewtonIterations = 8;
	// A 2 uA maximum KCL residual is below the device/model uncertainty by
	// orders of magnitude and corresponds to sub-millivolt output error. Solving
	// toward floating-point noise only burns iterations in a 192 kHz loop.
	static constexpr double NewtonTolerance = 2.0e-6;
	static constexpr double RealtimeFailureTolerance = 2.0e-2;
	static constexpr double MaximumNewtonStep = 10.0;
	static constexpr double MaximumJunctionStep = 0.30;
	static constexpr double PositiveRegulatedRail = 25.0;
	static constexpr double InputCapacitance = 6.4e-6;
	static constexpr double InputBiasResistance = 68000.0;
	static constexpr double Q1CollectorResistance = 6800.0;
	static constexpr double FeedbackGroundResistance = 100.0;
	static constexpr double FeedbackResistance = 5600.0;
	static constexpr double FeedbackCapacitance = 1200.0e-12;
	static constexpr double Q2EmitterResistance = 47.0;
	static constexpr double Q2CollectorEmitterResistance = 22000.0;
	static constexpr double DriverMillerCapacitance = 100.0e-12;
	static constexpr double CoreLossResistance = 42000.0;
	static constexpr double TransformerTurnsRatio = 1.55;
	static constexpr double SecondaryWindingResistance = 2.7;
	static constexpr double MagnetizingInductance = 8.0;
	static constexpr double BaseEmitterShuntResistance = 820.0;
	static constexpr double OutputMillerCapacitance = 0.01e-6;
	static constexpr double OutputEmitterResistance = 0.50;
	static constexpr double BleederResistance = 270.0;
	static constexpr DeviceParameters SmallSignalNpn{
		4.0e-14, 0.02585, 180.0 / 181.0, 4.0 / 5.0};
	// Reduced Ebers-Moll fit to the low-current/output region of the offline
	// DTG110B Gummel-Poon proxy. Copying that model's IS into a single-exponential
	// device left the transformer-driven class-B crossover about 10--20 dB too
	// dirty on ordinary high notes. These effective values reproduce the
	// reference's 0.05 V H2--H4 within the checked tolerance; they are not claimed
	// manufacturer parameters for the house-numbered 120725 device.
	static constexpr DeviceParameters GermaniumPowerDevice{
		5.0e-3, 0.036, 90.0 / 91.0, 3.0 / 4.0};

	static constexpr double LoadDcResistance = 12.8;
	static constexpr double LoadVoiceCoilInductance = 0.00055;
	static constexpr double LoadMovingMass = 0.018;
	static constexpr double LoadCompliance = 0.0002502;
	static constexpr double LoadMechanicalResistance = 2.12;
	static constexpr double LoadForceFactor = 10.50;

	double sampleRate_ = 192000.0;
	std::array<double, UnknownCount> unknowns_{};
	CapacitorCompanion inputCapacitor_;
	CapacitorCompanion feedbackCapacitor_;
	CapacitorCompanion driverMillerCapacitor_;
	CapacitorCompanion upperMillerCapacitor_;
	CapacitorCompanion lowerMillerCapacitor_;
	double primaryFlux_{};
	double operatingMagnetizingCurrent_{};
	Eigen::Matrix3d loadTransition_ = Eigen::Matrix3d::Identity();
	Eigen::Vector3d loadInput_ = Eigen::Vector3d::Zero();
	Eigen::Vector3d loadState_ = Eigen::Vector3d::Zero();
	double previousLoadVoltage_{};
	std::size_t solverFailures_{};
	int maximumIterationsUsed_{};
	double maximumRejectedResidual_{};
	mutable RealtimeDeviceLinearizations realtimeDevices_{};
	std::uint64_t totalIterations_{};
	std::uint64_t processedSamples_{};
};

} // namespace tfdsp
