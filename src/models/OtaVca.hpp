#pragma once

#include <algorithm>
#include <cmath>

namespace tfdsp
{

// Quasi-static differential-pair OTA. The interface deliberately uses circuit
// units: differential input voltage, control current, output-node voltage, and
// output current. Envelope voltages, Rack scaling, coupling capacitors, and
// application-specific loads belong in wrappers such as Tb303Vca.
//
// For a matched pair the exact large-signal law is
//
//   i_out = i_abc * tanh(v_diff / (2 V_T)),
//
// with g_m = i_abc / (2 V_T) at small signal. Efficiency accounts for base
// current and current-mirror transfer loss. Optional mismatch and Early-effect
// terms exist for reference/sensitivity tests; the production default leaves
// them disabled because no sufficiently broad original-BA662 dataset exists.
class OtaVcaCore
{
public:
	struct Configuration
	{
		double thermalVoltage{0.02585};
		double currentTransferEfficiency{0.85};
		double inputOffsetVolts{};
		double mirrorImbalance{};
		double earlyVoltage{100.0};
		bool modelEarlyEffect{};
	};

	OtaVcaCore() = default;
	explicit OtaVcaCore(Configuration configuration)
		: _configuration(configuration)
	{
	}

	void Configure(Configuration configuration)
	{
		_configuration = configuration;
	}

	const Configuration& GetConfiguration() const
	{
		return _configuration;
	}

	double SmallSignalTransconductance(double controlCurrentAmps) const
	{
		if (!std::isfinite(controlCurrentAmps) || controlCurrentAmps <= 0.0)
			return 0.0;
		const double thermalVoltage = SafeThermalVoltage();
		return SafeEfficiency() * std::min(controlCurrentAmps, MaximumControlCurrent) /
			(2.0 * thermalVoltage);
	}

	double ProcessCurrent(double differentialInputVolts,
		double controlCurrentAmps, double outputNodeVolts = 0.0) const
	{
		if (!std::isfinite(differentialInputVolts) ||
			!std::isfinite(controlCurrentAmps) ||
			!std::isfinite(outputNodeVolts) || controlCurrentAmps <= 0.0)
			return 0.0;

		const double control = SafeEfficiency() *
			std::clamp(controlCurrentAmps, 0.0, MaximumControlCurrent);
		const double argument = (differentialInputVolts +
			SafeInputOffset()) / (2.0 * SafeThermalVoltage());
		double output = control * std::tanh(argument);

		// A fixed mirror error produces control-current feedthrough and DC
		// offset even at zero differential input. Keep it separately switchable
		// so tests can distinguish this from the differential-pair distortion.
		const double imbalance = std::isfinite(_configuration.mirrorImbalance) ?
			std::clamp(_configuration.mirrorImbalance, -0.1, 0.1) : 0.0;
		output += control * imbalance;

		if (_configuration.modelEarlyEffect)
		{
			const double earlyVoltage = std::isfinite(_configuration.earlyVoltage) ?
				std::max(std::abs(_configuration.earlyVoltage), 10.0) : 100.0;
			output *= std::clamp(1.0 + outputNodeVolts / earlyVoltage, 0.5, 1.5);
		}
		return std::clamp(output, -MaximumControlCurrent,
			MaximumControlCurrent);
	}

	void Reset() {}

	static constexpr double MaximumControlCurrent = 2.0e-3;

private:
	Configuration _configuration{};

	double SafeThermalVoltage() const
	{
		return std::isfinite(_configuration.thermalVoltage) ?
			std::clamp(_configuration.thermalVoltage, 0.020, 0.035) : 0.02585;
	}

	double SafeEfficiency() const
	{
		return std::isfinite(_configuration.currentTransferEfficiency) ?
			std::clamp(_configuration.currentTransferEfficiency, 0.1, 1.2) : 0.85;
	}

	double SafeInputOffset() const
	{
		return std::isfinite(_configuration.inputOffsetVolts) ?
			std::clamp(_configuration.inputOffsetVolts, -0.010, 0.010) : 0.0;
	}
};

} // namespace tfdsp
