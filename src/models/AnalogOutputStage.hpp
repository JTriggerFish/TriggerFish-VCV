#pragma once

#include <cmath>

namespace tfdsp
{

// Rack-facing approximation of transistor/output-buffer voltage compliance.
// It is exactly linear below the knee, has unity slope at the join, and
// approaches the rail smoothly. Processing it before decimation keeps the
// harmonics produced by overload inside the oversampled signal path.
class AnalogOutputStage
{
public:
	static double Process(double voltage)
	{
		if (!std::isfinite(voltage))
			return 0.0;
		const double magnitude = std::abs(voltage);
		if (magnitude <= KneeVolts)
			return voltage;
		const double curved = KneeVolts + HeadroomVolts * std::tanh(
			(magnitude - KneeVolts) / HeadroomVolts);
		return std::copysign(curved, voltage);
	}

	// Catch interpolation/decimation overshoot beyond the modeled rail without
	// introducing a second hard boundary. This is linear through the 11 V join
	// and approaches the 12 V numerical limit algebraically.
	static double ProcessSafety(double voltage)
	{
		if (!std::isfinite(voltage))
			return 0.0;
		const double magnitude = std::abs(voltage);
		if (magnitude <= RailVolts)
			return voltage;
		const double excess = (magnitude - RailVolts) / SafetyHeadroomVolts;
		const double curved = RailVolts + SafetyHeadroomVolts * excess /
			std::sqrt(1.0 + excess * excess);
		return std::copysign(curved, voltage);
	}

	static constexpr double KneeVolts = 8.0;
	static constexpr double RailVolts = 11.0;
	static constexpr double SafetyLimitVolts = 12.0;

private:
	static constexpr double HeadroomVolts = RailVolts - KneeVolts;
	static constexpr double SafetyHeadroomVolts = SafetyLimitVolts - RailVolts;
};

} // namespace tfdsp
