#pragma once

#include <cmath>

namespace tfdsp
{

// TB-303 transistor/output-buffer voltage compliance, calibrated in Rack
// volts. It is exactly linear below the knee, has unity slope at the join, and
// approaches its circuit rail smoothly. Generic Rack cable headroom is handled
// separately by RackOutputAdapter.
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

	static constexpr double KneeVolts = 8.0;
	static constexpr double RailVolts = 11.0;

private:
	static constexpr double HeadroomVolts = RailVolts - KneeVolts;
};

} // namespace tfdsp
