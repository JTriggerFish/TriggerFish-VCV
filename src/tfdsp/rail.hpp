#pragma once

#include <cmath>

namespace tfdsp
{

// Converts an internal circuit-domain output to a Rack cable voltage without
// a hard boundary. The oversampled stage reserves headroom for IIR decimator
// overshoot. The post-decimation stage is normally an identity operation and
// guarantees that the final cable voltage remains below the Rack/Eurorack
// protected-rail convention. A model whose own oversampled circuit compliance
// already stays below 11.5 V can skip ProcessOversampled, but every Rack-facing
// path should use ProcessPostDecimation after its final processing stage.
class RackOutputAdapter
{
public:
	static double ProcessOversampled(double voltage)
	{
		return SoftRail(voltage, OversampledKneeVolts,
			OversampledLimitVolts);
	}

	static double ProcessPostDecimation(double voltage)
	{
		return SoftRail(voltage, OversampledLimitVolts, CableLimitVolts);
	}

	// Rack defines ±10 V as full scale. The half-volt margin keeps normal full-
	// scale signals linear, while a one-volt transition gives overload enough
	// curvature to avoid a bright, near-hard boundary.
	static constexpr double OversampledKneeVolts = 10.5;
	static constexpr double OversampledLimitVolts = 11.5;
	static constexpr double CableLimitVolts = 11.7;

private:
	static double SoftRail(double voltage, double knee, double limit)
	{
		if (!std::isfinite(voltage))
			return 0.0;
		const double magnitude = std::abs(voltage);
		if (magnitude <= knee)
			return voltage;
		const double headroom = limit - knee;
		const double excess = (magnitude - knee) / headroom;
		const double curved = knee + headroom * excess /
			std::hypot(1.0, excess);
		return std::copysign(curved, voltage);
	}
};

} // namespace tfdsp
