#pragma once

#include <cstdint>
#include <cstring>

namespace tfdsp
{
	/** Rack-compatible 2^x approximation with at most 6e-6 relative error.
	 * The polynomial matches Rack's exp2_taylor5, whose coefficients are
	 * credited there to Andy Simper.
	 * The caller must provide a finite exponent in the normal float range.
	 */
	inline float Exp2Taylor5(float value)
	{
		const float biased = value + 127.0f;
		const std::int32_t exponent = static_cast<std::int32_t>(biased);
		const float fraction = biased - static_cast<float>(exponent);
		const std::uint32_t integerBits = static_cast<std::uint32_t>(exponent) << 23;
		float integerPart;
		std::memcpy(&integerPart, &integerBits, sizeof(integerPart));
		const float polynomial = 1.0f + fraction * (0.69315169353961f + fraction *
			(0.2401595990753f + fraction * (0.055817908652f + fraction *
			(0.008991698010f + fraction * 0.001879100722f))));
		return integerPart * polynomial;
	}
}
