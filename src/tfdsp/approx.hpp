#pragma once

#include <algorithm>
#include <cmath>
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

	/** High-accuracy bounded tanh approximation for audio transfer functions.
	 * The [7/6] Padé form is effectively exact around zero; clamping its tiny
	 * high-argument overshoot keeps the result physical. Across all finite inputs
	 * its maximum absolute error versus std::tanh is below 1.1e-4.
	 */
	inline double TanhPade76(double value)
	{
		if (std::isnan(value))
			return 0.0;
		if (value <= -5.0)
			return -1.0;
		if (value >= 5.0)
			return 1.0;
		const double square = value * value;
		const double numerator = value * (135135.0 + square *
			(17325.0 + square * (378.0 + square)));
		const double denominator = 135135.0 + square *
			(62370.0 + square * (3150.0 + 28.0 * square));
		return std::clamp(numerator / denominator, -1.0, 1.0);
	}

	/** Accurate tan(x) on the ARP filter's prewarp domain.
	 *
	 * The 4072 runs at 2x or 4x and limits its host-domain cutoff to 0.45 of
	 * Nyquist, so the prewarp argument is bounded by 0.225*pi. The degree-21
	 * expansion avoids a general libm call while remaining below 2e-8 relative
	 * error on that interval.
	 */
	inline double TanPrewarp(double value)
	{
		const double square = value * value;
		return value * (1.0 + square * (1.0 / 3.0 + square *
			(2.0 / 15.0 + square * (17.0 / 315.0 + square *
				(62.0 / 2835.0 + square * (1382.0 / 155925.0 + square *
					(21844.0 / 6081075.0 + square *
						(929569.0 / 638512875.0 + square *
							(6404582.0 / 10854718875.0 + square *
								(443861162.0 / 1856156927625.0 + square *
									(18888466084.0 / 194896477400625.0)))))))))));
	}

	/** sin(2*pi*phase) with phase reduction to [-pi/2, pi/2].
	 * The degree-13 odd expansion stays below 7e-10 absolute error and is used
	 * for slow control oscillators where a general argument-reducing libm call
	 * would dominate the actual state update.
	 */
	inline double SinTwoPi(double phase)
	{
		phase -= std::floor(phase);
		if (phase > 0.5)
			phase -= 1.0;
		if (phase > 0.25)
			phase = 0.5 - phase;
		else if (phase < -0.25)
			phase = -0.5 - phase;
		const double angle = 6.283185307179586476925286766559 * phase;
		const double square = angle * angle;
		return angle * (1.0 + square * (-1.0 / 6.0 + square *
			(1.0 / 120.0 + square * (-1.0 / 5040.0 + square *
				(1.0 / 362880.0 + square * (-1.0 / 39916800.0 +
					square / 6227020800.0))))));
	}

	/** Approximate x^-1.3 over the pickup model's physical radial domain.
	 * A degree-six Chebyshev fit on the normalized mantissa keeps relative error
	 * below 5.6e-5 for 0.5 <= x < 2. Values outside that calibrated interval use
	 * the standard library rather than extrapolating the polynomial.
	 */
	inline double PowNegative1p3(double value)
	{
		if (!(value >= 0.5 && value < 2.0))
			return std::pow(value, -1.3);
		constexpr double TwoNegative1p3 = 0.40612619817811774;
		constexpr double coefficients[]{
			1.5875060990310508,
			-0.7040845424791434,
			0.13865337134021724,
			-0.02613043565624296,
			0.004816829894491764,
			-0.0008607565352371927,
			0.00015508195847975295};
		const bool upperOctave = value >= 1.0;
		const double mantissa = upperOctave ? 0.5 * value : value;
		const double coordinate = 4.0 * mantissa - 3.0;
		double next = 0.0;
		double following = 0.0;
		for (int index = 6; index >= 1; --index)
		{
			const double current = 2.0 * coordinate * next - following +
				coefficients[index];
			following = next;
			next = current;
		}
		const double result = coordinate * next - following + coefficients[0];
		return upperOctave ? TwoNegative1p3 * result : result;
	}
}
