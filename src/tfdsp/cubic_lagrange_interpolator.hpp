#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace tfdsp {

// Centred four-point Lagrange interpolation over a causal delay-line reader.
// The callback receives an integer distance behind the current write position.
template <typename IntegerReader>
float ReadCubicLagrange(float delaySamples, const std::size_t capacity,
                        IntegerReader readInteger) noexcept {
  if (capacity < 5)
    return 0.f;
  if (!std::isfinite(delaySamples))
    delaySamples = 2.f;
  delaySamples = std::clamp(delaySamples, 2.f,
                            static_cast<float>(capacity - 3));
  const auto integer = static_cast<std::size_t>(std::floor(delaySamples));
  const float mu = delaySamples - static_cast<float>(integer);
  const std::array<float, 4> coefficients{
      -mu * (mu - 1.f) * (mu - 2.f) / 6.f,
      (mu + 1.f) * (mu - 1.f) * (mu - 2.f) / 2.f,
      -(mu + 1.f) * mu * (mu - 2.f) / 2.f,
      (mu + 1.f) * mu * (mu - 1.f) / 6.f};
  constexpr std::array<int, 4> offsets{-1, 0, 1, 2};
  float value = 0.f;
  for (std::size_t tap = 0; tap < coefficients.size(); ++tap)
    value += coefficients[tap] * readInteger(static_cast<std::size_t>(
        static_cast<int>(integer) + offsets[tap]));
  return value;
}

} // namespace tfdsp
