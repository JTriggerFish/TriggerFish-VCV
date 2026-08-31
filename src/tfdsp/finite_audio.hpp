#pragma once

#include <cmath>
#include <limits>

namespace tfdsp {

inline float FiniteNormalOrZero(const float value) noexcept {
  return std::isfinite(value) &&
                 std::abs(value) >= std::numeric_limits<float>::min()
             ? value
             : 0.f;
}

inline double FiniteNormalOrZero(const double value) noexcept {
  return std::isfinite(value) &&
                 std::abs(value) >= std::numeric_limits<double>::min()
             ? value
             : 0.0;
}

} // namespace tfdsp
