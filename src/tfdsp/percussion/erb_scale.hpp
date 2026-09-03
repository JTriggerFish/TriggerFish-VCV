#pragma once

#include <cmath>

namespace tfdsp::percussion {

inline float ErbRate(const float frequencyHz) noexcept {
  return 21.4f * std::log10(1.f + .00437f * frequencyHz);
}

inline float InverseErbRate(const float erbRate) noexcept {
  return (std::pow(10.f, erbRate / 21.4f) - 1.f) / .00437f;
}

} // namespace tfdsp::percussion
