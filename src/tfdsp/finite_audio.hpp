#pragma once

#include <cstdint>
#include <cstring>
#include <limits>

namespace tfdsp {

static_assert(sizeof(float) == sizeof(std::uint32_t) &&
                  std::numeric_limits<float>::is_iec559,
              "audio float classification requires IEEE-754 binary32");
static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                  std::numeric_limits<double>::is_iec559,
              "audio double classification requires IEEE-754 binary64");

inline float FiniteNormalOrZero(const float value) noexcept {
  std::uint32_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  const std::uint32_t magnitude = bits & 0x7fffffffu;
  return magnitude >= 0x00800000u && magnitude < 0x7f800000u ? value : 0.f;
}

inline double FiniteNormalOrZero(const double value) noexcept {
  std::uint64_t bits{};
  std::memcpy(&bits, &value, sizeof(bits));
  const std::uint64_t magnitude = bits & 0x7fffffffffffffffull;
  return magnitude >= 0x0010000000000000ull &&
                 magnitude < 0x7ff0000000000000ull
             ? value
             : 0.0;
}

} // namespace tfdsp
