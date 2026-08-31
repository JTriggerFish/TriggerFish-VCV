#pragma once

#include <cmath>
#include <iostream>
#include <string_view>

namespace percussion_test {

inline int failures{};

inline void Check(const bool condition, const std::string_view message) {
  if (condition)
    return;
  std::cerr << "FAIL: " << message << '\n';
  ++failures;
}

inline void CheckNear(const double actual, const double expected,
                      const double tolerance, const std::string_view message) {
  if (std::abs(actual - expected) <= tolerance)
    return;
  std::cerr << "FAIL: " << message << " (actual " << actual << ", expected "
            << expected << " +/- " << tolerance << ")\n";
  ++failures;
}

inline float Sine(const std::size_t sample, const float frequencyHz,
                  const float sampleRate) noexcept {
  constexpr double TwoPi = 6.283185307179586476925286766559;
  return static_cast<float>(std::sin(
      TwoPi * frequencyHz * static_cast<double>(sample) / sampleRate));
}

} // namespace percussion_test
