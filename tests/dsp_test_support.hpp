#pragma once

#include <cmath>
#include <iostream>
#include <string_view>

// Small shared assertions for standalone numerical tests.
namespace dsp_test {
inline int failures{};

inline void Check(bool condition, std::string_view message) {
  if (condition)
    return;
  std::cerr << "FAIL: " << message << '\n';
  ++failures;
}

inline void CheckNear(double actual, double expected, double tolerance,
                      std::string_view message) {
  if (std::abs(actual - expected) <= tolerance)
    return;
  std::cerr << "FAIL: " << message << " (actual " << actual << ", expected "
            << expected << " +/- " << tolerance << ")\n";
  ++failures;
}
} // namespace dsp_test
