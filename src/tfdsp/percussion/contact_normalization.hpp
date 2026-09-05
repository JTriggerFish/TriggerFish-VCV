#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>

namespace tfdsp::percussion::contact_normalization {

// Reciprocal area of sin(pi*(k+1)/(N+1)). Body ports carry impulse
// increments, whereas the source's direct-audio port retains peak units.
inline float HalfSineImpulseScale(const std::size_t count) noexcept {
  constexpr double Pi = 3.14159265358979323846;
  return static_cast<float>(std::tan(Pi / (2.0 * (count + 1.0))));
}

// Coherent chirps use envelope area, not signed waveform area (which can
// vanish). This bounds their impulse while retaining frequency selectivity.
inline float ChirpImpulseScale(const std::size_t count,
                              const float decayNepers) noexcept {
  constexpr double Pi = 3.14159265358979323846;
  const double denominator = static_cast<double>(count - 1);
  const auto ratio = std::exp(std::complex<double>(-decayNepers, Pi) /
                              denominator);
  const double area = std::imag((1.0 - std::pow(ratio, count)) /
                                (1.0 - ratio));
  return static_cast<float>(1.0 / std::max(area, 1.e-12));
}

// Exact sum of squared smoothstep attack samples; O(1), including long
// brush gestures. Summation formulae avoid an event-time rendering pass.
inline double SmoothstepSquaredSum(const std::size_t count) noexcept {
  if (count == 0) return 0.0;
  const double n = static_cast<double>(count);
  const double n2 = n * n;
  const double n4 = n2 * n2;
  const double s4 = n * (n + 1) * (2*n + 1) * (3*n2 + 3*n - 1) / 30;
  const double s5 = n2 * (n + 1) * (n + 1) * (2*n2 + 2*n - 1) / 12;
  const double s6 = n * (n + 1) * (2*n + 1) *
      (3*n4 + 6*n2*n - 3*n + 1) / 42;
  return 9*s4/n4 - 12*s5/(n4*n) + 4*s6/(n4*n2);
}

inline float NoiseImpulseScale(const std::size_t attack,
                              const std::size_t hold,
                              const std::size_t decay) noexcept {
  double energy = SmoothstepSquaredSum(attack) + hold;
  if (decay > 0) {
    const double step = 2.0 * std::log(1.e-4) / decay;
    energy += std::exp(step) * std::expm1(step * decay) / std::expm1(step);
  }
  // Bipolar uniform noise has variance 1/3. Tilt remains an explicit
  // spectral filter; do not adapt this scale to the realized random seed.
  return static_cast<float>(std::sqrt(3.0 / std::max(energy, 1.e-12)));
}

} // namespace tfdsp::percussion::contact_normalization
