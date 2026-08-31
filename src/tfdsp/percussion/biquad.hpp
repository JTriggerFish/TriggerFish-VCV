#pragma once

#include <cmath>

namespace tfdsp::percussion {

struct BiquadCoefficients {
  float b0{1.f};
  float b1{};
  float b2{};
  float a1{};
  float a2{};
};

// Transposed-direct-form-II biquad. Invalid or unstable coefficient sets fall
// back to identity instead of placing unsafe state on the audio thread.
class Biquad {
public:
  void Reset() noexcept { state1_ = state2_ = 0.f; }

  void SetCoefficients(const BiquadCoefficients &coefficients) noexcept {
    coefficients_ = IsStable(coefficients) ? coefficients : BiquadCoefficients{};
  }

  float Process(float input) noexcept {
    if (!std::isfinite(input))
      input = 0.f;
    const float output = coefficients_.b0 * input + state1_;
    state1_ = coefficients_.b1 * input - coefficients_.a1 * output + state2_;
    state2_ = coefficients_.b2 * input - coefficients_.a2 * output;
    if (!std::isfinite(output)) {
      Reset();
      return 0.f;
    }
    return output;
  }

private:
  static bool IsStable(const BiquadCoefficients &value) noexcept {
    const bool finite = std::isfinite(value.b0) && std::isfinite(value.b1) &&
        std::isfinite(value.b2) && std::isfinite(value.a1) &&
        std::isfinite(value.a2);
    return finite && std::abs(value.a2) < 1.f &&
        1.f + value.a1 + value.a2 > 0.f &&
        1.f - value.a1 + value.a2 > 0.f;
  }

  BiquadCoefficients coefficients_{};
  float state1_{};
  float state2_{};
};

} // namespace tfdsp::percussion
