#pragma once

#include "tfdsp/finite_audio.hpp"

#include <cmath>
#include <cstddef>

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
    targetCoefficients_ = coefficients_;
    coefficientStep_ = {};
    transitionSamplesRemaining_ = 0;
  }

  void SetTargetCoefficients(const BiquadCoefficients &coefficients,
                             const std::size_t transitionSamples) noexcept {
    const auto safe = IsStable(coefficients) ? coefficients
                                             : BiquadCoefficients{};
    if (transitionSamples == 0) {
      SetCoefficients(safe);
      return;
    }
    targetCoefficients_ = safe;
    const float inverseSamples = 1.f / static_cast<float>(transitionSamples);
    coefficientStep_ = Difference(targetCoefficients_, coefficients_,
                                  inverseSamples);
    transitionSamplesRemaining_ = transitionSamples;
  }

  float Process(float input) noexcept {
    AdvanceCoefficients();
    input = tfdsp::FiniteNormalOrZero(input);
    const float output = coefficients_.b0 * input + state1_;
    state1_ = coefficients_.b1 * input - coefficients_.a1 * output + state2_;
    state2_ = coefficients_.b2 * input - coefficients_.a2 * output;
    if (!std::isfinite(output)) {
      Reset();
      return 0.f;
    }
    state1_ = tfdsp::FiniteNormalOrZero(state1_);
    state2_ = tfdsp::FiniteNormalOrZero(state2_);
    return tfdsp::FiniteNormalOrZero(output);
  }

private:
  static BiquadCoefficients Difference(const BiquadCoefficients &target,
                                       const BiquadCoefficients &current,
                                       const float scale) noexcept {
    return {(target.b0 - current.b0) * scale,
            (target.b1 - current.b1) * scale,
            (target.b2 - current.b2) * scale,
            (target.a1 - current.a1) * scale,
            (target.a2 - current.a2) * scale};
  }

  void AdvanceCoefficients() noexcept {
    if (transitionSamplesRemaining_ == 0)
      return;
    coefficients_.b0 += coefficientStep_.b0;
    coefficients_.b1 += coefficientStep_.b1;
    coefficients_.b2 += coefficientStep_.b2;
    coefficients_.a1 += coefficientStep_.a1;
    coefficients_.a2 += coefficientStep_.a2;
    if (--transitionSamplesRemaining_ == 0)
      coefficients_ = targetCoefficients_;
  }

  static bool IsStable(const BiquadCoefficients &value) noexcept {
    const bool finite = std::isfinite(value.b0) && std::isfinite(value.b1) &&
        std::isfinite(value.b2) && std::isfinite(value.a1) &&
        std::isfinite(value.a2);
    return finite && std::abs(value.a2) < 1.f &&
        1.f + value.a1 + value.a2 > 0.f &&
        1.f - value.a1 + value.a2 > 0.f;
  }

  BiquadCoefficients coefficients_{};
  BiquadCoefficients targetCoefficients_{};
  BiquadCoefficients coefficientStep_{};
  float state1_{};
  float state2_{};
  std::size_t transitionSamplesRemaining_{};
};

} // namespace tfdsp::percussion
