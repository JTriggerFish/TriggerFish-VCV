#pragma once

#include "finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace tfdsp {

// Structure-of-arrays decay filter for a fixed FDN frame. Stable coefficients
// are cached per line and independent state recurrence remains SIMD-friendly.
template <std::size_t LineCount> class MultibandDecayFilterBank {
public:
  using Frame = std::array<float, LineCount>;

  void Prepare(const double sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument("decay-filter sample rate must be positive");
    constexpr float Pi = 3.14159265358979323846f;
    constexpr float LowCrossoverHz = 220.f;
    const float rate = static_cast<float>(sampleRate);
    const float highCrossoverHz = std::min(3'500.f, .20f * rate);
    lowAlpha_ = 1.f - std::exp(-2.f * Pi * LowCrossoverHz / rate);
    highAlpha_ = 1.f - std::exp(-2.f * Pi * highCrossoverHz / rate);
    Reset();
  }

  void Reset() noexcept {
    lowState_.fill(0.f);
    belowHighState_.fill(0.f);
    lastPathSeconds_.fill(-1.f);
    lowGain_.fill(1.f);
    midGain_.fill(1.f);
    highGain_.fill(1.f);
    lastLowT60_ = lastMidT60_ = lastHighT60_ = 0.f;
  }

  Frame Process(const Frame &input, const Frame &pathSeconds,
                const float lowT60, const float midT60,
                const float highT60) noexcept {
    const bool t60Changed = lowT60 != lastLowT60_ || midT60 != lastMidT60_ ||
                            highT60 != lastHighT60_;
    Frame output{};
    for (std::size_t line = 0; line < LineCount; ++line) {
      if (t60Changed || pathSeconds[line] != lastPathSeconds_[line])
        UpdateGain(line, pathSeconds[line], lowT60, midT60, highT60);
      const float safeInput = FiniteNormalOrZero(input[line]);
      lowState_[line] += lowAlpha_ * (safeInput - lowState_[line]);
      belowHighState_[line] +=
          highAlpha_ * (safeInput - belowHighState_[line]);
      lowState_[line] = FiniteNormalOrZero(lowState_[line]);
      belowHighState_[line] = FiniteNormalOrZero(belowHighState_[line]);
      const float low = lowState_[line];
      const float mid = belowHighState_[line] - low;
      const float high = safeInput - belowHighState_[line];
      output[line] = FiniteNormalOrZero(lowGain_[line] * low +
                                        midGain_[line] * mid +
                                        highGain_[line] * high);
    }
    lastLowT60_ = lowT60;
    lastMidT60_ = midT60;
    lastHighT60_ = highT60;
    return output;
  }

private:
  void UpdateGain(const std::size_t line, const float pathSeconds,
                  const float lowT60, const float midT60,
                  const float highT60) noexcept {
    lastPathSeconds_[line] = pathSeconds;
    lowGain_[line] = GainForT60(pathSeconds, lowT60);
    midGain_[line] = GainForT60(pathSeconds, midT60);
    highGain_[line] = GainForT60(pathSeconds, highT60);
  }

  static float GainForT60(const float pathSeconds,
                          const float t60) noexcept {
    if (t60 == std::numeric_limits<float>::infinity())
      return 1.f;
    if (!std::isfinite(pathSeconds) || !std::isfinite(t60))
      return 0.f;
    return std::pow(10.f, -3.f * std::max(pathSeconds, 0.f) /
                              std::max(t60, 1.e-4f));
  }

  Frame lowState_{};
  Frame belowHighState_{};
  Frame lastPathSeconds_{};
  Frame lowGain_{};
  Frame midGain_{};
  Frame highGain_{};
  float lowAlpha_{};
  float highAlpha_{};
  float lastLowT60_{};
  float lastMidT60_{};
  float lastHighT60_{};
};

} // namespace tfdsp
