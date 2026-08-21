#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tfdsp {

// A complementary three-band attenuation filter for an FDN delay segment.
// When all three gains are one, low + mid + high reconstructs the input
// exactly. Different gain-per-traversal values therefore set band-dependent
// T60 without adding an unrelated tone filter to the loop.
class MultibandDecayFilter {
public:
  void Prepare(const double sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument("decay-filter sample rate must be positive");
    sampleRate_ = static_cast<float>(sampleRate);
    constexpr float Pi = 3.14159265358979323846f;
    constexpr float LowCrossoverHz = 220.f;
    const float highCrossoverHz = std::min(3'500.f, 0.20f * sampleRate_);
    lowAlpha_ =
        1.f - std::exp(-2.f * Pi * LowCrossoverHz / sampleRate_);
    highAlpha_ =
        1.f - std::exp(-2.f * Pi * highCrossoverHz / sampleRate_);
    Reset();
  }

  void Reset() noexcept {
    lowState_ = 0.f;
    belowHighState_ = 0.f;
    lastPathSeconds_ = -1.f;
  }

  float Process(const float input, const float pathSeconds,
                const float lowT60, const float midT60,
                const float highT60) noexcept {
    UpdateGains(pathSeconds, lowT60, midT60, highT60);
    lowState_ += lowAlpha_ * (input - lowState_);
    belowHighState_ += highAlpha_ * (input - belowHighState_);
    const float low = lowState_;
    const float mid = belowHighState_ - lowState_;
    const float high = input - belowHighState_;
    return lowGain_ * low + midGain_ * mid + highGain_ * high;
  }

private:
  void UpdateGains(const float pathSeconds, const float lowT60,
                   const float midT60, const float highT60) noexcept {
    if (pathSeconds == lastPathSeconds_ && lowT60 == lastLowT60_ &&
        midT60 == lastMidT60_ && highT60 == lastHighT60_)
      return;
    lastPathSeconds_ = pathSeconds;
    lastLowT60_ = lowT60;
    lastMidT60_ = midT60;
    lastHighT60_ = highT60;
    lowGain_ = GainForT60(pathSeconds, lowT60);
    midGain_ = GainForT60(pathSeconds, midT60);
    highGain_ = GainForT60(pathSeconds, highT60);
  }

  static float GainForT60(const float pathSeconds,
                          const float t60) noexcept {
    if (!std::isfinite(t60))
      return 1.f;
    return std::pow(10.f, -3.f * std::max(pathSeconds, 0.f) /
                              std::max(t60, 1.e-4f));
  }

  float sampleRate_{48'000.f};
  float lowAlpha_{};
  float highAlpha_{};
  float lowState_{};
  float belowHighState_{};
  float lastPathSeconds_{-1.f};
  float lastLowT60_{};
  float lastMidT60_{};
  float lastHighT60_{};
  float lowGain_{1.f};
  float midGain_{1.f};
  float highGain_{1.f};
};

} // namespace tfdsp
