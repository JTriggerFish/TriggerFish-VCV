#pragma once

#include "static_fractional_delay.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {

// Fractional-delay Schroeder allpass. The configured delay is static and its
// Thiran interpolator keeps the ideal allpass magnitude in recursive use.
class SchroederAllpass {
public:
  void Prepare(const float maximumDelaySamples, const float delaySamples,
               const float feedbackGain) {
    delay_.Prepare(maximumDelaySamples, delaySamples);
    SetFeedbackGain(feedbackGain);
  }

  void Reset() noexcept { delay_.Reset(); }

  void SetFeedbackGain(float feedbackGain) noexcept {
    if (!std::isfinite(feedbackGain))
      feedbackGain = 0.f;
    feedbackGain_ = std::clamp(feedbackGain, -0.999f, 0.999f);
  }

  float Process(const float input) noexcept {
    const float delayed = delay_.Read();
    const float output = delayed - feedbackGain_ * input;
    delay_.Push(input + feedbackGain_ * output);
    return std::isfinite(output) ? output : 0.f;
  }

  float DelaySamples() const noexcept { return delay_.DelaySamples(); }

private:
  StaticFractionalDelay delay_{};
  float feedbackGain_{};
};

} // namespace tfdsp::percussion
