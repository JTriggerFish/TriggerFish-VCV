#pragma once

#include "static_fractional_delay.hpp"
#include "tfdsp/finite_audio.hpp"

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
    const float safeInput = tfdsp::FiniteNormalOrZero(input);
    const float delayed = delay_.Read();
    const float output = tfdsp::FiniteNormalOrZero(
        delayed - feedbackGain_ * safeInput);
    delay_.Push(safeInput + feedbackGain_ * output);
    return output;
  }

  float DelaySamples() const noexcept { return delay_.DelaySamples(); }

  // A non-zero Schroeder coefficient creates an immediate feed-through term.
  // The configured delay remains its nominal group-delay scale, not latency.
  float MinimumPropagationSamples() const noexcept {
    return feedbackGain_ == 0.f ? delay_.DelaySamples() : 0.f;
  }

private:
  StaticFractionalDelay delay_{};
  float feedbackGain_{};
};

} // namespace tfdsp::percussion
