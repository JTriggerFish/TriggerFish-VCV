#pragma once

#include "fractional_delay_line.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace tfdsp::percussion {

// Static integer delay followed by a first-order Thiran allpass. The allpass
// preserves magnitude in a feedback loop; SetDelaySamples() resets its state
// and is therefore a prepare-time operation rather than an automatable control.
class StaticFractionalDelay {
public:
  static constexpr float MinimumDelaySamples = 1.f;

  void Prepare(const float maximumDelaySamples, const float delaySamples) {
    if (!std::isfinite(maximumDelaySamples) || maximumDelaySamples < 2.f)
      throw std::invalid_argument("static delay capacity is too small");
    line_.Prepare(std::max(maximumDelaySamples,
                           FractionalDelayLine::MinimumSincDelaySamples));
    maximumDelaySamples_ = maximumDelaySamples;
    SetDelaySamples(delaySamples);
  }

  void SetDelaySamples(float delaySamples) {
    if (!std::isfinite(delaySamples))
      throw std::invalid_argument("static delay must be finite");
    delaySamples = std::clamp(delaySamples, MinimumDelaySamples,
                              maximumDelaySamples_);
    integerDelaySamples_ =
        static_cast<std::size_t>(std::floor(delaySamples));
    const float fraction = delaySamples -
                           static_cast<float>(integerDelaySamples_);
    bypassAllpass_ = fraction < 1.e-6f;
    allpassCoefficient_ = bypassAllpass_ ? 0.f : (1.f - fraction) / (1.f + fraction);
    delaySamples_ = delaySamples;
    Reset();
  }

  void Reset() noexcept {
    line_.Reset();
    previousInput_ = previousOutput_ = 0.f;
  }

  float Read() noexcept {
    const float input = line_.ReadInteger(integerDelaySamples_);
    if (bypassAllpass_)
      return input;
    const float output = allpassCoefficient_ * input + previousInput_ -
                         allpassCoefficient_ * previousOutput_;
    previousInput_ = input;
    previousOutput_ = output;
    return std::isfinite(output) ? output : 0.f;
  }

  void Push(const float input) noexcept { line_.Push(input); }

  float Process(const float input) noexcept {
    const float output = Read();
    Push(input);
    return output;
  }

  float DelaySamples() const noexcept { return delaySamples_; }

private:
  FractionalDelayLine line_{};
  std::size_t integerDelaySamples_{1};
  float maximumDelaySamples_{2.f};
  float delaySamples_{1.f};
  float allpassCoefficient_{};
  float previousInput_{};
  float previousOutput_{};
  bool bypassAllpass_{true};
};

} // namespace tfdsp::percussion
