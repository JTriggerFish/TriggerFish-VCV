#pragma once

#include "fractional_delay_line.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace tfdsp::percussion {

// Static magnitude-preserving observation delay. Unlike a feedback delay it
// supports exact zero delay, which is required for a direct microphone path.
class ObservationDelay {
public:
  void Prepare(const float maximumDelaySamples, const float delaySamples) {
    if (!std::isfinite(maximumDelaySamples) || maximumDelaySamples < 1.f)
      throw std::invalid_argument("observation delay capacity is too small");
    maximumDelaySamples_ = maximumDelaySamples;
    line_.Prepare(std::max(maximumDelaySamples,
                           FractionalDelayLine::MinimumSincDelaySamples));
    SetStaticDelaySamples(delaySamples);
  }

  void Reset() noexcept {
    line_.Reset();
    previousInput_ = previousOutput_ = 0.f;
  }

  void SetStaticDelaySamples(float delaySamples) {
    if (!std::isfinite(delaySamples))
      throw std::invalid_argument("observation delay must be finite");
    delaySamples_ = std::clamp(delaySamples, 0.f, maximumDelaySamples_);
    integerDelaySamples_ = static_cast<std::size_t>(std::floor(delaySamples_));
    const float fraction = delaySamples_ -
                           static_cast<float>(integerDelaySamples_);
    bypassAllpass_ = fraction < 1.e-6f;
    allpassCoefficient_ = bypassAllpass_ ? 0.f :
        (1.f - fraction) / (1.f + fraction);
    Reset();
  }

  float Process(float input) noexcept {
    input = tfdsp::FiniteNormalOrZero(input);
    const float delayed = integerDelaySamples_ == 0
                              ? input
                              : line_.ReadInteger(integerDelaySamples_);
    float output = delayed;
    if (!bypassAllpass_) {
      output = allpassCoefficient_ * delayed + previousInput_ -
               allpassCoefficient_ * previousOutput_;
      previousInput_ = delayed;
      previousOutput_ = tfdsp::FiniteNormalOrZero(output);
    }
    line_.Push(input);
    return tfdsp::FiniteNormalOrZero(output);
  }

  float DelaySamples() const noexcept { return delaySamples_; }

private:
  FractionalDelayLine line_{};
  float maximumDelaySamples_{1.f};
  float delaySamples_{};
  float allpassCoefficient_{};
  float previousInput_{};
  float previousOutput_{};
  std::size_t integerDelaySamples_{};
  bool bypassAllpass_{true};
};

} // namespace tfdsp::percussion
