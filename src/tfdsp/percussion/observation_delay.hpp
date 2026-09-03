#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace tfdsp::percussion {

// Static magnitude-preserving observation delay. Unlike a feedback delay it
// supports exact zero delay, which is required for a direct microphone path.
class ObservationDelay {
public:
  // Production observation paths fit inline (10 ms at the supported 384 kHz
  // ceiling). Longer experimental paths retain an explicit dynamic fallback.
  static constexpr std::size_t InlineDelaySamples = 3840;

  void Prepare(const float maximumDelaySamples, const float delaySamples) {
    if (!std::isfinite(maximumDelaySamples) || maximumDelaySamples < 1.f)
      throw std::invalid_argument("observation delay capacity is too small");
    maximumDelaySamples_ = maximumDelaySamples;
    bufferSize_ = static_cast<std::size_t>(
        std::ceil(maximumDelaySamples)) + 1;
    if (bufferSize_ > inlineBuffer_.size())
      extendedBuffer_.assign(bufferSize_, 0.f);
    else
      extendedBuffer_.clear();
    SetStaticDelaySamples(delaySamples);
  }

  void Reset() noexcept {
    if (UsesExtendedBuffer())
      std::fill(extendedBuffer_.begin(), extendedBuffer_.end(), 0.f);
    else
      inlineBuffer_.fill(0.f);
    writeIndex_ = 0;
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
                              : ReadInteger(integerDelaySamples_);
    float output = delayed;
    if (!bypassAllpass_) {
      output = allpassCoefficient_ * delayed + previousInput_ -
               allpassCoefficient_ * previousOutput_;
      previousInput_ = delayed;
      previousOutput_ = tfdsp::FiniteNormalOrZero(output);
    }
    Sample(writeIndex_) = input;
    if (++writeIndex_ == bufferSize_) writeIndex_ = 0;
    return tfdsp::FiniteNormalOrZero(output);
  }

  float DelaySamples() const noexcept { return delaySamples_; }

private:
  float ReadInteger(const std::size_t delaySamples) const noexcept {
    const auto index = writeIndex_ >= delaySamples
                           ? writeIndex_ - delaySamples
                           : writeIndex_ + bufferSize_ - delaySamples;
    return Sample(index);
  }

  bool UsesExtendedBuffer() const noexcept {
    return bufferSize_ > inlineBuffer_.size();
  }

  float &Sample(const std::size_t index) noexcept {
    return UsesExtendedBuffer() ? extendedBuffer_[index]
                                : inlineBuffer_[index];
  }

  float Sample(const std::size_t index) const noexcept {
    return UsesExtendedBuffer() ? extendedBuffer_[index]
                                : inlineBuffer_[index];
  }

  std::array<float, InlineDelaySamples + 1> inlineBuffer_{};
  std::vector<float> extendedBuffer_{};
  std::size_t bufferSize_{2};
  std::size_t writeIndex_{};
  float maximumDelaySamples_{1.f};
  float delaySamples_{};
  float allpassCoefficient_{};
  float previousInput_{};
  float previousOutput_{};
  std::size_t integerDelaySamples_{};
  bool bypassAllpass_{true};
};

} // namespace tfdsp::percussion
