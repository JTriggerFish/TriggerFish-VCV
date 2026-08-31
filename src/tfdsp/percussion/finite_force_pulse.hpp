#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace tfdsp::percussion {

// Peak-normalized half-sine force pulse. Duration controls compliance while
// amplitude controls peak force; retriggering starts a fresh contact.
class FiniteForcePulse {
public:
  void Prepare(const float sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("force-pulse sample rate must be positive");
    sampleRate_ = sampleRate;
    Reset();
  }

  void Reset() noexcept {
    sample_ = 0;
    sampleCount_ = 0;
    amplitude_ = 0.f;
    sine_ = 0.f;
    cosine_ = 1.f;
  }

  void Trigger(float durationSeconds, float amplitude) noexcept {
    if (!std::isfinite(durationSeconds))
      durationSeconds = 0.f;
    if (!std::isfinite(amplitude))
      amplitude = 0.f;
    sampleCount_ = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::lround(
               std::clamp(durationSeconds, 0.f, 1.f) * sampleRate_)));
    amplitude_ = std::clamp(tfdsp::FiniteNormalOrZero(amplitude), 0.f, 16.f);
    if (amplitude_ == 0.f) {
      sampleCount_ = 0;
      return;
    }
    sample_ = 0;
    constexpr double Pi = 3.1415926535897932384626433832795;
    const double step = Pi / (static_cast<double>(sampleCount_) + 1.0);
    rotationSine_ = static_cast<float>(std::sin(step));
    rotationCosine_ = static_cast<float>(std::cos(step));
    sine_ = rotationSine_;
    cosine_ = rotationCosine_;
  }

  float Process() noexcept {
    if (sample_ >= sampleCount_)
      return 0.f;
    const float output = amplitude_ * sine_;
    const float nextSine = sine_ * rotationCosine_ + cosine_ * rotationSine_;
    cosine_ = cosine_ * rotationCosine_ - sine_ * rotationSine_;
    sine_ = nextSine;
    ++sample_;
    return output;
  }

  bool Active() const noexcept { return sample_ < sampleCount_; }

private:
  float sampleRate_{48000.f};
  float amplitude_{};
  float sine_{};
  float cosine_{1.f};
  float rotationSine_{};
  float rotationCosine_{1.f};
  std::size_t sample_{};
  std::size_t sampleCount_{};
};

} // namespace tfdsp::percussion
