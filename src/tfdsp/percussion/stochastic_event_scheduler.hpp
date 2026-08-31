#pragma once

#include "deterministic_random.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace tfdsp::percussion {

// Continuous-time Poisson scheduler sampled on the audio grid. Random work is
// performed only when an event is scheduled, not on every quiet sample.
class StochasticEventScheduler {
public:
  void Prepare(const float sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("contact scheduler rate must be positive");
    sampleRate_ = sampleRate;
    Reset(1);
  }

  void Reset(const std::uint32_t seed) noexcept {
    random_.Seed(seed);
    samplesUntilEvent_ = 0.f;
  }

  void SetRateHz(const float rateHz) noexcept {
    const float next = std::clamp(
        std::isfinite(rateHz) ? rateHz : 0.f, 0.f, .45f * sampleRate_);
    if (next > 0.f && rateHz_ > 0.f)
      samplesUntilEvent_ *= rateHz_ / next;
    else if (next > 0.f)
      samplesUntilEvent_ = 0.f;
    rateHz_ = next;
    if (rateHz_ == 0.f)
      samplesUntilEvent_ = 0.f;
  }

  bool Process() noexcept {
    if (rateHz_ == 0.f)
      return false;
    if (samplesUntilEvent_ > 0.f) {
      samplesUntilEvent_ -= 1.f;
      return false;
    }
    ScheduleNext();
    return true;
  }

private:
  void ScheduleNext() noexcept {
    const float uniform = std::max(random_.Uniform(), 1.f / 16777216.f);
    samplesUntilEvent_ = -std::log(uniform) * sampleRate_ / rateHz_;
  }

  DeterministicRandom random_{};
  float sampleRate_{48000.f};
  float rateHz_{};
  float samplesUntilEvent_{};
};

} // namespace tfdsp::percussion
