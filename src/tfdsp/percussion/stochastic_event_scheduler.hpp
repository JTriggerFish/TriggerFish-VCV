#pragma once

#include "deterministic_random.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
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
    rateHz_ = 0.f;
    samplesUntilEvent_ = 0.f;
  }

  void SetRateHz(const float rateHz) noexcept {
    const float next = std::clamp(
        std::isfinite(rateHz) ? rateHz : 0.f, 0.f, .45f * sampleRate_);
    if (next > 0.f && rateHz_ > 0.f) {
      samplesUntilEvent_ *= rateHz_ / next;
    } else if (next > 0.f) {
      rateHz_ = next;
      samplesUntilEvent_ = NextIntervalSamples();
      return;
    }
    rateHz_ = next;
    if (rateHz_ == 0.f)
      samplesUntilEvent_ = 0.f;
  }

  std::size_t Process() noexcept {
    if (rateHz_ == 0.f)
      return 0;
    samplesUntilEvent_ -= 1.f;
    std::size_t events = 0;
    while (samplesUntilEvent_ <= 0.f && events < MaximumEventsPerSample) {
      ++events;
      samplesUntilEvent_ += NextIntervalSamples();
    }
    return events;
  }

private:
  static constexpr std::size_t MaximumEventsPerSample = 16;

  float NextIntervalSamples() noexcept {
    const float uniform = std::max(random_.Uniform(), 1.f / 16777216.f);
    return -std::log(uniform) * sampleRate_ / rateHz_;
  }

  DeterministicRandom random_{};
  float sampleRate_{48000.f};
  float rateHz_{};
  float samplesUntilEvent_{};
};

} // namespace tfdsp::percussion
