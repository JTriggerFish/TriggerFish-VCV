#pragma once

#include "modulated_fractional_delay.hpp"
#include "tfdsp/smooth_random_modulator.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace tfdsp::percussion {

// Signal-independent, band-limited delay motion for decorrelating successive
// dispersion-loop circulations without the periodic signature of an LFO.
class SlowModulatedDelay {
public:
  void Prepare(const float sampleRate, const float maximumDelaySamples,
               const std::uint32_t seed) {
    delay_.Prepare(maximumDelaySamples);
    modulator_.Prepare(sampleRate, rateHz_, seed);
    maximumDelaySamples_ = maximumDelaySamples;
    Reset();
  }

  void Reset() noexcept {
    delay_.Reset();
    modulator_.Reset();
  }

  // Static configuration resets the random trajectory. Live depth automation
  // will use a separate smoothed control and must not repeatedly call this.
  void SetStaticParameters(float centreDelaySamples, float depthSamples,
                           float rateHz, const float sampleRate,
                           const std::uint32_t seed) noexcept {
    centreDelaySamples_ = std::clamp(
        std::isfinite(centreDelaySamples) ? centreDelaySamples :
            ModulatedFractionalDelay::MinimumDelaySamples,
        ModulatedFractionalDelay::MinimumDelaySamples, maximumDelaySamples_);
    const float margin = std::min(
        centreDelaySamples_ - ModulatedFractionalDelay::MinimumDelaySamples,
        maximumDelaySamples_ - centreDelaySamples_);
    depthSamples_ = std::clamp(std::isfinite(depthSamples) ? depthSamples : 0.f,
                               0.f, std::max(0.f, margin));
    rateHz_ = std::clamp(std::isfinite(rateHz) ? rateHz : .1f, .01f, 20.f);
    modulator_.Prepare(sampleRate, rateHz_, seed);
  }

  float Process(const float input) noexcept {
    const float tap = centreDelaySamples_ + depthSamples_ * modulator_.Next();
    return delay_.Process(input, tap);
  }

  float CentreDelaySamples() const noexcept { return centreDelaySamples_; }

private:
  ModulatedFractionalDelay delay_{};
  tfdsp::SmoothRandomModulator modulator_{};
  float maximumDelaySamples_{64.f};
  float centreDelaySamples_{16.f};
  float depthSamples_{};
  float rateHz_{.1f};
};

} // namespace tfdsp::percussion
