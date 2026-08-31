#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace tfdsp::percussion {

struct ThreeBandDecayTimes {
  float lowSeconds{2.f};
  float middleSeconds{1.f};
  float highSeconds{.5f};
};

// Complementary low/middle/high loop-loss filter. Equal gains reconstruct the
// input exactly, while T60 parameters map directly to per-circulation gains.
class ThreeBandDecayFilter {
public:
  void Prepare(const float sampleRate, const float lowCrossoverHz,
               const float highCrossoverHz) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("decay-filter sample rate must be positive");
    sampleRate_ = sampleRate;
    lowCoefficient_ = Pole(lowCrossoverHz);
    highCoefficient_ = Pole(std::max(lowCrossoverHz, highCrossoverHz));
    Reset();
  }

  void Reset() noexcept { lowState_ = belowHighState_ = 0.f; }

  void SetDecayTimes(const float pathSeconds,
                     const ThreeBandDecayTimes &times) noexcept {
    lowGain_ = GainForT60(pathSeconds, times.lowSeconds);
    middleGain_ = GainForT60(pathSeconds, times.middleSeconds);
    highGain_ = GainForT60(pathSeconds, times.highSeconds);
  }

  float Process(float input) noexcept {
    if (!std::isfinite(input)) {
      Reset();
      return 0.f;
    }
    lowState_ += lowCoefficient_ * (input - lowState_);
    belowHighState_ += highCoefficient_ * (input - belowHighState_);
    lowState_ = tfdsp::FiniteNormalOrZero(lowState_);
    belowHighState_ = tfdsp::FiniteNormalOrZero(belowHighState_);
    const float low = lowState_;
    const float middle = belowHighState_ - lowState_;
    const float high = input - belowHighState_;
    const float output =
        lowGain_ * low + middleGain_ * middle + highGain_ * high;
    if (!std::isfinite(output)) {
      Reset();
      return 0.f;
    }
    return tfdsp::FiniteNormalOrZero(output);
  }

  static float GainForT60(const float pathSeconds,
                          const float t60Seconds) noexcept {
    if (!std::isfinite(pathSeconds) || std::isnan(t60Seconds))
      return 0.f;
    if (t60Seconds == std::numeric_limits<float>::infinity())
      return 1.f;
    return std::pow(10.f, -3.f * std::max(0.f, pathSeconds) /
                              std::max(1.e-4f, t60Seconds));
  }

private:
  float Pole(float frequencyHz) const noexcept {
    frequencyHz = std::clamp(std::isfinite(frequencyHz) ? frequencyHz : 0.f,
                             0.f, .49f * sampleRate_);
    return 1.f - std::exp(-6.283185307179586f * frequencyHz / sampleRate_);
  }

  float sampleRate_{48000.f};
  float lowCoefficient_{};
  float highCoefficient_{};
  float lowState_{};
  float belowHighState_{};
  float lowGain_{1.f};
  float middleGain_{1.f};
  float highGain_{1.f};
};

} // namespace tfdsp::percussion
