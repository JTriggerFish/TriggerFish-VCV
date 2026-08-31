#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tfdsp::percussion {

// Complementary one-pole shelf tilt. Zero dB is an exact wire, negative tilt
// darkens the source, and positive tilt emphasizes its upper band.
class SpectralTiltFilter {
public:
  void Prepare(const float sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("tilt-filter sample rate must be positive");
    sampleRate_ = sampleRate;
    SetTilt(0.f, 3000.f);
    Reset();
  }

  void Reset() noexcept { lowState_ = 0.f; }

  void SetTilt(float tiltDb, float pivotHz) noexcept {
    tiltDb = std::clamp(std::isfinite(tiltDb) ? tiltDb : 0.f, -24.f, 24.f);
    pivotHz = std::clamp(std::isfinite(pivotHz) ? pivotHz : 3000.f,
                         20.f, .45f * sampleRate_);
    coefficient_ = 1.f - std::exp(-6.283185307179586f * pivotHz / sampleRate_);
    lowGain_ = std::pow(10.f, -tiltDb / 40.f);
    highGain_ = std::pow(10.f, tiltDb / 40.f);
  }

  float Process(const float input) noexcept {
    const float safeInput = tfdsp::FiniteNormalOrZero(input);
    lowState_ += coefficient_ * (safeInput - lowState_);
    lowState_ = tfdsp::FiniteNormalOrZero(lowState_);
    const float high = safeInput - lowState_;
    return tfdsp::FiniteNormalOrZero(
        lowGain_ * lowState_ + highGain_ * high);
  }

private:
  float sampleRate_{48000.f};
  float coefficient_{};
  float lowState_{};
  float lowGain_{1.f};
  float highGain_{1.f};
};

} // namespace tfdsp::percussion
