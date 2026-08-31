#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tfdsp::percussion {

// Multipliers applied inside a resonator recurrence. The ordered band gains
// describe the usual increasingly strong high-frequency damping of a mute.
struct PassiveConstraintGains {
  float broadband{1.f};
  float low{1.f};
  float middle{1.f};
  float high{1.f};
};

inline PassiveConstraintGains
SanitizePassiveConstraint(PassiveConstraintGains gains) noexcept {
  gains.broadband = std::clamp(tfdsp::FiniteNormalOrZero(gains.broadband), 0.f,
                               1.f);
  gains.low = std::clamp(tfdsp::FiniteNormalOrZero(gains.low), 0.f, 1.f);
  gains.middle = std::clamp(tfdsp::FiniteNormalOrZero(gains.middle), 0.f,
                            gains.low);
  gains.high = std::clamp(tfdsp::FiniteNormalOrZero(gains.high), 0.f,
                          gains.middle);
  return gains;
}

// Asymmetric gain-domain smoothing for a passive constraint. Increasing loss
// follows the attack time; releasing the object follows the release time.
class DynamicLossController {
public:
  void Prepare(const float sampleRate, const float attackSeconds,
               const float releaseSeconds) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("loss-controller sample rate must be positive");
    sampleRate_ = sampleRate;
    SetTimes(attackSeconds, releaseSeconds);
    Reset();
  }

  void SetTimes(const float attackSeconds, const float releaseSeconds) noexcept {
    attackCoefficient_ = Coefficient(attackSeconds);
    releaseCoefficient_ = Coefficient(releaseSeconds);
  }

  void Reset(PassiveConstraintGains gains = {}) noexcept {
    current_ = target_ = SanitizePassiveConstraint(gains);
  }

  void SetTarget(const PassiveConstraintGains gains) noexcept {
    target_ = SanitizePassiveConstraint(gains);
  }

  PassiveConstraintGains Process() noexcept {
    current_.broadband = Follow(current_.broadband, target_.broadband);
    current_.low = Follow(current_.low, target_.low);
    current_.middle = Follow(current_.middle, target_.middle);
    current_.high = Follow(current_.high, target_.high);
    current_ = SanitizePassiveConstraint(current_);
    return current_;
  }

  PassiveConstraintGains Current() const noexcept { return current_; }

private:
  float Coefficient(float seconds) const noexcept {
    seconds = std::clamp(std::isfinite(seconds) ? seconds : 0.f, 0.f, 10.f);
    return seconds == 0.f ? 1.f : -std::expm1(-1.f / (seconds * sampleRate_));
  }

  float Follow(const float current, const float target) const noexcept {
    const float coefficient = target < current ? attackCoefficient_
                                               : releaseCoefficient_;
    return tfdsp::FiniteNormalOrZero(current + coefficient * (target - current));
  }

  PassiveConstraintGains current_{};
  PassiveConstraintGains target_{};
  float sampleRate_{48000.f};
  float attackCoefficient_{1.f};
  float releaseCoefficient_{1.f};
};

} // namespace tfdsp::percussion
