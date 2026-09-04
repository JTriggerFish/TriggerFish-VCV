#pragma once

#include "passive_constraint.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tfdsp::percussion {

// Per-sample attenuation for modal recurrences. These gains are separate from
// PassiveConstraintGains, whose attenuation is defined per delay traversal.
struct ModalDampingGains {
  float broadband{1.f};
  float low{1.f};
  float middle{1.f};
  float high{1.f};
};

// Converts traversal-domain passive losses into smoothly changing per-sample
// modal damping. Conversion happens only when the target changes.
class ModalConstraintController {
public:
  void Prepare(float sampleRate, float referenceTraversalSeconds,
               float attackSeconds, float releaseSeconds) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("modal constraint sample rate must be positive");
    sampleRate_ = sampleRate;
    referenceSamples_ = std::max(1.f, referenceTraversalSeconds * sampleRate);
    attackCoefficient_ = SmoothingCoefficient(attackSeconds);
    releaseCoefficient_ = SmoothingCoefficient(releaseSeconds);
    Reset();
  }

  void Reset(const PassiveConstraintGains gains = {}) noexcept {
    SetTarget(gains);
    current_ = target_;
  }

  void SetTarget(const PassiveConstraintGains gains) noexcept {
    const auto safe = SanitizePassiveConstraint(gains);
    target_ = {PerSample(safe.broadband), PerSample(safe.low),
               PerSample(safe.middle), PerSample(safe.high)};
  }

  ModalDampingGains Process() noexcept {
    current_.broadband = Follow(current_.broadband, target_.broadband);
    current_.low = Follow(current_.low, target_.low);
    current_.middle = Follow(current_.middle, target_.middle);
    current_.high = Follow(current_.high, target_.high);
    return current_;
  }

private:
  float PerSample(const float traversalGain) const noexcept {
    return traversalGain <= 0.f
        ? 0.f
        : std::exp(std::log(traversalGain) / referenceSamples_);
  }

  float SmoothingCoefficient(float seconds) const noexcept {
    seconds = std::clamp(std::isfinite(seconds) ? seconds : 0.f, 0.f, 10.f);
    return seconds == 0.f ? 1.f : -std::expm1(-1.f / (seconds * sampleRate_));
  }

  float Follow(const float current, const float target) const noexcept {
    const float coefficient = target < current ? attackCoefficient_
                                               : releaseCoefficient_;
    return tfdsp::FiniteNormalOrZero(current + coefficient * (target - current));
  }

  ModalDampingGains current_{};
  ModalDampingGains target_{};
  float sampleRate_{48000.f};
  float referenceSamples_{48.f};
  float attackCoefficient_{1.f};
  float releaseCoefficient_{1.f};
};

} // namespace tfdsp::percussion
