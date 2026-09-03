#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tfdsp::percussion {

struct StrikeEnergyEnvelopeParameters {
  float releaseSeconds{.12f};
  float capacity{2.f};
  float tensionOctaves{.12f};
};

// Bounded, event-derived strike history used to modulate a body's tension.
// This is a control envelope, not a measurement of resonator-state energy.
class StrikeEnergyEnvelope {
public:
  void Prepare(const float sampleRate,
               const StrikeEnergyEnvelopeParameters &parameters) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("strike-energy sample rate must be positive");
    capacity_ = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters.capacity), .001f, 64.f);
    tensionOctaves_ = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters.tensionOctaves), -.5f, 1.f);
    const float release = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters.releaseSeconds), .001f, 30.f);
    releaseCoefficient_ = std::exp(-1.f / (release * sampleRate));
    Reset();
  }

  void Reset() noexcept { value_ = 0.f; }

  void InjectStrike(const float normalizedEnergy) noexcept {
    value_ = std::min(capacity_, value_ + std::clamp(
        tfdsp::FiniteNormalOrZero(normalizedEnergy), 0.f, capacity_));
  }

  float Process() noexcept {
    value_ = tfdsp::FiniteNormalOrZero(value_ * releaseCoefficient_);
    return value_;
  }

  float TensionScale() const noexcept {
    const float normalized = std::sqrt(std::clamp(
        value_ / capacity_, 0.f, 1.f));
    return std::exp2(tensionOctaves_ * normalized);
  }

  float Value() const noexcept { return value_; }

private:
  float value_{};
  float releaseCoefficient_{};
  float capacity_{2.f};
  float tensionOctaves_{.12f};
};

} // namespace tfdsp::percussion
