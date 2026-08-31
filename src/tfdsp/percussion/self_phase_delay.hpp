#pragma once

#include "modulated_fractional_delay.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tfdsp::percussion {

struct SelfPhaseDelayParameters {
  float centreDelaySamples{12.f};
  float maximumExcursionSamples{2.f};
  float drive{1.f};
  float toneHz{5000.f};
  float envelopeReleaseSeconds{.01f};
  float normalization{.5f};
};

// Signed audio controls a bounded moving delay, creating correlated sidebands.
// The envelope only controls partial level normalization; it is not the PM source.
class SelfPhaseDelay {
public:
  void Prepare(const float sampleRate, const float maximumDelaySamples) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("self-phase sample rate must be positive");
    sampleRate_ = sampleRate;
    maximumDelaySamples_ = maximumDelaySamples;
    delay_.Prepare(maximumDelaySamples);
    SetParameters({});
    Reset();
  }

  void Reset() noexcept {
    delay_.Reset();
    toneState_ = envelope_ = 0.f;
  }

  void SetParameters(const SelfPhaseDelayParameters &parameters) noexcept {
    centreDelaySamples_ = std::clamp(FiniteOr(parameters.centreDelaySamples, 12.f),
        ModulatedFractionalDelay::MinimumDelaySamples, maximumDelaySamples_);
    const float margin = std::min(
        centreDelaySamples_ - ModulatedFractionalDelay::MinimumDelaySamples,
        maximumDelaySamples_ - centreDelaySamples_);
    excursionSamples_ = std::clamp(
        FiniteOr(parameters.maximumExcursionSamples, 0.f), 0.f, std::max(0.f, margin));
    drive_ = std::clamp(FiniteOr(parameters.drive, 0.f), 0.f, 32.f);
    normalization_ = std::clamp(FiniteOr(parameters.normalization, 0.f), 0.f, 1.f);
    const float toneHz = std::clamp(FiniteOr(parameters.toneHz, 0.f), 0.f,
                                    .45f * sampleRate_);
    toneCoefficient_ = 1.f - std::exp(-6.283185307179586f * toneHz / sampleRate_);
    const float release = std::clamp(
        FiniteOr(parameters.envelopeReleaseSeconds, .01f), 1.f / sampleRate_, 1.f);
    envelopeRelease_ = std::exp(-1.f / (release * sampleRate_));
  }

  float Process(float input) noexcept {
    if (!std::isfinite(input))
      input = 0.f;
    toneState_ += toneCoefficient_ * (input - toneState_);
    envelope_ = std::max(std::abs(input), envelopeRelease_ * envelope_);
    const float denominator = std::pow(1.e-5f + envelope_, normalization_);
    const float control = toneState_ / denominator;
    const float offset = excursionSamples_ * std::tanh(drive_ * control);
    return delay_.Process(input, centreDelaySamples_ + offset);
  }

  float CentreDelaySamples() const noexcept { return centreDelaySamples_; }

private:
  static float FiniteOr(const float value, const float fallback) noexcept {
    return std::isfinite(value) ? value : fallback;
  }

  ModulatedFractionalDelay delay_{};
  float sampleRate_{48000.f};
  float maximumDelaySamples_{64.f};
  float centreDelaySamples_{12.f};
  float excursionSamples_{2.f};
  float drive_{1.f};
  float normalization_{.5f};
  float toneCoefficient_{};
  float envelopeRelease_{};
  float toneState_{};
  float envelope_{};
};

} // namespace tfdsp::percussion
