#pragma once

#include "modulated_fractional_delay.hpp"
#include "tfdsp/finite_audio.hpp"

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

// Audio-rate nonlinear moving-delay core. Delay dimensions are expressed at
// this core's rate; the public oversampled wrapper performs host conversion.
class SelfPhaseDelayCore {
public:
  void Prepare(const float sampleRate, const float maximumDelaySamples) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("self-phase core rate must be positive");
    sampleRate_ = sampleRate;
    maximumDelaySamples_ = maximumDelaySamples;
    delay_.Prepare(maximumDelaySamples);
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
    // Above eight the bounded delay is almost a hard switch and the 2x path no
    // longer converges adequately to the 4x reference. Keep the declared
    // production range inside the measured nonlinear bandwidth.
    drive_ = std::clamp(FiniteOr(parameters.drive, 0.f), 0.f, 8.f);
    normalization_ = std::clamp(FiniteOr(parameters.normalization, 0.f), 0.f, 1.f);
    const float toneHz = std::clamp(FiniteOr(parameters.toneHz, 0.f), 0.f,
                                    .45f * sampleRate_);
    toneCoefficient_ = 1.f - std::exp(-6.283185307179586f * toneHz / sampleRate_);
    const float release = std::clamp(
        FiniteOr(parameters.envelopeReleaseSeconds, .01f), 1.f / sampleRate_, 1.f);
    envelopeRelease_ = std::exp(-1.f / (release * sampleRate_));
  }

  float Process(float input) noexcept {
    input = tfdsp::FiniteNormalOrZero(input);
    toneState_ += toneCoefficient_ * (input - toneState_);
    envelope_ = std::max(std::abs(input), envelopeRelease_ * envelope_);
    toneState_ = tfdsp::FiniteNormalOrZero(toneState_);
    envelope_ = tfdsp::FiniteNormalOrZero(envelope_);
    const float scale = std::pow(1.e-5f + envelope_, normalization_);
    const float offset = excursionSamples_ * std::tanh(drive_ * toneState_ / scale);
    return delay_.Process(input, centreDelaySamples_ + offset);
  }

  float CentreDelaySamples() const noexcept { return centreDelaySamples_; }

private:
  static float FiniteOr(const float value, const float fallback) noexcept {
    return std::isfinite(value) ? value : fallback;
  }

  ModulatedFractionalDelay delay_{};
  float sampleRate_{96000.f};
  float maximumDelaySamples_{128.f};
  float centreDelaySamples_{24.f};
  float excursionSamples_{4.f};
  float drive_{1.f};
  float normalization_{.5f};
  float toneCoefficient_{};
  float envelopeRelease_{};
  float toneState_{};
  float envelope_{};
};

} // namespace tfdsp::percussion
