#pragma once

#include "deterministic_random.hpp"
#include "spectral_tilt_filter.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace tfdsp::percussion {

struct EnvelopedNoiseBurstParameters {
  float attackSeconds{.0002f};
  float holdSeconds{.001f};
  float decaySeconds{.008f};
  float amplitude{1.f};
  float tiltDb{0.f};
  float tiltPivotHz{3000.f};
  std::uint32_t seed{1};
};

// Seeded white noise with a smooth attack, optional hold, and exponential
// decay. Spectral shaping is deliberately separate from this source primitive.
class EnvelopedNoiseBurst {
public:
  void Prepare(const float sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("noise-burst sample rate must be positive");
    sampleRate_ = sampleRate;
    tilt_.Prepare(sampleRate);
    Reset();
  }

  void Reset() noexcept {
    sample_ = attackSamples_ = holdSamples_ = decaySamples_ = 0;
    decayEnvelope_ = 1.f;
    tilt_.Reset();
  }

  void Trigger(const EnvelopedNoiseBurstParameters &parameters) noexcept {
    Reset();
    attackSamples_ = TimeToSamples(parameters.attackSeconds);
    holdSamples_ = TimeToSamples(parameters.holdSeconds);
    decaySamples_ = TimeToSamples(parameters.decaySeconds);
    if (attackSamples_ + holdSamples_ + decaySamples_ == 0)
      decaySamples_ = 1;
    amplitude_ = std::max(0.f, FiniteOr(parameters.amplitude, 0.f));
    tilt_.SetTilt(parameters.tiltDb, parameters.tiltPivotHz);
    random_.Seed(parameters.seed);
    constexpr float EndLevel = 1.e-4f;
    decayMultiplier_ = std::pow(
        EndLevel, 1.f / static_cast<float>(
            std::max<std::size_t>(1, decaySamples_)));
  }

  float Process() noexcept {
    if (!Active())
      return 0.f;
    float envelope = 1.f;
    if (sample_ < attackSamples_) {
      const float position = static_cast<float>(sample_ + 1) /
                             static_cast<float>(attackSamples_);
      envelope = position * position * (3.f - 2.f * position);
    } else if (sample_ >= attackSamples_ + holdSamples_) {
      decayEnvelope_ *= decayMultiplier_;
      envelope = decayEnvelope_;
    }
    ++sample_;
    return amplitude_ * envelope * tilt_.Process(random_.Bipolar());
  }

  bool Active() const noexcept {
    return sample_ < attackSamples_ + holdSamples_ + decaySamples_;
  }

private:
  static float FiniteOr(const float value, const float fallback) noexcept {
    return std::isfinite(value) ? value : fallback;
  }

  std::size_t TimeToSamples(const float seconds) const noexcept {
    const float bounded = std::clamp(FiniteOr(seconds, 0.f), 0.f, 10.f);
    return static_cast<std::size_t>(std::lround(bounded * sampleRate_));
  }

  DeterministicRandom random_{};
  SpectralTiltFilter tilt_{};
  float sampleRate_{48000.f};
  float amplitude_{1.f};
  float decayEnvelope_{1.f};
  float decayMultiplier_{};
  std::size_t sample_{};
  std::size_t attackSamples_{};
  std::size_t holdSamples_{};
  std::size_t decaySamples_{};
};

} // namespace tfdsp::percussion
