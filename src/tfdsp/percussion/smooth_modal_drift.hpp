#pragma once

#include "deterministic_random.hpp"
#include "tfdsp/finite_audio.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace tfdsp::percussion {

// Independent, C2-continuous random frequency trajectories. Depth is a peak
// fractional frequency deviation, not linewidth or a change to the pole radius.
template <std::size_t Count> class SmoothModalDrift {
public:
  // Absolute-Hz movement for narrow, irregular ringing. Unlike fractional
  // drift it neither suppresses clean modes nor widens the treble in proportion
  // to pitch. Boundary compression keeps every instantaneous pole frequency safe.
  void PrepareHz(float sampleRate, const std::array<float, Count> &frequency,
                 float depthHz, float knotsPerSecond, std::uint32_t seed) {
    std::array<float, Count> amount{};
    amount.fill(1.f);
    Prepare(sampleRate, frequency, amount, 0.f, knotsPerSecond, seed);
    const float depth = std::clamp(tfdsp::FiniteNormalOrZero(depthHz), 0.f, 100.f);
    enabled_ = depth > 0.f;
    for (std::size_t i = 0; i < Count; ++i) {
      const float centre = std::clamp(tfdsp::FiniteNormalOrZero(frequency[i]),
                                     0.f, .499f * sampleRate);
      const float excursion = std::min({depth, centre, .499f * sampleRate-centre});
      angleDepth_[i] = 6.28318530718f * excursion / sampleRate;
    }
  }

  void Prepare(float sampleRate, const std::array<float, Count> &frequency,
               const std::array<float, Count> &amount, float depthPercent,
               float knotsPerSecond, std::uint32_t seed) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("drift sample rate must be positive");
    const float depth = .01f * std::clamp(
        tfdsp::FiniteNormalOrZero(depthPercent), 0.f, 10.f);
    increment_ = std::min(.25f, std::clamp(
        tfdsp::FiniteNormalOrZero(knotsPerSecond), .1f, 40.f) / sampleRate);
    seed_ = seed;
    enabled_ = depth > 0.f;
    for (std::size_t i = 0; i < Count; ++i) {
      // Keep instantaneous centres inside the positive/Nyquist boundaries.
      const float centre = std::clamp(tfdsp::FiniteNormalOrZero(frequency[i]),
                                      0.f, .499f * sampleRate);
      const float excursion = std::min(centre * depth,
          std::max(0.f, .499f * sampleRate - centre));
      angleDepth_[i] = 6.28318530718f * excursion / sampleRate * std::clamp(
          tfdsp::FiniteNormalOrZero(amount[i]), 0.f, 1.f);
    }
    Reset();
  }

  void Reset() noexcept {
    random_.Seed(seed_);
    for (std::size_t i = 0; i < Count; ++i) {
      phase_[i] = random_.Uniform();
      from_[i] = random_.Bipolar();
      to_[i] = random_.Bipolar();
    }
  }

  bool Enabled() const noexcept { return enabled_; }

  float NextAngle(std::size_t i) noexcept {
    phase_[i] += increment_;
    if (phase_[i] >= 1.f) {
      phase_[i] -= 1.f;
      from_[i] = to_[i];
      to_[i] = random_.Bipolar();
    }
    const float t = phase_[i];
    const float weight = t*t*t * (10.f + t * (-15.f + 6.f*t));
    return angleDepth_[i] * (from_[i] + weight * (to_[i] - from_[i]));
  }

  // Cayley rotation has unit length algebraically. The tan approximation is
  // accurate over our bounded +/-0.31 rad range; no per-sample trig is needed.
  void RotateCoefficients(std::size_t i, float &cosine, float &sine) noexcept {
    const float angle = NextAngle(i);
    const float a2 = angle * angle;
    const float tangent = .5f * angle *
        (1.f + a2 * (1.f/12.f + a2 * (1.f/120.f + a2 * (17.f/20160.f))));
    const float inverse = 1.f / (1.f + tangent*tangent);
    const float c = (1.f - tangent*tangent) * inverse;
    const float s = 2.f * tangent * inverse;
    const float prior = cosine;
    cosine = c * prior - s * sine;
    sine = s * prior + c * sine;
  }

private:
  std::array<float, Count> phase_{}, from_{}, to_{}, angleDepth_{};
  DeterministicRandom random_{};
  std::uint32_t seed_{};
  float increment_{};
  bool enabled_{};
};

} // namespace tfdsp::percussion
