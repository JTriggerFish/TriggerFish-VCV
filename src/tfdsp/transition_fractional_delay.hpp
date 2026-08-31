#pragma once

#include "cubic_lagrange_interpolator.hpp"
#include "finite_audio.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace tfdsp {
class TransitionFractionalDelay {
public:
  void Prepare(const std::size_t capacity, const float transitionIncrement) {
    if (capacity < 8)
      throw std::invalid_argument("VFM delay capacity is too small");
    buffer_.assign(capacity, 0.f);
    transitionIncrement_ = transitionIncrement;
    Reset();
  }
  void Reset() noexcept {
    std::fill(buffer_.begin(), buffer_.end(), 0.f);
    writeIndex_ = 0;
    current_ = from_ = to_ = pending_ = 1;
    phase_ = 1.f;
    initialized_ = false;
  }
  void SetTarget(const std::size_t delaySamples) noexcept {
    const std::size_t target =
        std::clamp(delaySamples, std::size_t{1}, buffer_.size() - 1);
    if (!initialized_) {
      current_ = from_ = to_ = pending_ = target;
      initialized_ = true;
      return;
    }
    pending_ = target;
    if (phase_ < 1.f || pending_ == current_)
      return;
    from_ = current_;
    to_ = pending_;
    phase_ = 0.f;
  }
  float Process(const float input, const float normalizedModulation,
                const float relativeDepth,
                const float maximumDepthSamples) noexcept {
    float output = ReadModulated(current_, normalizedModulation,
                                 relativeDepth, maximumDepthSamples);
    if (phase_ < 1.f) {
      const float fade = phase_ * phase_ * (3.f - 2.f * phase_);
      const float from = ReadModulated(from_, normalizedModulation,
                                       relativeDepth, maximumDepthSamples);
      const float to = ReadModulated(to_, normalizedModulation,
                                     relativeDepth, maximumDepthSamples);
      output = from + fade * (to - from);
      phase_ = std::min(1.f, phase_ + transitionIncrement_);
      if (phase_ >= 1.f) {
        current_ = to_;
        if (pending_ != current_) {
          from_ = current_;
          to_ = pending_;
          phase_ = 0.f;
        }
      }
    }
    buffer_[writeIndex_] = FiniteNormalOrZero(input);
    if (++writeIndex_ == buffer_.size())
      writeIndex_ = 0;
    return output;
  }
  std::size_t EffectiveSamples() const noexcept {
    return phase_ < 1.f ? std::max(from_, to_) : current_;
  }

private:
  float ReadInteger(const std::size_t samples) const noexcept {
    const std::size_t distance =
        std::clamp(samples, std::size_t{1}, buffer_.size() - 1);
    const std::size_t index = writeIndex_ >= distance
                                  ? writeIndex_ - distance
                                  : writeIndex_ + buffer_.size() - distance;
    return buffer_[index];
  }

  float Read(float delaySamples) const noexcept {
    const float nearest = std::round(delaySamples);
    if (std::abs(delaySamples - nearest) < 1.e-6f)
      return ReadInteger(static_cast<std::size_t>(std::max(1.f, nearest)));

    // A centred four-point Lagrange read needs one older and one newer
    // neighbour around the integer tap. Modulation is disabled for base taps
    // too short to maintain this causal two-sample lower bound.
    return ReadCubicLagrange(delaySamples, buffer_.size(),
                            [&](const std::size_t distance) {
      const auto index = writeIndex_ >= distance
                             ? writeIndex_ - distance
                             : writeIndex_ + buffer_.size() - distance;
      return buffer_[index];
    });
  }

  float ReadModulated(const std::size_t baseSamples,
                      const float normalizedModulation,
                      const float relativeDepth,
                      const float maximumDepthSamples) const noexcept {
    const float base = static_cast<float>(baseSamples);
    const float causalDepth = std::max(0.f, base - 2.f);
    const float relative = std::max(
        0.f, std::isfinite(relativeDepth) ? relativeDepth : 0.f);
    const float maximum = std::max(
        0.f, std::isfinite(maximumDepthSamples) ? maximumDepthSamples : 0.f);
    const float depth = std::min(
        {relative * base, maximum, causalDepth});
    const float modulation = std::clamp(
        std::isfinite(normalizedModulation) ? normalizedModulation : 0.f,
        -1.f, 1.f);
    return Read(base + modulation * depth);
  }

  std::vector<float> buffer_{};
  std::size_t writeIndex_{};
  std::size_t current_{1};
  std::size_t from_{1};
  std::size_t to_{1};
  std::size_t pending_{1};
  float phase_{1.f};
  float transitionIncrement_{1.f};
  bool initialized_{};
};

} // namespace tfdsp
