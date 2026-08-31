#pragma once

#include "cubic_lagrange_interpolator.hpp"
#include "finite_audio.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace tfdsp {
class CubicFractionalDelay {
public:
  void Prepare(const std::size_t capacity) {
    if (capacity < 8)
      throw std::invalid_argument("fractional delay capacity is too small");
    buffer_.assign(capacity, 0.f);
    writeIndex_ = 0;
  }
  void Reset() noexcept {
    std::fill(buffer_.begin(), buffer_.end(), 0.f);
    writeIndex_ = 0;
  }
  float Read(float delaySamples) const noexcept {
    if (buffer_.empty())
      return 0.f;
    return ReadCubicLagrange(delaySamples, buffer_.size(),
                            [&](const std::size_t distance) {
      const auto index = writeIndex_ >= distance
                             ? writeIndex_ - distance
                             : writeIndex_ + buffer_.size() - distance;
      return buffer_[index];
    });
  }
  void Push(const float value) noexcept {
    if (buffer_.empty())
      return;
    buffer_[writeIndex_] = FiniteNormalOrZero(value);
    if (++writeIndex_ == buffer_.size())
      writeIndex_ = 0;
  }

private:
  std::vector<float> buffer_{};
  std::size_t writeIndex_{};
};

class LiveFractionalDelay {
public:
  void Prepare(const std::size_t capacity) {
    if (capacity < 8)
      throw std::invalid_argument("room pre-delay capacity is too small");
    buffer_.assign(capacity, 0.f);
    writeIndex_ = 0;
  }

  void Reset() noexcept {
    std::fill(buffer_.begin(), buffer_.end(), 0.f);
    writeIndex_ = 0;
  }

  float Read(float delaySamples, const float liveSample) const noexcept {
    if (buffer_.empty())
      return 0.f;
    delaySamples = std::max(0.f, std::isfinite(delaySamples) ? delaySamples : 0.f);
    const float safeLiveSample = FiniteNormalOrZero(liveSample);
    if (delaySamples < 2.f) {
      const auto history = [&](const std::size_t distance) {
        const auto index = writeIndex_ >= distance
                               ? writeIndex_ - distance
                               : writeIndex_ + buffer_.size() - distance;
        return buffer_[index];
      };
      if (delaySamples < 1.f)
        return safeLiveSample + delaySamples * (history(1) - safeLiveSample);
      const float fraction = delaySamples - 1.f;
      return history(1) + fraction * (history(2) - history(1));
    }
    return ReadCubicLagrange(delaySamples, buffer_.size(),
                            [&](const std::size_t distance) {
      const auto index = writeIndex_ >= distance
                             ? writeIndex_ - distance
                             : writeIndex_ + buffer_.size() - distance;
      return buffer_[index];
    });
  }

  void Push(const float value) noexcept {
    buffer_[writeIndex_] = FiniteNormalOrZero(value);
    if (++writeIndex_ == buffer_.size())
      writeIndex_ = 0;
  }

private:
  std::vector<float> buffer_{};
  std::size_t writeIndex_{};
};

} // namespace tfdsp
