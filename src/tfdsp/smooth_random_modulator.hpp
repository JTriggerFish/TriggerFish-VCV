#pragma once

#include <algorithm>
#include <cstdint>

namespace tfdsp {

// Deterministic, band-limited random motion for decorrelating delay lines.
// Smoothstep interpolation gives zero slope at every randomly scheduled turn,
// avoiding the corners and strong periodicity of a conventional LFO.
class SmoothRandomModulator {
public:
  void Prepare(const double sampleRate, const float rate,
               const std::uint32_t seed) noexcept {
    sampleRate_ = static_cast<float>(sampleRate);
    rate_ = std::max(rate, 0.01f);
    seed_ = seed != 0 ? seed : 0x9e3779b9u;
    Reset();
  }

  void Reset() noexcept {
    state_ = seed_;
    current_ = Random();
    target_ = Random();
    phase_ = 0.f;
    Schedule();
  }

  float Next() noexcept {
    if (phase_ >= 1.f) {
      current_ = target_;
      target_ = Random();
      phase_ = 0.f;
      Schedule();
    }
    const float fade = phase_ * phase_ * (3.f - 2.f * phase_);
    const float value = current_ + fade * (target_ - current_);
    phase_ += increment_;
    return value;
  }

private:
  float Random() noexcept {
    state_ ^= state_ << 13;
    state_ ^= state_ >> 17;
    state_ ^= state_ << 5;
    return 2.f * static_cast<float>(state_ & 0x00ffffffu) /
               static_cast<float>(0x00ffffffu) -
           1.f;
  }

  void Schedule() noexcept {
    const float duration = (0.35f + 0.30f * (0.5f * Random() + 0.5f)) / rate_;
    increment_ = 1.f / std::max(1.f, duration * sampleRate_);
  }

  float sampleRate_{48'000.f};
  float rate_{0.1f};
  float current_{};
  float target_{};
  float phase_{};
  float increment_{1.f};
  std::uint32_t seed_{0x9e3779b9u};
  std::uint32_t state_{seed_};
};

} // namespace tfdsp
