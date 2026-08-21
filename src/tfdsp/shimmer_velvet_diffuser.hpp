#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace tfdsp {

// A four-channel, two-stage feed-forward velvet diffuser used exclusively by
// the shimmer layer. Every transform is orthonormal and every stage contains
// only pure delays, so the static operator is paraunitary: it decorrelates the
// pitch-shifted buses without changing their total energy or adding poles to
// the main room network.
class ShimmerVelvetDiffuser {
public:
  static constexpr std::size_t BusCount = 4;
  static constexpr std::size_t StageCount = 2;
  using Frame = std::array<float, BusCount>;

private:
  static constexpr float MaximumRoomScale = 2.25f;
  static constexpr std::array<std::array<float, BusCount>, StageCount>
      DelayMs{{{{2.9f, 4.7f, 7.1f, 10.3f}},
               {{3.7f, 6.1f, 8.9f, 13.1f}}}};

  class TransitionDelay {
  public:
    void Prepare(const std::size_t capacity, const float transitionIncrement) {
      if (capacity < 3)
        throw std::invalid_argument("shimmer diffuser delay is too small");
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

    void SetTarget(const std::size_t samples) noexcept {
      const auto target =
          std::clamp(samples, std::size_t{1}, buffer_.size() - 1);
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

    float Process(const float input) noexcept {
      float output = Read(current_);
      if (phase_ < 1.f) {
        const float fade = phase_ * phase_ * (3.f - 2.f * phase_);
        output = Read(from_) + fade * (Read(to_) - Read(from_));
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
      buffer_[writeIndex_] = std::isfinite(input) ? input : 0.f;
      if (++writeIndex_ == buffer_.size())
        writeIndex_ = 0;
      return output;
    }

  private:
    float Read(const std::size_t samples) const noexcept {
      return buffer_[(writeIndex_ + buffer_.size() - samples) % buffer_.size()];
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

  std::array<std::array<TransitionDelay, BusCount>, StageCount> delays_{};
  float sampleRate_{48'000.f};
  float lastRoomScale_{-1.f};

  static Frame Transform(const Frame &input, const std::size_t transform) {
    constexpr float Half = 0.5f;
    Frame mixed{{Half * (input[0] + input[1] + input[2] + input[3]),
                 Half * (input[0] - input[1] + input[2] - input[3]),
                 Half * (input[0] + input[1] - input[2] - input[3]),
                 Half * (input[0] - input[1] - input[2] + input[3])}};
    constexpr std::array<std::array<std::size_t, BusCount>, StageCount + 1>
        Permutation{{{{2, 0, 3, 1}}, {{1, 3, 0, 2}}, {{3, 2, 1, 0}}}};
    constexpr std::array<std::array<float, BusCount>, StageCount + 1> Sign{{
        {{1.f, -1.f, 1.f, 1.f}},
        {{-1.f, 1.f, 1.f, -1.f}},
        {{1.f, 1.f, -1.f, 1.f}},
    }};
    Frame result{};
    for (std::size_t bus = 0; bus < BusCount; ++bus)
      result[bus] = Sign[transform][bus] *
                    mixed[Permutation[transform][bus]];
    return result;
  }

  void UpdateRoomScale(const float requestedScale) noexcept {
    const float scale =
        std::clamp(std::isfinite(requestedScale) ? requestedScale : 1.f,
                   0.35f, MaximumRoomScale);
    if (std::abs(scale - lastRoomScale_) < 0.002f)
      return;
    lastRoomScale_ = scale;
    for (std::size_t stage = 0; stage < StageCount; ++stage)
      for (std::size_t bus = 0; bus < BusCount; ++bus) {
        const float samples =
            DelayMs[stage][bus] * 0.001f * scale * sampleRate_;
        delays_[stage][bus].SetTarget(static_cast<std::size_t>(
            std::max(1.f, std::round(samples))));
      }
  }

public:
  void Prepare(const double sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument(
          "shimmer diffuser sample rate must be positive");
    sampleRate_ = static_cast<float>(sampleRate);
    constexpr float TransitionSeconds = 0.030f;
    const float increment =
        1.f / std::max(1.f, TransitionSeconds * sampleRate_);
    for (std::size_t stage = 0; stage < StageCount; ++stage)
      for (std::size_t bus = 0; bus < BusCount; ++bus) {
        const auto capacity = static_cast<std::size_t>(std::ceil(
                                  DelayMs[stage][bus] * 0.001f *
                                  MaximumRoomScale * sampleRate_)) +
                              3;
        delays_[stage][bus].Prepare(capacity, increment);
      }
    Reset();
  }

  void Reset() noexcept {
    for (auto &stage : delays_)
      for (auto &delay : stage)
        delay.Reset();
    lastRoomScale_ = -1.f;
  }

  Frame Process(Frame frame, const float roomScale) noexcept {
    UpdateRoomScale(roomScale);
    frame = Transform(frame, 0);
    for (std::size_t stage = 0; stage < StageCount; ++stage) {
      for (std::size_t bus = 0; bus < BusCount; ++bus)
        frame[bus] = delays_[stage][bus].Process(frame[bus]);
      frame = Transform(frame, stage + 1);
    }
    return frame;
  }
};

} // namespace tfdsp
