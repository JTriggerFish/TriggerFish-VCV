#pragma once

#include "late_reverb_coefficients.hpp"
#include "multiband_decay_filter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace tfdsp {

// Factorized velvet feedback matrix U2 D2(z) U1 D1(z) U0. All transforms are
// fixed, fully mixed Hadamard butterflies with signed permutations. Diffusion
// changes the temporal span of the pure-delay banks, not the scattering
// matrices. At any static setting the complete operator is paraunitary.
class VelvetFeedbackMatrix {
public:
  static constexpr std::size_t LineCount =
      late_reverb_coefficients::LineCount;
  using Frame = std::array<float, LineCount>;

private:
  class TransitionIntegerDelay {
  public:
    void Prepare(const std::size_t capacity, const float transitionIncrement) {
      if (capacity < 3)
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
    std::size_t EffectiveSamples() const noexcept {
      return phase_ < 1.f ? std::max(from_, to_) : current_;
    }

  private:
    float Read(const std::size_t samples) const noexcept {
      const std::size_t index =
          (writeIndex_ + buffer_.size() - samples) % buffer_.size();
      return buffer_[index];
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

  std::array<std::array<TransitionIntegerDelay, LineCount>,
             late_reverb_coefficients::VelvetStageCount>
      delays_{};
  std::array<std::array<MultibandDecayFilter, LineCount>,
             late_reverb_coefficients::VelvetStageCount>
      decayFilters_{};
  float sampleRate_{48'000.f};
  float lastDiffusion_{-1.f};
  float lastRoomScale_{-1.f};

  static void ButterflyLayer(Frame &frame, const std::size_t stride,
                             const float cosine, const float sine) noexcept {
    for (std::size_t base = 0; base < LineCount; base += 2 * stride)
      for (std::size_t offset = 0; offset < stride; ++offset) {
        const std::size_t first = base + offset;
        const std::size_t second = first + stride;
        const float left = frame[first];
        const float right = frame[second];
        frame[first] = cosine * left + sine * right;
        frame[second] = sine * left - cosine * right;
      }
  }

  static Frame Transform(Frame frame, const std::size_t transform) noexcept {
    constexpr float InverseRootTwo = 0.7071067811865475f;
    for (std::size_t stride = 1; stride < LineCount; stride *= 2)
      ButterflyLayer(frame, stride, InverseRootTwo, InverseRootTwo);
    Frame permuted{};
    for (std::size_t line = 0; line < LineCount; ++line)
      permuted[line] =
          late_reverb_coefficients::TransformSign[transform][line] *
          frame[late_reverb_coefficients::TransformPermutation[transform]
                                                        [line]];
    return permuted;
  }

  void UpdateDelayTargets(const float control,
                          const float requestedRoomScale) noexcept {
    const float x = std::clamp(control, 0.f, 1.f);
    const float roomScale =
        std::clamp(std::isfinite(requestedRoomScale) ? requestedRoomScale : 1.f,
                   0.35f, 2.25f);
    if (std::abs(x - lastDiffusion_) < 0.002f &&
        std::abs(roomScale - lastRoomScale_) < 0.002f)
      return;
    lastDiffusion_ = x;
    lastRoomScale_ = roomScale;
    const float smooth = x * x * (3.f - 2.f * x);
    const float scale = (0.20f + 1.30f * smooth) * roomScale;
    for (std::size_t stage = 0; stage < delays_.size(); ++stage)
      for (std::size_t line = 0; line < LineCount; ++line) {
        const float seconds =
            late_reverb_coefficients::VelvetDelayMs[stage][line] * scale *
            0.001f;
        delays_[stage][line].SetTarget(static_cast<std::size_t>(
            std::max(1.f, std::round(seconds * sampleRate_))));
      }
  }

public:
  void Prepare(const double sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument("VFM sample rate must be positive");
    sampleRate_ = static_cast<float>(sampleRate);
    constexpr float MaximumDiffusionScale = 1.5f;
    constexpr float MaximumRoomScale = 2.25f;
    // Size can move every sample under CV. A 100 ms fixed-tap crossfade keeps
    // the succession of discrete VFM target changes below audible pitch and
    // tremolo sidebands while matching the geometric ER transition time.
    constexpr float TransitionSeconds = 0.100f;
    const float transitionIncrement =
        1.f / std::max(1.f, TransitionSeconds * sampleRate_);
    for (std::size_t stage = 0; stage < delays_.size(); ++stage)
      for (std::size_t line = 0; line < LineCount; ++line) {
        const auto capacity = static_cast<std::size_t>(std::ceil(
                                  late_reverb_coefficients::
                                          VelvetDelayMs[stage][line] *
                                      MaximumDiffusionScale *
                                      MaximumRoomScale * 0.001f *
                                      sampleRate_)) +
                              3;
        delays_[stage][line].Prepare(capacity, transitionIncrement);
        decayFilters_[stage][line].Prepare(sampleRate);
      }
    Reset();
  }
  void Reset() noexcept {
    for (auto &stage : delays_)
      for (auto &delay : stage)
        delay.Reset();
    for (auto &stage : decayFilters_)
      for (auto &filter : stage)
        filter.Reset();
    lastDiffusion_ = -1.f;
    lastRoomScale_ = -1.f;
  }

  // Lossless reference path, used to verify paraunitarity independently of
  // the reverb's attenuation filters.
  Frame Process(Frame frame, const float diffusion) noexcept {
    return Process(frame, diffusion, 1.f);
  }

  Frame Process(Frame frame, const float diffusion,
                const float roomScale) noexcept {
    constexpr float infinity = std::numeric_limits<float>::infinity();
    return Process(frame, diffusion, infinity, infinity, infinity, roomScale);
  }

  // Apply gain-per-sample attenuation to every actual delay in the VFM.  This
  // is the cascaded form of homogeneous FFDN loss: the decay control changes
  // attenuation only, never the feedback topology.
  Frame Process(Frame frame, const float diffusion, const float midT60,
                const float highT60, const float lowT60,
                const float roomScale = 1.f) noexcept {
    UpdateDelayTargets(diffusion, roomScale);
    frame = Transform(frame, 0);
    for (std::size_t stage = 0; stage < delays_.size(); ++stage) {
      for (std::size_t line = 0; line < LineCount; ++line) {
        frame[line] = delays_[stage][line].Process(frame[line]);
        const float seconds =
            static_cast<float>(delays_[stage][line].EffectiveSamples()) /
            sampleRate_;
        frame[line] = decayFilters_[stage][line].Process(
            frame[line], seconds, lowT60, midT60, highT60);
      }
      frame = Transform(frame, stage + 1);
    }
    return frame;
  }
};

} // namespace tfdsp
