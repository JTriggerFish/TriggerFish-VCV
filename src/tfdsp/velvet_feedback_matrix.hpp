#pragma once

#include "cubic_lagrange_interpolator.hpp"
#include "finite_audio.hpp"
#include "transition_fractional_delay.hpp"
#include "late_reverb_coefficient_sets.hpp"
#include "multiband_decay_filter.hpp"
#include "smooth_random_modulator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace tfdsp {

// Factorized velvet feedback matrix U2 D2(z) U1 D1(z) U0. Every transform uses
// orthogonal variable-angle butterflies followed by signed permutations.
// Diffusion controls both the butterfly angle and the temporal span of the
// pure-delay banks. With modulation at zero the complete operator is
// paraunitary; bounded fractional modulation deliberately makes it time
// varying to decorrelate the late field.
class VelvetFeedbackMatrix {
public:
  static constexpr std::size_t LineCount =
      late_reverb_coefficients::LineCount;
  using Frame = std::array<float, LineCount>;

private:
  std::array<std::array<TransitionFractionalDelay, LineCount>,
             late_reverb_coefficients::VelvetStageCount>
      delays_{};
  std::array<std::array<SmoothRandomModulator, LineCount>,
             late_reverb_coefficients::VelvetStageCount>
      modulators_{};
  std::array<std::array<MultibandDecayFilter, LineCount>,
             late_reverb_coefficients::VelvetStageCount>
      decayFilters_{};
  float sampleRate_{48'000.f};
  float lastDiffusion_{-1.f};
  float lastRoomScale_{-1.f};
  float currentMixAngle_{0.7853981633974483f};
  float targetMixAngle_{0.7853981633974483f};
  float mixCosine_{0.7071067811865475f};
  float mixSine_{0.7071067811865475f};
  float mixAngleSmoothingCoefficient_{1.f};
  std::size_t rotationNormalizationCountdown_{256};
  bool mixAngleInitialized_{};
  LateReverbFlavour transitionFrom_{DefaultLateReverbFlavour};
  LateReverbFlavour transitionTarget_{DefaultLateReverbFlavour};
  float flavourTransitionPhase_{1.f};
  float flavourTransitionIncrement_{1.f};

  // The short first bank gets the smaller relative and absolute excursion.
  // Both limits are applied per base tap; the fractional reader also reserves
  // a causal two-sample margin for exceptionally short low-Diffusion taps.
  static constexpr std::array<float,
                              late_reverb_coefficients::VelvetStageCount>
      ModulationDepthRatio{0.08f, 0.12f};
  static constexpr std::array<float,
                              late_reverb_coefficients::VelvetStageCount>
      MaximumModulationSeconds{0.00035f, 0.00065f};

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

  static Frame TransformFor(Frame frame, const std::size_t transform,
                            const LateReverbFlavour flavour,
                            const float cosine,
                            const float sine) noexcept {
    for (std::size_t stride = 1; stride < LineCount; stride *= 2)
      ButterflyLayer(frame, stride, cosine, sine);
    const auto &signs = LateSigns(flavour);
    const auto &permutations = LatePermutations(flavour);
    Frame permuted{};
    for (std::size_t line = 0; line < LineCount; ++line)
      permuted[line] = signs[transform][line] *
                       frame[permutations[transform][line]];
    return permuted;
  }

  Frame Transform(const Frame &frame, const std::size_t transform) const
      noexcept {
    if (flavourTransitionPhase_ >= 1.f ||
        transitionFrom_ == transitionTarget_)
      return TransformFor(frame, transform, transitionTarget_,
                          mixCosine_, mixSine_);
    const float phase = flavourTransitionPhase_;
    const float fade = phase * phase * (3.f - 2.f * phase);
    const auto from = TransformFor(frame, transform, transitionFrom_,
                                   mixCosine_, mixSine_);
    const auto to = TransformFor(frame, transform, transitionTarget_,
                                 mixCosine_, mixSine_);
    Frame result{};
    for (std::size_t line = 0; line < LineCount; ++line)
      result[line] = from[line] + fade * (to[line] - from[line]);
    return result;
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
    constexpr float MinimumMixAngle = 0.1963495408493621f; // pi / 16
    constexpr float MaximumMixAngle = 0.7853981633974483f; // pi / 4
    targetMixAngle_ =
        MinimumMixAngle + smooth * (MaximumMixAngle - MinimumMixAngle);
    if (!mixAngleInitialized_) {
      currentMixAngle_ = targetMixAngle_;
      mixCosine_ = std::cos(currentMixAngle_);
      mixSine_ = std::sin(currentMixAngle_);
      mixAngleInitialized_ = true;
    }
    const float scale = (0.20f + 1.30f * smooth) * roomScale;
    const auto &velvetDelays =
        LateVelvetDelayMilliseconds(transitionTarget_);
    for (std::size_t stage = 0; stage < delays_.size(); ++stage)
      for (std::size_t line = 0; line < LineCount; ++line) {
        const float seconds = velvetDelays[stage][line] * scale * 0.001f;
        delays_[stage][line].SetTarget(static_cast<std::size_t>(
            std::max(1.f, std::round(seconds * sampleRate_))));
      }
  }

  void AdvanceMixRotation() noexcept {
    const float previous = currentMixAngle_;
    currentMixAngle_ += mixAngleSmoothingCoefficient_ *
                        (targetMixAngle_ - currentMixAngle_);
    const float delta = currentMixAngle_ - previous;
    if (delta == 0.f)
      return;
    const float square = delta * delta;
    const float deltaSine = delta * (1.f - square / 6.f);
    const float deltaCosine = 1.f - .5f * square;
    const float nextCosine =
        mixCosine_ * deltaCosine - mixSine_ * deltaSine;
    mixSine_ = mixSine_ * deltaCosine + mixCosine_ * deltaSine;
    mixCosine_ = nextCosine;
    if (--rotationNormalizationCountdown_ == 0) {
      const float inverse = 1.f / std::sqrt(
          std::max(1.e-20f, mixCosine_ * mixCosine_ + mixSine_ * mixSine_));
      mixCosine_ *= inverse;
      mixSine_ *= inverse;
      rotationNormalizationCountdown_ = 256;
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
    flavourTransitionIncrement_ = transitionIncrement;
    mixAngleSmoothingCoefficient_ =
        1.f - std::exp(-1.f / (TransitionSeconds * sampleRate_));
    const auto &baseDelays =
        LateVelvetDelayMilliseconds(LateReverbFlavour::Base);
    const auto &optimizedDelays =
        LateVelvetDelayMilliseconds(LateReverbFlavour::Optimized);
    for (std::size_t stage = 0; stage < delays_.size(); ++stage)
      for (std::size_t line = 0; line < LineCount; ++line) {
        const float maximumDelay =
            std::max(baseDelays[stage][line], optimizedDelays[stage][line]);
        const auto capacity = static_cast<std::size_t>(std::ceil(
                                  maximumDelay * MaximumDiffusionScale *
                                      MaximumRoomScale * 0.001f * sampleRate_ +
                                  MaximumModulationSeconds[stage] *
                                      sampleRate_)) +
                              4;
        delays_[stage][line].Prepare(capacity, transitionIncrement);
        decayFilters_[stage][line].Prepare(sampleRate);
        const std::size_t rateIndex =
            (5 * line + 7 * stage) % late_reverb_coefficients::LineCount;
        const float rateMultiplier = stage == 0 ? 1.37f : 1.11f;
        modulators_[stage][line].Prepare(
            sampleRate,
            rateMultiplier *
                late_reverb_coefficients::ModulationRateHz[rateIndex],
            0xa511e9b3u + static_cast<std::uint32_t>(stage) * 0x63d83595u +
                static_cast<std::uint32_t>(line) * 0x9e3779b9u);
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
    for (auto &stage : modulators_)
      for (auto &modulator : stage)
        modulator.Reset();
    lastDiffusion_ = -1.f;
    lastRoomScale_ = -1.f;
    currentMixAngle_ = targetMixAngle_ = 0.7853981633974483f;
    mixCosine_ = mixSine_ = 0.7071067811865475f;
    rotationNormalizationCountdown_ = 256;
    mixAngleInitialized_ = false;
    transitionFrom_ = transitionTarget_;
    flavourTransitionPhase_ = 1.f;
  }

  void SetFlavour(const LateReverbFlavour flavour) noexcept {
    if (flavour == transitionTarget_)
      return;
    if (flavourTransitionPhase_ >= 1.f)
      transitionFrom_ = transitionTarget_;
    else if (flavourTransitionPhase_ >= 0.5f)
      transitionFrom_ = transitionTarget_;
    transitionTarget_ = flavour;
    flavourTransitionPhase_ = 0.f;
    lastDiffusion_ = -1.f;
    lastRoomScale_ = -1.f;
  }

  void SetFlavourImmediate(const LateReverbFlavour flavour) noexcept {
    transitionFrom_ = transitionTarget_ = flavour;
    flavourTransitionPhase_ = 1.f;
    lastDiffusion_ = -1.f;
    lastRoomScale_ = -1.f;
  }

  LateReverbFlavour Flavour() const noexcept { return transitionTarget_; }

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
                const float roomScale = 1.f,
                const float modulationAmount = 0.f) noexcept {
    UpdateDelayTargets(diffusion, roomScale);
    const float safeModulationAmount =
        std::clamp(std::isfinite(modulationAmount) ? modulationAmount : 0.f,
                   0.f, 1.f);
    AdvanceMixRotation();
    std::array<Frame, late_reverb_coefficients::VelvetStageCount>
        modulation{};
    for (std::size_t stage = 0; stage < modulation.size(); ++stage) {
      // Remove coherent delay motion across a complete bank. Rescaling only
      // when needed retains that zero mean while keeping every offset bounded.
      float mean = 0.f;
      for (std::size_t line = 0; line < LineCount; ++line) {
        modulation[stage][line] = modulators_[stage][line].Next();
        mean += modulation[stage][line];
      }
      mean /= static_cast<float>(LineCount);
      float maximum = 1.f;
      for (float &value : modulation[stage]) {
        value -= mean;
        maximum = std::max(maximum, std::abs(value));
      }
      for (float &value : modulation[stage])
        value = safeModulationAmount * value / maximum;
    }
    frame = Transform(frame, 0);
    for (std::size_t stage = 0; stage < delays_.size(); ++stage) {
      for (std::size_t line = 0; line < LineCount; ++line) {
        frame[line] = delays_[stage][line].Process(
            frame[line], modulation[stage][line],
            ModulationDepthRatio[stage],
            MaximumModulationSeconds[stage] * sampleRate_);
        const float seconds =
            static_cast<float>(delays_[stage][line].EffectiveSamples()) /
            sampleRate_;
        frame[line] = decayFilters_[stage][line].Process(
            frame[line], seconds, lowT60, midT60, highT60);
      }
      frame = Transform(frame, stage + 1);
    }
    if (flavourTransitionPhase_ < 1.f) {
      flavourTransitionPhase_ =
          std::min(1.f,
                   flavourTransitionPhase_ + flavourTransitionIncrement_);
      if (flavourTransitionPhase_ >= 1.f)
        transitionFrom_ = transitionTarget_;
    }
    return frame;
  }
};

} // namespace tfdsp
