#pragma once

#include "late_reverb_coefficient_sets.hpp"
#include "multiband_decay_filter.hpp"
#include "reverb_defaults.hpp"
#include "smooth_random_modulator.hpp"
#include "velvet_feedback_matrix.hpp"
#include "windowed_pitch_shifter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace tfdsp {

struct LateReverbControls {
  float decay{reverb_defaults::Decay};
  float damping{reverb_defaults::Damping};
  float diffusion{reverb_defaults::Diffusion};
  float modulation{reverb_defaults::Modulation};
  float shimmer{reverb_defaults::Shimmer};
  std::array<float, 3> roomDimensionsMetres{
      reverb_defaults::RoomDimensionsMetres};
};

class LateReverb {
public:
  static constexpr std::size_t FdnLineCount =
      late_reverb_coefficients::LineCount;
  static constexpr std::size_t WallCount = late_reverb_coefficients::WallCount;
  using StereoFrame = std::array<float, 2>;
  using LineFrame = std::array<float, FdnLineCount>;
  using WallFrame = std::array<float, WallCount>;

private:
  static constexpr float Pi = 3.14159265358979323846f;
  static constexpr float SpeedOfSound = 343.f;
  static constexpr float MaximumFdnModulationSeconds = 0.000500f;
  // Shimmer adds a bounded octave return to four orthogonal coordinates inside
  // the main velvet loop without replacing the unshifted room feedback. The
  // remaining twelve coordinates are unaltered.
  static constexpr float MaximumShimmerLoopGain = 0.85f;
  static constexpr float MainDelayTransitionSeconds = 0.050f;
  static constexpr std::size_t ControlUpdateInterval = 64;

  class FractionalDelay {
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
      delaySamples =
          std::clamp(delaySamples, 2.f, static_cast<float>(buffer_.size() - 3));
      const auto integer = static_cast<std::size_t>(std::floor(delaySamples));
      const float mu = delaySamples - static_cast<float>(integer);
      const std::array<float, 4> c{-mu * (mu - 1.f) * (mu - 2.f) / 6.f,
                                   (mu + 1.f) * (mu - 1.f) * (mu - 2.f) / 2.f,
                                   -(mu + 1.f) * mu * (mu - 2.f) / 2.f,
                                   (mu + 1.f) * mu * (mu - 1.f) / 6.f};
      constexpr std::array<int, 4> offsets{-1, 0, 1, 2};
      float value = 0.f;
      for (std::size_t tap = 0; tap < c.size(); ++tap) {
        const auto distance =
            static_cast<std::size_t>(static_cast<int>(integer) + offsets[tap]);
        const auto index = writeIndex_ >= distance
                               ? writeIndex_ - distance
                               : writeIndex_ + buffer_.size() - distance;
        value += c[tap] * buffer_[index];
      }
      return value;
    }
    void Push(const float value) noexcept {
      if (buffer_.empty())
        return;
      buffer_[writeIndex_] = std::isfinite(value) ? value : 0.f;
      if (++writeIndex_ == buffer_.size())
        writeIndex_ = 0;
    }

  private:
    std::vector<float> buffer_{};
    std::size_t writeIndex_{};
  };

  double sampleRate_{48'000.0};
  std::array<FractionalDelay, FdnLineCount> mainDelays_{};
  VelvetFeedbackMatrix feedbackMatrix_{};
  std::array<SmoothRandomModulator, FdnLineCount> modulators_{};
  std::array<MultibandDecayFilter, FdnLineCount> mainDecayFilters_{};
  std::array<WindowedPitchShifter, late_reverb_coefficients::ShimmerBusCount>
      shimmerShifters_{};
  std::array<float, late_reverb_coefficients::ShimmerBusCount>
      shimmerHighpassState_{};
  std::array<std::array<float, late_reverb_coefficients::ShimmerBusCount>, 2>
      shimmerLowpassState_{};
  float shimmerHighpassAlpha_{};
  float shimmerLowpassAlpha_{};
  float currentMeanDelaySamples_{480.f};
  float fromMeanDelaySamples_{480.f};
  float toMeanDelaySamples_{480.f};
  float pendingMeanDelaySamples_{480.f};
  float mainDelayPhase_{1.f};
  float mainDelayIncrement_{1.f};
  bool mainDelayInitialized_{};
  LateMainDelayRatios mainRatioFrom_ =
      LateMainRatios(DefaultLateReverbFlavour);
  LateMainDelayRatios mainRatioTarget_ =
      LateMainRatios(DefaultLateReverbFlavour);
  float mainRatioTransitionPhase_{1.f};
  float mainRatioTransitionIncrement_{1.f};
  LateReverbFlavour flavour_{DefaultLateReverbFlavour};
  std::size_t controlCountdown_{};
  float controlRoomScale_{1.f};
  float controlModulationDepth_{};
  float controlModulationAmount_{};
  float controlDecaySeconds_{1.f};
  float controlHighT60_{1.f};
  float controlLowT60_{1.f};
  float controlDiffusion_{1.f};
  float controlShimmerAmount_{};

  static float ClampControl(const float value) noexcept {
    return std::clamp(std::isfinite(value) ? value : 0.f, 0.f, 1.f);
  }
  static float Smoothstep(const float value) noexcept {
    const float x = ClampControl(value);
    return x * x * (3.f - 2.f * x);
  }

  void UpdateControlState(const LateReverbControls &controls) noexcept {
    const float meanFreeTime =
        MeanFreeTimeSeconds(controls.roomDimensionsMetres);
    UpdateMeanDelay(meanFreeTime * static_cast<float>(sampleRate_));
    controlRoomScale_ =
        meanFreeTime /
        MeanFreeTimeSeconds(reverb_defaults::RoomDimensionsMetres);
    const float modulation = ClampControl(controls.modulation);
    controlModulationAmount_ = modulation * modulation;
    controlModulationDepth_ = controlModulationAmount_ *
                              MaximumFdnModulationSeconds *
                              static_cast<float>(sampleRate_);
    const float damping = Smoothstep(controls.damping);
    controlDecaySeconds_ =
        0.25f * std::exp(ClampControl(controls.decay) * std::log(32.f));
    controlHighT60_ = controlDecaySeconds_ * std::pow(0.20f, damping);
    controlLowT60_ = controlDecaySeconds_ * std::pow(0.72f, damping);
    controlDiffusion_ = ClampControl(controls.diffusion);
    controlShimmerAmount_ = Smoothstep(controls.shimmer);
  }

  void PollControlState(const LateReverbControls &controls) noexcept {
    if (controlCountdown_ == 0) {
      UpdateControlState(controls);
      controlCountdown_ = ControlUpdateInterval;
    }
    --controlCountdown_;
  }

  void UpdateMeanDelay(const float targetSamples) noexcept {
    if (!mainDelayInitialized_) {
      currentMeanDelaySamples_ = fromMeanDelaySamples_ = toMeanDelaySamples_ =
          pendingMeanDelaySamples_ = targetSamples;
      mainDelayPhase_ = 1.f;
      mainDelayInitialized_ = true;
      return;
    }
    pendingMeanDelaySamples_ = targetSamples;
    if (mainDelayPhase_ < 1.f ||
        std::abs(targetSamples - currentMeanDelaySamples_) < 0.25f)
      return;
    fromMeanDelaySamples_ = currentMeanDelaySamples_;
    toMeanDelaySamples_ = targetSamples;
    mainDelayPhase_ = 0.f;
  }

  float ReadMainDelayForRatio(const std::size_t line, const float delayRatio,
                              const float modulation) const noexcept {
    if (mainDelayPhase_ >= 1.f)
      return mainDelays_[line].Read(delayRatio * currentMeanDelaySamples_ +
                                    modulation);
    const float fade =
        mainDelayPhase_ * mainDelayPhase_ * (3.f - 2.f * mainDelayPhase_);
    const float from =
        mainDelays_[line].Read(delayRatio * fromMeanDelaySamples_ + modulation);
    const float to =
        mainDelays_[line].Read(delayRatio * toMeanDelaySamples_ + modulation);
    return from + fade * (to - from);
  }

  float ReadMainDelay(const std::size_t line,
                      const float modulation) const noexcept {
    if (mainRatioTransitionPhase_ >= 1.f)
      return ReadMainDelayForRatio(line, mainRatioTarget_[line], modulation);
    const float phase = mainRatioTransitionPhase_;
    const float fade = phase * phase * (3.f - 2.f * phase);
    const float from =
        ReadMainDelayForRatio(line, mainRatioFrom_[line], modulation);
    const float to =
        ReadMainDelayForRatio(line, mainRatioTarget_[line], modulation);
    return from + fade * (to - from);
  }

  static float
  MeanFreeTimeSeconds(const std::array<float, 3> &dimensions) noexcept {
    const float length = std::clamp(dimensions[0], 1.f, 40.f);
    const float width = std::clamp(dimensions[1], 1.f, 40.f);
    const float height = std::clamp(dimensions[2], 1.f, 40.f);
    const float volume = length * width * height;
    const float surface =
        2.f * (length * width + length * height + width * height);
    return std::clamp(4.f * volume / (surface * SpeedOfSound), 0.003f, 0.078f);
  }

  static std::array<float, FdnLineCount>
  ProjectWalls(const std::array<float, WallCount> &walls) noexcept {
    std::array<float, FdnLineCount> frame{};
    for (std::size_t line = 0; line < FdnLineCount; ++line)
      for (std::size_t wall = 0; wall < WallCount; ++wall)
        frame[line] += late_reverb_coefficients::InverseRootLineCount *
                       late_reverb_coefficients::WallProjection[line][wall] *
                       walls[wall];
    return frame;
  }

  LineFrame ProcessFdnLoop(const LineFrame &injection) noexcept {
    LineFrame delayed{};
    for (std::size_t line = 0; line < FdnLineCount; ++line)
      delayed[line] = ReadMainDelay(
          line, controlModulationDepth_ * modulators_[line].Next());

    const float attenuationMeanSamples =
        mainDelayPhase_ < 1.f
            ? std::max(fromMeanDelaySamples_, toMeanDelaySamples_)
            : currentMeanDelaySamples_;
    std::array<float, FdnLineCount> attenuated{};
    for (std::size_t line = 0; line < FdnLineCount; ++line) {
      const float attenuationRatio =
          mainRatioTransitionPhase_ < 1.f
              ? std::max(mainRatioFrom_[line], mainRatioTarget_[line])
              : mainRatioTarget_[line];
      const float pathSeconds = attenuationRatio * attenuationMeanSamples /
                                static_cast<float>(sampleRate_);
      attenuated[line] = mainDecayFilters_[line].Process(
          delayed[line], pathSeconds, controlLowT60_, controlDecaySeconds_,
          controlHighT60_);
    }

    // The octave branch belongs inside the production velvet loop. Projecting
    // a bounded return onto four orthonormal coordinates and then running the
    // complete 16-line VFM avoids a direct grainy output layer and a separate
    // recursive tank. At zero gain the shifters advance their filters, input
    // history and grain state but omit window evaluation and interpolated
    // reads, so automation remains immediate without paying the full render
    // cost for an inaudible branch.
    std::array<float, late_reverb_coefficients::ShimmerBusCount> shimmerBus{};
    std::array<float, late_reverb_coefficients::ShimmerBusCount> shiftedBus{};
    const bool renderShimmer = controlShimmerAmount_ > 1.e-6f;
    for (std::size_t bus = 0; bus < shimmerBus.size(); ++bus) {
      for (std::size_t line = 0; line < FdnLineCount; ++line)
        shimmerBus[bus] +=
            late_reverb_coefficients::ShimmerProjection[bus][line] *
            attenuated[line];
      shimmerHighpassState_[bus] +=
          shimmerHighpassAlpha_ *
          (shimmerBus[bus] - shimmerHighpassState_[bus]);
      const float highpassed = shimmerBus[bus] - shimmerHighpassState_[bus];
      float shifted = 0.f;
      if (renderShimmer)
        shifted = shimmerShifters_[bus].Process(highpassed);
      else
        shimmerShifters_[bus].Advance(highpassed);
      for (auto &stage : shimmerLowpassState_) {
        stage[bus] += shimmerLowpassAlpha_ * (shifted - stage[bus]);
        shifted = stage[bus];
      }
      shiftedBus[bus] = std::isfinite(shifted) ? shifted : 0.f;
    }

    const float shimmerGain =
        MaximumShimmerLoopGain * controlShimmerAmount_;
    if (shimmerGain > 1.e-6f) {
      for (std::size_t line = 0; line < FdnLineCount; ++line)
        for (std::size_t bus = 0; bus < shimmerBus.size(); ++bus)
          attenuated[line] +=
              late_reverb_coefficients::ShimmerProjection[bus][line] *
              (shimmerGain * shiftedBus[bus]);
    }

    auto feedback = feedbackMatrix_.Process(
        attenuated, controlDiffusion_, controlDecaySeconds_, controlHighT60_,
        controlLowT60_, controlRoomScale_, controlModulationAmount_);

    for (std::size_t line = 0; line < FdnLineCount; ++line)
      mainDelays_[line].Push(feedback[line] + 0.42f * injection[line]);

    if (mainDelayPhase_ < 1.f) {
      mainDelayPhase_ = std::min(1.f, mainDelayPhase_ + mainDelayIncrement_);
      if (mainDelayPhase_ >= 1.f) {
        currentMeanDelaySamples_ = toMeanDelaySamples_;
        if (std::abs(pendingMeanDelaySamples_ - currentMeanDelaySamples_) >=
            0.25f) {
          fromMeanDelaySamples_ = currentMeanDelaySamples_;
          toMeanDelaySamples_ = pendingMeanDelaySamples_;
          mainDelayPhase_ = 0.f;
        }
      }
    }
    if (mainRatioTransitionPhase_ < 1.f) {
      mainRatioTransitionPhase_ =
          std::min(1.f, mainRatioTransitionPhase_ +
                            mainRatioTransitionIncrement_);
      if (mainRatioTransitionPhase_ >= 1.f)
        mainRatioFrom_ = mainRatioTarget_;
    }
    for (const auto value : delayed)
      if (!std::isfinite(value)) {
        Reset();
        return {};
      }
    return delayed;
  }

  WallFrame ProcessWallLoop(const WallFrame &inputWalls) noexcept {
    const auto delayed = ProcessFdnLoop(ProjectWalls(inputWalls));
    WallFrame outputWalls{};
    for (std::size_t wall = 0; wall < WallCount; ++wall)
      for (std::size_t line = 0; line < FdnLineCount; ++line)
        outputWalls[wall] +=
            late_reverb_coefficients::InverseRootLineCount *
            late_reverb_coefficients::WallProjection[line][wall] *
            delayed[line];
    return outputWalls;
  }

public:
  LateReverb() { SetSampleRate(sampleRate_); }

  void SetSampleRate(const double sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument("late reverb sample rate must be positive");
    sampleRate_ = sampleRate;
    const auto &baseRatios = LateMainRatios(LateReverbFlavour::Base);
    const auto &optimizedRatios = LateMainRatios(LateReverbFlavour::Optimized);
    const float maximumRatio =
        std::max(*std::max_element(baseRatios.begin(), baseRatios.end()),
                 *std::max_element(optimizedRatios.begin(),
                                   optimizedRatios.end()));
    const auto mainCapacity = static_cast<std::size_t>(std::ceil(
        (0.078 * maximumRatio + MaximumFdnModulationSeconds) * sampleRate_ +
        8.0));
    for (std::size_t line = 0; line < FdnLineCount; ++line) {
      mainDelays_[line].Prepare(mainCapacity);
      mainDecayFilters_[line].Prepare(sampleRate_);
    }
    feedbackMatrix_.Prepare(sampleRate_);
    for (std::size_t line = 0; line < FdnLineCount; ++line)
      modulators_[line].Prepare(
          sampleRate_, late_reverb_coefficients::ModulationRateHz[line],
          0x51f15e5du + static_cast<std::uint32_t>(line) * 0x9e3779b9u);
    constexpr std::array<std::uint32_t,
                         late_reverb_coefficients::ShimmerBusCount>
        ShimmerSeeds{0x8f3a21d7u, 0xc4b82e39u, 0x71d934abu, 0xe29c6f15u};
    for (std::size_t bus = 0; bus < shimmerShifters_.size(); ++bus)
      shimmerShifters_[bus].Prepare(
          sampleRate_, 0.120f,
          static_cast<float>(bus) / static_cast<float>(shimmerShifters_.size()),
          ShimmerSeeds[bus]);
    shimmerHighpassAlpha_ =
        1.f - std::exp(-2.f * Pi * 250.f / static_cast<float>(sampleRate_));
    const float shimmerLowpassCutoff =
        std::min(6'500.f, 0.20f * static_cast<float>(sampleRate_));
    shimmerLowpassAlpha_ = 1.f -
                           std::exp(-2.f * Pi * shimmerLowpassCutoff /
                                    static_cast<float>(sampleRate_));
    mainDelayIncrement_ =
        1.f / std::max(1.f, MainDelayTransitionSeconds *
                                static_cast<float>(sampleRate_));
    mainRatioTransitionIncrement_ = mainDelayIncrement_;
    Reset();
  }

  void Reset() noexcept {
    for (auto &delay : mainDelays_)
      delay.Reset();
    feedbackMatrix_.Reset();
    for (auto &modulator : modulators_)
      modulator.Reset();
    for (auto &filter : mainDecayFilters_)
      filter.Reset();
    for (auto &shifter : shimmerShifters_)
      shifter.Reset();
    shimmerHighpassState_.fill(0.f);
    for (auto &stage : shimmerLowpassState_)
      stage.fill(0.f);
    currentMeanDelaySamples_ = fromMeanDelaySamples_ = toMeanDelaySamples_ =
        pendingMeanDelaySamples_ = static_cast<float>(0.010 * sampleRate_);
    mainDelayPhase_ = 1.f;
    mainDelayInitialized_ = false;
    mainRatioFrom_ = mainRatioTarget_ = LateMainRatios(flavour_);
    mainRatioTransitionPhase_ = 1.f;
    controlCountdown_ = 0;
    controlRoomScale_ = 1.f;
    controlModulationDepth_ = 0.f;
    controlModulationAmount_ = 0.f;
    controlDecaySeconds_ = controlHighT60_ = controlLowT60_ = 1.f;
    controlDiffusion_ = 1.f;
    controlShimmerAmount_ = 0.f;
  }

  double SampleRate() const noexcept { return sampleRate_; }

  void SetFlavour(const LateReverbFlavour flavour) noexcept {
    if (flavour == flavour_)
      return;
    if (mainRatioTransitionPhase_ >= 0.5f)
      mainRatioFrom_ = mainRatioTarget_;
    mainRatioTarget_ = LateMainRatios(flavour);
    mainRatioTransitionPhase_ = 0.f;
    flavour_ = flavour;
    feedbackMatrix_.SetFlavour(flavour);
  }

  // Select a coefficient set without a transition. This is intended for
  // restoring a saved flavour before the first audio sample; live changes use
  // SetFlavour() so the stored tail is crossfaded instead of reinterpreted.
  void SetFlavourImmediate(const LateReverbFlavour flavour) noexcept {
    flavour_ = flavour;
    mainRatioFrom_ = mainRatioTarget_ = LateMainRatios(flavour);
    mainRatioTransitionPhase_ = 1.f;
    feedbackMatrix_.SetFlavourImmediate(flavour);
  }

  LateReverbFlavour Flavour() const noexcept { return flavour_; }

  float IntrinsicFirstOutputSeconds(
      const std::array<float, 3> &roomDimensionsMetres) const noexcept {
    const auto &ratios = LateMainRatios(flavour_);
    return MeanFreeTimeSeconds(roomDimensionsMetres) *
           *std::min_element(ratios.begin(), ratios.end());
  }

  // Static wall-domain probe retained for differentiable-reference parity.
  // Production stereo rendering uses the fixed line-domain vectors below.
  WallFrame ProcessWallFrame(const WallFrame &inputWalls,
                             const LateReverbControls &controls) noexcept {
    PollControlState(controls);
    return ProcessWallLoop(inputWalls);
  }

  StereoFrame Process(const float input,
                      const LateReverbControls &controls) noexcept {
    PollControlState(controls);
    LineFrame injection{};
    const float safeInput = std::isfinite(input) ? input : 0.f;
    for (std::size_t line = 0; line < FdnLineCount; ++line)
      injection[line] =
          safeInput * late_reverb_coefficients::FixedInputVector[line];
    const auto delayed = ProcessFdnLoop(injection);
    StereoFrame stereo{};
    for (std::size_t channel = 0; channel < stereo.size(); ++channel)
      for (std::size_t line = 0; line < FdnLineCount; ++line)
        stereo[channel] +=
            late_reverb_coefficients::FixedStereoOutput[channel][line] *
            delayed[line];
    for (auto &value : stereo)
      value *= late_reverb_coefficients::TailOutputCalibration;

    if (!std::isfinite(stereo[0]) || !std::isfinite(stereo[1])) {
      Reset();
      return {};
    }
    return stereo;
  }
};

} // namespace tfdsp
