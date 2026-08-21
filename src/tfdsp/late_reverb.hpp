#pragma once

#include "late_reverb_coefficient_sets.hpp"
#include "multiband_decay_filter.hpp"
#include "reverb_defaults.hpp"
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
  float exciterPosition{reverb_defaults::Source[0]};
  std::array<float, 3> listener{reverb_defaults::Listener};
  std::array<float, 3> roomDimensionsMetres{
      reverb_defaults::RoomDimensionsMetres};
};

class LateReverb {
public:
  static constexpr std::size_t FdnLineCount =
      late_reverb_coefficients::LineCount;
  static constexpr std::size_t WallCount = late_reverb_coefficients::WallCount;
  static constexpr std::size_t MaximumSources = 8;
  using StereoFrame = std::array<float, 2>;
  using WallFrame = std::array<float, WallCount>;
  using SourceFrame = std::array<float, MaximumSources>;
  using SourcePosition = std::array<float, 3>;
  using SourcePositions = std::array<SourcePosition, MaximumSources>;

private:
  static constexpr float Pi = 3.14159265358979323846f;
  static constexpr float SpeedOfSound = 343.f;
  static constexpr float MaximumFdnModulationSeconds = 0.000250f;
  // Shimmer adds a bounded octave return to four orthogonal coordinates inside
  // the main velvet loop without replacing the unshifted room feedback. The
  // remaining twelve coordinates are unaltered.
  static constexpr float MaximumShimmerLoopGain = 0.85f;
  static constexpr float MainDelayTransitionSeconds = 0.050f;
  static constexpr float GeometryTransitionSeconds = 0.100f;
  static constexpr std::size_t GeometryUpdateInterval = 64;
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
        const auto index =
            (writeIndex_ + buffer_.size() - distance) % buffer_.size();
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

  // Fixed read heads crossfade when geometry changes. This changes spatial
  // delays without the pitch sweep produced by moving a read head.
  class TransitionDelay {
  public:
    void Prepare(const std::size_t capacity, const float increment) {
      delay_.Prepare(capacity);
      increment_ = increment;
      Reset();
    }
    void Reset() noexcept {
      delay_.Reset();
      current_ = from_ = to_ = pending_ = 2.f;
      phase_ = 1.f;
      initialized_ = false;
    }
    void SetTarget(const float samples) noexcept {
      const float target = std::max(samples, 2.f);
      if (!initialized_) {
        current_ = from_ = to_ = pending_ = target;
        initialized_ = true;
        return;
      }
      pending_ = target;
      if (phase_ < 1.f || std::abs(pending_ - current_) < 0.25f)
        return;
      from_ = current_;
      to_ = pending_;
      phase_ = 0.f;
    }
    float Process(const float input) noexcept {
      float output = 0.f;
      if (phase_ >= 1.f) {
        output = delay_.Read(current_);
      } else {
        const float fade = phase_ * phase_ * (3.f - 2.f * phase_);
        output =
            delay_.Read(from_) + fade * (delay_.Read(to_) - delay_.Read(from_));
        phase_ = std::min(1.f, phase_ + increment_);
        if (phase_ >= 1.f) {
          current_ = to_;
          if (std::abs(pending_ - current_) >= 0.25f) {
            from_ = current_;
            to_ = pending_;
            phase_ = 0.f;
          }
        }
      }
      delay_.Push(input);
      return output;
    }

  private:
    FractionalDelay delay_{};
    float current_{2.f};
    float from_{2.f};
    float to_{2.f};
    float pending_{2.f};
    float phase_{1.f};
    float increment_{1.f};
    bool initialized_{};
  };

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

  double sampleRate_{48'000.0};
  std::array<FractionalDelay, FdnLineCount> mainDelays_{};
  VelvetFeedbackMatrix feedbackMatrix_{};
  std::array<SmoothRandomModulator, FdnLineCount> modulators_{};
  std::array<MultibandDecayFilter, FdnLineCount> mainDecayFilters_{};
  std::array<std::array<TransitionDelay, WallCount>, MaximumSources>
      sourceWallDelays_{};
  std::array<TransitionDelay, WallCount> listenerWallDelays_{};
  std::array<std::array<float, WallCount>, MaximumSources> sourceWallWeights_{};
  std::array<float, WallCount> listenerWallWeights_{};
  std::array<std::array<float, WallCount>, 2> listenerDecoder_{};
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
  std::size_t geometryCountdown_{};
  std::size_t controlCountdown_{};
  float controlRoomScale_{1.f};
  float controlModulationDepth_{};
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
    const float modulationRange =
        std::max(0.f, (modulation - 0.35f) / 0.65f);
    controlModulationDepth_ =
        modulationRange * modulationRange * modulationRange *
        MaximumFdnModulationSeconds * static_cast<float>(sampleRate_);
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

  static std::array<float, WallCount>
  WallDistances(const SourcePosition &normalized,
                const std::array<float, 3> &dimensions) noexcept {
    const float x = ClampControl(normalized[0]) * dimensions[0];
    const float y = ClampControl(normalized[1]) * dimensions[1];
    const float z = ClampControl(normalized[2]) * dimensions[2];
    return {x, dimensions[0] - x, y, dimensions[1] - y, z, dimensions[2] - z};
  }

  static std::array<float, WallCount>
  WallWeights(const std::array<float, WallCount> &distance,
              const std::array<float, 3> &dimensions) noexcept {
    const std::array<float, WallCount> area{
        dimensions[1] * dimensions[2], dimensions[1] * dimensions[2],
        dimensions[0] * dimensions[2], dimensions[0] * dimensions[2],
        dimensions[0] * dimensions[1], dimensions[0] * dimensions[1]};
    std::array<float, WallCount> result{};
    float norm = 0.f;
    for (std::size_t wall = 0; wall < WallCount; ++wall) {
      const float solidAngleProxy =
          area[wall] /
          (distance[wall] * distance[wall] + area[wall] / Pi + 1.e-4f);
      result[wall] = std::sqrt(std::max(solidAngleProxy, 1.e-6f));
      norm += result[wall] * result[wall];
    }
    const float inverseNorm = 1.f / std::sqrt(std::max(norm, 1.e-12f));
    for (auto &value : result)
      value *= inverseNorm;
    return result;
  }

  void UpdateGeometry(const SourcePositions &positions,
                      const std::size_t sourceCount,
                      const LateReverbControls &controls) noexcept {
    std::array<float, 3> dimensions{};
    for (std::size_t axis = 0; axis < 3; ++axis)
      dimensions[axis] =
          std::clamp(std::isfinite(controls.roomDimensionsMetres[axis])
                         ? controls.roomDimensionsMetres[axis]
                         : 3.f,
                     1.f, 40.f);
    const float samplesPerMetre =
        static_cast<float>(sampleRate_) / SpeedOfSound;
    const std::size_t active = std::min(sourceCount, MaximumSources);
    for (std::size_t source = 0; source < active; ++source) {
      const auto distances = WallDistances(positions[source], dimensions);
      sourceWallWeights_[source] = WallWeights(distances, dimensions);
      const float nearest =
          *std::min_element(distances.begin(), distances.end());
      for (std::size_t wall = 0; wall < WallCount; ++wall)
        sourceWallDelays_[source][wall].SetTarget(
            2.f + (distances[wall] - nearest) * samplesPerMetre);
    }
    const auto listenerDistances = WallDistances(controls.listener, dimensions);
    listenerWallWeights_ = WallWeights(listenerDistances, dimensions);
    const float nearest =
        *std::min_element(listenerDistances.begin(), listenerDistances.end());
    for (std::size_t wall = 0; wall < WallCount; ++wall)
      listenerWallDelays_[wall].SetTarget(
          2.f + (listenerDistances[wall] - nearest) * samplesPerMetre);

    // Deterministic equal-energy stereo factorization. Mid follows the six
    // physical listener-wall weights. Side starts as a directional wall mode,
    // is projected orthogonal to Mid, then normalized. Consequently the L/R
    // decoder rows are orthonormal and exactly equal energy for every scene.
    constexpr std::array<float, WallCount> SidePattern{1.f,    -1.f,  0.35f,
                                                       -0.35f, 0.15f, -0.15f};
    std::array<float, WallCount> side{};
    float projection = 0.f;
    for (std::size_t wall = 0; wall < WallCount; ++wall) {
      side[wall] = listenerWallWeights_[wall] * SidePattern[wall];
      projection += side[wall] * listenerWallWeights_[wall];
    }
    float sideNormSquared = 0.f;
    for (std::size_t wall = 0; wall < WallCount; ++wall) {
      side[wall] -= projection * listenerWallWeights_[wall];
      sideNormSquared += side[wall] * side[wall];
    }
    if (sideNormSquared < 1.e-8f) {
      side = {1.f, -1.f, 0.f, 0.f, 0.f, 0.f};
      projection = 0.f;
      for (std::size_t wall = 0; wall < WallCount; ++wall)
        projection += side[wall] * listenerWallWeights_[wall];
      sideNormSquared = 0.f;
      for (std::size_t wall = 0; wall < WallCount; ++wall) {
        side[wall] -= projection * listenerWallWeights_[wall];
        sideNormSquared += side[wall] * side[wall];
      }
    }
    const float inverseSideNorm =
        1.f / std::sqrt(std::max(sideNormSquared, 1.e-12f));
    constexpr float InverseRootTwo = 0.7071067811865475f;
    for (std::size_t wall = 0; wall < WallCount; ++wall) {
      side[wall] *= inverseSideNorm;
      listenerDecoder_[0][wall] =
          InverseRootTwo * (listenerWallWeights_[wall] + side[wall]);
      listenerDecoder_[1][wall] =
          InverseRootTwo * (listenerWallWeights_[wall] - side[wall]);
    }
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

  WallFrame ProcessWallLoop(const WallFrame &inputWalls) noexcept {
    const auto injection = ProjectWalls(inputWalls);
    std::array<float, FdnLineCount> delayed{};
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
        controlLowT60_, controlRoomScale_);

    for (std::size_t line = 0; line < FdnLineCount; ++line)
      mainDelays_[line].Push(feedback[line] + 0.42f * injection[line]);

    WallFrame outputWalls{};
    for (std::size_t wall = 0; wall < WallCount; ++wall)
      for (std::size_t line = 0; line < FdnLineCount; ++line)
        outputWalls[wall] +=
            late_reverb_coefficients::InverseRootLineCount *
            late_reverb_coefficients::WallProjection[line][wall] *
            delayed[line];

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
    for (const auto value : outputWalls)
      if (!std::isfinite(value)) {
        Reset();
        return {};
      }
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
    const auto geometryCapacity =
        static_cast<std::size_t>(std::ceil(0.120 * sampleRate_)) + 8;
    const float geometryIncrement =
        1.f / std::max(1.f, GeometryTransitionSeconds *
                                static_cast<float>(sampleRate_));
    for (auto &source : sourceWallDelays_)
      for (auto &delay : source)
        delay.Prepare(geometryCapacity, geometryIncrement);
    for (auto &delay : listenerWallDelays_)
      delay.Prepare(geometryCapacity, geometryIncrement);
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
    for (auto &source : sourceWallDelays_)
      for (auto &delay : source)
        delay.Reset();
    for (auto &delay : listenerWallDelays_)
      delay.Reset();
    for (auto &shifter : shimmerShifters_)
      shifter.Reset();
    shimmerHighpassState_.fill(0.f);
    for (auto &stage : shimmerLowpassState_)
      stage.fill(0.f);
    sourceWallWeights_ = {};
    listenerWallWeights_ = {};
    listenerDecoder_ = {};
    currentMeanDelaySamples_ = fromMeanDelaySamples_ = toMeanDelaySamples_ =
        pendingMeanDelaySamples_ = static_cast<float>(0.010 * sampleRate_);
    mainDelayPhase_ = 1.f;
    mainDelayInitialized_ = false;
    mainRatioFrom_ = mainRatioTarget_ = LateMainRatios(flavour_);
    mainRatioTransitionPhase_ = 1.f;
    geometryCountdown_ = 0;
    controlCountdown_ = 0;
    controlRoomScale_ = 1.f;
    controlModulationDepth_ = 0.f;
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

  // Static wall-domain probe used by the differentiable-reference parity
  // test. It runs the exact production delay/loss/VFM loop while bypassing
  // only the position-dependent source and listener connections.
  WallFrame ProcessWallFrame(const WallFrame &inputWalls,
                             const LateReverbControls &controls) noexcept {
    PollControlState(controls);
    return ProcessWallLoop(inputWalls);
  }

  StereoFrame Process(const float input,
                      const LateReverbControls &controls) noexcept {
    SourceFrame sources{};
    sources[0] = input;
    SourcePositions positions{};
    positions[0] = {controls.exciterPosition, reverb_defaults::Source[1],
                    reverb_defaults::Source[2]};
    return Process(sources, positions, 1, controls);
  }

  StereoFrame Process(const SourceFrame &sources,
                      const SourcePositions &positions,
                      const std::size_t sourceCount,
                      const LateReverbControls &controls) noexcept {
    SourceFrame gains{};
    gains.fill(1.f);
    return Process(sources, positions, sourceCount, controls, gains);
  }

  StereoFrame Process(const SourceFrame &sources,
                      const SourcePositions &positions,
                      const std::size_t sourceCount,
                      const LateReverbControls &controls,
                      const SourceFrame &sourceGains) noexcept {
    PollControlState(controls);
    if (geometryCountdown_ == 0) {
      UpdateGeometry(positions, sourceCount, controls);
      geometryCountdown_ = GeometryUpdateInterval;
    }
    --geometryCountdown_;

    WallFrame inputWalls{};
    const std::size_t active = std::min(sourceCount, MaximumSources);
    for (std::size_t source = 0; source < MaximumSources; ++source) {
      const float input = source < active && std::isfinite(sources[source])
                              ? sources[source]
                              : 0.f;
      const float gain = source < active && std::isfinite(sourceGains[source])
                             ? std::max(0.f, sourceGains[source])
                             : 0.f;
      for (std::size_t wall = 0; wall < WallCount; ++wall) {
        const float connected = sourceWallDelays_[source][wall].Process(input);
        if (source < active)
          inputWalls[wall] +=
              gain * sourceWallWeights_[source][wall] * connected;
      }
    }
    const auto outputWalls = ProcessWallLoop(inputWalls);

    StereoFrame stereo{};
    for (std::size_t wall = 0; wall < WallCount; ++wall) {
      const float connected =
          listenerWallDelays_[wall].Process(outputWalls[wall]);
      stereo[0] += listenerDecoder_[0][wall] * connected;
      stereo[1] += listenerDecoder_[1][wall] * connected;
    }
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
