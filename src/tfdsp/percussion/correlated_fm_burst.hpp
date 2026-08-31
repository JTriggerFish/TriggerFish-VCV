#pragma once

#include "breakpoint_trajectory.hpp"
#include "deterministic_random.hpp"
#include "fixed_oversampling.hpp"
#include "radiation_filter.hpp"
#include "tfdsp/finite_audio.hpp"
#include "tfdsp/sampleRate.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>

namespace tfdsp::percussion {

constexpr std::size_t CorrelatedFmMaximumSegments = 8;

struct CorrelatedFmTrajectory {
  float initialValue{};
  std::array<TrajectorySegment, CorrelatedFmMaximumSegments> segments{};
  std::size_t segmentCount{};
};

struct CorrelatedFmBurstParameters {
  CorrelatedFmTrajectory amplitude{};
  CorrelatedFmTrajectory carrierFrequencyHz{};
  CorrelatedFmTrajectory frequencyDeviationHz{};
  float irregularCutoffHz{8000.f};
  float periodicModulatorHz{2000.f};
  float periodicMix{};
  float carrierPhaseCycles{};
  float modulatorPhaseCycles{};
  float pitchPerturbationCents{};
  float deviationPerturbation{};
  std::uint32_t seed{1};
  bool radiationEnabled{};
  RadiationFilterParameters radiation{};
};

// A finite sine body whose pitch and noisy FM depth share the same event.
// The complete generator runs above host rate and is anti-aliased on return;
// it can represent a compact drum body or supplement another resonator.
template <typename ResamplerType> class OversampledCorrelatedFmBurst {
public:
  static constexpr int OversamplingFactor = ResamplerType::ResamplingFactor;
  using Trajectory = BreakpointTrajectory<CorrelatedFmMaximumSegments>;

  OversampledCorrelatedFmBurst()
      : decimator_(FixedResamplerFactory<ResamplerType>::Create()) {}

  void Prepare(const float sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("FM burst sample rate must be positive");
    sampleRate_ = sampleRate;
    internalRate_ = sampleRate * OversamplingFactor;
    amplitude_.Prepare(internalRate_);
    carrierFrequency_.Prepare(internalRate_);
    frequencyDeviation_.Prepare(internalRate_);
    radiation_.Prepare(sampleRate, {});
    Reset();
  }

  void Reset() noexcept {
    amplitude_.Reset();
    carrierFrequency_.Reset();
    frequencyDeviation_.Reset();
    decimator_->Reset();
    radiation_.Reset();
    irregularState_ = carrierPhase_ = modulatorPhase_ = 0.f;
    drainFrames_ = 0;
    sourceWasActive_ = false;
  }

  void Trigger(const CorrelatedFmBurstParameters &parameters) noexcept {
    Reset();
    random_.Seed(parameters.seed);
    const float pitchRandom = random_.Bipolar();
    const float deviationRandom = random_.Bipolar();
    pitchScale_ = std::exp2(std::clamp(
        FiniteOr(parameters.pitchPerturbationCents, 0.f), -1200.f, 1200.f) *
        pitchRandom / 1200.f);
    deviationScale_ = std::max(0.f, 1.f + std::clamp(
        FiniteOr(parameters.deviationPerturbation, 0.f), 0.f, 1.f) *
        deviationRandom);
    irregularCoefficient_ = 1.f - std::exp(
        -6.283185307179586f *
        std::clamp(FiniteOr(parameters.irregularCutoffHz, 0.f), 0.f,
                   .45f * internalRate_) /
        internalRate_);
    periodicFrequencyHz_ = std::clamp(
        FiniteOr(parameters.periodicModulatorHz, 0.f), 0.f,
        .45f * internalRate_);
    periodicMix_ = std::clamp(FiniteOr(parameters.periodicMix, 0.f), 0.f, 1.f);
    carrierPhase_ = WrappedPhase(parameters.carrierPhaseCycles);
    modulatorPhase_ = WrappedPhase(parameters.modulatorPhaseCycles);
    radiationEnabled_ = parameters.radiationEnabled;
    radiation_.SetStaticParameters(parameters.radiation);
    Start(amplitude_, parameters.amplitude);
    Start(carrierFrequency_, parameters.carrierFrequencyHz);
    Start(frequencyDeviation_, parameters.frequencyDeviationHz);
    sourceWasActive_ = amplitude_.Active();
  }

  float Process() noexcept {
    if (!Active())
      return 0.f;
    Eigen::Array<double, OversamplingFactor, 1> frame;
    for (int phase = 0; phase < OversamplingFactor; ++phase)
      frame(phase) = amplitude_.Active() ? GenerateInternal() : 0.0;
    if (sourceWasActive_ && !amplitude_.Active()) {
      sourceWasActive_ = false;
      drainFrames_ = DecimatorDrainFrames;
    } else if (!sourceWasActive_ && drainFrames_ > 0) {
      --drainFrames_;
    }
    float output = static_cast<float>(decimator_->Downsample(frame));
    output = tfdsp::FiniteNormalOrZero(output);
    return radiationEnabled_ ? radiation_.Process(output) : output;
  }

  bool Active() const noexcept {
    return amplitude_.Active() || drainFrames_ > 0;
  }

private:
  static constexpr std::size_t DecimatorDrainFrames = 32;
  static constexpr float TwoPi = 6.2831853071795864769f;

  static float FiniteOr(const float value, const float fallback) noexcept {
    return std::isfinite(value) ? value : fallback;
  }

  static float WrappedPhase(const float cycles) noexcept {
    const float finite = FiniteOr(cycles, 0.f);
    return TwoPi * (finite - std::floor(finite));
  }

  static void Start(Trajectory &trajectory,
                    const CorrelatedFmTrajectory &parameters) noexcept {
    trajectory.Start(parameters.initialValue, parameters.segments,
                     parameters.segmentCount);
  }

  float GenerateInternal() noexcept {
    const float amplitude = amplitude_.Process();
    const float carrierHz = pitchScale_ * carrierFrequency_.Process();
    const float deviationHz = deviationScale_ * frequencyDeviation_.Process();
    irregularState_ += irregularCoefficient_ *
                       (random_.Bipolar() - irregularState_);
    irregularState_ = tfdsp::FiniteNormalOrZero(irregularState_);
    const float periodic = std::sin(modulatorPhase_);
    const float modulator =
        (1.f - periodicMix_) * irregularState_ + periodicMix_ * periodic;
    const float instantaneousHz = std::clamp(
        carrierHz + deviationHz * modulator, -.45f * internalRate_,
        .45f * internalRate_);
    const float output = amplitude * std::sin(carrierPhase_);
    AdvancePhase(carrierPhase_, instantaneousHz);
    AdvancePhase(modulatorPhase_, periodicFrequencyHz_);
    return tfdsp::FiniteNormalOrZero(output);
  }

  void AdvancePhase(float &phase, const float frequencyHz) const noexcept {
    phase += TwoPi * frequencyHz / internalRate_;
    if (phase >= TwoPi || phase < 0.f)
      phase -= TwoPi * std::floor(phase / TwoPi);
  }

  std::unique_ptr<ResamplerType> decimator_;
  Trajectory amplitude_{};
  Trajectory carrierFrequency_{};
  Trajectory frequencyDeviation_{};
  DeterministicRandom random_{};
  RadiationFilter radiation_{};
  float sampleRate_{48000.f};
  float internalRate_{96000.f};
  float pitchScale_{1.f};
  float deviationScale_{1.f};
  float irregularCoefficient_{};
  float irregularState_{};
  float periodicFrequencyHz_{};
  float periodicMix_{};
  float carrierPhase_{};
  float modulatorPhase_{};
  std::size_t drainFrames_{};
  bool sourceWasActive_{};
  bool radiationEnabled_{};
};

using CorrelatedFmBurst =
    OversampledCorrelatedFmBurst<tfdsp::X2Resampler_Order7>;
using CorrelatedFmBurstReference4x =
    OversampledCorrelatedFmBurst<tfdsp::X4Resampler_Order7>;

} // namespace tfdsp::percussion
