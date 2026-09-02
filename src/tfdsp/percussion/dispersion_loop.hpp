#pragma once

#include "schroeder_allpass.hpp"
#include "self_phase_delay.hpp"
#include "slow_modulated_delay.hpp"
#include "static_fractional_delay.hpp"
#include "tfdsp/finite_audio.hpp"
#include "three_band_decay_filter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace tfdsp::percussion {

struct DispersionLoopParameters {
  float baseDelaySamples{16.f};
  float slowDelaySamples{12.f};
  float slowDepthSamples{1.f};
  float slowRateHz{.2f};
  float firstAllpassDelaySamples{9.f};
  float firstAllpassGain{.55f};
  float secondAllpassDelaySamples{13.f};
  float secondAllpassGain{-.45f};
  float thirdAllpassDelaySamples{5.f};
  float thirdAllpassGain{.3f};
  float fourthAllpassDelaySamples{17.f};
  float fourthAllpassGain{-.25f};
  SelfPhaseDelayParameters selfPhase{};
  ThreeBandDecayTimes decay{.8f, .55f, .25f};
  float feedbackGain{.95f};
  float lowCrossoverHz{600.f};
  float highCrossoverHz{6000.f};
  std::uint32_t modulationSeed{1};
};

// Causal serial dispersion processor inside one explicit feedback loop. The
// returned tap is for analysis/resonator drive and is not a direct output mix.
class DispersionLoop {
public:
  void Prepare(const float sampleRate, const float maximumDelaySamples,
               const DispersionLoopParameters &parameters) {
    sampleRate_ = sampleRate;
    base_.Prepare(maximumDelaySamples, parameters.baseDelaySamples);
    slow_.Prepare(sampleRate, maximumDelaySamples, parameters.modulationSeed);
    firstAllpass_.Prepare(maximumDelaySamples,
                          parameters.firstAllpassDelaySamples,
                          parameters.firstAllpassGain);
    secondAllpass_.Prepare(maximumDelaySamples,
                           parameters.secondAllpassDelaySamples,
                           parameters.secondAllpassGain);
    thirdAllpass_.Prepare(maximumDelaySamples,
                          parameters.thirdAllpassDelaySamples,
                          parameters.thirdAllpassGain);
    fourthAllpass_.Prepare(maximumDelaySamples,
                           parameters.fourthAllpassDelaySamples,
                           parameters.fourthAllpassGain);
    selfPhase_.Prepare(sampleRate, maximumDelaySamples);
    loss_.Prepare(sampleRate, parameters.lowCrossoverHz,
                  parameters.highCrossoverHz);
    ConfigurePreparedParameters(parameters);
    Reset();
  }

  void Reset() noexcept {
    base_.Reset();
    slow_.Reset();
    firstAllpass_.Reset();
    secondAllpass_.Reset();
    thirdAllpass_.Reset();
    fourthAllpass_.Reset();
    selfPhase_.Reset();
    loss_.Reset();
  }

  float Process(float bodyDrive) noexcept {
    return Process(bodyDrive, {});
  }

  void SetSelfPhaseParameters(
      const SelfPhaseDelayParameters &parameters) noexcept {
    selfPhase_.SetParameters(parameters);
  }

  float Process(float bodyDrive,
                const PassiveConstraintGains constraint) noexcept {
    bodyDrive = tfdsp::FiniteNormalOrZero(bodyDrive);
    float circulating = base_.Read();
    circulating = slow_.Process(circulating);
    circulating = firstAllpass_.Process(circulating);
    circulating = secondAllpass_.Process(circulating);
    circulating = thirdAllpass_.Process(circulating);
    circulating = fourthAllpass_.Process(circulating);
    circulating = selfPhase_.Process(circulating);
    const float feedback = feedbackGain_ * loss_.Process(circulating, constraint);
    base_.Push(bodyDrive + feedback);
    return tfdsp::FiniteNormalOrZero(circulating);
  }

  float MinimumPropagationSamples() const noexcept {
    return minimumPropagationSamples_;
  }

  float NominalPropagationSamples() const noexcept {
    return nominalPropagationSamples_;
  }

private:
  void ConfigurePreparedParameters(
      const DispersionLoopParameters &parameters) noexcept {
    slow_.SetStaticParameters(parameters.slowDelaySamples,
                              parameters.slowDepthSamples,
                              parameters.slowRateHz, sampleRate_,
                              parameters.modulationSeed);
    firstAllpass_.SetFeedbackGain(parameters.firstAllpassGain);
    secondAllpass_.SetFeedbackGain(parameters.secondAllpassGain);
    thirdAllpass_.SetFeedbackGain(parameters.thirdAllpassGain);
    fourthAllpass_.SetFeedbackGain(parameters.fourthAllpassGain);
    selfPhase_.SetParameters(parameters.selfPhase);
    feedbackGain_ = std::clamp(
        std::isfinite(parameters.feedbackGain) ? parameters.feedbackGain : 0.f,
        0.f, .9999f);
    minimumPropagationSamples_ =
        base_.DelaySamples() + slow_.CentreDelaySamples() +
        firstAllpass_.DelaySamples() + secondAllpass_.DelaySamples() +
        thirdAllpass_.DelaySamples() + fourthAllpass_.DelaySamples() +
        selfPhase_.CentreDelaySamples();
    nominalPropagationSamples_ = minimumPropagationSamples_ +
        static_cast<float>(selfPhase_.NominalLatencySamples());
    loss_.SetDecayTimes(nominalPropagationSamples_ / sampleRate_,
                        parameters.decay);
  }

  StaticFractionalDelay base_{};
  SlowModulatedDelay slow_{};
  SchroederAllpass firstAllpass_{};
  SchroederAllpass secondAllpass_{};
  SchroederAllpass thirdAllpass_{};
  SchroederAllpass fourthAllpass_{};
  SelfPhaseDelay selfPhase_{};
  ThreeBandDecayFilter loss_{};
  float sampleRate_{48000.f};
  float feedbackGain_{};
  float minimumPropagationSamples_{};
  float nominalPropagationSamples_{};
};

} // namespace tfdsp::percussion
