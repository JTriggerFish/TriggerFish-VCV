#pragma once

#include "schroeder_allpass.hpp"
#include "self_phase_delay.hpp"
#include "slow_modulated_delay.hpp"
#include "static_fractional_delay.hpp"
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
    selfPhase_.Reset();
    loss_.Reset();
  }

  float Process(float bodyDrive) noexcept {
    if (!std::isfinite(bodyDrive))
      bodyDrive = 0.f;
    float circulating = base_.Read();
    circulating = slow_.Process(circulating);
    circulating = firstAllpass_.Process(circulating);
    circulating = secondAllpass_.Process(circulating);
    circulating = selfPhase_.Process(circulating);
    const float feedback = feedbackGain_ * loss_.Process(circulating);
    base_.Push(bodyDrive + feedback);
    return std::isfinite(circulating) ? circulating : 0.f;
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
    selfPhase_.SetParameters(parameters.selfPhase);
    feedbackGain_ = std::clamp(
        std::isfinite(parameters.feedbackGain) ? parameters.feedbackGain : 0.f,
        0.f, .999f);
    const float pathSamples = base_.DelaySamples() + slow_.CentreDelaySamples() +
        firstAllpass_.DelaySamples() + secondAllpass_.DelaySamples() +
        selfPhase_.CentreDelaySamples();
    loss_.SetDecayTimes(pathSamples / sampleRate_, parameters.decay);
  }

  StaticFractionalDelay base_{};
  SlowModulatedDelay slow_{};
  SchroederAllpass firstAllpass_{};
  SchroederAllpass secondAllpass_{};
  SelfPhaseDelay selfPhase_{};
  ThreeBandDecayFilter loss_{};
  float sampleRate_{48000.f};
  float feedbackGain_{};
};

} // namespace tfdsp::percussion
