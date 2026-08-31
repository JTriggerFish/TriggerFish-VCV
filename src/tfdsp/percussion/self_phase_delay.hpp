#pragma once

#include "self_phase_delay_core.hpp"
#include "tfdsp/finite_audio.hpp"
#include "tfdsp/sampleRate.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>

namespace tfdsp::percussion {

template <typename ResamplerType> struct SelfPhaseResamplerFactory;

template <> struct SelfPhaseResamplerFactory<tfdsp::DummyResampler> {
  static auto Create() { return tfdsp::CreateDummyResampler(); }
};

template <> struct SelfPhaseResamplerFactory<tfdsp::X2Resampler_Order7> {
  static auto Create() { return tfdsp::CreateX2Resampler_Chebychev7(); }
};

template <> struct SelfPhaseResamplerFactory<tfdsp::X4Resampler_Order7> {
  static auto Create() { return tfdsp::CreateX4Resampler_Cheby7(); }
};

// The production stage uses 2x oversampling. The template remains public so
// tests and offline tools can render the identical nonlinear core at 4x.
template <typename ResamplerType> class OversampledSelfPhaseDelay {
public:
  static constexpr int OversamplingFactor = ResamplerType::ResamplingFactor;

  OversampledSelfPhaseDelay()
      : interpolator_(SelfPhaseResamplerFactory<ResamplerType>::Create()),
        decimator_(SelfPhaseResamplerFactory<ResamplerType>::Create()) {}

  void Prepare(const float sampleRate, const float maximumDelaySamples) {
    core_.Prepare(sampleRate * OversamplingFactor,
                  maximumDelaySamples * OversamplingFactor);
    MeasureResamplingLatency(maximumDelaySamples);
    ConfigureCore(parameters_);
    Reset();
  }

  void Reset() noexcept {
    interpolator_->Reset();
    decimator_->Reset();
    core_.Reset();
  }

  void SetParameters(const SelfPhaseDelayParameters &parameters) noexcept {
    parameters_ = parameters;
    ConfigureCore(parameters);
  }

  std::size_t NominalLatencySamples() const noexcept {
    return resamplingLatencySamples_;
  }

private:
  void ConfigureCore(const SelfPhaseDelayParameters &parameters) noexcept {
    auto scaled = parameters;
    scaled.centreDelaySamples =
        (std::isfinite(parameters.centreDelaySamples)
             ? parameters.centreDelaySamples
             : SelfPhaseDelayParameters{}.centreDelaySamples) *
        OversamplingFactor;
    scaled.maximumExcursionSamples =
        (std::isfinite(parameters.maximumExcursionSamples)
             ? parameters.maximumExcursionSamples
             : 0.f) *
        OversamplingFactor;
    core_.SetParameters(scaled);
  }

  void MeasureResamplingLatency(const float maximumDelaySamples) noexcept {
    SelfPhaseDelayParameters linear;
    linear.centreDelaySamples = std::clamp(12.f, 4.f, maximumDelaySamples);
    linear.maximumExcursionSamples = 0.f;
    linear.drive = 0.f;
    ConfigureCore(linear);
    Reset();
    float peak = 0.f;
    std::size_t peakSample = 0;
    const auto count = static_cast<std::size_t>(maximumDelaySamples) + 64;
    for (std::size_t sample = 0; sample < count; ++sample) {
      const float output = Process(sample == 0 ? 1.f : 0.f);
      if (std::abs(output) > peak) {
        peak = std::abs(output);
        peakSample = sample;
      }
    }
    const auto centre = static_cast<std::size_t>(std::lround(
        linear.centreDelaySamples));
    resamplingLatencySamples_ = peakSample > centre ? peakSample - centre : 0;
  }

public:
  float Process(float input) noexcept {
    input = tfdsp::FiniteNormalOrZero(input);
    const auto upsampled = interpolator_->Upsample(input);
    Eigen::Array<double, OversamplingFactor, 1> processed;
    for (int phase = 0; phase < OversamplingFactor; ++phase)
      processed(phase) = core_.Process(static_cast<float>(upsampled(phase)));
    const float output = static_cast<float>(decimator_->Downsample(processed));
    return tfdsp::FiniteNormalOrZero(output);
  }

  float CentreDelaySamples() const noexcept {
    return core_.CentreDelaySamples() / OversamplingFactor;
  }

private:
  std::unique_ptr<ResamplerType> interpolator_;
  std::unique_ptr<ResamplerType> decimator_;
  SelfPhaseDelayCore core_{};
  SelfPhaseDelayParameters parameters_{};
  std::size_t resamplingLatencySamples_{};
};

using SelfPhaseDelay = OversampledSelfPhaseDelay<tfdsp::X2Resampler_Order7>;
using SelfPhaseDelayReference1x =
    OversampledSelfPhaseDelay<tfdsp::DummyResampler>;
using SelfPhaseDelayReference4x =
    OversampledSelfPhaseDelay<tfdsp::X4Resampler_Order7>;

} // namespace tfdsp::percussion
