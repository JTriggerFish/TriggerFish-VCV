#pragma once

#include "deterministic_random.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace tfdsp::percussion {

struct MicroContactBurstParameters {
  float durationSeconds{.012f};
  float densityHz{8000.f};
  float microDecaySeconds{.0005f};
  float brightness{.7f};
  float amplitude{1.f};
  std::uint32_t seed{1};
};

// Bounded stochastic micro-contact cluster. Random impacts maintain a short
// envelope over filtered noise, so dense settings fuse instead of exposing
// independent grains or adding an unrelated noise layer.
class MicroContactBurst {
public:
  void Prepare(const float sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("micro-contact sample rate must be positive");
    sampleRate_ = sampleRate;
    highpassCoefficient_ =
        std::exp(-6.283185307179586f * 40.f / sampleRate_);
    Reset();
  }

  void Reset() noexcept {
    sample_ = sampleCount_ = 0;
    contactEnvelope_ = lowpassState_ = previousLowpass_ = 0.f;
    highpassState_ = 0.f;
    windowSine_ = 0.f;
    windowCosine_ = 1.f;
  }

  void Trigger(const MicroContactBurstParameters &parameters) noexcept {
    Reset();
    const float duration = FiniteOr(parameters.durationSeconds, 0.f);
    sampleCount_ = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::lround(
               std::clamp(duration, 0.f, 10.f) * sampleRate_)));
    eventProbability_ = std::clamp(
        FiniteOr(parameters.densityHz, 0.f) / sampleRate_, 0.f, .5f);
    const float decaySeconds = std::clamp(
        FiniteOr(parameters.microDecaySeconds, 0.f), 1.f / sampleRate_, .1f);
    contactDecay_ = std::exp(-1.f / (decaySeconds * sampleRate_));
    const float brightness = std::clamp(FiniteOr(parameters.brightness, 0.f), 0.f, 1.f);
    const float cutoffHz = 800.f * std::pow(.45f * sampleRate_ / 800.f, brightness);
    lowpassCoefficient_ = std::exp(-6.283185307179586f * cutoffHz / sampleRate_);
    const auto whiteNoiseVariance = [](const float pole) noexcept {
      const float feed = 1.f - pole;
      return feed / (2.f - feed);
    };
    constexpr float ReferenceBrightness = .7f;
    const float referenceCutoff = 800.f *
        std::pow(.45f * sampleRate_ / 800.f, ReferenceBrightness);
    const float referencePole =
        std::exp(-6.283185307179586f * referenceCutoff / sampleRate_);
    colourGain_ = std::sqrt(
        whiteNoiseVariance(referencePole) /
        std::max(whiteNoiseVariance(lowpassCoefficient_), 1.e-12f));
    amplitude_ = std::clamp(
        tfdsp::FiniteNormalOrZero(FiniteOr(parameters.amplitude, 0.f)),
        0.f, 16.f);
    if (amplitude_ == 0.f) {
      sampleCount_ = 0;
      return;
    }
    random_.Seed(parameters.seed);
    constexpr double Pi = 3.1415926535897932384626433832795;
    const double windowStep = Pi / (static_cast<double>(sampleCount_) + 1.0);
    windowRotationSine_ = static_cast<float>(std::sin(windowStep));
    windowRotationCosine_ = static_cast<float>(std::cos(windowStep));
    windowSine_ = windowRotationSine_;
    windowCosine_ = windowRotationCosine_;
  }

  float Process() noexcept {
    if (sample_ >= sampleCount_)
      return 0.f;
    if (random_.Uniform() < eventProbability_)
      contactEnvelope_ = std::max(contactEnvelope_, .5f + .5f * random_.Uniform());
    else
      contactEnvelope_ *= contactDecay_;
    contactEnvelope_ = tfdsp::FiniteNormalOrZero(contactEnvelope_);

    const float noise = random_.Bipolar();
    lowpassState_ = (1.f - lowpassCoefficient_) * noise +
                    lowpassCoefficient_ * lowpassState_;
    highpassState_ = highpassCoefficient_ *
        (highpassState_ + lowpassState_ - previousLowpass_);
    lowpassState_ = tfdsp::FiniteNormalOrZero(lowpassState_);
    highpassState_ = tfdsp::FiniteNormalOrZero(highpassState_);
    previousLowpass_ = lowpassState_;

    const float output = .5f * colourGain_ * amplitude_ * windowSine_ *
                         contactEnvelope_ * highpassState_;
    const float nextWindowSine = windowSine_ * windowRotationCosine_ +
                                 windowCosine_ * windowRotationSine_;
    windowCosine_ = windowCosine_ * windowRotationCosine_ -
                    windowSine_ * windowRotationSine_;
    windowSine_ = nextWindowSine;
    ++sample_;
    return tfdsp::FiniteNormalOrZero(output);
  }

  bool Active() const noexcept { return sample_ < sampleCount_; }

private:
  static float FiniteOr(const float value, const float fallback) noexcept {
    return std::isfinite(value) ? value : fallback;
  }

  DeterministicRandom random_{};
  float sampleRate_{48000.f};
  float eventProbability_{};
  float contactDecay_{};
  float lowpassCoefficient_{};
  float colourGain_{1.f};
  float highpassCoefficient_{};
  float contactEnvelope_{};
  float lowpassState_{};
  float previousLowpass_{};
  float highpassState_{};
  float amplitude_{1.f};
  float windowSine_{};
  float windowCosine_{1.f};
  float windowRotationSine_{};
  float windowRotationCosine_{1.f};
  std::size_t sample_{};
  std::size_t sampleCount_{};
};

} // namespace tfdsp::percussion
