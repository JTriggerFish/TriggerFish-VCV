#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace tfdsp::percussion {

struct TonalContactChirpParameters {
  float durationSeconds{.004f};
  float startFrequencyHz{4000.f};
  float endFrequencyHz{1800.f};
  float amplitude{1.f};
  float decayNepers{2.f};
};

// Windowed exponential-frequency chirp for coherent tip articulation. Its
// endpoints are exactly zero, preventing a separate truncation click.
class TonalContactChirp {
public:
  void Prepare(const float sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("contact-chirp sample rate must be positive");
    sampleRate_ = sampleRate;
    Reset();
  }

  void Reset() noexcept {
    sample_ = sampleCount_ = 0;
    phase_ = 0.0;
    windowSine_ = 0.f;
    windowCosine_ = 1.f;
    decay_ = 1.f;
  }

  void Trigger(const TonalContactChirpParameters &parameters) noexcept {
    const float duration = FiniteOr(parameters.durationSeconds, 0.f);
    sampleCount_ = std::max<std::size_t>(
        3, static_cast<std::size_t>(std::lround(
               std::clamp(duration, 0.f, 1.f) * sampleRate_)));
    startHz_ = BoundFrequency(parameters.startFrequencyHz);
    endHz_ = BoundFrequency(parameters.endFrequencyHz);
    amplitude_ = std::clamp(FiniteOr(parameters.amplitude, 0.f), 0.f, 16.f);
    if (amplitude_ == 0.f) {
      sampleCount_ = 0;
      return;
    }
    decayNepers_ = std::clamp(FiniteOr(parameters.decayNepers, 0.f), 0.f, 20.f);
    sample_ = 0;
    phase_ = 0.0;
    constexpr double Pi = 3.1415926535897932384626433832795;
    const double denominator = static_cast<double>(sampleCount_ - 1);
    const double windowStep = Pi / denominator;
    windowRotationSine_ = static_cast<float>(std::sin(windowStep));
    windowRotationCosine_ = static_cast<float>(std::cos(windowStep));
    windowSine_ = 0.f;
    windowCosine_ = 1.f;
    currentHz_ = startHz_;
    frequencyMultiplier_ = static_cast<float>(
        std::pow(endHz_ / startHz_, 1.0 / denominator));
    decay_ = 1.f;
    decayMultiplier_ = static_cast<float>(
        std::exp(-decayNepers_ / denominator));
  }

  float Process() noexcept {
    if (sample_ >= sampleCount_)
      return 0.f;
    constexpr double Pi = 3.1415926535897932384626433832795;
    constexpr double TwoPi = 2.0 * Pi;
    const float output = amplitude_ * windowSine_ * decay_ *
                         static_cast<float>(std::sin(phase_));
    phase_ = std::remainder(phase_ + TwoPi * currentHz_ / sampleRate_, TwoPi);
    currentHz_ *= frequencyMultiplier_;
    decay_ *= decayMultiplier_;
    const float nextWindowSine = windowSine_ * windowRotationCosine_ +
                                 windowCosine_ * windowRotationSine_;
    windowCosine_ = windowCosine_ * windowRotationCosine_ -
                    windowSine_ * windowRotationSine_;
    windowSine_ = nextWindowSine;
    ++sample_;
    return output;
  }

  bool Active() const noexcept { return sample_ < sampleCount_; }

private:
  static float FiniteOr(const float value, const float fallback) noexcept {
    return std::isfinite(value) ? value : fallback;
  }

  float BoundFrequency(const float frequencyHz) const noexcept {
    return std::clamp(FiniteOr(frequencyHz, 1.f), 1.f, .45f * sampleRate_);
  }

  float sampleRate_{48000.f};
  float startHz_{4000.f};
  float endHz_{1800.f};
  float amplitude_{1.f};
  float decayNepers_{2.f};
  float currentHz_{4000.f};
  float frequencyMultiplier_{1.f};
  float decay_{1.f};
  float decayMultiplier_{1.f};
  float windowSine_{};
  float windowCosine_{1.f};
  float windowRotationSine_{};
  float windowRotationCosine_{1.f};
  double phase_{};
  std::size_t sample_{};
  std::size_t sampleCount_{};
};

} // namespace tfdsp::percussion
