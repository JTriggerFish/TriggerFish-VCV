#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace tfdsp::percussion {

struct QuadratureSample {
  float cosine{};
  float sine{};
};

// Phase-continuous quadrature oscillator. Frequency changes update one complex
// rotation; the audio-rate step itself uses four multiplies and two additions.
class QuadratureOscillator {
public:
  void Prepare(const float sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("oscillator sample rate must be positive");
    sampleRate_ = sampleRate;
    SetFrequencyHz(0.f);
    Reset();
  }

  void Reset() noexcept {
    cosine_ = 1.f;
    sine_ = 0.f;
    samplesUntilNormalization_ = NormalizationPeriod;
  }

  void SetFrequencyHz(float frequencyHz) noexcept {
    if (!std::isfinite(frequencyHz))
      frequencyHz = 0.f;
    frequencyHz = std::clamp(frequencyHz, -.49f * sampleRate_,
                             .49f * sampleRate_);
    if (frequencyHz == frequencyHz_)
      return;
    constexpr double TwoPi = 6.283185307179586476925286766559;
    const double step = TwoPi * frequencyHz / sampleRate_;
    rotationCosine_ = static_cast<float>(std::cos(step));
    rotationSine_ = static_cast<float>(std::sin(step));
    frequencyHz_ = frequencyHz;
  }

  QuadratureSample Process() noexcept {
    const QuadratureSample output{cosine_, sine_};
    const float nextCosine = cosine_ * rotationCosine_ -
                             sine_ * rotationSine_;
    sine_ = sine_ * rotationCosine_ + cosine_ * rotationSine_;
    cosine_ = nextCosine;
    if (--samplesUntilNormalization_ == 0)
      Normalize();
    return output;
  }

private:
  void Normalize() noexcept {
    const float inverseMagnitude = 1.f / std::sqrt(
        std::max(1.e-20f, cosine_ * cosine_ + sine_ * sine_));
    cosine_ *= inverseMagnitude;
    sine_ *= inverseMagnitude;
    samplesUntilNormalization_ = NormalizationPeriod;
  }

  static constexpr std::size_t NormalizationPeriod = 64;
  float sampleRate_{48000.f};
  float frequencyHz_{};
  float cosine_{1.f};
  float sine_{};
  float rotationCosine_{1.f};
  float rotationSine_{};
  std::size_t samplesUntilNormalization_{NormalizationPeriod};
};

} // namespace tfdsp::percussion
