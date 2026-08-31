#pragma once

#include "biquad.hpp"
#include "biquad_design.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace tfdsp::percussion {

// Fourth-order source/output guards for an SSB translation. The source band
// is narrowed so translated positive-frequency content reaches neither DC nor
// Nyquist; the output guard removes transition-band residue.
class TranslationBandLimiter {
public:
  void Prepare(const float sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("translation limiter rate must be positive");
    sampleRate_ = sampleRate;
    coefficientTransitionSamples_ = static_cast<std::size_t>(
        std::max(1.f, .002f * sampleRate_));
    const float guard = GuardHz();
    Configure(outputHighpass_, true, guard, 0);
    Configure(outputLowpass_, false, .5f * sampleRate_ - guard, 0);
    configuredLowerHz_ = configuredUpperHz_ =
        std::numeric_limits<float>::quiet_NaN();
    initialized_ = false;
    SetShiftHz(0.f);
    Reset();
  }

  void Reset() noexcept {
    for (auto *bank : {&sourceHighpass_, &sourceLowpass_,
                       &outputHighpass_, &outputLowpass_})
      for (auto &filter : *bank)
        filter.Reset();
  }

  void SetShiftHz(const float shiftHz) noexcept {
    const float safeShift = std::isfinite(shiftHz) ? shiftHz : 0.f;
    const float guard = GuardHz();
    const float transition = std::max(200.f, .01f * sampleRate_);
    const float nyquist = .5f * sampleRate_;
    const float quantum = std::max(8.f, .0005f * sampleRate_);
    const float requestedLower = std::max(guard, -safeShift + transition);
    const float requestedUpper = std::min(
        nyquist - guard, nyquist - safeShift - transition);
    // Quantize toward the rejected band. Slowly moving CV then redesigns only
    // when its safe boundary crosses a small frequency cell, never per sample.
    const float lower = quantum * std::ceil(requestedLower / quantum);
    const float upper = quantum * std::floor(requestedUpper / quantum);
    const std::size_t smoothing =
        initialized_ ? coefficientTransitionSamples_ : 0;
    if (lower != configuredLowerHz_) {
      Configure(sourceHighpass_, true, lower, smoothing);
      configuredLowerHz_ = lower;
    }
    if (upper != configuredUpperHz_) {
      Configure(sourceLowpass_, false, upper, smoothing);
      configuredUpperHz_ = upper;
    }
    initialized_ = true;
  }

  float ProcessSource(float input) noexcept {
    return ProcessBank(sourceLowpass_, ProcessBank(sourceHighpass_, input));
  }

  float ProcessOutput(float input) noexcept {
    return ProcessBank(outputLowpass_, ProcessBank(outputHighpass_, input));
  }

private:
  using Bank = std::array<Biquad, 2>;

  float GuardHz() const noexcept {
    return std::max(80.f, .003f * sampleRate_);
  }

  void Configure(Bank &bank, const bool highpass, const float frequency,
                 const std::size_t transitionSamples) noexcept {
    constexpr float Q1 = .541196100146197f;
    constexpr float Q2 = 1.306562964876377f;
    const auto design = [&](const float q) {
      return highpass ? biquad_design::Highpass(frequency, q, sampleRate_)
                      : biquad_design::Lowpass(frequency, q, sampleRate_);
    };
    bank[0].SetTargetCoefficients(design(Q1), transitionSamples);
    bank[1].SetTargetCoefficients(design(Q2), transitionSamples);
  }

  static float ProcessBank(Bank &bank, float input) noexcept {
    for (auto &filter : bank)
      input = filter.Process(input);
    return input;
  }

  Bank sourceHighpass_{};
  Bank sourceLowpass_{};
  Bank outputHighpass_{};
  Bank outputLowpass_{};
  float sampleRate_{48000.f};
  float configuredLowerHz_{};
  float configuredUpperHz_{};
  std::size_t coefficientTransitionSamples_{};
  bool initialized_{};
};

} // namespace tfdsp::percussion
