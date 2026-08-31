#pragma once

#include "fir_hilbert_transform.hpp"
#include "quadrature_oscillator.hpp"
#include "translation_band_limiter.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {

// Single-sideband frequency shifter. Positive and negative offsets share the
// same phase-continuous oscillator; zero offset is the exactly delayed input.
class FrequencyShifter {
public:
  static constexpr std::size_t LatencySamples =
      FirHilbertTransform::LatencySamples;

  void Prepare(const float sampleRate) {
    oscillator_.Prepare(sampleRate);
    limiter_.Prepare(sampleRate);
    limiterShiftHz_ = 0.f;
    mixStep_ = 1.f / std::max(1.f, .002f * sampleRate);
    Reset();
  }

  void Reset() noexcept {
    hilbert_.Reset();
    oscillator_.Reset();
    limiter_.Reset();
    limiterMix_ = 0.f;
  }

  void SetShiftHz(const float shiftHz) noexcept {
    const float safeShift = std::isfinite(shiftHz) ? shiftHz : 0.f;
    oscillator_.SetFrequencyHz(safeShift);
    // Hosts commonly resend unchanged parameters for every sample. Avoid eight
    // redundant trigonometric biquad designs in that normal real-time path.
    if (safeShift != limiterShiftHz_) {
      limiter_.SetShiftHz(safeShift);
      limiterShiftHz_ = safeShift;
    }
    limiterTarget_ = std::clamp((std::abs(safeShift) - 80.f) / 160.f, 0.f, 1.f);
  }

  float Process(float input) noexcept {
    if (!std::isfinite(input))
      input = 0.f;
    if (limiterMix_ < limiterTarget_)
      limiterMix_ = std::min(limiterTarget_, limiterMix_ + mixStep_);
    else
      limiterMix_ = std::max(limiterTarget_, limiterMix_ - mixStep_);
    const float limitedInput = limiter_.ProcessSource(input);
    const auto analytic = hilbert_.Process(
        input + limiterMix_ * (limitedInput - input));
    const auto rotation = oscillator_.Process();
    const float shifted = analytic.real * rotation.cosine -
                          analytic.quadrature * rotation.sine;
    const float limited = limiter_.ProcessOutput(shifted);
    const float output = shifted + limiterMix_ * (limited - shifted);
    return std::isfinite(output) ? output : 0.f;
  }

private:
  FirHilbertTransform hilbert_{};
  QuadratureOscillator oscillator_{};
  TranslationBandLimiter limiter_{};
  float limiterMix_{};
  float limiterTarget_{};
  float limiterShiftHz_{};
  float mixStep_{1.f};
};

} // namespace tfdsp::percussion
