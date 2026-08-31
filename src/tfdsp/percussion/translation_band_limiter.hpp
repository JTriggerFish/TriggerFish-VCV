#pragma once

#include "biquad.hpp"
#include "biquad_design.hpp"

#include <array>
#include <cmath>
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
    constexpr float Q1 = .541196100146197f;
    constexpr float Q2 = 1.306562964876377f;
    const float guard = std::max(80.f, .003f * sampleRate_);
    const float transition = std::max(200.f, .01f * sampleRate_);
    const float nyquist = .5f * sampleRate_;
    const float lower = std::max(guard, -shiftHz + transition);
    const float upper = std::min(nyquist - guard,
                                 nyquist - shiftHz - transition);
    Configure(sourceHighpass_, true, lower, Q1, Q2);
    Configure(sourceLowpass_, false, upper, Q1, Q2);
    Configure(outputHighpass_, true, guard, Q1, Q2);
    Configure(outputLowpass_, false, nyquist - guard, Q1, Q2);
  }

  float ProcessSource(float input) noexcept {
    return ProcessBank(sourceLowpass_, ProcessBank(sourceHighpass_, input));
  }

  float ProcessOutput(float input) noexcept {
    return ProcessBank(outputLowpass_, ProcessBank(outputHighpass_, input));
  }

private:
  using Bank = std::array<Biquad, 2>;

  void Configure(Bank &bank, const bool highpass, const float frequency,
                 const float q1, const float q2) noexcept {
    const auto design = [&](const float q) {
      return highpass ? biquad_design::Highpass(frequency, q, sampleRate_)
                      : biquad_design::Lowpass(frequency, q, sampleRate_);
    };
    bank[0].SetCoefficients(design(q1));
    bank[1].SetCoefficients(design(q2));
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
};

} // namespace tfdsp::percussion
