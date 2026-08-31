#pragma once

#include "fractional_delay_line.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {

// High-quality moving delay. The caller owns the modulation law; this class
// bounds the tap causally and performs a continuous 12-tap sinc interpolation.
class ModulatedFractionalDelay {
public:
  static constexpr float MinimumDelaySamples =
      FractionalDelayLine::MinimumSincDelaySamples;

  void Prepare(const float maximumDelaySamples) {
    line_.Prepare(maximumDelaySamples);
  }

  void Reset() noexcept { line_.Reset(); }

  float Process(const float input, float delaySamples) noexcept {
    if (!std::isfinite(delaySamples))
      delaySamples = MinimumDelaySamples;
    const float output = line_.ReadSinc(delaySamples);
    line_.Push(input);
    return tfdsp::FiniteNormalOrZero(output);
  }

  float MaximumDelaySamples() const noexcept {
    return line_.MaximumDelaySamples();
  }

private:
  FractionalDelayLine line_{};
};

} // namespace tfdsp::percussion
