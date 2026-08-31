#pragma once

#include "biquad.hpp"
#include "biquad_design.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tfdsp::percussion {

struct RadiationFilterParameters {
  float lowCutHz{40.f};
  float lowCutQ{.707f};
  float colourFrequencyHz{5000.f};
  float colourGainDb{0.f};
  float colourQ{.8f};
  float highCutHz{18000.f};
  float highCutQ{.707f};
  float outputGain{1.f};
};

// Static object-radiation EQ: high-pass, fitted colour peak, then low-pass.
// Live location morphing is performed by a later state-preserving wrapper.
class RadiationFilter {
public:
  void Prepare(const float sampleRate,
               const RadiationFilterParameters &parameters) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("radiation sample rate must be positive");
    sampleRate_ = sampleRate;
    SetStaticParameters(parameters);
    Reset();
  }

  void Reset() noexcept {
    highpass_.Reset();
    colour_.Reset();
    lowpass_.Reset();
  }

  void SetStaticParameters(const RadiationFilterParameters &parameters) noexcept {
    highpass_.SetCoefficients(biquad_design::Highpass(
        parameters.lowCutHz, parameters.lowCutQ, sampleRate_));
    colour_.SetCoefficients(biquad_design::Peaking(
        parameters.colourFrequencyHz, parameters.colourQ,
        parameters.colourGainDb, sampleRate_));
    lowpass_.SetCoefficients(biquad_design::Lowpass(
        parameters.highCutHz, parameters.highCutQ, sampleRate_));
    outputGain_ = std::clamp(
        std::isfinite(parameters.outputGain) ? parameters.outputGain : 0.f,
        0.f, 16.f);
  }

  float Process(const float input) noexcept {
    const float output = outputGain_ *
        lowpass_.Process(colour_.Process(highpass_.Process(input)));
    if (!std::isfinite(output)) {
      Reset();
      return 0.f;
    }
    return output;
  }

private:
  Biquad highpass_{};
  Biquad colour_{};
  Biquad lowpass_{};
  float sampleRate_{48000.f};
  float outputGain_{1.f};
};

} // namespace tfdsp::percussion
