#pragma once

#include "observation_delay.hpp"
#include "radiation_filter.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace tfdsp::percussion {

struct ObservationPathParameters {
  float gain{1.f};
  float delaySeconds{};
  bool invertPolarity{};
  bool radiationEnabled{true};
  RadiationFilterParameters radiation{};
};

// One fitted microphone/listening view over explicit source-state taps. It
// changes observation only; it never feeds energy back into an instrument.
template <std::size_t SourceCount> class ObservationModel {
public:
  using SourceFrame = std::array<float, SourceCount>;
  using Parameters = std::array<ObservationPathParameters, SourceCount>;

  void Prepare(const float sampleRate, const float maximumDelaySeconds,
               const Parameters &parameters) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("observation sample rate must be positive");
    sampleRate_ = sampleRate;
    maximumDelaySeconds_ = std::clamp(
        std::isfinite(maximumDelaySeconds) ? maximumDelaySeconds : 0.f,
        1.f / sampleRate_, 10.f);
    const float capacity = maximumDelaySeconds_ * sampleRate_;
    for (std::size_t source = 0; source < SourceCount; ++source) {
      delays_[source].Prepare(capacity, 0.f);
      filters_[source].Prepare(sampleRate, parameters[source].radiation);
    }
    SetStaticParameters(parameters);
    Reset();
  }

  void Reset() noexcept {
    for (std::size_t source = 0; source < SourceCount; ++source) {
      delays_[source].Reset();
      filters_[source].Reset();
    }
  }

  void SetStaticParameters(const Parameters &parameters) {
    for (std::size_t source = 0; source < SourceCount; ++source) {
      parameters_[source] = parameters[source];
      const float magnitude = std::clamp(
          tfdsp::FiniteNormalOrZero(parameters[source].gain), 0.f, 16.f);
      gains_[source] = parameters[source].invertPolarity ? -magnitude : magnitude;
      const float seconds = std::clamp(
          std::isfinite(parameters[source].delaySeconds)
              ? parameters[source].delaySeconds : 0.f,
          0.f, maximumDelaySeconds_);
      delays_[source].SetStaticDelaySamples(seconds * sampleRate_);
      filters_[source].SetStaticParameters(parameters[source].radiation);
    }
  }

  float Process(const SourceFrame &sources) noexcept {
    float output = 0.f;
    for (std::size_t source = 0; source < SourceCount; ++source) {
      float path = delays_[source].Process(sources[source]);
      if (parameters_[source].radiationEnabled)
        path = filters_[source].Process(path);
      output += gains_[source] * path;
    }
    return tfdsp::FiniteNormalOrZero(output);
  }

private:
  std::array<ObservationDelay, SourceCount> delays_{};
  std::array<RadiationFilter, SourceCount> filters_{};
  Parameters parameters_{};
  std::array<float, SourceCount> gains_{};
  float sampleRate_{48000.f};
  float maximumDelaySeconds_{.1f};
};

} // namespace tfdsp::percussion
