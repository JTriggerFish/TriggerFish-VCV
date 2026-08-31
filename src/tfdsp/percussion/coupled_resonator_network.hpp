#pragma once

#include "orthogonal_mixer.hpp"
#include "static_fractional_delay.hpp"
#include "three_band_decay_filter.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace tfdsp::percussion {

struct ResonatorLineParameters {
  float delaySamples{32.f};
  ThreeBandDecayTimes decay{};
  float inputGain{1.f};
  float outputGain{1.f};
};

// Wet-only fractional feedback-comb bank with optional energy-preserving
// coupling. At zero coupling it reduces exactly to independent parallel combs.
template <std::size_t LineCount> class CoupledResonatorNetwork {
public:
  using Parameters = std::array<ResonatorLineParameters, LineCount>;
  using Frame = std::array<float, LineCount>;

  void Prepare(const float sampleRate, const float maximumDelaySamples,
               const Parameters &parameters, const float lowCrossoverHz,
               const float highCrossoverHz) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("resonator sample rate must be positive");
    sampleRate_ = sampleRate;
    for (std::size_t line = 0; line < LineCount; ++line) {
      delays_[line].Prepare(maximumDelaySamples,
                            parameters[line].delaySamples);
      losses_[line].Prepare(sampleRate, lowCrossoverHz, highCrossoverHz);
    }
    SetStaticParameters(parameters);
    Reset();
  }

  void Reset() noexcept {
    for (std::size_t line = 0; line < LineCount; ++line) {
      delays_[line].Reset();
      losses_[line].Reset();
    }
  }

  void SetCoupling(const float coupling) noexcept {
    mixer_.SetCoupling(coupling);
  }

  // Static retuning clears delay state. Live pitch motion requires a separate
  // state-preserving transition and must not call this configuration method.
  void SetStaticParameters(const Parameters &parameters) {
    for (std::size_t line = 0; line < LineCount; ++line) {
      delays_[line].SetDelaySamples(parameters[line].delaySamples);
      losses_[line].Reset();
      parameters_[line] = parameters[line];
      parameters_[line].delaySamples = delays_[line].DelaySamples();
      if (!std::isfinite(parameters_[line].inputGain))
        parameters_[line].inputGain = 0.f;
      if (!std::isfinite(parameters_[line].outputGain))
        parameters_[line].outputGain = 0.f;
      const float pathSeconds = parameters_[line].delaySamples / sampleRate_;
      losses_[line].SetDecayTimes(pathSeconds, parameters_[line].decay);
    }
  }

  float Process(float input) noexcept {
    if (!std::isfinite(input))
      input = 0.f;
    Frame wet{};
    float output = 0.f;
    for (std::size_t line = 0; line < LineCount; ++line) {
      wet[line] = delays_[line].Read();
      output += parameters_[line].outputGain * wet[line];
    }
    const Frame feedback = mixer_.Process(wet);
    for (std::size_t line = 0; line < LineCount; ++line) {
      const float drive = parameters_[line].inputGain * input;
      delays_[line].Push(drive + losses_[line].Process(feedback[line]));
    }
    return std::isfinite(output) ? output : 0.f;
  }

private:
  std::array<StaticFractionalDelay, LineCount> delays_{};
  std::array<ThreeBandDecayFilter, LineCount> losses_{};
  OrthogonalMixer<LineCount> mixer_{};
  Parameters parameters_{};
  float sampleRate_{48000.f};
};

} // namespace tfdsp::percussion
