#pragma once

#include "deterministic_random.hpp"
#include "modal_constraint.hpp"
#include "three_band_decay_filter.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace tfdsp::percussion {

struct TurbulentResidualParameters {
  std::array<float, 3> gain{};
  ThreeBandDecayTimes decay{.8f, .55f, .3f};
  float injectionSeconds{.02f};
  float lowCrossoverHz{700.f};
  float highCrossoverHz{6500.f};
  std::uint32_t seed{0x54555242u};
};

// Continuous stochastic presentation of unresolved nonlinear mode coupling.
// Squared propagation energy fills three passive reservoirs; one correlated,
// complementary noise field reads them without exposing grains or transients.
class TurbulentResidual {
public:
  void Prepare(const float sampleRate,
               const TurbulentResidualParameters &parameters) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("turbulent-residual sample rate must be positive");
    sampleRate_ = sampleRate;
    Configure(parameters);
    Reset();
  }

  void Reset() noexcept {
    driveLow_ = driveBelowHigh_ = 0.f;
    noiseLow_ = noiseBelowHigh_ = 0.f;
    energy_ = {};
    random_.Seed(seed_);
  }

  void Seed(const std::uint32_t seed) noexcept {
    random_.Seed(seed == 0 ? seed_ : seed);
  }

  float Process(float drive,
                const ModalDampingGains constraint = {}) noexcept {
    const auto driveBands = Split(tfdsp::FiniteNormalOrZero(drive),
                                  driveLow_, driveBelowHigh_);
    const auto noiseBands = Split(1.7320508075688772f * random_.Bipolar(),
                                  noiseLow_, noiseBelowHigh_);
    const std::array<float, 3> loss{
        Unit(constraint.broadband) * Unit(constraint.low),
        Unit(constraint.broadband) * Unit(constraint.middle),
        Unit(constraint.broadband) * Unit(constraint.high)};
    float output = 0.f;
    for (std::size_t band = 0; band < energy_.size(); ++band) {
      energy_[band] = decayMultiplier_[band] * loss[band] * loss[band] *
              energy_[band] +
          injectionGain_ * driveBands[band] * driveBands[band];
      energy_[band] = std::clamp(
          tfdsp::FiniteNormalOrZero(energy_[band]), 0.f, 16.f);
      output += gain_[band] * std::sqrt(energy_[band]) * noiseBands[band];
    }
    return tfdsp::FiniteNormalOrZero(output);
  }

  float StoredEnergy() const noexcept {
    return energy_[0] + energy_[1] + energy_[2];
  }

private:
  static float Unit(const float value) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(value), 0.f, 1.f);
  }

  std::array<float, 3> Split(const float input, float &low,
                             float &belowHigh) const noexcept {
    low += lowCoefficient_ * (input - low);
    belowHigh += highCoefficient_ * (input - belowHigh);
    low = tfdsp::FiniteNormalOrZero(low);
    belowHigh = tfdsp::FiniteNormalOrZero(belowHigh);
    return {low, belowHigh - low, input - belowHigh};
  }

  void Configure(const TurbulentResidualParameters &parameters) noexcept {
    for (std::size_t band = 0; band < gain_.size(); ++band)
      gain_[band] = std::clamp(
          tfdsp::FiniteNormalOrZero(parameters.gain[band]), 0.f, 4.f);
    const float injectionSeconds = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters.injectionSeconds), .0001f, 1.f);
    injectionGain_ = 1.f / (sampleRate_ * injectionSeconds);
    const std::array<float, 3> decay{
        parameters.decay.lowSeconds, parameters.decay.middleSeconds,
        parameters.decay.highSeconds};
    for (std::size_t band = 0; band < decay.size(); ++band) {
      const float seconds = std::clamp(
          tfdsp::FiniteNormalOrZero(decay[band]), .005f, 30.f);
      decayMultiplier_[band] = std::pow(10.f, -6.f / (sampleRate_ * seconds));
    }
    lowCoefficient_ = Pole(parameters.lowCrossoverHz);
    highCoefficient_ = Pole(
        std::max(parameters.lowCrossoverHz, parameters.highCrossoverHz));
    seed_ = parameters.seed;
  }

  float Pole(float frequencyHz) const noexcept {
    frequencyHz = std::clamp(tfdsp::FiniteNormalOrZero(frequencyHz),
                             0.f, .49f * sampleRate_);
    return 1.f - std::exp(-6.283185307179586f * frequencyHz / sampleRate_);
  }

  DeterministicRandom random_{};
  std::array<float, 3> gain_{};
  std::array<float, 3> decayMultiplier_{};
  std::array<float, 3> energy_{};
  float sampleRate_{48000.f};
  float injectionGain_{};
  float lowCoefficient_{};
  float highCoefficient_{};
  float driveLow_{};
  float driveBelowHigh_{};
  float noiseLow_{};
  float noiseBelowHigh_{};
  std::uint32_t seed_{0x54555242u};
};

} // namespace tfdsp::percussion
