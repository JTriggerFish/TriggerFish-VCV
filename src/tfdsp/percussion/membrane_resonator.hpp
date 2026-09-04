#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace tfdsp::percussion {

struct MembraneModeParameters {
  float frequencyHz{100.f};
  float decaySeconds{1.f};
  float inputGain{1.f};
  float outputGain{1.f};
  float centerProjection{1.f};
  float edgeProjection{1.f};
};

// A shared, dynamically tensioned modal membrane. Coefficients are refreshed
// at a short control interval while quadrature state remains continuous.
template <std::size_t ModeCount> class MembraneResonator {
public:
  using Parameters = std::array<MembraneModeParameters, ModeCount>;
  using Drive = std::array<float, ModeCount>;
  struct PreparedParameters {
    std::array<float, ModeCount> frequencies{};
    std::array<float, ModeCount> radii{};
    std::array<float, ModeCount> cosines{};
    std::array<float, ModeCount> sines{};
    std::array<float, ModeCount> inputGains{};
    std::array<float, ModeCount> outputGains{};
    std::array<float, ModeCount> center{};
    std::array<float, ModeCount> edge{};
    float sampleRate{48000.f};
    float maximumStoredEnergy{1.f};
  };

  void Prepare(const float sampleRate, const Parameters &parameters,
               const float maximumStoredEnergy = 1.f) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("membrane sample rate must be positive");
    LoadPrepared(PrepareParameters(sampleRate, parameters,
                                   maximumStoredEnergy));
  }

  static PreparedParameters PrepareParameters(
      const float sampleRate, const Parameters &parameters,
      const float maximumStoredEnergy = 1.f) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("membrane sample rate must be positive");
    PreparedParameters result;
    result.sampleRate = sampleRate;
    float driveNormSquared = 0.f;
    float outputNormSquared = 0.f;
    for (std::size_t mode = 0; mode < ModeCount; ++mode) {
      const auto &source = parameters[mode];
      result.frequencies[mode] = std::clamp(
          tfdsp::FiniteNormalOrZero(source.frequencyHz), 1.f,
          .45f * sampleRate);
      const float decay = std::clamp(
          tfdsp::FiniteNormalOrZero(source.decaySeconds), .002f, 30.f);
      result.radii[mode] = std::exp(std::log(.001f) / (decay * sampleRate));
      result.inputGains[mode] = SafeGain(source.inputGain, 64.f);
      result.outputGains[mode] = SafeGain(source.outputGain, 16.f);
      result.center[mode] = SafeGain(source.centerProjection, 4.f);
      result.edge[mode] = SafeGain(source.edgeProjection, 4.f);
      const float projection = std::max(std::abs(result.center[mode]),
                                        std::abs(result.edge[mode]));
      driveNormSquared += result.inputGains[mode] * result.inputGains[mode] *
                          projection * projection;
      outputNormSquared += result.outputGains[mode] * result.outputGains[mode];
    }
    // Keep loudness and stored energy independent of mode count. The two
    // constants define the bank's calibrated force and observation units.
    const float driveScale = .25f / std::sqrt(std::max(1.e-12f,
                                                       driveNormSquared));
    const float outputScale = 2.f / std::sqrt(std::max(1.e-12f,
                                                        outputNormSquared));
    constexpr float TwoPi = 6.28318530717958647692f;
    for (std::size_t mode = 0; mode < ModeCount; ++mode) {
      result.inputGains[mode] *= driveScale;
      result.outputGains[mode] *= outputScale;
      const float angle = TwoPi * result.frequencies[mode] / sampleRate;
      result.cosines[mode] = std::cos(angle);
      result.sines[mode] = std::sin(angle);
    }
    result.maximumStoredEnergy = std::clamp(
        tfdsp::FiniteNormalOrZero(maximumStoredEnergy), .001f, 64.f);
    return result;
  }

  void LoadPrepared(const PreparedParameters &prepared) noexcept {
    sampleRate_ = prepared.sampleRate;
    frequencies_ = prepared.frequencies;
    radii_ = prepared.radii;
    cosines_ = prepared.cosines;
    sines_ = prepared.sines;
    inputGains_ = prepared.inputGains;
    outputGains_ = prepared.outputGains;
    center_ = prepared.center;
    edge_ = prepared.edge;
    maximumStoredEnergy_ = prepared.maximumStoredEnergy;
    tensionScale_ = 1.f;
    Reset();
  }

  void Reset() noexcept {
    real_.fill(0.f);
    imaginary_.fill(0.f);
    updateCountdown_ = 0;
  }

  Drive Project(const float force, const float location) const noexcept {
    Drive result{};
    const float safeForce = tfdsp::FiniteNormalOrZero(force);
    const float position = std::clamp(
        tfdsp::FiniteNormalOrZero(location), 0.f, 1.f);
    for (std::size_t mode = 0; mode < ModeCount; ++mode)
      result[mode] = safeForce *
          (center_[mode] + position * (edge_[mode] - center_[mode]));
    return result;
  }

  float Process(const Drive &drive, const float tensionScale) noexcept {
    UpdateCoefficients(tensionScale);
    std::array<float, ModeCount> forces{};
    float baseEnergy = 0.f;
    float driveEnergy = 0.f;
    float crossEnergy = 0.f;
    for (std::size_t mode = 0; mode < ModeCount; ++mode) {
      const float priorReal = real_[mode];
      const float priorImaginary = imaginary_[mode];
      forces[mode] = inputGains_[mode] *
          tfdsp::FiniteNormalOrZero(drive[mode]);
      real_[mode] = tfdsp::FiniteNormalOrZero(
          radii_[mode] * (cosines_[mode] * priorReal -
                          sines_[mode] * priorImaginary));
      imaginary_[mode] = tfdsp::FiniteNormalOrZero(
          radii_[mode] * (sines_[mode] * priorReal +
                          cosines_[mode] * priorImaginary));
      baseEnergy += real_[mode] * real_[mode] +
                    imaginary_[mode] * imaginary_[mode];
      crossEnergy += real_[mode] * forces[mode];
      driveEnergy += forces[mode] * forces[mode];
    }
    const float driveScale = AvailableDriveScale(
        baseEnergy, crossEnergy, driveEnergy);
    float output = 0.f;
    for (std::size_t mode = 0; mode < ModeCount; ++mode) {
      real_[mode] = tfdsp::FiniteNormalOrZero(
          real_[mode] + driveScale * forces[mode]);
      output += outputGains_[mode] * real_[mode];
    }
    return tfdsp::FiniteNormalOrZero(output);
  }

  float StoredEnergy() const noexcept {
    float result = 0.f;
    for (std::size_t mode = 0; mode < ModeCount; ++mode)
      result += real_[mode] * real_[mode] +
                imaginary_[mode] * imaginary_[mode];
    return tfdsp::FiniteNormalOrZero(result);
  }

private:
  static float SafeGain(const float value, const float maximum) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(value), -maximum, maximum);
  }

  float AvailableDriveScale(const float baseEnergy, const float crossEnergy,
                            const float driveEnergy) const noexcept {
    const float proposed = baseEnergy + 2.f * crossEnergy + driveEnergy;
    if (proposed <= maximumStoredEnergy_) return 1.f;
    if (driveEnergy <= 1.e-20f) return 0.f;
    const float discriminant = crossEnergy * crossEnergy + driveEnergy *
        std::max(0.f, maximumStoredEnergy_ - baseEnergy);
    return std::clamp((-crossEnergy + std::sqrt(discriminant)) /
                          driveEnergy,
                      0.f, 1.f);
  }

  void UpdateCoefficients(float tensionScale) noexcept {
    tensionScale = std::clamp(
        std::isfinite(tensionScale) ? tensionScale : 1.f, .5f, 2.f);
    if (updateCountdown_ != 0) {
      --updateCountdown_;
      return;
    }
    updateCountdown_ = CoefficientInterval - 1;
    if (std::abs(tensionScale - tensionScale_) < 1.e-6f) return;
    tensionScale_ = tensionScale;
    constexpr float TwoPi = 6.28318530717958647692f;
    for (std::size_t mode = 0; mode < ModeCount; ++mode) {
      const float frequency = std::min(
          .45f * sampleRate_, frequencies_[mode] * tensionScale_);
      const float angle = TwoPi * frequency / sampleRate_;
      cosines_[mode] = std::cos(angle);
      sines_[mode] = std::sin(angle);
    }
  }

  static constexpr std::size_t CoefficientInterval = 16;
  std::array<float, ModeCount> real_{};
  std::array<float, ModeCount> imaginary_{};
  std::array<float, ModeCount> frequencies_{};
  std::array<float, ModeCount> radii_{};
  std::array<float, ModeCount> cosines_{};
  std::array<float, ModeCount> sines_{};
  std::array<float, ModeCount> inputGains_{};
  std::array<float, ModeCount> outputGains_{};
  std::array<float, ModeCount> center_{};
  std::array<float, ModeCount> edge_{};
  float sampleRate_{48000.f};
  float tensionScale_{-1.f};
  float maximumStoredEnergy_{1.f};
  std::size_t updateCountdown_{};
};

} // namespace tfdsp::percussion
