#pragma once

#include "modal_constraint.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace tfdsp::percussion {

struct ModalModeParameters {
  float frequencyHz{1000.f};
  float decaySeconds{1.f};
  float inputGain{1.f};
  float outputGain{1.f};
  float inputPhaseRadians{};
};

// Fixed-capacity bank of independently placed damped complex modes. The bank
// owns one mono body state; excitation and observation projections do not
// duplicate that state.
template <std::size_t ModeCount> class ModalBank {
public:
  static_assert(ModeCount > 0, "a modal bank needs at least one mode");
  using Parameters = std::array<ModalModeParameters, ModeCount>;
  using Projection = std::array<float, ModeCount>;

  void Prepare(const float sampleRate, const Parameters &parameters,
               const float lowCrossoverHz, const float highCrossoverHz) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("modal bank sample rate must be positive");
    if (!(lowCrossoverHz > 0.f && highCrossoverHz > lowCrossoverHz))
      throw std::invalid_argument("modal bank crossovers must be ordered");
    sampleRate_ = sampleRate;
    lowCrossoverHz_ = lowCrossoverHz;
    highCrossoverHz_ = highCrossoverHz;
    unityProjection_.fill(1.f);
    excitationProjection_.fill(1.f);
    SetStaticParameters(parameters);
  }

  // Static structural changes intentionally clear stored body energy.
  void SetStaticParameters(const Parameters &parameters) noexcept {
    constexpr float TwoPi = 6.28318530717958647692f;
    for (std::size_t mode = 0; mode < ModeCount; ++mode) {
      const float frequency = std::clamp(
          tfdsp::FiniteNormalOrZero(parameters[mode].frequencyHz), 1.f,
          .49f * sampleRate_);
      const float decay = std::clamp(
          tfdsp::FiniteNormalOrZero(parameters[mode].decaySeconds), .001f, 60.f);
      const float angle = TwoPi * frequency / sampleRate_;
      cosine_[mode] = std::cos(angle);
      sine_[mode] = std::sin(angle);
      radius_[mode] = std::exp(std::log(.001f) / (decay * sampleRate_));
      inputGain_[mode] = SafeInputGain(parameters[mode].inputGain);
      outputGain_[mode] = SafeOutputGain(parameters[mode].outputGain);
      const float phase = std::clamp(
          tfdsp::FiniteNormalOrZero(parameters[mode].inputPhaseRadians),
          -3.14159265358979323846f, 3.14159265358979323846f);
      inputPhaseCosine_[mode] = std::cos(phase);
      inputPhaseSine_[mode] = std::sin(phase);
      band_[mode] = frequency < lowCrossoverHz_
          ? 0
          : static_cast<std::uint8_t>(frequency < highCrossoverHz_ ? 1 : 2);
      effectiveRadius_[mode] = radius_[mode];
    }
    damping_ = {};
    Reset();
  }

  void Reset() noexcept {
    real_.fill(0.f);
    imaginary_.fill(0.f);
  }

  float Process(const float input,
                const ModalDampingGains damping = {}) noexcept {
    return ProcessPrepared(input, unityProjection_, unityProjection_, damping);
  }

  void SetExcitationProjection(const Projection &projection) noexcept {
    for (std::size_t mode = 0; mode < ModeCount; ++mode)
      excitationProjection_[mode] = SafeProjection(projection[mode]);
  }

  float ProcessExcited(const float input,
                       const ModalDampingGains damping = {}) noexcept {
    return ProcessPrepared(
        input, excitationProjection_, unityProjection_, damping);
  }

  float ProcessExcited(const float input, const Projection &excitation,
                       const ModalDampingGains damping = {}) noexcept {
    return ProcessProjected(input, excitation, unityProjection_, damping);
  }

  float ProcessProjected(float input, const Projection &excitation,
                         const Projection &observation,
                         const ModalDampingGains damping = {}) noexcept {
    input = tfdsp::FiniteNormalOrZero(input);
    UpdateDamping(damping);
    for (std::size_t mode = 0; mode < ModeCount; ++mode) {
      const float priorReal = real_[mode];
      const float priorImaginary = imaginary_[mode];
      const float drive = inputGain_[mode] *
          SafeProjection(excitation[mode]) * input;
      real_[mode] = tfdsp::FiniteNormalOrZero(
          effectiveRadius_[mode] * (cosine_[mode] * priorReal -
                                    sine_[mode] * priorImaginary) +
          inputPhaseCosine_[mode] * drive);
      imaginary_[mode] = tfdsp::FiniteNormalOrZero(
          effectiveRadius_[mode] * (sine_[mode] * priorReal +
                                    cosine_[mode] * priorImaginary) +
          inputPhaseSine_[mode] * drive);
    }
    float output = 0.f;
    for (std::size_t mode = 0; mode < ModeCount; ++mode)
      output += outputGain_[mode] * SafeProjection(observation[mode]) *
          real_[mode];
    return tfdsp::FiniteNormalOrZero(output);
  }

private:
  static float SafeInputGain(const float gain) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(gain), -256.f, 256.f);
  }

  static float SafeOutputGain(const float gain) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(gain), -16.f, 16.f);
  }

  static float SafeProjection(const float projection) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(projection), -4.f, 4.f);
  }

  static float UnitGain(const float gain) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(gain), 0.f, 1.f);
  }

  void UpdateDamping(const ModalDampingGains damping) noexcept {
    const ModalDampingGains safe{
        UnitGain(damping.broadband), UnitGain(damping.low),
        UnitGain(damping.middle), UnitGain(damping.high)};
    if (safe.broadband == damping_.broadband && safe.low == damping_.low &&
        safe.middle == damping_.middle && safe.high == damping_.high)
      return;
    const std::array<float, 3> bandGain{
        safe.low, safe.middle, safe.high};
    for (std::size_t mode = 0; mode < ModeCount; ++mode)
      effectiveRadius_[mode] =
          radius_[mode] * safe.broadband * bandGain[band_[mode]];
    damping_ = safe;
  }

  float ProcessPrepared(float input, const Projection &excitation,
                        const Projection &observation,
                        const ModalDampingGains damping) noexcept {
    input = tfdsp::FiniteNormalOrZero(input);
    UpdateDamping(damping);
    for (std::size_t mode = 0; mode < ModeCount; ++mode) {
      const float priorReal = real_[mode];
      const float priorImaginary = imaginary_[mode];
      const float drive = inputGain_[mode] * excitation[mode] * input;
      real_[mode] = tfdsp::FiniteNormalOrZero(
          effectiveRadius_[mode] * (cosine_[mode] * priorReal -
                                    sine_[mode] * priorImaginary) +
          inputPhaseCosine_[mode] * drive);
      imaginary_[mode] = tfdsp::FiniteNormalOrZero(
          effectiveRadius_[mode] * (sine_[mode] * priorReal +
                                    cosine_[mode] * priorImaginary) +
          inputPhaseSine_[mode] * drive);
    }
    float output = 0.f;
    for (std::size_t mode = 0; mode < ModeCount; ++mode)
      output += outputGain_[mode] * observation[mode] * real_[mode];
    return tfdsp::FiniteNormalOrZero(output);
  }

  std::array<float, ModeCount> real_{};
  std::array<float, ModeCount> imaginary_{};
  std::array<float, ModeCount> cosine_{};
  std::array<float, ModeCount> sine_{};
  std::array<float, ModeCount> radius_{};
  std::array<float, ModeCount> effectiveRadius_{};
  std::array<float, ModeCount> inputGain_{};
  std::array<float, ModeCount> outputGain_{};
  std::array<float, ModeCount> inputPhaseCosine_{};
  std::array<float, ModeCount> inputPhaseSine_{};
  std::array<std::uint8_t, ModeCount> band_{};
  Projection unityProjection_{};
  Projection excitationProjection_{};
  ModalDampingGains damping_{};
  float sampleRate_{48000.f};
  float lowCrossoverHz_{700.f};
  float highCrossoverHz_{6500.f};
};

} // namespace tfdsp::percussion
