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
    secondaryExcitationProjection_.fill(1.f);
    SetStaticParameters(parameters);
  }

  // Static structural changes intentionally clear stored body energy.
  void SetStaticParameters(const Parameters &parameters) noexcept {
    constexpr float TwoPi = 6.28318530717958647692f;
    activeModeCount_ = 0;
    for (std::size_t source = 0; source < ModeCount; ++source) {
      const float inputGain = SafeInputGain(parameters[source].inputGain);
      const float outputGain = SafeOutputGain(parameters[source].outputGain);
      if (inputGain == 0.f || outputGain == 0.f)
        continue;
      const std::size_t mode = activeModeCount_++;
      sourceIndex_[mode] = source;
      const float frequency = std::clamp(
          tfdsp::FiniteNormalOrZero(parameters[source].frequencyHz), 1.f,
          .49f * sampleRate_);
      const float decay = std::clamp(
          tfdsp::FiniteNormalOrZero(parameters[source].decaySeconds), .001f,
          60.f);
      const float angle = TwoPi * frequency / sampleRate_;
      cosine_[mode] = std::cos(angle);
      sine_[mode] = std::sin(angle);
      radius_[mode] = std::exp(std::log(.001f) / (decay * sampleRate_));
      inputGain_[mode] = inputGain;
      outputGain_[mode] = outputGain;
      const float phase = std::clamp(
          tfdsp::FiniteNormalOrZero(parameters[source].inputPhaseRadians),
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
    return ProcessPreparedInputs<false>(
        input, unityProjection_, 0.f, unityProjection_,
        unityProjection_, damping);
  }

  void SetExcitationProjection(const Projection &projection) noexcept {
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode)
      excitationProjection_[mode] =
          SafeProjection(projection[sourceIndex_[mode]]);
  }

  void SetSecondaryExcitationProjection(
      const Projection &projection) noexcept {
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode)
      secondaryExcitationProjection_[mode] =
          SafeProjection(projection[sourceIndex_[mode]]);
  }

  std::size_t ActiveModeCount() const noexcept { return activeModeCount_; }

  float ProcessExcited(const float input,
                       const ModalDampingGains damping = {}) noexcept {
    return ProcessPreparedInputs<false>(
        input, excitationProjection_, 0.f, unityProjection_,
        unityProjection_, damping);
  }

  // Adds two independently projected forces to the same stored modal state.
  // This lets a local contact and a body-wide bloom coexist without a new
  // strike recolouring energy that is already circulating in the body.
  float ProcessExcitedPair(
      const float primaryInput, const float secondaryInput,
      const ModalDampingGains damping = {}) noexcept {
    return ProcessPreparedInputs<true>(
        primaryInput, excitationProjection_, secondaryInput,
        secondaryExcitationProjection_, unityProjection_, damping);
  }

  float ProcessSecondaryExcited(
      const float input, const ModalDampingGains damping = {}) noexcept {
    return ProcessPreparedInputs<false>(
        input, secondaryExcitationProjection_, 0.f, unityProjection_,
        unityProjection_, damping);
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
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode) {
      const std::size_t source = sourceIndex_[mode];
      const float priorReal = real_[mode];
      const float priorImaginary = imaginary_[mode];
      const float drive = inputGain_[mode] *
          SafeProjection(excitation[source]) * input;
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
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode)
      output += outputGain_[mode] *
          SafeProjection(observation[sourceIndex_[mode]]) *
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
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode)
      effectiveRadius_[mode] =
          radius_[mode] * safe.broadband * bandGain[band_[mode]];
    damping_ = safe;
  }

  template <bool HasSecondary>
  float ProcessPreparedInputs(
      float primaryInput, const Projection &primaryExcitation,
      float secondaryInput, const Projection &secondaryExcitation,
      const Projection &observation,
      const ModalDampingGains damping) noexcept {
    primaryInput = tfdsp::FiniteNormalOrZero(primaryInput);
    secondaryInput = tfdsp::FiniteNormalOrZero(secondaryInput);
    UpdateDamping(damping);
    if (activeModeCount_ == ModeCount)
      return ProcessPreparedFixed<ModeCount, HasSecondary>(
          primaryInput, primaryExcitation, secondaryInput,
          secondaryExcitation, observation);
    if constexpr (ModeCount == 2048) {
      switch (activeModeCount_) {
      case 64:
        return ProcessPreparedFixed<64, HasSecondary>(
            primaryInput, primaryExcitation, secondaryInput,
            secondaryExcitation, observation);
      case 128:
        return ProcessPreparedFixed<128, HasSecondary>(
            primaryInput, primaryExcitation, secondaryInput,
            secondaryExcitation, observation);
      case 256:
        return ProcessPreparedFixed<256, HasSecondary>(
            primaryInput, primaryExcitation, secondaryInput,
            secondaryExcitation, observation);
      case 512:
        return ProcessPreparedFixed<512, HasSecondary>(
            primaryInput, primaryExcitation, secondaryInput,
            secondaryExcitation, observation);
      case 768:
        return ProcessPreparedFixed<768, HasSecondary>(
            primaryInput, primaryExcitation, secondaryInput,
            secondaryExcitation, observation);
      case 1024:
        return ProcessPreparedFixed<1024, HasSecondary>(
            primaryInput, primaryExcitation, secondaryInput,
            secondaryExcitation, observation);
      case 1536:
        return ProcessPreparedFixed<1536, HasSecondary>(
            primaryInput, primaryExcitation, secondaryInput,
            secondaryExcitation, observation);
      case 2048:
        return ProcessPreparedFixed<2048, HasSecondary>(
            primaryInput, primaryExcitation, secondaryInput,
            secondaryExcitation, observation);
      default: break;
      }
    }
    return ProcessPreparedRuntime<HasSecondary>(
        primaryInput, primaryExcitation, secondaryInput,
        secondaryExcitation, observation, activeModeCount_);
  }

  template <std::size_t ActiveCount, bool HasSecondary>
  [[gnu::noinline]] float ProcessPreparedFixed(
      const float primaryInput, const Projection &primaryExcitation,
      const float secondaryInput, const Projection &secondaryExcitation,
      const Projection &observation) noexcept {
    for (std::size_t mode = 0; mode < ActiveCount; ++mode) {
      const float priorReal = real_[mode];
      const float priorImaginary = imaginary_[mode];
      float projectedInput = primaryExcitation[mode] * primaryInput;
      if constexpr (HasSecondary)
        projectedInput += secondaryExcitation[mode] * secondaryInput;
      const float drive = inputGain_[mode] * projectedInput;
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
    for (std::size_t mode = 0; mode < ActiveCount; ++mode)
      output += outputGain_[mode] * observation[mode] * real_[mode];
    return tfdsp::FiniteNormalOrZero(output);
  }

  template <bool HasSecondary>
  float ProcessPreparedRuntime(
      const float primaryInput, const Projection &primaryExcitation,
      const float secondaryInput, const Projection &secondaryExcitation,
      const Projection &observation,
      const std::size_t activeCount) noexcept {
    for (std::size_t mode = 0; mode < activeCount; ++mode) {
      const float priorReal = real_[mode];
      const float priorImaginary = imaginary_[mode];
      float projectedInput = primaryExcitation[mode] * primaryInput;
      if constexpr (HasSecondary)
        projectedInput += secondaryExcitation[mode] * secondaryInput;
      const float drive = inputGain_[mode] * projectedInput;
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
    for (std::size_t mode = 0; mode < activeCount; ++mode)
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
  std::array<std::size_t, ModeCount> sourceIndex_{};
  Projection unityProjection_{};
  Projection excitationProjection_{};
  Projection secondaryExcitationProjection_{};
  ModalDampingGains damping_{};
  float sampleRate_{48000.f};
  float lowCrossoverHz_{700.f};
  float highCrossoverHz_{6500.f};
  std::size_t activeModeCount_{ModeCount};
};

} // namespace tfdsp::percussion
