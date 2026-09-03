#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace tfdsp::percussion {

struct StochasticModalModeParameters {
  float frequencyHz{1000.f};
  float decaySeconds{1.f};
  float inputGain{1.f};
  float outputGain{1.f};
  float inputPhaseRadians{};
  float phaseBandwidthHz{};
  std::uint16_t packet{};
  float exchangeAmount{1.f};
};

struct StochasticModalFieldControls {
  float exchangeAngleRadians{};
  std::uint32_t seed{0x4649454cu};
};

template <std::size_t ModeCount> struct PreparedStochasticModalField {
  std::array<float, ModeCount> cosineCentre{};
  std::array<float, ModeCount> cosineSpread{};
  std::array<float, ModeCount> sineCentre{};
  std::array<float, ModeCount> sineSpread{};
  std::array<float, ModeCount> radius{};
  std::array<float, ModeCount> inputGain{};
  std::array<float, ModeCount> outputGain{};
  std::array<float, ModeCount> inputPhaseCosine{};
  std::array<float, ModeCount> inputPhaseSine{};
  std::array<float, ModeCount> exchangeCosine{};
  std::array<float, ModeCount> exchangeSine{};
  std::array<float, ModeCount> exchangeAmount{};
  std::array<std::uint32_t, ModeCount> sourceIndex{};
  std::array<std::uint16_t, ModeCount> packet{};
  std::array<std::uint8_t, ModeCount> band{};
  float sampleRate{48000.f};
  float lowCrossoverHz{700.f};
  float highCrossoverHz{6500.f};
  float maximumExchangeAngle{};
  std::uint32_t seed{0x4649454cu};
  std::uint32_t activeModeCount{};
};

namespace detail {

inline void NormalizeModalRotation(float &cosine, float &sine) noexcept {
  const float inverseLength = 1.f / std::sqrt(cosine * cosine + sine * sine);
  cosine *= inverseLength;
  sine *= inverseLength;
  if (std::abs(cosine) >= std::abs(sine)) {
    sine = std::copysign(
        std::sqrt(std::max(0.f, 1.f - cosine * cosine)), sine);
  } else {
    cosine = std::copysign(
        std::sqrt(std::max(0.f, 1.f - sine * sine)), cosine);
  }
}

inline float ModalInputGain(const float value) noexcept {
  return std::clamp(tfdsp::FiniteNormalOrZero(value), -256.f, 256.f);
}

inline float ModalOutputGain(const float value) noexcept {
  return std::clamp(tfdsp::FiniteNormalOrZero(value), -16.f, 16.f);
}

inline float UnitModalValue(const float value) noexcept {
  return std::clamp(tfdsp::FiniteNormalOrZero(value), 0.f, 1.f);
}

} // namespace detail

template <std::size_t ModeCount>
PreparedStochasticModalField<ModeCount> PrepareStochasticModalField(
    const float sampleRate,
    const std::array<StochasticModalModeParameters, ModeCount> &parameters,
    const StochasticModalFieldControls controls,
    const float lowCrossoverHz, const float highCrossoverHz) {
  if (!std::isfinite(sampleRate) || sampleRate < 1.f)
    throw std::invalid_argument("modal-field sample rate must be positive");
  if (!(lowCrossoverHz > 0.f && highCrossoverHz > lowCrossoverHz))
    throw std::invalid_argument("modal-field crossovers must be ordered");

  PreparedStochasticModalField<ModeCount> result{};
  result.sampleRate = sampleRate;
  result.lowCrossoverHz = lowCrossoverHz;
  result.highCrossoverHz = highCrossoverHz;
  result.maximumExchangeAngle = std::clamp(
      tfdsp::FiniteNormalOrZero(controls.exchangeAngleRadians), 0.f, .05f);
  result.seed = controls.seed;
  constexpr float TwoPi = 6.28318530717958647692f;
  for (std::size_t source = 0; source < ModeCount; ++source) {
    const float inputGain = detail::ModalInputGain(parameters[source].inputGain);
    const float outputGain =
        detail::ModalOutputGain(parameters[source].outputGain);
    if (inputGain == 0.f || outputGain == 0.f) continue;
    const std::size_t mode = result.activeModeCount++;
    result.sourceIndex[mode] = static_cast<std::uint32_t>(source);
    const float frequency = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters[source].frequencyHz), 1.f,
        .49f * sampleRate);
    const float decay = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters[source].decaySeconds), .001f,
        60.f);
    const float angle = TwoPi * frequency / sampleRate;
    const float oscillatorCosine = std::cos(angle);
    const float oscillatorSine = std::sin(angle);
    result.radius[mode] = std::exp(std::log(.001f) / (decay * sampleRate));
    result.inputGain[mode] = inputGain;
    result.outputGain[mode] = outputGain;
    const float phase = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters[source].inputPhaseRadians),
        -3.14159265358979323846f, 3.14159265358979323846f);
    result.inputPhaseCosine[mode] = std::cos(phase);
    result.inputPhaseSine[mode] = std::sin(phase);
    const float bandwidth = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters[source].phaseBandwidthHz), 0.f,
        .45f * sampleRate);
    const float phaseCosine = std::exp(
        -3.14159265358979323846f * bandwidth / sampleRate);
    const float phaseSine = std::sqrt(
        std::max(0.f, 1.f - phaseCosine * phaseCosine));
    float positiveCosine = phaseCosine * oscillatorCosine -
        phaseSine * oscillatorSine;
    float positiveSine = phaseSine * oscillatorCosine +
        phaseCosine * oscillatorSine;
    float negativeCosine = phaseCosine * oscillatorCosine +
        phaseSine * oscillatorSine;
    float negativeSine = phaseCosine * oscillatorSine -
        phaseSine * oscillatorCosine;
    detail::NormalizeModalRotation(positiveCosine, positiveSine);
    detail::NormalizeModalRotation(negativeCosine, negativeSine);
    result.cosineCentre[mode] = .5f * (positiveCosine + negativeCosine);
    result.cosineSpread[mode] = .5f * (positiveCosine - negativeCosine);
    result.sineCentre[mode] = .5f * (positiveSine + negativeSine);
    result.sineSpread[mode] = .5f * (positiveSine - negativeSine);
    result.packet[mode] = parameters[source].packet;
    result.exchangeAmount[mode] =
        detail::UnitModalValue(parameters[source].exchangeAmount);
    result.band[mode] = frequency < lowCrossoverHz ? 0 :
        static_cast<std::uint8_t>(frequency < highCrossoverHz ? 1 : 2);
  }
  result.exchangeCosine.fill(1.f);
  for (std::size_t first = 0; first + 1 < result.activeModeCount; ++first) {
    const auto distance = result.packet[first] > result.packet[first + 1]
        ? result.packet[first] - result.packet[first + 1]
        : result.packet[first + 1] - result.packet[first];
    if (distance > 1) continue;
    const float amount = std::sqrt(
        result.exchangeAmount[first] * result.exchangeAmount[first + 1]);
    const float exchangeAngle = result.maximumExchangeAngle * amount;
    result.exchangeCosine[first] = std::cos(exchangeAngle);
    result.exchangeSine[first] = std::sin(exchangeAngle);
  }
  return result;
}

} // namespace tfdsp::percussion
