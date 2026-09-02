#pragma once

#include "deterministic_random.hpp"
#include "modal_bank.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace tfdsp::percussion {

struct StatisticalModalCloudParameters {
  static constexpr std::size_t GainEnvelopePointCount = 33;
  static constexpr std::size_t DecayEnvelopePointCount = 6;
  float minimumFrequencyHz{650.f};
  float maximumFrequencyHz{18000.f};
  float frequencyWarp{1.f};
  float spacingJitter{.8f};
  float modeDensity{1.f};
  float lowDecaySeconds{3.f};
  float highDecaySeconds{.18f};
  float decayCurve{.8f};
  std::array<float, DecayEnvelopePointCount> decayEnvelopeOctaves{};
  float decaySpreadOctaves{.35f};
  float tiltDbPerOctave{-1.5f};
  std::array<float, GainEnvelopePointCount> gainEnvelopeDb{};
  float gainSpreadDb{4.f};
  float outputGain{1.f};
  std::uint32_t seed{0x43594d42u};
};

template <std::size_t PointCount>
inline float InterpolateEnvelope(const std::array<float, PointCount> &points,
    const float position) noexcept {
  static_assert(PointCount >= 2, "an envelope needs at least two points");
  constexpr std::size_t Last = PointCount - 1;
  const float coordinate = std::clamp(position, 0.f, 1.f) * Last;
  const auto left = static_cast<std::size_t>(coordinate);
  const auto right = std::min(left + 1, Last);
  const float amount = coordinate - static_cast<float>(left);
  return points[left] + amount * (points[right] - points[left]);
}

inline float ErbRate(const float frequencyHz) noexcept {
  return 21.4f * std::log10(1.f + .00437f * frequencyHz);
}

inline float InverseErbRate(const float erbRate) noexcept {
  return (std::pow(10.f, erbRate / 21.4f) - 1.f) / .00437f;
}

inline float GeometricMix(float first, float second,
                          const float amount) noexcept {
  first = std::max(first, .001f);
  second = std::max(second, .001f);
  return std::exp(std::log(first) + amount * (std::log(second) - std::log(first)));
}

// Builds a repeatable high-density modal population from a small set of smooth
// statistical controls. Individual frequencies are nuisance variables rather
// than calibration controls.
template <std::size_t ModeCount>
typename ModalBank<ModeCount>::Parameters MakeStatisticalModalCloud(
    const float sampleRate, const StatisticalModalCloudParameters &parameters) {
  using Bank = ModalBank<ModeCount>;
  typename Bank::Parameters result{};
  DeterministicRandom random;
  random.Seed(parameters.seed);

  const float minimum = std::clamp(parameters.minimumFrequencyHz, 20.f,
                                   .45f * sampleRate);
  const float maximum = std::clamp(parameters.maximumFrequencyHz,
                                   minimum + 1.f, .48f * sampleRate);
  const float minimumErb = ErbRate(minimum);
  const float maximumErb = ErbRate(maximum);
  const float warp = std::clamp(parameters.frequencyWarp, .25f, 4.f);
  const float jitter = std::clamp(parameters.spacingJitter, 0.f, .95f);
  std::array<float, ModeCount> rawGain{};

  for (std::size_t mode = 0; mode < ModeCount; ++mode) {
    const float cell = (static_cast<float>(mode) + .5f +
                        .48f * jitter * random.Bipolar()) /
        static_cast<float>(ModeCount);
    const float position = std::pow(std::clamp(cell, 0.f, 1.f), warp);
    const float frequency = InverseErbRate(
        minimumErb + position * (maximumErb - minimumErb));
    const float decayPosition = std::pow(
        std::clamp((frequency - minimum) / (maximum - minimum), 0.f, 1.f),
        std::clamp(parameters.decayCurve, .1f, 4.f));
    const float baseDecay = GeometricMix(parameters.lowDecaySeconds,
                                         parameters.highDecaySeconds,
                                         decayPosition);
    const float decayShape = std::clamp(
        InterpolateEnvelope(parameters.decayEnvelopeOctaves, position),
        -4.f, 4.f);
    const float decay = baseDecay * std::exp2(
        decayShape +
        std::clamp(parameters.decaySpreadOctaves, 0.f, 2.f) * random.Bipolar());
    const float octaves = std::log2(frequency / 1000.f);
    const float tiltGain = std::pow(
        10.f, parameters.tiltDbPerOctave * octaves / 20.f);
    const float envelopeGain = std::pow(
        10.f, std::clamp(InterpolateEnvelope(parameters.gainEnvelopeDb, position),
                         -24.f, 24.f) /
                  20.f);
    const float variation = std::pow(
        10.f, std::clamp(parameters.gainSpreadDb, 0.f, 18.f) *
                  random.Bipolar() / 20.f);
    rawGain[mode] = tiltGain * envelopeGain * variation;
    result[mode] = {frequency, decay, 1.f, 0.f, 0.f};
  }

  // Density activates an exact, nested subset. Raising it only adds modes;
  // existing modal frequencies and nuisance parameters never move.
  DeterministicRandom densityRandom;
  densityRandom.Seed(parameters.seed ^ 0x44454e53u);
  const float density = std::clamp(parameters.modeDensity, 0.f, 1.f);
  const float requestedCount = density * static_cast<float>(ModeCount);
  const std::size_t fullCount = std::min<std::size_t>(
      static_cast<std::size_t>(requestedCount), ModeCount);
  const float fractional = fullCount < ModeCount
      ? requestedCount - static_cast<float>(fullCount) : 0.f;
  std::array<std::pair<float, std::size_t>, ModeCount> priority{};
  for (std::size_t mode = 0; mode < ModeCount; ++mode)
    priority[mode] = {densityRandom.Uniform(), mode};
  std::sort(priority.begin(), priority.end(), [](const auto &left,
                                                  const auto &right) {
    return left.first < right.first ||
        (left.first == right.first && left.second < right.second);
  });
  std::array<float, ModeCount> activation{};
  for (std::size_t rank = 0; rank < fullCount; ++rank)
    activation[priority[rank].second] = 1.f;
  if (fractional > 0.f)
    activation[priority[fullCount].second] = fractional;
  float squaredGain = 0.f;
  for (std::size_t mode = 0; mode < ModeCount; ++mode) {
    rawGain[mode] *= activation[mode];
    squaredGain += rawGain[mode] * rawGain[mode];
  }

  const float normalization = parameters.outputGain /
      std::sqrt(std::max(squaredGain, 1.e-12f));
  random.Seed(parameters.seed ^ 0x4f555450u);
  constexpr float Pi = 3.14159265358979323846f;
  for (std::size_t mode = 0; mode < ModeCount; ++mode) {
    result[mode].outputGain = normalization * rawGain[mode];
    result[mode].inputPhaseRadians = Pi * random.Bipolar();
  }
  return result;
}

} // namespace tfdsp::percussion
