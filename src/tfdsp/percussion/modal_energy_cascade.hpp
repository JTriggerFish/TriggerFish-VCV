#pragma once

#include "deterministic_random.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace tfdsp::percussion {

struct ModalEnergyCascadeParameters {
  float rateOctavesPerSecond{};
  float energyDependence{};
  float phaseDiffusion{};
  std::uint32_t seed{0x43415343u};
};

// Passive, one-way transport between frequency-ordered modal packets. Members
// of each packet are supplied as one contiguous run; the stochastic-field
// preparation API validates that contract. Energy
// already stored in a lower packet moves only into its next higher neighbour;
// processing pairs from high to low prevents a sample from crossing multiple
// packet boundaries. The transfer itself is energy preserving.
template <std::size_t ModeCount> class ModalEnergyCascade {
public:
  void Prepare(const float sampleRate,
               const std::array<float, ModeCount> &frequencyHz,
               const std::array<float, ModeCount> &inputGain,
               const std::array<std::uint16_t, ModeCount> &packet,
               const std::size_t activeModeCount,
               const ModalEnergyCascadeParameters parameters) noexcept {
    sampleRate_ = std::max(sampleRate, 1.f);
    activeModeCount_ = std::min(activeModeCount, ModeCount);
    rateOctavesPerSecond_ = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters.rateOctavesPerSecond), 0.f, 32.f);
    energyDependence_ = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters.energyDependence), 0.f, 1.f);
    phaseDiffusion_ = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters.phaseDiffusion), 0.f, 1.f);
    seed_ = parameters.seed;
    BuildPackets(frequencyHz, inputGain, packet);
    Reset();
  }

  void Reset() noexcept {
    random_.Seed(seed_);
    lastTransferredEnergy_ = 0.f;
  }

  float Process(std::array<float, ModeCount> &real,
                std::array<float, ModeCount> &imaginary) noexcept {
    lastTransferredEnergy_ = 0.f;
    if (rateOctavesPerSecond_ <= 0.f || packetCount_ < 2)
      return 0.f;
    for (std::size_t packet = 0; packet < packetCount_; ++packet) {
      originalEnergy_[packet] = Energy(upward_[packet], real, imaginary);
      finalEnergy_[packet] = originalEnergy_[packet];
      receivedFraction_[packet] = 0.f;
    }
    for (std::size_t upper = packetCount_ - 1; upper > 0; --upper)
      TransferEnergy(upper - 1, upper);
    for (std::size_t packet = 0; packet < packetCount_; ++packet)
      ApplyPacketUpdate(packet, real, imaginary);
    return lastTransferredEnergy_;
  }

  float LastTransferredEnergy() const noexcept {
    return lastTransferredEnergy_;
  }

private:
  struct Packet {
    std::size_t begin{};
    std::size_t end{};
    float centreFrequencyHz{1000.f};
    float referenceEnergy{1.f};
    float inverseGapToNext{};
  };

  void BuildPackets(const std::array<float, ModeCount> &frequencyHz,
                    const std::array<float, ModeCount> &inputGain,
                    const std::array<std::uint16_t, ModeCount> &packet) noexcept {
    packetCount_ = 0;
    std::size_t begin = 0;
    while (begin < activeModeCount_ && packetCount_ < ModeCount) {
      const auto id = packet[begin];
      std::size_t end = begin + 1;
      while (end < activeModeCount_ && packet[end] == id) ++end;
      double weightedFrequency = 0.0;
      double weight = 0.0;
      for (std::size_t mode = begin; mode < end; ++mode) {
        const double gainSquared = static_cast<double>(inputGain[mode]) *
            inputGain[mode];
        weightedFrequency += gainSquared * std::max(frequencyHz[mode], 1.f);
        weight += gainSquared;
      }
      upward_[packetCount_++] = {
          begin, end,
          static_cast<float>(weightedFrequency / std::max(weight, 1.e-20)),
          static_cast<float>(std::max(weight, 1.e-20))};
      begin = end;
    }
    for (std::size_t index = 1; index < packetCount_; ++index) {
      const Packet value = upward_[index];
      std::size_t position = index;
      while (position > 0 &&
             upward_[position - 1].centreFrequencyHz >
                 value.centreFrequencyHz) {
        upward_[position] = upward_[position - 1];
        --position;
      }
      upward_[position] = value;
    }
    for (std::size_t lower = 0; lower + 1 < packetCount_; ++lower) {
      const float gapOctaves = std::max(
          std::log2(upward_[lower + 1].centreFrequencyHz /
                    upward_[lower].centreFrequencyHz),
          .02f);
      upward_[lower].inverseGapToNext = 1.f / gapOctaves;
    }
  }

  static float Energy(const Packet &packet,
                      const std::array<float, ModeCount> &real,
                      const std::array<float, ModeCount> &imaginary) noexcept {
    double result = 0.0;
    for (std::size_t mode = packet.begin; mode < packet.end; ++mode)
      result += static_cast<double>(real[mode]) * real[mode] +
          static_cast<double>(imaginary[mode]) * imaginary[mode];
    return static_cast<float>(result);
  }

  void TransferEnergy(const std::size_t lowerIndex,
                      const std::size_t upperIndex) noexcept {
    const Packet &lower = upward_[lowerIndex];
    const float lowerEnergy = finalEnergy_[lowerIndex];
    if (!(lowerEnergy > 1.e-20f)) return;
    const float normalizedEnergy = lowerEnergy /
        (lowerEnergy + lower.referenceEnergy);
    const float activation = 1.f + energyDependence_ *
        (normalizedEnergy - 1.f);
    const float exponent = rateOctavesPerSecond_ * activation *
        lower.inverseGapToNext / sampleRate_;
    const float fraction = TransferFraction(exponent);
    const float transferred = fraction * lowerEnergy;
    if (!(transferred > 0.f)) return;
    finalEnergy_[lowerIndex] -= transferred;
    finalEnergy_[upperIndex] += transferred;
    receivedFraction_[upperIndex] = fraction;
    lastTransferredEnergy_ += transferred;
  }

  static float TransferFraction(const float exponent) noexcept {
    // [2/2] Pade form of 1-exp(-x). Below the .25 safety ceiling its
    // absolute error is below 2e-6, without a transcendental in the audio loop.
    constexpr float CeilingExponent = .2876820724517809f;
    if (exponent >= CeilingExponent) return .25f;
    const float positive = std::max(exponent, 0.f);
    return positive /
        (1.f + .5f * positive + positive * positive / 12.f);
  }

  void ApplyPacketUpdate(const std::size_t packetIndex,
                         std::array<float, ModeCount> &real,
                         std::array<float, ModeCount> &imaginary) noexcept {
    const Packet &packet = upward_[packetIndex];
    const float original = originalEnergy_[packetIndex];
    const float target = finalEnergy_[packetIndex];
    if (!(target > 1.e-20f)) return;
    if (!(original > 1.e-20f)) {
      SeedSilentPacket(packet, target, real, imaginary);
    } else {
      const float scale = std::sqrt(target / original);
      for (std::size_t mode = packet.begin; mode < packet.end; ++mode) {
        real[mode] *= scale;
        imaginary[mode] *= scale;
      }
    }
    DiffusePacket(packet, receivedFraction_[packetIndex], real, imaginary);
  }

  void SeedSilentPacket(const Packet &packet, const float energy,
                        std::array<float, ModeCount> &real,
                        std::array<float, ModeCount> &imaginary) noexcept {
    const float perMode = std::sqrt(
        energy / static_cast<float>(std::max<std::size_t>(
                     packet.end - packet.begin, 1)));
    for (std::size_t mode = packet.begin; mode < packet.end; ++mode) {
      real[mode] = random_.Bipolar() >= 0.f ? perMode : -perMode;
      imaginary[mode] = 0.f;
    }
  }

  void DiffusePacket(const Packet &packet, const float transferFraction,
                     std::array<float, ModeCount> &real,
                     std::array<float, ModeCount> &imaginary) noexcept {
    const float angle = 1.57079632679489661923f * phaseDiffusion_ *
        std::sqrt(transferFraction);
    if (!(angle > 0.f)) return;
    const float squared = angle * angle;
    // The maximum angle is pi/4. This fifth-order sine approximation remains
    // within 4e-5 there; deriving cosine restores exact unit magnitude.
    const float sine = angle *
        (1.f - squared / 6.f + squared * squared / 120.f);
    const float cosine = std::sqrt(std::max(0.f, 1.f - sine * sine));
    for (std::size_t mode = packet.begin; mode < packet.end; ++mode) {
      const float signedSine = random_.Bipolar() >= 0.f ? sine : -sine;
      const float oldReal = real[mode];
      real[mode] = cosine * oldReal - signedSine * imaginary[mode];
      imaginary[mode] = signedSine * oldReal + cosine * imaginary[mode];
    }
  }

  std::array<Packet, ModeCount> upward_{};
  std::array<float, ModeCount> originalEnergy_{};
  std::array<float, ModeCount> finalEnergy_{};
  std::array<float, ModeCount> receivedFraction_{};
  DeterministicRandom random_{};
  float sampleRate_{48000.f};
  float rateOctavesPerSecond_{};
  float energyDependence_{};
  float phaseDiffusion_{};
  float lastTransferredEnergy_{};
  std::uint32_t seed_{0x43415343u};
  std::size_t activeModeCount_{};
  std::size_t packetCount_{};
};

} // namespace tfdsp::percussion
