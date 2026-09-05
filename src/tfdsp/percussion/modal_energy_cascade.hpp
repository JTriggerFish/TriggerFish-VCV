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
  float energyAcceleration{};
  float phaseDiffusion{};
  std::uint32_t seed{0x43415343u};
};

// Passive, one-way transport between frequency-ordered modal packets. Members
// of each packet are supplied as one contiguous run; the stochastic-field
// preparation API validates that contract. A fixed half-octave transport
// stencil is interpolated onto the painted packets, so adding intermediate
// anchors refines the spectrum without adding serial transport stages. All
// transfers use the states measured at the start of the sample and preserve
// their total energy.
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
    energyAcceleration_ = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters.energyAcceleration), 0.f, 1.f);
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
    const float activation = EnergyActivation();
    for (std::size_t source = 0; source + 1 < packetCount_; ++source)
      TransferEnergy(source, activation);
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
  };

  struct Destination {
    std::size_t packet{};
    float weight{};
  };

  struct Route {
    Destination lower{};
    Destination upper{};
    float inverseDistanceOctaves{};
  };

  static constexpr float TransportStepOctaves = .5f;

  void BuildPackets(const std::array<float, ModeCount> &frequencyHz,
                    const std::array<float, ModeCount> &inputGain,
                    const std::array<std::uint16_t, ModeCount> &packet) noexcept {
    packetCount_ = 0;
    totalReferenceEnergy_ = 0.f;
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
          static_cast<float>(weightedFrequency / std::max(weight, 1.e-20))};
      totalReferenceEnergy_ += static_cast<float>(weight);
      for (std::size_t mode = begin; mode < end; ++mode) {
        const float gain = inputGain[mode];
        seedEnergyWeight_[mode] = gain * gain /
            static_cast<float>(std::max(weight, 1.e-20));
      }
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
    BuildRoutes();
  }

  void BuildRoutes() noexcept {
    if (packetCount_ < 2) return;
    const float highest = std::log2(
        std::max(upward_[packetCount_ - 1].centreFrequencyHz, 1.f));
    for (std::size_t source = 0; source + 1 < packetCount_; ++source) {
      const float origin = std::log2(
          std::max(upward_[source].centreFrequencyHz, 1.f));
      const float target = std::min(origin + TransportStepOctaves, highest);
      std::size_t above = source + 1;
      while (above + 1 < packetCount_ &&
             std::log2(std::max(upward_[above].centreFrequencyHz, 1.f)) <
                 target)
        ++above;
      const float abovePosition = std::log2(
          std::max(upward_[above].centreFrequencyHz, 1.f));
      if (abovePosition <= target) {
        routes_[source] = {{above, 1.f}, {},
            1.f / std::max(target - origin, .02f)};
        continue;
      }
      const std::size_t below = above - 1;
      const float belowPosition = std::log2(
          std::max(upward_[below].centreFrequencyHz, 1.f));
      const float upperWeight = std::clamp(
          (target - belowPosition) /
              std::max(abovePosition - belowPosition, .02f),
          0.f, 1.f);
      routes_[source] = {{below, 1.f - upperWeight}, {above, upperWeight},
          1.f / std::max(target - origin, .02f)};
    }
  }

  float EnergyActivation() const noexcept {
    double totalEnergy = 0.0;
    for (std::size_t packet = 0; packet < packetCount_; ++packet)
      totalEnergy += originalEnergy_[packet];
    const float normalizedEnergy = static_cast<float>(totalEnergy /
        (totalEnergy + std::max(totalReferenceEnergy_, 1.e-20f)));
    // The rate control is the baseline transport speed. Energy acceleration is
    // an orthogonal increase above that rate; field-wide normalization keeps
    // it independent of how many handles and sidebands represent the body.
    // One unit spans a useful 1x..8x range. The previous 1x..2x mapping was
    // too weak to separate a fast high-energy bloom from a slow late tail.
    return 1.f + 7.f * energyAcceleration_ * normalizedEnergy;
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

  void TransferEnergy(const std::size_t source,
                      const float activation) noexcept {
    const float sourceEnergy = originalEnergy_[source];
    if (!(sourceEnergy > 1.e-20f)) return;
    const Route &route = routes_[source];
    const float exponent = rateOctavesPerSecond_ * activation *
        route.inverseDistanceOctaves / sampleRate_;
    const float fraction = TransferFraction(exponent);
    const float eventEnergy = fraction * sourceEnergy;
    if (!(eventEnergy > 0.f)) return;
    finalEnergy_[source] -= eventEnergy;
    AddArrival(source, route.lower, eventEnergy, fraction);
    AddArrival(source, route.upper, eventEnergy, fraction);
  }

  void AddArrival(const std::size_t source, const Destination destination,
                  const float eventEnergy, const float fraction) noexcept {
    if (!(destination.weight > 0.f)) return;
    const float arrival = destination.weight * eventEnergy;
    finalEnergy_[destination.packet] += arrival;
    if (destination.packet == source) return;
    receivedFraction_[destination.packet] = std::max(
        receivedFraction_[destination.packet],
        destination.weight * fraction);
    lastTransferredEnergy_ += arrival;
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
    for (std::size_t mode = packet.begin; mode < packet.end; ++mode) {
      const float magnitude = std::sqrt(
          energy * seedEnergyWeight_[mode]);
      real[mode] = random_.Bipolar() >= 0.f ? magnitude : -magnitude;
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
  std::array<Route, ModeCount> routes_{};
  std::array<float, ModeCount> originalEnergy_{};
  std::array<float, ModeCount> finalEnergy_{};
  std::array<float, ModeCount> receivedFraction_{};
  std::array<float, ModeCount> seedEnergyWeight_{};
  DeterministicRandom random_{};
  float sampleRate_{48000.f};
  float rateOctavesPerSecond_{};
  float energyAcceleration_{};
  float phaseDiffusion_{};
  float totalReferenceEnergy_{1.f};
  float lastTransferredEnergy_{};
  std::uint32_t seed_{0x43415343u};
  std::size_t activeModeCount_{};
  std::size_t packetCount_{};
};

} // namespace tfdsp::percussion
