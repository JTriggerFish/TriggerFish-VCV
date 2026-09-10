#pragma once

#include "deterministic_random.hpp"
#include "modal_spectral_diffusion.hpp"
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
  bool spectralDiffusion{};
  float energySensitivity{}; // Total-energy exponent for spectral diffusion only.
};

// Passive transport between frequency-ordered modal packets: the original
// one-way stencil or experimental nonlinear spectral diffusion. Members
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
    spectralDiffusion_ = parameters.spectralDiffusion;
    energySensitivity_ = std::clamp(
        tfdsp::FiniteNormalOrZero(parameters.energySensitivity), 0.f, 2.f);
    BuildPackets(frequencyHz, inputGain, packet);
    if (spectralDiffusion_) PrepareDiffusion(inputGain);
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
    if (spectralDiffusion_) {
      DiffuseEnergy();
    } else {
      const float activation = EnergyActivation();
      for (std::size_t source = 0; source + 1 < packetCount_; ++source)
        TransferEnergy(source, activation);
    }
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

  void PrepareDiffusion(const std::array<float, ModeCount> &inputGain) noexcept {
    std::array<float, ModeCount> centres{};
    std::array<double, ModeCount> weights{};
    for (std::size_t i = 0; i < packetCount_; ++i) {
      centres[i] = upward_[i].centreFrequencyHz;
      for (std::size_t mode = upward_[i].begin; mode < upward_[i].end; ++mode)
        weights[i] += double(inputGain[mode]) * inputGain[mode];
    }
    diffusion_.Prepare(centres, weights, packetCount_, sampleRate_);
  }

  void DiffuseEnergy() noexcept {
    diffusion_.Process(originalEnergy_, finalEnergy_, totalReferenceEnergy_,
        double(rateOctavesPerSecond_) / sampleRate_, energyAcceleration_, energySensitivity_);
    for (std::size_t i = 0; i < packetCount_; ++i) {
      const float arrival = std::max(0.f, finalEnergy_[i] - originalEnergy_[i]);
      receivedFraction_[i] = arrival / std::max(finalEnergy_[i], 1.e-20f);
      // Net relocated energy, not the sum of gross internal face fluxes.
      lastTransferredEnergy_ += arrival;
    }
  }

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
          weight > 0 ? static_cast<float>(weightedFrequency / weight)
                     : std::max(frequencyHz[begin], 1.f)};
      totalReferenceEnergy_ += static_cast<float>(weight);
      for (std::size_t mode = begin; mode < end; ++mode) {
        const double gain = inputGain[mode];
        seedEnergyWeight_[mode] = weight > 0
            ? static_cast<float>(gain * gain / weight) : 0.f;
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
    if (!(target > 0.f)) {
      for (std::size_t mode = packet.begin; mode < packet.end; ++mode)
        real[mode] = imaginary[mode] = 0.f;
      return;
    }
    if (!(original > 0.f)) {
      SeedSilentPacket(packet, target, real, imaginary);
    } else {
      // Preserve quiet stored states and their phase; only zero energy needs
      // seeding. Double division avoids overflow for a tiny receiving state.
      const float scale = static_cast<float>(std::sqrt(double(target) / original));
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
    // Legacy arrivals are at most pi/4; their fifth-order sine approximation
    // stays within 4e-5. Diffusion may fill a silent cell in one step, requiring
    // angles up to pi/2 and the exact sine. Cosine restores unit magnitude.
    const float sine = transferFraction <= .25f ? angle *
        (1.f - squared / 6.f + squared * squared / 120.f) : std::sin(angle);
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
  ModalSpectralDiffusion<ModeCount> diffusion_{};
  bool spectralDiffusion_{};
  float sampleRate_{48000.f};
  float rateOctavesPerSecond_{};
  float energyAcceleration_{};
  float energySensitivity_{};
  float phaseDiffusion_{};
  float totalReferenceEnergy_{1.f};
  float lastTransferredEnergy_{};
  std::uint32_t seed_{0x43415343u};
  std::size_t activeModeCount_{};
  std::size_t packetCount_{};
};

} // namespace tfdsp::percussion
