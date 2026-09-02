#pragma once

#include "deterministic_random.hpp"
#include "modal_constraint.hpp"
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

// One bank spans coherent ridges, overlapping modal packets and stochastic
// wash. Random phase kicks and local Givens rotations preserve modal energy;
// only the declared pole radii and external constraints remove it.
template <std::size_t ModeCount> class StochasticModalField {
public:
  static_assert(ModeCount > 0, "a modal field needs at least one mode");
  using Parameters = std::array<StochasticModalModeParameters, ModeCount>;
  using Projection = std::array<float, ModeCount>;

  void Prepare(const float sampleRate, const Parameters &parameters,
               const StochasticModalFieldControls controls,
               const float lowCrossoverHz, const float highCrossoverHz) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("modal-field sample rate must be positive");
    if (!(lowCrossoverHz > 0.f && highCrossoverHz > lowCrossoverHz))
      throw std::invalid_argument("modal-field crossovers must be ordered");
    sampleRate_ = sampleRate;
    lowCrossoverHz_ = lowCrossoverHz;
    highCrossoverHz_ = highCrossoverHz;
    seed_ = controls.seed;
    SetExchangeAngle(controls.exchangeAngleRadians);
    excitationProjection_.fill(1.f);
    secondaryExcitationProjection_.fill(1.f);
    SetStaticParameters(parameters);
  }

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
      const float bandwidth = std::clamp(
          tfdsp::FiniteNormalOrZero(parameters[source].phaseBandwidthHz),
          0.f, .45f * sampleRate_);
      phaseCosine_[mode] = std::exp(-3.14159265358979323846f * bandwidth /
                                    sampleRate_);
      phaseSine_[mode] = std::sqrt(
          std::max(0.f, 1.f - phaseCosine_[mode] * phaseCosine_[mode]));
      packet_[mode] = parameters[source].packet;
      exchangeAmount_[mode] = Unit(parameters[source].exchangeAmount);
      band_[mode] = frequency < lowCrossoverHz_
          ? 0
          : static_cast<std::uint8_t>(frequency < highCrossoverHz_ ? 1 : 2);
      effectiveRadius_[mode] = radius_[mode];
    }
    PrepareExchange();
    damping_ = {};
    Reset();
  }

  void Reset() noexcept {
    real_.fill(0.f);
    imaginary_.fill(0.f);
    random_.Seed(seed_);
    oddExchange_ = false;
  }

  void SetExcitationProjection(const Projection &projection) noexcept {
    SetProjection(projection, excitationProjection_);
  }

  void SetSecondaryExcitationProjection(
      const Projection &projection) noexcept {
    SetProjection(projection, secondaryExcitationProjection_);
  }

  float ProcessExcitedPair(
      float primaryInput, float secondaryInput,
      const ModalDampingGains damping = {}) noexcept {
    primaryInput = tfdsp::FiniteNormalOrZero(primaryInput);
    secondaryInput = tfdsp::FiniteNormalOrZero(secondaryInput);
    UpdateDamping(damping);
    Propagate(primaryInput, secondaryInput);
    ExchangeNeighbours();
    return SumOutput();
  }

  std::size_t ActiveModeCount() const noexcept { return activeModeCount_; }

  double StoredEnergy() const noexcept {
    double result = 0.0;
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode)
      result += static_cast<double>(real_[mode]) * real_[mode] +
          static_cast<double>(imaginary_[mode]) * imaginary_[mode];
    return result;
  }

private:
  static float SafeInputGain(const float gain) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(gain), -256.f, 256.f);
  }

  static float SafeOutputGain(const float gain) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(gain), -16.f, 16.f);
  }

  static float SafeProjection(const float value) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(value), -4.f, 4.f);
  }

  static float Unit(const float value) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(value), 0.f, 1.f);
  }

  void SetExchangeAngle(float angle) noexcept {
    angle = std::clamp(tfdsp::FiniteNormalOrZero(angle), 0.f, .05f);
    maximumExchangeAngle_ = angle;
  }

  void PrepareExchange() noexcept {
    exchangeCosine_.fill(1.f);
    exchangeSine_.fill(0.f);
    for (std::size_t first = 0; first + 1 < activeModeCount_; ++first) {
      const float amount = std::sqrt(
          exchangeAmount_[first] * exchangeAmount_[first + 1]);
      const float angle = maximumExchangeAngle_ * amount;
      exchangeCosine_[first] = std::cos(angle);
      exchangeSine_[first] = std::sin(angle);
    }
  }

  void SetProjection(const Projection &source, Projection &destination) noexcept {
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode)
      destination[mode] = SafeProjection(source[sourceIndex_[mode]]);
  }

  void UpdateDamping(const ModalDampingGains damping) noexcept {
    const ModalDampingGains safe{
        Unit(damping.broadband), Unit(damping.low),
        Unit(damping.middle), Unit(damping.high)};
    if (safe.broadband == damping_.broadband && safe.low == damping_.low &&
        safe.middle == damping_.middle && safe.high == damping_.high)
      return;
    const std::array<float, 3> bandGain{safe.low, safe.middle, safe.high};
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode)
      effectiveRadius_[mode] =
          radius_[mode] * safe.broadband * bandGain[band_[mode]];
    damping_ = safe;
  }

  void Propagate(const float primaryInput,
                 const float secondaryInput) noexcept {
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode) {
      const float priorReal = real_[mode];
      const float priorImaginary = imaginary_[mode];
      float rotatedReal = cosine_[mode] * priorReal -
          sine_[mode] * priorImaginary;
      float rotatedImaginary = sine_[mode] * priorReal +
          cosine_[mode] * priorImaginary;
      if (phaseSine_[mode] > 0.f) {
        const float randomSine = (random_.NextBits() & 1u)
            ? phaseSine_[mode] : -phaseSine_[mode];
        const float diffusedReal = phaseCosine_[mode] * rotatedReal -
            randomSine * rotatedImaginary;
        rotatedImaginary = randomSine * rotatedReal +
            phaseCosine_[mode] * rotatedImaginary;
        rotatedReal = diffusedReal;
      }
      const float drive = inputGain_[mode] *
          (excitationProjection_[mode] * primaryInput +
           secondaryExcitationProjection_[mode] * secondaryInput);
      real_[mode] = tfdsp::FiniteNormalOrZero(
          effectiveRadius_[mode] * rotatedReal +
          inputPhaseCosine_[mode] * drive);
      imaginary_[mode] = tfdsp::FiniteNormalOrZero(
          effectiveRadius_[mode] * rotatedImaginary +
          inputPhaseSine_[mode] * drive);
    }
  }

  void ExchangeNeighbours() noexcept {
    if (maximumExchangeAngle_ == 0.f || activeModeCount_ < 2)
      return;
    const std::size_t offset = oddExchange_ ? 1 : 0;
    oddExchange_ = !oddExchange_;
    for (std::size_t first = offset; first + 1 < activeModeCount_; first += 2) {
      const std::size_t second = first + 1;
      const auto packetDistance = packet_[first] > packet_[second]
          ? packet_[first] - packet_[second]
          : packet_[second] - packet_[first];
      if (packetDistance > 1)
        continue;
      const float sine = (random_.NextBits() & 1u)
          ? exchangeSine_[first] : -exchangeSine_[first];
      RotatePair(real_[first], real_[second], exchangeCosine_[first], sine);
      RotatePair(imaginary_[first], imaginary_[second],
                 exchangeCosine_[first], sine);
    }
  }

  static void RotatePair(float &first, float &second, const float cosine,
                         const float sine) noexcept {
    const float oldFirst = first;
    first = tfdsp::FiniteNormalOrZero(
        cosine * oldFirst - sine * second);
    second = tfdsp::FiniteNormalOrZero(
        sine * oldFirst + cosine * second);
  }

  float SumOutput() const noexcept {
    float output = 0.f;
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode)
      output += outputGain_[mode] * real_[mode];
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
  std::array<float, ModeCount> phaseCosine_{};
  std::array<float, ModeCount> phaseSine_{};
  std::array<std::uint8_t, ModeCount> band_{};
  std::array<std::uint16_t, ModeCount> packet_{};
  std::array<float, ModeCount> exchangeAmount_{};
  std::array<float, ModeCount> exchangeCosine_{};
  std::array<float, ModeCount> exchangeSine_{};
  std::array<std::size_t, ModeCount> sourceIndex_{};
  Projection excitationProjection_{};
  Projection secondaryExcitationProjection_{};
  DeterministicRandom random_{};
  ModalDampingGains damping_{};
  float sampleRate_{48000.f};
  float lowCrossoverHz_{700.f};
  float highCrossoverHz_{6500.f};
  float maximumExchangeAngle_{};
  std::uint32_t seed_{0x4649454cu};
  std::size_t activeModeCount_{};
  bool oddExchange_{};
};

} // namespace tfdsp::percussion
