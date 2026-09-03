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
      const float oscillatorCosine = std::cos(angle);
      const float oscillatorSine = std::sin(angle);
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
      const float phaseCosine = std::exp(
          -3.14159265358979323846f * bandwidth / sampleRate_);
      const float phaseSine = std::sqrt(
          std::max(0.f, 1.f - phaseCosine * phaseCosine));
      // A stochastic phase kick follows the ordinary modal rotation. Compose
      // the two rotations here so the audio loop needs only one complex
      // multiply for either random sign.
      float positiveCosine = phaseCosine * oscillatorCosine -
          phaseSine * oscillatorSine;
      float positiveSine = phaseSine * oscillatorCosine +
          phaseCosine * oscillatorSine;
      float negativeCosine = phaseCosine * oscillatorCosine +
          phaseSine * oscillatorSine;
      float negativeSine = phaseCosine * oscillatorSine -
          phaseSine * oscillatorCosine;
      NormalizeRotation(positiveCosine, positiveSine);
      NormalizeRotation(negativeCosine, negativeSine);
      cosineCentre_[mode] = .5f * (positiveCosine + negativeCosine);
      cosineSpread_[mode] = .5f * (positiveCosine - negativeCosine);
      sineCentre_[mode] = .5f * (positiveSine + negativeSine);
      sineSpread_[mode] = .5f * (positiveSine - negativeSine);
      packet_[mode] = parameters[source].packet;
      exchangeAmount_[mode] = Unit(parameters[source].exchangeAmount);
      band_[mode] = frequency < lowCrossoverHz_
          ? 0
          : static_cast<std::uint8_t>(frequency < highCrossoverHz_ ? 1 : 2);
      effectiveRadius_[mode] = radius_[mode];
      primaryDriveGain_[mode] = inputGain * excitationProjection_[mode];
      secondaryDriveGain_[mode] =
          inputGain * secondaryExcitationProjection_[mode];
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
    SetProjection(projection, excitationProjection_, primaryDriveGain_);
  }

  void SetSecondaryExcitationProjection(
      const Projection &projection) noexcept {
    SetProjection(
        projection, secondaryExcitationProjection_, secondaryDriveGain_);
  }

  float ProcessExcitedPair(
      float primaryInput, float secondaryInput,
      const ModalDampingGains damping = {}) noexcept {
    primaryInput = tfdsp::FiniteNormalOrZero(primaryInput);
    secondaryInput = tfdsp::FiniteNormalOrZero(secondaryInput);
    UpdateDamping(damping);
    Propagate(primaryInput, secondaryInput);
    return ExchangeNeighboursAndSum();
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

  static void NormalizeRotation(float &cosine, float &sine) noexcept {
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

  void SetExchangeAngle(float angle) noexcept {
    angle = std::clamp(tfdsp::FiniteNormalOrZero(angle), 0.f, .05f);
    maximumExchangeAngle_ = angle;
  }

  void PrepareExchange() noexcept {
    exchangeCosine_.fill(1.f);
    exchangeSine_.fill(0.f);
    for (std::size_t first = 0; first + 1 < activeModeCount_; ++first) {
      const auto packetDistance = packet_[first] > packet_[first + 1]
          ? packet_[first] - packet_[first + 1]
          : packet_[first + 1] - packet_[first];
      if (packetDistance > 1)
        continue;
      const float amount = std::sqrt(
          exchangeAmount_[first] * exchangeAmount_[first + 1]);
      const float angle = maximumExchangeAngle_ * amount;
      exchangeCosine_[first] = std::cos(angle);
      exchangeSine_[first] = std::sin(angle);
    }
  }

  void SetProjection(const Projection &source, Projection &destination,
                     Projection &driveGain) noexcept {
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode) {
      destination[mode] = SafeProjection(source[sourceIndex_[mode]]);
      driveGain[mode] = inputGain_[mode] * destination[mode];
    }
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
    PrepareRandomSigns();
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode) {
      const float priorReal = real_[mode];
      const float priorImaginary = imaginary_[mode];
      const float cosine = cosineCentre_[mode] +
          randomSign_[mode] * cosineSpread_[mode];
      const float sine = sineCentre_[mode] +
          randomSign_[mode] * sineSpread_[mode];
      const float rotatedReal = cosine * priorReal - sine * priorImaginary;
      const float rotatedImaginary = sine * priorReal + cosine * priorImaginary;
      const float drive = primaryDriveGain_[mode] * primaryInput +
          secondaryDriveGain_[mode] * secondaryInput;
      real_[mode] = tfdsp::FiniteNormalOrZero(
          effectiveRadius_[mode] * rotatedReal +
          inputPhaseCosine_[mode] * drive);
      imaginary_[mode] = tfdsp::FiniteNormalOrZero(
          effectiveRadius_[mode] * rotatedImaginary +
          inputPhaseSine_[mode] * drive);
    }
  }

  void PrepareRandomSigns() noexcept {
    std::uint32_t randomBits = 0;
    unsigned randomBitsRemaining = 0;
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode) {
      if (randomBitsRemaining == 0) {
        randomBits = random_.NextBits();
        randomBitsRemaining = 32;
      }
      randomSign_[mode] = static_cast<float>(2 * (randomBits & 1u)) - 1.f;
      randomBits >>= 1;
      --randomBitsRemaining;
    }
  }

  float ExchangeNeighboursAndSum() noexcept {
    if (maximumExchangeAngle_ == 0.f || activeModeCount_ < 2)
      return SumOutput();
    const std::size_t offset = oddExchange_ ? 1 : 0;
    oddExchange_ = !oddExchange_;
    PrepareExchangeSigns(offset);
    float output = offset == 1 ? outputGain_[0] * real_[0] : 0.f;
    std::size_t first = offset;
    for (; first + 1 < activeModeCount_; first += 2) {
      const std::size_t second = first + 1;
      const float sine = randomSign_[first] * exchangeSine_[first];
      RotatePair(real_[first], real_[second], exchangeCosine_[first], sine);
      RotatePair(imaginary_[first], imaginary_[second],
                 exchangeCosine_[first], sine);
      output += outputGain_[first] * real_[first] +
          outputGain_[second] * real_[second];
    }
    if (first < activeModeCount_)
      output += outputGain_[first] * real_[first];
    return tfdsp::FiniteNormalOrZero(output);
  }

  void PrepareExchangeSigns(const std::size_t offset) noexcept {
    std::uint32_t randomBits = 0;
    unsigned randomBitsRemaining = 0;
    for (std::size_t first = offset; first + 1 < activeModeCount_; first += 2) {
      if (randomBitsRemaining == 0) {
        randomBits = random_.NextBits();
        randomBitsRemaining = 32;
      }
      randomSign_[first] = static_cast<float>(2 * (randomBits & 1u)) - 1.f;
      randomBits >>= 1;
      --randomBitsRemaining;
    }
  }

  static void RotatePair(float &first, float &second, const float cosine,
                         const float sine) noexcept {
    const float oldFirst = first;
    first = cosine * oldFirst - sine * second;
    second = sine * oldFirst + cosine * second;
  }

  float SumOutput() const noexcept {
    float output = 0.f;
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode)
      output += outputGain_[mode] * real_[mode];
    return tfdsp::FiniteNormalOrZero(output);
  }

  std::array<float, ModeCount> real_{};
  std::array<float, ModeCount> imaginary_{};
  std::array<float, ModeCount> cosineCentre_{};
  std::array<float, ModeCount> cosineSpread_{};
  std::array<float, ModeCount> sineCentre_{};
  std::array<float, ModeCount> sineSpread_{};
  std::array<float, ModeCount> randomSign_{};
  std::array<float, ModeCount> radius_{};
  std::array<float, ModeCount> effectiveRadius_{};
  std::array<float, ModeCount> inputGain_{};
  std::array<float, ModeCount> primaryDriveGain_{};
  std::array<float, ModeCount> secondaryDriveGain_{};
  std::array<float, ModeCount> outputGain_{};
  std::array<float, ModeCount> inputPhaseCosine_{};
  std::array<float, ModeCount> inputPhaseSine_{};
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
