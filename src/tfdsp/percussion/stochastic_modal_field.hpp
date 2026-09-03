#pragma once

#include "deterministic_random.hpp"
#include "modal_constraint.hpp"
#include "stochastic_modal_field_parameters.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace tfdsp::percussion {

// One bank spans coherent ridges, overlapping modal packets and stochastic
// wash. Random phase kicks and local Givens rotations preserve modal energy;
// only the declared pole radii and external constraints remove it.
template <std::size_t ModeCount> class StochasticModalField {
public:
  static_assert(ModeCount > 0, "a modal field needs at least one mode");
  using Parameters = std::array<StochasticModalModeParameters, ModeCount>;
  using Projection = std::array<float, ModeCount>;
  using PreparedParameters = PreparedStochasticModalField<ModeCount>;

  void Prepare(const float sampleRate, const Parameters &parameters,
               const StochasticModalFieldControls controls,
               const float lowCrossoverHz, const float highCrossoverHz) {
    LoadPrepared(PrepareStochasticModalField(
        sampleRate, parameters, controls, lowCrossoverHz, highCrossoverHz));
  }

  void SetStaticParameters(const Parameters &parameters) {
    const auto primaryProjection = excitationProjection_;
    const auto secondaryProjection = secondaryExcitationProjection_;
    LoadPrepared(PrepareStochasticModalField(
        sampleRate_, parameters,
        {maximumExchangeAngle_, seed_}, lowCrossoverHz_, highCrossoverHz_));
    SetProjection(
        primaryProjection, excitationProjection_, primaryDriveGain_);
    SetProjection(secondaryProjection, secondaryExcitationProjection_,
                  secondaryDriveGain_);
  }

  void LoadPrepared(const PreparedParameters &prepared) noexcept {
    sampleRate_ = prepared.sampleRate;
    lowCrossoverHz_ = prepared.lowCrossoverHz;
    highCrossoverHz_ = prepared.highCrossoverHz;
    maximumExchangeAngle_ = prepared.maximumExchangeAngle;
    seed_ = prepared.seed;
    activeModeCount_ = std::min<std::size_t>(
        prepared.activeModeCount, ModeCount);
    cosineCentre_ = prepared.cosineCentre;
    cosineSpread_ = prepared.cosineSpread;
    sineCentre_ = prepared.sineCentre;
    sineSpread_ = prepared.sineSpread;
    radius_ = prepared.radius;
    inputGain_ = prepared.inputGain;
    outputGain_ = prepared.outputGain;
    inputPhaseCosine_ = prepared.inputPhaseCosine;
    inputPhaseSine_ = prepared.inputPhaseSine;
    exchangeCosine_ = prepared.exchangeCosine;
    exchangeSine_ = prepared.exchangeSine;
    exchangeAmount_ = prepared.exchangeAmount;
    sourceIndex_ = prepared.sourceIndex;
    packet_ = prepared.packet;
    band_ = prepared.band;
    excitationProjection_.fill(1.f);
    secondaryExcitationProjection_.fill(1.f);
    for (std::size_t mode = 0; mode < activeModeCount_; ++mode) {
      effectiveRadius_[mode] = radius_[mode];
      primaryDriveGain_[mode] = inputGain_[mode];
      secondaryDriveGain_[mode] = inputGain_[mode];
    }
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
  static float SafeProjection(const float value) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(value), -4.f, 4.f);
  }

  static float Unit(const float value) noexcept {
    return std::clamp(tfdsp::FiniteNormalOrZero(value), 0.f, 1.f);
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
  std::array<std::uint32_t, ModeCount> sourceIndex_{};
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
