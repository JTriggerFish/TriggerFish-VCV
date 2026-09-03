#include "wire_rack.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tfdsp::percussion {
namespace {

float Safe(float value, float fallback, float low, float high) noexcept {
  return std::clamp(std::isfinite(value) ? value : fallback, low, high);
}

} // namespace

WireRackPreparedParameters
PrepareWireRackParameters(const float sampleRate,
                          const WireRackParameters &source) {
  if (!std::isfinite(sampleRate) || sampleRate < 1.f)
    throw std::invalid_argument("wire-rack sample rate must be positive");
  WireRackPreparedParameters result;
  result.sampleRate = sampleRate;
  const float low =
      Safe(source.minimumFrequencyHz, 900.f, 40.f, .4f * sampleRate);
  const float high =
      Safe(source.maximumFrequencyHz, 15500.f, low + 1.f, .48f * sampleRate);
  const float density = Safe(source.density, .8f, 0.f, 1.f);
  result.activeModeCount =
      std::clamp<std::size_t>(static_cast<std::size_t>(std::lround(
                                  8.f + density * (WireRackModeCount - 8.f))),
                              8, WireRackModeCount);
  const float decay = Safe(source.decaySeconds, .16f, .008f, 3.f);
  const float decayTilt = Safe(source.decayTilt, .7f, -1.f, 1.f);
  DeterministicRandom random;
  random.Seed(source.seed);
  float outputNormSquared = 0.f;
  constexpr float TwoPi = 6.28318530717958647692f;
  for (std::size_t mode = 0; mode < result.activeModeCount; ++mode) {
    const float position =
        (static_cast<float>(mode) + .5f + .38f * random.Bipolar()) /
        static_cast<float>(result.activeModeCount);
    const float frequency =
        low * std::pow(high / low, std::clamp(position, 0.f, 1.f));
    const float angle = TwoPi * frequency / sampleRate;
    const float modeDecay =
        decay * std::exp2(-2.f * decayTilt * std::clamp(position, 0.f, 1.f));
    const float phase = 3.14159265358979323846f * random.Bipolar();
    const float gain =
        std::exp2(-.35f * position) * (.75f + .5f * random.Uniform());
    result.cosine[mode] = std::cos(angle);
    result.sine[mode] = std::sin(angle);
    result.radius[mode] = std::exp(std::log(.001f) / (modeDecay * sampleRate));
    result.inputPhaseCosine[mode] = std::cos(phase);
    result.inputPhaseSine[mode] = std::sin(phase);
    result.modeOutputGain[mode] = gain;
    outputNormSquared += gain * gain;
  }
  const float outputScale =
      1.25f / std::sqrt(std::max(outputNormSquared, 1.e-12f));
  for (std::size_t mode = 0; mode < result.activeModeCount; ++mode)
    result.modeOutputGain[mode] *= outputScale;
  const float motionHz =
      Safe(source.motionHighpassHz, 140.f, 5.f, .45f * sampleRate);
  result.motionCoefficient = 1.f - std::exp(-TwoPi * motionHz / sampleRate);
  const float attack = Safe(source.attackSeconds, .002f, .0001f, .08f);
  result.attackCoefficient = std::exp(-1.f / (attack * sampleRate));
  const float release = Safe(source.releaseSeconds, .018f, .0005f, 1.f);
  result.releaseCoefficient = std::exp(-1.f / (release * sampleRate));
  result.sensitivity = Safe(source.sensitivity, 9.f, 0.f, 64.f);
  result.threshold = Safe(source.threshold, .004f, 0.f, 1.f);
  result.noiseTiltDb = 16.f * (Safe(source.brightness, .62f, 0.f, 1.f) - .5f);
  result.noiseMix = Safe(source.noiseMix, .6f, 0.f, 2.f);
  result.modalMix = Safe(source.modalMix, .75f, 0.f, 2.f);
  result.outputLevel = Safe(source.outputGain, .42f, 0.f, 4.f);
  result.maximumModalEnergy = Safe(source.maximumModalEnergy, 1.f, .001f, 16.f);
  result.seed = source.seed;
  return result;
}

void WireRack::Prepare(const float sampleRate,
                       const WireRackParameters &parameters) {
  Prepare(PrepareWireRackParameters(sampleRate, parameters));
}

void WireRack::Prepare(const WireRackPreparedParameters &parameters) {
  parameters_ = parameters;
  tilt_.Prepare(parameters.sampleRate);
  tilt_.SetTilt(parameters.noiseTiltDb, 3500.f);
  Reset();
}

void WireRack::Reset() noexcept {
  real_.fill(0.f);
  imaginary_.fill(0.f);
  motionLowpass_ = contactEnvelope_ = 0.f;
  tilt_.Reset();
  random_.Seed(parameters_.seed);
}

void WireRack::Seed(const std::uint32_t seed) noexcept {
  random_.Seed(seed == 0 ? parameters_.seed : seed);
}

float WireRack::Process(float bodyMotion) noexcept {
  bodyMotion = tfdsp::FiniteNormalOrZero(bodyMotion);
  motionLowpass_ +=
      parameters_.motionCoefficient * (bodyMotion - motionLowpass_);
  motionLowpass_ = tfdsp::FiniteNormalOrZero(motionLowpass_);
  const float motion = bodyMotion - motionLowpass_;
  const float target = std::max(0.f, std::abs(motion) - parameters_.threshold);
  contactEnvelope_ = target > contactEnvelope_
                         ? parameters_.attackCoefficient * contactEnvelope_ +
                               (1.f - parameters_.attackCoefficient) * target
                         : parameters_.releaseCoefficient * contactEnvelope_;
  contactEnvelope_ = tfdsp::FiniteNormalOrZero(contactEnvelope_);
  const float linearContact =
      std::clamp(parameters_.sensitivity * contactEnvelope_, 0.f, 1.f);
  // Wire contact grows with both the number of touching strands and their
  // individual force. Squaring the normalized follower gives a smooth onset
  // without introducing a second trigger or a delayed noise burst.
  const float contact = linearContact * linearContact;
  const float noise = tilt_.Process(1.7320508075688772f * random_.Bipolar());
  const float drive = contact * noise;

  float baseEnergy = 0.f;
  float crossEnergy = 0.f;
  float driveEnergy = 0.f;
  const float inputGain =
      1.f / std::sqrt(static_cast<float>(parameters_.activeModeCount));
  for (std::size_t mode = 0; mode < parameters_.activeModeCount; ++mode) {
    const float priorReal = real_[mode];
    const float priorImaginary = imaginary_[mode];
    real_[mode] =
        parameters_.radius[mode] * (parameters_.cosine[mode] * priorReal -
                                    parameters_.sine[mode] * priorImaginary);
    imaginary_[mode] =
        parameters_.radius[mode] * (parameters_.sine[mode] * priorReal +
                                    parameters_.cosine[mode] * priorImaginary);
    const float force = inputGain * drive;
    baseEnergy +=
        real_[mode] * real_[mode] + imaginary_[mode] * imaginary_[mode];
    crossEnergy +=
        force * (parameters_.inputPhaseCosine[mode] * real_[mode] +
                 parameters_.inputPhaseSine[mode] * imaginary_[mode]);
    driveEnergy += force * force;
  }
  const float driveScale =
      AvailableDriveScale(baseEnergy, crossEnergy, driveEnergy);
  float modal = 0.f;
  for (std::size_t mode = 0; mode < parameters_.activeModeCount; ++mode) {
    const float force = driveScale * inputGain * drive;
    real_[mode] = tfdsp::FiniteNormalOrZero(
        real_[mode] + parameters_.inputPhaseCosine[mode] * force);
    imaginary_[mode] = tfdsp::FiniteNormalOrZero(
        imaginary_[mode] + parameters_.inputPhaseSine[mode] * force);
    modal += parameters_.modeOutputGain[mode] * real_[mode];
  }
  return tfdsp::FiniteNormalOrZero(
      parameters_.outputLevel *
      (parameters_.noiseMix * drive + parameters_.modalMix * modal));
}

float WireRack::StoredEnergy() const noexcept {
  float result = 0.f;
  for (std::size_t mode = 0; mode < parameters_.activeModeCount; ++mode)
    result += real_[mode] * real_[mode] + imaginary_[mode] * imaginary_[mode];
  return tfdsp::FiniteNormalOrZero(result);
}

float WireRack::AvailableDriveScale(const float baseEnergy,
                                    const float crossEnergy,
                                    const float driveEnergy) const noexcept {
  const float proposed = baseEnergy + 2.f * crossEnergy + driveEnergy;
  if (proposed <= parameters_.maximumModalEnergy)
    return 1.f;
  if (driveEnergy <= 1.e-20f)
    return 0.f;
  const float discriminant =
      crossEnergy * crossEnergy +
      driveEnergy * std::max(0.f, parameters_.maximumModalEnergy - baseEnergy);
  return std::clamp((-crossEnergy + std::sqrt(discriminant)) / driveEnergy, 0.f,
                    1.f);
}

} // namespace tfdsp::percussion
