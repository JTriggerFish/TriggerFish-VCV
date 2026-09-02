#include "crash_cymbal.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace tfdsp::percussion {
namespace {

float Unit(const float value) noexcept {
  return std::clamp(tfdsp::FiniteNormalOrZero(value), 0.f, 1.f);
}

float Mix(const float first, const float second, const float amount) noexcept {
  return first + Unit(amount) * (second - first);
}

template <std::size_t Count>
std::array<float, Count> InterpolateProjection(
    const float location, const std::array<float, Count> &bell,
    const std::array<float, Count> &bow,
    const std::array<float, Count> &edge) noexcept {
  const bool outer = location >= .5f;
  const float amount = outer ? 2.f * location - 1.f : 2.f * location;
  const auto &first = outer ? bow : bell;
  const auto &second = outer ? edge : bow;
  std::array<float, Count> result{};
  for (std::size_t mode = 0; mode < Count; ++mode)
    result[mode] = Mix(first[mode], second[mode], amount);
  return result;
}

} // namespace

void CrashCymbal::Prepare(const float sampleRate,
                          const CrashCymbalParameters &parameters) {
  if (!std::isfinite(sampleRate) || sampleRate < 1.f)
    throw std::invalid_argument("crash sample rate must be positive");
  sampleRate_ = sampleRate;
  parameters_ = parameters;
  contact_.Prepare(sampleRate);
  const float maximumDelay = std::max(512.f, .025f * sampleRate);
  dispersion_.Prepare(sampleRate, maximumDelay, parameters.dispersion);
  turbulence_.Prepare(sampleRate, parameters.turbulence);
  sparseModes_.Prepare(sampleRate, parameters.sparseModes, 700.f, 6500.f);
  denseModes_.Prepare(sampleRate, parameters.denseModes, 700.f, 6500.f);
  observation_.Prepare(sampleRate, .01f, parameters.observation);
  sparseBloomGain_ = std::clamp(
      tfdsp::FiniteNormalOrZero(parameters.fit.sparseBloomGain), 0.f, 2.f);
  denseBypassGain_ = std::clamp(
      tfdsp::FiniteNormalOrZero(parameters.fit.bodyBypassGain), 0.f, 2.f);
  sparseEnabled_ = parameters.observation[1].gain > 0.f;
  denseEnabled_ = parameters.observation[2].gain > 0.f;
  turbulenceEnabled_ = std::any_of(
      parameters.turbulence.gain.begin(), parameters.turbulence.gain.end(),
      [](const float gain) { return gain > 0.f; });
  delayConstraint_.Prepare(sampleRate, .003f, .08f);
  modalConstraint_.Prepare(sampleRate, .001f, .003f, .08f);
  Reset();
}

void CrashCymbal::Reset() noexcept {
  contact_.Reset();
  dispersion_.Reset();
  turbulence_.Reset();
  sparseModes_.Reset();
  denseModes_.Reset();
  observation_.Reset();
  delayConstraint_.Reset();
  modalConstraint_.Reset();
  SetLocation(1.f);
  bodyDriveScale_ = 1.f;
  denseDriveScale_ = 1.f;
  denseVelocityLoss_ = 1.f;
}

void CrashCymbal::Trigger(const CrashCymbalHit &hit) noexcept {
  const float strength = Unit(hit.strength);
  SetLocation(hit.location);
  contact_.Trigger(ContactParameters(hit));
  turbulence_.Seed(hit.seed ^ 0x54555242u);
  const float directGamma = std::clamp(parameters_.fit.strengthGamma, .25f, 4.f);
  const float bodyGamma =
      std::clamp(parameters_.fit.bodyStrengthGamma, .1f, 4.f);
  const float denseGamma =
      std::clamp(parameters_.fit.denseStrengthGamma, .1f, 4.f);
  bodyDriveScale_ = strength > 0.f
      ? std::pow(strength, bodyGamma - directGamma)
      : 0.f;
  denseDriveScale_ = strength > 0.f
      ? std::pow(strength, denseGamma - bodyGamma)
      : 0.f;
  const float velocityLoss = std::clamp(
      tfdsp::FiniteNormalOrZero(
          parameters_.fit.denseVelocityLossNepersPerSecond), 0.f, 8.f);
  denseVelocityLoss_ = std::exp(-velocityLoss * strength / sampleRate_);
  SetBloomDrive(strength);
}

CrashCymbalFrame CrashCymbal::ProcessFrame() noexcept {
  const auto delayLoss = delayConstraint_.Process();
  const auto modalLoss = modalConstraint_.Process();
  const auto contact = contact_.Process();
  const float bodyDrive = bodyDriveScale_ * contact.bodyDrive;
  const float bloom = dispersion_.Process(bodyDrive, delayLoss);
  const float sparseDrive = bodyDrive + sparseBloomGain_ * bloom;
  const float denseDrive = denseDriveScale_ * (
      bloom + denseBypassGain_ * bodyDrive);
  auto denseLoss = modalLoss;
  denseLoss.broadband *= denseVelocityLoss_;
  const float sparse = sparseEnabled_
      ? sparseModes_.ProcessExcited(sparseDrive, modalLoss)
      : 0.f;
  const float modalDense = denseEnabled_
      ? denseModes_.ProcessExcited(denseDrive, denseLoss)
      : 0.f;
  const float turbulent = denseEnabled_ && turbulenceEnabled_
      ? turbulence_.Process(bloom, denseLoss)
      : 0.f;
  const float dense = tfdsp::FiniteNormalOrZero(modalDense + turbulent);
  const ObservationModel<3>::SourceFrame sources{
      contact.directRadiation, sparse, dense};
  return {
      contact.directRadiation,
      bloom,
      sparse,
      dense,
      tfdsp::FiniteNormalOrZero(
          parameters_.fit.outputGain * observation_.Process(sources))};
}

float CrashCymbal::Process() noexcept {
  return ProcessFrame().output;
}

void CrashCymbal::SetMute(const float amount) noexcept {
  const float mute = Unit(amount);
  const float broadband = std::exp(-3.2f * mute);
  const PassiveConstraintGains target{
      broadband, std::exp(-2.2f * mute),
      std::exp(-4.2f * mute), std::exp(-7.5f * mute)};
  delayConstraint_.SetTarget(target);
  modalConstraint_.SetTarget(target);
}

float CrashCymbal::MinimumBodyDelaySamples() const noexcept {
  return dispersion_.MinimumPropagationSamples();
}

ContactExciterParameters CrashCymbal::ContactParameters(
    const CrashCymbalHit &hit) const noexcept {
  const float strength = Unit(hit.strength);
  const float location = Unit(hit.location);
  const float hardness = Unit(hit.hardness);
  const float amplitude = std::pow(
      strength, std::clamp(parameters_.fit.strengthGamma, .25f, 4.f));
  const float bell = Unit(1.f - 2.f * location);
  const float edge = Unit(2.f * location - 1.f);
  const auto &fit = parameters_.fit;
  const float durationScale = std::clamp(
      tfdsp::FiniteNormalOrZero(fit.contactDurationScale), .25f, 4.f);
  const float noiseDurationScale = std::clamp(
      tfdsp::FiniteNormalOrZero(fit.contactNoiseDurationScale), .25f, 4.f);
  const float microDurationScale = std::clamp(
      tfdsp::FiniteNormalOrZero(fit.contactMicroDurationScale), .25f, 4.f);

  ContactExciterParameters result;
  result.pulseDurationSeconds =
      durationScale * Mix(.0035f, .00055f, hardness);
  result.pulseAmplitude = amplitude * Mix(.7f, 1.15f, hardness) *
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactPulseGain), 0.f, 4.f);
  result.chirp.durationSeconds =
      durationScale * Mix(.006f, .0015f, hardness);
  const float chirpFrequencyScale = std::clamp(
      tfdsp::FiniteNormalOrZero(fit.contactChirpFrequencyScale), .25f, 4.f);
  result.chirp.startFrequencyHz = Mix(2800.f, 9200.f, hardness) *
      (1.f + .35f * bell) * chirpFrequencyScale;
  result.chirp.endFrequencyHz = Mix(1300.f, 4100.f, hardness) *
      (1.f + .25f * bell) * chirpFrequencyScale;
  result.chirp.amplitude = amplitude * hardness *
      Mix(.75f, .35f, edge) * (1.f + .5f * bell) *
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactChirpGain), 0.f, 4.f);
  result.chirp.decayNepers = 2.4f;
  result.noise.attackSeconds = .00015f;
  result.noise.holdSeconds =
      noiseDurationScale * Mix(.0004f, .0015f, edge);
  result.noise.decaySeconds =
      noiseDurationScale * Mix(.006f, .018f, edge);
  result.noise.amplitude = amplitude * Mix(.25f, .7f, edge) *
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactNoiseGain), 0.f, 4.f);
  result.noise.tiltDb = Mix(-4.f, 1.f, hardness) +
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactNoiseTiltDb), -18.f, 18.f);
  result.noise.tiltPivotHz = 4200.f;
  result.noise.seed = hit.seed ^ 0x4e4f4953u;
  result.microContacts.durationSeconds =
      microDurationScale * Mix(.004f, .014f, edge);
  result.microContacts.densityHz = Mix(5500.f, 15000.f, hardness) *
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactMicroDensityScale),
                 .25f, 4.f);
  result.microContacts.microDecaySeconds = Mix(.0009f, .00025f, hardness);
  result.microContacts.brightness = Mix(.45f, .95f, hardness);
  result.microContacts.amplitude = amplitude * Mix(.18f, .5f, edge) *
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactMicroGain), 0.f, 4.f);
  result.microContacts.seed = hit.seed ^ 0x4d494352u;
  result.routing = {.02f, .95f, .8f, .28f, .32f, .72f, .25f, .85f};
  return result;
}

void CrashCymbal::SetLocation(const float value) noexcept {
  const float location = Unit(value);
  sparseProjection_ = InterpolateProjection(
      location, parameters_.sparseBellProjection,
      parameters_.sparseBowProjection, parameters_.sparseEdgeProjection);
  denseProjection_ = InterpolateProjection(
      location, parameters_.denseBellProjection,
      parameters_.denseBowProjection, parameters_.denseEdgeProjection);
  sparseModes_.SetExcitationProjection(sparseProjection_);
  denseModes_.SetExcitationProjection(denseProjection_);
}

void CrashCymbal::SetBloomDrive(const float strength) noexcept {
  auto selfPhase = parameters_.dispersion.selfPhase;
  selfPhase.drive *= .18f + .82f * std::sqrt(Unit(strength));
  selfPhase.maximumExcursionSamples *= .35f + .65f * Unit(strength);
  dispersion_.SetSelfPhaseParameters(selfPhase);
}

} // namespace tfdsp::percussion
