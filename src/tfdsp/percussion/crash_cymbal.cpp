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

template <std::size_t Count, typename Parameters>
void ApplyVelocityColour(std::array<float, Count> &projection,
                         const Parameters &modes, const float strength,
                         const float brightnessDbPerOctave) noexcept {
  const float tilt = std::clamp(brightnessDbPerOctave, 0.f, 12.f) *
      (Unit(strength) - .8f);
  for (std::size_t mode = 0; mode < Count; ++mode) {
    const float octaves = std::max(
        std::log2(std::max(modes[mode].frequencyHz, 700.f) / 700.f), 0.f);
    const float gain = std::pow(10.f, tilt * octaves / 20.f);
    projection[mode] *= std::clamp(gain, .08f, 3.f);
  }
}

template <typename Bank>
float ProcessModalBody(Bank &bank, const bool enabled,
                       const float contact, const float bloom,
                       const ModalDampingGains damping) noexcept {
  if (!enabled)
    return 0.f;
  return contact != 0.f
      ? bank.ProcessExcitedPair(contact, bloom, damping)
      : bank.ProcessSecondaryExcited(bloom, damping);
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
  denseExtensionModes_.Prepare(
      sampleRate, parameters.denseExtensionModes, 700.f, 6500.f);
  observation_.Prepare(sampleRate, .01f, parameters.observation);
  sparseBloomGain_ = std::clamp(
      tfdsp::FiniteNormalOrZero(parameters.fit.sparseBloomGain), 0.f, 2.f);
  denseBypassGain_ = std::clamp(
      tfdsp::FiniteNormalOrZero(parameters.fit.bodyBypassGain), 0.f, 2.f);
  sparseEnabled_ = parameters.observation[1].gain > 0.f;
  denseEnabled_ = parameters.observation[2].gain > 0.f;
  denseExtensionEnabled_ = denseExtensionModes_.ActiveModeCount() > 0;
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
  denseExtensionModes_.Reset();
  observation_.Reset();
  delayConstraint_.Reset();
  modalConstraint_.Reset();
  SetExcitationProjection(1.f, .8f);
  sparseModes_.SetSecondaryExcitationProjection(
      parameters_.sparseBowProjection);
  denseModes_.SetSecondaryExcitationProjection(
      parameters_.denseBowProjection);
  denseExtensionModes_.SetSecondaryExcitationProjection(
      parameters_.denseExtensionBowProjection);
  bodyDriveScale_ = 1.f;
}

void CrashCymbal::Trigger(const CrashCymbalHit &hit) noexcept {
  const float strength = Unit(hit.strength);
  if (strength <= 0.f)
    return;
  SetExcitationProjection(hit.location, strength);
  contact_.Trigger(ContactParameters(hit));
  const float directGamma = std::clamp(parameters_.fit.strengthGamma, .25f, 4.f);
  const float bodyGamma =
      std::clamp(parameters_.fit.bodyStrengthGamma, .1f, 4.f);
  bodyDriveScale_ = std::pow(strength, bodyGamma - directGamma);
}

CrashCymbalFrame CrashCymbal::ProcessFrame() noexcept {
  const auto delayLoss = delayConstraint_.Process();
  const auto modalLoss = modalConstraint_.Process();
  const auto contact = contact_.Process();
  const float bodyDrive = bodyDriveScale_ * contact.bodyDrive;
  const float bloom = dispersion_.Process(bodyDrive, delayLoss);
  const float sparseBloom = sparseBloomGain_ * bloom;
  const float sparse = ProcessModalBody(
      sparseModes_, sparseEnabled_, bodyDrive, sparseBloom, modalLoss);
  const float denseContact = denseBypassGain_ * bodyDrive;
  const float modalDense = ProcessModalBody(
      denseModes_, denseEnabled_, denseContact, bloom, modalLoss);
  const float modalDenseExtension = ProcessModalBody(
      denseExtensionModes_, denseEnabled_ && denseExtensionEnabled_,
      denseContact, bloom, modalLoss);
  const float turbulent = turbulenceEnabled_
      ? turbulence_.Process(bloom, modalLoss)
      : 0.f;
  const float dense = tfdsp::FiniteNormalOrZero(
      modalDense + modalDenseExtension);
  const ObservationModel<4>::SourceFrame sources{
      contact.directRadiation, sparse, dense, turbulent};
  return {
      contact.directRadiation,
      bloom,
      sparse,
      dense,
      turbulent,
      tfdsp::FiniteNormalOrZero(dense + turbulent),
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
  const float implement = Unit(hit.implement);
  const float spread = Unit(hit.contactSpread);
  const float brush = Unit(1.f - 2.f * implement);
  const float stick = Unit(2.f * implement - 1.f);
  const float mallet = 1.f - brush - stick;
  const auto implementMix = [&](const float brushValue,
                                const float malletValue,
                                const float stickValue) noexcept {
    return brush * brushValue + mallet * malletValue + stick * stickValue;
  };
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
  const float implementWidth = implementMix(1.2f, 2.2f, .75f);
  // Faster collisions are shorter and transfer disproportionately more energy
  // to fine contact structure. The .8 anchor preserves the calibrated nominal
  // strike while still giving velocity a timbral, rather than merely scalar,
  // response.
  const float velocityOffset = strength - .8f;
  const float velocityDuration = std::exp2(-.55f * velocityOffset);
  const float velocityNoise = std::exp2(1.2f * velocityOffset);
  const float velocityMicro = std::exp2(1.5f * velocityOffset);
  const float velocityChirp = std::exp2(.45f * velocityOffset);
  // Even a brush tap distributes contact across a bundle of bristles. Spread
  // extends that baseline gesture; it never collapses the brush into a stick-
  // like impulse.
  const float brushSpread = .3f + .7f * spread;

  ContactExciterParameters result;
  result.pulseDurationSeconds =
      velocityDuration * implementWidth * durationScale *
      Mix(.0035f, .00055f, hardness);
  result.pulseAmplitude = amplitude * Mix(.7f, 1.15f, hardness) *
      implementMix(0.f, .85f, 1.08f) *
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactPulseGain), 0.f, 4.f);
  result.chirp.durationSeconds =
      velocityDuration * implementWidth * durationScale *
      Mix(.006f, .0015f, hardness);
  const float chirpFrequencyScale = std::clamp(
      tfdsp::FiniteNormalOrZero(fit.contactChirpFrequencyScale), .25f, 4.f);
  result.chirp.startFrequencyHz = Mix(2800.f, 9200.f, hardness) *
      velocityChirp * (1.f + .35f * bell) * chirpFrequencyScale;
  result.chirp.endFrequencyHz = Mix(1300.f, 4100.f, hardness) *
      velocityChirp * (1.f + .25f * bell) * chirpFrequencyScale;
  result.chirp.amplitude = amplitude * hardness *
      Mix(.75f, .35f, edge) * (1.f + .5f * bell) *
      implementMix(0.f, .08f, 1.f) *
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactChirpGain), 0.f, 4.f);
  result.chirp.decayNepers = 2.4f;
  result.noise.attackSeconds = noiseDurationScale *
      implementMix(Mix(.003f, .009f, brushSpread), .00035f, .00015f);
  result.noise.holdSeconds =
      noiseDurationScale *
      implementMix(Mix(8.f, 35.f, brushSpread), 1.5f, .85f) *
      Mix(.0004f, .0015f, edge);
  result.noise.decaySeconds =
      noiseDurationScale *
      implementMix(Mix(10.f, 28.f, brushSpread), 1.5f, .85f) *
      Mix(.006f, .018f, edge);
  result.noise.amplitude = velocityNoise * amplitude * Mix(.25f, .7f, edge) *
      implementMix(1.1f, .7f, .88f) *
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactNoiseGain), 0.f, 4.f);
  result.noise.tiltDb = 7.f * velocityOffset + Mix(-4.f, 1.f, hardness) +
      implementMix(-1.f, -6.f, 1.f) +
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactNoiseTiltDb), -18.f, 18.f);
  result.noise.tiltPivotHz = 4200.f;
  result.noise.seed = hit.seed ^ 0x4e4f4953u;
  result.microContacts.durationSeconds =
      microDurationScale *
      implementMix(Mix(8.f, 36.f, brushSpread), 1.2f, .8f) *
      Mix(.004f, .014f, edge);
  result.microContacts.densityHz = implementMix(
      Mix(1200.f, 6500.f, hardness) * Mix(.75f, 1.25f, spread),
      .25f * Mix(5500.f, 15000.f, hardness),
      1.1f * Mix(5500.f, 15000.f, hardness)) *
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactMicroDensityScale),
                 .25f, 4.f);
  result.microContacts.microDecaySeconds = implementMix(
      Mix(.0014f, .00055f, hardness),
      Mix(.0009f, .00025f, hardness),
      Mix(.0009f, .00025f, hardness));
  result.microContacts.brightness = implementMix(.65f, .25f, 1.f) *
      Mix(.7f, 1.f, hardness);
  result.microContacts.amplitude = velocityMicro * amplitude *
      Mix(.18f, .5f, edge) *
      implementMix(1.4f, .3f, 1.f) *
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactMicroGain), 0.f, 4.f);
  result.microContacts.seed = hit.seed ^ 0x4d494352u;
  result.routing = {
      implementMix(0.f, .02f, .02f),
      implementMix(0.f, .95f, .95f),
      implementMix(0.f, .2f, .8f),
      implementMix(0.f, .12f, .28f),
      implementMix(.5f, .24f, .32f),
      implementMix(.35f, .5f, .72f),
      implementMix(1.f, .22f, .25f),
      implementMix(.7f, .42f, .85f)};
  return result;
}

void CrashCymbal::SetExcitationProjection(const float value,
                                          const float strength) noexcept {
  const float location = Unit(value);
  sparseProjection_ = InterpolateProjection(
      location, parameters_.sparseBellProjection,
      parameters_.sparseBowProjection, parameters_.sparseEdgeProjection);
  denseProjection_ = InterpolateProjection(
      location, parameters_.denseBellProjection,
      parameters_.denseBowProjection, parameters_.denseEdgeProjection);
  denseExtensionProjection_ = InterpolateProjection(
      location, parameters_.denseExtensionBellProjection,
      parameters_.denseExtensionBowProjection,
      parameters_.denseExtensionEdgeProjection);
  const float velocityBrightness =
      parameters_.fit.velocityBrightnessDbPerOctave;
  ApplyVelocityColour(sparseProjection_, parameters_.sparseModes, strength,
                      velocityBrightness);
  ApplyVelocityColour(denseProjection_, parameters_.denseModes, strength,
                      velocityBrightness);
  ApplyVelocityColour(denseExtensionProjection_,
                      parameters_.denseExtensionModes, strength,
                      velocityBrightness);
  sparseModes_.SetExcitationProjection(sparseProjection_);
  denseModes_.SetExcitationProjection(denseProjection_);
  denseExtensionModes_.SetExcitationProjection(denseExtensionProjection_);
}

} // namespace tfdsp::percussion
