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

} // namespace

void CrashCymbal::Prepare(const float sampleRate,
                          const CrashCymbalParameters &parameters) {
  PrepareComponents(sampleRate, parameters);
  modalField_.Prepare(sampleRate, parameters.modalField,
                      parameters.modalFieldControls, 700.f, 6500.f);
  Reset();
}

void CrashCymbal::Prepare(
    const CrashCymbalPreparedParameters &prepared) {
  PrepareComponents(prepared.sampleRate, prepared.parameters);
  modalField_.LoadPrepared(prepared.modalField);
  Reset();
}

void CrashCymbal::PrepareComponents(
    const float sampleRate, const CrashCymbalParameters &parameters) {
  if (!std::isfinite(sampleRate) || sampleRate < 1.f)
    throw std::invalid_argument("crash sample rate must be positive");
  sampleRate_ = sampleRate;
  parameters_ = parameters;
  contact_.Prepare(sampleRate);
  const float maximumDelay = std::max(512.f, .025f * sampleRate);
  dispersion_.Prepare(sampleRate, maximumDelay, parameters.dispersion);
  observation_.Prepare(sampleRate, .01f, parameters.observation);
  bloomBodyGain_ = std::clamp(
      tfdsp::FiniteNormalOrZero(parameters.fit.bloomBodyGain), 0.f, 2.f);
  routing_ = parameters.routing;
  delayConstraint_.Prepare(sampleRate, .003f, .08f);
  modalConstraint_.Prepare(sampleRate, .001f, .003f, .08f);
}

void CrashCymbal::Reset() noexcept {
  contact_.Reset();
  dispersion_.Reset();
  modalField_.Reset();
  observation_.Reset();
  delayConstraint_.Reset();
  modalConstraint_.Reset();
  SetExcitationProjection(1.f, .8f);
  modalField_.SetSecondaryExcitationProjection(parameters_.fieldBowProjection);
  bloomDriveScale_ = 1.f;
  bodyDriveScale_ = 1.f;
}

void CrashCymbal::Trigger(const CrashCymbalHit &hit) noexcept {
  const float strength = Unit(hit.strength);
  if (strength <= 0.f)
    return;
  SetExcitationProjection(hit.location, strength);
  contact_.Trigger(ContactParameters(hit));
  const float implement = Unit(hit.implement);
  const float brush = Unit(1.f - 2.f * implement);
  const float stick = Unit(2.f * implement - 1.f);
  const float mallet = 1.f - brush - stick;
  bloomDriveScale_ = .075f * brush + .55f * mallet + stick;
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
  const float bloom = dispersion_.Process(
      routing_.Get(MetallicPlateRoute::ContactToDispersion) *
          bloomDriveScale_ * bodyDrive,
      delayLoss);
  const float bodyBloom =
      routing_.Get(MetallicPlateRoute::DispersionToBody) *
      bloomBodyGain_ * bloom;
  const float body = modalField_.ProcessExcitedPair(
      routing_.Get(MetallicPlateRoute::ContactToBody) * bodyDrive,
      bodyBloom, modalLoss);
  const ObservationModel<2>::SourceFrame sources{
      routing_.Get(MetallicPlateRoute::ContactToObservation) *
          contact.directRadiation,
      routing_.Get(MetallicPlateRoute::BodyToObservation) * body};
  return {
      contact.directRadiation,
      bloom,
      body,
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
  // Reserve more of the knob's travel for the bright end, where the cymbal
  // body becomes increasingly sensitive to small bandwidth changes.
  const float brushCharacter = hardness * hardness;
  // Fine bristles sustain more overlapping contact energy. Keep perceived
  // strike strength broadly stable while retaining that longer, darker shape.
  const float brushLevel = Mix(
      .62f, 1.f, brushCharacter * brushCharacter);

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
  const float contactTiltDb = implementMix(
      Mix(-10.f, 4.f, brushCharacter),
      Mix(-10.f, -5.f, hardness),
      Mix(-3.f, 2.f, hardness));
  result.noise.tiltDb = 7.f * velocityOffset + contactTiltDb +
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactNoiseTiltDb), -18.f, 18.f);
  result.noise.tiltPivotHz = 4200.f;
  result.noise.seed = hit.seed ^ 0x4e4f4953u;
  result.microContacts.durationSeconds =
      microDurationScale *
      implementMix(Mix(8.f, 36.f, brushSpread), 1.2f, .8f) *
      Mix(.004f, .014f, edge);
  // Keep a brush's event realization stable while stiffness moves. Modulating
  // the stochastic event probability changes which contacts exist and makes a
  // continuous knob produce unrelated, non-monotonic gestures.
  result.microContacts.densityHz = implementMix(
      9000.f * Mix(.75f, 1.25f, spread),
      .25f * Mix(5500.f, 15000.f, hardness),
      1.1f * Mix(5500.f, 15000.f, hardness)) *
      std::clamp(tfdsp::FiniteNormalOrZero(fit.contactMicroDensityScale),
                 .25f, 4.f);
  result.microContacts.microDecaySeconds = implementMix(
      Mix(.0018f, .0005f, brushCharacter),
      Mix(.0009f, .00025f, hardness),
      Mix(.0009f, .00025f, hardness));
  result.microContacts.brightness = implementMix(
      Mix(.2f, .85f, brushCharacter),
      .25f * Mix(.7f, 1.f, hardness),
      Mix(.7f, 1.f, hardness));
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
      // A brush contains many more contacts than a stick. Its routing is
      // energy-normalized here while its long stochastic gesture is retained.
      implementMix(.035f * brushLevel, .24f, .32f),
      implementMix(.14f * brushLevel, .5f, .72f),
      implementMix(.07f * brushLevel, .22f, .25f),
      implementMix(.28f * brushLevel, .42f, .85f)};
  return result;
}

void CrashCymbal::SetExcitationProjection(const float value,
                                          const float strength) noexcept {
  const float location = Unit(value);
  const float velocityBrightness =
      parameters_.fit.velocityBrightnessDbPerOctave;
  fieldProjection_ = InterpolateProjection(
      location, parameters_.fieldBellProjection,
      parameters_.fieldBowProjection, parameters_.fieldEdgeProjection);
  ApplyVelocityColour(fieldProjection_, parameters_.modalField, strength,
                      velocityBrightness);
  modalField_.SetExcitationProjection(fieldProjection_);
}

} // namespace tfdsp::percussion
