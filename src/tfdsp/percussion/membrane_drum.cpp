#include "membrane_drum.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {
namespace {

float Unit(const float value, const float fallback = 0.f) noexcept {
  return std::clamp(std::isfinite(value) ? value : fallback, 0.f, 1.f);
}

ContactExciterParameters ContactForHit(
    ContactExciterParameters result, const MembraneDrumHit &hit) noexcept {
  const float hardness = Unit(hit.hardness, .5f);
  const float implement = Unit(hit.implement, 1.f);
  const float spread = Unit(hit.contactSpread, .25f);
  const float brush = std::clamp(1.f - 2.f * implement, 0.f, 1.f);
  const float mallet = std::clamp(1.f - 2.f * std::abs(implement - .5f), 0.f, 1.f);
  const float stick = std::clamp(2.f * implement - 1.f, 0.f, 1.f);

  result.pulseDurationSeconds *= .45f * stick + 1.7f * mallet +
      (2.f + 3.f * spread) * brush;
  result.pulseAmplitude *= .95f * stick + .65f * mallet + .08f * brush;
  result.chirp.durationSeconds *= .7f + .8f * (1.f - hardness);
  result.chirp.amplitude *= stick * (.25f + 1.2f * hardness) +
      .18f * mallet + .015f * brush;
  result.chirp.startFrequencyHz *= .65f + .75f * hardness;
  result.chirp.endFrequencyHz *= .72f + .5f * hardness;
  result.noise.decaySeconds *= .65f + .9f * spread + 1.7f * brush;
  result.noise.amplitude *= .3f * stick + .45f * mallet +
      (1.1f + .7f * spread) * brush;
  result.noise.tiltDb += 7.f * (hardness - .5f) - 5.f * brush;
  result.microContacts.durationSeconds *= .7f + 3.2f * brush + spread;
  result.microContacts.densityHz *= .7f + 1.2f * brush + .6f * spread;
  result.microContacts.microDecaySeconds *= 1.5f - .8f * hardness;
  result.microContacts.brightness = Unit(
      result.microContacts.brightness + .35f * (hardness - .5f));
  result.microContacts.amplitude *= .18f * stick + .42f * mallet +
      (1.25f + .5f * spread) * brush;
  result.noise.seed = hit.seed ^ 0x9e3779b9u;
  result.microContacts.seed = hit.seed ^ 0x85ebca6bu;
  return result;
}

} // namespace

void MembraneDrum::EventVoice::Prepare(const float sampleRate) {
  contact.Prepare(sampleRate);
  fm.Prepare(sampleRate);
  stealStep = 1.f / std::max(1.f, .001f * sampleRate);
  activityRelease = std::exp(-1.f / std::max(1.f, .015f * sampleRate));
}

void MembraneDrum::EventVoice::Reset() noexcept {
  contact.Reset();
  fm.Reset();
  location = .5f;
  amplitude = contactLevel = fmLevel = activity = 0.f;
  stealDirect = stealBody = stealFm = stealGain = 0.f;
  last = {};
  generation = 0;
}

void MembraneDrum::EventVoice::Trigger(
    const MembraneDrumParameters &parameters, const MembraneDrumHit &hit,
    const std::uint64_t eventGeneration) noexcept {
  if (Active()) {
    stealDirect = last.contactDirect;
    stealBody = last.contactBody;
    stealFm = last.fm;
    stealGain = 1.f;
  } else {
    stealDirect = stealBody = stealFm = stealGain = 0.f;
  }
  const float strength = Unit(hit.strength, .8f);
  auto contactParameters = ContactForHit(parameters.contact, hit);
  auto fmParameters = parameters.fm;
  fmParameters.seed = hit.seed;
  fmParameters.frequencyDeviationHz.initialValue *= .35f + .9f * strength;
  contact.Trigger(contactParameters);
  fm.Trigger(fmParameters);
  location = Unit(hit.location, .5f);
  amplitude = std::pow(strength, .72f);
  contactLevel = parameters.contactLevel;
  fmLevel = parameters.fmLevel * (.25f + .75f * strength);
  generation = eventGeneration;
}

MembraneDrum::EventVoice::Sample MembraneDrum::EventVoice::Process() noexcept {
  const auto contactSample = contact.Process();
  Sample result{contactLevel * amplitude * contactSample.directRadiation,
                contactLevel * amplitude * contactSample.bodyDrive,
                fmLevel * amplitude * fm.Process()};
  if (stealGain > 0.f) {
    result.contactDirect += stealGain * stealDirect;
    result.contactBody += stealGain * stealBody;
    result.fm += stealGain * stealFm;
    stealGain = std::max(0.f, stealGain - stealStep);
  }
  last = result;
  const float peak = std::max({std::abs(result.contactDirect),
                               std::abs(result.contactBody),
                               std::abs(result.fm)});
  activity = std::max(peak, activityRelease * activity);
  return result;
}

bool MembraneDrum::EventVoice::Active() const noexcept {
  return contact.Active() || fm.Active() || stealGain > 0.f;
}

void MembraneDrum::Prepare(const float sampleRate,
                           const MembraneDrumParameters &parameters) {
  PrepareComponents(sampleRate, parameters);
  membrane_.Prepare(sampleRate, parameters.membrane,
                    parameters.maximumModalEnergy);
  Reset();
}

void MembraneDrum::Prepare(const MembraneDrumPreparedParameters &prepared) {
  PrepareComponents(prepared.sampleRate, prepared.parameters);
  membrane_.LoadPrepared(prepared.membrane);
  Reset();
}

void MembraneDrum::PrepareComponents(
    const float sampleRate, const MembraneDrumParameters &parameters) {
  parameters_ = parameters;
  for (auto &voice : voices_) voice.Prepare(sampleRate);
  strikeEnergy_.Prepare(sampleRate, parameters.strikeEnergy);
  observation_.Prepare(sampleRate, .01f, parameters.observation);
  equalizer_.Prepare(sampleRate, parameters.equalizer);
  directMixer_.SetGains({
      parameters.routing.Get(MembraneDrumRoute::ContactToDirect),
      parameters.routing.Get(MembraneDrumRoute::FmToDirect)});
  bodyMixer_.SetGains({
      parameters.routing.Get(MembraneDrumRoute::ContactToBody),
      parameters.routing.Get(MembraneDrumRoute::FmToBody)});
}

void MembraneDrum::Reset() noexcept {
  for (auto &voice : voices_) voice.Reset();
  membrane_.Reset();
  strikeEnergy_.Reset();
  observation_.Reset();
  equalizer_.Reset();
  generation_ = 0;
}

MembraneDrum::EventVoice &MembraneDrum::NextVoice() noexcept {
  const auto inactive = std::find_if(voices_.begin(), voices_.end(),
      [](const auto &voice) { return !voice.Active(); });
  if (inactive != voices_.end()) return *inactive;
  return *std::min_element(voices_.begin(), voices_.end(),
      [](const auto &left, const auto &right) {
        if (left.Activity() != right.Activity())
          return left.Activity() < right.Activity();
        return left.generation < right.generation;
      });
}

void MembraneDrum::Trigger(const MembraneDrumHit &hit) noexcept {
  const float strength = Unit(hit.strength, .8f);
  if (strength <= 0.f) return;
  strikeEnergy_.InjectStrike(strength * strength);
  NextVoice().Trigger(parameters_, hit, ++generation_);
}

MembraneDrumFrame MembraneDrum::ProcessFrame() noexcept {
  MembraneResonator<MembraneModeCount>::Drive membraneDrive{};
  float direct = 0.f;
  for (auto &voice : voices_) {
    const auto event = voice.Process();
    direct += directMixer_.Process({event.contactDirect, event.fm});
    const float bodyForce = bodyMixer_.Process({event.contactBody, event.fm});
    const auto projected = membrane_.Project(bodyForce, voice.location);
    for (std::size_t mode = 0; mode < membraneDrive.size(); ++mode)
      membraneDrive[mode] += projected[mode];
  }
  strikeEnergy_.Process();
  const float body = parameters_.routing.Get(
      MembraneDrumRoute::BodyToObservation) *
      membrane_.Process(membraneDrive, strikeEnergy_.TensionScale());
  const float observed = observation_.Process({direct, body});
  const float output = tfdsp::FiniteNormalOrZero(
      parameters_.outputGain * equalizer_.Process(observed));
  return {direct, body, output};
}

float MembraneDrum::Process() noexcept { return ProcessFrame().output; }

} // namespace tfdsp::percussion
