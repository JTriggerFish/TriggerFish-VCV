#include "crash_cymbal.hpp"

#include <stdexcept>

namespace tfdsp::percussion {
namespace {

float Unit(const float value) noexcept {
  return std::clamp(tfdsp::FiniteNormalOrZero(value), 0.f, 1.f);
}

float Mix(const float first, const float second, const float amount) noexcept {
  return first + Unit(amount) * (second - first);
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
  resonators_.Prepare(sampleRate, maximumDelay, parameters.resonators,
                      700.f, 6500.f);
  resonators_.SetCoupling(parameters.fit.resonatorCoupling);
  ResonatorSubmix<CrashResonatorCount,
                  CrashResonatorBusCount>::Weights weights{};
  for (std::size_t bus = 0; bus < CrashResonatorBusCount; ++bus) {
    for (std::size_t line = 0; line < CrashResonatorCount; ++line)
      weights[bus][line] = line % CrashResonatorBusCount == bus ? 1.f : 0.f;
    shifters_[bus].Prepare(sampleRate);
    shifters_[bus].SetShiftHz(parameters.resonatorShiftHz[bus]);
  }
  submix_.SetWeights(weights);
  observation_.Prepare(sampleRate, .01f, parameters.observation);
  constraint_.Prepare(sampleRate, .003f, .08f);
  Reset();
}

void CrashCymbal::Reset() noexcept {
  contact_.Reset();
  dispersion_.Reset();
  resonators_.Reset();
  for (auto &shifter : shifters_)
    shifter.Reset();
  observation_.Reset();
  constraint_.Reset();
  projection_ = parameters_.edgeProjection;
}

void CrashCymbal::Trigger(const CrashCymbalHit &hit) noexcept {
  const float strength = Unit(hit.strength);
  projection_ = LocationProjection(hit.location);
  contact_.Trigger(ContactParameters(hit));
  SetBloomDrive(strength);
}

float CrashCymbal::Process() noexcept {
  const auto loss = constraint_.Process();
  const auto contact = contact_.Process();
  const float bloom = dispersion_.Process(contact.bodyDrive, loss);
  const float networkDrive = bloom +
      parameters_.fit.bodyBypassGain * contact.bodyDrive;
  const auto body = resonators_.ProcessProjected(networkDrive, projection_, loss);
  const auto buses = submix_.Process(body.lines);
  float shiftedBody = 0.f;
  for (std::size_t bus = 0; bus < CrashResonatorBusCount; ++bus)
    shiftedBody += shifters_[bus].Process(buses[bus]);
  const ObservationModel<2>::SourceFrame sources{
      contact.directRadiation, shiftedBody};
  return tfdsp::FiniteNormalOrZero(
      parameters_.fit.outputGain * observation_.Process(sources));
}

void CrashCymbal::SetMute(const float amount) noexcept {
  const float mute = Unit(amount);
  const float broadband = std::exp(-3.2f * mute);
  constraint_.SetTarget({broadband, std::exp(-2.2f * mute),
                         std::exp(-4.2f * mute), std::exp(-7.5f * mute)});
}

float CrashCymbal::MinimumBodyDelaySamples() const noexcept {
  float shortest = parameters_.resonators.front().delaySamples;
  for (const auto &line : parameters_.resonators)
    shortest = std::min(shortest, line.delaySamples);
  return dispersion_.MinimumPropagationSamples() + shortest;
}

ContactExciterParameters CrashCymbal::ContactParameters(
    const CrashCymbalHit &hit) const noexcept {
  const float strength = Unit(hit.strength);
  const float location = Unit(hit.location);
  const float hardness = Unit(hit.hardness);
  const float amplitude = std::pow(strength,
      std::clamp(parameters_.fit.strengthGamma, .25f, 4.f));
  const float bell = Unit(1.f - 2.f * location);
  const float edge = Unit(2.f * location - 1.f);

  ContactExciterParameters result;
  result.pulseDurationSeconds = Mix(.0035f, .00055f, hardness);
  result.pulseAmplitude = amplitude * Mix(.7f, 1.15f, hardness);
  result.chirp.durationSeconds = Mix(.006f, .0015f, hardness);
  result.chirp.startFrequencyHz = Mix(2800.f, 9200.f, hardness) *
      (1.f + .35f * bell);
  result.chirp.endFrequencyHz = Mix(1300.f, 4100.f, hardness) *
      (1.f + .25f * bell);
  result.chirp.amplitude = amplitude * hardness *
      Mix(.75f, .35f, edge) * (1.f + .5f * bell);
  result.chirp.decayNepers = 2.4f;
  result.noise.attackSeconds = .00015f;
  result.noise.holdSeconds = Mix(.0004f, .0015f, edge);
  result.noise.decaySeconds = Mix(.006f, .018f, edge);
  result.noise.amplitude = amplitude * Mix(.25f, .7f, edge);
  result.noise.tiltDb = Mix(-4.f, 1.f, hardness);
  result.noise.tiltPivotHz = 4200.f;
  result.noise.seed = hit.seed ^ 0x4e4f4953u;
  result.microContacts.durationSeconds = Mix(.004f, .014f, edge);
  result.microContacts.densityHz = Mix(5500.f, 15000.f, hardness);
  result.microContacts.microDecaySeconds = Mix(.0009f, .00025f, hardness);
  result.microContacts.brightness = Mix(.45f, .95f, hardness);
  result.microContacts.amplitude = amplitude * Mix(.18f, .5f, edge);
  result.microContacts.seed = hit.seed ^ 0x4d494352u;
  result.routing = {.02f, .95f, .8f, .28f, .32f, .72f, .25f, .85f};
  return result;
}

CrashCymbal::Projection CrashCymbal::LocationProjection(
    const float value) const noexcept {
  const float location = Unit(value);
  const bool outer = location >= .5f;
  const float amount = outer ? 2.f * location - 1.f : 2.f * location;
  const Projection &first = outer ? parameters_.bowProjection
                                  : parameters_.bellProjection;
  const Projection &second = outer ? parameters_.edgeProjection
                                   : parameters_.bowProjection;
  Projection result{};
  for (std::size_t line = 0; line < result.size(); ++line)
    result[line] = Mix(first[line], second[line], amount);
  return result;
}

void CrashCymbal::SetBloomDrive(const float strength) noexcept {
  auto selfPhase = parameters_.dispersion.selfPhase;
  selfPhase.drive *= .18f + .82f * std::sqrt(Unit(strength));
  selfPhase.maximumExcursionSamples *= .35f + .65f * Unit(strength);
  dispersion_.SetSelfPhaseParameters(selfPhase);
}

} // namespace tfdsp::percussion
