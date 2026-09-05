#include "membrane_drum_parameters.hpp"

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace tfdsp::percussion {
namespace {

using Curve = TrajectoryCurve;

float Safe(const float value, const float fallback,
           const float low, const float high) noexcept {
  return std::clamp(std::isfinite(value) ? value : fallback, low, high);
}

CorrelatedFmBurstParameters MakeFm(const MembraneDrumControls &controls,
                                   const float fundamental) noexcept {
  CorrelatedFmBurstParameters result;
  const float decay = Safe(controls.fmDecaySeconds, .07f, .003f, 1.f);
  const float drop = Safe(controls.pitchDropOctaves, .28f, 0.f, 3.f);
  result.amplitude.initialValue = 0.f;
  result.amplitude.segments[0] = {1.f, .0004f, Curve::Linear};
  result.amplitude.segments[1] = {.0001f, decay, Curve::Geometric};
  result.amplitude.segments[2] = {0.f, .001f, Curve::Linear};
  result.amplitude.segmentCount = 3;
  result.carrierFrequencyHz.initialValue = fundamental * std::exp2(drop);
  result.carrierFrequencyHz.segments[0] = {
      fundamental, std::min(.5f, .7f * decay), Curve::Geometric};
  result.carrierFrequencyHz.segmentCount = 1;
  result.frequencyDeviationHz.initialValue = Safe(
      controls.fmDepthHz, 260.f, 0.f, 8000.f);
  result.frequencyDeviationHz.segments[0] = {
      0.f, std::min(.5f, .8f * decay), Curve::Geometric};
  result.frequencyDeviationHz.segmentCount = 1;
  result.irregularCutoffHz = 2600.f;
  result.periodicModulatorHz = 1.58f * fundamental;
  result.periodicMix = .18f;
  result.pitchPerturbationCents = 2.f;
  result.deviationPerturbation = .03f;
  return result;
}

} // namespace

bool MembraneDrumRouting::Enabled(const MembraneDrumRoute route) const noexcept {
  return enabled[static_cast<std::size_t>(route)];
}

void MembraneDrumRouting::SetEnabled(const std::size_t index,
                                     const bool value) noexcept {
  if (index < enabled.size()) enabled[index] = value;
}

MembraneDrumParameters DefaultMembraneDrumParameters(
    const MembraneDrumControls &source) noexcept {
  MembraneDrumParameters result;
  const float fundamental = Safe(source.fundamentalHz, 105.f, 25.f, 500.f);
  const float decay = Safe(source.decaySeconds, 1.15f, .03f, 8.f);
  const float decayTilt = Safe(source.decayTilt, .55f, -1.f, 1.f);
  const float inharmonicity = Safe(source.inharmonicity, .35f, 0.f, 1.f);
  const float brightness = Safe(source.bodyBrightness, .55f, 0.f, 1.f);
  // First axisymmetric and split circular-membrane families, then a smooth
  // constructive interpolation toward a less strictly membrane-like body.
  constexpr std::array<float, MembraneModeCount> membraneRatios{
      1.f, 1.593f, 2.136f, 2.296f, 2.653f, 2.918f, 3.156f, 3.501f,
      3.600f, 3.652f, 4.060f, 4.154f, 4.291f, 4.601f, 4.831f, 5.136f};
  for (std::size_t mode = 0; mode < MembraneModeCount; ++mode) {
    const float harmonic = static_cast<float>(mode + 1);
    const float ratio = harmonic + inharmonicity *
        (membraneRatios[mode] - harmonic);
    const float normalized = static_cast<float>(mode) /
        static_cast<float>(MembraneModeCount - 1);
    auto &target = result.membrane[mode];
    target.frequencyHz = fundamental * ratio;
    target.decaySeconds = decay * std::exp2(
        -2.f * decayTilt * normalized);
    target.inputGain = (.45f + 1.1f * brightness * normalized) /
        std::sqrt(1.f + static_cast<float>(mode));
    target.outputGain = (1.f - .55f * normalized) /
        std::sqrt(1.f + static_cast<float>(mode));
    const float alternating = (mode & 1u) ? -1.f : 1.f;
    target.centerProjection = (mode == 0 ? 1.f : .65f * alternating) *
        (1.f - .45f * normalized);
    target.edgeProjection = .18f + .82f * normalized;
  }

  result.strikeEnergy.releaseSeconds = Safe(
      source.tensionDecaySeconds, .13f, .005f, 4.f);
  result.strikeEnergy.capacity = 2.5f;
  result.strikeEnergy.tensionOctaves = Safe(
      source.tensionOctaves, .11f, -.5f, .8f);
  result.fm = MakeFm(source, fundamental);

  result.contact.pulseDurationSeconds = Safe(
      source.contactDurationSeconds, .004f, .0002f, .08f);
  result.contact.pulseAmplitude = 1.f;
  result.contact.chirp.durationSeconds = .004f;
  result.contact.chirp.startFrequencyHz = 5200.f;
  result.contact.chirp.endFrequencyHz = 1700.f;
  result.contact.chirp.amplitude = .5f;
  result.contact.noise.attackSeconds = .0001f;
  result.contact.noise.holdSeconds = .0004f;
  result.contact.noise.decaySeconds = .012f;
  result.contact.noise.amplitude = .45f;
  result.contact.noise.tiltDb = 12.f *
      (Safe(source.contactBrightness, .58f, 0.f, 1.f) - .5f);
  result.contact.microContacts.durationSeconds = .018f;
  result.contact.microContacts.densityHz = 9000.f;
  result.contact.microContacts.microDecaySeconds = .00045f;
  result.contact.microContacts.brightness = Safe(
      source.contactBrightness, .58f, 0.f, 1.f);
  result.contact.microContacts.amplitude = .35f;
  result.contactDirectLevel = Safe(
      source.contactDirectLevel, .245f, 0.f, 4.f);
  result.contactBodyLevel = Safe(source.contactBodyLevel, .7f, 0.f, 4.f);
  result.fmDirectLevel = Safe(source.fmDirectLevel, .0144f, 0.f, 3.f);
  result.fmBodyLevel = Safe(source.fmBodyLevel, .081f, 0.f, 3.f);

  result.observation[0].gain = Safe(source.directLevel, .9f, 0.f, 4.f);
  result.observation[0].delaySeconds = Safe(
      source.directDelaySeconds, 0.f, 0.f, .01f);
  result.observation[0].radiationEnabled = false;
  result.observation[1].gain = Safe(source.bodyLevel, 3.f, 0.f, 4.f);
  result.observation[1].radiationEnabled = false;

  result.equalizer.mode = source.equalizerMode;
  result.equalizer.lowCutHz = Safe(source.lowCutHz, 24.f, 5.f, 1000.f);
  result.equalizer.highCutHz = Safe(
      source.highCutHz, 18000.f, 500.f, 30000.f);
  result.equalizer.radiation.lowCutHz = result.equalizer.lowCutHz;
  result.equalizer.radiation.highCutHz = result.equalizer.highCutHz;
  result.equalizer.radiation.colourFrequencyHz = Safe(
      source.colourFrequencyHz, 2800.f, 40.f, 20000.f);
  result.equalizer.radiation.colourGainDb = Safe(
      source.colourGainDb, 0.f, -24.f, 24.f);
  result.equalizer.radiation.colourQ = .7f;
  result.equalizer.outputGain = 1.f;
  result.outputGain = Safe(source.outputGain, .5f, 0.f, 4.f);
  return result;
}

MembraneDrumPreparedParameters PrepareMembraneDrumParameters(
    const float sampleRate, const MembraneDrumParameters &parameters) {
  MembraneDrumPreparedParameters result;
  result.parameters = parameters;
  result.membrane = MembraneResonator<MembraneModeCount>::PrepareParameters(
      sampleRate, parameters.membrane);
  result.sampleRate = sampleRate;
  return result;
}

} // namespace tfdsp::percussion
