#include "snare_drum.hpp"

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {

float SnareDrumRouting::Get(const SnareDrumRoute route) const noexcept {
  return gains[static_cast<std::size_t>(route)];
}

void SnareDrumRouting::Set(const std::size_t index, const float gain) noexcept {
  if (index < gains.size())
    gains[index] = std::clamp(tfdsp::FiniteNormalOrZero(gain), 0.f, 2.f);
}

SnareDrumParameters DefaultSnareDrumParameters() noexcept {
  SnareDrumParameters result;
  MembraneDrumControls body;
  body.fundamentalHz = 185.f;
  body.decaySeconds = .12f;
  body.decayTilt = .45f;
  body.inharmonicity = .72f;
  body.bodyBrightness = .68f;
  body.tensionOctaves = .08f;
  body.tensionDecaySeconds = .09f;
  body.contactLevel = .78f;
  body.contactDurationSeconds = .0022f;
  body.contactBrightness = .38f;
  body.fmLevel = .08f;
  body.fmDepthHz = 180.f;
  body.fmDecaySeconds = .035f;
  body.pitchDropOctaves = .12f;
  result.membrane = DefaultMembraneDrumParameters(body);
  result.membrane.directVelocityExponent = 2.59f;
  result.membrane.bodyVelocityExponent = 2.23f;
  result.membrane.velocitySaturation = 3.64f;
  // Keep the safety ceiling inactive for normal hits. The generic ceiling is
  // intentionally conservative, but it compressed all snare velocities to
  // nearly the same stored body energy.
  result.membrane.maximumModalEnergy = 16.f;
  result.membrane.membrane[5].frequencyHz = 675.f;
  result.membrane.membrane[5].decaySeconds = 1.8f;
  result.membrane.membrane[5].outputGain = .3f;
  result.membrane.observation[0].gain = 1.f;
  result.membrane.observation[1].gain = 1.f;
  result.membrane.equalizer.mode = ObservationEqualizerMode::Bypass;
  result.membrane.outputGain = 1.f;
  result.membrane.routing.gains = {.35f, 1.f, .05f, .32f, 1.f};

  result.wires.sensitivity = 4.f;
  result.wires.attackSeconds = .045f;
  result.wires.releaseSeconds = .03f;
  result.wires.threshold = .05f;
  result.wires.motionHighpassHz = 1000.f;
  result.wires.minimumFrequencyHz = 520.f;
  result.wires.maximumFrequencyHz = 9000.f;
  result.wires.decaySeconds = .025f;
  result.wires.decayTilt = .78f;
  result.wires.density = .9f;
  result.wires.brightness = .05f;
  result.wires.noiseMix = .42f;
  result.wires.modalMix = .7f;
  result.wires.outputGain = 2.5f;

  result.observation[0].gain = .06f;
  result.observation[1].gain = .95f;
  result.observation[2].gain = 4.f;
  for (auto &path : result.observation)
    path.radiationEnabled = false;
  result.equalizer.mode = ObservationEqualizerMode::Radiation;
  result.equalizer.radiation.lowCutHz = 38.f;
  result.equalizer.radiation.highCutHz = 10000.f;
  result.equalizer.radiation.colourFrequencyHz = 210.f;
  result.equalizer.radiation.colourGainDb = 0.f;
  result.equalizer.radiation.colourQ = .7f;
  result.equalizer.outputGain = 1.f;
  result.outputGain = .1f;
  return result;
}

SnareDrumPreparedParameters
PrepareSnareDrumParameters(const float sampleRate,
                           const SnareDrumParameters &source) {
  SnareDrumPreparedParameters result;
  auto membrane = source.membrane;
  membrane.routing.gains = {source.routing.Get(SnareDrumRoute::ContactToDirect),
                            source.routing.Get(SnareDrumRoute::ContactToBody),
                            source.routing.Get(SnareDrumRoute::FmToDirect),
                            source.routing.Get(SnareDrumRoute::FmToBody), 1.f};
  result.membrane = PrepareMembraneDrumParameters(sampleRate, membrane);
  result.wires = PrepareWireRackParameters(sampleRate, source.wires);
  result.observation = source.observation;
  result.equalizer = source.equalizer;
  result.routing = source.routing;
  result.outputGain =
      std::clamp(tfdsp::FiniteNormalOrZero(source.outputGain), 0.f, 4.f);
  result.sampleRate = sampleRate;
  return result;
}

void SnareDrum::Prepare(const float sampleRate,
                        const SnareDrumParameters &parameters) {
  Prepare(PrepareSnareDrumParameters(sampleRate, parameters));
}

void SnareDrum::Prepare(const SnareDrumPreparedParameters &prepared) {
  membrane_.Prepare(prepared.membrane);
  wires_.Prepare(prepared.wires);
  observation_.Prepare(prepared.sampleRate, .05f, prepared.observation);
  equalizer_.Prepare(prepared.sampleRate, prepared.equalizer);
  routing_ = prepared.routing;
  outputGain_ = prepared.outputGain;
  Reset();
}

void SnareDrum::Reset() noexcept {
  membrane_.Reset();
  wires_.Reset();
  observation_.Reset();
  equalizer_.Reset();
}

void SnareDrum::Trigger(const MembraneDrumHit &hit) noexcept {
  membrane_.Trigger(hit);
  wires_.Seed(hit.seed ^ 0xa511e9b3u);
}

SnareDrumFrame SnareDrum::ProcessFrame() noexcept {
  const auto membrane = membrane_.ProcessSources();
  const float wireDrive =
      routing_.Get(SnareDrumRoute::BodyToWires) * membrane.body;
  const float wires = wires_.Process(wireDrive);
  const float observed = observation_.Process(
      {membrane.direct,
       routing_.Get(SnareDrumRoute::BodyToObservation) * membrane.body,
       routing_.Get(SnareDrumRoute::WiresToObservation) * wires});
  const float output =
      tfdsp::FiniteNormalOrZero(outputGain_ * equalizer_.Process(observed));
  return {membrane.direct, membrane.body, wires, output};
}

float SnareDrum::Process() noexcept { return ProcessFrame().output; }

} // namespace tfdsp::percussion
