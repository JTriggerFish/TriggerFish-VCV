#include "snare_drum.hpp"

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {

bool SnareDrumRouting::Enabled(const SnareDrumRoute route) const noexcept {
  return enabled[static_cast<std::size_t>(route)];
}

void SnareDrumRouting::SetEnabled(const std::size_t index,
                                  const bool value) noexcept {
  if (index < enabled.size()) enabled[index] = value;
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
  body.contactDirectLevel = .273f;
  body.contactBodyLevel = .78f;
  body.contactDurationSeconds = .0022f;
  body.contactBrightness = .25f;
  body.fmDirectLevel = .004f;
  body.fmBodyLevel = .0256f;
  body.fmDepthHz = 180.f;
  body.fmDecaySeconds = .035f;
  body.pitchDropOctaves = .12f;
  result.membrane = DefaultMembraneDrumParameters(body);
  result.membrane.membrane[5].frequencyHz = 675.f;
  result.membrane.membrane[5].decaySeconds = 1.8f;
  result.membrane.membrane[5].outputGain = .3f;
  result.membrane.observation[0].gain = 1.f;
  result.membrane.observation[1].gain = 1.f;
  result.membrane.equalizer.mode = ObservationEqualizerMode::Bypass;
  result.membrane.outputGain = 1.f;

  result.wires.sensitivity = 1.5f;
  result.wires.attackSeconds = .0015f;
  result.wires.releaseSeconds = .08f;
  result.wires.threshold = .0015f;
  result.wires.motionHighpassHz = 140.f;
  result.wires.minimumFrequencyHz = 250.f;
  result.wires.maximumFrequencyHz = 14000.f;
  result.wires.decaySeconds = .35f;
  result.wires.decayTilt = .8f;
  result.wires.density = .9f;
  result.wires.brightness = 0.f;
  result.wires.noiseMix = .6f;
  result.wires.modalMix = .45f;

  result.observation[0].gain = .24f;
  result.observation[1].gain = 2.85f;
  result.observation[2].gain = 1.38f;
  for (auto &path : result.observation)
    path.radiationEnabled = false;
  result.equalizer.mode = ObservationEqualizerMode::Radiation;
  result.equalizer.radiation.lowCutHz = 38.f;
  result.equalizer.radiation.highCutHz = 8000.f;
  result.equalizer.radiation.colourFrequencyHz = 210.f;
  result.equalizer.radiation.colourGainDb = 0.f;
  result.equalizer.radiation.colourQ = .7f;
  result.equalizer.outputGain = 1.f;
  result.outputGain = .5f;
  return result;
}

SnareDrumPreparedParameters
PrepareSnareDrumParameters(const float sampleRate,
                           const SnareDrumParameters &source) {
  SnareDrumPreparedParameters result;
  auto membrane = source.membrane;
  membrane.routing.SetEnabled(
      static_cast<std::size_t>(MembraneDrumRoute::ContactToDirect),
      source.routing.Enabled(SnareDrumRoute::ContactToDirect));
  membrane.routing.SetEnabled(
      static_cast<std::size_t>(MembraneDrumRoute::ContactToBody),
      source.routing.Enabled(SnareDrumRoute::ContactToBody));
  membrane.routing.SetEnabled(
      static_cast<std::size_t>(MembraneDrumRoute::FmToDirect),
      source.routing.Enabled(SnareDrumRoute::FmToDirect));
  membrane.routing.SetEnabled(
      static_cast<std::size_t>(MembraneDrumRoute::FmToBody),
      source.routing.Enabled(SnareDrumRoute::FmToBody));
  membrane.routing.SetEnabled(
      static_cast<std::size_t>(MembraneDrumRoute::BodyToObservation), true);
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
      (routing_.Enabled(SnareDrumRoute::BodyToWires) ? 1.f : 0.f) *
      membrane.body;
  const float wires = wires_.Process(wireDrive);
  const float observed = observation_.Process(
      {membrane.direct,
       (routing_.Enabled(SnareDrumRoute::BodyToObservation) ? 1.f : 0.f) *
           membrane.body,
       (routing_.Enabled(SnareDrumRoute::WiresToObservation) ? 1.f : 0.f) *
           wires});
  const float output =
      tfdsp::FiniteNormalOrZero(outputGain_ * equalizer_.Process(observed));
  return {membrane.direct, membrane.body, wires, output};
}

float SnareDrum::Process() noexcept { return ProcessFrame().output; }

} // namespace tfdsp::percussion
