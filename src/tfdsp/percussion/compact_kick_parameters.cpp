#include "compact_kick_parameters.hpp"

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {
namespace {

using Curve = TrajectoryCurve;

float Finite(const float value, const float fallback) noexcept {
  return std::isfinite(value) ? value : fallback;
}

CorrelatedFmBurstParameters Burst(const float fundamentalHz,
                                  const float dropOctaves,
                                  const float pitchSeconds,
                                  const float decaySeconds,
                                  const float deviationHz,
                                  const float deviationSeconds,
                                  const float roughnessHz) noexcept {
  CorrelatedFmBurstParameters result;
  result.amplitude.initialValue = 0.f;
  result.amplitude.segments[0] = {1.f, .00035f, Curve::Linear};
  result.amplitude.segments[1] = {.0001f, decaySeconds, Curve::Geometric};
  result.amplitude.segments[2] = {0.f, .001f, Curve::Linear};
  result.amplitude.segmentCount = 3;
  result.carrierFrequencyHz.initialValue =
      fundamentalHz * std::exp2(dropOctaves);
  result.carrierFrequencyHz.segments[0] = {
      fundamentalHz, pitchSeconds, Curve::Geometric,
  };
  result.carrierFrequencyHz.segmentCount = 1;
  result.frequencyDeviationHz.initialValue = deviationHz;
  result.frequencyDeviationHz.segments[0] = {
      std::max(1.f, .02f * deviationHz), deviationSeconds, Curve::Geometric,
  };
  result.frequencyDeviationHz.segments[1] = {0.f, .002f, Curve::Linear};
  result.frequencyDeviationHz.segmentCount = 2;
  result.irregularCutoffHz = roughnessHz;
  result.periodicModulatorHz = 1.7f * fundamentalHz;
  result.periodicMix = .12f;
  result.pitchPerturbationCents = 2.f;
  result.deviationPerturbation = .025f;
  return result;
}

} // namespace

bool CompactKickRouting::Enabled(const CompactKickRoute route) const noexcept {
  return enabled[static_cast<std::size_t>(route)];
}

void CompactKickRouting::SetEnabled(const std::size_t index,
                                    const bool value) noexcept {
  if (index < enabled.size()) enabled[index] = value;
}

CompactKickParameters DefaultCompactKickParameters(
    const CompactKickControls &source) noexcept {
  CompactKickParameters result;
  const float fundamental = std::clamp(
      Finite(source.fundamentalHz, 52.f), 20.f, 180.f);
  const float drop = std::clamp(
      Finite(source.pitchDropOctaves, 1.8f), 0.f, 4.f);
  const float pitchSeconds = std::clamp(
      Finite(source.pitchDecaySeconds, .055f), .002f, .5f);
  const float decaySeconds = std::clamp(
      Finite(source.bodyDecaySeconds, .38f), .03f, 4.f);
  const float deviation = std::clamp(
      Finite(source.fmDepthHz, 720.f), 0.f, 6000.f);
  const float deviationSeconds = std::clamp(
      Finite(source.fmDecaySeconds, .035f), .002f, .5f);
  const float roughness = std::clamp(
      Finite(source.fmRoughnessHz, 4200.f), 20.f, 20000.f);
  const float ratio = std::clamp(
      Finite(source.secondaryRatio, 1.52f), .5f, 4.f);

  result.primary = Burst(fundamental, drop, pitchSeconds, decaySeconds,
                         deviation, deviationSeconds, roughness);
  result.secondary = Burst(
      ratio * fundamental, .72f * drop, .7f * pitchSeconds,
      .52f * decaySeconds, 1.2f * deviation, .7f * deviationSeconds,
      std::min(20000.f, 1.35f * roughness));
  result.secondary.carrierPhaseCycles = .19f;
  result.secondary.modulatorPhaseCycles = .31f;

  result.click.attackSeconds = .0001f;
  result.click.holdSeconds = .0003f;
  result.click.decaySeconds = std::clamp(
      Finite(source.clickDecaySeconds, .009f), .001f, .08f);
  result.click.tiltDb = std::clamp(
      Finite(source.clickTiltDb, 3.f), -12.f, 12.f);
  result.click.tiltPivotHz = 3000.f;

  auto &observation = result.observation[0];
  observation.radiationEnabled = true;
  observation.radiation.lowCutHz = std::clamp(
      Finite(source.lowCutHz, 18.f), 5.f, 500.f);
  observation.radiation.colourFrequencyHz = 95.f;
  observation.radiation.colourGainDb = 1.5f;
  observation.radiation.colourQ = .7f;
  observation.radiation.highCutHz = std::clamp(
      Finite(source.highCutHz, 15500.f), 1000.f, 22000.f);

  result.secondaryLevel = std::clamp(
      Finite(source.secondaryLevel, .32f), 0.f, 2.f);
  result.primaryLevel = std::clamp(Finite(source.primaryLevel, 1.f), 0.f, 2.f);
  result.clickLevel = std::clamp(Finite(source.clickLevel, .16f), 0.f, 2.f);
  result.outputGain = std::clamp(Finite(source.outputGain, .7f), 0.f, 4.f);
  return result;
}

} // namespace tfdsp::percussion
