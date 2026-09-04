#include "compact_kick.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {
namespace {

float Bounded(const float value, const float fallback,
              const float low, const float high) noexcept {
  return std::clamp(std::isfinite(value) ? value : fallback, low, high);
}

void ScaleDeviation(CorrelatedFmTrajectory &trajectory,
                    const float scale) noexcept {
  trajectory.initialValue *= scale;
  for (std::size_t index = 0; index < trajectory.segmentCount; ++index)
    trajectory.segments[index].targetValue *= scale;
}

} // namespace

void CompactKick::EventVoice::Prepare(const float sampleRate) {
  primary.Prepare(sampleRate);
  secondary.Prepare(sampleRate);
  click.Prepare(sampleRate);
  stealStep = 1.f / std::max(1.f, .001f * sampleRate);
  activityRelease = std::exp(-1.f / std::max(1.f, .01f * sampleRate));
}

void CompactKick::EventVoice::Reset() noexcept {
  primary.Reset();
  secondary.Reset();
  click.Reset();
  amplitude = 0.f;
  lastPrimary = lastSecondary = lastClick = 0.f;
  stealPrimary = stealSecondary = stealClick = 0.f;
  stealGain = activity = 0.f;
  generation = 0;
}

void CompactKick::EventVoice::Trigger(
    const CompactKickParameters &parameters, const CompactKickHit &hit,
    const std::uint64_t eventGeneration) noexcept {
  if (Active()) {
    stealPrimary = lastPrimary;
    stealSecondary = lastSecondary;
    stealClick = lastClick;
    stealGain = 1.f;
  } else {
    stealPrimary = stealSecondary = stealClick = 0.f;
    stealGain = 0.f;
  }
  const float strength = Bounded(hit.strength, .8f, 0.f, 1.f);
  const float hardness = Bounded(hit.hardness, .5f, 0.f, 1.f);
  auto primaryParameters = parameters.primary;
  auto secondaryParameters = parameters.secondary;
  auto clickParameters = parameters.click;
  const float pitchScale = std::exp2(.28f * (strength - .8f));
  primaryParameters.carrierFrequencyHz.initialValue *= pitchScale;
  secondaryParameters.carrierFrequencyHz.initialValue *= pitchScale;
  ScaleDeviation(primaryParameters.frequencyDeviationHz, .35f + .8f * strength);
  ScaleDeviation(secondaryParameters.frequencyDeviationHz,
                 .25f + .55f * strength + .8f * hardness);
  primaryParameters.seed = hit.seed;
  secondaryParameters.seed = hit.seed;
  clickParameters.seed = hit.seed ^ 0x9e3779b9u;
  clickParameters.amplitude = .3f + .7f * hardness;
  primary.Trigger(primaryParameters);
  secondary.Trigger(secondaryParameters);
  click.Trigger(clickParameters);
  amplitude = strength;
  primaryLevel = parameters.primaryLevel;
  secondaryLevel = parameters.secondaryLevel;
  clickLevel = parameters.clickLevel * amplitude * (.2f + .8f * strength) *
      (.05f + 2.5f * hardness);
  generation = eventGeneration;
}

CompactKickFrame CompactKick::EventVoice::Process() noexcept {
  float primarySample = primaryLevel * amplitude * primary.Process();
  float secondarySample = secondaryLevel * amplitude * secondary.Process();
  float clickSample = clickLevel * click.Process();
  if (stealGain > 0.f) {
    primarySample += stealGain * stealPrimary;
    secondarySample += stealGain * stealSecondary;
    clickSample += stealGain * stealClick;
    stealGain = std::max(0.f, stealGain - stealStep);
  }
  lastPrimary = primarySample;
  lastSecondary = secondarySample;
  lastClick = clickSample;
  const float peak = std::max({std::abs(primarySample),
                               std::abs(secondarySample),
                               std::abs(clickSample)});
  activity = std::max(peak, activityRelease * activity);
  return {primarySample, secondarySample, clickSample, 0.f};
}

bool CompactKick::EventVoice::Active() const noexcept {
  return primary.Active() || secondary.Active() || click.Active() ||
      stealGain > 0.f;
}

float CompactKick::EventVoice::Activity() const noexcept {
  return activity;
}

void CompactKick::Prepare(const float sampleRate,
                          const CompactKickParameters &parameters) {
  parameters_ = parameters;
  for (auto &voice : voices_) voice.Prepare(sampleRate);
  observation_.Prepare(sampleRate, .01f, parameters.observation);
  sourceMixer_.SetGains(parameters.routing.gains);
  Reset();
}

void CompactKick::Reset() noexcept {
  for (auto &voice : voices_) voice.Reset();
  observation_.Reset();
  generation_ = 0;
}

CompactKick::EventVoice &CompactKick::NextVoice() noexcept {
  const auto inactive = std::find_if(
      voices_.begin(), voices_.end(), [](const auto &voice) {
        return !voice.Active();
      });
  if (inactive != voices_.end()) return *inactive;
  return *std::min_element(
      voices_.begin(), voices_.end(), [](const auto &left, const auto &right) {
        if (left.Activity() != right.Activity())
          return left.Activity() < right.Activity();
        return left.generation < right.generation;
      });
}

void CompactKick::Trigger(const CompactKickHit &hit) noexcept {
  if (Bounded(hit.strength, .8f, 0.f, 1.f) <= 0.f) return;
  NextVoice().Trigger(parameters_, hit, ++generation_);
}

CompactKickFrame CompactKick::ProcessFrame() noexcept {
  CompactKickFrame frame{};
  for (auto &voice : voices_) {
    const auto event = voice.Process();
    frame.primary += event.primary;
    frame.secondary += event.secondary;
    frame.click += event.click;
  }
  const float mix = sourceMixer_.Process(
      {frame.primary, frame.secondary, frame.click});
  frame.output = parameters_.outputGain * observation_.Process({mix});
  frame.output = tfdsp::FiniteNormalOrZero(frame.output);
  return frame;
}

float CompactKick::Process() noexcept {
  return ProcessFrame().output;
}

} // namespace tfdsp::percussion
