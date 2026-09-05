#include "kick_voice_parameters.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {
namespace {
float Safe(float value, float fallback, float low, float high) noexcept {
  return std::clamp(std::isfinite(value) ? value : fallback, low, high);
}
} // namespace

std::array<KickModeControl, MembraneModeCount> DefaultKickModes() noexcept {
  // Factory data, not a constraint on the sounding bank. Editor templates can
  // replace this entire list; no runtime harmonic/membrane interpolation remains.
  return {{{55.f, 0.f, 1.f, .18f}, {87.615f, -4.03f, -.63f, .235f},
           {117.48f, -6.57f, .611f, .289f}, {126.28f, -7.19f, -.592f, .344f},
           {145.915f, -8.44f, .572f, .399f}, {160.49f, -9.27f, -.552f, .453f},
           {173.58f, -9.95f, .533f, .508f}, {192.555f, -10.85f, -.514f, .563f},
           {198.f, -11.09f, .494f, .617f}, {200.86f, -11.21f, -.475f, .672f},
           {223.3f, -12.13f, .455f, .727f}, {228.47f, -12.33f, -.435f, .781f},
           {236.005f, -12.61f, .416f, .836f}, {253.055f, -13.21f, -.397f, .891f},
           {265.705f, -13.63f, .377f, .945f}, {282.48f, -14.15f, -.357f, 1.f}}};
}

MembraneDrumParameters
DefaultKickVoiceParameters(const KickVoiceControls &source) noexcept {
  MembraneDrumControls body;
  body.fundamentalHz = Safe(source.thumpPitchHz, 28.f, 20.f, 180.f);
  body.decaySeconds = source.resonanceDecaySeconds;
  body.decayTilt = source.resonanceDecayTilt;
  body.tensionOctaves = source.tensionOctaves;
  body.tensionDecaySeconds = source.tensionRecoverySeconds;
  body.contactDirectLevel = source.contactLevel;
  body.contactBodyLevel = 1.f;
  body.contactDurationSeconds = source.contactWidthSeconds;
  body.contactBrightness = source.contactColour;
  body.contactNoiseLevel = source.contactNoiseLevel;
  body.contactNoiseDecaySeconds =
      4.f / 3.f * Safe(source.contactNoiseDecaySeconds, .205f, .00075f, .75f);
  body.fmDirectLevel = source.thumpLevel;
  body.fmBodyLevel = 0.f;
  body.fmDepthHz = 0.f;
  body.fmPitchDecaySeconds = source.thumpPitchFallSeconds;
  body.pitchDropOctaves = source.thumpPitchDropOctaves;
  body.directLevel = 1.f;
  body.bodyLevel = 1.f;
  body.equalizerMode = ObservationEqualizerMode::Bypass;
  body.outputGain = Safe(source.outputGain, .25f, 0.f, 1.f);
  auto result = DefaultMembraneDrumParameters(body);
  // T60 means -60 dB amplitude. The finite source closes after -80 dB.
  result.fm.amplitude.segments[1].durationSeconds =
      4.f / 3.f * Safe(source.thumpDecaySeconds, .306f, .005f, 3.f);
  result.fm.carrierFrequencyHz.initialValue =
      body.fundamentalHz *
      std::exp2(Safe(source.thumpPitchDropOctaves, 1.47f, 0.f, 4.f));
  result.fm.carrierFrequencyHz.segments[0].targetValue = body.fundamentalHz;
  result.fm.pitchPerturbationCents =
      0.f; // stable thump; contact supplies variation
  result.contact.noise.decaySeconds = body.contactNoiseDecaySeconds;
  const float decay = Safe(source.resonanceDecaySeconds, .6f, .03f, 8.f);
  const float slope = Safe(source.resonanceDecayTilt, .5f, -1.f, 1.f);
  for (std::size_t i = 0; i < result.membrane.size(); ++i) {
    const auto &control = source.modes[i];
    auto &mode = result.membrane[i];
    mode.frequencyHz = Safe(control.frequencyHz, 55.f, 20.f, 15000.f);
    const float level = Safe(control.levelDb, -72.f, -72.f, 6.f);
    // Share prominence between passive drive and observation. Subtract the
    // displayed off-floor so both approach zero continuously at -72 dB;
    // deleting a quiet mode must not suddenly redistribute a full drive slot.
    const float prominence = std::max(
        0.f, std::pow(10.f, level / 20.f) - std::pow(10.f, -72.f / 20.f));
    mode.inputGain = mode.outputGain = std::sqrt(prominence);
    mode.centerProjection = Safe(control.centreCoupling, 1.f, -1.f, 1.f);
    mode.edgeProjection = Safe(control.edgeCoupling, 1.f, -1.f, 1.f);
    // One frequency-based damping law, independent of slot/order/mode count.
    mode.decaySeconds = std::clamp(
        decay * std::pow(mode.frequencyHz / 100.f, -slope), .002f, 30.f);
  }
  result.observation[1].gain = Safe(source.resonanceLevel, 4.72f, 0.f, 12.f);
  result.fmDirectLevel = Safe(source.thumpLevel, 1.77f, 0.f, 4.f);
  ApplyKickRouting(result, {});
  return result;
}

void ApplyKickRouting(MembraneDrumParameters &parameters,
                      const KickVoiceRouting &routing) noexcept {
  parameters.routing.enabled = {routing.enabled[0], true,
                                routing.enabled[1], false, routing.enabled[2]};
}

} // namespace tfdsp::percussion
