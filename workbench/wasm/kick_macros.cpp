#include "kick_macros.hpp"

#include <array>
#include <cmath>

namespace tfworkbench {
namespace {

using Scale = ParameterScale;

const std::array<ParameterDescriptor, KickParameterCount> Descriptors{{
    {"model_level_db", "Model level", "dB", -60.f, 0.f, -12.f},
    {"primary_level", "Primary level", "x", 0.f, 2.f, 1.f},
    {"fundamental_hz", "Fundamental", "Hz", 25.f, 120.f, 52.f,
     Scale::Logarithmic},
    {"pitch_drop_octaves", "Pitch drop", "oct", 0.f, 4.f, 1.8f},
    {"pitch_decay_seconds", "Pitch decay", "s", .002f, .3f, .055f,
     Scale::Logarithmic},
    {"body_decay_seconds", "Body decay", "s", .03f, 2.f, .38f,
     Scale::Logarithmic},
    {"fm_depth_hz", "FM depth", "Hz", 0.f, 4000.f, 720.f},
    {"fm_decay_seconds", "FM decay", "s", .002f, .2f, .035f,
     Scale::Logarithmic},
    {"fm_roughness_hz", "FM roughness", "Hz", 100.f, 16000.f, 4200.f,
     Scale::Logarithmic},
    {"secondary_ratio", "Second body ratio", "x", .5f, 3.f, 1.52f},
    {"secondary_level", "Second body level", "x", 0.f, 1.5f, .32f},
    {"click_level", "Click level", "x", 0.f, 1.5f, .16f},
    {"click_decay_seconds", "Click decay", "s", .001f, .08f, .009f,
     Scale::Logarithmic},
    {"click_tilt_db", "Click tilt", "dB/oct", -12.f, 12.f, 3.f},
    {"low_cut_hz", "Low cut", "Hz", 5.f, 500.f, 18.f,
     Scale::Logarithmic},
    {"high_cut_hz", "High cut", "Hz", 1000.f, 22000.f, 15500.f,
     Scale::Logarithmic},
}};

std::size_t Index(const KickParameter parameter) noexcept {
  return static_cast<std::size_t>(parameter);
}

} // namespace

const ParameterDescriptor &KickParameterDescription(
    const std::size_t index) noexcept {
  return Descriptors[index < Descriptors.size() ? index : 0];
}

KickParameterValues DefaultKickParameters() noexcept {
  KickParameterValues result{};
  for (std::size_t index = 0; index < result.size(); ++index)
    result[index] = Descriptors[index].defaultValue;
  return result;
}

tfdsp::percussion::CompactKickParameters ApplyKickParameters(
    const KickParameterValues &values) noexcept {
  tfdsp::percussion::CompactKickControls controls;
  controls.primaryLevel = values[Index(KickParameter::PrimaryLevel)];
  controls.fundamentalHz = values[Index(KickParameter::FundamentalHz)];
  controls.pitchDropOctaves = values[Index(KickParameter::PitchDropOctaves)];
  controls.pitchDecaySeconds = values[Index(KickParameter::PitchDecaySeconds)];
  controls.bodyDecaySeconds = values[Index(KickParameter::BodyDecaySeconds)];
  controls.fmDepthHz = values[Index(KickParameter::FmDepthHz)];
  controls.fmDecaySeconds = values[Index(KickParameter::FmDecaySeconds)];
  controls.fmRoughnessHz = values[Index(KickParameter::FmRoughnessHz)];
  controls.secondaryRatio = values[Index(KickParameter::SecondaryRatio)];
  controls.secondaryLevel = values[Index(KickParameter::SecondaryLevel)];
  controls.clickLevel = values[Index(KickParameter::ClickLevel)];
  controls.clickDecaySeconds = values[Index(KickParameter::ClickDecaySeconds)];
  controls.clickTiltDb = values[Index(KickParameter::ClickTiltDb)];
  controls.lowCutHz = values[Index(KickParameter::LowCutHz)];
  controls.highCutHz = values[Index(KickParameter::HighCutHz)];
  controls.outputGain = std::pow(
      10.f, values[Index(KickParameter::ModelLevelDb)] / 20.f);
  return tfdsp::percussion::DefaultCompactKickParameters(controls);
}

} // namespace tfworkbench
