#pragma once

#include "parameter_descriptor.hpp"
#include "kick_mode_macros.hpp"
#include "tfdsp/percussion/kick_voice_parameters.hpp"
#include <array>
#include <cstddef>

namespace tfworkbench {
enum class KickParameter : std::size_t {
  ModelLevelDb,
  ContactLevel,
  ContactWidth,
  ContactColour,
  ContactNoise,
  ContactNoiseDecay,
  ThumpLevel,
  ThumpPitch,
  ThumpDrop,
  ThumpFall,
  ThumpDecay,
  ResonanceLevel,
  ResonanceDecay,
  ResonanceDecayTilt,
  TensionOctaves,
  TensionRecovery,
  EqualizerMode,
  LowCutHz,
  HighCutHz,
  ColourFrequency,
  ColourGain,
  Band1Frequency,
  Band1Gain,
  Band2Frequency,
  Band2Gain,
  Band3Frequency,
  Band3Gain,
  Band4Frequency,
  Band4Gain,
  Count
};
inline constexpr std::size_t KickParameterCount =
    static_cast<std::size_t>(KickParameter::Count) + KickModeParameterCount;
using KickParameterValues = std::array<float, KickParameterCount>;
const ParameterDescriptor &KickParameterDescription(std::size_t index) noexcept;
KickParameterValues DefaultKickParameters() noexcept;
tfdsp::percussion::MembraneDrumParameters
ApplyKickParameters(const KickParameterValues &values) noexcept;
} // namespace tfworkbench
