#pragma once

#include "parameter_descriptor.hpp"
#include "tfdsp/percussion/compact_kick_parameters.hpp"

#include <array>
#include <cstddef>

namespace tfworkbench {

enum class KickParameter : std::size_t {
  ModelLevelDb,
  PrimaryLevel,
  FundamentalHz,
  PitchDropOctaves,
  PitchDecaySeconds,
  BodyDecaySeconds,
  FmDepthHz,
  FmDecaySeconds,
  FmRoughnessHz,
  SecondaryRatio,
  SecondaryLevel,
  ClickLevel,
  ClickDecaySeconds,
  ClickTiltDb,
  LowCutHz,
  HighCutHz,
  Count
};

inline constexpr std::size_t KickParameterCount =
    static_cast<std::size_t>(KickParameter::Count);
using KickParameterValues = std::array<float, KickParameterCount>;

const ParameterDescriptor &KickParameterDescription(std::size_t index) noexcept;
KickParameterValues DefaultKickParameters() noexcept;
tfdsp::percussion::CompactKickParameters ApplyKickParameters(
    const KickParameterValues &values) noexcept;

} // namespace tfworkbench
