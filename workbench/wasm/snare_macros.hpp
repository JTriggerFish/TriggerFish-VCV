#pragma once

#include "membrane_macros.hpp"
#include "parameter_descriptor.hpp"
#include "tfdsp/percussion/snare_drum.hpp"

#include <array>
#include <cstddef>

namespace tfworkbench {

enum class SnareParameter : std::size_t {
  WireLevel = MembraneParameterCount,
  WireSensitivity, WireThreshold, WireAttackSeconds, WireReleaseSeconds,
  WireMinimumHz, WireMaximumHz, WireDecaySeconds, WireDecayTilt,
  WireDensity, WireBrightness, WireNoiseMix, WireModalMix,
  RingFrequencyHz, RingDecaySeconds, RingLevel, Count
};

inline constexpr std::size_t SnareParameterCount =
    static_cast<std::size_t>(SnareParameter::Count);
using SnareParameterValues = std::array<float, SnareParameterCount>;

const ParameterDescriptor &SnareParameterDescription(
    std::size_t index) noexcept;
SnareParameterValues DefaultSnareParameters() noexcept;
tfdsp::percussion::SnareDrumParameters ApplySnareParameters(
    const SnareParameterValues &values) noexcept;

} // namespace tfworkbench
