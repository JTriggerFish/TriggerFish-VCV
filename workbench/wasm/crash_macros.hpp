#pragma once

#include "tfdsp/percussion/crash_cymbal_parameters.hpp"

#include <array>
#include <cstddef>

namespace tfworkbench {

enum class CrashMacro : std::size_t {
  ModelLevelDb,
  ImpactToneNoise,
  ImpactWidth,
  BloomAmount,
  BloomDevelopment,
  BodyToneWash,
  BodyBrightness,
  LowDecay,
  MiddleDecay,
  HighDecay,
  Count
};

inline constexpr std::size_t CrashMacroCount =
    static_cast<std::size_t>(CrashMacro::Count);

struct CrashMacroDescriptor {
  const char *name;
  const char *unit;
  float minimum;
  float maximum;
  float defaultValue;
};

using CrashMacroValues = std::array<float, CrashMacroCount>;

const CrashMacroDescriptor &CrashMacroDescription(std::size_t index) noexcept;
CrashMacroValues DefaultCrashMacros() noexcept;
tfdsp::percussion::CrashCymbalFitParameters ApplyCrashMacros(
    const tfdsp::percussion::CrashCymbalFitParameters &base,
    const CrashMacroValues &values) noexcept;

} // namespace tfworkbench
