#pragma once

#include "parameter_descriptor.hpp"
#include "tfdsp/percussion/crash_cymbal_parameters.hpp"

#include <array>
#include <cstddef>
#include <string>

namespace tfworkbench {

inline constexpr std::size_t ResolvedModePointCount =
    tfdsp::percussion::CrashModalAnchorCapacity;
inline constexpr std::size_t BodyDecayCurvePointCount =
    tfdsp::percussion::CrashBodyDecayPointCount;
inline constexpr std::size_t BodyDecayInteriorPointCount =
    tfdsp::percussion::CrashBodyDecayInteriorPointCount;

using CrashMacroScale = ParameterScale;

enum class CrashMacro : std::size_t {
  ModelLevelDb,
  ImpactToneNoise,
  ImpactWidth,
  BloomRate,
  BloomEnergyAcceleration,
  BloomPhaseDiffusion,
  BodyBrightness,
  BodyExcitationCentre,
  FieldTurbulence,
  FieldTurbulenceSlope,
  FieldTurbulenceCentre,
  FieldPacketSpread,
  FieldSatelliteDensity,
  FieldPhaseBandwidth,
  FieldExchange,
  BodyExcitation,
  FieldGain,
  DirectGain,
  DirectRadiationEnabled,
  DirectLowCut,
  DirectColourFrequency,
  DirectColourGain,
  DirectHighCut,
  BodyRadiationEnabled,
  BodyLowCut,
  BodyColourFrequency,
  BodyColourGain,
  BodyHighCut,
  BodyDecayFrequencyFirst,
  BodyDecaySecondsFirst =
      BodyDecayFrequencyFirst + BodyDecayInteriorPointCount,
  BodyDecayActiveFirst =
      BodyDecaySecondsFirst + BodyDecayCurvePointCount,
  ResolvedFrequencyFirst =
      BodyDecayActiveFirst + BodyDecayInteriorPointCount,
  ResolvedLevelFirst = ResolvedFrequencyFirst + ResolvedModePointCount,
  ResolvedTurbulenceFirst = ResolvedLevelFirst + ResolvedModePointCount,
  ImpactChirpPitch = ResolvedTurbulenceFirst + ResolvedModePointCount,
  ImpactNoiseTilt,
  ImpactMicroDensity,
  VelocityBrightness,
  BodyTune,
  Count
};

inline constexpr std::size_t CrashMacroCount =
    static_cast<std::size_t>(CrashMacro::Count);

inline constexpr std::size_t ActiveCrashMacroCount =
    CrashMacroCount;

inline constexpr std::array<std::size_t, ActiveCrashMacroCount>
BuildActiveCrashMacroIndices() noexcept {
  std::array<std::size_t, ActiveCrashMacroCount> result{};
  for (std::size_t index = 0; index < CrashMacroCount; ++index)
    result[index] = index;
  return result;
}

inline constexpr auto ActiveCrashMacroIndices =
    BuildActiveCrashMacroIndices();

using CrashMacroDescriptor = ParameterDescriptor;

using CrashMacroValues = std::array<float, CrashMacroCount>;

const CrashMacroDescriptor &CrashMacroDescription(std::size_t index) noexcept;
const CrashMacroDescriptor &ActiveCrashMacroDescription(
    std::size_t index) noexcept;
CrashMacroValues DefaultCrashMacros() noexcept;
tfdsp::percussion::CrashCymbalFitParameters MetallicWorkbenchBaseFit() noexcept;
tfdsp::percussion::CrashCymbalFitParameters ApplyCrashMacros(
    const tfdsp::percussion::CrashCymbalFitParameters &base,
    const CrashMacroValues &values) noexcept;

} // namespace tfworkbench
