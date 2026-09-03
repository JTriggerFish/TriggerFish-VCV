#pragma once

#include "parameter_descriptor.hpp"
#include "tfdsp/percussion/crash_cymbal_parameters.hpp"

#include <array>
#include <cstddef>
#include <string>

namespace tfworkbench {

inline constexpr std::size_t ResolvedModePointCount =
    tfdsp::percussion::CrashSparseModeCount;
inline constexpr std::size_t DenseWashCurvePointCount = 8;
inline constexpr std::size_t TurbulenceCurvePointCount = 3;
inline constexpr std::size_t BodyDecayCurvePointCount =
    tfdsp::percussion::CrashBodyDecayPointCount;

using CrashMacroScale = ParameterScale;

enum class CrashMacro : std::size_t {
  ModelLevelDb,
  ImpactToneNoise,
  ImpactWidth,
  BloomLevel,
  BloomNonlinearity,
  BloomDevelopment,
  BloomDiffusion,
  BodyToneWash,
  BodyBrightness,
  UnifiedBodyEnabled,
  FieldTurbulence,
  FieldPacketSpread,
  FieldPhaseBandwidth,
  FieldExchange,
  FieldGain,
  DirectGain,
  TurbulenceEnabled,
  TurbulenceAmount,
  TurbulencePersistence,
  TurbulenceFrequencyFirst,
  TurbulenceLevelFirst =
      TurbulenceFrequencyFirst + TurbulenceCurvePointCount,
  DirectRadiationEnabled =
      TurbulenceLevelFirst + TurbulenceCurvePointCount,
  DirectLowCut,
  DirectLowCutQ,
  DirectColourFrequency,
  DirectColourGain,
  DirectColourQ,
  DirectHighCut,
  DirectHighCutQ,
  SparseRadiationEnabled,
  SparseLowCut,
  SparseLowCutQ,
  SparseColourFrequency,
  SparseColourGain,
  SparseColourQ,
  SparseHighCut,
  SparseHighCutQ,
  DenseRadiationEnabled,
  DenseLowCut,
  DenseLowCutQ,
  DenseColourFrequency,
  DenseColourGain,
  DenseColourQ,
  DenseHighCut,
  DenseHighCutQ,
  DenseMinimumFrequency,
  DenseMaximumFrequency,
  DenseModeDensity,
  DenseSpacingJitter,
  DenseDecaySpread,
  DenseGainSpread,
  DenseWashFrequencyFirst,
  DenseWashLevelFirst =
      DenseWashFrequencyFirst + DenseWashCurvePointCount,
  BodyDecayFrequencyFirst =
      DenseWashLevelFirst + DenseWashCurvePointCount,
  BodyDecaySecondsFirst =
      BodyDecayFrequencyFirst + BodyDecayCurvePointCount,
  BodyDecayActiveFirst =
      BodyDecaySecondsFirst + BodyDecayCurvePointCount,
  ResolvedFrequencyFirst =
      BodyDecayActiveFirst + BodyDecayCurvePointCount,
  ResolvedLevelFirst = ResolvedFrequencyFirst + ResolvedModePointCount,
  ResolvedTurbulenceFirst = ResolvedLevelFirst + ResolvedModePointCount,
  Count = ResolvedTurbulenceFirst + ResolvedModePointCount
};

inline constexpr std::size_t CrashMacroCount =
    static_cast<std::size_t>(CrashMacro::Count);

// The compatibility API retains the historical positional surface so old
// fitting tools can still be migrated. Production recipes expose only values
// consumed by the active unified metallic model.
inline constexpr bool IsActiveCrashMacro(const std::size_t index) noexcept {
  return index <= static_cast<std::size_t>(CrashMacro::BloomDiffusion) ||
      index == static_cast<std::size_t>(CrashMacro::BodyBrightness) ||
      (index >= static_cast<std::size_t>(CrashMacro::FieldTurbulence) &&
       index <= static_cast<std::size_t>(CrashMacro::DirectGain)) ||
      (index >= static_cast<std::size_t>(CrashMacro::DirectRadiationEnabled) &&
       index <= static_cast<std::size_t>(CrashMacro::DirectHighCutQ)) ||
      (index >= static_cast<std::size_t>(CrashMacro::DenseRadiationEnabled) &&
       index <= static_cast<std::size_t>(CrashMacro::DenseHighCutQ)) ||
      index >= static_cast<std::size_t>(CrashMacro::BodyDecayFrequencyFirst);
}

inline constexpr std::size_t CountActiveCrashMacros() noexcept {
  std::size_t result = 0;
  for (std::size_t index = 0; index < CrashMacroCount; ++index)
    result += IsActiveCrashMacro(index) ? 1 : 0;
  return result;
}

inline constexpr std::size_t ActiveCrashMacroCount =
    CountActiveCrashMacros();

inline constexpr std::array<std::size_t, ActiveCrashMacroCount>
BuildActiveCrashMacroIndices() noexcept {
  std::array<std::size_t, ActiveCrashMacroCount> result{};
  std::size_t target = 0;
  for (std::size_t source = 0; source < CrashMacroCount; ++source) {
    if (IsActiveCrashMacro(source)) result[target++] = source;
  }
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
tfdsp::percussion::CrashCymbalFitParameters CrashWorkbenchBaseFit() noexcept;
tfdsp::percussion::CrashCymbalFitParameters ApplyCrashMacros(
    const tfdsp::percussion::CrashCymbalFitParameters &base,
    const CrashMacroValues &values) noexcept;

} // namespace tfworkbench
