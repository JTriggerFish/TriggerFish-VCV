#pragma once

#include "tfdsp/percussion/crash_cymbal_parameters.hpp"

#include <array>
#include <cstddef>
#include <string>

namespace tfworkbench {

inline constexpr std::size_t ResolvedModePointCount = 12;
inline constexpr std::size_t DenseWashCurvePointCount = 8;
inline constexpr std::size_t TurbulenceCurvePointCount = 3;
inline constexpr std::size_t BodyDecayCurvePointCount =
    tfdsp::percussion::CrashBodyDecayPointCount;

enum class CrashMacroScale : int { Linear, Logarithmic, Boolean };

enum class CrashMacro : std::size_t {
  ModelLevelDb,
  ImpactToneNoise,
  ImpactWidth,
  BloomAmount,
  BloomDevelopment,
  BodyToneWash,
  BodyBrightness,
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
  ResolvedFrequencyFirst =
      BodyDecaySecondsFirst + BodyDecayCurvePointCount,
  ResolvedLevelFirst = ResolvedFrequencyFirst + ResolvedModePointCount,
  Count = ResolvedLevelFirst + ResolvedModePointCount
};

inline constexpr std::size_t CrashMacroCount =
    static_cast<std::size_t>(CrashMacro::Count);

struct CrashMacroDescriptor {
  std::string key;
  std::string name;
  std::string unit;
  float minimum;
  float maximum;
  float defaultValue;
  CrashMacroScale scale{CrashMacroScale::Linear};
};

using CrashMacroValues = std::array<float, CrashMacroCount>;

const CrashMacroDescriptor &CrashMacroDescription(std::size_t index) noexcept;
CrashMacroValues DefaultCrashMacros() noexcept;
tfdsp::percussion::CrashCymbalFitParameters ApplyCrashMacros(
    const tfdsp::percussion::CrashCymbalFitParameters &base,
    const CrashMacroValues &values) noexcept;

} // namespace tfworkbench
