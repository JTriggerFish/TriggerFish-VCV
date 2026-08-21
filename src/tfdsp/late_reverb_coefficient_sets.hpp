#pragma once

#include "late_reverb_coefficients.hpp"
#include "late_reverb_optimized_coefficients.hpp"

#include <array>
#include <cstddef>

namespace tfdsp {

enum class LateReverbFlavour { Base = 0, Optimized = 1 };
inline constexpr LateReverbFlavour DefaultLateReverbFlavour =
    LateReverbFlavour::Optimized;

static_assert(late_reverb_coefficients::LineCount ==
                  late_reverb_optimized_coefficients::LineCount,
              "late-reverb flavours must have the same line count");
static_assert(late_reverb_coefficients::VelvetStageCount ==
                  late_reverb_optimized_coefficients::VelvetStageCount,
              "late-reverb flavours must have the same velvet topology");

using LateMainDelayRatios =
    std::array<float, late_reverb_coefficients::LineCount>;
using LateVelvetDelays = std::array<
    std::array<float, late_reverb_coefficients::LineCount>,
    late_reverb_coefficients::VelvetStageCount>;
using LateTransformPermutations = std::array<
    std::array<std::size_t, late_reverb_coefficients::LineCount>,
    late_reverb_coefficients::VelvetStageCount + 1>;
using LateTransformSigns = std::array<
    std::array<float, late_reverb_coefficients::LineCount>,
    late_reverb_coefficients::VelvetStageCount + 1>;

inline const LateMainDelayRatios &
LateMainRatios(const LateReverbFlavour flavour) noexcept {
  return flavour == LateReverbFlavour::Optimized
             ? late_reverb_optimized_coefficients::MainDelayRatio
             : late_reverb_coefficients::MainDelayRatio;
}

inline const LateVelvetDelays &
LateVelvetDelayMilliseconds(const LateReverbFlavour flavour) noexcept {
  return flavour == LateReverbFlavour::Optimized
             ? late_reverb_optimized_coefficients::VelvetDelayMs
             : late_reverb_coefficients::VelvetDelayMs;
}

inline const LateTransformPermutations &
LatePermutations(const LateReverbFlavour flavour) noexcept {
  return flavour == LateReverbFlavour::Optimized
             ? late_reverb_optimized_coefficients::TransformPermutation
             : late_reverb_coefficients::TransformPermutation;
}

inline const LateTransformSigns &
LateSigns(const LateReverbFlavour flavour) noexcept {
  return flavour == LateReverbFlavour::Optimized
             ? late_reverb_optimized_coefficients::TransformSign
             : late_reverb_coefficients::TransformSign;
}

} // namespace tfdsp
