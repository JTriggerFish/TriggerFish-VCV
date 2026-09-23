#pragma once

#include "early_reflection_path_renderer.hpp"
#include "finite_audio.hpp"
#include "fixed_radix2_fft.hpp"

namespace tfdsp {
template <typename Sample = float, std::size_t PartitionSize = 128,
          std::size_t MaximumSources = EarlyReflectionMaximumSources>
class EarlyReflectionConvolver {
#include "early_reflection_convolver/core.inl"
#include "early_reflection_convolver/public_api.inl"
};

} // namespace tfdsp
