#pragma once
#include "parameter_descriptor.hpp"
#include <cstddef>

namespace tfworkbench {
inline constexpr std::size_t KickModeParameterCount = 16 * 4;
const ParameterDescriptor &KickModeDescription(std::size_t index) noexcept;
}
