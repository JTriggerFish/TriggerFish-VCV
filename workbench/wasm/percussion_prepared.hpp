#pragma once

#include "tfdsp/percussion/compact_kick_parameters.hpp"
#include "tfdsp/percussion/crash_cymbal_parameters.hpp"
#include "tfdsp/percussion/membrane_drum_parameters.hpp"

#include <cstdint>
#include <type_traits>

namespace tfworkbench {

inline constexpr std::uint32_t PreparedRecipeMagic = 0x54465052u;
inline constexpr std::uint32_t PreparedRecipeVersion = 1;

struct PreparedRecipeHeader {
  std::uint32_t magic{PreparedRecipeMagic};
  std::uint32_t version{PreparedRecipeVersion};
  std::uint32_t recipe{};
  std::uint32_t byteSize{};
  float sampleRate{48000.f};
};

struct PreparedMetallicRecipe {
  PreparedRecipeHeader header{};
  tfdsp::percussion::CrashCymbalPreparedParameters parameters{};
};

struct PreparedKickRecipe {
  PreparedRecipeHeader header{};
  tfdsp::percussion::CompactKickParameters parameters{};
};

struct PreparedMembraneRecipe {
  PreparedRecipeHeader header{};
  tfdsp::percussion::MembraneDrumPreparedParameters parameters{};
};

static_assert(std::is_trivially_copyable_v<PreparedMetallicRecipe>);
static_assert(std::is_trivially_copyable_v<PreparedKickRecipe>);
static_assert(std::is_trivially_copyable_v<PreparedMembraneRecipe>);

} // namespace tfworkbench
