#pragma once

#include "crash_macros.hpp"
#include "kick_macros.hpp"
#include "tfdsp/percussion/compact_kick.hpp"
#include "tfdsp/percussion/crash_cymbal.hpp"

#include <cstddef>
#include <cstdint>

namespace tfworkbench::detail {

enum class Recipe : std::uint32_t { MetallicPlate, CompactKick, Count };

struct Session {
  tfdsp::percussion::CrashCymbal cymbal{};
  tfdsp::percussion::CompactKick kick{};
  tfdsp::percussion::CrashCymbalFitParameters crashBase{};
  CrashMacroValues crashValues{};
  KickParameterValues kickValues{};
  tfdsp::percussion::MetallicPlateRouting cymbalRouting{};
  tfdsp::percussion::CompactKickRouting kickRouting{};
  float sampleRate{48000.f};
  Recipe recipe{Recipe::MetallicPlate};
  bool active{};
};

Session *Find(std::uint32_t handle) noexcept;
std::uint32_t Allocate(std::uint32_t recipe, float sampleRate) noexcept;
void Prepare(Session &session);
float Process(Session &session) noexcept;

const ParameterDescriptor *Description(
    const Session &session, std::size_t index) noexcept;
std::size_t ParameterCount(const Session &session) noexcept;

} // namespace tfworkbench::detail
