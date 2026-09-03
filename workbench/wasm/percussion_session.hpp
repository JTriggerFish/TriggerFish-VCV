#pragma once

#include "crash_macros.hpp"
#include "kick_macros.hpp"
#include "membrane_macros.hpp"
#include "tfdsp/percussion/compact_kick.hpp"
#include "tfdsp/percussion/crash_cymbal.hpp"
#include "tfdsp/percussion/membrane_drum.hpp"

#include <cstddef>
#include <cstdint>

namespace tfworkbench::detail {

enum class Recipe : std::uint32_t {
  MetallicPlate, CompactKick, MembraneDrum, Count
};

struct Session {
  tfdsp::percussion::CrashCymbal cymbal{};
  tfdsp::percussion::CompactKick kick{};
  tfdsp::percussion::MembraneDrum membrane{};
  tfdsp::percussion::CrashCymbalFitParameters crashBase{};
  CrashMacroValues crashValues{};
  KickParameterValues kickValues{};
  MembraneParameterValues membraneValues{};
  tfdsp::percussion::MetallicPlateRouting cymbalRouting{};
  tfdsp::percussion::CompactKickRouting kickRouting{};
  tfdsp::percussion::MembraneDrumRouting membraneRouting{};
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
