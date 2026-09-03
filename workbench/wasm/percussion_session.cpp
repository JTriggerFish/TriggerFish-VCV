#include "percussion_session.hpp"

#include <array>
#include <cmath>

namespace tfworkbench::detail {
namespace {

constexpr std::size_t MaximumSessions = 4;
std::array<Session, MaximumSessions> sessions{};

bool ValidRecipeAndRate(const std::uint32_t recipe,
                        const float sampleRate) noexcept {
  return recipe < static_cast<std::uint32_t>(Recipe::Count) &&
      std::isfinite(sampleRate) && sampleRate >= 8000.f &&
      sampleRate <= 384000.f;
}

void Initialize(Session &session, const Recipe recipe,
                const float sampleRate) noexcept {
  session.recipe = recipe;
  session.sampleRate = sampleRate;
  session.crashBase = CrashWorkbenchBaseFit();
  session.crashValues = DefaultCrashMacros();
  session.kickValues = DefaultKickParameters();
  session.cymbalRouting = {};
  session.kickRouting = {};
}

} // namespace

Session *Find(const std::uint32_t handle) noexcept {
  if (handle == 0 || handle > sessions.size()) return nullptr;
  auto &session = sessions[handle - 1];
  return session.active ? &session : nullptr;
}

std::uint32_t Allocate(const std::uint32_t recipe,
                       const float sampleRate) noexcept {
  if (!ValidRecipeAndRate(recipe, sampleRate)) return 0;
  for (std::size_t index = 0; index < sessions.size(); ++index) {
    auto &session = sessions[index];
    if (session.active) continue;
    Initialize(session, static_cast<Recipe>(recipe), sampleRate);
    session.active = true;
    return static_cast<std::uint32_t>(index + 1);
  }
  return 0;
}

const ParameterDescriptor *Description(
    const Session &session, const std::size_t index) noexcept {
  if (session.recipe == Recipe::MetallicPlate) {
    return index < ActiveCrashMacroCount
        ? &ActiveCrashMacroDescription(index) : nullptr;
  }
  return index < session.kickValues.size()
      ? &KickParameterDescription(index) : nullptr;
}

std::size_t ParameterCount(const Session &session) noexcept {
  return session.recipe == Recipe::MetallicPlate
      ? ActiveCrashMacroCount : session.kickValues.size();
}

void Prepare(Session &session) {
  if (session.recipe == Recipe::MetallicPlate) {
    auto parameters = tfdsp::percussion::DefaultCrashCymbalParameters(
        session.sampleRate, ApplyCrashMacros(
            session.crashBase, session.crashValues));
    parameters.routing = session.cymbalRouting;
    session.cymbal.Prepare(session.sampleRate, parameters);
    return;
  }
  auto parameters = ApplyKickParameters(session.kickValues);
  parameters.routing = session.kickRouting;
  session.kick.Prepare(session.sampleRate, parameters);
}

float Process(Session &session) noexcept {
  return session.recipe == Recipe::MetallicPlate
      ? session.cymbal.Process() : session.kick.Process();
}

} // namespace tfworkbench::detail
