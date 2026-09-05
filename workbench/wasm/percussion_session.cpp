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
  session.crashBase = MetallicWorkbenchBaseFit();
  session.crashValues = DefaultCrashMacros();
  session.kickValues = DefaultKickParameters();
  session.membraneValues = DefaultMembraneParameters();
  session.snareValues = DefaultSnareParameters();
  session.cymbalRouting = {};
  session.kickRouting = {};
  session.membraneRouting = {};
  session.snareRouting = {};
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
  switch (session.recipe) {
  case Recipe::MetallicPlate:
    return index < ActiveCrashMacroCount
        ? &ActiveCrashMacroDescription(index) : nullptr;
  case Recipe::Kick:
    return index < session.kickValues.size()
        ? &KickParameterDescription(index) : nullptr;
  case Recipe::MembraneDrum:
    return index < session.membraneValues.size()
        ? &MembraneParameterDescription(index) : nullptr;
  case Recipe::SnareDrum:
    return index < session.snareValues.size()
        ? &SnareParameterDescription(index) : nullptr;
  default: return nullptr;
  }
}

std::size_t ParameterCount(const Session &session) noexcept {
  switch (session.recipe) {
  case Recipe::MetallicPlate: return ActiveCrashMacroCount;
  case Recipe::Kick: return session.kickValues.size();
  case Recipe::MembraneDrum: return session.membraneValues.size();
  case Recipe::SnareDrum: return session.snareValues.size();
  default: return 0;
  }
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
  if (session.recipe == Recipe::Kick) {
    auto parameters = ApplyKickParameters(session.kickValues);
    tfdsp::percussion::ApplyKickRouting(parameters, session.kickRouting);
    session.kick.Prepare(session.sampleRate, parameters);
    return;
  }
  if (session.recipe == Recipe::MembraneDrum) {
    auto parameters = ApplyMembraneParameters(session.membraneValues);
    parameters.routing = session.membraneRouting;
    session.membrane.Prepare(session.sampleRate, parameters);
    return;
  }
  auto parameters = ApplySnareParameters(session.snareValues);
  parameters.routing = session.snareRouting;
  session.snare.Prepare(session.sampleRate, parameters);
}

float Process(Session &session) noexcept {
  switch (session.recipe) {
  case Recipe::MetallicPlate: return session.cymbal.Process();
  case Recipe::Kick: return session.kick.Process();
  case Recipe::MembraneDrum: return session.membrane.Process();
  case Recipe::SnareDrum: return session.snare.Process();
  default: return 0.f;
  }
}

} // namespace tfworkbench::detail
