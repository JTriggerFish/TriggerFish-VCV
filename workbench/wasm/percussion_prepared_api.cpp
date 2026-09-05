#include "percussion_api.hpp"

#include "percussion_prepared.hpp"
#include "percussion_session.hpp"

#include <cmath>
#include <cstring>

using tfworkbench::detail::Find;
using tfworkbench::detail::Recipe;

extern "C" {

std::uint32_t tf_percussion_prepared_size(
    const std::uint32_t handle) noexcept {
  const auto *session = Find(handle);
  if (!session) return 0;
  switch (session->recipe) {
  case Recipe::MetallicPlate:
    return sizeof(tfworkbench::PreparedMetallicRecipe);
  case Recipe::Kick: return sizeof(tfworkbench::PreparedKickRecipe);
  case Recipe::MembraneDrum:
    return sizeof(tfworkbench::PreparedMembraneRecipe);
  case Recipe::SnareDrum: return sizeof(tfworkbench::PreparedSnareRecipe);
  default: return 0;
  }
}

int tf_percussion_export_prepared(
    const std::uint32_t handle, void *destination,
    const std::uint32_t byteSize) noexcept {
  const auto *session = Find(handle);
  if (!session || !destination ||
      byteSize != tf_percussion_prepared_size(handle)) return 0;
  try {
    if (session->recipe == Recipe::MetallicPlate) {
      tfworkbench::PreparedMetallicRecipe prepared;
      prepared.header.recipe = static_cast<std::uint32_t>(session->recipe);
      prepared.header.byteSize = sizeof(prepared);
      prepared.header.sampleRate = session->sampleRate;
      auto parameters = tfdsp::percussion::DefaultCrashCymbalParameters(
          session->sampleRate, tfworkbench::ApplyCrashMacros(
              session->crashBase, session->crashValues));
      parameters.routing = session->cymbalRouting;
      prepared.parameters = tfdsp::percussion::PrepareCrashCymbalParameters(
          session->sampleRate, parameters);
      std::memcpy(destination, &prepared, sizeof(prepared));
    } else if (session->recipe == Recipe::Kick) {
      tfworkbench::PreparedKickRecipe prepared;
      prepared.header.recipe = static_cast<std::uint32_t>(session->recipe);
      prepared.header.byteSize = sizeof(prepared);
      prepared.header.sampleRate = session->sampleRate;
      auto parameters = tfworkbench::ApplyKickParameters(session->kickValues);
      tfdsp::percussion::ApplyKickRouting(parameters, session->kickRouting);
      prepared.parameters = tfdsp::percussion::PrepareMembraneDrumParameters(
          session->sampleRate, parameters);
      std::memcpy(destination, &prepared, sizeof(prepared));
    } else if (session->recipe == Recipe::MembraneDrum) {
      tfworkbench::PreparedMembraneRecipe prepared;
      prepared.header.recipe = static_cast<std::uint32_t>(session->recipe);
      prepared.header.byteSize = sizeof(prepared);
      prepared.header.sampleRate = session->sampleRate;
      auto parameters = tfworkbench::ApplyMembraneParameters(
          session->membraneValues);
      parameters.routing = session->membraneRouting;
      prepared.parameters = tfdsp::percussion::PrepareMembraneDrumParameters(
          session->sampleRate, parameters);
      std::memcpy(destination, &prepared, sizeof(prepared));
    } else {
      tfworkbench::PreparedSnareRecipe prepared;
      prepared.header.recipe = static_cast<std::uint32_t>(session->recipe);
      prepared.header.byteSize = sizeof(prepared);
      prepared.header.sampleRate = session->sampleRate;
      auto parameters = tfworkbench::ApplySnareParameters(
          session->snareValues);
      parameters.routing = session->snareRouting;
      prepared.parameters = tfdsp::percussion::PrepareSnareDrumParameters(
          session->sampleRate, parameters);
      std::memcpy(destination, &prepared, sizeof(prepared));
    }
    return 1;
  } catch (...) {
    return 0;
  }
}

int tf_percussion_apply_prepared(
    const std::uint32_t handle, const void *source,
    const std::uint32_t byteSize) noexcept {
  auto *session = Find(handle);
  if (!session || !source ||
      byteSize < sizeof(tfworkbench::PreparedRecipeHeader)) return 0;
  tfworkbench::PreparedRecipeHeader header;
  std::memcpy(&header, source, sizeof(header));
  if (header.magic != tfworkbench::PreparedRecipeMagic ||
      header.version != tfworkbench::PreparedRecipeVersion ||
      header.recipe != static_cast<std::uint32_t>(session->recipe) ||
      header.byteSize != byteSize ||
      !std::isfinite(header.sampleRate) ||
      std::abs(header.sampleRate - session->sampleRate) > .01f) return 0;
  try {
    if (session->recipe == Recipe::MetallicPlate &&
        byteSize == sizeof(tfworkbench::PreparedMetallicRecipe)) {
      tfworkbench::PreparedMetallicRecipe prepared;
      std::memcpy(&prepared, source, sizeof(prepared));
      if (!std::isfinite(prepared.parameters.sampleRate) ||
          !std::isfinite(prepared.parameters.modalField.sampleRate) ||
          std::abs(prepared.parameters.sampleRate - session->sampleRate) >
              .01f ||
          std::abs(prepared.parameters.modalField.sampleRate -
                   session->sampleRate) > .01f) return 0;
      session->cymbal.Prepare(prepared.parameters);
      return 1;
    }
    if (session->recipe == Recipe::Kick &&
        byteSize == sizeof(tfworkbench::PreparedKickRecipe)) {
      tfworkbench::PreparedKickRecipe prepared;
      std::memcpy(&prepared, source, sizeof(prepared));
      if (!std::isfinite(prepared.parameters.sampleRate) ||
          !std::isfinite(prepared.parameters.membrane.sampleRate) ||
          std::abs(prepared.parameters.sampleRate - session->sampleRate) > .01f ||
          std::abs(prepared.parameters.membrane.sampleRate - session->sampleRate) > .01f) return 0;
      session->kick.Prepare(prepared.parameters);
      return 1;
    }
    if (session->recipe == Recipe::MembraneDrum &&
        byteSize == sizeof(tfworkbench::PreparedMembraneRecipe)) {
      tfworkbench::PreparedMembraneRecipe prepared;
      std::memcpy(&prepared, source, sizeof(prepared));
      if (!std::isfinite(prepared.parameters.sampleRate) ||
          !std::isfinite(prepared.parameters.membrane.sampleRate) ||
          std::abs(prepared.parameters.sampleRate - session->sampleRate) >
              .01f ||
          std::abs(prepared.parameters.membrane.sampleRate -
                   session->sampleRate) > .01f) return 0;
      session->membrane.Prepare(prepared.parameters);
      return 1;
    }
    if (session->recipe == Recipe::SnareDrum &&
        byteSize == sizeof(tfworkbench::PreparedSnareRecipe)) {
      tfworkbench::PreparedSnareRecipe prepared;
      std::memcpy(&prepared, source, sizeof(prepared));
      if (!std::isfinite(prepared.parameters.sampleRate) ||
          !std::isfinite(prepared.parameters.membrane.sampleRate) ||
          !std::isfinite(prepared.parameters.wires.sampleRate) ||
          std::abs(prepared.parameters.sampleRate - session->sampleRate) > .01f ||
          std::abs(prepared.parameters.membrane.sampleRate -
                   session->sampleRate) > .01f ||
          std::abs(prepared.parameters.wires.sampleRate -
                   session->sampleRate) > .01f) return 0;
      session->snare.Prepare(prepared.parameters);
      return 1;
    }
  } catch (...) {
    return 0;
  }
  return 0;
}

} // extern "C"
