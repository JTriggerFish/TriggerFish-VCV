#include "percussion_api.hpp"

#include "crash_macros.hpp"
#include "percussion_session.hpp"

#include <algorithm>
#include <array>
#include <cmath>

using tfworkbench::detail::Allocate;
using tfworkbench::detail::Description;
using tfworkbench::detail::Find;
using tfworkbench::detail::ParameterCount;
using tfworkbench::detail::Prepare;
using tfworkbench::detail::Process;
using tfworkbench::detail::Recipe;

extern "C" {

std::uint32_t tf_percussion_api_version() noexcept { return 1; }
std::uint32_t tf_percussion_recipe_count() noexcept {
  return static_cast<std::uint32_t>(Recipe::Count);
}

const char *tf_percussion_recipe_key(const std::uint32_t recipe) noexcept {
  static constexpr std::array keys{
      "metal.cymbal.v1", "drum.kick-fm.v1", "drum.membrane.v1",
      "drum.snare.v1"};
  return recipe < keys.size() ? keys[recipe] : nullptr;
}

const char *tf_percussion_recipe_name(const std::uint32_t recipe) noexcept {
  static constexpr std::array names{
      "Metallic plate", "Compact FM kick", "Membrane drum", "Snare drum"};
  return recipe < names.size() ? names[recipe] : nullptr;
}

std::uint32_t tf_percussion_create(const std::uint32_t recipe,
                                   const float sampleRate) noexcept {
  const auto handle = tf_percussion_create_unprepared(recipe, sampleRate);
  auto *session = Find(handle);
  if (!session) return 0;
  try {
    Prepare(*session);
    return handle;
  } catch (...) {
    tf_percussion_destroy(handle);
    return 0;
  }
}

std::uint32_t tf_percussion_create_unprepared(
    const std::uint32_t recipe, const float sampleRate) noexcept {
  return Allocate(recipe, sampleRate);
}

void tf_percussion_destroy(const std::uint32_t handle) noexcept {
  if (auto *session = Find(handle)) {
    switch (session->recipe) {
    case Recipe::MetallicPlate: session->cymbal.Reset(); break;
    case Recipe::CompactKick: session->kick.Reset(); break;
    case Recipe::MembraneDrum: session->membrane.Reset(); break;
    case Recipe::SnareDrum: session->snare.Reset(); break;
    default: break;
    }
    session->active = false;
  }
}

std::uint32_t tf_percussion_recipe(const std::uint32_t handle) noexcept {
  const auto *session = Find(handle);
  return session ? static_cast<std::uint32_t>(session->recipe) : UINT32_MAX;
}

int tf_percussion_reset(const std::uint32_t handle) noexcept {
  auto *session = Find(handle);
  if (!session) return 0;
  switch (session->recipe) {
  case Recipe::MetallicPlate: session->cymbal.Reset(); break;
  case Recipe::CompactKick: session->kick.Reset(); break;
  case Recipe::MembraneDrum: session->membrane.Reset(); break;
  case Recipe::SnareDrum: session->snare.Reset(); break;
  default: return 0;
  }
  return 1;
}

int tf_percussion_trigger(
    const std::uint32_t handle, const float strength, const float location,
    const float hardness, const float implement, const float contactSpread,
    const std::uint32_t seed) noexcept {
  auto *session = Find(handle);
  if (!session) return 0;
  switch (session->recipe) {
  case Recipe::MetallicPlate:
    session->cymbal.Trigger(
        {strength, location, hardness, seed, implement, contactSpread});
    break;
  case Recipe::CompactKick:
    session->kick.Trigger({strength, hardness, seed});
    break;
  case Recipe::MembraneDrum:
    session->membrane.Trigger(
        {strength, location, hardness, implement, contactSpread, seed});
    break;
  case Recipe::SnareDrum:
    session->snare.Trigger(
        {strength, location, hardness, implement, contactSpread, seed});
    break;
  default: return 0;
  }
  return 1;
}

int tf_percussion_set_mute(const std::uint32_t handle,
                           const float amount) noexcept {
  auto *session = Find(handle);
  if (!session || !std::isfinite(amount)) return 0;
  if (session->recipe == Recipe::MetallicPlate)
    session->cymbal.SetMute(amount);
  return 1;
}

int tf_percussion_process(const std::uint32_t handle, float *output,
                          const std::uint32_t frames) noexcept {
  auto *session = Find(handle);
  if (!session || (!output && frames != 0)) return 0;
  for (std::uint32_t frame = 0; frame < frames; ++frame)
    output[frame] = Process(*session);
  return 1;
}

std::uint32_t tf_percussion_parameter_count(
    const std::uint32_t handle) noexcept {
  const auto *session = Find(handle);
  return session ? static_cast<std::uint32_t>(ParameterCount(*session)) : 0;
}

const char *tf_percussion_parameter_key(
    const std::uint32_t handle, const std::uint32_t index) noexcept {
  const auto *session = Find(handle);
  const auto *descriptor = session ? Description(*session, index) : nullptr;
  return descriptor ? descriptor->key.c_str() : nullptr;
}

const char *tf_percussion_parameter_name(
    const std::uint32_t handle, const std::uint32_t index) noexcept {
  const auto *session = Find(handle);
  const auto *descriptor = session ? Description(*session, index) : nullptr;
  return descriptor ? descriptor->name.c_str() : nullptr;
}

const char *tf_percussion_parameter_unit(
    const std::uint32_t handle, const std::uint32_t index) noexcept {
  const auto *session = Find(handle);
  const auto *descriptor = session ? Description(*session, index) : nullptr;
  return descriptor ? descriptor->unit.c_str() : nullptr;
}

int tf_percussion_parameter_scale(const std::uint32_t handle,
                                  const std::uint32_t index) noexcept {
  const auto *session = Find(handle);
  const auto *descriptor = session ? Description(*session, index) : nullptr;
  return descriptor ? static_cast<int>(descriptor->scale) : -1;
}

float tf_percussion_parameter_minimum(
    const std::uint32_t handle, const std::uint32_t index) noexcept {
  const auto *session = Find(handle);
  const auto *descriptor = session ? Description(*session, index) : nullptr;
  return descriptor ? descriptor->minimum : 0.f;
}

float tf_percussion_parameter_maximum(
    const std::uint32_t handle, const std::uint32_t index) noexcept {
  const auto *session = Find(handle);
  const auto *descriptor = session ? Description(*session, index) : nullptr;
  return descriptor ? descriptor->maximum : 0.f;
}

float tf_percussion_parameter_default(
    const std::uint32_t handle, const std::uint32_t index) noexcept {
  const auto *session = Find(handle);
  const auto *descriptor = session ? Description(*session, index) : nullptr;
  return descriptor ? descriptor->defaultValue : 0.f;
}

float tf_percussion_parameter_get(const std::uint32_t handle,
                                  const std::uint32_t index) noexcept {
  const auto *session = Find(handle);
  if (!session || index >= ParameterCount(*session)) return 0.f;
  switch (session->recipe) {
  case Recipe::MetallicPlate:
    return session->crashValues[tfworkbench::ActiveCrashMacroIndices[index]];
  case Recipe::CompactKick: return session->kickValues[index];
  case Recipe::MembraneDrum: return session->membraneValues[index];
  case Recipe::SnareDrum: return session->snareValues[index];
  default: return 0.f;
  }
}

int tf_percussion_parameter_set(const std::uint32_t handle,
                                const std::uint32_t index,
                                const float value) noexcept {
  auto *session = Find(handle);
  const auto *descriptor = session ? Description(*session, index) : nullptr;
  if (!descriptor || !std::isfinite(value) ||
      descriptor->scale == tfworkbench::ParameterScale::Choice &&
          value != std::round(value)) return 0;
  const float bounded = std::clamp(
      value, descriptor->minimum, descriptor->maximum);
  switch (session->recipe) {
  case Recipe::MetallicPlate:
    session->crashValues[tfworkbench::ActiveCrashMacroIndices[index]] = bounded;
    break;
  case Recipe::CompactKick:
    session->kickValues[index] = bounded;
    break;
  case Recipe::MembraneDrum:
    session->membraneValues[index] = bounded;
    break;
  case Recipe::SnareDrum:
    session->snareValues[index] = bounded;
    break;
  default: return 0;
  }
  return 1;
}

int tf_percussion_commit(const std::uint32_t handle) noexcept {
  auto *session = Find(handle);
  if (!session) return 0;
  try {
    Prepare(*session);
    return 1;
  } catch (...) {
    return 0;
  }
}

std::uint32_t tf_percussion_route_count(
    const std::uint32_t handle) noexcept {
  const auto *session = Find(handle);
  if (!session) return 0;
  switch (session->recipe) {
  case Recipe::MetallicPlate:
    return static_cast<std::uint32_t>(session->cymbalRouting.enabled.size());
  case Recipe::CompactKick:
    return static_cast<std::uint32_t>(session->kickRouting.enabled.size());
  case Recipe::MembraneDrum:
    return static_cast<std::uint32_t>(session->membraneRouting.enabled.size());
  case Recipe::SnareDrum:
    return static_cast<std::uint32_t>(session->snareRouting.enabled.size());
  default: return 0;
  }
}

int tf_percussion_route_enabled(const std::uint32_t handle,
                                const std::uint32_t index) noexcept {
  const auto *session = Find(handle);
  if (!session || index >= tf_percussion_route_count(handle)) return 0;
  switch (session->recipe) {
  case Recipe::MetallicPlate: return session->cymbalRouting.enabled[index];
  case Recipe::CompactKick: return session->kickRouting.enabled[index];
  case Recipe::MembraneDrum: return session->membraneRouting.enabled[index];
  case Recipe::SnareDrum: return session->snareRouting.enabled[index];
  default: return 0;
  }
}

int tf_percussion_route_enable(const std::uint32_t handle,
                               const std::uint32_t index,
                               const int enabled) noexcept {
  auto *session = Find(handle);
  if (!session || index >= tf_percussion_route_count(handle) ||
      (enabled != 0 && enabled != 1)) return 0;
  switch (session->recipe) {
  case Recipe::MetallicPlate:
    session->cymbalRouting.SetEnabled(index, enabled != 0);
    break;
  case Recipe::CompactKick:
    session->kickRouting.SetEnabled(index, enabled != 0);
    break;
  case Recipe::MembraneDrum:
    session->membraneRouting.SetEnabled(index, enabled != 0);
    break;
  case Recipe::SnareDrum:
    session->snareRouting.SetEnabled(index, enabled != 0);
    break;
  default: return 0;
  }
  return 1;
}

} // extern "C"

namespace tfworkbench::detail {

float LegacyCrashParameterGet(const std::uint32_t handle,
                              const std::uint32_t index) noexcept {
  const auto *session = Find(handle);
  if (!session || session->recipe != Recipe::MetallicPlate ||
      index >= session->crashValues.size()) return 0.f;
  return session->crashValues[index];
}

int LegacyCrashParameterSet(const std::uint32_t handle,
                            const std::uint32_t index,
                            const float value) noexcept {
  auto *session = Find(handle);
  if (!session || session->recipe != Recipe::MetallicPlate ||
      index >= session->crashValues.size() || !std::isfinite(value)) return 0;
  const auto &descriptor = CrashMacroDescription(index);
  session->crashValues[index] = std::clamp(
      value, descriptor.minimum, descriptor.maximum);
  return 1;
}

} // namespace tfworkbench::detail
