#include "crash_api.hpp"

#include "crash_macros.hpp"
#include "percussion_api.hpp"
#include "tfdsp/percussion/metallic_plate_routing.hpp"

extern "C" {

std::uint32_t tf_crash_api_version() noexcept {
  return tf_percussion_api_version();
}

std::uint32_t tf_crash_create(const float sampleRate) noexcept {
  return tf_percussion_create(0, sampleRate);
}

void tf_crash_destroy(const std::uint32_t handle) noexcept {
  tf_percussion_destroy(handle);
}

int tf_crash_reset(const std::uint32_t handle) noexcept {
  return tf_percussion_reset(handle);
}

int tf_crash_trigger(
    const std::uint32_t handle, const float strength, const float location,
    const float hardness, const float implement, const float contactSpread,
    const std::uint32_t seed) noexcept {
  return tf_percussion_trigger(handle, strength, location, hardness, implement,
                               contactSpread, seed);
}

int tf_crash_set_mute(const std::uint32_t handle,
                      const float amount) noexcept {
  return tf_percussion_set_mute(handle, amount);
}

int tf_crash_process(const std::uint32_t handle, float *output,
                     const std::uint32_t frames) noexcept {
  return tf_percussion_process(handle, output, frames);
}

std::uint32_t tf_crash_macro_count() noexcept {
  return static_cast<std::uint32_t>(tfworkbench::CrashMacroCount);
}

const char *tf_crash_macro_key(const std::uint32_t index) noexcept {
  return index < tfworkbench::CrashMacroCount
      ? tfworkbench::CrashMacroDescription(index).key.c_str() : nullptr;
}

const char *tf_crash_macro_name(const std::uint32_t index) noexcept {
  return index < tfworkbench::CrashMacroCount
      ? tfworkbench::CrashMacroDescription(index).name.c_str() : nullptr;
}

const char *tf_crash_macro_unit(const std::uint32_t index) noexcept {
  return index < tfworkbench::CrashMacroCount
      ? tfworkbench::CrashMacroDescription(index).unit.c_str() : nullptr;
}

int tf_crash_macro_scale(const std::uint32_t index) noexcept {
  return index < tfworkbench::CrashMacroCount
      ? static_cast<int>(tfworkbench::CrashMacroDescription(index).scale) : -1;
}

float tf_crash_macro_minimum(const std::uint32_t index) noexcept {
  return index < tfworkbench::CrashMacroCount
      ? tfworkbench::CrashMacroDescription(index).minimum : 0.f;
}

float tf_crash_macro_maximum(const std::uint32_t index) noexcept {
  return index < tfworkbench::CrashMacroCount
      ? tfworkbench::CrashMacroDescription(index).maximum : 0.f;
}

float tf_crash_macro_default(const std::uint32_t index) noexcept {
  return index < tfworkbench::CrashMacroCount
      ? tfworkbench::CrashMacroDescription(index).defaultValue : 0.f;
}

float tf_crash_macro_get(const std::uint32_t handle,
                         const std::uint32_t index) noexcept {
  return tfworkbench::detail::LegacyCrashParameterGet(handle, index);
}

int tf_crash_macro_set(const std::uint32_t handle, const std::uint32_t index,
                       const float value) noexcept {
  return tfworkbench::detail::LegacyCrashParameterSet(handle, index, value);
}

int tf_crash_macro_commit(const std::uint32_t handle) noexcept {
  return tf_percussion_commit(handle);
}

std::uint32_t tf_crash_route_count() noexcept {
  return static_cast<std::uint32_t>(
      tfdsp::percussion::MetallicPlateRouting::Count);
}

int tf_crash_route_enabled(const std::uint32_t handle,
                           const std::uint32_t index) noexcept {
  return tf_percussion_route_enabled(handle, index);
}

int tf_crash_route_enable(const std::uint32_t handle,
                          const std::uint32_t index,
                          const int enabled) noexcept {
  return tf_percussion_route_enable(handle, index, enabled);
}

} // extern "C"
