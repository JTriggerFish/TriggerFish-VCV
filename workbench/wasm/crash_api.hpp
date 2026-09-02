#pragma once

#include <cstdint>

// Stable C boundary shared by the offline Wasm renderer and the future audio
// worklet. Handles are local to one Wasm instance and calls are single-threaded.
extern "C" {

std::uint32_t tf_crash_api_version() noexcept;
std::uint32_t tf_crash_create(float sampleRate) noexcept;
void tf_crash_destroy(std::uint32_t handle) noexcept;
int tf_crash_reset(std::uint32_t handle) noexcept;
int tf_crash_trigger(std::uint32_t handle, float strength, float location,
                     float hardness, std::uint32_t seed) noexcept;
int tf_crash_set_mute(std::uint32_t handle, float amount) noexcept;
int tf_crash_process(std::uint32_t handle, float *output,
                     std::uint32_t frames) noexcept;
std::uint32_t tf_crash_macro_count() noexcept;
const char *tf_crash_macro_name(std::uint32_t index) noexcept;
const char *tf_crash_macro_unit(std::uint32_t index) noexcept;
float tf_crash_macro_minimum(std::uint32_t index) noexcept;
float tf_crash_macro_maximum(std::uint32_t index) noexcept;
float tf_crash_macro_default(std::uint32_t index) noexcept;
float tf_crash_macro_get(std::uint32_t handle, std::uint32_t index) noexcept;
int tf_crash_macro_set(std::uint32_t handle, std::uint32_t index,
                       float value) noexcept;
int tf_crash_macro_commit(std::uint32_t handle) noexcept;

}
