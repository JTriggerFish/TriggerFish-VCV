#pragma once

#include <cstdint>

extern "C" {

std::uint32_t tf_percussion_api_version() noexcept;
std::uint32_t tf_percussion_recipe_count() noexcept;
const char *tf_percussion_recipe_key(std::uint32_t recipe) noexcept;
const char *tf_percussion_recipe_name(std::uint32_t recipe) noexcept;

std::uint32_t tf_percussion_create(std::uint32_t recipe,
                                   float sampleRate) noexcept;
std::uint32_t tf_percussion_create_unprepared(std::uint32_t recipe,
                                              float sampleRate) noexcept;
void tf_percussion_destroy(std::uint32_t handle) noexcept;
std::uint32_t tf_percussion_recipe(std::uint32_t handle) noexcept;
int tf_percussion_reset(std::uint32_t handle) noexcept;
int tf_percussion_trigger(std::uint32_t handle, float strength,
                          float location, float hardness, float implement,
                          float contactSpread, std::uint32_t seed) noexcept;
int tf_percussion_set_mute(std::uint32_t handle, float amount) noexcept;
int tf_percussion_process(std::uint32_t handle, float *output,
                          std::uint32_t frames) noexcept;

std::uint32_t tf_percussion_parameter_count(std::uint32_t handle) noexcept;
const char *tf_percussion_parameter_key(std::uint32_t handle,
                                        std::uint32_t index) noexcept;
const char *tf_percussion_parameter_name(std::uint32_t handle,
                                         std::uint32_t index) noexcept;
const char *tf_percussion_parameter_unit(std::uint32_t handle,
                                         std::uint32_t index) noexcept;
int tf_percussion_parameter_scale(std::uint32_t handle,
                                  std::uint32_t index) noexcept;
float tf_percussion_parameter_minimum(std::uint32_t handle,
                                      std::uint32_t index) noexcept;
float tf_percussion_parameter_maximum(std::uint32_t handle,
                                      std::uint32_t index) noexcept;
float tf_percussion_parameter_default(std::uint32_t handle,
                                      std::uint32_t index) noexcept;
float tf_percussion_parameter_get(std::uint32_t handle,
                                  std::uint32_t index) noexcept;
int tf_percussion_parameter_set(std::uint32_t handle, std::uint32_t index,
                                float value) noexcept;
int tf_percussion_commit(std::uint32_t handle) noexcept;
std::uint32_t tf_percussion_prepared_size(std::uint32_t handle) noexcept;
int tf_percussion_export_prepared(std::uint32_t handle, void *destination,
                                  std::uint32_t byteSize) noexcept;
int tf_percussion_apply_prepared(std::uint32_t handle, const void *source,
                                 std::uint32_t byteSize) noexcept;

std::uint32_t tf_percussion_route_count(std::uint32_t handle) noexcept;
float tf_percussion_route_get(std::uint32_t handle,
                              std::uint32_t index) noexcept;
int tf_percussion_route_set(std::uint32_t handle, std::uint32_t index,
                            float gain) noexcept;

}

namespace tfworkbench::detail {

float LegacyCrashParameterGet(std::uint32_t handle,
                              std::uint32_t index) noexcept;
int LegacyCrashParameterSet(std::uint32_t handle, std::uint32_t index,
                            float value) noexcept;

} // namespace tfworkbench::detail
