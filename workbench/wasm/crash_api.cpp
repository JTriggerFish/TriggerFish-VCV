#include "crash_api.hpp"
#include "crash_macros.hpp"

#include "tfdsp/percussion/crash_cymbal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace {

constexpr std::size_t MaximumSessions = 4;

struct Session {
  tfdsp::percussion::CrashCymbal cymbal{};
  tfdsp::percussion::CrashCymbalFitParameters baseFit{};
  tfworkbench::CrashMacroValues macros{};
  float sampleRate{48000.f};
  bool active{};
};

std::array<Session, MaximumSessions> sessions{};

Session *Find(const std::uint32_t handle) noexcept {
  if (handle == 0 || handle > sessions.size())
    return nullptr;
  auto &session = sessions[handle - 1];
  return session.active ? &session : nullptr;
}

} // namespace

extern "C" {

std::uint32_t tf_crash_api_version() noexcept {
  return 12;
}

std::uint32_t tf_crash_create(const float sampleRate) noexcept {
  if (!std::isfinite(sampleRate) || sampleRate < 8000.f || sampleRate > 384000.f)
    return 0;
  for (std::size_t index = 0; index < sessions.size(); ++index) {
    auto &session = sessions[index];
    if (session.active)
      continue;
    session.baseFit = tfworkbench::CrashWorkbenchBaseFit();
    session.macros = tfworkbench::DefaultCrashMacros();
    session.sampleRate = sampleRate;
    session.cymbal.Prepare(
        sampleRate, tfdsp::percussion::DefaultCrashCymbalParameters(
            sampleRate, tfworkbench::ApplyCrashMacros(
                session.baseFit, session.macros)));
    session.active = true;
    return static_cast<std::uint32_t>(index + 1);
  }
  return 0;
}

void tf_crash_destroy(const std::uint32_t handle) noexcept {
  if (auto *session = Find(handle)) {
    session->cymbal.Reset();
    session->active = false;
  }
}

int tf_crash_reset(const std::uint32_t handle) noexcept {
  auto *session = Find(handle);
  if (!session)
    return 0;
  session->cymbal.Reset();
  return 1;
}

int tf_crash_trigger(const std::uint32_t handle, const float strength,
                     const float location, const float hardness,
                     const float implement, const float contactSpread,
                     const std::uint32_t seed) noexcept {
  auto *session = Find(handle);
  if (!session)
    return 0;
  session->cymbal.Trigger(
      {strength, location, hardness, seed, implement, contactSpread});
  return 1;
}

int tf_crash_set_mute(const std::uint32_t handle, const float amount) noexcept {
  auto *session = Find(handle);
  if (!session)
    return 0;
  session->cymbal.SetMute(amount);
  return 1;
}

int tf_crash_process(const std::uint32_t handle, float *output,
                     const std::uint32_t frames) noexcept {
  auto *session = Find(handle);
  if (!session || (!output && frames != 0))
    return 0;
  for (std::uint32_t frame = 0; frame < frames; ++frame)
    output[frame] = session->cymbal.Process();
  return 1;
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
  auto *session = Find(handle);
  return session && index < session->macros.size()
      ? session->macros[index] : 0.f;
}

int tf_crash_macro_set(const std::uint32_t handle, const std::uint32_t index,
                       const float value) noexcept {
  auto *session = Find(handle);
  if (!session || index >= session->macros.size() || !std::isfinite(value))
    return 0;
  const auto &descriptor = tfworkbench::CrashMacroDescription(index);
  session->macros[index] = std::clamp(value, descriptor.minimum,
                                     descriptor.maximum);
  return 1;
}

int tf_crash_macro_commit(const std::uint32_t handle) noexcept {
  auto *session = Find(handle);
  if (!session)
    return 0;
  const auto fit = tfworkbench::ApplyCrashMacros(
      session->baseFit, session->macros);
  session->cymbal.Prepare(
      session->sampleRate,
      tfdsp::percussion::DefaultCrashCymbalParameters(session->sampleRate, fit));
  return 1;
}

} // extern "C"
