#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace tfui {

enum class CursorTravelCurve { Linear, Smoothstep };

inline std::int64_t arrangementCursorGroup(double relativeClockBeat,
                                           int clocksPerPulse) noexcept {
  const double group =
      std::floor(std::max(0.0, relativeClockBeat) /
                 static_cast<double>(std::max(1, clocksPerPulse)));
  if (group >= static_cast<double>(std::numeric_limits<std::int64_t>::max()))
    return std::numeric_limits<std::int64_t>::max();
  return static_cast<std::int64_t>(group);
}

constexpr double CursorHeadSeconds = 0.055;
constexpr double CursorTailSeconds = 0.32;
constexpr double ExecutionBloomExpansionSeconds = 0.18;

inline float exponentialDecay(double age, double timeConstant) noexcept {
  if (!(timeConstant > 0.0))
    return 0.f;
  return static_cast<float>(std::exp(-std::max(0.0, age) / timeConstant));
}

inline double cursorTravelDuration(double pulseInterval) noexcept {
  // Keep the beam visible for several display frames, while still scaling it
  // with the musical interval so fast patterns do not turn into a solid bar.
  return std::clamp(0.75 * std::max(0.0, pulseInterval), 0.035, 0.12);
}

inline double cursorMotionTailDuration(double pulseInterval) noexcept {
  return std::clamp(1.1 * std::max(0.0, pulseInterval), 0.075, 0.30);
}

inline double cursorBloomExpansionDuration(double pulseInterval) noexcept {
  return std::clamp(0.85 * std::max(0.0, pulseInterval), 0.070, 0.20);
}

inline double cursorBloomTailDuration(double pulseInterval) noexcept {
  return std::clamp(1.4 * std::max(0.0, pulseInterval), 0.15, 0.42);
}

inline float cursorBloomExpansion(double age, double duration) noexcept {
  if (!(duration > 0.0))
    return 1.f;
  const float progress =
      static_cast<float>(std::clamp(age / duration, 0.0, 1.0));
  return progress * progress * (3.f - 2.f * progress);
}

inline float cursorTravelPosition(CursorTravelCurve curve,
                                  float progress) noexcept {
  progress = std::clamp(progress, 0.f, 1.f);
  if (curve == CursorTravelCurve::Smoothstep)
    return progress * progress * (3.f - 2.f * progress);
  return progress;
}

inline float cursorMotionIntensity(double age, double duration,
                                   double tailDuration) noexcept {
  if (!(duration > 0.0) || !(tailDuration > 0.0))
    return 0.f;
  const float progress =
      static_cast<float>(std::clamp(age / duration, 0.0, 1.0));
  if (progress < 1.f)
    return 0.65f + 0.35f * std::sin(3.14159265358979323846f * progress);
  return 0.65f * exponentialDecay(age - duration, tailDuration);
}

inline float accumulatedTail(float previousEnergy, double age, float impulse,
                             double tailDuration = CursorTailSeconds) noexcept {
  return std::clamp(
      previousEnergy * exponentialDecay(age, tailDuration) + impulse, 0.f, 1.f);
}

} // namespace tfui
