#pragma once

#include "breakpoint_trajectory.hpp"

namespace tfdsp::percussion {

// Coherent FM audio -> modal impulse increments. Normalize envelope area
// relative to its peak, preserving the user's amplitude trajectory level.
// O(segment count), with no calibration sample rate or fitted gain constant.
template <typename Trajectory>
float TrajectoryImpulseScale(const Trajectory &trajectory,
                             const float sampleRate) noexcept {
  double previous = std::abs(tfdsp::FiniteNormalOrZero(trajectory.initialValue));
  double peak = previous;
  double area = 0.0;
  const auto count = std::min(trajectory.segmentCount, trajectory.segments.size());
  for (std::size_t index = 0; index < count; ++index) {
    const auto &segment = trajectory.segments[index];
    const double next = std::abs(tfdsp::FiniteNormalOrZero(segment.targetValue));
    const double duration = std::clamp(
        tfdsp::FiniteNormalOrZero(segment.durationSeconds), 0.f, 60.f);
    double mean = .5 * (previous + next);
    if (segment.curve == TrajectoryCurve::Geometric && previous > 0 &&
        next > 0 && std::abs(next - previous) > 1.e-12)
      mean = (next - previous) / std::log(next / previous);
    area += duration * mean;
    peak = std::max(peak, next);
    previous = next;
  }
  return peak > 0 ? static_cast<float>(
      peak / std::max(area * sampleRate, peak)) : 0.f;
}

} // namespace tfdsp::percussion
