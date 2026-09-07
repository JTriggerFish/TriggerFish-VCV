#pragma once

#include "correlated_fm_burst.hpp"
#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {

// log amplitude = -ln(1000) * ((1-shape)*u + shape*u*u), u=t/T60.
// Shape 1 has a rounded full-level shoulder and an accelerating decay. This
// shapes a source, never clips its waveform or rescales stored modal energy.
inline void ShapeKickThump(CorrelatedFmTrajectory &amplitude,
                          float t60, float hold, float shape) noexcept {
  const auto safe = [](float value, float low, float high) {
    return std::clamp(std::isfinite(value) ? value : low, low, high);
  };
  t60 = safe(t60, .005f, 3.f);
  hold = safe(hold, 0.f, .08f);
  shape = safe(shape, 0.f, 1.f);
  const double linear = 1. - shape;
  const double end = (8. / 3.) /
      (std::sqrt(linear * linear + 16. / 3. * shape) + linear);
  auto &decay = amplitude.segments[1];
  decay.durationSeconds = static_cast<float>(t60 * end);
  decay.logCurvature = static_cast<float>(.75 * shape * end * end);
  if (hold > 0.f) {
    amplitude.segments[3] = amplitude.segments[2];
    amplitude.segments[2] = decay;
    amplitude.segments[1] = {1.f, hold, TrajectoryCurve::Linear};
    amplitude.segmentCount = 4;
  }
}

} // namespace tfdsp::percussion
