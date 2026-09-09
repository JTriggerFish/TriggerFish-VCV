#pragma once
#include "tfdsp/finite_audio.hpp"
#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {
struct TurbulenceResponse {
  float intensity{}, diffuseEnergy{}, exchangeAmount{};
};

// Classic is retained for existing patches. Relaxed intensity has no knee at
// one: only the energy partition approaches a unit upper bound.
inline TurbulenceResponse EvaluateTurbulence(float frequency, float level,
    float slope, float centre, float local, bool relaxed) noexcept {
  const auto finite = [](float x) { return tfdsp::FiniteNormalOrZero(x); };
  frequency = std::clamp(finite(frequency), 1.f, 384000.f);
  centre = std::clamp(finite(centre), 1.f, 384000.f);
  slope = std::clamp(finite(slope), -1.f, 1.f);
  local = std::clamp(finite(local), 0.f, 2.f);
  const float octaves = std::log2(frequency / centre);
  if (!relaxed) {
    const float spectral = std::clamp(
        std::clamp(finite(level), 0.f, 1.f) + slope * octaves, 0.f, 1.f);
    const float intensity = std::clamp(spectral * local, 0.f, 1.f);
    return {intensity, .9f * intensity * intensity, intensity * intensity};
  }
  // At the fixed 1-kHz pivot, legacy centre/level combinations can reach 4000.
  const float intensity = std::clamp(finite(level), 0.f, 4000.f) *
      std::exp2(slope * octaves) * local;
  // 90% satellite energy at intensity one; smoothly tends to 100%.
  const float fraction = -std::expm1(-2.302585092994f * intensity * intensity);
  return {intensity, fraction, fraction};
}
} // namespace tfdsp::percussion
