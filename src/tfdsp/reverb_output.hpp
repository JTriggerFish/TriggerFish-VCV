#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace tfdsp {

inline float ReverbControlSmoothstep(const float value) noexcept {
  const float limited =
      std::clamp(std::isfinite(value) ? value : 0.f, 0.f, 1.f);
  return limited * limited * (3.f - 2.f * limited);
}

inline float ReverbDecibelsToGain(const float decibels) noexcept {
  if (!std::isfinite(decibels) || decibels <= -60.f)
    return 0.f;
  return std::pow(10.f, std::clamp(decibels, -60.f, 6.f) / 20.f);
}

struct ReverbOutputGains {
  float dry{};
  float wet{};
};

inline ReverbOutputGains
CalculateReverbOutputGains(const float mixControl,
                           const float outputLevelDb) noexcept {
  constexpr float HalfPi = 1.57079632679489661923f;
  const float mix = ReverbControlSmoothstep(mixControl);
  const float level = ReverbDecibelsToGain(outputLevelDb);
  return {level * std::cos(HalfPi * mix),
          level * std::sin(HalfPi * mix)};
}

inline std::array<float, 2>
MixReverbOutput(const float dry, const std::array<float, 2> &wet,
                const ReverbOutputGains &gains) noexcept {
  return {gains.dry * dry + gains.wet * wet[0],
          gains.dry * dry + gains.wet * wet[1]};
}

inline std::array<float, 2>
MixReverbOutput(const float dry, const std::array<float, 2> &wet,
                const float mixControl, const float outputLevelDb) noexcept {
  return MixReverbOutput(
      dry, wet, CalculateReverbOutputGains(mixControl, outputLevelDb));
}

} // namespace tfdsp
