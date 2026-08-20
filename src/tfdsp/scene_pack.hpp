#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace tfdsp {

inline constexpr std::size_t ScenePackLocalSourceCount = 4;
inline constexpr std::size_t ScenePackMaximumSourceCount = 8;

struct ScenePackSource {
  float audio{};
  float x{5.0f};
  float y{5.0f};
  float z{5.0f};
};

struct ScenePackInput {
  std::array<ScenePackSource, ScenePackMaximumSourceCount> bus{};
  std::size_t busCount{};
  std::array<ScenePackSource, ScenePackLocalSourceCount> local{};
  std::array<bool, ScenePackLocalSourceCount> localConnected{};
};

struct ScenePackOutput {
  std::array<ScenePackSource, ScenePackMaximumSourceCount> sources{};
  std::size_t sourceCount{};
};

inline float SanitizeSceneAudio(const float value) noexcept {
  return std::isfinite(value) ? value : 0.0f;
}

inline float SanitizeScenePosition(const float value) noexcept {
  return std::clamp(std::isfinite(value) ? value : 5.0f, 0.0f, 10.0f);
}

inline ScenePackSource SanitizeSceneSource(ScenePackSource source) noexcept {
  source.audio = SanitizeSceneAudio(source.audio);
  source.x = SanitizeScenePosition(source.x);
  source.y = SanitizeScenePosition(source.y);
  source.z = SanitizeScenePosition(source.z);
  return source;
}

inline ScenePackOutput PackSceneSources(const ScenePackInput &input) noexcept {
  ScenePackOutput output;
  const std::size_t busCount =
      std::min(input.busCount, ScenePackMaximumSourceCount);
  for (; output.sourceCount < busCount; ++output.sourceCount)
    output.sources[output.sourceCount] =
        SanitizeSceneSource(input.bus[output.sourceCount]);
  for (std::size_t lane = 0;
       lane < ScenePackLocalSourceCount &&
       output.sourceCount < ScenePackMaximumSourceCount;
       ++lane) {
    if (!input.localConnected[lane])
      continue;
    output.sources[output.sourceCount++] =
        SanitizeSceneSource(input.local[lane]);
  }
  return output;
}

} // namespace tfdsp
