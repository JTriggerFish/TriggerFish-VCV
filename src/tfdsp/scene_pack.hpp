#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace tfdsp {

inline constexpr std::size_t ScenePackInputCount = 4;
inline constexpr std::size_t ScenePackMaximumSourceCount = 8;

// Each physical input may itself be polyphonic. Connected channels are
// concatenated in input order, retaining channel order within each input.
// This makes the four-jack utility chainable without a separate scene bus.
struct ScenePackInput {
  std::array<std::array<float, ScenePackMaximumSourceCount>,
             ScenePackInputCount>
      inputs{};
  std::array<std::size_t, ScenePackInputCount> channelCounts{};
};

struct ScenePackOutput {
  std::array<float, ScenePackMaximumSourceCount> audio{};
  std::size_t sourceCount{};
};

inline float SanitizeSceneAudio(const float value) noexcept {
  return std::isfinite(value) ? value : 0.0f;
}

inline ScenePackOutput PackSceneSources(const ScenePackInput &input) noexcept {
  ScenePackOutput output;
  for (std::size_t port = 0;
       port < ScenePackInputCount &&
       output.sourceCount < ScenePackMaximumSourceCount;
       ++port) {
    const std::size_t channels = std::min(
        input.channelCounts[port], ScenePackMaximumSourceCount);
    for (std::size_t channel = 0;
         channel < channels &&
         output.sourceCount < ScenePackMaximumSourceCount;
         ++channel)
      output.audio[output.sourceCount++] =
          SanitizeSceneAudio(input.inputs[port][channel]);
  }
  return output;
}

} // namespace tfdsp
