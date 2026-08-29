#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace tfdsp {

inline constexpr std::size_t EventSignalCount = 5;
inline constexpr std::size_t EventVoiceLimit = 16;

// A performance-event bundle keeps Pitch, Gate, Trigger, Velocity, and Accent
// channel-aligned. Merging whole bundles avoids the invalid pitch-voltage
// addition produced by stacking independent polyphonic cables at a Rack input.
struct EventBundle {
  std::array<std::array<float, EventVoiceLimit>, EventSignalCount> signals{};
  std::size_t voiceCount{};
};

inline EventBundle MergeEventBundles(const EventBundle &first,
                                     const EventBundle &second) noexcept {
  EventBundle merged;
  const auto append = [&](const EventBundle &source) {
    const std::size_t count = std::min(source.voiceCount, EventVoiceLimit);
    for (std::size_t voice = 0;
         voice < count && merged.voiceCount < EventVoiceLimit; ++voice) {
      for (std::size_t signal = 0; signal < EventSignalCount; ++signal) {
        const float value = source.signals[signal][voice];
        merged.signals[signal][merged.voiceCount] =
            std::isfinite(value) ? value : 0.f;
      }
      ++merged.voiceCount;
    }
  };
  append(first);
  append(second);
  return merged;
}

} // namespace tfdsp
