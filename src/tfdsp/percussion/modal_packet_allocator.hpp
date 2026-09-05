#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace tfdsp::percussion {

struct ModalPacketRequest {
  float centreErb{};
  float spreadErb{};
  bool active{};
};

template <std::size_t HandleCapacity>
struct ModalPacketAllocation {
  std::array<std::size_t, HandleCapacity> sidebandPairs{};
  std::size_t activeHandleCount{};
  std::size_t stateCount{};
};

// Allocates a fixed real-time state pool during parameter preparation. Painted
// centres always win one state; the remaining budget follows useful ERB
// coverage, with overlap reducing redundant satellites rather than removing
// deliberately close centre handles.
template <std::size_t HandleCapacity>
ModalPacketAllocation<HandleCapacity> AllocateModalPackets(
    const std::array<ModalPacketRequest, HandleCapacity> &requests,
    const std::size_t stateCapacity, const float density) noexcept {
  ModalPacketAllocation<HandleCapacity> result;
  std::array<float, HandleCapacity> desiredPairs{};
  float desiredTotal = 0.f;
  for (std::size_t handle = 0; handle < HandleCapacity; ++handle) {
    if (!requests[handle].active) continue;
    ++result.activeHandleCount;
    float nearest = 1.e6f;
    for (std::size_t other = 0; other < HandleCapacity; ++other) {
      if (other == handle || !requests[other].active) continue;
      nearest = std::min(nearest, std::abs(
          requests[handle].centreErb - requests[other].centreErb));
    }
    const float spread = std::max(requests[handle].spreadErb, 0.f);
    const float uniqueCoverage = std::clamp(
        nearest / std::max(2.f * spread, .25f), 0.f, 1.f);
    const float overlapScale = .5f + .5f * uniqueCoverage;
    desiredPairs[handle] = 4.f * std::clamp(density, 0.f, 1.f) *
        spread * overlapScale;
    desiredTotal += desiredPairs[handle];
  }

  result.stateCount = std::min(result.activeHandleCount, stateCapacity);
  if (result.activeHandleCount >= stateCapacity || !(desiredTotal > 0.f))
    return result;
  const std::size_t pairCapacity =
      (stateCapacity - result.activeHandleCount) / 2;
  const float scale = std::min(
      1.f, static_cast<float>(pairCapacity) / desiredTotal);
  std::array<float, HandleCapacity> remainder{};
  std::size_t allocatedPairs = 0;
  for (std::size_t handle = 0; handle < HandleCapacity; ++handle) {
    const float target = desiredPairs[handle] * scale;
    const auto pairs = static_cast<std::size_t>(std::floor(target));
    result.sidebandPairs[handle] = pairs;
    remainder[handle] = target - static_cast<float>(pairs);
    allocatedPairs += pairs;
  }
  while (allocatedPairs < pairCapacity) {
    const auto best = std::max_element(remainder.begin(), remainder.end());
    if (best == remainder.end() || !(*best > 0.f)) break;
    const std::size_t handle = static_cast<std::size_t>(
        std::distance(remainder.begin(), best));
    ++result.sidebandPairs[handle];
    *best = 0.f;
    ++allocatedPairs;
  }
  result.stateCount += 2 * allocatedPairs;
  return result;
}

} // namespace tfdsp::percussion
