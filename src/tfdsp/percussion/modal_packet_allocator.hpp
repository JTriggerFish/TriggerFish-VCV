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
  float allocationWeight{1.f};
};

template <std::size_t HandleCapacity>
struct ModalPacketAllocation {
  std::array<std::size_t, HandleCapacity> sidebandPairs{};
  std::size_t activeHandleCount{};
  std::size_t stateCount{};
};

// Allocates a fixed real-time state pool during parameter preparation. Painted
// centres reserve one state (two for paired rings). Density chooses the budget independently
// of packet width; sqrt(width) and visible local weights distribute that budget.
// Narrow packets can therefore contain many stable beating modes without blur.
template <std::size_t HandleCapacity>
ModalPacketAllocation<HandleCapacity> AllocateModalPackets(
    const std::array<ModalPacketRequest, HandleCapacity> &requests,
    const std::size_t stateCapacity, const float density,
    const bool pairedCentres = false) noexcept {
  ModalPacketAllocation<HandleCapacity> result;
  std::array<float, HandleCapacity> desiredPairs{};
  float desiredTotal = 0.f;
  for (std::size_t handle = 0; handle < HandleCapacity; ++handle) {
    if (!requests[handle].active) continue;
    ++result.activeHandleCount;
    const float spread = requests[handle].spreadErb;
    const float weight = requests[handle].allocationWeight;
    desiredPairs[handle] = std::isfinite(spread) && spread > 0.f &&
        std::isfinite(weight) && weight > 0.f
        ? std::sqrt(spread) * std::min(weight, 4.f) : 0.f;
    desiredTotal += desiredPairs[handle];
  }

  const std::size_t centres = result.activeHandleCount * (pairedCentres ? 2 : 1);
  result.stateCount = std::min(centres, stateCapacity);
  if (centres >= stateCapacity || !(desiredTotal > 0.f))
    return result;
  const float safeDensity = std::isfinite(density) ? std::clamp(density, 0.f, 1.f) : 0.f;
  const std::size_t pairCapacity = static_cast<std::size_t>(std::round(
      float((stateCapacity - centres) / 2) * safeDensity));
  const float scale = static_cast<float>(pairCapacity) / desiredTotal;
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
