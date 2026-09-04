#pragma once

#include <array>
#include <cstddef>

namespace tfdsp::percussion {

enum class MetallicPlateRoute : std::size_t {
  ContactToBody,
  ContactToObservation,
  BodyToObservation,
  Count
};

// Prepared on/off switches for the bounded metallic-plate recipe. Signal
// levels belong to the connected modules, never to graph connections.
struct MetallicPlateRouting {
  static constexpr std::size_t Count =
      static_cast<std::size_t>(MetallicPlateRoute::Count);

  bool Enabled(const MetallicPlateRoute route) const noexcept {
    return enabled[static_cast<std::size_t>(route)];
  }

  void SetEnabled(const std::size_t index, const bool value) noexcept {
    if (index < enabled.size()) enabled[index] = value;
  }

  std::array<bool, Count> enabled{true, true, true};
};

} // namespace tfdsp::percussion
