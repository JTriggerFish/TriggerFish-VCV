#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cstddef>

namespace tfdsp::percussion {

enum class MetallicPlateRoute : std::size_t {
  ContactToDispersion,
  ContactToBody,
  DispersionToBody,
  ContactToObservation,
  BodyToObservation,
  Count
};

// Prepared gains for the bounded metallic-plate recipe. The patch compiler
// owns graph validation; the voice only applies its fixed, direct call sites.
struct MetallicPlateRouting {
  static constexpr std::size_t Count =
      static_cast<std::size_t>(MetallicPlateRoute::Count);

  float Get(const MetallicPlateRoute route) const noexcept {
    return gains[static_cast<std::size_t>(route)];
  }

  void Set(const std::size_t index, const float gain) noexcept {
    if (index < gains.size())
      gains[index] = std::clamp(tfdsp::FiniteNormalOrZero(gain), 0.f, 2.f);
  }

  std::array<float, Count> gains{1.f, 1.f, 1.f, 1.f, 1.f};
};

} // namespace tfdsp::percussion
