#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace tfdsp::percussion {

// Fixed linear observation of resonator-line outputs. It exposes logical
// groups for filtering or frequency shifting without creating another body.
template <std::size_t LineCount, std::size_t BusCount> class ResonatorSubmix {
public:
  using LineFrame = std::array<float, LineCount>;
  using BusFrame = std::array<float, BusCount>;
  using Weights = std::array<LineFrame, BusCount>;

  void SetWeights(const Weights &weights) noexcept {
    for (std::size_t bus = 0; bus < BusCount; ++bus)
      for (std::size_t line = 0; line < LineCount; ++line)
        weights_[bus][line] = std::clamp(
            tfdsp::FiniteNormalOrZero(weights[bus][line]), -16.f, 16.f);
  }

  BusFrame Process(const LineFrame &lines) const noexcept {
    BusFrame buses{};
    for (std::size_t bus = 0; bus < BusCount; ++bus) {
      float sum = 0.f;
      for (std::size_t line = 0; line < LineCount; ++line)
        sum += weights_[bus][line] * tfdsp::FiniteNormalOrZero(lines[line]);
      buses[bus] = tfdsp::FiniteNormalOrZero(sum);
    }
    return buses;
  }

private:
  Weights weights_{};
};

} // namespace tfdsp::percussion
