#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cstddef>

namespace tfdsp::percussion {

// A bounded, allocation-free mixer for compiled percussion graphs.
template <std::size_t SourceCount> class FixedMixer {
public:
  using Frame = std::array<float, SourceCount>;

  void SetGains(const Frame &gains) noexcept {
    for (std::size_t source = 0; source < SourceCount; ++source)
      gains_[source] = std::clamp(
          tfdsp::FiniteNormalOrZero(gains[source]), -16.f, 16.f);
  }

  float Process(const Frame &sources) const noexcept {
    float output = 0.f;
    for (std::size_t source = 0; source < SourceCount; ++source)
      output += gains_[source] * tfdsp::FiniteNormalOrZero(sources[source]);
    return tfdsp::FiniteNormalOrZero(output);
  }

private:
  Frame gains_{};
};

} // namespace tfdsp::percussion
