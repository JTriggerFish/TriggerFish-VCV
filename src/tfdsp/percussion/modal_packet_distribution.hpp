#pragma once

#include "erb_scale.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace tfdsp::percussion {

enum class ModalPacketDistribution { Scattered, Even, Doublets, PairedRing, PlateCloud };

// An explicit constructive law, not a physical claim about doublet splitting.
// Zero tilt preserves a shared Hz gap; .5 doubles rate over two octaves.
inline float RingBeatRate(const float centre, const float baseHz,
                         const float tilt) noexcept {
  return std::clamp(baseHz * std::pow(std::max(centre, 1.f) / 125.f,
                                    std::clamp(tilt, -1.f, 1.f)), 0.f, 80.f);
}

inline bool HasPairedRing(const ModalPacketDistribution layout,
                          const float splitHz) noexcept {
  return (layout == ModalPacketDistribution::PairedRing ||
          layout == ModalPacketDistribution::PlateCloud) &&
      std::isfinite(splitHz) && splitHz > 0.f;
}

// A complete doublet keeps the same sum of squared excitation weights at
// every depth. An odd, unpaired sideband retains its original weight.
inline float DoubletWeightScale(std::size_t pair, std::size_t count,
                                float depth) noexcept {
  if ((count & 1u) && pair + 1 == count) return 1.f;
  depth = std::clamp(std::isfinite(depth) ? depth : 0.f, 0.f, 1.f);
  const float scale = std::sqrt(2.f / (1.f + depth * depth));
  return (pair & 1u) ? scale : depth * scale;
}

// Keep the painted pitch at the midpoint and compress the pair at boundaries.
inline float RingHalfSplit(const float centre, const float splitHz,
                           const float maximumHz) noexcept {
  return std::min(.5f * std::clamp(splitHz, 0.f, 80.f),
      std::max(0.f, std::min(centre - 1.f, maximumHz - centre)));
}

// Nested low-discrepancy positions keep existing frequencies stable as the
// shared pool grows. Unlike a count-dependent regular grid they do not retune
// every oscillator on a one-state allocation change.
inline float PacketRadius(std::size_t index) noexcept {
  if (index == 0) return 1.f;
  float result = 0.f, place = .5f;
  while (index) {
    result += place * float(index & 1u);
    index >>= 1u;
    place *= .5f;
  }
  return result;
}

inline float PacketSideFrequency(const float centre, const float spreadErb,
    const std::size_t pair, const float side, const ModalPacketDistribution layout,
    const float splitHz, const float jitter, const float maximumHz) noexcept {
  if (layout == ModalPacketDistribution::PlateCloud) {
    // A finite collection of stable resonances, approximately uniform per Hz.
    // Width is the local ERB bandwidth times the existing spread control.
    // Compress each side at boundaries instead of piling modes on a clamp.
    const float width = spreadErb * 24.7f * (1.f + .00437f * centre);
    const float room = side < 0 ? centre - 1.f : maximumHz - centre;
    return centre + side * std::min(width, std::max(0.f, room)) *
        PacketRadius(pair) * jitter;
  }
  const bool doublets = layout == ModalPacketDistribution::Doublets;
  const float radius = doublets ? .05f + .9f * PacketRadius(pair / 2 + 1)
      : PacketRadius(pair);
  const float scatter = layout == ModalPacketDistribution::Scattered ||
      layout == ModalPacketDistribution::PairedRing ? jitter : 1.f;
  // Compress the support at boundaries rather than stacking clipped modes on
  // DC/Nyquist. An ERB packet is symmetric only away from those boundaries.
  const float room = side < 0 ? ErbRate(centre) - ErbRate(1.f)
                             : ErbRate(maximumHz) - ErbRate(centre);
  float frequency = InverseErbRate(ErbRate(centre) + side *
      std::min(spreadErb, std::max(0.f, room)) * radius * scatter);
  if (doublets) {
    const float halfSplit = std::min(.5f * std::max(0.f, splitHz),
        .45f * std::min(frequency - 1.f, maximumHz - frequency));
    frequency += (pair & 1u) ? halfSplit : -halfSplit;
  }
  return std::clamp(frequency, 1.f, maximumHz);
}

} // namespace tfdsp::percussion
