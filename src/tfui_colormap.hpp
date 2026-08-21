#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace tfui {

struct HeatmapRgb {
  float red = 0.f;
  float green = 0.f;
  float blue = 0.f;
};

template <std::size_t Size>
HeatmapRgb sampleHeatmap(const std::array<HeatmapRgb, Size> &heatmap,
                         float intensity) noexcept {
  static_assert(Size >= 2, "a heatmap requires at least two samples");
  intensity = std::clamp(intensity, 0.f, 1.f);
  const float position = intensity * static_cast<float>(Size - 1);
  const auto lower = static_cast<std::size_t>(std::floor(position));
  const auto upper = std::min(lower + 1, Size - 1);
  const float mix = position - static_cast<float>(lower);
  const auto interpolate = [mix](float left, float right) {
    return left + (right - left) * mix;
  };
  return {interpolate(heatmap[lower].red, heatmap[upper].red),
          interpolate(heatmap[lower].green, heatmap[upper].green),
          interpolate(heatmap[lower].blue, heatmap[upper].blue)};
}

// Uniform samples 0, 17, ... 255 from Matplotlib's canonical 256-entry magma
// table. UI code supplies only normalized intensity; replacing this one table
// changes the complete editor colour mapping without changing rendering logic.
// Source: matplotlib/lib/matplotlib/_cm_listed.py (_magma_data).
inline constexpr std::array<HeatmapRgb, 16> MagmaHeatmap{{
    {0.001462f, 0.000466f, 0.013866f},
    {0.043830f, 0.033830f, 0.141886f},
    {0.123833f, 0.067295f, 0.295879f},
    {0.225302f, 0.060445f, 0.431742f},
    {0.335308f, 0.078236f, 0.491024f},
    {0.439062f, 0.120298f, 0.506555f},
    {0.550287f, 0.161158f, 0.505719f},
    {0.658483f, 0.196027f, 0.490253f},
    {0.767398f, 0.233705f, 0.457755f},
    {0.863320f, 0.283729f, 0.412403f},
    {0.944006f, 0.377643f, 0.365136f},
    {0.981000f, 0.498428f, 0.369734f},
    {0.994738f, 0.624350f, 0.427397f},
    {0.997228f, 0.747981f, 0.516859f},
    {0.993545f, 0.862859f, 0.619299f},
    {0.987053f, 0.991438f, 0.749504f},
}};

// This alias is the only palette choice made by the program editor.
inline constexpr const auto &ProgramEditorHeatmap = MagmaHeatmap;

} // namespace tfui
