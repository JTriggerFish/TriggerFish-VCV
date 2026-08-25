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

enum class HeatmapPalette {
  Magma = 0,
  Inferno,
  Plasma,
  Viridis,
  Cividis,
  CrtGreen,
  CrtBlue,
  CrtYellow,
  CrtRed,
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

// Uniform samples 0, 17, ... 255 from Matplotlib's canonical 256-entry
// perceptually uniform tables in matplotlib/lib/matplotlib/_cm_listed.py.
inline constexpr std::array<HeatmapRgb, 16> MagmaHeatmap{{
    {0.001462f, 0.000466f, 0.013866f},
    {0.043830f, 0.033830f, 0.141886f},
    {0.123833f, 0.067295f, 0.295879f},
    {0.232077f, 0.059889f, 0.437695f},
    {0.341482f, 0.080564f, 0.492631f},
    {0.445163f, 0.122724f, 0.506901f},
    {0.550287f, 0.161158f, 0.505719f},
    {0.658483f, 0.196027f, 0.490253f},
    {0.767398f, 0.233705f, 0.457755f},
    {0.868793f, 0.287728f, 0.409303f},
    {0.944006f, 0.377643f, 0.365136f},
    {0.981000f, 0.498428f, 0.369734f},
    {0.994738f, 0.624350f, 0.427397f},
    {0.997228f, 0.747981f, 0.516859f},
    {0.993170f, 0.870024f, 0.626189f},
    {0.987053f, 0.991438f, 0.749504f},
}};

inline constexpr std::array<HeatmapRgb, 16> InfernoHeatmap{{
    {0.001462f, 0.000466f, 0.013866f},
    {0.046915f, 0.030324f, 0.150164f},
    {0.142378f, 0.046242f, 0.308553f},
    {0.258234f, 0.038571f, 0.406485f},
    {0.366529f, 0.071579f, 0.431994f},
    {0.472328f, 0.110547f, 0.428334f},
    {0.578304f, 0.148039f, 0.404411f},
    {0.682656f, 0.189501f, 0.360757f},
    {0.780517f, 0.243327f, 0.299523f},
    {0.865006f, 0.316822f, 0.226055f},
    {0.929644f, 0.411479f, 0.145367f},
    {0.970919f, 0.522853f, 0.058367f},
    {0.987622f, 0.645320f, 0.039886f},
    {0.978806f, 0.774545f, 0.176037f},
    {0.950018f, 0.903409f, 0.380271f},
    {0.988362f, 0.998364f, 0.644924f},
}};

inline constexpr std::array<HeatmapRgb, 16> PlasmaHeatmap{{
    {0.050383f, 0.029803f, 0.527975f},
    {0.200445f, 0.017902f, 0.593364f},
    {0.312543f, 0.008239f, 0.635700f},
    {0.417642f, 0.000564f, 0.658390f},
    {0.517933f, 0.021563f, 0.654109f},
    {0.610667f, 0.090204f, 0.619951f},
    {0.692840f, 0.165141f, 0.564522f},
    {0.764193f, 0.240396f, 0.502126f},
    {0.826588f, 0.315714f, 0.441316f},
    {0.881443f, 0.392529f, 0.383229f},
    {0.928329f, 0.472975f, 0.326067f},
    {0.965024f, 0.559118f, 0.268513f},
    {0.988260f, 0.652325f, 0.211364f},
    {0.994141f, 0.753137f, 0.161404f},
    {0.977995f, 0.861432f, 0.142808f},
    {0.940015f, 0.975158f, 0.131326f},
}};

inline constexpr std::array<HeatmapRgb, 16> ViridisHeatmap{{
    {0.267004f, 0.004874f, 0.329415f},
    {0.282656f, 0.100196f, 0.422160f},
    {0.277134f, 0.185228f, 0.489898f},
    {0.253935f, 0.265254f, 0.529983f},
    {0.221989f, 0.339161f, 0.548752f},
    {0.190631f, 0.407061f, 0.556089f},
    {0.163625f, 0.471133f, 0.558148f},
    {0.139147f, 0.533812f, 0.555298f},
    {0.120565f, 0.596422f, 0.543611f},
    {0.134692f, 0.658636f, 0.517649f},
    {0.208030f, 0.718701f, 0.472873f},
    {0.327796f, 0.773980f, 0.406640f},
    {0.477504f, 0.821444f, 0.318195f},
    {0.647257f, 0.858400f, 0.209861f},
    {0.824940f, 0.884720f, 0.106217f},
    {0.993248f, 0.906157f, 0.143936f},
}};

inline constexpr std::array<HeatmapRgb, 16> CividisHeatmap{{
    {0.000000f, 0.135112f, 0.304751f},
    {0.000000f, 0.181610f, 0.421859f},
    {0.117612f, 0.225935f, 0.434308f},
    {0.208926f, 0.272546f, 0.424809f},
    {0.279411f, 0.318677f, 0.423031f},
    {0.342246f, 0.364939f, 0.428559f},
    {0.401418f, 0.411790f, 0.440708f},
    {0.458366f, 0.459552f, 0.460457f},
    {0.517920f, 0.508454f, 0.472707f},
    {0.582087f, 0.558670f, 0.468118f},
    {0.648222f, 0.610553f, 0.454801f},
    {0.716177f, 0.664384f, 0.432386f},
    {0.785965f, 0.720438f, 0.399613f},
    {0.857809f, 0.778969f, 0.353259f},
    {0.932180f, 0.840159f, 0.285880f},
    {0.995737f, 0.909344f, 0.217772f},
}};

// Monochrome CRT ramps begin at an unlit black screen, rise through the
// saturated phosphor colour, and become paler near maximum excitation. The
// fourth-power bloom term confines that whitening to the bright end.
inline constexpr std::array<HeatmapRgb, 16>
makeCrtHeatmap(const HeatmapRgb phosphor, const HeatmapRgb highlight) noexcept {
  std::array<HeatmapRgb, 16> result{};
  for (std::size_t index = 0; index < result.size(); ++index) {
    const float intensity =
        static_cast<float>(index) / static_cast<float>(result.size() - 1);
    const float emission =
        intensity * intensity * (3.f - 2.f * intensity);
    const float squared = intensity * intensity;
    const float bloom = squared * squared;
    const auto channel = [=](const float base, const float bright) {
      return emission * (base + bloom * (bright - base));
    };
    result[index] = {channel(phosphor.red, highlight.red),
                     channel(phosphor.green, highlight.green),
                     channel(phosphor.blue, highlight.blue)};
  }
  return result;
}

inline constexpr auto CrtGreenHeatmap =
    makeCrtHeatmap({0.05f, 0.92f, 0.18f}, {0.72f, 1.f, 0.76f});
inline constexpr auto CrtBlueHeatmap =
    makeCrtHeatmap({0.05f, 0.42f, 1.f}, {0.70f, 0.88f, 1.f});
inline constexpr auto CrtYellowHeatmap =
    makeCrtHeatmap({1.f, 0.72f, 0.03f}, {1.f, 0.98f, 0.68f});
inline constexpr auto CrtRedHeatmap =
    makeCrtHeatmap({1.f, 0.08f, 0.03f}, {1.f, 0.72f, 0.62f});

inline constexpr HeatmapPalette heatmapPaletteFromInt(int value) noexcept {
  switch (value) {
  case static_cast<int>(HeatmapPalette::Inferno):
    return HeatmapPalette::Inferno;
  case static_cast<int>(HeatmapPalette::Plasma):
    return HeatmapPalette::Plasma;
  case static_cast<int>(HeatmapPalette::Viridis):
    return HeatmapPalette::Viridis;
  case static_cast<int>(HeatmapPalette::Cividis):
    return HeatmapPalette::Cividis;
  case static_cast<int>(HeatmapPalette::CrtGreen):
    return HeatmapPalette::CrtGreen;
  case static_cast<int>(HeatmapPalette::CrtBlue):
    return HeatmapPalette::CrtBlue;
  case static_cast<int>(HeatmapPalette::CrtYellow):
    return HeatmapPalette::CrtYellow;
  case static_cast<int>(HeatmapPalette::CrtRed):
    return HeatmapPalette::CrtRed;
  default:
    return HeatmapPalette::Magma;
  }
}

inline constexpr const char *
heatmapPaletteName(const HeatmapPalette palette) noexcept {
  switch (palette) {
  case HeatmapPalette::Inferno:
    return "Inferno";
  case HeatmapPalette::Plasma:
    return "Plasma";
  case HeatmapPalette::Viridis:
    return "Viridis";
  case HeatmapPalette::Cividis:
    return "Cividis";
  case HeatmapPalette::CrtGreen:
    return "CRT Green";
  case HeatmapPalette::CrtBlue:
    return "CRT Blue";
  case HeatmapPalette::CrtYellow:
    return "CRT Yellow";
  case HeatmapPalette::CrtRed:
    return "CRT Red";
  case HeatmapPalette::Magma:
  default:
    return "Magma";
  }
}

inline HeatmapRgb sampleHeatmap(HeatmapPalette palette,
                                float intensity) noexcept {
  switch (palette) {
  case HeatmapPalette::Inferno:
    return sampleHeatmap(InfernoHeatmap, intensity);
  case HeatmapPalette::Plasma:
    return sampleHeatmap(PlasmaHeatmap, intensity);
  case HeatmapPalette::Viridis:
    return sampleHeatmap(ViridisHeatmap, intensity);
  case HeatmapPalette::Cividis:
    return sampleHeatmap(CividisHeatmap, intensity);
  case HeatmapPalette::CrtGreen:
    return sampleHeatmap(CrtGreenHeatmap, intensity);
  case HeatmapPalette::CrtBlue:
    return sampleHeatmap(CrtBlueHeatmap, intensity);
  case HeatmapPalette::CrtYellow:
    return sampleHeatmap(CrtYellowHeatmap, intensity);
  case HeatmapPalette::CrtRed:
    return sampleHeatmap(CrtRedHeatmap, intensity);
  case HeatmapPalette::Magma:
  default:
    return sampleHeatmap(MagmaHeatmap, intensity);
  }
}

// Backward-compatible alias for callers that want the default palette.
inline constexpr const auto &ProgramEditorHeatmap = MagmaHeatmap;

} // namespace tfui
