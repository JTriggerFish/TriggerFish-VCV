#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

namespace tfdsp::late_reverb_coefficients {

// Production heuristic VFM design. The main-delay ratios come from the prime
// FDN16 sequence published by Fagerstrom et al., normalized to unit mean so
// room geometry can set their absolute mean free time. Two sparse paraunitary
// Hadamard delay/mix stages supply the short scattering time scale and follow
// the same room-scale ratio without exposing a bare short-delay FDN at minimum
// Diffusion.
inline constexpr std::uint32_t DesignVersion = 6;
inline constexpr std::size_t LineCount = 16;
inline constexpr std::size_t WallCount = 6;
inline constexpr std::size_t VelvetStageCount = 2;

inline constexpr std::array<float, LineCount> MainDelayRatio{
    0.5668773f, 0.6296269f, 0.6799816f, 0.7295616f,
    0.7899872f, 0.8504128f, 0.9201346f, 0.9681652f,
    1.0258794f, 1.0773961f, 1.1362723f, 1.2001840f,
    1.2629336f, 1.3280074f, 1.3864962f, 1.4480839f};
inline constexpr std::array<std::array<float, LineCount>, VelvetStageCount>
    VelvetDelayMs{{
        {{0.125f, 0.875f, 1.375f, 2.0f, 2.75f, 3.25f, 3.875f, 4.5f,
          5.25f, 5.75f, 6.375f, 7.125f, 7.75f, 8.25f, 9.125f, 9.75f}},
        {{0.375f, 1.0f, 1.625f, 2.375f, 2.875f, 3.625f, 4.125f, 4.875f,
          5.375f, 6.125f, 6.75f, 7.375f, 8.0f, 8.75f, 9.375f, 10.125f}}
    }};
inline constexpr std::array<std::array<std::size_t, LineCount>,
                            VelvetStageCount + 1>
    TransformPermutation{{
        {{11, 1, 15, 0, 6, 3, 4, 13, 8, 10, 14, 12, 9, 7, 2, 5}},
        {{9, 10, 12, 4, 7, 5, 13, 1, 3, 6, 0, 8, 2, 15, 11, 14}},
        {{5, 12, 1, 9, 14, 3, 10, 0, 7, 15, 4, 11, 2, 13, 8, 6}}
    }};
inline constexpr std::array<std::array<float, LineCount>,
                            VelvetStageCount + 1>
    TransformSign{{
    {{-1.f, -1.f, -1.f, -1.f, -1.f, -1.f, 1.f, -1.f,
      -1.f, -1.f, 1.f, -1.f, -1.f, -1.f, 1.f, -1.f}},
    {{1.f, -1.f, -1.f, 1.f, -1.f, 1.f, -1.f, 1.f,
      1.f, 1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f}},
    {{1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, -1.f,
      1.f, -1.f, -1.f, 1.f, -1.f, 1.f, 1.f, -1.f}}
}};

inline constexpr float InverseRootLineCount = 0.25f;
constexpr float WalshSign(std::size_t row, std::size_t mask) noexcept {
  std::size_t bits = row & mask;
  bool odd = false;
  while (bits != 0) {
    odd = !odd;
    bits &= bits - 1;
  }
  return odd ? -1.f : 1.f;
}
constexpr auto MakeWallProjection() noexcept {
  std::array<std::array<float, WallCount>, LineCount> result{};
  for (std::size_t row = 0; row < LineCount; ++row)
    for (std::size_t wall = 0; wall < WallCount; ++wall)
      result[row][wall] = WalshSign(row, wall);
  return result;
}
inline constexpr auto WallProjection = MakeWallProjection();

inline constexpr std::size_t ShimmerBusCount = 4;
constexpr auto MakeShimmerProjection() noexcept {
  std::array<std::array<float, LineCount>, ShimmerBusCount> result{};
  for (std::size_t bus = 0; bus < ShimmerBusCount; ++bus)
    for (std::size_t line = 0; line < LineCount; ++line)
      result[bus][line] = InverseRootLineCount * WalshSign(line, 1u << bus);
  return result;
}
inline constexpr auto ShimmerProjection = MakeShimmerProjection();

inline constexpr std::array<float, LineCount> ModulationRateHz{
    0.061f, 0.067f, 0.073f, 0.079f, 0.089f, 0.097f, 0.107f, 0.117f,
    0.127f, 0.139f, 0.151f, 0.163f, 0.179f, 0.193f, 0.211f, 0.229f};
inline constexpr float UnitHandoffPower = 2.5e-5f;
inline constexpr float TailOutputCalibration = 1.30f;

} // namespace tfdsp::late_reverb_coefficients
