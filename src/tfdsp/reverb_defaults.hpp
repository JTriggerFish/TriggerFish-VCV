#pragma once

#include <algorithm>
#include <array>
#include <cstddef>

namespace tfdsp::reverb_defaults {

// Complete acoustic and placement snapshots for the Rack context menu.
// Frequency controls use the module's normalized logarithmic coordinates.
// Output level is omitted deliberately so applying a preset never overrides
// the user's master trim.
struct ReverbPreset {
  float space;
  float aspect;
  float preDelay;
  std::array<float, 3> source;
  std::array<float, 3> listener;
  float decay;
  float damping;
  float diffusion;
  float modulation;
  float shimmer;
  float width;
  float balance;
  float lowCut;
  float highCut;
  float mix;
};

inline constexpr float ProgressiveSourceX(const std::size_t source,
                                          const std::size_t sourceCount) {
  if (sourceCount <= 1)
    return 0.5f;
  const float halfSpan = std::min(
      0.35f, 0.12f + 0.04f * static_cast<float>(sourceCount - 2));
  const float position =
      static_cast<float>(std::min(source, sourceCount - 1)) /
      static_cast<float>(sourceCount - 1);
  return 0.5f + halfSpan * (2.f * position - 1.f);
}

// Approximately 9.35 x 12.51 x 5.24 metres, with a clear wet onset and a
// dense, gently dark tail. This is also the new-module factory default.
inline constexpr ReverbPreset MediumHall{
    0.60f,
    0.50f,
    0.048f, // 12 ms
    {{0.50f, 0.35f, 0.42f}},
    {{0.50f, 0.682f, 0.45f}},
    0.64f, // 2.30 s
    0.28f,
    0.82f,
    0.30f,
    0.f,
    0.6130368569f, // 100% displayed width
    0.50f,
    0.f,          // 20 Hz
    0.8294822177f, // 12 kHz
    0.35f};

// Approximately 3.74 x 4.76 x 2.90 metres. Early energy and a short, damped
// tail put an instrument in a pleasant real room without sounding effected.
inline constexpr ReverbPreset SmallRoom{
    0.25f,
    0.50f,
    0.f,
    {{0.38f, 0.35f, 0.42f}},
    {{0.55f, 0.68f, 0.45f}},
    0.34f, // 0.81 s
    0.40f,
    0.62f,
    0.15f,
    0.f,
    0.55f, // approximately 86% displayed width
    0.40f,
    0.1771838201f, // 40 Hz
    0.7334515827f, // 9 kHz
    0.22f};

// A large evolving texture whose pre-delay, filtering, and conservative Mix
// retain a clear dry attack despite the long modulated shimmer tail.
inline constexpr ReverbPreset Superlush{
    0.78f,
    0.55f,
    0.112f, // 28 ms
    {{0.50f, 0.32f, 0.42f}},
    {{0.50f, 0.68f, 0.45f}},
    0.84f, // 4.59 s
    0.26f,
    1.f,
    1.f,
    0.65f,
    0.75f, // approximately 127% displayed width
    0.75f,
    0.4580135308f, // 120 Hz
    0.6941346395f, // 8 kHz
    0.40f};

// The single authoritative factory baseline shared by the Rack module, DSP
// engines, native response tests, and generated smoke patches.
inline constexpr float Space = MediumHall.space;
inline constexpr float Aspect = MediumHall.aspect;
inline constexpr float PreDelay = MediumHall.preDelay;
inline constexpr auto Source = MediumHall.source;
inline constexpr auto Listener = MediumHall.listener;
inline constexpr float Decay = MediumHall.decay;
inline constexpr float Damping = MediumHall.damping;
inline constexpr float Diffusion = MediumHall.diffusion;
inline constexpr float Modulation = MediumHall.modulation;
inline constexpr float Shimmer = MediumHall.shimmer;
inline constexpr float Width = MediumHall.width;
inline constexpr float Balance = MediumHall.balance;
inline constexpr float LowCut = MediumHall.lowCut;
inline constexpr float HighCut = MediumHall.highCut;
inline constexpr float Mix = MediumHall.mix;
inline constexpr float OutputLevelDb = 0.f;

// RoomReverb::MakeRoom() evaluated at the factory geometry. Keeping the
// standalone late engine on this baseline prevents its tests from silently
// exercising an unrelated room.
inline constexpr std::array<float, 3> RoomDimensionsMetres{
    {9.35009902f, 12.51345042f, 5.23644647f}};

} // namespace tfdsp::reverb_defaults
