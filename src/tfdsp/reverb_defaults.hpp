#pragma once

#include <array>

namespace tfdsp::reverb_defaults {

// The single authoritative factory baseline shared by the Rack module, DSP
// engines, and native response tests. All normalized controls are in [0, 1].
inline constexpr float Space = 0.5f;
inline constexpr float Aspect = 0.5f;
// Retained only as the serialized default for legacy parameter slot 2. Room
// height is now derived from Space and is not a DSP control.
inline constexpr float Height = 0.28f;
inline constexpr float PreDelay = 0.f;
inline constexpr std::array<float, 3> Source{{0.5f, 0.35f, 0.42f}};
inline constexpr std::array<float, 3> Listener{{0.5f, 0.682f, 0.45f}};
inline constexpr float Decay = 0.55f;
inline constexpr float Damping = 0.18f;
inline constexpr float Diffusion = 0.75f;
inline constexpr float Modulation = 0.4f;
inline constexpr float Shimmer = 0.f;
inline constexpr float Width = 0.6130368569f;
inline constexpr float EarlyLevelDb = 0.f;
inline constexpr float TailLevelDb = 0.f;
inline constexpr float LowCut = 0.f;
inline constexpr float HighCut = 0.9039693650f;
inline constexpr float Mix = 0.35f;
inline constexpr float OutputLevelDb = 0.f;

// RoomReverb::MakeRoom() evaluated at the factory geometry. Keeping the
// standalone late engine on this baseline prevents its tests from silently
// exercising an unrelated room.
inline constexpr std::array<float, 3> RoomDimensionsMetres{
    {7.09929574f, 9.35414347f, 4.38178046f}};

} // namespace tfdsp::reverb_defaults
