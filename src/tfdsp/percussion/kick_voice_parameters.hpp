#pragma once

#include "membrane_drum_parameters.hpp"

#include <array>

namespace tfdsp::percussion {

struct KickModeControl {
  float frequencyHz{55.f};
  float levelDb{-72.f}; // relative observation; -72 disables this mode
  float centreCoupling{1.f};
  float edgeCoupling{1.f};
};

std::array<KickModeControl, MembraneModeCount> DefaultKickModes() noexcept;

// Functional kick surface. Levels are observation gains, not extra drive
// stages. One unit contact drives the passive body; thump never drives it.
struct KickVoiceControls {
  float contactLevel{.4f};
  float contactWidthSeconds{.011f};
  float contactColour{.33f};
  float contactNoiseLevel{.57f};
  float contactNoiseDecaySeconds{.205f}; // base T60, before implement shaping
  float thumpLevel{1.77f};
  float thumpPitchHz{28.f};
  float thumpPitchDropOctaves{1.47f};
  float thumpPitchFallSeconds{.059f};
  float thumpDecaySeconds{.306f}; // amplitude T60
  float resonanceLevel{4.72f};
  float resonanceDecaySeconds{.6f}; // T60 at 100 Hz
  float resonanceDecayTilt{.5f}; // T60 octaves lost per frequency octave
  std::array<KickModeControl, MembraneModeCount> modes{DefaultKickModes()};
  float tensionOctaves{.056f};
  float tensionRecoverySeconds{.012f};
  float outputGain{.25f};
};

struct KickVoiceRouting {
  std::array<bool, 3> enabled{true, true, true}; // contact, thump, resonance
  void SetEnabled(std::size_t index, bool value) noexcept {
    if (index < enabled.size())
      enabled[index] = value;
  }
};

MembraneDrumParameters
DefaultKickVoiceParameters(const KickVoiceControls &controls = {}) noexcept;
void ApplyKickRouting(MembraneDrumParameters &parameters,
                      const KickVoiceRouting &routing) noexcept;

} // namespace tfdsp::percussion
