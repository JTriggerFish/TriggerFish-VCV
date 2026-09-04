#pragma once

#include "contact_exciter.hpp"
#include "correlated_fm_burst.hpp"
#include "membrane_resonator.hpp"
#include "observation_equalizer.hpp"
#include "observation_model.hpp"
#include "strike_energy_envelope.hpp"

#include <array>
#include <cstddef>

namespace tfdsp::percussion {

inline constexpr std::size_t MembraneModeCount = 16;

enum class MembraneDrumRoute : std::size_t {
  ContactToDirect,
  ContactToBody,
  FmToDirect,
  FmToBody,
  BodyToObservation,
  Count
};

struct MembraneDrumRouting {
  static constexpr std::size_t Count =
      static_cast<std::size_t>(MembraneDrumRoute::Count);
  bool Enabled(MembraneDrumRoute route) const noexcept;
  void SetEnabled(std::size_t index, bool value) noexcept;
  std::array<bool, Count> enabled{true, true, true, true, true};
};

struct MembraneDrumControls {
  float fundamentalHz{105.f};
  float decaySeconds{1.15f};
  float decayTilt{.55f};
  float inharmonicity{.35f};
  float bodyBrightness{.55f};
  float tensionOctaves{.11f};
  float tensionDecaySeconds{.13f};
  float contactDirectLevel{.245f};
  float contactBodyLevel{.7f};
  float contactDurationSeconds{.004f};
  float contactBrightness{.58f};
  float fmDirectLevel{.0144f};
  float fmBodyLevel{.081f};
  float fmDepthHz{260.f};
  float fmDecaySeconds{.07f};
  float pitchDropOctaves{.28f};
  float directLevel{.3f};
  float bodyLevel{1.f};
  float directDelaySeconds{};
  ObservationEqualizerMode equalizerMode{ObservationEqualizerMode::Radiation};
  float lowCutHz{24.f};
  float highCutHz{18000.f};
  float colourFrequencyHz{2800.f};
  float colourGainDb{};
  float outputGain{.08f};
};

struct MembraneDrumParameters {
  ContactExciterParameters contact{};
  CorrelatedFmBurstParameters fm{};
  MembraneResonator<MembraneModeCount>::Parameters membrane{};
  StrikeEnergyEnvelopeParameters strikeEnergy{};
  ObservationModel<2>::Parameters observation{};
  ObservationEqualizerParameters equalizer{};
  MembraneDrumRouting routing{};
  float contactDirectLevel{.245f};
  float contactBodyLevel{.7f};
  float fmDirectLevel{.0144f};
  float fmBodyLevel{.081f};
  float outputGain{.08f};
  float maximumModalEnergy{64.f};
};

struct MembraneDrumPreparedParameters {
  MembraneDrumParameters parameters{};
  MembraneResonator<MembraneModeCount>::PreparedParameters membrane{};
  float sampleRate{48000.f};
};

MembraneDrumParameters DefaultMembraneDrumParameters(
    const MembraneDrumControls &controls = {}) noexcept;
MembraneDrumPreparedParameters
PrepareMembraneDrumParameters(float sampleRate,
                              const MembraneDrumParameters &parameters);

} // namespace tfdsp::percussion
