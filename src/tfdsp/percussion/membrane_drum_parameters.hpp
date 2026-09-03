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
  float Get(MembraneDrumRoute route) const noexcept;
  void Set(std::size_t index, float gain) noexcept;
  std::array<float, Count> gains{.35f, 1.f, .08f, .45f, 1.f};
};

struct MembraneDrumControls {
  float fundamentalHz{105.f};
  float decaySeconds{1.15f};
  float decayTilt{.55f};
  float inharmonicity{.35f};
  float bodyBrightness{.55f};
  float tensionOctaves{.11f};
  float tensionDecaySeconds{.13f};
  float contactLevel{.7f};
  float contactDurationSeconds{.004f};
  float contactBrightness{.58f};
  float fmLevel{.18f};
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
  float outputGain{.32f};
};

struct MembraneDrumParameters {
  ContactExciterParameters contact{};
  CorrelatedFmBurstParameters fm{};
  MembraneResonator<MembraneModeCount>::Parameters membrane{};
  StrikeEnergyEnvelopeParameters strikeEnergy{};
  ObservationModel<2>::Parameters observation{};
  ObservationEqualizerParameters equalizer{};
  MembraneDrumRouting routing{};
  float contactLevel{.7f};
  float directVelocityExponent{.72f};
  float bodyVelocityExponent{.72f};
  float velocitySaturation{};
  float fmLevel{.18f};
  float outputGain{.32f};
  float maximumModalEnergy{1.f};
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
