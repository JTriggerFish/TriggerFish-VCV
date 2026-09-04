#pragma once

#include "correlated_fm_burst.hpp"
#include "enveloped_noise_burst.hpp"
#include "observation_model.hpp"

#include <array>
#include <cstddef>

namespace tfdsp::percussion {

enum class CompactKickRoute : std::size_t {
  PrimaryToMix,
  SecondaryToMix,
  ClickToMix,
  Count
};

struct CompactKickRouting {
  static constexpr std::size_t Count =
      static_cast<std::size_t>(CompactKickRoute::Count);
  bool Enabled(CompactKickRoute route) const noexcept;
  void SetEnabled(std::size_t index, bool value) noexcept;
  std::array<bool, Count> enabled{true, true, true};
};

struct CompactKickControls {
  float primaryLevel{1.f};
  float fundamentalHz{52.f};
  float pitchDropOctaves{1.8f};
  float pitchDecaySeconds{.055f};
  float bodyDecaySeconds{.38f};
  float fmDepthHz{720.f};
  float fmDecaySeconds{.035f};
  float fmRoughnessHz{4200.f};
  float secondaryRatio{1.52f};
  float secondaryLevel{.32f};
  float clickLevel{.16f};
  float clickDecaySeconds{.009f};
  float clickTiltDb{3.f};
  float lowCutHz{18.f};
  float highCutHz{15500.f};
  float outputGain{.7f};
};

struct CompactKickParameters {
  CorrelatedFmBurstParameters primary{};
  CorrelatedFmBurstParameters secondary{};
  EnvelopedNoiseBurstParameters click{};
  ObservationModel<1>::Parameters observation{};
  CompactKickRouting routing{};
  float primaryLevel{1.f};
  float secondaryLevel{.32f};
  float clickLevel{.16f};
  float outputGain{.7f};
};

CompactKickParameters DefaultCompactKickParameters(
    const CompactKickControls &controls = {}) noexcept;

} // namespace tfdsp::percussion
