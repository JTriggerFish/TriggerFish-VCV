#pragma once

#include "crash_cymbal_parameters.hpp"
#include "frequency_shifter.hpp"
#include "passive_constraint.hpp"
#include "resonator_submix.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

namespace tfdsp::percussion {

struct CrashCymbalHit {
  float strength{.8f};
  float location{1.f};
  float hardness{.65f};
  std::uint32_t seed{1};
};

// Mono crash body. Contact is the only dry source; its body port drives one
// dispersion loop, which drives one wet-only coupled resonator network.
class CrashCymbal {
public:
  void Prepare(float sampleRate, const CrashCymbalParameters &parameters);
  void Reset() noexcept;
  void Trigger(const CrashCymbalHit &hit) noexcept;
  float Process() noexcept;
  void SetMute(float amount) noexcept;

  float MinimumBodyDelaySamples() const noexcept;

private:
  using Projection = std::array<float, CrashResonatorCount>;

  ContactExciterParameters ContactParameters(
      const CrashCymbalHit &hit) const noexcept;
  Projection LocationProjection(float location) const noexcept;
  void SetBloomDrive(float strength) noexcept;

  ContactExciter contact_{};
  DispersionLoop dispersion_{};
  CrashResonators resonators_{};
  ResonatorSubmix<CrashResonatorCount, CrashResonatorBusCount> submix_{};
  std::array<FrequencyShifter, CrashResonatorBusCount> shifters_{};
  ObservationModel<2> observation_{};
  DynamicLossController constraint_{};
  CrashCymbalParameters parameters_{};
  Projection projection_{};
  float sampleRate_{48000.f};
};

} // namespace tfdsp::percussion
