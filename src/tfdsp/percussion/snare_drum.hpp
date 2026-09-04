#pragma once

#include "membrane_drum.hpp"
#include "observation_equalizer.hpp"
#include "observation_model.hpp"
#include "wire_rack.hpp"

#include <array>
#include <cstddef>

namespace tfdsp::percussion {

enum class SnareDrumRoute : std::size_t {
  ContactToDirect, ContactToBody, FmToDirect, FmToBody,
  BodyToObservation, BodyToWires, WiresToObservation, Count
};

struct SnareDrumRouting {
  static constexpr std::size_t Count =
      static_cast<std::size_t>(SnareDrumRoute::Count);
  bool Enabled(SnareDrumRoute route) const noexcept;
  void SetEnabled(std::size_t index, bool value) noexcept;
  std::array<bool, Count> enabled{true, true, true, true, true, true, true};
};

struct SnareDrumParameters {
  MembraneDrumParameters membrane{};
  WireRackParameters wires{};
  ObservationModel<3>::Parameters observation{};
  ObservationEqualizerParameters equalizer{};
  SnareDrumRouting routing{};
  float outputGain{.2f};
};

struct SnareDrumPreparedParameters {
  MembraneDrumPreparedParameters membrane{};
  WireRackPreparedParameters wires{};
  ObservationModel<3>::Parameters observation{};
  ObservationEqualizerParameters equalizer{};
  SnareDrumRouting routing{};
  float outputGain{.2f};
  float sampleRate{48000.f};
};

struct SnareDrumFrame {
  float direct{};
  float body{};
  float wires{};
  float output{};
};

SnareDrumParameters DefaultSnareDrumParameters() noexcept;
SnareDrumPreparedParameters PrepareSnareDrumParameters(
    float sampleRate, const SnareDrumParameters &parameters);

class SnareDrum {
public:
  void Prepare(float sampleRate, const SnareDrumParameters &parameters);
  void Prepare(const SnareDrumPreparedParameters &prepared);
  void Reset() noexcept;
  void Trigger(const MembraneDrumHit &hit) noexcept;
  SnareDrumFrame ProcessFrame() noexcept;
  float Process() noexcept;
  float MembraneEnergy() const noexcept { return membrane_.ModalEnergy(); }
  float WireEnergy() const noexcept { return wires_.StoredEnergy(); }

private:
  MembraneDrum membrane_{};
  WireRack wires_{};
  ObservationModel<3> observation_{};
  ObservationEqualizer equalizer_{};
  SnareDrumRouting routing_{};
  float outputGain_{.2f};
};

} // namespace tfdsp::percussion
