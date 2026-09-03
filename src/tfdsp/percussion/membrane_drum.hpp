#pragma once

#include "fixed_mixer.hpp"
#include "membrane_drum_parameters.hpp"

#include <array>
#include <cstdint>

namespace tfdsp::percussion {

struct MembraneDrumHit {
  float strength{.8f};
  float location{.5f};
  float hardness{.5f};
  float implement{1.f};
  float contactSpread{.25f};
  std::uint32_t seed{1};
};

struct MembraneDrumFrame {
  float direct{};
  float body{};
  float output{};
};

class MembraneDrum {
public:
  static constexpr std::size_t VoiceCount = 8;

  void Prepare(float sampleRate, const MembraneDrumParameters &parameters);
  void Prepare(const MembraneDrumPreparedParameters &prepared);
  void Reset() noexcept;
  void Trigger(const MembraneDrumHit &hit) noexcept;
  MembraneDrumFrame ProcessFrame() noexcept;
  float Process() noexcept;
  float StrikeEnergy() const noexcept { return strikeEnergy_.Value(); }
  float ModalEnergy() const noexcept { return membrane_.StoredEnergy(); }

private:
  struct EventVoice {
    struct Sample { float contactDirect{}, contactBody{}, fm{}; };
    void Prepare(float sampleRate);
    void Reset() noexcept;
    void Trigger(const MembraneDrumParameters &parameters,
                 const MembraneDrumHit &hit, std::uint64_t generation) noexcept;
    Sample Process() noexcept;
    bool Active() const noexcept;
    float Activity() const noexcept { return activity; }

    ContactExciter contact{};
    CorrelatedFmBurst fm{};
    float location{.5f};
    float amplitude{};
    float contactLevel{};
    float fmLevel{};
    float activity{};
    float activityRelease{};
    float stealDirect{};
    float stealBody{};
    float stealFm{};
    float stealGain{};
    float stealStep{};
    Sample last{};
    std::uint64_t generation{};
  };

  EventVoice &NextVoice() noexcept;
  void PrepareComponents(float sampleRate,
                         const MembraneDrumParameters &parameters);

  std::array<EventVoice, VoiceCount> voices_{};
  MembraneResonator<MembraneModeCount> membrane_{};
  StrikeEnergyEnvelope strikeEnergy_{};
  ObservationModel<2> observation_{};
  ObservationEqualizer equalizer_{};
  FixedMixer<2> directMixer_{};
  FixedMixer<2> bodyMixer_{};
  MembraneDrumParameters parameters_{};
  std::uint64_t generation_{};
};

} // namespace tfdsp::percussion
