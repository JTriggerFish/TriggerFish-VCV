#pragma once

#include "compact_kick_parameters.hpp"
#include "fixed_mixer.hpp"
#include "tfdsp/finite_audio.hpp"

#include <array>
#include <cstdint>

namespace tfdsp::percussion {

struct CompactKickHit {
  float strength{.8f};
  float hardness{.5f};
  std::uint32_t seed{1};
};

struct CompactKickFrame {
  float primary{};
  float secondary{};
  float click{};
  float output{};
};

class CompactKick {
public:
  static constexpr std::size_t VoiceCount = 8;

  void Prepare(float sampleRate, const CompactKickParameters &parameters);
  void Reset() noexcept;
  void Trigger(const CompactKickHit &hit) noexcept;
  CompactKickFrame ProcessFrame() noexcept;
  float Process() noexcept;

private:
  struct EventVoice {
    void Prepare(float sampleRate);
    void Reset() noexcept;
    void Trigger(const CompactKickParameters &parameters,
                 const CompactKickHit &hit, std::uint64_t generation) noexcept;
    CompactKickFrame Process() noexcept;
    bool Active() const noexcept;
    float Activity() const noexcept;

    CorrelatedFmBurst primary{};
    CorrelatedFmBurst secondary{};
    EnvelopedNoiseBurst click{};
    float amplitude{};
    float primaryLevel{1.f};
    float secondaryLevel{};
    float clickLevel{};
    float lastPrimary{};
    float lastSecondary{};
    float lastClick{};
    float stealPrimary{};
    float stealSecondary{};
    float stealClick{};
    float stealGain{};
    float stealStep{1.f};
    float activity{};
    float activityRelease{};
    std::uint64_t generation{};
  };

  EventVoice &NextVoice() noexcept;

  std::array<EventVoice, VoiceCount> voices_{};
  ObservationModel<1> observation_{};
  FixedMixer<3> sourceMixer_{};
  CompactKickParameters parameters_{};
  std::uint64_t generation_{};
};

} // namespace tfdsp::percussion
