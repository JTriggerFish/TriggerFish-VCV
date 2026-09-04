#pragma once

#include "deterministic_random.hpp"
#include "spectral_tilt_filter.hpp"
#include "tfdsp/finite_audio.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace tfdsp::percussion {

inline constexpr std::size_t WireRackModeCount = 48;

struct WireRackParameters {
  float sensitivity{9.f};
  float threshold{.004f};
  float motionHighpassHz{140.f};
  float attackSeconds{.002f};
  float releaseSeconds{.018f};
  float minimumFrequencyHz{900.f};
  float maximumFrequencyHz{15500.f};
  float decaySeconds{.16f};
  float decayTilt{.7f};
  float density{.8f};
  float brightness{.62f};
  float noiseMix{.6f};
  float modalMix{.75f};
  float outputGain{.42f};
  float maximumModalEnergy{16.f};
  std::uint32_t seed{0x57495245u};
};

struct WireRackPreparedParameters {
  std::array<float, WireRackModeCount> cosine{};
  std::array<float, WireRackModeCount> sine{};
  std::array<float, WireRackModeCount> radius{};
  std::array<float, WireRackModeCount> inputPhaseCosine{};
  std::array<float, WireRackModeCount> inputPhaseSine{};
  std::array<float, WireRackModeCount> modeOutputGain{};
  float sampleRate{48000.f};
  float motionCoefficient{};
  float attackCoefficient{};
  float releaseCoefficient{};
  float sensitivity{9.f};
  float threshold{.004f};
  float noiseTiltDb{};
  float noiseMix{.6f};
  float modalMix{.75f};
  float outputLevel{.42f};
  float maximumModalEnergy{16.f};
  std::size_t activeModeCount{WireRackModeCount};
  std::uint32_t seed{0x57495245u};
};

WireRackPreparedParameters PrepareWireRackParameters(
    float sampleRate, const WireRackParameters &parameters);

// A compact snare-wire interaction driven continuously by body motion. A
// motion high-pass rejects static displacement, while a contact
// envelope drives correlated noise and a dense, normalized wire-mode bank.
class WireRack {
public:
  void Prepare(float sampleRate, const WireRackParameters &parameters);
  void Prepare(const WireRackPreparedParameters &parameters);
  void Reset() noexcept;
  void Seed(std::uint32_t seed) noexcept;
  float Process(float bodyMotion) noexcept;
  float StoredEnergy() const noexcept;
  float ContactAmount() const noexcept { return contactEnvelope_; }

private:
  float AvailableDriveScale(float baseEnergy, float crossEnergy,
                            float driveEnergy) const noexcept;

  std::array<float, WireRackModeCount> real_{};
  std::array<float, WireRackModeCount> imaginary_{};
  WireRackPreparedParameters parameters_{};
  DeterministicRandom random_{};
  SpectralTiltFilter tilt_{};
  float motionLowpass_{};
  float contactEnvelope_{};
};

} // namespace tfdsp::percussion
