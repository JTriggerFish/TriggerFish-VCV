#pragma once

#include <cstdint>

namespace tfdsp::percussion {

// Small deterministic generator for repeatable contact detail. Structural
// tuning never depends on this stream; instruments reseed it per hit.
class DeterministicRandom {
public:
  void Seed(const std::uint32_t seed) noexcept {
    state_ = seed == 0 ? DefaultSeed : seed;
  }

  std::uint32_t NextBits() noexcept {
    std::uint32_t value = state_;
    value ^= value << 13;
    value ^= value >> 17;
    value ^= value << 5;
    state_ = value;
    return value;
  }

  float Uniform() noexcept {
    return static_cast<float>(NextBits() >> 8) * (1.f / 16777216.f);
  }

  float Bipolar() noexcept { return 2.f * Uniform() - 1.f; }

private:
  static constexpr std::uint32_t DefaultSeed = 0x6d2b79f5u;
  std::uint32_t state_{DefaultSeed};
};

} // namespace tfdsp::percussion
