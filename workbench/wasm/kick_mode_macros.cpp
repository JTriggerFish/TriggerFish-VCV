#include "kick_mode_macros.hpp"
#include "tfdsp/percussion/kick_voice_parameters.hpp"
#include <array>
#include <cstdio>

namespace tfworkbench {
namespace {
struct ModeDescriptions {
  std::array<std::array<char, 48>, KickModeParameterCount> keys{};
  std::array<std::array<char, 48>, KickModeParameterCount> names{};
  std::array<ParameterDescriptor, KickModeParameterCount> values{};

  ModeDescriptions() {
    const auto defaults = tfdsp::percussion::DefaultKickModes();
    constexpr const char *suffix[] = {"frequency", "level"};
    constexpr const char *labels[] = {"frequency", "prominence"};
    for (std::size_t mode = 0; mode < defaults.size(); ++mode) {
      const auto &source = defaults[mode];
      const float initial[] = {source.frequencyHz, source.levelDb};
      for (std::size_t field = 0; field < 2; ++field) {
        const auto index = 2 * mode + field;
        std::snprintf(keys[index].data(), 48, "resonance_%s_%zu", suffix[field], mode);
        std::snprintf(names[index].data(), 48, "Mode %zu %s", mode + 1, labels[field]);
        values[index] = {keys[index].data(), names[index].data(),
          field == 0 ? "Hz" : "dB",
          field == 0 ? 20.f : -72.f,
          field == 0 ? 15000.f : 6.f, initial[field],
          field == 0 ? ParameterScale::Logarithmic : ParameterScale::Linear};
      }
    }
  }
};
const ModeDescriptions Modes;
}
const ParameterDescriptor &KickModeDescription(std::size_t index) noexcept {
  return Modes.values[index < Modes.values.size() ? index : 0];
}
}
