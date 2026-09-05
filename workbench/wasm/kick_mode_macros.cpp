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
    constexpr const char *suffix[] = {"frequency", "level", "centre", "edge"};
    constexpr const char *labels[] = {"frequency", "prominence", "centre coupling", "edge coupling"};
    for (std::size_t mode = 0; mode < defaults.size(); ++mode) {
      const auto &source = defaults[mode];
      const float initial[] = {source.frequencyHz, source.levelDb, source.centreCoupling, source.edgeCoupling};
      for (std::size_t field = 0; field < 4; ++field) {
        const auto index = 4 * mode + field;
        std::snprintf(keys[index].data(), 48, "resonance_%s_%zu", suffix[field], mode);
        std::snprintf(names[index].data(), 48, "Mode %zu %s", mode + 1, labels[field]);
        values[index] = {keys[index].data(), names[index].data(),
          field == 0 ? "Hz" : field == 1 ? "dB" : "x",
          field == 0 ? 20.f : field == 1 ? -72.f : -1.f,
          field == 0 ? 15000.f : field == 1 ? 6.f : 1.f, initial[field],
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
