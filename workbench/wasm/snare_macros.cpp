#include "snare_macros.hpp"

#include <algorithm>
#include <array>

namespace tfworkbench {
namespace {

using Scale = ParameterScale;

const std::array<ParameterDescriptor,
    SnareParameterCount - MembraneParameterCount> WireDescriptors{{
    {"wire_level", "Wire level", "x", 0.f, 3.f, .35f},
    {"wire_sensitivity", "Wire sensitivity", "x", 0.f, 32.f, 9.f},
    {"wire_threshold", "Wire threshold", "", 0.f, .08f, .012f},
    {"wire_attack_seconds", "Wire engagement", "s", .0001f, .05f, .006f,
     Scale::Logarithmic},
    {"wire_release_seconds", "Contact release", "s", .0005f, .2f, .01f,
     Scale::Logarithmic},
    {"wire_minimum_hz", "Wire low edge", "Hz", 80.f, 6000.f, 520.f,
     Scale::Logarithmic},
    {"wire_maximum_hz", "Wire high edge", "Hz", 2000.f, 22000.f, 9000.f,
     Scale::Logarithmic},
    {"wire_decay_seconds", "Wire decay", "s", .008f, 2.f, .06f,
     Scale::Logarithmic},
    {"wire_decay_tilt", "Wire decay tilt", "", -1.f, 1.f, .78f},
    {"wire_density", "Wire density", "", 0.f, 1.f, .9f},
    {"wire_brightness", "Wire brightness", "", 0.f, 1.f, .05f},
    {"wire_noise_mix", "Wire noise", "x", 0.f, 2.f, .42f},
    {"wire_modal_mix", "Wire modes", "x", 0.f, 2.f, .7f},
    {"ring_frequency_hz", "Persistent ring", "Hz", 100.f, 2000.f, 675.f,
     Scale::Logarithmic},
    {"ring_decay_seconds", "Ring decay", "s", .05f, 3.f, 1.8f,
     Scale::Logarithmic},
    {"ring_level", "Ring level", "x", 0.f, 3.f, .3f},
}};

std::size_t Index(const SnareParameter parameter) noexcept {
  return static_cast<std::size_t>(parameter);
}

void Set(SnareParameterValues &values, const MembraneParameter parameter,
         const float value) noexcept {
  values[static_cast<std::size_t>(parameter)] = value;
}

SnareParameterValues MakeDefaultValues() noexcept {
  SnareParameterValues result{};
  const auto membrane = DefaultMembraneParameters();
  std::copy(membrane.begin(), membrane.end(), result.begin());
  for (std::size_t index = MembraneParameterCount;
       index < result.size(); ++index)
    result[index] = WireDescriptors[index - MembraneParameterCount].defaultValue;
  using P = MembraneParameter;
  Set(result, P::FundamentalHz, 185.f);
  Set(result, P::DecaySeconds, .12f);
  Set(result, P::DecayTilt, .45f);
  Set(result, P::Inharmonicity, .72f);
  Set(result, P::BodyBrightness, .68f);
  Set(result, P::TensionOctaves, .08f);
  Set(result, P::TensionDecaySeconds, .09f);
  Set(result, P::ContactLevel, .78f);
  Set(result, P::ContactDurationSeconds, .0022f);
  Set(result, P::ContactBrightness, .38f);
  Set(result, P::FmLevel, .08f);
  Set(result, P::FmDepthHz, 180.f);
  Set(result, P::FmDecaySeconds, .035f);
  Set(result, P::PitchDropOctaves, .12f);
  Set(result, P::DirectLevel, .06f);
  Set(result, P::BodyLevel, .95f);
  Set(result, P::LowCutHz, 38.f);
  Set(result, P::HighCutHz, 10000.f);
  Set(result, P::ColourFrequencyHz, 210.f);
  Set(result, P::ColourGainDb, 0.f);
  return result;
}

const std::array<ParameterDescriptor, SnareParameterCount> &Descriptors() {
  static const auto result = [] {
    std::array<ParameterDescriptor, SnareParameterCount> descriptions{};
    for (std::size_t index = 0; index < MembraneParameterCount; ++index)
      descriptions[index] = MembraneParameterDescription(index);
    std::copy(WireDescriptors.begin(), WireDescriptors.end(),
              descriptions.begin() + MembraneParameterCount);
    const auto defaults = MakeDefaultValues();
    for (std::size_t index = 0; index < descriptions.size(); ++index)
      descriptions[index].defaultValue = defaults[index];
    return descriptions;
  }();
  return result;
}

} // namespace

const ParameterDescriptor &SnareParameterDescription(
    const std::size_t index) noexcept {
  const auto &descriptors = Descriptors();
  return descriptors[index < descriptors.size() ? index : 0];
}

SnareParameterValues DefaultSnareParameters() noexcept {
  return MakeDefaultValues();
}

tfdsp::percussion::SnareDrumParameters ApplySnareParameters(
    const SnareParameterValues &values) noexcept {
  MembraneParameterValues membraneValues{};
  std::copy_n(values.begin(), membraneValues.size(), membraneValues.begin());
  auto membrane = ApplyMembraneParameters(membraneValues);
  auto result = tfdsp::percussion::DefaultSnareDrumParameters();
  result.membrane = membrane;
  result.observation[0] = membrane.observation[0];
  result.observation[1] = membrane.observation[1];
  result.observation[2].gain = values[Index(SnareParameter::WireLevel)];
  result.observation[2].radiationEnabled = false;
  result.equalizer = membrane.equalizer;
  result.outputGain = membrane.outputGain;
  result.membrane.observation[0].gain = 1.f;
  result.membrane.observation[1].gain = 1.f;
  result.membrane.equalizer.mode =
      tfdsp::percussion::ObservationEqualizerMode::Bypass;
  result.membrane.outputGain = 1.f;
  result.wires.sensitivity = values[Index(SnareParameter::WireSensitivity)];
  result.wires.threshold = values[Index(SnareParameter::WireThreshold)];
  result.wires.attackSeconds = values[Index(SnareParameter::WireAttackSeconds)];
  result.wires.releaseSeconds = values[Index(SnareParameter::WireReleaseSeconds)];
  result.wires.minimumFrequencyHz = values[Index(SnareParameter::WireMinimumHz)];
  result.wires.maximumFrequencyHz = values[Index(SnareParameter::WireMaximumHz)];
  result.wires.decaySeconds = values[Index(SnareParameter::WireDecaySeconds)];
  result.wires.decayTilt = values[Index(SnareParameter::WireDecayTilt)];
  result.wires.density = values[Index(SnareParameter::WireDensity)];
  result.wires.brightness = values[Index(SnareParameter::WireBrightness)];
  result.wires.noiseMix = values[Index(SnareParameter::WireNoiseMix)];
  result.wires.modalMix = values[Index(SnareParameter::WireModalMix)];
  constexpr std::size_t RingMode = 5;
  result.membrane.membrane[RingMode].frequencyHz =
      values[Index(SnareParameter::RingFrequencyHz)];
  result.membrane.membrane[RingMode].decaySeconds =
      values[Index(SnareParameter::RingDecaySeconds)];
  result.membrane.membrane[RingMode].outputGain =
      values[Index(SnareParameter::RingLevel)];
  return result;
}

} // namespace tfworkbench
