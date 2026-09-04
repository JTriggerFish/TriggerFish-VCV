#include "membrane_macros.hpp"

#include <array>
#include <cmath>

namespace tfworkbench {
namespace {

using Scale = ParameterScale;

const std::array<ParameterDescriptor, MembraneParameterCount> Descriptors{{
    {"model_level_db", "Model level", "dB", -60.f, 0.f, -22.f},
    {"fundamental_hz", "Fundamental", "Hz", 25.f, 500.f, 105.f,
     Scale::Logarithmic},
    {"decay_seconds", "Body decay", "s", .03f, 8.f, 1.15f, Scale::Logarithmic},
    {"decay_tilt", "Decay tilt", "", -1.f, 1.f, .55f},
    {"inharmonicity", "Membrane character", "", 0.f, 1.f, .35f},
    {"body_brightness", "Body brightness", "", 0.f, 1.f, .55f},
    {"tension_octaves", "Energy pitch lift", "oct", -.25f, .6f, .11f},
    {"tension_decay_seconds", "Tension recovery", "s", .005f, 2.f, .13f,
     Scale::Logarithmic},
    {"contact_level", "Contact level", "x", 0.f, 3.f, .7f},
    {"contact_duration_seconds", "Contact width", "s", .0002f, .08f, .004f,
     Scale::Logarithmic},
    {"contact_brightness", "Contact brightness", "", 0.f, 1.f, .58f},
    {"direct_velocity_exponent", "Direct velocity curve", "exp", 1.f, 3.f,
     1.f},
    {"body_velocity_exponent", "Body velocity curve", "exp", 1.f, 3.f, 1.f},
    {"fm_level", "FM supplement", "x", 0.f, 2.f, .18f},
    {"fm_depth_hz", "FM depth", "Hz", 0.f, 8000.f, 260.f},
    {"fm_decay_seconds", "FM decay", "s", .003f, 1.f, .07f, Scale::Logarithmic},
    {"pitch_drop_octaves", "FM pitch drop", "oct", 0.f, 3.f, .28f},
    {"direct_level", "Direct level", "x", 0.f, 3.f, .3f},
    {"body_level", "Body level", "x", 0.f, 3.f, 1.f},
    {"direct_delay_ms", "Direct delay", "ms", 0.f, 10.f, 0.f},
    {"equalizer_mode", "Output EQ", "mode", 0.f, 2.f, 1.f, Scale::Choice},
    {"low_cut_hz", "High-pass", "Hz", 5.f, 1000.f, 24.f, Scale::Logarithmic},
    {"high_cut_hz", "Low-pass", "Hz", 500.f, 22000.f, 18000.f,
     Scale::Logarithmic},
    {"colour_frequency_hz", "Radiation colour", "Hz", 40.f, 20000.f, 2800.f,
     Scale::Logarithmic},
    {"colour_gain_db", "Radiation colour gain", "dB", -24.f, 24.f, 0.f},
    {"band_1_frequency_hz", "Band 1 frequency", "Hz", 30.f, 500.f, 90.f,
     Scale::Logarithmic},
    {"band_1_gain_db", "Band 1 gain", "dB", -24.f, 24.f, 0.f},
    {"band_2_frequency_hz", "Band 2 frequency", "Hz", 100.f, 2000.f, 350.f,
     Scale::Logarithmic},
    {"band_2_gain_db", "Band 2 gain", "dB", -24.f, 24.f, 0.f},
    {"band_3_frequency_hz", "Band 3 frequency", "Hz", 400.f, 8000.f, 1800.f,
     Scale::Logarithmic},
    {"band_3_gain_db", "Band 3 gain", "dB", -24.f, 24.f, 0.f},
    {"band_4_frequency_hz", "Band 4 frequency", "Hz", 1500.f, 20000.f, 7500.f,
     Scale::Logarithmic},
    {"band_4_gain_db", "Band 4 gain", "dB", -24.f, 24.f, 0.f},
}};

std::size_t Index(const MembraneParameter parameter) noexcept {
  return static_cast<std::size_t>(parameter);
}

} // namespace

const ParameterDescriptor &
MembraneParameterDescription(const std::size_t index) noexcept {
  return Descriptors[index < Descriptors.size() ? index : 0];
}

MembraneParameterValues DefaultMembraneParameters() noexcept {
  MembraneParameterValues result{};
  for (std::size_t index = 0; index < result.size(); ++index)
    result[index] = Descriptors[index].defaultValue;
  return result;
}

tfdsp::percussion::MembraneDrumParameters
ApplyMembraneParameters(const MembraneParameterValues &values) noexcept {
  using P = MembraneParameter;
  tfdsp::percussion::MembraneDrumControls controls;
  controls.fundamentalHz = values[Index(P::FundamentalHz)];
  controls.decaySeconds = values[Index(P::DecaySeconds)];
  controls.decayTilt = values[Index(P::DecayTilt)];
  controls.inharmonicity = values[Index(P::Inharmonicity)];
  controls.bodyBrightness = values[Index(P::BodyBrightness)];
  controls.tensionOctaves = values[Index(P::TensionOctaves)];
  controls.tensionDecaySeconds = values[Index(P::TensionDecaySeconds)];
  controls.contactLevel = values[Index(P::ContactLevel)];
  controls.contactDurationSeconds = values[Index(P::ContactDurationSeconds)];
  controls.contactBrightness = values[Index(P::ContactBrightness)];
  controls.fmLevel = values[Index(P::FmLevel)];
  controls.fmDepthHz = values[Index(P::FmDepthHz)];
  controls.fmDecaySeconds = values[Index(P::FmDecaySeconds)];
  controls.pitchDropOctaves = values[Index(P::PitchDropOctaves)];
  controls.directLevel = values[Index(P::DirectLevel)];
  controls.bodyLevel = values[Index(P::BodyLevel)];
  controls.directDelaySeconds = .001f * values[Index(P::DirectDelayMs)];
  controls.equalizerMode =
      static_cast<tfdsp::percussion::ObservationEqualizerMode>(
          static_cast<int>(std::lround(values[Index(P::EqualizerMode)])));
  controls.lowCutHz = values[Index(P::LowCutHz)];
  controls.highCutHz = values[Index(P::HighCutHz)];
  controls.colourFrequencyHz = values[Index(P::ColourFrequencyHz)];
  controls.colourGainDb = values[Index(P::ColourGainDb)];
  controls.outputGain = std::pow(10.f, values[Index(P::ModelLevelDb)] / 20.f);
  auto result = tfdsp::percussion::DefaultMembraneDrumParameters(controls);
  result.directVelocityExponent = values[Index(P::DirectVelocityExponent)];
  result.bodyVelocityExponent = values[Index(P::BodyVelocityExponent)];
  constexpr std::array<P, 4> frequencies{
      P::Band1FrequencyHz, P::Band2FrequencyHz, P::Band3FrequencyHz,
      P::Band4FrequencyHz};
  constexpr std::array<P, 4> gains{P::Band1GainDb, P::Band2GainDb,
                                   P::Band3GainDb, P::Band4GainDb};
  for (std::size_t band = 0; band < 4; ++band) {
    result.equalizer.bands[band].frequencyHz = values[Index(frequencies[band])];
    result.equalizer.bands[band].gainDb = values[Index(gains[band])];
  }
  return result;
}

} // namespace tfworkbench
