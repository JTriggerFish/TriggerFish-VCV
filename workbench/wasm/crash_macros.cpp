#include "crash_macros.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <utility>

namespace tfworkbench {
namespace {

using tfdsp::percussion::CrashCymbalFitParameters;
constexpr float PiOverTwo = 1.57079632679489661923f;
constexpr float DefaultImpactWidth = .65f;

constexpr std::size_t Index(const CrashMacro macro) noexcept {
  return static_cast<std::size_t>(macro);
}

CrashMacroDescriptor Linear(std::string key, std::string name,
                            std::string unit, const float minimum,
                            const float maximum, const float value) {
  return {std::move(key), std::move(name), std::move(unit), minimum, maximum,
          value, CrashMacroScale::Linear};
}

CrashMacroDescriptor Logarithmic(std::string key, std::string name,
                                 std::string unit, const float minimum,
                                 const float maximum, const float value) {
  auto result = Linear(std::move(key), std::move(name), std::move(unit),
                       minimum, maximum, value);
  result.scale = CrashMacroScale::Logarithmic;
  return result;
}

std::array<CrashMacroDescriptor, CrashMacroCount> BuildDescriptors() {
  std::array<CrashMacroDescriptor, CrashMacroCount> result{};
  const CrashCymbalFitParameters fit = MetallicWorkbenchBaseFit();
  const auto set = [&](const CrashMacro macro, CrashMacroDescriptor descriptor) {
    result[Index(macro)] = std::move(descriptor);
  };
  set(CrashMacro::ModelLevelDb,
      Linear("model_level_db", "Model level", "dB", -60.f, 0.f, -6.f));
  set(CrashMacro::ImpactToneNoise,
      Linear("impact_tone_noise", "Impact: ping to noise", "", 0.f, 1.f, .9f));
  set(CrashMacro::ImpactWidth,
      Logarithmic("impact_width", "Contact width", "x", .25f, 4.f,
                  DefaultImpactWidth));
  set(CrashMacro::BloomRate,
      Linear("bloom_rate", "Upward cascade rate", "oct/s", 0.f, 16.f,
             fit.bloomRateOctavesPerSecond));
  set(CrashMacro::BloomEnergyAcceleration,
      Linear("bloom_energy_acceleration", "Energy acceleration", "", 0.f,
             1.f, fit.bloomEnergyAcceleration));
  set(CrashMacro::BloomPhaseDiffusion,
      Linear("bloom_phase_diffusion", "Transfer diffusion", "", 0.f, 1.f,
             fit.bloomPhaseDiffusion));
  set(CrashMacro::BodyBrightness,
      Linear("body_brightness", "Initial excitation tilt", "dB/oct", -72.f, 24.f,
             fit.bodyTiltDbPerOctave));
  set(CrashMacro::BodyExcitationCentre,
      Logarithmic("body_excitation_centre", "Excitation centre", "Hz",
                  40.f, 15000.f, fit.bodyExcitationCentreHz));
  set(CrashMacro::FieldTurbulence,
      Linear("field_turbulence", "Turbulence", "", 0.f, 1.f,
             fit.fieldTurbulence));
  set(CrashMacro::FieldTurbulenceSlope,
      Linear("field_turbulence_slope", "Turbulence slope", "/oct", -1.f,
             1.f, fit.fieldTurbulenceSlopePerOctave));
  set(CrashMacro::FieldTurbulenceCentre,
      Logarithmic("field_turbulence_centre", "Turbulence centre", "Hz",
                  40.f, 15000.f, fit.fieldTurbulenceCentreHz));
  set(CrashMacro::FieldPacketSpread,
      Linear("field_packet_spread", "Packet spread", "ERB", 0.f, 12.f,
             fit.fieldPacketSpreadErb));
  set(CrashMacro::FieldSatelliteDensity,
      Linear("field_satellite_density", "Satellite density", "", 0.f, 1.f,
             fit.fieldSatelliteDensity));
  set(CrashMacro::FieldPhaseBandwidth,
      Linear("field_phase_bandwidth", "Phase diffusion", "ERB", 0.f, 4.f,
             fit.fieldPhaseBandwidthErb));
  set(CrashMacro::FieldExchange,
      Linear("field_exchange", "Neighbour exchange", "", 0.f, 1.f,
             fit.fieldExchange));
  set(CrashMacro::BodyExcitation,
      Logarithmic("body_excitation", "Body excitation", "x", .001f, 4.f,
                  fit.bodyExcitationGain));
  set(CrashMacro::FieldGain,
      Linear("field_gain", "Body observation level", "x", 0.f, 4.f,
             fit.fieldGain));
  set(CrashMacro::DirectGain,
      Linear("direct_gain", "Contact presence", "", 0.f, 2.f,
             fit.directGain));

  const auto radiation = [&](const CrashMacro enabled, const CrashMacro low,
                             const CrashMacro frequency, const CrashMacro gain,
                             const CrashMacro high, const std::string &prefix,
                             const bool enabledValue, const float lowValue,
                             const float frequencyValue, const float gainValue,
                             const float highValue) {
    set(enabled, {prefix + "_radiation_enabled", "Enable radiation EQ", "",
                  0.f, 1.f, enabledValue ? 1.f : 0.f,
                  CrashMacroScale::Boolean});
    set(low, Logarithmic(prefix + "_low_cut", "High-pass", "Hz", 10.f,
                         1000.f, lowValue));
    set(frequency, Logarithmic(prefix + "_colour_frequency",
                               "Colour frequency", "Hz", 100.f, 18000.f,
                               frequencyValue));
    set(gain, Linear(prefix + "_colour_gain", "Colour gain", "dB", -18.f,
                     18.f, gainValue));
    set(high, Logarithmic(prefix + "_high_cut", "Low-pass", "Hz", 1000.f,
                          22000.f, highValue));
  };
  radiation(CrashMacro::DirectRadiationEnabled, CrashMacro::DirectLowCut,
            CrashMacro::DirectColourFrequency, CrashMacro::DirectColourGain,
            CrashMacro::DirectHighCut, "direct", fit.directRadiationEnabled,
            fit.directLowCutHz, fit.directColourFrequencyHz,
            fit.directColourGainDb, fit.directHighCutHz);
  radiation(CrashMacro::BodyRadiationEnabled, CrashMacro::BodyLowCut,
            CrashMacro::BodyColourFrequency, CrashMacro::BodyColourGain,
            CrashMacro::BodyHighCut, "body", fit.bodyRadiationEnabled,
            fit.bodyLowCutHz, fit.bodyColourFrequencyHz,
            fit.bodyColourGainDb, fit.bodyHighCutHz);

  for (std::size_t interior = 0; interior < BodyDecayInteriorPointCount;
       ++interior) {
    const std::size_t point = interior + 1;
    result[Index(CrashMacro::BodyDecayFrequencyFirst) + interior] = Logarithmic(
        "body_decay_frequency_" + std::to_string(point),
        "Decay centre " + std::to_string(point + 1), "Hz", 40.f, 15000.f,
        fit.bodyDecayFrequencyHz[interior]);
  }
  for (std::size_t point = 0; point < BodyDecayCurvePointCount; ++point) {
    result[Index(CrashMacro::BodyDecaySecondsFirst) + point] = Logarithmic(
        "body_decay_seconds_" + std::to_string(point),
        "Modal T60 " + std::to_string(point + 1), "s", .02f, 30.f,
        fit.bodyDecaySeconds[point]);
  }
  for (std::size_t interior = 0; interior < BodyDecayInteriorPointCount;
       ++interior) {
    const std::size_t point = interior + 1;
    result[Index(CrashMacro::BodyDecayActiveFirst) + interior] = {
        "body_decay_active_" + std::to_string(point),
        "Enable T60 knot " + std::to_string(point + 1), "", 0.f, 1.f,
        fit.bodyDecayActive[interior] ? 1.f : 0.f, CrashMacroScale::Boolean};
  }

  for (std::size_t point = 0; point < ResolvedModePointCount; ++point) {
    result[Index(CrashMacro::ResolvedFrequencyFirst) + point] = Logarithmic(
        "resolved_frequency_" + std::to_string(point),
        "Resolved mode " + std::to_string(point + 1), "Hz", 40.f, 15000.f,
        fit.sparseFrequencyHz[point]);
    const float levelDb = 20.f * std::log10(
        std::max(fit.sparseAmplitude[point], 1.e-8f));
    result[Index(CrashMacro::ResolvedLevelFirst) + point] = Linear(
        "resolved_level_" + std::to_string(point),
        "Mode prominence " + std::to_string(point + 1), "dB", -72.f, 6.f,
        std::max(levelDb, -72.f));
    result[Index(CrashMacro::ResolvedTurbulenceFirst) + point] = Linear(
        "resolved_turbulence_" + std::to_string(point),
        "Turbulence response " + std::to_string(point + 1), "x", 0.f, 2.f,
        fit.fieldTurbulenceScale[point]);
  }
  set(CrashMacro::ImpactChirpPitch,
      Logarithmic("impact_chirp_pitch", "Ping pitch", "x", .05f, 4.f,
                  fit.contactChirpFrequencyScale));
  set(CrashMacro::ImpactNoiseTilt,
      Linear("impact_noise_tilt", "Impact noise tilt", "dB/oct", -18.f,
             18.f, fit.contactNoiseTiltDb));
  set(CrashMacro::ImpactMicroDensity,
      Logarithmic("impact_micro_density", "Micro-contact density", "x",
                  .25f, 4.f, fit.contactMicroDensityScale));
  set(CrashMacro::VelocityBrightness,
      Linear("velocity_brightness", "Velocity brightness", "dB/oct", 0.f,
             12.f, fit.velocityBrightnessDbPerOctave));
  set(CrashMacro::BodyTune,
      Logarithmic("body_tune", "Body tune", "x", .5f, 2.f,
                  fit.sparseTune));
  return result;
}

const auto Descriptors = BuildDescriptors();

float ValueAt(const CrashMacroValues &values, const std::size_t index) noexcept {
  const auto &descriptor = Descriptors[index];
  const float value = std::isfinite(values[index])
      ? values[index] : descriptor.defaultValue;
  return std::clamp(value, descriptor.minimum, descriptor.maximum);
}

float Value(const CrashMacroValues &values, const CrashMacro macro) noexcept {
  return ValueAt(values, Index(macro));
}

void ApplyResolvedPaint(CrashCymbalFitParameters &fit,
                        const CrashMacroValues &values) noexcept {
  for (std::size_t mode = 0; mode < ResolvedModePointCount; ++mode) {
    fit.sparseFrequencyHz[mode] = ValueAt(
        values, Index(CrashMacro::ResolvedFrequencyFirst) + mode);
    const float level = ValueAt(
        values, Index(CrashMacro::ResolvedLevelFirst) + mode);
    fit.sparseAmplitude[mode] = level <= -71.999f
        ? 0.f : std::pow(10.f, level / 20.f);
    fit.fieldTurbulenceScale[mode] = ValueAt(
        values, Index(CrashMacro::ResolvedTurbulenceFirst) + mode);
  }
}

void ApplyBodyDecay(CrashCymbalFitParameters &fit,
                    const CrashMacroValues &values) noexcept {
  for (std::size_t interior = 0; interior < BodyDecayInteriorPointCount;
       ++interior) {
    fit.bodyDecayFrequencyHz[interior] = ValueAt(
        values, Index(CrashMacro::BodyDecayFrequencyFirst) + interior);
    fit.bodyDecayActive[interior] = ValueAt(
        values, Index(CrashMacro::BodyDecayActiveFirst) + interior) >= .5f;
  }
  for (std::size_t point = 0; point < BodyDecayCurvePointCount; ++point) {
    fit.bodyDecaySeconds[point] = ValueAt(
        values, Index(CrashMacro::BodyDecaySecondsFirst) + point);
  }
}

} // namespace

const CrashMacroDescriptor &CrashMacroDescription(
    const std::size_t index) noexcept {
  return Descriptors[std::min(index, Descriptors.size() - 1)];
}

const CrashMacroDescriptor &ActiveCrashMacroDescription(
    const std::size_t index) noexcept {
  const std::size_t bounded = std::min(
      index, ActiveCrashMacroIndices.size() - 1);
  return Descriptors[ActiveCrashMacroIndices[bounded]];
}

CrashMacroValues DefaultCrashMacros() noexcept {
  CrashMacroValues result{};
  for (std::size_t index = 0; index < result.size(); ++index)
    result[index] = Descriptors[index].defaultValue;
  return result;
}

CrashCymbalFitParameters MetallicWorkbenchBaseFit() noexcept {
  // Instrument presets start from the documented DSP defaults. No fitted
  // crash, gong, ride, or hat state is allowed to leak into another preset.
  CrashCymbalFitParameters fit{};
  return fit;
}

CrashCymbalFitParameters ApplyCrashMacros(
    const CrashCymbalFitParameters &base,
    const CrashMacroValues &values) noexcept {
  auto fit = base;
  fit.outputGain = std::pow(10.f, Value(values, CrashMacro::ModelLevelDb) / 20.f);

  const float impact = Value(values, CrashMacro::ImpactToneNoise);
  const float tonal = std::cos(PiOverTwo * impact);
  const float noisy = std::sin(PiOverTwo * impact);
  fit.contactPulseGain = .55f + .3f * (1.f - impact);
  fit.contactChirpGain = 1.5f * tonal;
  fit.contactNoiseGain = noisy;
  fit.contactMicroGain = noisy;
  const float width = Value(values, CrashMacro::ImpactWidth);
  fit.contactDurationScale = width;
  fit.contactNoiseDurationScale = width;
  fit.contactMicroDurationScale = width;
  fit.contactChirpFrequencyScale =
      Value(values, CrashMacro::ImpactChirpPitch);
  fit.contactNoiseTiltDb = Value(values, CrashMacro::ImpactNoiseTilt);
  fit.contactMicroDensityScale =
      Value(values, CrashMacro::ImpactMicroDensity);
  fit.velocityBrightnessDbPerOctave =
      Value(values, CrashMacro::VelocityBrightness);

  fit.bloomRateOctavesPerSecond = Value(values, CrashMacro::BloomRate);
  fit.bloomEnergyAcceleration = Value(
      values, CrashMacro::BloomEnergyAcceleration);
  fit.bloomPhaseDiffusion = Value(values, CrashMacro::BloomPhaseDiffusion);

  fit.bodyTiltDbPerOctave = Value(values, CrashMacro::BodyBrightness);
  fit.bodyExcitationCentreHz = Value(
      values, CrashMacro::BodyExcitationCentre);
  fit.fieldTurbulence = Value(values, CrashMacro::FieldTurbulence);
  fit.fieldTurbulenceSlopePerOctave = Value(
      values, CrashMacro::FieldTurbulenceSlope);
  fit.fieldTurbulenceCentreHz = Value(
      values, CrashMacro::FieldTurbulenceCentre);
  fit.fieldPacketSpreadErb = Value(values, CrashMacro::FieldPacketSpread);
  fit.fieldSatelliteDensity = Value(
      values, CrashMacro::FieldSatelliteDensity);
  fit.fieldPhaseBandwidthErb =
      Value(values, CrashMacro::FieldPhaseBandwidth);
  fit.fieldExchange = Value(values, CrashMacro::FieldExchange);
  fit.bodyExcitationGain = Value(values, CrashMacro::BodyExcitation);
  fit.fieldGain = Value(values, CrashMacro::FieldGain);
  fit.directGain = Value(values, CrashMacro::DirectGain);
  fit.sparseTune = Value(values, CrashMacro::BodyTune);

  ApplyResolvedPaint(fit, values);
  ApplyBodyDecay(fit, values);

  fit.directRadiationEnabled =
      Value(values, CrashMacro::DirectRadiationEnabled) >= .5f;
  fit.directLowCutHz = Value(values, CrashMacro::DirectLowCut);
  fit.directColourFrequencyHz = Value(values, CrashMacro::DirectColourFrequency);
  fit.directColourGainDb = Value(values, CrashMacro::DirectColourGain);
  fit.directHighCutHz = Value(values, CrashMacro::DirectHighCut);
  fit.bodyRadiationEnabled =
      Value(values, CrashMacro::BodyRadiationEnabled) >= .5f;
  fit.bodyLowCutHz = Value(values, CrashMacro::BodyLowCut);
  fit.bodyColourFrequencyHz = Value(values, CrashMacro::BodyColourFrequency);
  fit.bodyColourGainDb = Value(values, CrashMacro::BodyColourGain);
  fit.bodyHighCutHz = Value(values, CrashMacro::BodyHighCut);
  return fit;
}

} // namespace tfworkbench
