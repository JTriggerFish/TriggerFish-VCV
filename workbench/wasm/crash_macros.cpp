#include "crash_macros.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <utility>

namespace tfworkbench {
namespace {

using tfdsp::percussion::CrashCymbalFitParameters;
using tfdsp::percussion::CrashModalMinimumFrequencyHz;
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
      Linear("bloom_rate", "Diffusion strength", "", 0.f, 16.f,
             fit.bloomRateOctavesPerSecond));
  set(CrashMacro::BloomEnergyAcceleration,
      Linear("bloom_energy_acceleration", "Concentration dependence", "", 0.f,
             1.f, fit.bloomEnergyAcceleration));
  set(CrashMacro::BloomEnergySensitivity,
      Linear("bloom_energy_sensitivity", "Energy sensitivity", "", 0.f,
             2.f, fit.bloomEnergySensitivity));
  set(CrashMacro::BodyBrightness,
      Linear("body_brightness", "Initial excitation tilt", "dB/oct", -72.f, 24.f,
             fit.bodyTiltDbPerOctave));
  set(CrashMacro::BodyExcitationCentre,
      Logarithmic("body_excitation_centre", "Excitation centre", "Hz",
                  CrashModalMinimumFrequencyHz, 15000.f, fit.bodyExcitationCentreHz));
  set(CrashMacro::FieldTurbulence,
      Linear("field_turbulence", "Noisiness at 1 kHz", "", 0.f, 4000.f,
             fit.fieldTurbulence));
  set(CrashMacro::FieldTurbulenceSlope,
      Linear("field_turbulence_slope", "Noisiness slope", "/oct", -1.f,
             1.f, fit.fieldTurbulenceSlopePerOctave));
  set(CrashMacro::FieldPacketSpread,
      Linear("field_packet_spread", "Packet spread", "ERB", 0.f, 12.f,
             fit.fieldPacketSpreadErb));
  set(CrashMacro::FieldSatelliteDensity,
      Linear("field_satellite_density", "Satellite density", "", 0.f, 1.f,
             fit.fieldSatelliteDensity));
  set(CrashMacro::FieldPhaseBandwidth,
      Linear("field_phase_bandwidth", "Phase blur", "ERB", 0.f, 4.f,
             fit.fieldPhaseBandwidthErb));
  set(CrashMacro::FieldPhaseTilt,
      Linear("field_phase_tilt", "Blur tilt (1 kHz pivot)", "oct/oct", -2.f, 2.f,
             fit.fieldPhaseTilt));
  set(CrashMacro::FieldDistribution,
      Linear("field_distribution", "Ring character", "", 0.f, 4.f,
             float(fit.fieldDistribution)));
  set(CrashMacro::FieldDoubletSplit,
      Linear("field_doublet_split", "Beat rate", "Hz", 0.f, 80.f,
                  fit.fieldDoubletSplitHz));
  set(CrashMacro::FieldBeatDepth,
      Linear("field_beat_depth", "Beat depth", "", 0.f, 1.f, fit.fieldBeatDepth));
  set(CrashMacro::FieldBeatRateTilt,
      Linear("field_beat_rate_tilt", "Beat rate tilt", "oct/oct", -1.f, 1.f,
             fit.fieldBeatRateTilt));
  set(CrashMacro::FieldWanderHz,
      Linear("field_wander_hz", "Pitch wander", "Hz",
             0.f, 100.f, fit.fieldWanderDepthHz));
  set(CrashMacro::FieldWanderRate,
      Logarithmic("field_wander_rate", "Wander speed", "changes/s",
                  .1f, 40.f, fit.fieldWanderKnotsPerSecond));
  set(CrashMacro::FieldMotionDepth,
      Linear("field_motion_depth", "Ridge movement", "rad", 0.f, 3.f, 0.f));
  set(CrashMacro::FieldMotionRate,
      Logarithmic("field_motion_rate", "Movement speed", "changes/s", .1f, 200.f, 40.f));
  set(CrashMacro::FieldMotionSharing,
      Linear("field_motion_sharing", "Packet sharing", "", 0.f, 1.f, .5f));
  set(CrashMacro::BodyExcitation,
      Logarithmic("body_excitation", "Body excitation", "x", .001f, 4.f,
                  fit.bodyExcitationGain));
  set(CrashMacro::FieldGain,
      Linear("field_gain", "Body observation level", "x", 0.f, 4.f,
             fit.fieldGain));
  set(CrashMacro::DirectGain,
      Linear("direct_gain", "Contact presence", "", 0.f, 2.f,
             fit.directGain));

  set(CrashMacro::OutputEqEnabled,
      {"output_eq_enabled", "Enable final EQ", "", 0.f, 1.f,
       fit.outputEqEnabled ? 1.f : 0.f, CrashMacroScale::Boolean});
  set(CrashMacro::OutputLowCut,
      Logarithmic("output_low_cut", "High-pass", "Hz", 10.f, 1000.f,
                  fit.outputLowCutHz));
  set(CrashMacro::OutputColourFrequency,
      Logarithmic("output_colour_frequency", "Colour frequency", "Hz",
                  100.f, 18000.f, fit.outputColourFrequencyHz));
  set(CrashMacro::OutputColourGain,
      Linear("output_colour_gain", "Colour gain", "dB", -18.f, 18.f,
             fit.outputColourGainDb));
  set(CrashMacro::OutputHighCut,
      Logarithmic("output_high_cut", "Low-pass", "Hz", 1000.f, 22000.f,
                  fit.outputHighCutHz));

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
    result[Index(CrashMacro::ResolvedAllocationFirst) + point] = Linear(
        "resolved_allocation_" + std::to_string(point),
        "Sideband allocation " + std::to_string(point + 1), "x", 0.f, 4.f,
        fit.fieldAllocationWeight[point]);
    result[Index(CrashMacro::ResolvedFrequencyFirst) + point] = Logarithmic(
        "resolved_frequency_" + std::to_string(point),
        "Resolved mode " + std::to_string(point + 1), "Hz", CrashModalMinimumFrequencyHz, 15000.f,
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
    fit.fieldAllocationWeight[mode] = ValueAt(
        values, Index(CrashMacro::ResolvedAllocationFirst) + mode);
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
  fit.fieldDistribution = tfdsp::percussion::ModalPacketDistribution::PairedRing;
  fit.bloomSpectralDiffusion = true;
  fit.fieldRelaxedTurbulence = true;
  fit.bloomPhaseDiffusion = 0.f;
  fit.fieldExchange = 0.f;
  fit.bloomEnergyAcceleration = 1.f;
  fit.bloomEnergySensitivity = 2.f;
  return fit;
}

CrashCymbalFitParameters ApplyCrashMacros(
    const CrashCymbalFitParameters &base,
    const CrashMacroValues &values) noexcept {
  auto fit = base;
  // This recipe has one transfer law and no secondary random energy exchange.
  fit.bloomSpectralDiffusion = true;
  fit.fieldRelaxedTurbulence = true;
  fit.bloomPhaseDiffusion = 0.f;
  fit.fieldExchange = 0.f;
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
  fit.bloomEnergySensitivity = Value(values, CrashMacro::BloomEnergySensitivity);

  fit.bodyTiltDbPerOctave = Value(values, CrashMacro::BodyBrightness);
  fit.bodyExcitationCentreHz = Value(
      values, CrashMacro::BodyExcitationCentre);
  fit.fieldTurbulence = Value(values, CrashMacro::FieldTurbulence);
  fit.fieldTurbulenceSlopePerOctave = Value(
      values, CrashMacro::FieldTurbulenceSlope);
  fit.fieldPacketSpreadErb = Value(values, CrashMacro::FieldPacketSpread);
  fit.fieldSatelliteDensity = Value(
      values, CrashMacro::FieldSatelliteDensity);
  fit.fieldPhaseBandwidthErb =
      Value(values, CrashMacro::FieldPhaseBandwidth);
  fit.fieldPhaseTilt = Value(values, CrashMacro::FieldPhaseTilt);
  fit.fieldDistribution = static_cast<tfdsp::percussion::ModalPacketDistribution>(
      int(std::round(Value(values, CrashMacro::FieldDistribution))));
  fit.fieldDoubletSplitHz = Value(values, CrashMacro::FieldDoubletSplit);
  fit.fieldBeatDepth = Value(values, CrashMacro::FieldBeatDepth);
  fit.fieldBeatRateTilt = Value(values, CrashMacro::FieldBeatRateTilt);
  fit.fieldWanderDepthHz = Value(values, CrashMacro::FieldWanderHz);
  fit.fieldWanderKnotsPerSecond = Value(values, CrashMacro::FieldWanderRate);
  fit.fieldMotion = {Value(values, CrashMacro::FieldMotionDepth),
                    Value(values, CrashMacro::FieldMotionRate),
                    Value(values, CrashMacro::FieldMotionSharing)};
  fit.bodyExcitationGain = Value(values, CrashMacro::BodyExcitation);
  fit.fieldGain = Value(values, CrashMacro::FieldGain);
  fit.directGain = Value(values, CrashMacro::DirectGain);
  fit.sparseTune = Value(values, CrashMacro::BodyTune);

  ApplyResolvedPaint(fit, values);
  ApplyBodyDecay(fit, values);

  fit.outputEqEnabled =
      Value(values, CrashMacro::OutputEqEnabled) >= .5f;
  fit.outputLowCutHz = Value(values, CrashMacro::OutputLowCut);
  fit.outputColourFrequencyHz = Value(values, CrashMacro::OutputColourFrequency);
  fit.outputColourGainDb = Value(values, CrashMacro::OutputColourGain);
  fit.outputHighCutHz = Value(values, CrashMacro::OutputHighCut);
  return fit;
}

} // namespace tfworkbench
