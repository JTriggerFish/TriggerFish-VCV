#include "crash_macros.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <utility>

namespace tfworkbench {
namespace {

using tfdsp::percussion::CrashCymbalFitParameters;
using tfdsp::percussion::ErbRate;
using tfdsp::percussion::InverseErbRate;
constexpr float PiOverTwo = 1.57079632679489661923f;
constexpr float DefaultBloomDevelopment = .44875992f;
constexpr std::array<float, DenseWashCurvePointCount>
    StartingDenseWashFrequencyHz{
        120.f, 500.f, 1500.f, 3000.f, 6000.f, 10000.f, 16000.f, 21000.f};
constexpr std::array<float, DenseWashCurvePointCount> StartingDenseWashLevelDb{
    -7.f, -4.5f, -2.5f, 0.f, 2.5f, 2.f, -4.f, -7.f};
constexpr std::array<float, BodyDecayCurvePointCount> StartingDecayFrequencyHz{
    150.f, 500.f, 1500.f, 6000.f, 16000.f};
constexpr std::array<float, BodyDecayCurvePointCount> StartingDecaySeconds{
    3.5f, 2.8f, 2.5f, 4.5f, .35f};

struct CurvePoint {
  float frequency;
  float level;
};

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
  const CrashCymbalFitParameters fit;
  const auto set = [&](const CrashMacro macro, CrashMacroDescriptor descriptor) {
    result[Index(macro)] = std::move(descriptor);
  };
  set(CrashMacro::ModelLevelDb,
      Linear("model_level_db", "Model level", "dB", -60.f, 12.f, -20.f));
  set(CrashMacro::ImpactToneNoise,
      Linear("impact_tone_noise", "Impact: ping to noise", "", 0.f, 1.f, .9f));
  set(CrashMacro::ImpactWidth,
      Logarithmic("impact_width", "Contact width", "x", .25f, 4.f, .65f));
  set(CrashMacro::BloomAmount,
      Linear("bloom_amount", "Bloom amount", "", 0.f, 1.f, .35f));
  set(CrashMacro::BloomDevelopment,
      Linear("bloom_development", "Bloom development", "", 0.f, 1.f,
             DefaultBloomDevelopment));
  set(CrashMacro::BodyToneWash,
      Linear("body_tone_wash", "Resolved to wash", "", 0.f, 1.f,
             1.f));
  set(CrashMacro::BodyBrightness,
      Linear("body_brightness", "Global wash tilt", "dB/oct", -8.f, 5.f,
             -1.f));
  set(CrashMacro::DirectGain,
      Linear("direct_gain", "Contact presence", "", 0.f, 2.f,
             .95f));

  set(CrashMacro::TurbulenceEnabled,
      {"turbulence_enabled", "Enable turbulent residual", "", 0.f, 1.f, 1.f,
       CrashMacroScale::Boolean});
  const float turbulenceAmount = std::cbrt(
      fit.turbulenceGain[0] * fit.turbulenceGain[1] * fit.turbulenceGain[2]);
  set(CrashMacro::TurbulenceAmount,
      Logarithmic("turbulence_amount", "Output level", "", .001f, 2.f,
                  turbulenceAmount));
  set(CrashMacro::TurbulencePersistence,
      Logarithmic("turbulence_persistence", "Persistence", "x", .25f, 4.f,
                  fit.turbulencePersistence));
  float turbulenceMeanDb = 0.f;
  std::array<float, TurbulenceCurvePointCount> turbulenceDb{};
  for (std::size_t point = 0; point < turbulenceDb.size(); ++point) {
    turbulenceDb[point] = 20.f * std::log10(
        std::max(fit.turbulenceGain[point], 1.e-6f));
    turbulenceMeanDb += turbulenceDb[point] / turbulenceDb.size();
  }
  for (std::size_t point = 0; point < TurbulenceCurvePointCount; ++point) {
    result[Index(CrashMacro::TurbulenceFrequencyFirst) + point] = Logarithmic(
        "turbulence_frequency_" + std::to_string(point),
        "Turbulence centre " + std::to_string(point + 1), "Hz",
        40.f, 20000.f, fit.turbulenceFrequencyHz[point]);
    result[Index(CrashMacro::TurbulenceLevelFirst) + point] = Linear(
        "turbulence_level_" + std::to_string(point),
        "Turbulence colour " + std::to_string(point + 1), "dB",
        -18.f, 18.f, turbulenceDb[point] - turbulenceMeanDb);
  }

  const auto radiation = [&](const CrashMacro enabled, const CrashMacro low,
                             const CrashMacro lowQ, const CrashMacro frequency,
                             const CrashMacro gain, const CrashMacro colourQ,
                             const CrashMacro high, const CrashMacro highQ,
                             const std::string &prefix, const bool enabledValue,
                             const float lowValue, const float lowQValue,
                             const float frequencyValue, const float gainValue,
                             const float colourQValue, const float highValue,
                             const float highQValue) {
    set(enabled, {prefix + "_radiation_enabled", "Enable radiation EQ", "",
                  0.f, 1.f, enabledValue ? 1.f : 0.f,
                  CrashMacroScale::Boolean});
    set(low, Logarithmic(prefix + "_low_cut", "High-pass", "Hz", 10.f,
                         1000.f, lowValue));
    set(lowQ, Logarithmic(prefix + "_low_cut_q", "High-pass Q", "", .25f,
                          4.f, lowQValue));
    set(frequency, Logarithmic(prefix + "_colour_frequency",
                               "Colour frequency", "Hz", 100.f, 18000.f,
                               frequencyValue));
    set(gain, Linear(prefix + "_colour_gain", "Colour gain", "dB", -18.f,
                     18.f, gainValue));
    set(colourQ, Logarithmic(prefix + "_colour_q", "Colour Q", "", .25f,
                             12.f, colourQValue));
    set(high, Logarithmic(prefix + "_high_cut", "Low-pass", "Hz", 1000.f,
                          22000.f, highValue));
    set(highQ, Logarithmic(prefix + "_high_cut_q", "Low-pass Q", "", .25f,
                           4.f, highQValue));
  };
  radiation(CrashMacro::DirectRadiationEnabled, CrashMacro::DirectLowCut,
            CrashMacro::DirectLowCutQ, CrashMacro::DirectColourFrequency,
            CrashMacro::DirectColourGain, CrashMacro::DirectColourQ,
            CrashMacro::DirectHighCut, CrashMacro::DirectHighCutQ, "direct",
            fit.directRadiationEnabled, fit.directLowCutHz, fit.directLowCutQ,
            7000.f, 4.f, .7f, fit.directHighCutHz, fit.directHighCutQ);
  radiation(CrashMacro::SparseRadiationEnabled, CrashMacro::SparseLowCut,
            CrashMacro::SparseLowCutQ, CrashMacro::SparseColourFrequency,
            CrashMacro::SparseColourGain, CrashMacro::SparseColourQ,
            CrashMacro::SparseHighCut, CrashMacro::SparseHighCutQ, "sparse",
            fit.sparseRadiationEnabled, fit.sparseLowCutHz, fit.sparseLowCutQ,
            fit.colourFrequencyHz, fit.colourGainDb, fit.sparseColourQ,
            fit.highCutHz, fit.sparseHighCutQ);
  radiation(CrashMacro::DenseRadiationEnabled, CrashMacro::DenseLowCut,
            CrashMacro::DenseLowCutQ, CrashMacro::DenseColourFrequency,
            CrashMacro::DenseColourGain, CrashMacro::DenseColourQ,
            CrashMacro::DenseHighCut, CrashMacro::DenseHighCutQ, "dense",
            fit.denseRadiationEnabled, fit.denseLowCutHz, fit.denseLowCutQ,
            fit.denseColourFrequencyHz, fit.denseColourGainDb, fit.denseColourQ,
            fit.denseHighCutHz, fit.denseHighCutQ);

  set(CrashMacro::DenseMinimumFrequency,
      Logarithmic("dense_minimum_frequency", "Lower frequency", "Hz", 40.f,
                  2000.f, fit.denseMinimumFrequencyHz));
  set(CrashMacro::DenseMaximumFrequency,
      Logarithmic("dense_maximum_frequency", "Upper frequency", "Hz", 4000.f,
                  22000.f, fit.denseMaximumFrequencyHz));
  set(CrashMacro::DenseModeDensity,
      Logarithmic("dense_mode_density", "Mode density", "x", 1.f / 32.f,
                  2.f, 2.f));
  set(CrashMacro::DenseSpacingJitter,
      Linear("dense_spacing_jitter", "Spacing irregularity", "", 0.f, .95f,
             fit.denseSpacingJitter));
  set(CrashMacro::DenseDecaySpread,
      Linear("dense_decay_spread", "Decay spread", "oct", 0.f, 2.f,
             fit.denseDecaySpreadOctaves));
  set(CrashMacro::DenseGainSpread,
      Linear("dense_gain_spread", "Level spread", "dB", 0.f, 18.f,
             fit.denseGainSpreadDb));

  for (std::size_t point = 0; point < DenseWashCurvePointCount; ++point) {
    result[Index(CrashMacro::DenseWashFrequencyFirst) + point] = Logarithmic(
        "dense_wash_frequency_" + std::to_string(point),
        "Wash colour centre " + std::to_string(point + 1), "Hz",
        40.f, 22000.f, StartingDenseWashFrequencyHz[point]);
    result[Index(CrashMacro::DenseWashLevelFirst) + point] = Linear(
        "dense_wash_level_" + std::to_string(point),
        "Wash colour " + std::to_string(point + 1), "dB", -24.f, 24.f,
        StartingDenseWashLevelDb[point]);
  }

  for (std::size_t point = 0; point < BodyDecayCurvePointCount; ++point) {
    result[Index(CrashMacro::BodyDecayFrequencyFirst) + point] = Logarithmic(
        "body_decay_frequency_" + std::to_string(point),
        "Decay centre " + std::to_string(point + 1), "Hz", 40.f, 20000.f,
        StartingDecayFrequencyHz[point]);
    result[Index(CrashMacro::BodyDecaySecondsFirst) + point] = Logarithmic(
        "body_decay_seconds_" + std::to_string(point),
        "Body T60 " + std::to_string(point + 1), "s", .02f, 20.f,
        StartingDecaySeconds[point]);
  }

  for (std::size_t point = 0; point < ResolvedModePointCount; ++point) {
    result[Index(CrashMacro::ResolvedFrequencyFirst) + point] = Logarithmic(
        "resolved_frequency_" + std::to_string(point),
        "Resolved mode " + std::to_string(point + 1), "Hz", 40.f, 22000.f,
        fit.sparseFrequencyHz[point]);
    result[Index(CrashMacro::ResolvedLevelFirst) + point] = Linear(
        "resolved_level_" + std::to_string(point),
        "Mode level " + std::to_string(point + 1), "dB", -24.f, 24.f,
        0.f);
  }
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

template <std::size_t Count>
std::array<CurvePoint, Count> ReadCurve(
    const CrashMacroValues &values, const CrashMacro frequencyFirst,
    const CrashMacro levelFirst) noexcept {
  std::array<CurvePoint, Count> points{};
  for (std::size_t point = 0; point < Count; ++point) {
    points[point] = {
        ValueAt(values, Index(frequencyFirst) + point),
        ValueAt(values, Index(levelFirst) + point)};
  }
  std::sort(points.begin(), points.end(), [](const auto &left,
                                              const auto &right) {
    return left.frequency < right.frequency;
  });
  return points;
}

template <std::size_t Count>
float SmoothCurve(const std::array<CurvePoint, Count> &points,
                  const float frequency) noexcept {
  if (frequency <= points.front().frequency)
    return points.front().level;
  for (std::size_t right = 1; right < Count; ++right) {
    if (frequency > points[right].frequency)
      continue;
    const float leftFrequency = std::max(points[right - 1].frequency, 1.f);
    const float rightFrequency = std::max(points[right].frequency,
                                          leftFrequency + 1.f);
    const float linear = std::clamp(
        std::log(frequency / leftFrequency) /
            std::log(rightFrequency / leftFrequency),
        0.f, 1.f);
    const float smooth = .5f - .5f * std::cos(2.f * PiOverTwo * linear);
    return points[right - 1].level + smooth *
        (points[right].level - points[right - 1].level);
  }
  return points.back().level;
}

void ApplyResolvedPaint(CrashCymbalFitParameters &fit,
                        const CrashMacroValues &values) noexcept {
  const auto points = ReadCurve<ResolvedModePointCount>(
      values, CrashMacro::ResolvedFrequencyFirst, CrashMacro::ResolvedLevelFirst);
  float meanDb = 0.f;
  for (const auto &point : points)
    meanDb += point.level / points.size();
  for (std::size_t mode = 0; mode < points.size(); ++mode) {
    fit.sparseFrequencyHz[mode] = points[mode].frequency;
    fit.sparseAmplitude[mode] *=
        std::pow(10.f, (points[mode].level - meanDb) / 20.f);
  }
}

void ApplyDenseWashPaint(CrashCymbalFitParameters &fit,
                         const CrashMacroValues &values) noexcept {
  const auto points = ReadCurve<DenseWashCurvePointCount>(
      values, CrashMacro::DenseWashFrequencyFirst,
      CrashMacro::DenseWashLevelFirst);
  const float minimumErb = ErbRate(fit.denseMinimumFrequencyHz);
  const float maximumErb = ErbRate(fit.denseMaximumFrequencyHz);
  constexpr std::size_t Count =
      tfdsp::percussion::StatisticalModalCloudParameters::GainEnvelopePointCount;
  std::array<float, Count> sampledCurve{};
  float meanDb = 0.f;
  for (std::size_t point = 0; point < Count; ++point) {
    const float position = static_cast<float>(point) /
        static_cast<float>(Count - 1);
    const float frequency = InverseErbRate(
        minimumErb + position * (maximumErb - minimumErb));
    sampledCurve[point] = SmoothCurve(points, frequency);
    meanDb += sampledCurve[point] / static_cast<float>(Count);
  }
  for (std::size_t point = 0; point < Count; ++point)
    fit.denseGainEnvelopeDb[point] += sampledCurve[point] - meanDb;
}

void ApplyBodyDecay(CrashCymbalFitParameters &fit,
                    const CrashMacroValues &values) noexcept {
  const auto points = ReadCurve<BodyDecayCurvePointCount>(
      values, CrashMacro::BodyDecayFrequencyFirst,
      CrashMacro::BodyDecaySecondsFirst);
  for (std::size_t point = 0; point < points.size(); ++point) {
    fit.bodyDecayFrequencyHz[point] = points[point].frequency;
    fit.bodyDecaySeconds[point] = points[point].level;
  }
}

void ApplyTurbulence(CrashCymbalFitParameters &fit,
                     const CrashMacroValues &values) noexcept {
  const auto points = ReadCurve<TurbulenceCurvePointCount>(
      values, CrashMacro::TurbulenceFrequencyFirst,
      CrashMacro::TurbulenceLevelFirst);
  float meanDb = 0.f;
  for (const auto &point : points)
    meanDb += point.level / points.size();
  const float amount = Value(values, CrashMacro::TurbulenceAmount);
  const bool enabled = Value(values, CrashMacro::TurbulenceEnabled) >= .5f;
  for (std::size_t point = 0; point < points.size(); ++point) {
    fit.turbulenceFrequencyHz[point] = points[point].frequency;
    fit.turbulenceGain[point] = enabled
        ? amount * std::pow(10.f, (points[point].level - meanDb) / 20.f)
        : 0.f;
  }
  fit.turbulencePersistence =
      Value(values, CrashMacro::TurbulencePersistence);
}

} // namespace

const CrashMacroDescriptor &CrashMacroDescription(
    const std::size_t index) noexcept {
  return Descriptors[std::min(index, Descriptors.size() - 1)];
}

CrashMacroValues DefaultCrashMacros() noexcept {
  CrashMacroValues result{};
  for (std::size_t index = 0; index < result.size(); ++index)
    result[index] = Descriptors[index].defaultValue;
  return result;
}

CrashCymbalFitParameters ApplyCrashMacros(
    const CrashCymbalFitParameters &base,
    const CrashMacroValues &values) noexcept {
  auto fit = base;
  fit.outputGain = base.outputGain *
      std::pow(10.f, Value(values, CrashMacro::ModelLevelDb) / 20.f);

  const float impact = Value(values, CrashMacro::ImpactToneNoise);
  const float tonal = 1.41421356237f * std::cos(PiOverTwo * impact);
  const float noisy = 1.41421356237f * std::sin(PiOverTwo * impact);
  fit.contactPulseGain = base.contactPulseGain * (1.15f - .3f * impact);
  fit.contactChirpGain = base.contactChirpGain * tonal;
  fit.contactNoiseGain = base.contactNoiseGain * noisy;
  fit.contactMicroGain = base.contactMicroGain * noisy;
  const float width = Value(values, CrashMacro::ImpactWidth);
  fit.contactDurationScale = base.contactDurationScale * width;
  fit.contactNoiseDurationScale = base.contactNoiseDurationScale * width;
  fit.contactMicroDurationScale = base.contactMicroDurationScale * width;

  const float bloom = Value(values, CrashMacro::BloomAmount);
  fit.dispersionDrive = 8.f * bloom;
  fit.dispersionExcursionSamples = 16.f * std::pow(bloom, 1.80708707f);
  const float development = Value(values, CrashMacro::BloomDevelopment);
  fit.dispersionFeedback = 1.f - std::exp(
      std::log(.03f) + development * (std::log(.00025f) - std::log(.03f)));
  const float developmentScale =
      std::exp2(2.f * (development - DefaultBloomDevelopment));
  fit.dispersionLowDecaySeconds =
      base.dispersionLowDecaySeconds * developmentScale;
  fit.dispersionMiddleDecaySeconds =
      base.dispersionMiddleDecaySeconds * developmentScale;
  fit.dispersionHighDecaySeconds =
      base.dispersionHighDecaySeconds * developmentScale;

  const float body = PiOverTwo * Value(values, CrashMacro::BodyToneWash);
  const float bodyGain = std::hypot(base.sparseGain, base.denseGain);
  fit.sparseGain = bodyGain * std::cos(body);
  fit.denseGain = bodyGain * std::sin(body);
  fit.denseTiltDbPerOctave = Value(values, CrashMacro::BodyBrightness);
  fit.directGain = Value(values, CrashMacro::DirectGain);

  fit.denseMinimumFrequencyHz = Value(values, CrashMacro::DenseMinimumFrequency);
  fit.denseMaximumFrequencyHz = Value(values, CrashMacro::DenseMaximumFrequency);
  fit.denseModeDensity = Value(values, CrashMacro::DenseModeDensity);
  fit.denseSpacingJitter = Value(values, CrashMacro::DenseSpacingJitter);
  fit.denseDecaySpreadOctaves = Value(values, CrashMacro::DenseDecaySpread);
  fit.denseGainSpreadDb = Value(values, CrashMacro::DenseGainSpread);
  ApplyResolvedPaint(fit, values);
  ApplyDenseWashPaint(fit, values);
  ApplyBodyDecay(fit, values);
  ApplyTurbulence(fit, values);

  fit.directRadiationEnabled =
      Value(values, CrashMacro::DirectRadiationEnabled) >= .5f;
  fit.directLowCutHz = Value(values, CrashMacro::DirectLowCut);
  fit.directLowCutQ = Value(values, CrashMacro::DirectLowCutQ);
  fit.directColourFrequencyHz = Value(values, CrashMacro::DirectColourFrequency);
  fit.directColourGainDb = Value(values, CrashMacro::DirectColourGain);
  fit.directColourQ = Value(values, CrashMacro::DirectColourQ);
  fit.directHighCutHz = Value(values, CrashMacro::DirectHighCut);
  fit.directHighCutQ = Value(values, CrashMacro::DirectHighCutQ);
  fit.sparseRadiationEnabled =
      Value(values, CrashMacro::SparseRadiationEnabled) >= .5f;
  fit.sparseLowCutHz = Value(values, CrashMacro::SparseLowCut);
  fit.sparseLowCutQ = Value(values, CrashMacro::SparseLowCutQ);
  fit.colourFrequencyHz = Value(values, CrashMacro::SparseColourFrequency);
  fit.colourGainDb = Value(values, CrashMacro::SparseColourGain);
  fit.sparseColourQ = Value(values, CrashMacro::SparseColourQ);
  fit.highCutHz = Value(values, CrashMacro::SparseHighCut);
  fit.sparseHighCutQ = Value(values, CrashMacro::SparseHighCutQ);
  fit.denseRadiationEnabled =
      Value(values, CrashMacro::DenseRadiationEnabled) >= .5f;
  fit.denseLowCutHz = Value(values, CrashMacro::DenseLowCut);
  fit.denseLowCutQ = Value(values, CrashMacro::DenseLowCutQ);
  fit.denseColourFrequencyHz = Value(values, CrashMacro::DenseColourFrequency);
  fit.denseColourGainDb = Value(values, CrashMacro::DenseColourGain);
  fit.denseColourQ = Value(values, CrashMacro::DenseColourQ);
  fit.denseHighCutHz = Value(values, CrashMacro::DenseHighCut);
  fit.denseHighCutQ = Value(values, CrashMacro::DenseHighCutQ);
  return fit;
}

} // namespace tfworkbench
