#include "crash_macros.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace tfworkbench {
namespace {

constexpr float PiOverTwo = 1.57079632679489661923f;
constexpr float DefaultBloomDevelopment = .44875992f;

constexpr std::array<CrashMacroDescriptor, CrashMacroCount> Descriptors{{
    {"Model level", "dB", -60.f, 12.f, -36.f},
    {"Impact tone / noise", "", 0.f, 1.f, .5f},
    {"Impact width", "", 0.f, 1.f, .5f},
    {"Bloom amount", "", 0.f, 1.f, .35f},
    {"Bloom development", "", 0.f, 1.f, DefaultBloomDevelopment},
    {"Body tone / wash", "", 0.f, 1.f, .68554716f},
    {"Body brightness", "", 0.f, 1.f, 7.f / 13.f},
    {"Low decay", "", 0.f, 1.f, .5f},
    {"Middle decay", "", 0.f, 1.f, .5f},
    {"High decay", "", 0.f, 1.f, .5f},
}};

float Value(const CrashMacroValues &values, const CrashMacro macro) noexcept {
  const auto index = static_cast<std::size_t>(macro);
  const auto &descriptor = Descriptors[index];
  const float value = std::isfinite(values[index])
      ? values[index] : descriptor.defaultValue;
  return std::clamp(value, descriptor.minimum, descriptor.maximum);
}

float DecayScale(const float value) noexcept {
  return std::exp2(4.f * (value - .5f));
}

float BandScale(const float position, const float low, const float middle,
                const float high) noexcept {
  return position < .5f
      ? std::exp(std::log(low) + 2.f * position * (std::log(middle) - std::log(low)))
      : std::exp(std::log(middle) + (2.f * position - 1.f) *
          (std::log(high) - std::log(middle)));
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

tfdsp::percussion::CrashCymbalFitParameters ApplyCrashMacros(
    const tfdsp::percussion::CrashCymbalFitParameters &base,
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

  const float width =
      std::exp2(4.f * (Value(values, CrashMacro::ImpactWidth) - .5f));
  fit.contactDurationScale = base.contactDurationScale * width;
  fit.contactNoiseDurationScale = base.contactNoiseDurationScale * width;
  fit.contactMicroDurationScale = base.contactMicroDurationScale * width;

  const float bloom = Value(values, CrashMacro::BloomAmount);
  fit.dispersionDrive = 8.f * bloom;
  fit.dispersionExcursionSamples = 16.f * std::pow(bloom, 1.80708707f);
  const float development = Value(values, CrashMacro::BloomDevelopment);
  fit.dispersionFeedback = 1.f -
      std::exp(std::log(.03f) + development * (std::log(.00025f) - std::log(.03f)));
  const float developmentScale =
      std::exp2(2.f * (development - DefaultBloomDevelopment));
  fit.dispersionLowDecaySeconds = base.dispersionLowDecaySeconds * developmentScale;
  fit.dispersionMiddleDecaySeconds =
      base.dispersionMiddleDecaySeconds * developmentScale;
  fit.dispersionHighDecaySeconds =
      base.dispersionHighDecaySeconds * developmentScale;

  const float body = PiOverTwo * Value(values, CrashMacro::BodyToneWash);
  const float bodyGain = std::hypot(base.sparseGain, base.denseGain);
  fit.sparseGain = bodyGain * std::cos(body);
  fit.denseGain = bodyGain * std::sin(body);
  fit.denseTiltDbPerOctave =
      -8.f + 13.f * Value(values, CrashMacro::BodyBrightness);

  const float low = DecayScale(Value(values, CrashMacro::LowDecay));
  const float middle = DecayScale(Value(values, CrashMacro::MiddleDecay));
  const float high = DecayScale(Value(values, CrashMacro::HighDecay));
  const float minimumLog = std::log2(500.f);
  const float maximumLog = std::log2(18000.f);
  for (std::size_t mode = 0; mode < fit.sparseDecaySeconds.size(); ++mode) {
    const float position = std::clamp(
        (std::log2(std::max(base.sparseFrequencyHz[mode], 500.f)) - minimumLog) /
            (maximumLog - minimumLog),
        0.f, 1.f);
    fit.sparseDecaySeconds[mode] = base.sparseDecaySeconds[mode] *
        BandScale(position, low, middle, high);
  }
  for (std::size_t point = 0; point < fit.denseDecayEnvelopeOctaves.size(); ++point) {
    const float position = static_cast<float>(point) /
        static_cast<float>(fit.denseDecayEnvelopeOctaves.size() - 1);
    fit.denseDecayEnvelopeOctaves[point] =
        base.denseDecayEnvelopeOctaves[point] +
        std::log2(BandScale(position, low, middle, high));
  }
  fit.turbulenceLowDecaySeconds = base.turbulenceLowDecaySeconds * low;
  fit.turbulenceMiddleDecaySeconds = base.turbulenceMiddleDecaySeconds * middle;
  fit.turbulenceHighDecaySeconds = base.turbulenceHighDecaySeconds * high;
  return fit;
}

} // namespace tfworkbench
