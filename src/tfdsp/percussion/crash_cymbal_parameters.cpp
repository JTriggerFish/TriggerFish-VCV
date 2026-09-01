#include "crash_cymbal_parameters.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace tfdsp::percussion {
namespace {

constexpr std::array<float, CrashResonatorCount> FundamentalsHz{
    183.f, 197.f, 211.f, 229.f, 247.f, 269.f, 283.f, 307.f,
    331.f, 359.f, 379.f, 409.f, 431.f, 463.f, 499.f, 541.f,
    577.f, 619.f, 659.f, 709.f, 761.f, 821.f, 887.f, 953.f,
    1031.f, 1117.f, 1213.f, 1321.f, 1451.f, 1583.f, 1741.f, 1931.f};

float Positive(const float value, const float fallback) noexcept {
  return std::isfinite(value) && value > 0.f ? value : fallback;
}

CrashResonators::Parameters Resonators(
    const float sampleRate, const CrashCymbalFitParameters &fit) {
  CrashResonators::Parameters result{};
  const float tune = std::clamp(Positive(fit.resonanceTune, 1.f), .5f, 2.f);
  for (std::size_t line = 0; line < result.size(); ++line) {
    const float frequency = FundamentalsHz[line];
    const float frequencyRatio = frequency / FundamentalsHz.front();
    result[line].delaySamples = sampleRate / (frequency * tune);
    result[line].decay = {
        13.f * std::pow(frequencyRatio, -.55f) *
            Positive(fit.lowDecayScale, 1.f),
        8.f * std::pow(frequencyRatio, -.55f) *
            Positive(fit.middleDecayScale, 1.f),
        3.2f * std::pow(frequencyRatio, -.6f) *
            Positive(fit.highDecayScale, 1.f)};
    result[line].inputGain = 1.f;
    result[line].outputGain = 1.f / static_cast<float>(CrashResonatorCount);
  }
  return result;
}

void SetProjections(CrashCymbalParameters &parameters) noexcept {
  constexpr float Pi = 3.14159265358979323846f;
  for (std::size_t line = 0; line < CrashResonatorCount; ++line) {
    const float position = static_cast<float>(line) /
        static_cast<float>(CrashResonatorCount - 1);
    const float firstSign = line % 2 == 0 ? 1.f : -1.f;
    const float secondSign = (line * 5 + 1) % 7 < 3 ? -1.f : 1.f;
    parameters.bellProjection[line] = firstSign * (.12f + .88f * position);
    parameters.bowProjection[line] = secondSign *
        (.55f + .42f * std::abs(std::sin(Pi * (2.7f * position + .13f))));
    parameters.edgeProjection[line] = firstSign * secondSign *
        (.68f + .3f * std::abs(std::sin(Pi * (4.1f * position + .31f))));
  }
}

DispersionLoopParameters Dispersion(
    const float sampleRate, const CrashCymbalFitParameters &fit) {
  const float scale = sampleRate / 48000.f;
  DispersionLoopParameters result;
  result.baseDelaySamples = 11.f * scale;
  result.slowDelaySamples = 8.f * scale;
  result.slowDepthSamples = 1.3f * scale;
  result.slowRateHz = .17f;
  result.firstAllpassDelaySamples = 7.f * scale;
  result.firstAllpassGain = .61f;
  result.secondAllpassDelaySamples = 13.f * scale;
  result.secondAllpassGain = -.53f;
  result.selfPhase.centreDelaySamples = 9.f * scale;
  result.selfPhase.maximumExcursionSamples =
      std::clamp(fit.dispersionExcursionSamples, 0.f, 6.f) * scale;
  result.selfPhase.drive = std::clamp(fit.dispersionDrive, 0.f, 8.f);
  result.selfPhase.toneHz = 7600.f;
  result.selfPhase.envelopeReleaseSeconds = .018f;
  result.selfPhase.normalization = .55f;
  result.decay = {
      std::clamp(fit.dispersionLowDecaySeconds, .05f, 5.f),
      std::clamp(fit.dispersionMiddleDecaySeconds, .05f, 5.f),
      std::clamp(fit.dispersionHighDecaySeconds, .05f, 5.f)};
  result.feedbackGain = std::clamp(fit.dispersionFeedback, 0.f, .995f);
  result.lowCrossoverHz = 700.f;
  result.highCrossoverHz = 6500.f;
  result.modulationSeed = 0x43524153u;
  return result;
}

ObservationModel<2>::Parameters Observation(
    const CrashCymbalFitParameters &fit) {
  ObservationModel<2>::Parameters result{};
  result[0].gain = std::clamp(fit.directGain, 0.f, 4.f);
  result[0].radiation.lowCutHz = 180.f;
  result[0].radiation.colourFrequencyHz = 7200.f;
  result[0].radiation.colourGainDb = 1.f;
  result[0].radiation.highCutHz = 20000.f;
  result[1].gain = std::clamp(fit.bodyGain, 0.f, 4.f);
  result[1].radiation.lowCutHz = 90.f;
  result[1].radiation.colourFrequencyHz =
      std::clamp(fit.colourFrequencyHz, 100.f, 18000.f);
  result[1].radiation.colourGainDb =
      std::clamp(fit.colourGainDb, -18.f, 18.f);
  result[1].radiation.highCutHz =
      std::clamp(fit.highCutHz, 1000.f, 22000.f);
  return result;
}

} // namespace

CrashCymbalParameters DefaultCrashCymbalParameters(
    const float sampleRate, const CrashCymbalFitParameters &fit) {
  CrashCymbalParameters result;
  result.fit = fit;
  result.resonators = Resonators(sampleRate, fit);
  SetProjections(result);
  constexpr std::array<float, CrashResonatorBusCount> BaseShiftsHz{
      0.f, 17.f, -23.f, 31.f};
  const float shiftScale = std::clamp(fit.resonatorShiftScale, 0.f, 4.f);
  for (std::size_t bus = 0; bus < CrashResonatorBusCount; ++bus)
    result.resonatorShiftHz[bus] = shiftScale * BaseShiftsHz[bus];
  result.dispersion = Dispersion(sampleRate, fit);
  result.observation = Observation(fit);
  return result;
}

} // namespace tfdsp::percussion
