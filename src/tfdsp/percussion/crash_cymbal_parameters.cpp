#include "crash_cymbal_parameters.hpp"

#include "statistical_modal_cloud.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace tfdsp::percussion {
namespace {

// Contact body drive and generalized modal force use different arbitrary units.
// ModalBank itself applies no hidden decay-dependent gain normalization.
constexpr float ModalDriveGain = .05f;

float Positive(const float value, const float fallback) noexcept {
  return std::isfinite(value) && value > 0.f ? value : fallback;
}

template <std::size_t Count>
float LogEnvelope(const std::array<float, Count> &frequencies,
                  const std::array<float, Count> &values,
                  const float frequency) noexcept {
  std::array<std::pair<float, float>, Count> points{};
  for (std::size_t point = 0; point < Count; ++point) {
    points[point] = {Positive(frequencies[point], 1.f),
                     Positive(values[point], 1.f)};
  }
  std::sort(points.begin(), points.end(), [](const auto &left,
                                              const auto &right) {
    return left.first < right.first;
  });
  const float target = Positive(frequency, points.front().first);
  if (target <= points.front().first)
    return points.front().second;
  for (std::size_t right = 1; right < Count; ++right) {
    if (target > points[right].first)
      continue;
    const float leftFrequency = points[right - 1].first;
    const float rightFrequency = std::max(
        points[right].first, leftFrequency + 1.f);
    const float amount = std::clamp(
        std::log(target / leftFrequency) /
            std::log(rightFrequency / leftFrequency),
        0.f, 1.f);
    const float left = points[right - 1].second;
    const float rightValue = points[right].second;
    return std::exp(std::log(left) + amount * (std::log(rightValue) -
                                                std::log(left)));
  }
  return points.back().second;
}

float BodyDecay(const CrashCymbalFitParameters &fit,
                const float frequency) noexcept {
  return std::clamp(LogEnvelope(fit.bodyDecayFrequencyHz,
                                fit.bodyDecaySeconds, frequency),
                    .01f, 30.f);
}

CrashSparseModes::Parameters SparseModes(
    const float sampleRate, const CrashCymbalFitParameters &fit) noexcept {
  CrashSparseModes::Parameters result{};
  const float tune = std::clamp(Positive(fit.sparseTune, 1.f), .5f, 2.f);
  float squaredGain = 0.f;
  for (const float amplitude : fit.sparseAmplitude) {
    const float gain = std::clamp(Positive(amplitude, 0.f), 0.f, 8.f);
    squaredGain += gain * gain;
  }
  const float normalization = 1.f / std::sqrt(std::max(squaredGain, 1.e-12f));
  for (std::size_t mode = 0; mode < result.size(); ++mode) {
    const float frequency = std::clamp(
        Positive(fit.sparseFrequencyHz[mode], 1000.f) * tune,
        20.f, .48f * sampleRate);
    result[mode] = {
        frequency,
        std::clamp(BodyDecay(fit, frequency) *
                       std::clamp(Positive(fit.sparseDecayRatio[mode], 1.f),
                                  .5f, 2.f),
                   .01f, 30.f),
        ModalDriveGain,
        normalization * std::clamp(Positive(fit.sparseAmplitude[mode], 0.f),
                                   0.f, 8.f),
        std::clamp(fit.sparsePhaseRadians[mode],
                   -3.14159265358979323846f, 3.14159265358979323846f)};
  }
  return result;
}

CrashDenseModes::Parameters DenseModes(
    const float sampleRate, const CrashCymbalFitParameters &fit,
    const bool extension) {
  StatisticalModalCloudParameters cloud;
  cloud.minimumFrequencyHz = fit.denseMinimumFrequencyHz;
  cloud.maximumFrequencyHz = fit.denseMaximumFrequencyHz;
  cloud.frequencyWarp = fit.denseFrequencyWarp;
  cloud.spacingJitter = fit.denseSpacingJitter;
  const float requestedBanks = std::clamp(fit.denseModeDensity, 0.f, 2.f);
  const float bankDensity = extension
      ? std::max(requestedBanks - 1.f, 0.f)
      : std::min(requestedBanks, 1.f);
  const float totalModes = requestedBanks * CrashDenseModeCount;
  const float bankModes = bankDensity * CrashDenseModeCount;
  cloud.modeDensity = bankDensity;
  cloud.lowDecaySeconds = 1.f;
  cloud.highDecaySeconds = 1.f;
  cloud.decayCurve = 1.f;
  const float minimum = std::clamp(fit.denseMinimumFrequencyHz, 20.f,
                                   .45f * sampleRate);
  const float maximum = std::clamp(fit.denseMaximumFrequencyHz,
                                   minimum + 1.f, .48f * sampleRate);
  const float minimumErb = ErbRate(minimum);
  const float maximumErb = ErbRate(maximum);
  for (std::size_t point = 0; point < cloud.decayEnvelopeOctaves.size(); ++point) {
    const float position = static_cast<float>(point) /
        static_cast<float>(cloud.decayEnvelopeOctaves.size() - 1);
    const float frequency = InverseErbRate(
        minimumErb + position * (maximumErb - minimumErb));
    cloud.decayEnvelopeOctaves[point] = std::log2(BodyDecay(fit, frequency));
  }
  cloud.decaySpreadOctaves = fit.denseDecaySpreadOctaves;
  cloud.tiltDbPerOctave = fit.denseTiltDbPerOctave;
  cloud.gainEnvelopeDb = fit.denseGainEnvelopeDb;
  cloud.gainSpreadDb = fit.denseGainSpreadDb;
  cloud.outputGain = totalModes > 0.f
      ? std::sqrt(bankModes / totalModes) : 0.f;
  cloud.seed = extension
      ? fit.denseModeSeed ^ 0x4558544eu : fit.denseModeSeed;
  auto result = MakeStatisticalModalCloud<CrashDenseModeCount>(sampleRate, cloud);
  // The painted cloud envelope is a normalized excitation-energy
  // distribution. Radiation and branch level remain observation concerns.
  // For independent linear modes this factorization preserves the transfer
  // product while giving the stored state the intended physical semantics.
  for (auto &mode : result) {
    mode.inputGain = ModalDriveGain * mode.outputGain;
    mode.outputGain = mode.inputGain != 0.f ? 1.f : 0.f;
  }
  return result;
}

template <std::size_t Count, typename Parameters>
void SetLocationProjections(const Parameters &modes,
                            std::array<float, Count> &bell,
                            std::array<float, Count> &bow,
                            std::array<float, Count> &edge) noexcept {
  for (std::size_t mode = 0; mode < Count; ++mode) {
    const float frequency = modes[mode].frequencyHz;
    const float position = static_cast<float>(mode) /
        static_cast<float>(Count - 1);
    const float logDistance = std::log2(std::max(frequency, 1.f) / 4200.f);
    const float bellFocus = std::exp(-.5f * logDistance * logDistance);
    bell[mode] = .25f + 1.15f * bellFocus;
    bow[mode] = .72f + .22f *
        std::abs(std::sin(3.14159265358979323846f * (2.7f * position + .13f)));
    edge[mode] = .65f + .55f * std::sqrt(position);
  }
}

DispersionLoopParameters Dispersion(
    const float sampleRate, const CrashCymbalFitParameters &fit) noexcept {
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
  result.thirdAllpassDelaySamples = 5.f * scale;
  result.thirdAllpassGain = .3f;
  result.fourthAllpassDelaySamples = 17.f * scale;
  result.fourthAllpassGain = -.25f;
  result.selfPhase.centreDelaySamples = 9.f * scale;
  result.selfPhase.maximumExcursionSamples =
      std::clamp(fit.dispersionExcursionSamples, 0.f, 16.f) * scale;
  result.selfPhase.drive = std::clamp(fit.dispersionDrive, 0.f, 8.f);
  result.selfPhase.toneHz = 7600.f;
  result.selfPhase.envelopeReleaseSeconds = .018f;
  result.selfPhase.normalization = .55f;
  result.decay = {
      std::clamp(fit.dispersionLowDecaySeconds, .05f, 5.f),
      std::clamp(fit.dispersionMiddleDecaySeconds, .05f, 5.f),
      std::clamp(fit.dispersionHighDecaySeconds, .05f, 5.f)};
  result.feedbackGain = std::clamp(fit.dispersionFeedback, 0.f, .9999f);
  result.lowCrossoverHz = 700.f;
  result.highCrossoverHz = 6500.f;
  result.modulationSeed = 0x43524153u;
  return result;
}

TurbulentResidualParameters Turbulence(
    const CrashCymbalFitParameters &fit) noexcept {
  TurbulentResidualParameters result;
  for (std::size_t band = 0; band < result.gain.size(); ++band)
    result.gain[band] = std::clamp(fit.turbulenceGain[band], 0.f, 4.f);
  std::array<float, 3> frequencies{};
  frequencies[0] = std::clamp(fit.turbulenceFrequencyHz[0], 40.f, 8000.f);
  frequencies[1] = std::clamp(fit.turbulenceFrequencyHz[1],
                              frequencies[0] * 1.05f, 14000.f);
  frequencies[2] = std::clamp(fit.turbulenceFrequencyHz[2],
                              frequencies[1] * 1.05f, 20000.f);
  const float persistence = std::clamp(
      Positive(fit.turbulencePersistence, 1.f), .25f, 4.f);
  result.decay = {
      persistence * BodyDecay(fit, frequencies[0]),
      persistence * BodyDecay(fit, frequencies[1]),
      persistence * BodyDecay(fit, frequencies[2])};
  result.lowCrossoverHz = std::sqrt(frequencies[0] * frequencies[1]);
  result.highCrossoverHz = std::sqrt(frequencies[1] * frequencies[2]);
  result.seed = 0x43524153u;
  return result;
}

ObservationModel<4>::Parameters Observation(
    const CrashCymbalFitParameters &fit) noexcept {
  ObservationModel<4>::Parameters result{};
  result[0].gain = std::clamp(fit.directGain, 0.f, 4.f);
  result[0].radiationEnabled = fit.directRadiationEnabled;
  result[0].radiation.lowCutHz =
      std::clamp(fit.directLowCutHz, 10.f, 1000.f);
  result[0].radiation.lowCutQ = std::clamp(fit.directLowCutQ, .25f, 4.f);
  result[0].radiation.colourFrequencyHz =
      std::clamp(fit.directColourFrequencyHz, 100.f, 18000.f);
  result[0].radiation.colourGainDb =
      std::clamp(fit.directColourGainDb, -18.f, 18.f);
  result[0].radiation.colourQ = std::clamp(fit.directColourQ, .25f, 12.f);
  result[0].radiation.highCutHz =
      std::clamp(fit.directHighCutHz, 1000.f, 22000.f);
  result[0].radiation.highCutQ = std::clamp(fit.directHighCutQ, .25f, 4.f);
  result[1].gain = std::clamp(fit.sparseGain, 0.f, 4.f);
  result[1].radiationEnabled = fit.sparseRadiationEnabled;
  result[1].radiation.lowCutHz =
      std::clamp(fit.sparseLowCutHz, 10.f, 1000.f);
  result[1].radiation.lowCutQ = std::clamp(fit.sparseLowCutQ, .25f, 4.f);
  result[1].radiation.colourFrequencyHz =
      std::clamp(fit.colourFrequencyHz, 100.f, 18000.f);
  result[1].radiation.colourGainDb =
      std::clamp(fit.colourGainDb, -18.f, 18.f);
  result[1].radiation.colourQ = std::clamp(fit.sparseColourQ, .25f, 12.f);
  result[1].radiation.highCutHz =
      std::clamp(fit.highCutHz, 1000.f, 22000.f);
  result[1].radiation.highCutQ = std::clamp(fit.sparseHighCutQ, .25f, 4.f);
  result[2].gain = std::clamp(fit.denseGain, 0.f, 4.f);
  result[2].radiationEnabled = fit.denseRadiationEnabled;
  result[2].radiation.lowCutHz =
      std::clamp(fit.denseLowCutHz, 10.f, 1000.f);
  result[2].radiation.lowCutQ = std::clamp(fit.denseLowCutQ, .25f, 4.f);
  result[2].radiation.colourFrequencyHz =
      std::clamp(fit.denseColourFrequencyHz, 100.f, 18000.f);
  result[2].radiation.colourGainDb =
      std::clamp(fit.denseColourGainDb, -18.f, 18.f);
  result[2].radiation.colourQ = std::clamp(fit.denseColourQ, .25f, 12.f);
  result[2].radiation.highCutHz =
      std::clamp(fit.denseHighCutHz, 1000.f, 22000.f);
  result[2].radiation.highCutQ = std::clamp(fit.denseHighCutQ, .25f, 4.f);
  // Turbulence is generated and observed independently from the dense modal
  // cloud. It starts with the same radiation curve, but owns separate filter
  // state and is never scaled by the resolved-to-wash balance.
  result[3] = result[2];
  result[3].gain = 1.f;
  return result;
}

} // namespace

CrashCymbalParameters DefaultCrashCymbalParameters(
    const float sampleRate, const CrashCymbalFitParameters &fit) {
  CrashCymbalParameters result;
  result.fit = fit;
  result.sparseModes = SparseModes(sampleRate, fit);
  result.denseModes = DenseModes(sampleRate, fit, false);
  result.denseExtensionModes = DenseModes(sampleRate, fit, true);
  SetLocationProjections<CrashSparseModeCount>(
      result.sparseModes, result.sparseBellProjection,
      result.sparseBowProjection, result.sparseEdgeProjection);
  SetLocationProjections<CrashDenseModeCount>(
      result.denseModes, result.denseBellProjection,
      result.denseBowProjection, result.denseEdgeProjection);
  SetLocationProjections<CrashDenseModeCount>(
      result.denseExtensionModes, result.denseExtensionBellProjection,
      result.denseExtensionBowProjection,
      result.denseExtensionEdgeProjection);
  result.dispersion = Dispersion(sampleRate, fit);
  result.turbulence = Turbulence(fit);
  result.observation = Observation(fit);
  return result;
}

} // namespace tfdsp::percussion
