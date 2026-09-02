#include "crash_cymbal_parameters.hpp"

#include "statistical_modal_cloud.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion {
namespace {

// Contact body drive and generalized modal force use different arbitrary units.
// ModalBank itself applies no hidden decay-dependent gain normalization.
constexpr float ModalDriveGain = .05f;

float Positive(const float value, const float fallback) noexcept {
  return std::isfinite(value) && value > 0.f ? value : fallback;
}

CrashSparseModes::Parameters SparseModes(
    const float sampleRate, const CrashCymbalFitParameters &fit) noexcept {
  CrashSparseModes::Parameters result{};
  const float tune = std::clamp(Positive(fit.sparseTune, 1.f), .5f, 2.f);
  const float decayScale =
      std::clamp(Positive(fit.sparseDecayScale, 1.f), .1f, 4.f);
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
        std::clamp(Positive(fit.sparseDecaySeconds[mode], 1.f) * decayScale,
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
    const float sampleRate, const CrashCymbalFitParameters &fit) {
  StatisticalModalCloudParameters cloud;
  cloud.minimumFrequencyHz = fit.denseMinimumFrequencyHz;
  cloud.maximumFrequencyHz = fit.denseMaximumFrequencyHz;
  cloud.frequencyWarp = fit.denseFrequencyWarp;
  cloud.spacingJitter = fit.denseSpacingJitter;
  cloud.lowDecaySeconds = fit.denseLowDecaySeconds;
  cloud.highDecaySeconds = fit.denseHighDecaySeconds;
  cloud.decayCurve = fit.denseDecayCurve;
  cloud.decayEnvelopeOctaves = fit.denseDecayEnvelopeOctaves;
  cloud.decaySpreadOctaves = fit.denseDecaySpreadOctaves;
  cloud.tiltDbPerOctave = fit.denseTiltDbPerOctave;
  cloud.gainEnvelopeDb = fit.denseGainEnvelopeDb;
  cloud.gainSpreadDb = fit.denseGainSpreadDb;
  cloud.outputGain = 1.f;
  cloud.seed = fit.denseModeSeed;
  auto result = MakeStatisticalModalCloud<CrashDenseModeCount>(sampleRate, cloud);
  for (auto &mode : result)
    mode.inputGain = ModalDriveGain;
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
  result.gain = {
      std::clamp(fit.turbulenceLowGain, 0.f, 4.f),
      std::clamp(fit.turbulenceMiddleGain, 0.f, 4.f),
      std::clamp(fit.turbulenceHighGain, 0.f, 4.f)};
  result.decay = {
      std::clamp(fit.turbulenceLowDecaySeconds, .005f, 30.f),
      std::clamp(fit.turbulenceMiddleDecaySeconds, .005f, 30.f),
      std::clamp(fit.turbulenceHighDecaySeconds, .005f, 30.f)};
  result.seed = 0x43524153u;
  return result;
}

ObservationModel<3>::Parameters Observation(
    const CrashCymbalFitParameters &fit) noexcept {
  ObservationModel<3>::Parameters result{};
  result[0].gain = std::clamp(fit.directGain, 0.f, 4.f);
  result[0].radiation.lowCutHz = 180.f;
  result[0].radiation.colourFrequencyHz = 7200.f;
  result[0].radiation.colourGainDb = 1.f;
  result[0].radiation.highCutHz = 20000.f;
  result[1].gain = std::clamp(fit.sparseGain, 0.f, 4.f);
  result[1].radiation.lowCutHz = 90.f;
  result[1].radiation.colourFrequencyHz =
      std::clamp(fit.colourFrequencyHz, 100.f, 18000.f);
  result[1].radiation.colourGainDb =
      std::clamp(fit.colourGainDb, -18.f, 18.f);
  result[1].radiation.highCutHz =
      std::clamp(fit.highCutHz, 1000.f, 22000.f);
  result[2].gain = std::clamp(fit.denseGain, 0.f, 4.f);
  result[2].radiation.lowCutHz = 140.f;
  result[2].radiation.colourFrequencyHz = 7200.f;
  result[2].radiation.colourGainDb = .5f;
  result[2].radiation.highCutHz =
      std::clamp(fit.highCutHz, 1000.f, 22000.f);
  return result;
}

} // namespace

CrashCymbalParameters DefaultCrashCymbalParameters(
    const float sampleRate, const CrashCymbalFitParameters &fit) {
  CrashCymbalParameters result;
  result.fit = fit;
  result.sparseModes = SparseModes(sampleRate, fit);
  result.denseModes = DenseModes(sampleRate, fit);
  SetLocationProjections<CrashSparseModeCount>(
      result.sparseModes, result.sparseBellProjection,
      result.sparseBowProjection, result.sparseEdgeProjection);
  SetLocationProjections<CrashDenseModeCount>(
      result.denseModes, result.denseBellProjection,
      result.denseBowProjection, result.denseEdgeProjection);
  result.dispersion = Dispersion(sampleRate, fit);
  result.turbulence = Turbulence(fit);
  result.observation = Observation(fit);
  return result;
}

} // namespace tfdsp::percussion
