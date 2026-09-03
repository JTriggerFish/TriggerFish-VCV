#include "crash_cymbal_parameters.hpp"

#include "erb_scale.hpp"

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

class BodyDecayEnvelope {
public:
  BodyDecayEnvelope(const float sampleRate,
                    const CrashCymbalFitParameters &fit) noexcept {
    const float nyquist = .5f * sampleRate;
    for (std::size_t index = 0; index < CrashBodyDecayPointCount; ++index) {
      if (index != 0 && index + 1 != CrashBodyDecayPointCount &&
          !fit.bodyDecayActive[index])
        continue;
      const float frequency = index == 0 ? 0.f :
          index + 1 == CrashBodyDecayPointCount ? nyquist :
          std::clamp(tfdsp::FiniteNormalOrZero(
                         fit.bodyDecayFrequencyHz[index]),
                     1.f, std::max(1.f, nyquist - 1.f));
      const float seconds = std::clamp(
          Positive(fit.bodyDecaySeconds[index], 1.f), .01f, 30.f);
      points_[count_++] = {ErbRate(frequency), std::log(seconds)};
    }
    std::sort(points_.begin(), points_.begin() + count_,
              [](const auto &left, const auto &right) {
                return left.erbRate < right.erbRate;
              });
  }

  float At(const float frequency) const noexcept {
    const float target = ErbRate(std::max(
        0.f, tfdsp::FiniteNormalOrZero(frequency)));
    if (target <= points_[0].erbRate)
      return std::exp(points_[0].logSeconds);
    for (std::size_t right = 1; right < count_; ++right) {
      if (target > points_[right].erbRate)
        continue;
      const auto &left = points_[right - 1];
      const float denominator = points_[right].erbRate - left.erbRate;
      const float amount = denominator > 1.e-6f
          ? std::clamp((target - left.erbRate) / denominator, 0.f, 1.f)
          : 1.f;
      return std::exp(left.logSeconds + amount *
          (points_[right].logSeconds - left.logSeconds));
    }
    return std::exp(points_[count_ - 1].logSeconds);
  }

private:
  struct Point {
    float erbRate{};
    float logSeconds{};
  };
  std::array<Point, CrashBodyDecayPointCount> points_{};
  std::size_t count_{};
};

float ErbBandwidth(const float frequencyHz) noexcept {
  return 24.7f * (1.f + .00437f * frequencyHz);
}

CrashModalField::Parameters ModalField(
    const float sampleRate, const CrashCymbalFitParameters &fit,
    const BodyDecayEnvelope &decay) noexcept {
  struct Anchor {
    float frequencyHz;
    float decayRatio;
    float amplitude;
    float phaseRadians;
    float turbulenceScale;
  };
  std::array<Anchor, CrashSparseModeCount> anchors{};
  for (std::size_t index = 0; index < anchors.size(); ++index) {
    anchors[index] = {
        fit.sparseFrequencyHz[index], fit.sparseDecayRatio[index],
        fit.sparseAmplitude[index], fit.sparsePhaseRadians[index],
        fit.fieldTurbulenceScale[index]};
  }
  std::sort(anchors.begin(), anchors.end(), [](const auto &left,
                                                const auto &right) {
    return left.frequencyHz < right.frequencyHz;
  });
  CrashModalField::Parameters result{};
  DeterministicRandom random;
  random.Seed(fit.denseModeSeed ^ 0x4649454cu);
  const float globalTurbulence = std::clamp(fit.fieldTurbulence, 0.f, 1.f);
  std::array<float, CrashSparseModeCount> anchorGains{};
  float anchorSquaredGain = 0.f;
  for (std::size_t anchor = 0; anchor < anchorGains.size(); ++anchor) {
    const float frequency = Positive(anchors[anchor].frequencyHz, 1000.f);
    const float tilt = std::pow(
        10.f, fit.denseTiltDbPerOctave * std::log2(frequency / 4000.f) / 20.f);
    const float gain = std::clamp(
        Positive(anchors[anchor].amplitude, 0.f) * tilt, 0.f, 8.f);
    anchorGains[anchor] = gain;
    anchorSquaredGain += gain * gain;
  }
  const float anchorNormalization =
      1.f / std::sqrt(std::max(anchorSquaredGain, 1.e-12f));
  constexpr float Pi = 3.14159265358979323846f;
  constexpr std::size_t PairCount = (CrashPacketModeCount - 1) / 2;

  for (std::size_t anchor = 0; anchor < CrashSparseModeCount; ++anchor) {
    const float turbulence = std::clamp(
        globalTurbulence * anchors[anchor].turbulenceScale, 0.f, 1.f);
    const float diffuseEnergy = .9f * turbulence * turbulence;
    const float coreWeight = std::sqrt(1.f - diffuseEnergy);
    const float satelliteWeight = std::sqrt(
        diffuseEnergy / static_cast<float>(CrashPacketModeCount - 1));
    const float spreadErb = turbulence * std::clamp(
        fit.fieldPacketSpreadErb, 0.f, 12.f);
    const float bandwidthErb = turbulence * turbulence * std::clamp(
        fit.fieldPhaseBandwidthErb, 0.f, 4.f);
    const float centre = std::clamp(
        Positive(anchors[anchor].frequencyHz, 1000.f) *
            std::clamp(Positive(fit.sparseTune, 1.f), .5f, 2.f),
        20.f, .48f * sampleRate);
    const float anchorGain = anchorNormalization * anchorGains[anchor];
    const float decayRatio = std::clamp(
        Positive(anchors[anchor].decayRatio, 1.f), .5f, 2.f);
    const auto makeMode = [&](const std::size_t slot, const float frequency,
                              const float weight, const float phase) {
      const std::size_t index = anchor * CrashPacketModeCount + slot;
      const float safeFrequency = std::clamp(frequency, 20.f,
                                              .48f * sampleRate);
      result[index] = {
          safeFrequency,
          std::clamp(decay.At(safeFrequency) * decayRatio, .01f, 30.f),
          ModalDriveGain * anchorGain * weight,
          1.f,
          phase,
          bandwidthErb * ErbBandwidth(safeFrequency) *
              (slot == 0 ? .35f : 1.f),
          static_cast<std::uint16_t>(anchor),
          turbulence * turbulence};
    };

    makeMode(0, centre, coreWeight,
             std::clamp(anchors[anchor].phaseRadians, -Pi, Pi));
    const float centreErb = ErbRate(centre);
    for (std::size_t pair = 0; pair < PairCount; ++pair) {
      const float radius = std::pow(
          (static_cast<float>(pair) + 1.f) / static_cast<float>(PairCount),
          .72f);
      const float jitter = .85f + .3f * random.Uniform();
      const float offset = spreadErb * radius * jitter;
      const float low = InverseErbRate(std::max(ErbRate(20.f),
                                                centreErb - offset));
      const float high = InverseErbRate(std::min(ErbRate(.48f * sampleRate),
                                                 centreErb + offset));
      makeMode(1 + 2 * pair, low, satelliteWeight, Pi * random.Bipolar());
      makeMode(2 + 2 * pair, high, satelliteWeight, Pi * random.Bipolar());
    }
    auto first = result.begin() + anchor * CrashPacketModeCount;
    std::sort(first, first + CrashPacketModeCount,
              [](const auto &left, const auto &right) {
                return left.frequencyHz < right.frequencyHz;
              });
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
  const float diffusion = std::clamp(fit.dispersionDiffusion, 0.f, 1.f);
  result.firstAllpassGain = diffusion * .61f;
  result.secondAllpassDelaySamples = 13.f * scale;
  result.secondAllpassGain = diffusion * -.53f;
  result.thirdAllpassDelaySamples = 5.f * scale;
  result.thirdAllpassGain = diffusion * .3f;
  result.fourthAllpassDelaySamples = 17.f * scale;
  result.fourthAllpassGain = diffusion * -.25f;
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

ObservationModel<2>::Parameters Observation(
    const CrashCymbalFitParameters &fit) noexcept {
  ObservationModel<2>::Parameters result{};
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
  result[1].gain = std::clamp(fit.fieldGain, 0.f, 4.f);
  result[1].radiationEnabled = fit.denseRadiationEnabled;
  result[1].radiation.lowCutHz =
      std::clamp(fit.denseLowCutHz, 10.f, 1000.f);
  result[1].radiation.lowCutQ = std::clamp(fit.denseLowCutQ, .25f, 4.f);
  result[1].radiation.colourFrequencyHz =
      std::clamp(fit.denseColourFrequencyHz, 100.f, 18000.f);
  result[1].radiation.colourGainDb =
      std::clamp(fit.denseColourGainDb, -18.f, 18.f);
  result[1].radiation.colourQ = std::clamp(fit.denseColourQ, .25f, 12.f);
  result[1].radiation.highCutHz =
      std::clamp(fit.denseHighCutHz, 1000.f, 22000.f);
  result[1].radiation.highCutQ = std::clamp(fit.denseHighCutQ, .25f, 4.f);
  return result;
}

} // namespace

CrashCymbalParameters DefaultCrashCymbalParameters(
    const float sampleRate, const CrashCymbalFitParameters &fit) {
  CrashCymbalParameters result;
  result.fit = fit;
  const BodyDecayEnvelope decay(sampleRate, fit);
  result.modalField = ModalField(sampleRate, fit, decay);
  SetLocationProjections<CrashModalFieldModeCount>(
      result.modalField, result.fieldBellProjection,
      result.fieldBowProjection, result.fieldEdgeProjection);
  result.modalFieldControls.exchangeAngleRadians =
      .012f * std::clamp(fit.fieldExchange, 0.f, 1.f);
  result.modalFieldControls.seed = fit.denseModeSeed ^ 0x4649454cu;
  result.dispersion = Dispersion(sampleRate, fit);
  result.observation = Observation(fit);
  return result;
}

CrashCymbalPreparedParameters PrepareCrashCymbalParameters(
    const float sampleRate, const CrashCymbalParameters &parameters) {
  CrashCymbalPreparedParameters result;
  result.parameters = parameters;
  result.modalField = PrepareStochasticModalField(
      sampleRate, parameters.modalField, parameters.modalFieldControls,
      700.f, 6500.f);
  result.sampleRate = sampleRate;
  return result;
}

} // namespace tfdsp::percussion
