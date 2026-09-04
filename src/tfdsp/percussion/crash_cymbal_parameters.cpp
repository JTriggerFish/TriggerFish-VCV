#include "crash_cymbal_parameters.hpp"

#include "erb_scale.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace tfdsp::percussion {
namespace {

constexpr std::uint32_t ModalFieldSeed = 0x43524153u;

float Positive(const float value, const float fallback) noexcept {
  return std::isfinite(value) && value > 0.f ? value : fallback;
}

class BodyDecayEnvelope {
public:
  BodyDecayEnvelope(const float sampleRate,
                    const CrashCymbalFitParameters &fit) noexcept {
    const float nyquist = .5f * sampleRate;
    points_[count_++] = {ErbRate(0.f), std::log(std::clamp(
        Positive(fit.bodyDecaySeconds.front(), 1.f), .01f, 30.f))};
    for (std::size_t index = 0; index < CrashBodyDecayInteriorPointCount;
         ++index) {
      if (!fit.bodyDecayActive[index]) continue;
      const float frequency = std::clamp(tfdsp::FiniteNormalOrZero(
          fit.bodyDecayFrequencyHz[index]), 1.f, std::max(1.f, nyquist - 1.f));
      const float seconds = std::clamp(
          Positive(fit.bodyDecaySeconds[index + 1], 1.f), .01f, 30.f);
      points_[count_++] = {ErbRate(frequency), std::log(seconds)};
    }
    points_[count_++] = {ErbRate(nyquist), std::log(std::clamp(
        Positive(fit.bodyDecaySeconds.back(), 1.f), .01f, 30.f))};
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
    float amplitude;
    float turbulenceScale;
  };
  std::array<Anchor, CrashSparseModeCount> anchors{};
  for (std::size_t index = 0; index < anchors.size(); ++index) {
    anchors[index] = {
        fit.sparseFrequencyHz[index], fit.sparseAmplitude[index],
        fit.fieldTurbulenceScale[index]};
  }
  std::sort(anchors.begin(), anchors.end(), [](const auto &left,
                                                const auto &right) {
    return left.frequencyHz < right.frequencyHz;
  });
  CrashModalField::Parameters result{};
  DeterministicRandom random;
  random.Seed(ModalFieldSeed ^ 0x4649454cu);
  const float globalTurbulence = std::clamp(fit.fieldTurbulence, 0.f, 1.f);
  const float turbulenceSlope = std::clamp(
      fit.fieldTurbulenceSlopePerOctave, -1.f, 1.f);
  const float turbulenceCentre = std::clamp(
      Positive(fit.fieldTurbulenceCentreHz, 4000.f), 20.f, .48f * sampleRate);
  const float bodyTiltCentre = std::clamp(
      Positive(fit.bodyTiltCentreHz, 4000.f), 20.f, .48f * sampleRate);
  std::array<float, CrashSparseModeCount> anchorGains{};
  std::array<float, CrashSparseModeCount> anchorOutputGains{};
  float anchorSquaredGain = 0.f;
  for (std::size_t anchor = 0; anchor < anchorGains.size(); ++anchor) {
    const float frequency = Positive(anchors[anchor].frequencyHz, 1000.f);
    const float tilt = std::pow(
        10.f, std::clamp(fit.bodyTiltDbPerOctave, -24.f, 24.f) *
            std::log2(frequency / bodyTiltCentre) / 20.f);
    const float level = std::clamp(
        Positive(anchors[anchor].amplitude, 0.f), 0.f, 8.f);
    const float observationGain = std::sqrt(level);
    const float gain = observationGain * tilt;
    anchorGains[anchor] = gain;
    anchorOutputGains[anchor] = observationGain;
    anchorSquaredGain += gain * gain;
  }
  const float anchorNormalization =
      1.f / std::sqrt(std::max(anchorSquaredGain, 1.e-12f));
  constexpr float Pi = 3.14159265358979323846f;
  constexpr std::size_t PairCount = (CrashPacketModeCount - 1) / 2;

  for (std::size_t anchor = 0; anchor < CrashSparseModeCount; ++anchor) {
    const float spectralTurbulence = std::clamp(
        globalTurbulence + turbulenceSlope * std::log2(
            std::max(anchors[anchor].frequencyHz, 20.f) / turbulenceCentre),
        0.f, 1.f);
    const float turbulence = std::clamp(
        spectralTurbulence * anchors[anchor].turbulenceScale, 0.f, 1.f);
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
    const float anchorOutputGain = anchorOutputGains[anchor];
    const auto makeMode = [&](const std::size_t slot, const float frequency,
                              const float weight, const float phase) {
      const std::size_t index = anchor * CrashPacketModeCount + slot;
      const float safeFrequency = std::clamp(frequency, 20.f,
                                              .48f * sampleRate);
      result[index] = {
          safeFrequency,
          std::clamp(decay.At(safeFrequency), .01f, 30.f),
          anchorGain * weight,
          anchorOutputGain,
          phase,
          bandwidthErb * ErbBandwidth(safeFrequency) *
              (slot == 0 ? .35f : 1.f),
          static_cast<std::uint16_t>(anchor),
          turbulence * turbulence};
    };

    makeMode(0, centre, coreWeight, 0.f);
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
    const float baseBell = .25f + 1.15f * bellFocus;
    const float baseBow = .72f + .22f *
        std::abs(std::sin(3.14159265358979323846f * (2.7f * position + .13f)));
    const float baseEdge = .65f + .55f * std::sqrt(position);
    bell[mode] = baseBell;
    bow[mode] = baseBow;
    edge[mode] = baseEdge;
  }
}

ObservationModel<2>::Parameters Observation(
    const CrashCymbalFitParameters &fit) noexcept {
  ObservationModel<2>::Parameters result{};
  result[0].gain = std::clamp(fit.directGain, 0.f, 4.f);
  result[0].radiationEnabled = fit.directRadiationEnabled;
  result[0].radiation.lowCutHz =
      std::clamp(fit.directLowCutHz, 10.f, 1000.f);
  result[0].radiation.lowCutQ = .70710678f;
  result[0].radiation.colourFrequencyHz =
      std::clamp(fit.directColourFrequencyHz, 100.f, 18000.f);
  result[0].radiation.colourGainDb =
      std::clamp(fit.directColourGainDb, -18.f, 18.f);
  result[0].radiation.colourQ = .8f;
  result[0].radiation.highCutHz =
      std::clamp(fit.directHighCutHz, 1000.f, 22000.f);
  result[0].radiation.highCutQ = .70710678f;
  result[1].gain = std::clamp(fit.fieldGain, 0.f, 4.f);
  result[1].radiationEnabled = fit.bodyRadiationEnabled;
  result[1].radiation.lowCutHz =
      std::clamp(fit.bodyLowCutHz, 10.f, 1000.f);
  result[1].radiation.lowCutQ = .70710678f;
  result[1].radiation.colourFrequencyHz =
      std::clamp(fit.bodyColourFrequencyHz, 100.f, 18000.f);
  result[1].radiation.colourGainDb =
      std::clamp(fit.bodyColourGainDb, -18.f, 18.f);
  result[1].radiation.colourQ = .8f;
  result[1].radiation.highCutHz =
      std::clamp(fit.bodyHighCutHz, 1000.f, 22000.f);
  result[1].radiation.highCutQ = .70710678f;
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
  result.modalFieldControls.seed = ModalFieldSeed ^ 0x4649454cu;
  result.modalFieldControls.cascade = {
      std::clamp(fit.bloomRateOctavesPerSecond, 0.f, 32.f),
      std::clamp(fit.bloomEnergyDependence, 0.f, 1.f),
      std::clamp(fit.bloomPhaseDiffusion, 0.f, 1.f),
      ModalFieldSeed ^ 0x43415343u};
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
