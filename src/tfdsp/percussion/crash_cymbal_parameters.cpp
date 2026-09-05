#include "crash_cymbal_parameters.hpp"

#include "erb_scale.hpp"
#include "modal_packet_allocator.hpp"

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
    const float maximumFrequency = std::min(
        CrashModalMaximumFrequencyHz, .48f * sampleRate);
    points_[count_++] = {ErbRate(CrashModalMinimumFrequencyHz),
        std::log(std::clamp(
            Positive(fit.bodyDecaySeconds.front(), 1.f), .02f, 30.f))};
    for (std::size_t index = 0; index < CrashBodyDecayInteriorPointCount;
         ++index) {
      if (!fit.bodyDecayActive[index]) continue;
      const float frequency = std::clamp(tfdsp::FiniteNormalOrZero(
          fit.bodyDecayFrequencyHz[index]), CrashModalMinimumFrequencyHz,
          maximumFrequency);
      const float seconds = std::clamp(
          Positive(fit.bodyDecaySeconds[index + 1], 1.f), .02f, 30.f);
      points_[count_++] = {ErbRate(frequency), std::log(seconds)};
    }
    points_[count_++] = {ErbRate(maximumFrequency), std::log(std::clamp(
        Positive(fit.bodyDecaySeconds.back(), 1.f), .02f, 30.f))};
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

struct ModalAnchor {
  float frequencyHz{};
  float amplitude{};
  float turbulence{};
  float spreadErb{};
};

float NestedRadius(std::size_t index) noexcept {
  if (index == 0) return 1.f;
  float result = 0.f;
  float place = .5f;
  while (index > 0) {
    result += place * static_cast<float>(index & 1u);
    index >>= 1u;
    place *= .5f;
  }
  return result;
}

float ExcitationTiltGain(const float frequencyHz, const float centreHz,
                         const float tiltDbPerOctave) noexcept {
  const float ratio = std::max(frequencyHz, 20.f) /
      std::max(centreHz, 20.f);
  // Smooth shelving knee: flat below the centre and asymptotically equal to
  // tiltDbPerOctave above it. Unlike a pivoted power law, the centre remains
  // meaningful after the complete excitation vector is energy-normalized.
  constexpr float FortyLog10Two = 12.041199826559248f;
  const float exponent = std::clamp(tiltDbPerOctave, -72.f, 24.f) /
      FortyLog10Two;
  return std::pow(1.f + ratio * ratio, exponent);
}

std::size_t BuildActiveAnchors(
    const CrashCymbalFitParameters &fit,
    std::array<ModalAnchor, CrashModalAnchorCapacity> &anchors) noexcept {
  std::size_t count = 0;
  for (std::size_t index = 0; index < anchors.size(); ++index) {
    if (!(fit.sparseAmplitude[index] > 0.f)) continue;
    anchors[count++] = {
        fit.sparseFrequencyHz[index], fit.sparseAmplitude[index],
        fit.fieldTurbulenceScale[index], 0.f};
  }
  std::sort(anchors.begin(), anchors.begin() + count,
            [](const auto &left, const auto &right) {
              return left.frequencyHz < right.frequencyHz;
            });
  return count;
}

CrashModalField::Parameters ModalField(
    const float sampleRate, const CrashCymbalFitParameters &fit,
    const BodyDecayEnvelope &decay) noexcept {
  CrashModalField::Parameters result{};
  for (auto &mode : result) {
    mode.inputGain = 0.f;
    mode.outputGain = 0.f;
  }
  std::array<ModalAnchor, CrashModalAnchorCapacity> anchors{};
  const std::size_t anchorCount = BuildActiveAnchors(fit, anchors);
  if (anchorCount == 0) return result;
  DeterministicRandom random;
  random.Seed(ModalFieldSeed ^ 0x4649454cu);
  const float globalTurbulence = std::clamp(fit.fieldTurbulence, 0.f, 1.f);
  const float turbulenceSlope = std::clamp(
      fit.fieldTurbulenceSlopePerOctave, -1.f, 1.f);
  const float turbulenceCentre = std::clamp(
      Positive(fit.fieldTurbulenceCentreHz, 4000.f), 20.f, .48f * sampleRate);
  std::array<float, CrashModalAnchorCapacity> anchorGains{};
  std::array<float, CrashModalAnchorCapacity> anchorOutputGains{};
  std::array<ModalPacketRequest, CrashModalAnchorCapacity> requests{};
  float anchorSquaredGain = 0.f;
  const float excitationCentreHz = std::clamp(
      Positive(fit.bodyExcitationCentreHz, 1000.f), 40.f,
      .48f * sampleRate);
  for (std::size_t anchor = 0; anchor < anchorCount; ++anchor) {
    const float frequency = Positive(anchors[anchor].frequencyHz, 1000.f);
    const float tilt = ExcitationTiltGain(
        frequency, excitationCentreHz, fit.bodyTiltDbPerOctave);
    const float level = std::clamp(
        Positive(anchors[anchor].amplitude, 0.f), 0.f, 8.f);
    // Painted levels describe observation prominence. The excitation curve is
    // deliberately independent: a unit-norm spatial input distributes the
    // contact impulse. Delivered energy also depends on time and phase.
    anchorGains[anchor] = tilt;
    anchorOutputGains[anchor] = level;
    anchorSquaredGain += tilt * tilt;
    const float spectralTurbulence = std::clamp(
        globalTurbulence + turbulenceSlope * std::log2(
            std::max(frequency, 20.f) / turbulenceCentre), 0.f, 1.f);
    anchors[anchor].turbulence = std::clamp(
        spectralTurbulence * anchors[anchor].turbulence, 0.f, 1.f);
    anchors[anchor].spreadErb = anchors[anchor].turbulence * std::clamp(
        fit.fieldPacketSpreadErb, 0.f, 12.f);
    requests[anchor] = {
        ErbRate(frequency), anchors[anchor].spreadErb, true};
  }
  const float anchorNormalization =
      1.f / std::sqrt(std::max(anchorSquaredGain, 1.e-12f));
  const auto allocation = AllocateModalPackets(
      requests, CrashModalFieldModeCount, fit.fieldSatelliteDensity);
  constexpr float Pi = 3.14159265358979323846f;
  std::size_t modeIndex = 0;
  for (std::size_t anchor = 0; anchor < anchorCount; ++anchor) {
    const float turbulence = anchors[anchor].turbulence;
    const std::size_t pairCount = allocation.sidebandPairs[anchor];
    const float diffuseEnergy = pairCount > 0
        ? .9f * turbulence * turbulence : 0.f;
    const float coreWeight = std::sqrt(1.f - diffuseEnergy);
    const float satelliteWeight = pairCount > 0 ? std::sqrt(
        diffuseEnergy / static_cast<float>(2 * pairCount)) : 0.f;
    const float spreadErb = anchors[anchor].spreadErb;
    const float bandwidthErb = turbulence * turbulence * std::clamp(
        fit.fieldPhaseBandwidthErb, 0.f, 4.f);
    const float centre = std::clamp(
        Positive(anchors[anchor].frequencyHz, 1000.f) *
            std::clamp(Positive(fit.sparseTune, 1.f), .5f, 2.f),
        20.f, .48f * sampleRate);
    const float anchorGain = anchorNormalization * anchorGains[anchor];
    // Bars are actual observation amplitudes, not a second input budget.
    // Replicating a packet's observation over its states preserves expected
    // incoherent power because its excitation weights have unit squared sum.
    // Do not divide again by handle count or renormalize a sounding tail.
    const float anchorOutputGain = anchorOutputGains[anchor];
    const auto makeMode = [&](const float frequency, const float weight,
                              const float phase, const float bandwidthScale) {
      const float safeFrequency = std::clamp(frequency, 20.f,
                                              .48f * sampleRate);
      result[modeIndex++] = {
          safeFrequency,
          std::clamp(decay.At(safeFrequency), .02f, 30.f),
          anchorGain * weight,
          anchorOutputGain,
          phase,
          bandwidthErb * ErbBandwidth(safeFrequency) * bandwidthScale,
          static_cast<std::uint16_t>(anchor),
          turbulence * turbulence};
    };

    makeMode(centre, coreWeight, 0.f, .35f);
    const float centreErb = ErbRate(centre);
    const std::size_t packetBegin = modeIndex - 1;
    for (std::size_t pair = 0; pair < pairCount; ++pair) {
      const float jitter = .92f + .08f * random.Uniform();
      const float offset = spreadErb * NestedRadius(pair) * jitter;
      const float low = InverseErbRate(std::max(ErbRate(20.f),
                                                centreErb - offset));
      const float high = InverseErbRate(std::min(ErbRate(.48f * sampleRate),
                                                 centreErb + offset));
      makeMode(low, satelliteWeight, Pi * random.Bipolar(), 1.f);
      makeMode(high, satelliteWeight, Pi * random.Bipolar(), 1.f);
    }
    auto first = result.begin() + packetBegin;
    std::sort(first, result.begin() + modeIndex,
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
    const float position = std::clamp(
        std::log2(std::max(frequency, 20.f) / 20.f) /
            std::log2(22000.f / 20.f),
        0.f, 1.f);
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
      std::clamp(fit.bloomEnergyAcceleration, 0.f, 1.f),
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
