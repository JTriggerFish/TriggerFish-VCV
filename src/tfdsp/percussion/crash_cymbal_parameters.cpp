#include "crash_cymbal_parameters.hpp"

#include "erb_scale.hpp"
#include "modal_packet_allocator.hpp"
#include "modal_spectral_diffusion.hpp"
#include "turbulence_profile.hpp"

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
    points_[count_++] = {ErbRate(CrashDecayMinimumFrequencyHz),
        std::log(std::clamp(
            Positive(fit.bodyDecaySeconds.front(), 1.f), .02f, 30.f))};
    for (std::size_t index = 0; index < CrashBodyDecayInteriorPointCount;
         ++index) {
      if (!fit.bodyDecayActive[index]) continue;
      const float frequency = std::clamp(tfdsp::FiniteNormalOrZero(
          fit.bodyDecayFrequencyHz[index]), CrashDecayMinimumFrequencyHz,
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
  float diffuseEnergy{};
  float exchangeAmount{};
  float allocationWeight{1.f};
};

float ExcitationTiltGain(const float frequencyHz, const float centreHz,
                         const float tiltDbPerOctave) noexcept {
  const float ratio = std::max(frequencyHz, CrashModalMinimumFrequencyHz) /
      std::max(centreHz, CrashModalMinimumFrequencyHz);
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
        Positive(fit.sparseFrequencyHz[index], 1000.f), fit.sparseAmplitude[index],
        fit.fieldTurbulenceScale[index], 0.f, 0.f, 0.f,
        fit.fieldAllocationWeight[index]};
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
  const float turbulenceSlope = std::clamp(
      fit.fieldTurbulenceSlopePerOctave, -1.f, 1.f);
  constexpr float turbulenceCentre = 1000.f; // Noisiness is defined at 1 kHz.
  std::array<float, CrashModalAnchorCapacity> anchorGains{};
  std::array<float, CrashModalAnchorCapacity> anchorOutputGains{};
  std::array<ModalPacketRequest, CrashModalAnchorCapacity> requests{};
  double anchorSquaredGain = 0.0;
  const float excitationCentreHz = std::clamp(
      Positive(fit.bodyExcitationCentreHz, 1000.f), CrashModalMinimumFrequencyHz,
      .48f * sampleRate);
  // The diffusion state is energy per frequency cell. Sample the excitation
  // density with those same quadrature masses; equal energy per handle would
  // invent large density spikes wherever the editor has closely spaced bars.
  ModalSpectralDiffusion<CrashModalAnchorCapacity> excitationGrid;
  if (fit.bloomSpectralDiffusion) {
    std::array<float, CrashModalAnchorCapacity> centres{}, weights{};
    for (std::size_t i = 0; i < anchorCount; ++i) {
      centres[i] = std::clamp(anchors[i].frequencyHz *
          std::clamp(Positive(fit.sparseTune, 1.f), .5f, 2.f),
          CrashModalMinimumFrequencyHz, .48f * sampleRate);
      weights[i] = 1.f;
    }
    excitationGrid.Prepare(centres, weights, anchorCount, sampleRate);
  }
  for (std::size_t anchor = 0; anchor < anchorCount; ++anchor) {
    const float frequency = Positive(anchors[anchor].frequencyHz, 1000.f);
    const float tilt = ExcitationTiltGain(
        frequency, excitationCentreHz, fit.bodyTiltDbPerOctave);
    const float level = std::clamp(
        Positive(anchors[anchor].amplitude, 0.f), 0.f, 8.f);
    // Painted levels describe observation prominence. The excitation curve is
    // deliberately independent: a unit-norm spatial input distributes the
    // contact impulse. Delivered energy also depends on time and phase.
    const float cellAmplitude = fit.bloomSpectralDiffusion
        ? float(std::sqrt(excitationGrid.ExcitationCellWeight(anchor))) : 1.f;
    anchorGains[anchor] = tilt * cellAmplitude;
    anchorOutputGains[anchor] = level;
    anchorSquaredGain += double(anchorGains[anchor]) * anchorGains[anchor];
    const auto response = EvaluateTurbulence(frequency, fit.fieldTurbulence,
        turbulenceSlope, turbulenceCentre, anchors[anchor].turbulence,
        fit.fieldRelaxedTurbulence);
    anchors[anchor].turbulence = response.intensity;
    anchors[anchor].diffuseEnergy = response.diffuseEnergy;
    anchors[anchor].exchangeAmount = response.exchangeAmount;
    anchors[anchor].spreadErb = anchors[anchor].turbulence * std::clamp(
        fit.fieldPacketSpreadErb, 0.f, 12.f);
    requests[anchor] = {
        ErbRate(frequency), anchors[anchor].spreadErb, true,
        anchors[anchor].allocationWeight};
  }
  // A very dark shelf must redistribute energy, not secretly turn the body
  // down when all its unnormalized weights fall below an arbitrary floor.
  const double anchorNormalization = anchorSquaredGain > 0
      ? 1.0 / std::sqrt(anchorSquaredGain) : 0.0;
  const bool pairedRing = HasPairedRing(fit.fieldDistribution, fit.fieldDoubletSplitHz);
  const auto allocation = AllocateModalPackets(
      requests, CrashModalFieldModeCount, fit.fieldSatelliteDensity, pairedRing);
  constexpr float Pi = 3.14159265358979323846f;
  std::size_t modeIndex = 0;
  for (std::size_t anchor = 0; anchor < anchorCount; ++anchor) {
    // Another packet gaining states must not randomize this packet's phases.
    random.Seed(ModalFieldSeed ^ (0x9e3779b9u * std::uint32_t(anchor + 1)));
    const float turbulence = anchors[anchor].turbulence;
    const std::size_t pairCount = allocation.sidebandPairs[anchor];
    const float diffuseEnergy = pairCount > 0
        ? anchors[anchor].diffuseEnergy : 0.f;
    const float coreWeight = std::sqrt(1.f - diffuseEnergy);
    const float satelliteWeight = pairCount > 0 ? std::sqrt(
        diffuseEnergy / static_cast<float>(2 * pairCount)) : 0.f;
    const float spreadErb = anchors[anchor].spreadErb;
    const float bandwidthErb = turbulence * turbulence * std::clamp(
        fit.fieldPhaseBandwidthErb, 0.f, 4.f);
    const float centre = std::clamp(
        Positive(anchors[anchor].frequencyHz, 1000.f) *
            std::clamp(Positive(fit.sparseTune, 1.f), .5f, 2.f),
        CrashModalMinimumFrequencyHz, .48f * sampleRate);
    const float anchorGain = anchorNormalization * anchorGains[anchor];
    // Bars are actual observation amplitudes, not a second input budget.
    // Replicating a packet's observation over its states preserves expected
    // incoherent power because its excitation weights have unit squared sum.
    // Do not divide again by handle count or renormalize a sounding tail.
    const float anchorOutputGain = anchorOutputGains[anchor];
    const auto makeMode = [&](const float frequency, const float weight,
                              const float phase, const float bandwidthScale) {
      const float safeFrequency = std::clamp(frequency, CrashModalMinimumFrequencyHz,
                                              .48f * sampleRate);
      // A preparation-only colour control: phase coherence can vary without
      // changing packet allocation, excitation energy or ordinary damping.
      const float blurTilt = std::clamp(
          tfdsp::FiniteNormalOrZero(fit.fieldPhaseTilt), -2.f, 2.f);
      const float blurColour = std::pow(safeFrequency / 1000.f, blurTilt);
      result[modeIndex++] = {
          safeFrequency,
          std::clamp(decay.At(safeFrequency), .02f, 30.f),
          anchorGain * weight,
          anchorOutputGain,
          phase,
          bandwidthErb * ErbBandwidth(safeFrequency) * bandwidthScale * blurColour,
          static_cast<std::uint16_t>(anchor),
          anchors[anchor].exchangeAmount,
          centre};
    };

    const std::size_t packetBegin = modeIndex;
    if (pairedRing) {
      const float depth = std::clamp(fit.fieldBeatDepth, 0.f, 1.f);
      const float rate = RingBeatRate(centre, fit.fieldDoubletSplitHz, fit.fieldBeatRateTilt);
      const float halfSplit = depth > 0.f ? RingHalfSplit(centre, rate, .48f * sampleRate) : 0.f;
      // Orthogonal launch: complex coefficients sum to one and their squared
      // norms sum to one at every depth. Zero depth is one unsplit oscillator;
      // reserve its unused partner to avoid reallocating the surrounding cloud.
      const float normalization = 1.f / std::sqrt(1.f + depth * depth);
      const float angle = std::atan(depth);
      makeMode(centre - halfSplit, coreWeight * depth * normalization, Pi * .5f - angle, .35f);
      makeMode(centre + halfSplit, coreWeight * normalization, -angle, .35f);
    } else {
      makeMode(centre, coreWeight, 0.f, .35f);
    }
    for (std::size_t pair = 0; pair < pairCount; ++pair) {
      const float jitter = .92f + .08f * random.Uniform();
      const bool doublets = fit.fieldDistribution == ModalPacketDistribution::Doublets;
      const float pairRate = doublets ? RingBeatRate(centre,
          fit.fieldDoubletSplitHz, fit.fieldBeatRateTilt) : fit.fieldDoubletSplitHz;
      const float weight = satelliteWeight * (doublets
          ? DoubletWeightScale(pair, pairCount, fit.fieldBeatDepth) : 1.f);
      const auto sideFrequency = [&](float side) {
        return PacketSideFrequency(centre, spreadErb, pair, side,
            fit.fieldDistribution, pairRate, jitter, .48f * sampleRate);
      };
      const float low = sideFrequency(-1.f), high = sideFrequency(1.f);
      makeMode(low, weight, Pi * random.Bipolar(), 1.f);
      makeMode(high, weight, Pi * random.Bipolar(), 1.f);
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

RadiationFilterParameters OutputEq(
    const CrashCymbalFitParameters &fit) noexcept {
  RadiationFilterParameters result{};
  result.lowCutHz =
      std::clamp(fit.outputLowCutHz, 10.f, 1000.f);
  result.lowCutQ = .70710678f;
  result.colourFrequencyHz =
      std::clamp(fit.outputColourFrequencyHz, 100.f, 18000.f);
  result.colourGainDb =
      std::clamp(fit.outputColourGainDb, -18.f, 18.f);
  result.colourQ = .8f;
  result.highCutHz =
      std::clamp(fit.outputHighCutHz, 1000.f, 22000.f);
  result.highCutQ = .70710678f;
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
  result.modalFieldControls.driftDepthHz = fit.fieldWanderDepthHz;
  result.modalFieldControls.driftKnotsPerSecond = fit.fieldWanderKnotsPerSecond;
  result.modalFieldControls.motion = fit.fieldMotion;
  result.modalFieldControls.cascade = {
      std::clamp(fit.bloomRateOctavesPerSecond, 0.f, 32.f),
      std::clamp(fit.bloomEnergyAcceleration, 0.f, 1.f),
      std::clamp(fit.bloomPhaseDiffusion, 0.f, 1.f),
      ModalFieldSeed ^ 0x43415343u, fit.bloomSpectralDiffusion,
      std::clamp(fit.bloomEnergySensitivity, 0.f, 2.f)};
  result.outputEq = OutputEq(fit);
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
