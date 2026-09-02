#pragma once

#include "contact_exciter.hpp"
#include "dispersion_loop.hpp"
#include "modal_bank.hpp"
#include "observation_model.hpp"
#include "statistical_modal_cloud.hpp"
#include "turbulent_residual.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace tfdsp::percussion {

inline constexpr std::size_t CrashSparseModeCount = 12;
inline constexpr std::size_t CrashDenseModeCount = 512;
using CrashSparseModes = ModalBank<CrashSparseModeCount>;
using CrashDenseModes = ModalBank<CrashDenseModeCount>;

struct CrashCymbalFitParameters {
  std::array<float, CrashSparseModeCount> sparseFrequencyHz{
      522.f, 689.f, 1094.f, 1475.f, 2009.f, 2138.f,
      2573.f, 2753.f, 3589.f, 3923.f, 4428.f, 5707.f};
  std::array<float, CrashSparseModeCount> sparseDecaySeconds{
      5.5f, 6.5f, 4.8f, 3.8f, 3.2f, 3.f,
      2.5f, 2.2f, 1.6f, 1.35f, 1.05f, .75f};
  std::array<float, CrashSparseModeCount> sparseAmplitude{
      .25f, .32f, .55f, .25f, .8f, .65f,
      .75f, .62f, .55f, .48f, .4f, .25f};
  std::array<float, CrashSparseModeCount> sparsePhaseRadians{};
  float sparseTune{1.f};
  float sparseDecayScale{1.f};
  float denseMinimumFrequencyHz{700.f};
  float denseMaximumFrequencyHz{18000.f};
  float denseFrequencyWarp{1.f};
  float denseSpacingJitter{.82f};
  float denseLowDecaySeconds{3.2f};
  float denseHighDecaySeconds{.22f};
  float denseDecayCurve{.75f};
  std::array<float,
             StatisticalModalCloudParameters::DecayEnvelopePointCount>
      denseDecayEnvelopeOctaves{};
  float denseDecaySpreadOctaves{.4f};
  float denseTiltDbPerOctave{-1.f};
  std::array<float,
             StatisticalModalCloudParameters::GainEnvelopePointCount>
      denseGainEnvelopeDb{};
  float denseGainSpreadDb{4.5f};
  std::uint32_t denseModeSeed{0x43524153u};
  float turbulenceLowGain{0.f};
  float turbulenceMiddleGain{0.f};
  float turbulenceHighGain{0.f};
  float turbulenceLowDecaySeconds{.8f};
  float turbulenceMiddleDecaySeconds{.55f};
  float turbulenceHighDecaySeconds{.3f};
  float dispersionFeedback{.9965f};
  float dispersionDrive{2.8f};
  float dispersionExcursionSamples{2.4f};
  float dispersionLowDecaySeconds{.9f};
  float dispersionMiddleDecaySeconds{.65f};
  float dispersionHighDecaySeconds{.42f};
  float contactDurationScale{1.f};
  float contactPulseGain{1.f};
  float contactChirpGain{1.f};
  float contactChirpFrequencyScale{1.f};
  float contactNoiseGain{1.f};
  float contactNoiseDurationScale{1.f};
  float contactNoiseTiltDb{0.f};
  float contactMicroGain{1.f};
  float contactMicroDurationScale{1.f};
  float contactMicroDensityScale{1.f};
  float directGain{.18f};
  float sparseGain{.35f};
  float denseGain{.65f};
  float sparseBloomGain{0.f};
  float bodyBypassGain{.06f};
  float outputGain{1.f};
  float colourFrequencyHz{5200.f};
  float colourGainDb{1.5f};
  float highCutHz{19000.f};
  float strengthGamma{1.15f};
  float bodyStrengthGamma{.8f};
  float denseStrengthGamma{.8f};
  float denseVelocityLossNepersPerSecond{0.f};
};

struct CrashCymbalParameters {
  CrashSparseModes::Parameters sparseModes{};
  CrashDenseModes::Parameters denseModes{};
  CrashSparseModes::Projection sparseBellProjection{};
  CrashSparseModes::Projection sparseBowProjection{};
  CrashSparseModes::Projection sparseEdgeProjection{};
  CrashDenseModes::Projection denseBellProjection{};
  CrashDenseModes::Projection denseBowProjection{};
  CrashDenseModes::Projection denseEdgeProjection{};
  DispersionLoopParameters dispersion{};
  TurbulentResidualParameters turbulence{};
  ObservationModel<3>::Parameters observation{};
  CrashCymbalFitParameters fit{};
};

CrashCymbalParameters DefaultCrashCymbalParameters(
    float sampleRate, const CrashCymbalFitParameters &fit = {});

} // namespace tfdsp::percussion
