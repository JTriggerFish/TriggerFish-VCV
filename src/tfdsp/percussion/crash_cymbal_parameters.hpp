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
inline constexpr std::size_t CrashDenseModeCount = 2048;
inline constexpr std::size_t CrashBodyDecayPointCount = 5;
using CrashSparseModes = ModalBank<CrashSparseModeCount>;
using CrashDenseModes = ModalBank<CrashDenseModeCount>;

struct CrashCymbalFitParameters {
  std::array<float, CrashSparseModeCount> sparseFrequencyHz{
      421.f, 522.f, 689.f, 1094.f, 1475.f, 2009.f,
      2138.f, 2573.f, 2753.f, 3589.f, 4428.f, 5707.f};
  std::array<float, CrashSparseModeCount> sparseDecayRatio{
      .7f, .7f, .7f, 1.25f, 1.f, 1.f,
      .95f, .85f, .8f, .7f, .7f, .7f};
  std::array<float, CrashSparseModeCount> sparseAmplitude{
      .35f, .15f, .15f, .55f, .25f, .7f,
      .65f, .7f, .55f, .5f, .4f, .25f};
  std::array<float, CrashSparseModeCount> sparsePhaseRadians{};
  float sparseTune{1.f};
  std::array<float, CrashBodyDecayPointCount> bodyDecayFrequencyHz{
      200.f, 500.f, 1500.f, 5000.f, 15000.f};
  std::array<float, CrashBodyDecayPointCount> bodyDecaySeconds{
      4.f, 4.f, 3.8f, 2.3f, 1.2f};
  float denseMinimumFrequencyHz{180.f};
  float denseMaximumFrequencyHz{18000.f};
  float denseFrequencyWarp{1.f};
  float denseSpacingJitter{.82f};
  // Active cloud density in 2048-mode banks. Values above one progressively
  // add an independent extension bank without changing the primary cloud.
  float denseModeDensity{1.f};
  float denseDecaySpreadOctaves{.15f};
  float denseTiltDbPerOctave{-1.f};
  std::array<float,
             StatisticalModalCloudParameters::GainEnvelopePointCount>
      denseGainEnvelopeDb{
          4.f, 1.375f, -1.25f, -3.875f, -6.5f, -7.25f,
          -5.5f, -3.75f, -2.f, -.25f, 1.125f, 2.4375f,
          3.75f, 5.0625f, 6.125f, 6.5625f, 7.f, 7.4375f,
          7.875f, 7.6875f, 7.25f,
          6.8125f, 6.375f, 5.9375f, 5.5f, 5.0625f, 4.625f,
          4.1875f, 3.75f, 3.3125f, 2.875f, 2.4375f, 2.f};
  float denseGainSpreadDb{2.f};
  std::uint32_t denseModeSeed{0x43524153u};
  std::array<float, 3> turbulenceFrequencyHz{350.f, 2200.f, 10000.f};
  std::array<float, 3> turbulenceGain{.14f, .1f, .06f};
  float turbulencePersistence{1.f};
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
  // Resolved modes hear both the immediate strike and the dispersed body.
  float sparseBloomGain{1.f};
  float bodyBypassGain{.06f};
  float outputGain{1.f};
  bool directRadiationEnabled{true};
  float directLowCutHz{40.f};
  float directLowCutQ{.707f};
  float directColourFrequencyHz{7200.f};
  float directColourGainDb{1.f};
  float directColourQ{.8f};
  float directHighCutHz{20000.f};
  float directHighCutQ{.707f};
  bool sparseRadiationEnabled{true};
  float sparseLowCutHz{40.f};
  float sparseLowCutQ{.707f};
  float colourFrequencyHz{5200.f};
  float colourGainDb{1.5f};
  float sparseColourQ{.8f};
  float highCutHz{19000.f};
  float sparseHighCutQ{.707f};
  bool denseRadiationEnabled{true};
  float denseLowCutHz{40.f};
  float denseLowCutQ{.707f};
  float denseColourFrequencyHz{7200.f};
  float denseColourGainDb{.5f};
  float denseColourQ{.8f};
  float denseHighCutHz{19000.f};
  float denseHighCutQ{.707f};
  float strengthGamma{1.15f};
  float bodyStrengthGamma{.8f};
  float velocityBrightnessDbPerOctave{4.f};
};

struct CrashCymbalParameters {
  CrashSparseModes::Parameters sparseModes{};
  CrashDenseModes::Parameters denseModes{};
  CrashDenseModes::Parameters denseExtensionModes{};
  CrashSparseModes::Projection sparseBellProjection{};
  CrashSparseModes::Projection sparseBowProjection{};
  CrashSparseModes::Projection sparseEdgeProjection{};
  CrashDenseModes::Projection denseBellProjection{};
  CrashDenseModes::Projection denseBowProjection{};
  CrashDenseModes::Projection denseEdgeProjection{};
  CrashDenseModes::Projection denseExtensionBellProjection{};
  CrashDenseModes::Projection denseExtensionBowProjection{};
  CrashDenseModes::Projection denseExtensionEdgeProjection{};
  DispersionLoopParameters dispersion{};
  TurbulentResidualParameters turbulence{};
  ObservationModel<4>::Parameters observation{};
  CrashCymbalFitParameters fit{};
};

CrashCymbalParameters DefaultCrashCymbalParameters(
    float sampleRate, const CrashCymbalFitParameters &fit = {});

} // namespace tfdsp::percussion
