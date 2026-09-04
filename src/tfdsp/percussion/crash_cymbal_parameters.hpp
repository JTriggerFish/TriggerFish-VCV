#pragma once

#include "contact_exciter.hpp"
#include "metallic_plate_routing.hpp"
#include "observation_model.hpp"
#include "stochastic_modal_field.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace tfdsp::percussion {

inline constexpr std::size_t CrashSparseModeCount = 24;
inline constexpr std::size_t CrashPacketModeCount = 17;
inline constexpr std::size_t CrashModalFieldModeCount =
    CrashSparseModeCount * CrashPacketModeCount;
inline constexpr std::size_t CrashBodyDecayPointCount = 8;
inline constexpr std::size_t CrashBodyDecayInteriorPointCount =
    CrashBodyDecayPointCount - 2;
using CrashModalField = StochasticModalField<CrashModalFieldModeCount>;

struct CrashCymbalFitParameters {
  std::array<float, CrashSparseModeCount> sparseFrequencyHz{
      421.f, 522.f, 689.f, 1094.f, 1475.f, 2009.f,
      2138.f, 2573.f, 2753.f, 3589.f, 4428.f, 5707.f,
      6500.f, 7350.f, 8250.f, 9200.f, 10250.f, 11350.f,
      12000.f, 12750.f, 13500.f, 14100.f, 14600.f, 15000.f};
  std::array<float, CrashSparseModeCount> sparseAmplitude{
      .35f, .15f, .15f, .55f, .25f, .7f,
      .65f, .7f, .55f, .5f, .4f, .25f,
      .32f, .38f, .44f, .5f, .56f, .6f,
      .62f, .6f, .54f, .46f, .35f, .22f};
  std::array<float, CrashSparseModeCount> fieldTurbulenceScale{
      1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f,
      1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f,
      1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f};
  float sparseTune{1.f};
  std::array<float, CrashBodyDecayInteriorPointCount> bodyDecayFrequencyHz{
      500.f, 1500.f, 5000.f, 8000.f, 12000.f, 16000.f};
  std::array<float, CrashBodyDecayPointCount> bodyDecaySeconds{
      4.f, 2.965782f, 2.372841f, 1.782984f,
      1.585960f, 1.431738f, 1.330840f, 1.2f};
  std::array<bool, CrashBodyDecayInteriorPointCount> bodyDecayActive{
      false, false, false, false, false, false};
  float bodyTiltDbPerOctave{-1.f};
  float bodyTiltCentreHz{4000.f};
  // Intrinsic upward transport of energy already stored in the modal body.
  float bloomRateOctavesPerSecond{2.f};
  float bloomEnergyDependence{.7f};
  float bloomPhaseDiffusion{.7f};
  // Visible gain between the contact body port and the nonlinear modal field.
  float bodyExcitationGain{.05f};
  float fieldGain{.73824115f};
  float fieldTurbulence{.65f};
  float fieldTurbulenceSlopePerOctave{0.f};
  float fieldTurbulenceCentreHz{4000.f};
  float fieldPacketSpreadErb{6.f};
  float fieldPhaseBandwidthErb{1.f};
  float fieldExchange{.35f};
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
  float outputGain{1.f};
  bool directRadiationEnabled{true};
  float directLowCutHz{40.f};
  float directColourFrequencyHz{7200.f};
  float directColourGainDb{1.f};
  float directHighCutHz{20000.f};
  bool bodyRadiationEnabled{true};
  float bodyLowCutHz{40.f};
  float bodyColourFrequencyHz{7200.f};
  float bodyColourGainDb{.5f};
  float bodyHighCutHz{19000.f};
  float velocityBrightnessDbPerOctave{4.f};
};

struct CrashCymbalParameters {
  CrashModalField::Parameters modalField{};
  CrashModalField::Projection fieldBellProjection{};
  CrashModalField::Projection fieldBowProjection{};
  CrashModalField::Projection fieldEdgeProjection{};
  StochasticModalFieldControls modalFieldControls{};
  ObservationModel<2>::Parameters observation{};
  MetallicPlateRouting routing{};
  CrashCymbalFitParameters fit{};
};

struct CrashCymbalPreparedParameters {
  CrashCymbalParameters parameters{};
  CrashModalField::PreparedParameters modalField{};
  float sampleRate{48000.f};
};

CrashCymbalParameters DefaultCrashCymbalParameters(
    float sampleRate, const CrashCymbalFitParameters &fit = {});
CrashCymbalPreparedParameters PrepareCrashCymbalParameters(
    float sampleRate, const CrashCymbalParameters &parameters);

} // namespace tfdsp::percussion
