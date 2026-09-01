#pragma once

#include "contact_exciter.hpp"
#include "coupled_resonator_network.hpp"
#include "dispersion_loop.hpp"
#include "observation_model.hpp"

#include <array>
#include <cstddef>

namespace tfdsp::percussion {

inline constexpr std::size_t CrashResonatorCount = 32;
inline constexpr std::size_t CrashResonatorBusCount = 4;
using CrashResonators = CoupledResonatorNetwork<CrashResonatorCount>;

struct CrashCymbalFitParameters {
  float resonanceTune{1.f};
  float lowDecayScale{1.f};
  float middleDecayScale{1.f};
  float highDecayScale{1.f};
  float resonatorCoupling{.38f};
  float resonatorShiftScale{1.f};
  float dispersionFeedback{.93f};
  float dispersionDrive{2.8f};
  float dispersionExcursionSamples{2.4f};
  float dispersionLowDecaySeconds{.9f};
  float dispersionMiddleDecaySeconds{.65f};
  float dispersionHighDecaySeconds{.42f};
  float directGain{.18f};
  float bodyGain{.12f};
  float bodyBypassGain{.06f};
  float outputGain{1.f};
  float colourFrequencyHz{5200.f};
  float colourGainDb{1.5f};
  float highCutHz{19000.f};
  float strengthGamma{1.15f};
};

struct CrashCymbalParameters {
  CrashResonators::Parameters resonators{};
  std::array<float, CrashResonatorCount> bellProjection{};
  std::array<float, CrashResonatorCount> bowProjection{};
  std::array<float, CrashResonatorCount> edgeProjection{};
  std::array<float, CrashResonatorBusCount> resonatorShiftHz{};
  DispersionLoopParameters dispersion{};
  ObservationModel<2>::Parameters observation{};
  CrashCymbalFitParameters fit{};
};

CrashCymbalParameters DefaultCrashCymbalParameters(
    float sampleRate, const CrashCymbalFitParameters &fit = {});

} // namespace tfdsp::percussion
