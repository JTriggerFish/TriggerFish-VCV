#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace tfdsp {
inline constexpr std::size_t EarlyReflectionOutputCount = 2;
inline constexpr std::size_t EarlyReflectionBandCount = 4;
inline constexpr std::size_t EarlyReflectionSurfaceCount = 4;
inline constexpr std::size_t EarlyReflectionMaximumSources = 8;
inline constexpr double EarlyReflectionPi =
    3.141592653589793238462643383279502884;

using EarlyReflectionVector = std::array<double, 3>;

enum class EarlyReflectionRoomFamily { Room, Chamber, Hall };

struct EarlyReflectionConfig {
  double sampleRate{48000.0};
  double speedOfSound{343.0};
  // Hard ceiling for direct-relative ER generation. The adaptive mixing-time
  // predictor normally selects a considerably shorter horizon.
  double responseTimeSeconds{0.150};
  double materialFilterTailSeconds{0.020};
  double minimumHandoffSeconds{0.020};
  double handoffOverlapSeconds{0.020};
  double analysisSafetySeconds{0.020};
  double analysisWindowSeconds{0.020};
  double analysisHopSeconds{0.0025};
  double analysisStableSeconds{0.015};
  double referenceDistanceMetres{1.0};
  std::array<double, 3> crossoverHz{250.0, 1000.0, 4000.0};

  void Validate() const {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument("sample rate must be positive and finite");
    if (!std::isfinite(speedOfSound) || speedOfSound <= 0.0)
      throw std::invalid_argument("speed of sound must be positive and finite");
    if (!std::isfinite(responseTimeSeconds) || responseTimeSeconds <= 0.0)
      throw std::invalid_argument(
          "ER response time must be positive and finite");
    if (!std::isfinite(materialFilterTailSeconds) ||
        materialFilterTailSeconds <= 0.0)
      throw std::invalid_argument(
          "material-filter tail must be positive and finite");
    if (!std::isfinite(minimumHandoffSeconds) || minimumHandoffSeconds <= 0.0 ||
        minimumHandoffSeconds > responseTimeSeconds)
      throw std::invalid_argument(
          "minimum handoff must lie inside the ER horizon");
    if (!std::isfinite(handoffOverlapSeconds) || handoffOverlapSeconds <= 0.0)
      throw std::invalid_argument(
          "handoff overlap must be positive and finite");
    if (!std::isfinite(analysisSafetySeconds) || analysisSafetySeconds < 0.0)
      throw std::invalid_argument(
          "analysis safety margin must be finite and non-negative");
    if (!std::isfinite(analysisWindowSeconds) || analysisWindowSeconds <= 0.0 ||
        !std::isfinite(analysisHopSeconds) || analysisHopSeconds <= 0.0 ||
        !std::isfinite(analysisStableSeconds) || analysisStableSeconds <= 0.0)
      throw std::invalid_argument(
          "mixing analysis durations must be positive and finite");
    if (!std::isfinite(referenceDistanceMetres) ||
        referenceDistanceMetres <= 0.0)
      throw std::invalid_argument(
          "reference distance must be positive and finite");
    double previous = 0.0;
    for (const double crossover : crossoverHz) {
      if (!std::isfinite(crossover) || crossover <= previous ||
          crossover >= 0.5 * sampleRate)
        throw std::invalid_argument(
            "ER crossovers must increase and remain below Nyquist");
      previous = crossover;
    }
  }
};

struct EarlyReflectionRoom {
  EarlyReflectionVector dimensionsMetres{};
  EarlyReflectionVector listenerPositionMetres{};

  void Validate() const {
    for (std::size_t axis = 0; axis < 3; ++axis) {
      const double dimension = dimensionsMetres[axis];
      if (!std::isfinite(dimension) || dimension < 1.0)
        throw std::invalid_argument(
            "each room dimension must be finite and at least one metre");
      const double listener = listenerPositionMetres[axis];
      if (!std::isfinite(listener) || listener <= 0.0 || listener >= dimension)
        throw std::invalid_argument("listener must lie inside the room");
    }
  }
};

struct EarlyReflectionSource {
  EarlyReflectionVector positionMetres{};
};

struct EarlyReflectionMixingPrediction {
  double roomVolumeCubicMetres{};
  double roomSurfaceSquareMetres{};
  double averageSeconds{};
  double conservativeSeconds{};
  double generationHorizonSeconds{};
};

struct EarlyReflectionSourceHandoff {
  double directPropagationSeconds{};
  double predictedStartSeconds{};
  double predictedEndSeconds{};
  double detectedMixingSeconds{};
  double finalStartSeconds{};
  double finalEndSeconds{};
  double normalizedEchoDensity{};
  double excessKurtosis{};
  // Untapered energy descriptors at the ER/late handoff. They let the room
  // engine set a geometry-dependent late send without fitting an external IR.
  std::array<double, EarlyReflectionBandCount> bandHandoffPower{};
  std::array<double, EarlyReflectionBandCount> bandEarlyEnergy{};
  double broadbandHandoffPower{};
  double broadbandEarlyEnergy{};
  std::size_t imagePathCount{};
  std::size_t analysisPathCount{};
  bool detectedFromResponse{};
};

struct EarlyReflectionMaterials {
  // Pressure-amplitude coefficients. Rows are floor, ceiling, side walls,
  // and front/rear walls; columns are four complementary frequency bands.
  std::array<std::array<double, EarlyReflectionBandCount>,
             EarlyReflectionSurfaceCount>
      reflectionAmplitudes{};
  std::array<double, EarlyReflectionBandCount> airAbsorptionDbPerMetre{};

  void Validate() const {
    for (const auto &surface : reflectionAmplitudes)
      for (const double amplitude : surface)
        if (!std::isfinite(amplitude) || amplitude < 0.0 || amplitude > 1.0)
          throw std::invalid_argument(
              "surface reflection amplitudes must be in [0, 1]");
    for (const double absorption : airAbsorptionDbPerMetre)
      if (!std::isfinite(absorption) || absorption < 0.0)
        throw std::invalid_argument(
            "air absorption must be finite and non-negative");
  }
};

struct EarlyReflectionPath {
  std::array<int, 3> imageIndex{};
  std::size_t sourceIndex{};
  std::size_t reflectionOrder{};
  std::array<std::uint16_t, EarlyReflectionSurfaceCount> surfaceCounts{};
  EarlyReflectionVector imagePositionMetres{};
  EarlyReflectionVector arrivalDirection{};
  double distanceMetres{};
  double propagationSeconds{};
  double excessDelaySeconds{};
  std::array<double, EarlyReflectionOutputCount> outputGains{};
  std::array<double, EarlyReflectionBandCount> bandGains{};
};

struct EarlyReflectionImpulseResponse {
  double sampleRate{};
  std::size_t sourceCount{};
  std::size_t imagePathCount{};
  std::size_t analysisPathCount{};
  EarlyReflectionMixingPrediction mixingPrediction{};
  std::vector<EarlyReflectionSourceHandoff> sourceHandoffs{};
  // Output-major layout: kernel(output, source) = output * sourceCount +
  // source.
  std::vector<std::vector<double>> kernels{};

  std::size_t Size() const noexcept {
    return kernels.empty() ? 0 : kernels.front().size();
  }

  std::size_t KernelIndex(const std::size_t output,
                          const std::size_t source) const noexcept {
    return output * sourceCount + source;
  }

  void Validate() const {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument("ER FIR sample rate is invalid");
    if (sourceCount == 0 || sourceCount > EarlyReflectionMaximumSources)
      throw std::invalid_argument("ER FIR source count is invalid");
    if (kernels.size() != EarlyReflectionOutputCount * sourceCount)
      throw std::invalid_argument("ER FIR must contain two kernels per source");
    if (sourceHandoffs.size() != sourceCount)
      throw std::invalid_argument("ER FIR must contain one handoff per source");
    if (!std::isfinite(mixingPrediction.generationHorizonSeconds) ||
        mixingPrediction.generationHorizonSeconds <= 0.0)
      throw std::invalid_argument("ER FIR mixing horizon is invalid");
    if (analysisPathCount < imagePathCount)
      throw std::invalid_argument(
          "ER FIR audible path count exceeds its analysis path count");
    for (const auto &handoff : sourceHandoffs) {
      if (!std::isfinite(handoff.directPropagationSeconds) ||
          handoff.directPropagationSeconds < 0.0 ||
          !std::isfinite(handoff.finalStartSeconds) ||
          !std::isfinite(handoff.finalEndSeconds) ||
          handoff.finalStartSeconds < 0.0 ||
          handoff.finalEndSeconds < handoff.finalStartSeconds ||
          handoff.finalEndSeconds >
              mixingPrediction.generationHorizonSeconds + 1.0e-12)
        throw std::invalid_argument("ER FIR source handoff is invalid");
      if (!std::isfinite(handoff.broadbandHandoffPower) ||
          handoff.broadbandHandoffPower < 0.0 ||
          !std::isfinite(handoff.broadbandEarlyEnergy) ||
          handoff.broadbandEarlyEnergy < 0.0)
        throw std::invalid_argument("ER FIR handoff energy is invalid");
      for (std::size_t band = 0; band < EarlyReflectionBandCount; ++band)
        if (!std::isfinite(handoff.bandHandoffPower[band]) ||
            handoff.bandHandoffPower[band] < 0.0 ||
            !std::isfinite(handoff.bandEarlyEnergy[band]) ||
            handoff.bandEarlyEnergy[band] < 0.0)
          throw std::invalid_argument("ER FIR band energy is invalid");
      if (handoff.analysisPathCount < handoff.imagePathCount)
        throw std::invalid_argument(
            "ER source audible path count exceeds its analysis path count");
    }
    const std::size_t size = Size();
    if (size == 0)
      throw std::invalid_argument("ER FIR cannot be empty");
    for (const auto &kernel : kernels) {
      if (kernel.size() != size)
        throw std::invalid_argument(
            "all ER FIR kernels must have equal length");
      for (const double sample : kernel)
        if (!std::isfinite(sample))
          throw std::invalid_argument("ER FIR coefficients must be finite");
    }
  }
};


} // namespace tfdsp
