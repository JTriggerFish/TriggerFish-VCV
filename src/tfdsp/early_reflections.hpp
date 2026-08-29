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

inline double EarlyReflectionSmoothstep(const double value) noexcept {
  const double amount = std::clamp(value, 0.0, 1.0);
  return amount * amount * (3.0 - 2.0 * amount);
}

inline std::array<double, 4>
EarlyReflectionLagrange4(const double fraction) noexcept {
  const double mu = std::clamp(fraction, 0.0, 1.0);
  return {-mu * (mu - 1.0) * (mu - 2.0) / 6.0,
          (mu + 1.0) * (mu - 1.0) * (mu - 2.0) / 2.0,
          -(mu + 1.0) * mu * (mu - 2.0) / 2.0,
          (mu + 1.0) * mu * (mu - 1.0) / 6.0};
}

inline std::array<double, 4>
EarlyReflectionCausalLagrange4(const double delaySamples) noexcept {
  const double delay = std::clamp(delaySamples, 0.0, 1.0);
  return {-(delay - 1.0) * (delay - 2.0) * (delay - 3.0) / 6.0,
          delay * (delay - 2.0) * (delay - 3.0) / 2.0,
          -delay * (delay - 1.0) * (delay - 3.0) / 2.0,
          delay * (delay - 1.0) * (delay - 2.0) / 6.0};
}

inline std::size_t
EarlyReflectionImageCountThroughOrder(const std::size_t maximumOrder) noexcept {
  return maximumOrder == 0 ? 0
                           : 2 * maximumOrder * (maximumOrder + 1) *
                                     (2 * maximumOrder + 1) / 3 +
                                 2 * maximumOrder;
}

namespace detail {

// These constructors are shared with the realtime room-reverb control path.
// Their callers must provide the bounded inputs documented by the checked
// public wrappers below; keeping validation out of these functions makes that
// realtime path structurally non-throwing.
inline EarlyReflectionRoom MakeEarlyReflectionRoomFromValidatedInputs(
    const double space, const EarlyReflectionRoomFamily family) noexcept {
  EarlyReflectionVector minimum{};
  EarlyReflectionVector maximum{};
  switch (family) {
  case EarlyReflectionRoomFamily::Room:
    minimum = {2.8, 3.5, 2.4};
    maximum = {18.0, 25.0, 8.0};
    break;
  case EarlyReflectionRoomFamily::Chamber:
    minimum = {2.5, 2.8, 2.4};
    maximum = {12.0, 15.0, 6.0};
    break;
  case EarlyReflectionRoomFamily::Hall:
    minimum = {5.0, 8.0, 3.0};
    maximum = {24.0, 35.0, 14.0};
    break;
  }
  EarlyReflectionRoom room;
  const double amount = EarlyReflectionSmoothstep(space);
  for (std::size_t axis = 0; axis < 3; ++axis)
    room.dimensionsMetres[axis] =
        std::exp(std::log(minimum[axis]) +
                 amount * (std::log(maximum[axis]) - std::log(minimum[axis])));
  const EarlyReflectionVector listenerFraction{0.5, 0.682, 0.45};
  for (std::size_t axis = 0; axis < 3; ++axis)
    room.listenerPositionMetres[axis] =
        listenerFraction[axis] * room.dimensionsMetres[axis];
  return room;
}

inline EarlyReflectionSource
MakeEarlyReflectionSourceFromValidatedInputs(
    const EarlyReflectionRoom &room,
    const EarlyReflectionVector &normalizedPosition) noexcept {
  EarlyReflectionSource source;
  for (std::size_t axis = 0; axis < 3; ++axis) {
    constexpr double margin = 0.001;
    const double position =
        std::clamp(normalizedPosition[axis], margin, 1.0 - margin);
    source.positionMetres[axis] = position * room.dimensionsMetres[axis];
  }
  return source;
}

} // namespace detail

inline EarlyReflectionRoom MakeEarlyReflectionRoom(
    const double space,
    const EarlyReflectionRoomFamily family = EarlyReflectionRoomFamily::Room) {
  if (!std::isfinite(space) || space < 0.0 || space > 1.0)
    throw std::invalid_argument("Space must be in [0, 1]");
  auto room =
      detail::MakeEarlyReflectionRoomFromValidatedInputs(space, family);
  room.Validate();
  return room;
}

inline EarlyReflectionSource
MakeEarlyReflectionSource(const EarlyReflectionRoom &room,
                          const EarlyReflectionVector &normalizedPosition) {
  room.Validate();
  for (std::size_t axis = 0; axis < 3; ++axis) {
    if (!std::isfinite(normalizedPosition[axis]) ||
        normalizedPosition[axis] < 0.0 || normalizedPosition[axis] > 1.0)
      throw std::invalid_argument(
          "normalized source coordinates must lie in [0, 1]");
  }
  return detail::MakeEarlyReflectionSourceFromValidatedInputs(
      room, normalizedPosition);
}

inline std::vector<EarlyReflectionSource>
MakeDefaultEarlyReflectionSources(const EarlyReflectionRoom &room,
                                  const std::size_t sourceCount = 1,
                                  const double distance = 0.5) {
  if (sourceCount == 0 || sourceCount > EarlyReflectionMaximumSources)
    throw std::invalid_argument("source count must be between 1 and 8");
  if (!std::isfinite(distance) || distance < 0.0 || distance > 1.0)
    throw std::invalid_argument("Distance must be in [0, 1]");
  const double distanceAmount = EarlyReflectionSmoothstep(distance);
  const double y = 0.52 + distanceAmount * (0.12 - 0.52);
  const double z = 0.45 + distanceAmount * (0.42 - 0.45);
  std::vector<EarlyReflectionSource> sources;
  sources.reserve(sourceCount);
  for (std::size_t source = 0; source < sourceCount; ++source) {
    const double x = sourceCount == 1
                         ? 0.5
                         : 0.30 + 0.40 * static_cast<double>(source) /
                                      static_cast<double>(sourceCount - 1);
    sources.push_back(MakeEarlyReflectionSource(room, {x, y, z}));
  }
  return sources;
}

namespace detail {

inline EarlyReflectionMaterials MakeEarlyReflectionMaterialsFromValidatedInput(
    const double damping) noexcept {
  constexpr std::array<std::array<double, EarlyReflectionBandCount>,
                       EarlyReflectionSurfaceCount>
      bright{{{{0.97, 0.94, 0.88, 0.78}},
              {{0.96, 0.91, 0.82, 0.68}},
              {{0.97, 0.93, 0.85, 0.72}},
              {{0.96, 0.92, 0.83, 0.69}}}};
  constexpr std::array<std::array<double, EarlyReflectionBandCount>,
                       EarlyReflectionSurfaceCount>
      damped{{{{0.92, 0.82, 0.60, 0.32}},
              {{0.90, 0.75, 0.45, 0.15}},
              {{0.92, 0.78, 0.50, 0.20}},
              {{0.90, 0.74, 0.42, 0.16}}}};
  EarlyReflectionMaterials materials;
  const double amount = EarlyReflectionSmoothstep(damping);
  for (std::size_t surface = 0; surface < EarlyReflectionSurfaceCount;
       ++surface)
    for (std::size_t band = 0; band < EarlyReflectionBandCount; ++band)
      materials.reflectionAmplitudes[surface][band] =
          std::exp(std::log(bright[surface][band]) +
                   amount * (std::log(damped[surface][band]) -
                             std::log(bright[surface][band])));
  materials.airAbsorptionDbPerMetre = {0.0, 0.001, 0.006, 0.020};
  return materials;
}

} // namespace detail

inline EarlyReflectionMaterials
MakeEarlyReflectionMaterials(const double damping) {
  if (!std::isfinite(damping) || damping < 0.0 || damping > 1.0)
    throw std::invalid_argument("Damping must be in [0, 1]");
  return detail::MakeEarlyReflectionMaterialsFromValidatedInput(damping);
}

inline double
EarlyReflectionDirectDistance(const EarlyReflectionRoom &room,
                              const EarlyReflectionSource &source) noexcept {
  EarlyReflectionVector difference{};
  for (std::size_t axis = 0; axis < difference.size(); ++axis)
    difference[axis] =
        source.positionMetres[axis] - room.listenerPositionMetres[axis];
  return std::sqrt(difference[0] * difference[0] +
                   difference[1] * difference[1] +
                   difference[2] * difference[2]);
}

inline EarlyReflectionMixingPrediction
PredictEarlyReflectionMixing(const EarlyReflectionConfig &config,
                             const EarlyReflectionRoom &room) {
  config.Validate();
  room.Validate();

  const double length = room.dimensionsMetres[0];
  const double width = room.dimensionsMetres[1];
  const double height = room.dimensionsMetres[2];
  const double volume = length * width * height;
  const double surface =
      2.0 * (length * width + length * height + width * height);
  const double average = 0.001 * (20.0 * volume / surface + 12.0);
  const double conservative = 0.001 * (0.0117 * volume + 50.1);

  EarlyReflectionMixingPrediction prediction;
  prediction.roomVolumeCubicMetres = volume;
  prediction.roomSurfaceSquareMetres = surface;
  prediction.averageSeconds = std::clamp(average, config.minimumHandoffSeconds,
                                         config.responseTimeSeconds);
  prediction.conservativeSeconds =
      std::clamp(std::max(conservative, prediction.averageSeconds +
                                            config.handoffOverlapSeconds),
                 prediction.averageSeconds, config.responseTimeSeconds);
  prediction.generationHorizonSeconds =
      std::min(config.responseTimeSeconds,
               prediction.conservativeSeconds + config.analysisSafetySeconds);
  return prediction;
}

namespace early_reflection_detail {
inline int Parity(const int value) noexcept { return std::abs(value) % 2; }

inline double ImageCoordinate(const int index, const double dimension,
                              const double source) noexcept {
  const int parity = Parity(index);
  const int lattice = (index + parity) / 2;
  return 2.0 * lattice * dimension + (1.0 - 2.0 * parity) * source;
}

inline std::pair<int, int> LowerUpperCounts(const int index) noexcept {
  const int absolute = std::abs(index);
  const int major = (absolute + 1) / 2;
  const int minor = absolute / 2;
  if (index < 0)
    return {major, minor};
  if (index > 0)
    return {minor, major};
  return {0, 0};
}

inline double VectorLength(const EarlyReflectionVector &vector) noexcept {
  return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] +
                   vector[2] * vector[2]);
}

inline void ValidateSources(const EarlyReflectionRoom &room,
                            const std::vector<EarlyReflectionSource> &sources) {
  if (sources.empty() || sources.size() > EarlyReflectionMaximumSources)
    throw std::invalid_argument(
        "the ER scene must contain between 1 and 8 sources");
  for (const auto &source : sources)
    for (std::size_t axis = 0; axis < 3; ++axis)
      if (!std::isfinite(source.positionMetres[axis]) ||
          source.positionMetres[axis] <= 0.0 ||
          source.positionMetres[axis] >= room.dimensionsMetres[axis])
        throw std::invalid_argument("every source must lie inside the room");
}

class TwoPoleLowPass {
  double coefficient_{};
  double state1_{};
  double state2_{};

public:
  void Prepare(const double cutoffHz, const double sampleRate) noexcept {
    coefficient_ =
        1.0 - std::exp(-2.0 * EarlyReflectionPi * cutoffHz / sampleRate);
    state1_ = 0.0;
    state2_ = 0.0;
  }

  double Process(const double input) noexcept {
    state1_ += coefficient_ * (input - state1_);
    state2_ += coefficient_ * (state1_ - state2_);
    return state2_;
  }
};

class ComplementaryBandFilter {
  std::array<TwoPoleLowPass, EarlyReflectionBandCount - 1> lowPasses_{};
  std::size_t band_{};

public:
  void Prepare(const EarlyReflectionConfig &config,
               const std::size_t band) noexcept {
    band_ = band;
    for (std::size_t crossover = 0; crossover < lowPasses_.size(); ++crossover)
      lowPasses_[crossover].Prepare(config.crossoverHz[crossover],
                                    config.sampleRate);
  }

  double Process(const double input) noexcept {
    double residual = input;
    for (std::size_t crossover = 0; crossover < lowPasses_.size();
         ++crossover) {
      const double low = lowPasses_[crossover].Process(residual);
      if (band_ == crossover)
        return low;
      residual -= low;
    }
    return residual;
  }
};
} // namespace early_reflection_detail

inline std::vector<EarlyReflectionPath> EnumerateEarlyReflectionPaths(
    const EarlyReflectionConfig &config, const EarlyReflectionRoom &room,
    const std::vector<EarlyReflectionSource> &sources,
    const EarlyReflectionMaterials &materials) {
  config.Validate();
  room.Validate();
  materials.Validate();
  early_reflection_detail::ValidateSources(room, sources);
  const auto mixing = PredictEarlyReflectionMixing(config, room);
  std::vector<double> directDistances;
  directDistances.reserve(sources.size());
  double maximumDistance = 0.0;
  for (const auto &source : sources) {
    const double directDistance = EarlyReflectionDirectDistance(room, source);
    directDistances.push_back(directDistance);
    maximumDistance = std::max(
        maximumDistance,
        directDistance + mixing.generationHorizonSeconds * config.speedOfSound);
  }
  std::array<int, 3> limits{};
  for (std::size_t axis = 0; axis < 3; ++axis)
    limits[axis] = static_cast<int>(std::ceil(maximumDistance /
                                              room.dimensionsMetres[axis])) +
                   2;

  std::vector<EarlyReflectionPath> paths;
  for (std::size_t sourceIndex = 0; sourceIndex < sources.size(); ++sourceIndex)
    for (int x = -limits[0]; x <= limits[0]; ++x)
      for (int y = -limits[1]; y <= limits[1]; ++y)
        for (int z = -limits[2]; z <= limits[2]; ++z) {
          if (x == 0 && y == 0 && z == 0)
            continue;
          const std::array<int, 3> imageIndex{x, y, z};
          EarlyReflectionVector imagePosition{};
          for (std::size_t axis = 0; axis < 3; ++axis)
            imagePosition[axis] = early_reflection_detail::ImageCoordinate(
                imageIndex[axis], room.dimensionsMetres[axis],
                sources[sourceIndex].positionMetres[axis]);
          EarlyReflectionVector pathVector{};
          for (std::size_t axis = 0; axis < 3; ++axis)
            pathVector[axis] =
                imagePosition[axis] - room.listenerPositionMetres[axis];
          const double distance =
              early_reflection_detail::VectorLength(pathVector);
          const double excessDistance = distance - directDistances[sourceIndex];
          if (excessDistance < -1.0e-9 ||
              excessDistance >
                  mixing.generationHorizonSeconds * config.speedOfSound)
            continue;

          EarlyReflectionPath path;
          path.imageIndex = imageIndex;
          path.sourceIndex = sourceIndex;
          path.reflectionOrder =
              static_cast<std::size_t>(std::abs(x) + std::abs(y) + std::abs(z));
          path.imagePositionMetres = imagePosition;
          path.distanceMetres = distance;
          path.propagationSeconds = distance / config.speedOfSound;
          path.excessDelaySeconds =
              std::max(0.0, excessDistance / config.speedOfSound);
          for (std::size_t axis = 0; axis < 3; ++axis)
            path.arrivalDirection[axis] = pathVector[axis] / distance;

          const auto xCounts = early_reflection_detail::LowerUpperCounts(x);
          const auto yCounts = early_reflection_detail::LowerUpperCounts(y);
          const auto zCounts = early_reflection_detail::LowerUpperCounts(z);
          path.surfaceCounts = {
              static_cast<std::uint16_t>(zCounts.first),
              static_cast<std::uint16_t>(zCounts.second),
              static_cast<std::uint16_t>(xCounts.first + xCounts.second),
              static_cast<std::uint16_t>(yCounts.first + yCounts.second)};

          const double horizontalLength =
              std::hypot(pathVector[0], pathVector[1]);
          const double lateral = horizontalLength > 1.0e-12
                                     ? std::clamp(pathVector[0] /
                                                      horizontalLength,
                                                  -1.0, 1.0)
                                     : 0.0;
          const double pan = 0.5 * (lateral + 1.0);
          path.outputGains = {std::cos(0.5 * EarlyReflectionPi * pan),
                              std::sin(0.5 * EarlyReflectionPi * pan)};

          const double spreading = config.referenceDistanceMetres / distance;
          for (std::size_t band = 0; band < EarlyReflectionBandCount; ++band) {
            double materialGain = 1.0;
            for (std::size_t surface = 0; surface < EarlyReflectionSurfaceCount;
                 ++surface)
              materialGain *=
                  std::pow(materials.reflectionAmplitudes[surface][band],
                           path.surfaceCounts[surface]);
            const double airGain = std::pow(
                10.0,
                -distance * materials.airAbsorptionDbPerMetre[band] / 20.0);
            path.bandGains[band] = spreading * airGain * materialGain;
          }
          paths.push_back(path);
        }

  std::sort(
      paths.begin(), paths.end(),
      [](const EarlyReflectionPath &left, const EarlyReflectionPath &right) {
        if (left.propagationSeconds != right.propagationSeconds)
          return left.propagationSeconds < right.propagationSeconds;
        if (left.sourceIndex != right.sourceIndex)
          return left.sourceIndex < right.sourceIndex;
        return left.imageIndex < right.imageIndex;
      });
  return paths;
}

inline std::size_t
MaximumEarlyReflectionImpulseSamples(const EarlyReflectionConfig &config,
                                     const EarlyReflectionRoom &room) {
  config.Validate();
  room.Validate();
  return static_cast<std::size_t>(
             std::ceil(config.responseTimeSeconds * config.sampleRate)) +
         4;
}

struct EarlyReflectionDensityStatistics {
  double normalizedEchoDensity{};
  double excessKurtosis{};
  bool valid{};
};

inline EarlyReflectionDensityStatistics CalculateEarlyReflectionDensity(
    const std::vector<double> &left, const std::vector<double> &right,
    const std::size_t firstSample, const std::size_t sampleCount) noexcept {
  EarlyReflectionDensityStatistics statistics;
  if (sampleCount < 4 || firstSample >= left.size() ||
      firstSample >= right.size() || firstSample + sampleCount > left.size() ||
      firstSample + sampleCount > right.size())
    return statistics;

  const double count = static_cast<double>(2 * sampleCount);
  double mean = 0.0;
  for (std::size_t sample = firstSample; sample < firstSample + sampleCount;
       ++sample)
    mean += left[sample] + right[sample];
  mean /= count;

  double secondMoment = 0.0;
  double fourthMoment = 0.0;
  for (std::size_t sample = firstSample; sample < firstSample + sampleCount;
       ++sample)
    for (const double value : {left[sample], right[sample]}) {
      const double centred = value - mean;
      const double squared = centred * centred;
      secondMoment += squared;
      fourthMoment += squared * squared;
    }
  secondMoment /= count;
  fourthMoment /= count;
  if (!std::isfinite(secondMoment) || secondMoment <= 1.0e-30)
    return statistics;

  const double deviation = std::sqrt(secondMoment);
  std::size_t outside = 0;
  for (std::size_t sample = firstSample; sample < firstSample + sampleCount;
       ++sample) {
    outside +=
        static_cast<std::size_t>(std::abs(left[sample] - mean) > deviation);
    outside +=
        static_cast<std::size_t>(std::abs(right[sample] - mean) > deviation);
  }
  constexpr double GaussianOutsideOneSigma = 0.31731050786291415;
  statistics.normalizedEchoDensity =
      static_cast<double>(outside) / (count * GaussianOutsideOneSigma);
  statistics.excessKurtosis =
      fourthMoment / (secondMoment * secondMoment) - 3.0;
  statistics.valid = std::isfinite(statistics.normalizedEchoDensity) &&
                     std::isfinite(statistics.excessKurtosis);
  return statistics;
}

inline void RefineEarlyReflectionHandoffs(
    const EarlyReflectionConfig &config,
    EarlyReflectionImpulseResponse &response) noexcept {
  const std::size_t windowSamples = std::max<std::size_t>(
      4, static_cast<std::size_t>(
             std::llround(config.analysisWindowSeconds * config.sampleRate)));
  const std::size_t hopSamples = std::max<std::size_t>(
      1, static_cast<std::size_t>(
             std::llround(config.analysisHopSeconds * config.sampleRate)));
  const std::size_t stableWindows = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::ceil(config.analysisStableSeconds /
                                            config.analysisHopSeconds)));

  for (std::size_t source = 0; source < response.sourceCount; ++source) {
    auto &handoff = response.sourceHandoffs[source];
    const auto &left = response.kernels[response.KernelIndex(0, source)];
    const auto &right = response.kernels[response.KernelIndex(1, source)];
    const double lowerRelative = 0.75 * handoff.predictedStartSeconds;
    const double upperRelative =
        std::min(response.mixingPrediction.generationHorizonSeconds,
                 1.25 * handoff.predictedEndSeconds);
    const auto firstCentre = static_cast<std::size_t>(std::max(
        0.0, std::ceil(lowerRelative * config.sampleRate)));
    const auto lastCentre = static_cast<std::size_t>(std::max(
        0.0, std::floor(upperRelative * config.sampleRate)));

    std::size_t consecutive = 0;
    double firstStableRelative = 0.0;
    for (std::size_t centre = firstCentre;
         centre <= lastCentre && centre < response.Size();
         centre += hopSamples) {
      const std::size_t halfWindow = windowSamples / 2;
      if (centre < halfWindow ||
          centre - halfWindow + windowSamples > response.Size())
        continue;
      const auto statistics = CalculateEarlyReflectionDensity(
          left, right, centre - halfWindow, windowSamples);
      const bool diffuse = statistics.valid &&
                           statistics.normalizedEchoDensity >= 0.9 &&
                           statistics.normalizedEchoDensity <= 1.1 &&
                           std::abs(statistics.excessKurtosis) < 0.5;
      if (!diffuse) {
        consecutive = 0;
        continue;
      }
      if (consecutive == 0)
        firstStableRelative = static_cast<double>(centre) / config.sampleRate;
      ++consecutive;
      if (consecutive < stableWindows)
        continue;

      handoff.detectedFromResponse = true;
      handoff.detectedMixingSeconds = firstStableRelative;
      handoff.normalizedEchoDensity = statistics.normalizedEchoDensity;
      handoff.excessKurtosis = statistics.excessKurtosis;
      handoff.finalStartSeconds =
          std::clamp(firstStableRelative, lowerRelative, upperRelative);
      handoff.finalEndSeconds = std::min(
          response.mixingPrediction.generationHorizonSeconds,
          std::max(handoff.predictedEndSeconds,
                   handoff.finalStartSeconds + config.handoffOverlapSeconds));
      break;
    }
  }
}

inline EarlyReflectionImpulseResponse BuildEarlyReflectionImpulseResponse(
    const EarlyReflectionConfig &config, const EarlyReflectionRoom &room,
    const std::vector<EarlyReflectionSource> &sources,
    const EarlyReflectionMaterials &materials,
    const std::vector<EarlyReflectionPath> *cachedGeometry = nullptr) {
  const auto mixing = PredictEarlyReflectionMixing(config, room);
  materials.Validate();
  early_reflection_detail::ValidateSources(room, sources);
  auto paths = cachedGeometry != nullptr
                   ? *cachedGeometry
                   : EnumerateEarlyReflectionPaths(config, room, sources,
                                                   materials);
  if (paths.empty())
    throw std::logic_error("ER geometry produced no reflected paths");
  for (auto &path : paths) {
    if (path.sourceIndex >= sources.size() ||
        path.excessDelaySeconds > mixing.generationHorizonSeconds + 1.0e-12)
      throw std::invalid_argument(
          "cached ER geometry does not match the requested scene");
    const double spreading =
        config.referenceDistanceMetres / path.distanceMetres;
    for (std::size_t band = 0; band < EarlyReflectionBandCount; ++band) {
      double materialGain = 1.0;
      for (std::size_t surface = 0; surface < EarlyReflectionSurfaceCount;
           ++surface)
        materialGain *= std::pow(materials.reflectionAmplitudes[surface][band],
                                 path.surfaceCounts[surface]);
      const double airGain =
          std::pow(10.0, -path.distanceMetres *
                             materials.airAbsorptionDbPerMetre[band] / 20.0);
      path.bandGains[band] = spreading * airGain * materialGain;
    }
  }

  double maximumDesiredDelaySamples = 0.0;
  for (const auto &path : paths) {
    const double delay = path.excessDelaySeconds * config.sampleRate;
    maximumDesiredDelaySamples = std::max(maximumDesiredDelaySamples, delay);
  }
  const std::size_t tapTrainSize =
      static_cast<std::size_t>(std::ceil(maximumDesiredDelaySamples)) + 3;
  const std::size_t filterTailSamples = std::max<std::size_t>(
      1, static_cast<std::size_t>(
             std::ceil(config.materialFilterTailSeconds * config.sampleRate)));
  const std::size_t responseSize = tapTrainSize + filterTailSamples;
  const std::size_t kernelCount = EarlyReflectionOutputCount * sources.size();
  std::vector<std::vector<double>> bandTrains(
      kernelCount * EarlyReflectionBandCount,
      std::vector<double>(responseSize, 0.0));

  for (const auto &path : paths) {
    const double delaySamples = path.excessDelaySeconds * config.sampleRate;
    const auto integerDelay =
        static_cast<std::size_t>(std::floor(delaySamples));
    const bool needsCausalInterpolator = integerDelay == 0;
    const auto coefficients = needsCausalInterpolator
                                  ? EarlyReflectionCausalLagrange4(delaySamples)
                                  : EarlyReflectionLagrange4(
                                        delaySamples -
                                        static_cast<double>(integerDelay));
    constexpr std::array<int, 4> centredOffsets{-1, 0, 1, 2};
    for (std::size_t output = 0; output < EarlyReflectionOutputCount;
         ++output) {
      const std::size_t kernelIndex =
          output * sources.size() + path.sourceIndex;
      for (std::size_t band = 0; band < EarlyReflectionBandCount; ++band) {
        auto &train = bandTrains[kernelIndex * EarlyReflectionBandCount + band];
        const double gain = path.outputGains[output] * path.bandGains[band];
        for (std::size_t tap = 0; tap < coefficients.size(); ++tap) {
          const std::size_t index =
              needsCausalInterpolator
                  ? tap
                  : static_cast<std::size_t>(
                        static_cast<std::ptrdiff_t>(integerDelay) +
                        centredOffsets[tap]);
          train[index] += gain * coefficients[tap];
        }
      }
    }
  }

  EarlyReflectionImpulseResponse response;
  response.sampleRate = config.sampleRate;
  response.sourceCount = sources.size();
  response.analysisPathCount = paths.size();
  response.mixingPrediction = mixing;
  response.sourceHandoffs.resize(sources.size());
  for (std::size_t source = 0; source < sources.size(); ++source) {
    auto &handoff = response.sourceHandoffs[source];
    handoff.directPropagationSeconds =
        EarlyReflectionDirectDistance(room, sources[source]) /
        config.speedOfSound;
    handoff.predictedStartSeconds = mixing.averageSeconds;
    handoff.predictedEndSeconds = mixing.conservativeSeconds;
    handoff.finalStartSeconds = mixing.averageSeconds;
    handoff.finalEndSeconds = mixing.conservativeSeconds;
  }
  for (const auto &path : paths)
    ++response.sourceHandoffs[path.sourceIndex].analysisPathCount;
  response.kernels.assign(kernelCount, std::vector<double>(responseSize, 0.0));
  for (std::size_t kernel = 0; kernel < kernelCount; ++kernel) {
    std::array<early_reflection_detail::ComplementaryBandFilter,
               EarlyReflectionBandCount>
        filters{};
    for (std::size_t band = 0; band < EarlyReflectionBandCount; ++band)
      filters[band].Prepare(config, band);
    for (std::size_t sample = 0; sample < responseSize; ++sample)
      for (std::size_t band = 0; band < EarlyReflectionBandCount; ++band)
        response.kernels[kernel][sample] += filters[band].Process(
            bandTrains[kernel * EarlyReflectionBandCount + band][sample]);
  }
  RefineEarlyReflectionHandoffs(config, response);

  // Measure the untapered response at the detected handoff. This is done on
  // the worker thread, where the extra filter passes are harmless, and keeps
  // the audio thread free of scene-analysis work. The same complementary
  // filters used to build the audible kernels provide full-band coverage.
  const std::size_t handoffWindowSamples = std::max<std::size_t>(
      16, static_cast<std::size_t>(
              std::llround(config.analysisWindowSeconds * config.sampleRate)));
  for (std::size_t source = 0; source < sources.size(); ++source) {
    auto &handoff = response.sourceHandoffs[source];
    const auto handoffCentre = static_cast<std::size_t>(std::clamp(
        std::llround(handoff.finalStartSeconds * config.sampleRate),
        0LL, static_cast<long long>(responseSize - 1)));
    const std::size_t windowFirst =
        handoffCentre > handoffWindowSamples / 2
            ? handoffCentre - handoffWindowSamples / 2
            : 0;
    const std::size_t windowLast =
        std::min(responseSize, windowFirst + handoffWindowSamples);
    constexpr std::size_t earlyFirst = 0;

    for (std::size_t output = 0; output < EarlyReflectionOutputCount;
         ++output) {
      const std::size_t kernel = response.KernelIndex(output, source);
      for (std::size_t sample = earlyFirst; sample < handoffCentre; ++sample) {
        const double value = response.kernels[kernel][sample];
        handoff.broadbandEarlyEnergy += value * value;
      }
      for (std::size_t sample = windowFirst; sample < windowLast; ++sample) {
        const double value = response.kernels[kernel][sample];
        handoff.broadbandHandoffPower += value * value;
      }

      for (std::size_t band = 0; band < EarlyReflectionBandCount; ++band) {
        early_reflection_detail::ComplementaryBandFilter filter;
        filter.Prepare(config, band);
        for (std::size_t sample = 0; sample < responseSize; ++sample) {
          const double value = filter.Process(
              bandTrains[kernel * EarlyReflectionBandCount + band][sample]);
          if (sample >= earlyFirst && sample < handoffCentre)
            handoff.bandEarlyEnergy[band] += value * value;
          if (sample >= windowFirst && sample < windowLast)
            handoff.bandHandoffPower[band] += value * value;
        }
      }
    }
    const double divisor = static_cast<double>(
        EarlyReflectionOutputCount * std::max<std::size_t>(
                                         1, windowLast - windowFirst));
    handoff.broadbandHandoffPower /= divisor;
    for (auto &power : handoff.bandHandoffPower)
      power /= divisor;
  }

  response.imagePathCount = 0;
  for (const auto &path : paths) {
    auto &handoff = response.sourceHandoffs[path.sourceIndex];
    if (path.excessDelaySeconds <= handoff.finalEndSeconds) {
      ++handoff.imagePathCount;
      ++response.imagePathCount;
    }
  }

  std::size_t finalResponseSize = 1;
  for (std::size_t source = 0; source < sources.size(); ++source) {
    const auto &handoff = response.sourceHandoffs[source];
    const double first = handoff.finalStartSeconds;
    const double duration = handoff.finalEndSeconds - first;
    for (std::size_t output = 0; output < EarlyReflectionOutputCount;
         ++output) {
      auto &kernel = response.kernels[response.KernelIndex(output, source)];
      for (std::size_t sample = 0; sample < kernel.size(); ++sample) {
        const double relativeSeconds =
            static_cast<double>(sample) / config.sampleRate;
        double gain = 1.0;
        if (relativeSeconds >= handoff.finalEndSeconds)
          gain = 0.0;
        else if (relativeSeconds > first && duration > 0.0)
          gain = std::cos(
              0.5 * EarlyReflectionPi *
              EarlyReflectionSmoothstep((relativeSeconds - first) / duration));
        kernel[sample] *= gain;
      }
    }
    const double finalStoredSample =
        handoff.finalEndSeconds * config.sampleRate;
    finalResponseSize = std::max(finalResponseSize,
                                 static_cast<std::size_t>(std::max(
                                     1.0, std::ceil(finalStoredSample) + 1.0)));
  }
  finalResponseSize = std::min(finalResponseSize, response.Size());
  for (auto &kernel : response.kernels)
    kernel.resize(finalResponseSize);
  response.Validate();
  return response;
}

template <typename Sample = float, std::size_t PartitionSize = 128,
          std::size_t MaximumSources = EarlyReflectionMaximumSources>
class EarlyReflectionConvolver {
  static_assert(std::is_floating_point_v<Sample>,
                "EarlyReflectionConvolver requires floating-point samples");
  static_assert(PartitionSize >= 16 &&
                    (PartitionSize & (PartitionSize - 1)) == 0,
                "ER convolution partition size must be a power of two");
  static_assert(MaximumSources >= 1 &&
                    MaximumSources <= EarlyReflectionMaximumSources,
                "ER convolver source capacity must lie in [1, 8]");
  static constexpr std::size_t FftSize = 2 * PartitionSize;
  static constexpr std::size_t MaximumKernelCount =
      EarlyReflectionOutputCount * MaximumSources;
  static constexpr std::size_t BankCount = 4;
  static constexpr std::size_t InvalidBank = BankCount;

public:
  using InputFrame = std::array<Sample, MaximumSources>;
  using OutputFrame = std::array<Sample, EarlyReflectionOutputCount>;
  // The direct-form head is sample synchronous. The partitioned tail retains
  // an internal PartitionSize latency and stores coefficients shifted left by
  // that amount, so the complete convolver has zero implementation latency.
  static constexpr std::size_t LatencySamples = 0;
  static constexpr std::size_t HeadSize = PartitionSize;

private:
  using Spectrum = std::array<std::complex<Sample>, FftSize>;

  struct KernelBank {
    std::size_t sourceCount{};
    std::size_t partitionCount{};
    std::size_t transitionSamples{};
    std::size_t sceneSequence{};
    std::array<std::array<Sample, HeadSize>, MaximumKernelCount> head{};
    std::array<std::size_t, MaximumKernelCount> headTapCount{};
    std::array<std::vector<Spectrum>, MaximumKernelCount> spectra{};
  };

  enum class BankState : std::uint8_t { Free, Preparing, Ready, Active };

  double sampleRate_{};
  std::size_t maximumImpulseSamples_{};
  std::size_t maximumPartitions_{};
  std::size_t blockPosition_{};
  std::size_t blockSourceCount_{};
  std::size_t historyWriteIndex_{};
  std::size_t currentBankIndex_{InvalidBank};
  std::size_t previousBankIndex_{InvalidBank};
  std::size_t transitionSamples_{};
  std::size_t transitionRemaining_{};
  std::size_t activeSceneSequence_{};
  std::size_t renderedSceneSequence_{};
  bool prepared_{};
  std::array<std::complex<Sample>, FftSize / 2> roots_{};
  std::array<std::size_t, FftSize> bitReversed_{};
  std::array<KernelBank, BankCount> banks_{};
  std::array<std::atomic<BankState>, BankCount> bankStates_{};
  std::atomic<int> readyBankIndex_{-1};
  std::array<std::vector<Spectrum>, MaximumSources> inputHistory_{};
  std::array<std::array<Sample, HeadSize>, MaximumSources> headHistory_{};
  std::size_t headWriteIndex_{};
  std::array<std::size_t, MaximumSources> sourceFlushBlocks_{};
  std::array<std::array<Sample, PartitionSize>, MaximumSources>
      previousInputBlock_{};
  std::array<std::array<Sample, PartitionSize>, MaximumSources> inputBlock_{};
  std::array<std::array<Sample, PartitionSize>, EarlyReflectionOutputCount>
      currentOutputBlock_{};
  std::array<std::array<Sample, PartitionSize>, EarlyReflectionOutputCount>
      previousOutputBlock_{};

  void PrepareFft() noexcept {
    for (std::size_t index = 0; index < roots_.size(); ++index) {
      const double angle = -2.0 * EarlyReflectionPi *
                           static_cast<double>(index) /
                           static_cast<double>(FftSize);
      roots_[index] = {static_cast<Sample>(std::cos(angle)),
                       static_cast<Sample>(std::sin(angle))};
    }
    for (std::size_t index = 0; index < FftSize; ++index) {
      std::size_t value = index;
      std::size_t reversed = 0;
      for (std::size_t bit = 1; bit < FftSize; bit <<= 1) {
        reversed = (reversed << 1) | (value & 1);
        value >>= 1;
      }
      bitReversed_[index] = reversed;
    }
  }

  void Transform(Spectrum &values, const bool inverse) const noexcept {
    for (std::size_t index = 0; index < FftSize; ++index)
      if (index < bitReversed_[index])
        std::swap(values[index], values[bitReversed_[index]]);
    for (std::size_t length = 2; length <= FftSize; length <<= 1) {
      const std::size_t half = length / 2;
      const std::size_t rootStep = FftSize / length;
      for (std::size_t start = 0; start < FftSize; start += length)
        for (std::size_t offset = 0; offset < half; ++offset) {
          const auto root = inverse ? std::conj(roots_[offset * rootStep])
                                    : roots_[offset * rootStep];
          const auto even = values[start + offset];
          const auto odd = values[start + offset + half] * root;
          values[start + offset] = even + odd;
          values[start + offset + half] = even - odd;
        }
    }
    if (inverse)
      for (auto &value : values)
        value /= static_cast<Sample>(FftSize);
  }

  void LoadBank(const EarlyReflectionImpulseResponse &response,
                const std::size_t bankIndex) {
    response.Validate();
    if (!prepared_)
      throw std::logic_error("prepare the ER convolver before loading an FIR");
    if (std::abs(response.sampleRate - sampleRate_) > 1.0e-9)
      throw std::invalid_argument(
          "ER FIR sample rate does not match the convolver");
    if (response.sourceCount > MaximumSources)
      throw std::invalid_argument(
          "ER FIR exceeds the convolver source capacity");
    if (response.Size() > maximumImpulseSamples_)
      throw std::invalid_argument(
          "ER FIR exceeds the prepared convolution capacity");
    auto &bank = banks_[bankIndex];
    bank.sourceCount = response.sourceCount;
    bank.partitionCount =
        response.Size() > HeadSize
            ? (response.Size() - HeadSize + PartitionSize - 1) / PartitionSize
            : 0;
    for (auto &head : bank.head)
      head.fill(Sample{});
    bank.headTapCount.fill(0);
    for (std::size_t output = 0; output < EarlyReflectionOutputCount; ++output)
      for (std::size_t source = 0; source < response.sourceCount; ++source) {
        const std::size_t kernelIndex = output * MaximumSources + source;
        const auto &impulse =
            response.kernels[response.KernelIndex(output, source)];
        for (std::size_t sample = 0;
             sample < std::min(HeadSize, impulse.size()); ++sample) {
          bank.head[kernelIndex][sample] =
              static_cast<Sample>(impulse[sample]);
          if (std::abs(impulse[sample]) > 1.0e-20)
            bank.headTapCount[kernelIndex] = sample + 1;
        }
        for (std::size_t partition = 0; partition < bank.partitionCount;
             ++partition) {
          auto &spectrum = bank.spectra[kernelIndex][partition];
          spectrum.fill({});
          for (std::size_t sample = 0; sample < PartitionSize; ++sample) {
            const std::size_t impulseIndex =
                HeadSize + partition * PartitionSize + sample;
            if (impulseIndex < impulse.size())
              spectrum[sample] = static_cast<Sample>(impulse[impulseIndex]);
          }
          Transform(spectrum, false);
        }
      }
  }

  OutputFrame RenderHead(const KernelBank &bank) const noexcept {
    OutputFrame output{};
    for (std::size_t channel = 0; channel < EarlyReflectionOutputCount;
         ++channel)
      for (std::size_t source = 0; source < bank.sourceCount; ++source) {
        const auto &kernel = bank.head[channel * MaximumSources + source];
        const std::size_t tapCount =
            bank.headTapCount[channel * MaximumSources + source];
        for (std::size_t tap = 0; tap < tapCount; ++tap) {
          const std::size_t historyIndex =
              (headWriteIndex_ + HeadSize - tap) % HeadSize;
          output[channel] += kernel[tap] * headHistory_[source][historyIndex];
        }
      }
    return output;
  }

  void RenderBank(const KernelBank &bank,
                  std::array<std::array<Sample, PartitionSize>,
                             EarlyReflectionOutputCount> &output) noexcept {
    for (std::size_t channel = 0; channel < EarlyReflectionOutputCount;
         ++channel) {
      Spectrum accumulated{};
      for (std::size_t source = 0; source < bank.sourceCount; ++source) {
        if (sourceFlushBlocks_[source] == 0)
          continue;
        const std::size_t kernel = channel * MaximumSources + source;
        for (std::size_t partition = 0; partition < bank.partitionCount;
             ++partition) {
          const std::size_t historyIndex =
              (historyWriteIndex_ + maximumPartitions_ - partition) %
              maximumPartitions_;
          const auto &input = inputHistory_[source][historyIndex];
          const auto &filter = bank.spectra[kernel][partition];
          for (std::size_t bin = 0; bin < FftSize; ++bin)
            accumulated[bin] += input[bin] * filter[bin];
        }
      }
      Transform(accumulated, true);
      for (std::size_t sample = 0; sample < PartitionSize; ++sample)
        output[channel][sample] = accumulated[PartitionSize + sample].real();
    }
  }

  void ReleasePreviousBank() noexcept {
    if (previousBankIndex_ == InvalidBank)
      return;
    bankStates_[previousBankIndex_].store(BankState::Free,
                                          std::memory_order_release);
    previousBankIndex_ = InvalidBank;
  }

  void AdoptLatestReadyBank() noexcept {
    if (previousBankIndex_ != InvalidBank)
      return;
    const int candidate =
        readyBankIndex_.exchange(-1, std::memory_order_acq_rel);
    if (candidate < 0)
      return;
    const auto candidateIndex = static_cast<std::size_t>(candidate);
    BankState expected = BankState::Ready;
    if (!bankStates_[candidateIndex].compare_exchange_strong(
            expected, BankState::Active, std::memory_order_acq_rel,
            std::memory_order_acquire))
      return;

    if (currentBankIndex_ == InvalidBank) {
      currentBankIndex_ = candidateIndex;
      activeSceneSequence_ = banks_[candidateIndex].sceneSequence;
      return;
    }
    previousBankIndex_ = currentBankIndex_;
    currentBankIndex_ = candidateIndex;
    activeSceneSequence_ = banks_[candidateIndex].sceneSequence;
    transitionSamples_ = banks_[candidateIndex].transitionSamples;
    transitionRemaining_ = transitionSamples_;
  }

  void ProcessBlock() noexcept {
    for (std::size_t source = 0; source < MaximumSources; ++source) {
      const bool hasCurrentInput = source < blockSourceCount_;
      if (hasCurrentInput)
        sourceFlushBlocks_[source] = maximumPartitions_ + 1;
      if (!hasCurrentInput && sourceFlushBlocks_[source] == 0) {
        previousInputBlock_[source].fill(Sample{});
        continue;
      }
      Spectrum spectrum{};
      for (std::size_t sample = 0; sample < PartitionSize; ++sample) {
        spectrum[sample] = previousInputBlock_[source][sample];
        spectrum[PartitionSize + sample] = inputBlock_[source][sample];
      }
      Transform(spectrum, false);
      inputHistory_[source][historyWriteIndex_] = spectrum;
      previousInputBlock_[source] = inputBlock_[source];
      if (!hasCurrentInput)
        --sourceFlushBlocks_[source];
    }
    blockSourceCount_ = 0;
    if (transitionRemaining_ == 0)
      ReleasePreviousBank();
    AdoptLatestReadyBank();
    if (currentBankIndex_ != InvalidBank)
      RenderBank(banks_[currentBankIndex_], currentOutputBlock_);
    else
      for (auto &output : currentOutputBlock_)
        output.fill(Sample{});
    if (previousBankIndex_ != InvalidBank && transitionRemaining_ > 0)
      RenderBank(banks_[previousBankIndex_], previousOutputBlock_);

    ++historyWriteIndex_;
    if (historyWriteIndex_ == maximumPartitions_)
      historyWriteIndex_ = 0;
  }

public:
  void Prepare(const double sampleRate,
               const std::size_t maximumImpulseSamples) {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument(
          "convolver sample rate must be positive and finite");
    if (maximumImpulseSamples == 0)
      throw std::invalid_argument("convolver FIR capacity must be non-zero");
    sampleRate_ = sampleRate;
    maximumImpulseSamples_ = maximumImpulseSamples;
    maximumPartitions_ = std::max<std::size_t>(
        1, maximumImpulseSamples > HeadSize
               ? (maximumImpulseSamples - HeadSize + PartitionSize - 1) /
                     PartitionSize
               : 1);
    PrepareFft();
    for (auto &bank : banks_)
      for (auto &kernel : bank.spectra)
        kernel.assign(maximumPartitions_, Spectrum{});
    for (auto &history : inputHistory_)
      history.assign(maximumPartitions_, Spectrum{});
    for (auto &state : bankStates_)
      state.store(BankState::Free, std::memory_order_relaxed);
    readyBankIndex_.store(-1, std::memory_order_relaxed);
    currentBankIndex_ = InvalidBank;
    previousBankIndex_ = InvalidBank;
    activeSceneSequence_ = 0;
    renderedSceneSequence_ = 0;
    prepared_ = true;
    Reset();
  }

  void Reset() noexcept {
    for (auto &history : inputHistory_)
      for (auto &spectrum : history)
        spectrum.fill({});
    for (auto &block : previousInputBlock_)
      block.fill(Sample{});
    for (auto &history : headHistory_)
      history.fill(Sample{});
    for (auto &block : inputBlock_)
      block.fill(Sample{});
    for (auto &block : currentOutputBlock_)
      block.fill(Sample{});
    for (auto &block : previousOutputBlock_)
      block.fill(Sample{});
    blockPosition_ = 0;
    blockSourceCount_ = 0;
    historyWriteIndex_ = 0;
    headWriteIndex_ = 0;
    sourceFlushBlocks_.fill(0);
    transitionSamples_ = 0;
    transitionRemaining_ = 0;
    renderedSceneSequence_ = activeSceneSequence_;
  }

  void SetImpulseResponse(const EarlyReflectionImpulseResponse &response,
                          const std::size_t sceneSequence = 0) {
    if (!prepared_)
      throw std::logic_error("prepare the ER convolver before loading an FIR");
    readyBankIndex_.store(-1, std::memory_order_relaxed);
    for (auto &state : bankStates_)
      state.store(BankState::Free, std::memory_order_relaxed);
    currentBankIndex_ = 0;
    previousBankIndex_ = InvalidBank;
    bankStates_[currentBankIndex_].store(BankState::Preparing,
                                         std::memory_order_relaxed);
    LoadBank(response, currentBankIndex_);
    banks_[currentBankIndex_].sceneSequence = sceneSequence;
    activeSceneSequence_ = sceneSequence;
    bankStates_[currentBankIndex_].store(BankState::Active,
                                         std::memory_order_release);
    Reset();
  }

  bool
  PrepareAndQueueImpulseResponse(const EarlyReflectionImpulseResponse &response,
                                 const double transitionSeconds = 0.100,
                                 const std::size_t sceneSequence = 0) {
    if (!std::isfinite(transitionSeconds) || transitionSeconds < 0.0)
      throw std::invalid_argument(
          "ER transition time must be finite and non-negative");
    if (!prepared_)
      throw std::logic_error("prepare the ER convolver before loading an FIR");

    std::size_t bankIndex = InvalidBank;
    for (std::size_t candidate = 0; candidate < BankCount; ++candidate) {
      BankState expected = BankState::Free;
      if (bankStates_[candidate].compare_exchange_strong(
              expected, BankState::Preparing, std::memory_order_acq_rel,
              std::memory_order_acquire)) {
        bankIndex = candidate;
        break;
      }
    }
    if (bankIndex == InvalidBank)
      return false;

    try {
      LoadBank(response, bankIndex);
      banks_[bankIndex].sceneSequence = sceneSequence;
      banks_[bankIndex].transitionSamples =
          std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(
                                       transitionSeconds * sampleRate_)));
    } catch (...) {
      bankStates_[bankIndex].store(BankState::Free, std::memory_order_release);
      throw;
    }
    bankStates_[bankIndex].store(BankState::Ready, std::memory_order_release);
    const int replaced = readyBankIndex_.exchange(static_cast<int>(bankIndex),
                                                  std::memory_order_acq_rel);
    if (replaced >= 0) {
      BankState expected = BankState::Ready;
      bankStates_[static_cast<std::size_t>(replaced)].compare_exchange_strong(
          expected, BankState::Free, std::memory_order_acq_rel,
          std::memory_order_acquire);
    }
    return true;
  }

  // Scene identity used for the most recently returned audio sample.
  std::size_t RenderedSceneSequence() const noexcept {
    return renderedSceneSequence_;
  }

  void
  TransitionToImpulseResponse(const EarlyReflectionImpulseResponse &response,
                              const double transitionSeconds = 0.100,
                              const std::size_t sceneSequence = 0) {
    if (currentBankIndex_ == InvalidBank) {
      SetImpulseResponse(response, sceneSequence);
      return;
    }
    if (!PrepareAndQueueImpulseResponse(response, transitionSeconds,
                                        sceneSequence))
      throw std::runtime_error("no free ER convolution bank is available");
  }

  OutputFrame Process(const InputFrame &input,
                      const std::size_t sourceCount) noexcept {
    if (!prepared_)
      return {};
    renderedSceneSequence_ = activeSceneSequence_;
    const std::size_t activeSources = std::min(sourceCount, MaximumSources);
    for (std::size_t source = 0; source < MaximumSources; ++source) {
      const Sample sample =
          source < activeSources && std::isfinite(input[source])
              ? input[source]
              : Sample{};
      inputBlock_[source][blockPosition_] = sample;
      headHistory_[source][headWriteIndex_] = sample;
    }
    blockSourceCount_ = std::max(blockSourceCount_, activeSources);

    const OutputFrame currentHead =
        currentBankIndex_ != InvalidBank ? RenderHead(banks_[currentBankIndex_])
                                         : OutputFrame{};
    const OutputFrame previousHead =
        previousBankIndex_ != InvalidBank
            ? RenderHead(banks_[previousBankIndex_])
            : OutputFrame{};
    OutputFrame output{};
    for (std::size_t channel = 0; channel < EarlyReflectionOutputCount;
         ++channel) {
      output[channel] =
          currentHead[channel] + currentOutputBlock_[channel][blockPosition_];
      if (transitionRemaining_ > 0) {
        const std::size_t completed = transitionSamples_ - transitionRemaining_;
        const double linear =
            transitionSamples_ <= 1
                ? 1.0
                : static_cast<double>(completed) /
                      static_cast<double>(transitionSamples_ - 1);
        const Sample amount =
            static_cast<Sample>(EarlyReflectionSmoothstep(linear));
        const Sample previous = previousHead[channel] +
                                previousOutputBlock_[channel][blockPosition_];
        output[channel] = previous +
            amount * (output[channel] -
                      previous);
      }
    }

    if (transitionRemaining_ > 0)
      --transitionRemaining_;
    if (++headWriteIndex_ == HeadSize)
      headWriteIndex_ = 0;
    ++blockPosition_;
    if (blockPosition_ == PartitionSize) {
      blockPosition_ = 0;
      ProcessBlock();
    }
    return output;
  }

  bool IsPrepared() const noexcept { return prepared_; }
  bool IsTransitioning() const noexcept {
    return transitionRemaining_ > 0 ||
           readyBankIndex_.load(std::memory_order_acquire) >= 0;
  }
  std::size_t MaximumImpulseSamples() const noexcept {
    return maximumImpulseSamples_;
  }
};
} // namespace tfdsp
