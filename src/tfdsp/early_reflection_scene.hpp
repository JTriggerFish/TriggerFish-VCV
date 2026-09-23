#pragma once

#include "early_reflection_types.hpp"

namespace tfdsp {
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
} // namespace early_reflection_detail

} // namespace tfdsp
