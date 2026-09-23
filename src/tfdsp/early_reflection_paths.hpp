#pragma once

#include "early_reflection_filters.hpp"

namespace tfdsp {
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


} // namespace tfdsp
