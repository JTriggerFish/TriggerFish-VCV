#pragma once

#include "early_reflection_paths.hpp"

namespace tfdsp {
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


} // namespace tfdsp
