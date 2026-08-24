#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>

#include "tfdsp/early_reflections.hpp"

namespace {
using Clock = std::chrono::steady_clock;
using Convolver = tfdsp::EarlyReflectionConvolver<float, 128, 8>;

void Run(const std::size_t sourceCount) {
  tfdsp::EarlyReflectionConfig config;
  config.responseTimeSeconds = 0.150;
  const auto room = tfdsp::MakeEarlyReflectionRoom(0.45);
  const auto sources =
      tfdsp::MakeDefaultEarlyReflectionSources(room, sourceCount, 0.5);
  const auto materials = tfdsp::MakeEarlyReflectionMaterials(0.5);

  const auto buildStart = Clock::now();
  const auto response = tfdsp::BuildEarlyReflectionImpulseResponse(
      config, room, sources, materials);
  const double buildMilliseconds =
      std::chrono::duration<double, std::milli>(Clock::now() - buildStart)
          .count();

  Convolver convolver;
  convolver.Prepare(config.sampleRate,
                    tfdsp::MaximumEarlyReflectionImpulseSamples(config, room));
  const auto fftStart = Clock::now();
  convolver.SetImpulseResponse(response);
  const double fftMilliseconds =
      std::chrono::duration<double, std::milli>(Clock::now() - fftStart)
          .count();

  constexpr std::size_t warmupSamples = 48000;
  constexpr std::size_t measuredSamples = 10 * 48000;
  float checksum = 0.0f;
  for (std::size_t sample = 0; sample < warmupSamples; ++sample) {
    Convolver::InputFrame input{};
    for (std::size_t source = 0; source < sourceCount; ++source)
      input[source] = 0.1f * static_cast<float>(source + 1);
    const auto output = convolver.Process(input, sourceCount);
    checksum += output[0] + output[1];
  }

  const auto start = Clock::now();
  for (std::size_t sample = 0; sample < measuredSamples; ++sample) {
    Convolver::InputFrame input{};
    for (std::size_t source = 0; source < sourceCount; ++source) {
      const double phase = 0.017 * static_cast<double>((source + 1) * sample);
      input[source] = static_cast<float>(0.1 * std::sin(phase));
    }
    const auto output = convolver.Process(input, sourceCount);
    checksum += output[0] + output[1];
  }
  const double seconds =
      std::chrono::duration<double>(Clock::now() - start).count();
  const double realtimeFraction =
      seconds / (static_cast<double>(measuredSamples) / config.sampleRate);

  std::cout << std::fixed << std::setprecision(3) << "sources: " << sourceCount
            << '\n'
            << "audible image paths: " << response.imagePathCount << '\n'
            << "analysis image paths: " << response.analysisPathCount << '\n'
            << "FIR samples: " << response.Size() << '\n'
            << "FIR generation ms: " << buildMilliseconds << '\n'
            << "partition FFT preparation ms: " << fftMilliseconds << '\n'
            << "worker update total ms: " << buildMilliseconds + fftMilliseconds
            << '\n'
            << "predicted handoff ms: "
            << 1000.0 * response.sourceHandoffs.front().predictedStartSeconds
            << " to "
            << 1000.0 * response.sourceHandoffs.front().predictedEndSeconds
            << '\n'
            << "single-core real-time fraction: " << realtimeFraction << '\n'
            << "checksum: " << checksum << "\n\n";
}
} // namespace

int main() {
  Run(1);
  Run(4);
  Run(8);
}
