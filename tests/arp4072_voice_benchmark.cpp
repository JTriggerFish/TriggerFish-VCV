#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>

#include "models/Arp4019Vca.hpp"
#include "models/Arp4072Filter.hpp"

namespace {
using Clock = std::chrono::steady_clock;
using Filter = tfdsp::Arp4072Filter<tfdsp::X4Resampler_Order7>;
using Vca = tfdsp::Arp4019Vca<tfdsp::X4Resampler_Order7>;

void Run(const double resonance, const double driveGain) {
  constexpr double SampleRate = 48'000.0;
  constexpr std::size_t WarmupSamples = 48'000;
  constexpr std::size_t MeasuredSamples = 10 * 48'000;
  Filter filter(tfdsp::CreateX4Resampler_Cheby7);
  Vca vca(tfdsp::CreateX4Resampler_Cheby7);
  filter.SetSampleRate(SampleRate);
  vca.SetSampleRate(SampleRate);

  double checksum = 0.0;
  std::size_t solverIterations = 0;
  const auto render = [&](const std::size_t first, const std::size_t count) {
    for (std::size_t offset = 0; offset < count; ++offset) {
      const std::size_t sample = first + offset;
      const double audio = 5.0 * std::sin(
          2.0 * 3.14159265358979323846 * 110.0 * sample / SampleRate);
      const auto output = filter.StepWithPostProcessorLogCutoff(
          audio, std::log2(1'200.0), resonance, driveGain, 0.0, 10.0,
          [&](const double filtered, const double linearCv,
              const double exponentialCv) {
            return vca.ProcessOversampled(filtered, linearCv, exponentialCv);
          });
      checksum += output.lowPass + output.postProcessed;
      solverIterations += static_cast<std::size_t>(filter.LastIterations());
    }
  };

  render(0, WarmupSamples);
  const auto start = Clock::now();
  render(WarmupSamples, MeasuredSamples);
  const double elapsed =
      std::chrono::duration<double>(Clock::now() - start).count();
  const double audioSeconds = static_cast<double>(MeasuredSamples) / SampleRate;
  std::cout << std::fixed << std::setprecision(4)
            << "resonance: " << resonance << ", drive: " << driveGain << '\n'
            << "single-core real-time fraction: " << elapsed / audioSeconds
            << '\n'
            << "nanoseconds per sample: "
            << 1.e9 * elapsed / static_cast<double>(MeasuredSamples) << '\n'
            << "last-internal-step solver iterations/sample: "
            << static_cast<double>(solverIterations) /
                   static_cast<double>(WarmupSamples + MeasuredSamples)
            << '\n'
            << "checksum: " << checksum << "\n\n";
}
} // namespace

int main() {
  Run(0.f, 1.5);
  Run(0.5, 3.3);
}
