#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>

#include "tfdsp/sampleRate.hpp"
#include "tfdsp/wavefolder.hpp"

namespace {
using Clock = std::chrono::steady_clock;
using Oscillator =
    tfdsp::WavefoldOscillator<tfdsp::X4Resampler_Order7>;

void Run(const bool renderOscillator) {
  constexpr double SampleRate = 48'000.0;
  constexpr std::size_t WarmupSamples = 48'000;
  constexpr std::size_t MeasuredSamples = 10 * 48'000;
  Oscillator oscillator(tfdsp::CreateX4Resampler_Cheby7);
  oscillator.SetSampleRate(SampleRate);
  // This matches the smoke patch's centre switch position.
  oscillator.SetCharacter(tfdsp::WavefolderCharacter::Hinge);

  double checksum = 0.0;
  const auto render = [&](const std::size_t count) {
    for (std::size_t sample = 0; sample < count; ++sample) {
      const auto output = oscillator.StepWithInput(
          164.81, 0.28, 0.50, 0.08, 0.0, false, renderOscillator, true);
      checksum += output.oscillator + output.folded;
    }
  };
  render(WarmupSamples);
  const auto start = Clock::now();
  render(MeasuredSamples);
  const double elapsed =
      std::chrono::duration<double>(Clock::now() - start).count();
  const double audioSeconds = static_cast<double>(MeasuredSamples) / SampleRate;
  std::cout << std::fixed << std::setprecision(4)
            << "raw oscillator output: "
            << (renderOscillator ? "on" : "off") << '\n'
            << "single-core real-time fraction: " << elapsed / audioSeconds
            << '\n'
            << "nanoseconds per sample: "
            << 1.e9 * elapsed / static_cast<double>(MeasuredSamples) << '\n'
            << "checksum: " << checksum << "\n\n";
}
} // namespace

int main() {
  Run(true);
  Run(false);
}
