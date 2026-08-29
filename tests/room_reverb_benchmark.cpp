#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>

#include "tfdsp/room_reverb.hpp"

namespace {
using Clock = std::chrono::steady_clock;

void Run(const float shimmer) {
  constexpr double SampleRate = 48'000.0;
  constexpr std::size_t SourceCount = 3;
  constexpr std::size_t WarmupSamples = 48'000;
  constexpr std::size_t MeasuredSamples = 10 * 48'000;

  tfdsp::RoomReverb reverb;
  reverb.SetSampleRate(SampleRate);
  tfdsp::RoomReverbControls controls;
  controls.space = 0.78f;
  controls.aspect = 0.55f;
  controls.preDelay = 0.112f;
  controls.decay = 0.84f;
  controls.damping = 0.26f;
  controls.diffusion = 1.f;
  controls.modulation = 1.f;
  controls.shimmer = shimmer;
  controls.width = 0.75f;
  controls.balance = 0.75f;

  tfdsp::RoomReverb::SourcePositions positions{};
  for (std::size_t source = 0; source < SourceCount; ++source) {
    positions[source] = tfdsp::reverb_defaults::Source;
    positions[source][0] = tfdsp::reverb_defaults::ProgressiveSourceX(
        source, SourceCount);
  }

  float checksum = 0.f;
  const auto render = [&](const std::size_t first, const std::size_t count) {
    for (std::size_t offset = 0; offset < count; ++offset) {
      const std::size_t sample = first + offset;
      tfdsp::RoomReverb::InputFrame input{};
      for (std::size_t source = 0; source < SourceCount; ++source)
        input[source] = static_cast<float>(
            0.1 * std::sin(0.017 * static_cast<double>((source + 1) * sample)));
      const auto output = reverb.Process(input, positions, SourceCount, controls);
      checksum += output.direct[0] + output.direct[1] + output.wet[0] +
                  output.wet[1];
    }
  };

  render(0, WarmupSamples);
  const auto start = Clock::now();
  render(WarmupSamples, MeasuredSamples);
  const double elapsed =
      std::chrono::duration<double>(Clock::now() - start).count();
  const double audioSeconds = static_cast<double>(MeasuredSamples) / SampleRate;
  std::cout << std::fixed << std::setprecision(4)
            << "shimmer: " << shimmer << '\n'
            << "single-core real-time fraction: " << elapsed / audioSeconds
            << '\n'
            << "nanoseconds per sample: "
            << 1.e9 * elapsed / static_cast<double>(MeasuredSamples) << '\n'
            << "checksum: " << checksum << "\n\n";
}
} // namespace

int main() {
  Run(0.f);
  Run(0.45f);
}
