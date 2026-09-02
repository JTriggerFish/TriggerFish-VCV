#include "crash_api.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

int main() {
  constexpr std::uint32_t FrameCount = 8192;
  const auto handle = tf_crash_create(48000.f);
  if (!handle || !tf_crash_trigger(handle, .8f, .8f, .65f, 17))
    return 1;
  std::vector<float> audio(FrameCount);
  if (!tf_crash_process(handle, audio.data(), FrameCount))
    return 1;
  tf_crash_destroy(handle);

  double peak = 0.0;
  double energy = 0.0;
  double absoluteSum = 0.0;
  double earlyEnergy = 0.0;
  for (std::size_t index = 0; index < audio.size(); ++index) {
    const double sample = audio[index];
    peak = std::max(peak, std::abs(sample));
    energy += sample * sample;
    absoluteSum += std::abs(sample);
    if (index < 1024)
      earlyEnergy += sample * sample;
  }
  std::cout.precision(17);
  std::cout << "{\"api\":1,\"frames\":" << FrameCount
            << ",\"peak\":" << peak << ",\"energy\":" << energy
            << ",\"absoluteSum\":" << absoluteSum
            << ",\"earlyEnergy\":" << earlyEnergy << "}\n";
}
