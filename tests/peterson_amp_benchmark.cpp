#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>

#include "models/PetersonPowerAmplifier.hpp"
#include "models/PetersonPreAmplifier.hpp"
#include "models/PetersonTonePreAmplifier.hpp"

namespace {
using Clock = std::chrono::steady_clock;
constexpr double Rate = 96'000.0;
constexpr std::size_t Warmup = 24'000;
constexpr std::size_t Samples = 10 * 96'000;
constexpr double TwoPi = 6.2831853071795864769;

template <typename Circuit, typename Step>
void Run(const char *name, Circuit &circuit, Step &&step) {
  double checksum = 0.0;
  for (std::size_t sample = 0; sample < Warmup; ++sample)
    checksum += step(sample);
  const auto start = Clock::now();
  for (std::size_t sample = 0; sample < Samples; ++sample)
    checksum += step(Warmup + sample);
  const double elapsed =
      std::chrono::duration<double>(Clock::now() - start).count();
  const double audioSeconds = static_cast<double>(Samples) / Rate;
  std::cout << std::left << std::setw(22) << name << std::right << std::fixed
            << std::setprecision(4) << elapsed / audioSeconds
            << " core, avg iterations " << circuit.AverageIterations()
            << ", checksum " << checksum << '\n';
}
} // namespace

int main() {
  tfdsp::PetersonPreAmplifier preamp;
  preamp.SetSampleRate(Rate);
  Run("input preamp", preamp, [&](const std::size_t sample) {
    return preamp.Step(0.05 * std::sin(TwoPi * 250.0 * sample / Rate)).voltage;
  });

  tfdsp::PetersonTonePreAmplifier tone;
  tone.SetSampleRate(Rate);
  Run("tone preamp", tone, [&](const std::size_t sample) {
    return tone.Step(0.30 * std::sin(TwoPi * 250.0 * sample / Rate), 0.5, 0.5)
        .voltage;
  });

  tfdsp::PetersonPowerAmplifier power;
  power.SetSampleRate(Rate);
  Run("power channel", power, [&](const std::size_t sample) {
    return power.Step(0.05 * std::sin(TwoPi * 250.0 * sample / Rate), 35.0)
        .voltage;
  });
}
