#include "tfdsp/percussion/contact_exciter.hpp"
#include "tfdsp/percussion/coupled_resonator_network.hpp"
#include "tfdsp/percussion/dispersion_loop.hpp"
#include "tfdsp/percussion/frequency_shifter.hpp"

#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string_view>

namespace {

constexpr std::size_t SampleCount = 2'000'000;
volatile float sink{};

template <typename Process>
void Measure(const std::string_view name, Process process) {
  const auto start = std::chrono::steady_clock::now();
  float sum = 0.f;
  for (std::size_t sample = 0; sample < SampleCount; ++sample)
    sum += process(sample);
  const auto stop = std::chrono::steady_clock::now();
  sink = sum;
  const double seconds = std::chrono::duration<double>(stop - start).count();
  std::cout << name << ": " << 1.e9 * seconds / SampleCount
            << " ns/sample\n";
}

void BenchmarkFrequencyShifter() {
  tfdsp::percussion::FrequencyShifter shifter;
  shifter.Prepare(48000.f);
  shifter.SetShiftHz(317.f);
  Measure("frequency shifter", [&](const std::size_t sample) {
    return shifter.Process(std::sin(.071f * static_cast<float>(sample)));
  });
}

void BenchmarkContactExciter() {
  tfdsp::percussion::ContactExciter exciter;
  exciter.Prepare(48000.f);
  tfdsp::percussion::ContactExciterParameters parameters;
  Measure("active contact exciter", [&](const std::size_t sample) {
    if (sample % 128 == 0)
      exciter.Trigger(parameters);
    const auto output = exciter.Process();
    return output.directRadiation + output.bodyDrive;
  });
}

void BenchmarkResonators() {
  using Network = tfdsp::percussion::CoupledResonatorNetwork<12>;
  Network::Parameters parameters{};
  for (std::size_t line = 0; line < parameters.size(); ++line) {
    parameters[line].delaySamples = 18.f + 3.7f * static_cast<float>(line);
    parameters[line].inputGain = .08f;
    parameters[line].outputGain = .08f;
    parameters[line].decay = {2.f, 1.f, .4f};
  }
  Network network;
  network.Prepare(48000.f, 128.f, parameters, 500.f, 5500.f);
  network.SetCoupling(.4f);
  Measure("12-line resonator", [&](const std::size_t sample) {
    return network.Process(sample % 24000 == 0 ? 1.f : 0.f);
  });
}

void BenchmarkDispersion() {
  tfdsp::percussion::DispersionLoop loop;
  tfdsp::percussion::DispersionLoopParameters parameters;
  loop.Prepare(48000.f, 128.f, parameters);
  Measure("dispersion loop", [&](const std::size_t sample) {
    return loop.Process(sample % 24000 == 0 ? 1.f : 0.f);
  });
}

} // namespace

int main() {
  BenchmarkContactExciter();
  BenchmarkFrequencyShifter();
  BenchmarkResonators();
  BenchmarkDispersion();
}
