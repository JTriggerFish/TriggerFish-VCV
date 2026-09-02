#include "tfdsp/percussion/contact_exciter.hpp"
#include "tfdsp/percussion/correlated_fm_burst.hpp"
#include "tfdsp/percussion/coupled_resonator_network.hpp"
#include "tfdsp/percussion/crash_cymbal.hpp"
#include "tfdsp/percussion/dispersion_loop.hpp"
#include "tfdsp/percussion/frequency_shifter.hpp"
#include "tfdsp/percussion/micro_contact_process.hpp"
#include "tfdsp/percussion/modal_bank.hpp"
#include "tfdsp/percussion/observation_model.hpp"
#include "tfdsp/percussion/statistical_modal_cloud.hpp"

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

  Measure("automated frequency shifter", [&](const std::size_t sample) {
    const float phase = static_cast<float>(sample % 96000) / 95999.f;
    shifter.SetShiftHz(-4000.f + 8000.f * phase);
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

void BenchmarkCorrelatedFm() {
  tfdsp::percussion::CorrelatedFmBurst burst;
  burst.Prepare(48000.f);
  tfdsp::percussion::CorrelatedFmBurstParameters parameters;
  parameters.amplitude.segments[0] = {1.f, .001f};
  parameters.amplitude.segments[1] = {0.f, .05f};
  parameters.amplitude.segmentCount = 2;
  parameters.carrierFrequencyHz.initialValue = 80.f;
  parameters.frequencyDeviationHz.initialValue = 2000.f;
  Measure("correlated FM burst", [&](const std::size_t sample) {
    if (sample % 4096 == 0)
      burst.Trigger(parameters);
    return burst.Process();
  });
}

void BenchmarkMicroContacts() {
  tfdsp::percussion::MicroContactProcess process;
  process.Prepare(48000.f);
  tfdsp::percussion::MicroContactProcessParameters parameters;
  parameters.densityHz = 8000.f;
  process.StartStream(parameters);
  Measure("micro-contact stream", [&](const std::size_t) {
    return process.Process();
  });
}

void BenchmarkObservation() {
  using Model = tfdsp::percussion::ObservationModel<4>;
  Model model;
  Model::Parameters parameters{};
  model.Prepare(48000.f, .1f, parameters);
  Measure("four-source observation", [&](const std::size_t sample) {
    const float input = std::sin(.071f * static_cast<float>(sample));
    return model.Process({input, .5f * input, -.25f * input, .1f * input});
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

void BenchmarkModalBanks() {
  using SparseBank = tfdsp::percussion::ModalBank<24>;
  using DenseBank = tfdsp::percussion::ModalBank<256>;
  using DenserBank = tfdsp::percussion::ModalBank<512>;
  using DensestBank = tfdsp::percussion::ModalBank<1024>;
  tfdsp::percussion::StatisticalModalCloudParameters cloud;
  const auto denseParameters =
      tfdsp::percussion::MakeStatisticalModalCloud<256>(48000.f, cloud);
  const auto denserParameters =
      tfdsp::percussion::MakeStatisticalModalCloud<512>(48000.f, cloud);
  const auto densestParameters =
      tfdsp::percussion::MakeStatisticalModalCloud<1024>(48000.f, cloud);
  const auto sparseParameters =
      tfdsp::percussion::MakeStatisticalModalCloud<24>(48000.f, cloud);
  SparseBank sparse;
  DenseBank dense;
  DenserBank denser;
  DensestBank densest;
  sparse.Prepare(48000.f, sparseParameters, 700.f, 6500.f);
  dense.Prepare(48000.f, denseParameters, 700.f, 6500.f);
  denser.Prepare(48000.f, denserParameters, 700.f, 6500.f);
  densest.Prepare(48000.f, densestParameters, 700.f, 6500.f);
  Measure("24-mode sparse bank", [&](const std::size_t sample) {
    return sparse.Process(sample % 24000 == 0 ? 1.f : 0.f);
  });
  Measure("256-mode dense cloud", [&](const std::size_t sample) {
    return dense.Process(sample % 24000 == 0 ? 1.f : 0.f);
  });
  Measure("512-mode dense cloud", [&](const std::size_t sample) {
    return denser.Process(sample % 24000 == 0 ? 1.f : 0.f);
  });
  Measure("1024-mode dense cloud", [&](const std::size_t sample) {
    return densest.Process(sample % 24000 == 0 ? 1.f : 0.f);
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

void BenchmarkCrashVariant(const std::string_view name,
                           tfdsp::percussion::CrashCymbalFitParameters fit) {
  tfdsp::percussion::CrashCymbal cymbal;
  cymbal.Prepare(
      48000.f, tfdsp::percussion::DefaultCrashCymbalParameters(48000.f, fit));
  const tfdsp::percussion::CrashCymbalHit hit{};
  Measure(name, [&](const std::size_t sample) {
    if (sample % 24000 == 0)
      cymbal.Trigger(hit);
    return cymbal.Process();
  });
}

void BenchmarkCrashCymbal() {
  tfdsp::percussion::CrashCymbalFitParameters fit;
  BenchmarkCrashVariant("crash without modal banks", [&] {
    auto result = fit;
    result.sparseGain = 0.f;
    result.denseGain = 0.f;
    return result;
  }());
  BenchmarkCrashVariant("crash with sparse bank", [&] {
    auto result = fit;
    result.denseGain = 0.f;
    return result;
  }());
  BenchmarkCrashVariant("crash with dense bank", [&] {
    auto result = fit;
    result.sparseGain = 0.f;
    return result;
  }());
  BenchmarkCrashVariant("complete crash", fit);
  fit.denseModeDensity = 2.f;
  BenchmarkCrashVariant("complete crash, 4096 wash modes", fit);
}

} // namespace

int main() {
  BenchmarkContactExciter();
  BenchmarkCorrelatedFm();
  BenchmarkMicroContacts();
  BenchmarkFrequencyShifter();
  BenchmarkResonators();
  BenchmarkModalBanks();
  BenchmarkDispersion();
  BenchmarkObservation();
  BenchmarkCrashCymbal();
}
