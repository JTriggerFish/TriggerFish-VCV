#include "percussion_test_support.hpp"

#include "tfdsp/percussion/micro_contact_process.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace {

using percussion_test::Check;

std::vector<float> RenderCluster(const std::uint32_t seed) {
  tfdsp::percussion::MicroContactProcess process;
  process.Prepare(48000.f);
  tfdsp::percussion::MicroContactProcessParameters parameters;
  parameters.densityHz = 5000.f;
  parameters.contactDurationSeconds = .0008f;
  parameters.releaseSeconds = .004f;
  parameters.seed = seed;
  process.TriggerCluster(parameters, .025f);
  std::vector<float> output;
  while (process.Active() && output.size() < 4096)
    output.push_back(process.Process());
  return output;
}

void TestFiniteCluster() {
  const auto first = RenderCluster(91);
  const auto repeated = RenderCluster(91);
  const auto different = RenderCluster(92);
  Check(first == repeated,
        "micro-contact cluster repeats exactly for a fixed seed");
  Check(first.size() == different.size() && first != different,
        "micro-contact seed changes detail without changing duration");
  double energy = 0.0;
  bool finite = true;
  for (const float sample : first) {
    energy += sample * sample;
    finite = finite && std::isfinite(sample);
  }
  Check(finite && energy > 1.e-5 && first.size() < 4096,
        "finite micro-contact cluster terminates with bounded energy");
}

void TestContinuousGate() {
  tfdsp::percussion::MicroContactProcess stream;
  stream.Prepare(48000.f);
  tfdsp::percussion::MicroContactProcessParameters parameters;
  parameters.densityHz = 8000.f;
  parameters.attackSeconds = .002f;
  parameters.releaseSeconds = .01f;
  stream.StartStream(parameters);
  for (std::size_t sample = 0; sample < 2400; ++sample)
    (void)stream.Process();
  Check(stream.Active(), "micro-contact stream remains active while held");
  stream.SetDensityHz(12000.f);
  stream.SetAmplitude(.5f);
  stream.SetTilt(6.f, 5000.f);
  bool automatedFinite = true;
  double automatedEnergy = 0.0;
  for (std::size_t sample = 0; sample < 2400; ++sample) {
    const float output = stream.Process();
    automatedFinite = automatedFinite && std::isfinite(output);
    automatedEnergy += output * output;
  }
  Check(automatedFinite && automatedEnergy > 1.e-5 && stream.Active(),
        "micro-contact controls automate without restarting the held stream");
  stream.StopStream();
  std::size_t releaseSamples = 0;
  while (stream.Active() && releaseSamples < 9600) {
    (void)stream.Process();
    ++releaseSamples;
  }
  Check(releaseSamples > 0 && releaseSamples < 9600,
        "micro-contact gate releases to exact inactivity");
  Check(stream.Process() == 0.f,
        "released micro-contact stream returns exact silence");
}

void TestSupportedSampleRates() {
  for (const float sampleRate :
       {44100.f, 48000.f, 88200.f, 96000.f, 192000.f}) {
    tfdsp::percussion::MicroContactProcess process;
    process.Prepare(sampleRate);
    tfdsp::percussion::MicroContactProcessParameters parameters;
    parameters.densityHz = 5000.f;
    parameters.contactDurationSeconds = .0008f;
    parameters.releaseSeconds = .004f;
    process.TriggerCluster(parameters, .025f);
    std::size_t samples = 0;
    double energy = 0.0;
    bool finite = true;
    while (process.Active() && samples < static_cast<std::size_t>(sampleRate)) {
      const float output = process.Process();
      energy += output * output;
      finite = finite && std::isfinite(output);
      ++samples;
    }
    Check(finite && energy > 1.e-5 && samples < sampleRate / 10,
          "micro-contact gesture is bounded at every supported rate");
  }
}

} // namespace

int main() {
  TestFiniteCluster();
  TestContinuousGate();
  TestSupportedSampleRates();
  if (percussion_test::failures == 0)
    std::cout << "All percussion micro-contact tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
