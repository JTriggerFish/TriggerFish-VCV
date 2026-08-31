#include "percussion_test_support.hpp"

#include "tfdsp/percussion/observation_delay.hpp"
#include "tfdsp/percussion/observation_model.hpp"

#include <array>
#include <cmath>
#include <limits>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

void TestZeroAndIntegerObservationDelay() {
  tfdsp::percussion::ObservationDelay delay;
  delay.Prepare(64.f, 0.f);
  for (std::size_t sample = 0; sample < 64; ++sample) {
    const float input = percussion_test::Sine(sample, 317.f, 48000.f);
    CheckNear(delay.Process(input), input, 0.0,
              "zero observation delay is an exact wire");
  }

  delay.SetStaticDelaySamples(7.f);
  for (std::size_t sample = 0; sample < 16; ++sample) {
    const float output = delay.Process(sample == 0 ? 1.f : 0.f);
    CheckNear(output, sample == 7 ? 1.0 : 0.0, 0.0,
              "integer observation delay returns at its declared sample");
  }
}

void TestObservationRouting() {
  using Model = tfdsp::percussion::ObservationModel<3>;
  Model::Parameters parameters{};
  for (auto &path : parameters) {
    path.gain = 0.f;
    path.radiationEnabled = false;
  }
  parameters[0].gain = 2.f;
  parameters[1].gain = .5f;
  parameters[1].invertPolarity = true;
  Model model;
  model.Prepare(48000.f, .01f, parameters);
  const Model::SourceFrame sources{.25f, 1.f, 100.f};
  CheckNear(model.Process(sources), 0.0, 1.e-7,
            "observation gains, isolation, and polarity are explicit");

  parameters[0].gain = std::numeric_limits<float>::infinity();
  parameters[1].gain = std::numeric_limits<float>::quiet_NaN();
  model.SetStaticParameters(parameters);
  Check(std::isfinite(model.Process(sources)),
        "observation model sanitizes non-finite path gains");
}

void TestObservationPathDelay() {
  using Model = tfdsp::percussion::ObservationModel<2>;
  Model::Parameters parameters{};
  for (auto &path : parameters) {
    path.gain = 1.f;
    path.radiationEnabled = false;
  }
  parameters[1].delaySeconds = 4.f / 48000.f;
  Model model;
  model.Prepare(48000.f, .01f, parameters);
  float immediate = 0.f;
  float delayed = 0.f;
  for (std::size_t sample = 0; sample < 8; ++sample) {
    const Model::SourceFrame sources{sample == 0 ? 1.f : 0.f,
                                     sample == 0 ? 1.f : 0.f};
    const float output = model.Process(sources);
    if (sample == 0)
      immediate = output;
    if (sample == 4)
      delayed = output;
  }
  CheckNear(immediate, 1.0, 0.0,
            "direct observation source has no forced latency");
  CheckNear(delayed, 1.0, 0.0,
            "relative observation delay is preserved per source");
}

void TestSupportedSampleRates() {
  for (const float sampleRate :
       {44100.f, 48000.f, 88200.f, 96000.f, 192000.f}) {
    tfdsp::percussion::ObservationDelay delay;
    delay.Prepare(.01f * sampleRate, 0.f);
    bool exact = true;
    for (std::size_t sample = 0; sample < 64; ++sample) {
      const float input = percussion_test::Sine(sample, 997.f, sampleRate);
      exact = exact && delay.Process(input) == input;
    }
    Check(exact, "zero observation delay remains exact at every supported rate");
  }
}

} // namespace

int main() {
  TestZeroAndIntegerObservationDelay();
  TestObservationRouting();
  TestObservationPathDelay();
  TestSupportedSampleRates();
  if (percussion_test::failures == 0)
    std::cout << "All percussion observation tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
