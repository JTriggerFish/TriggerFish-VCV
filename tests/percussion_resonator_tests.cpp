#include "percussion_test_support.hpp"

#include "tfdsp/percussion/coupled_resonator_network.hpp"
#include "tfdsp/percussion/orthogonal_mixer.hpp"
#include "tfdsp/percussion/resonator_submix.hpp"
#include "tfdsp/percussion/three_band_decay_filter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

template <std::size_t Size>
double Energy(const std::array<float, Size> &frame) {
  double energy = 0.0;
  for (const float value : frame)
    energy += value * value;
  return energy;
}

void TestOrthogonalMixer() {
  using Mixer = tfdsp::percussion::OrthogonalMixer<12>;
  Mixer mixer;
  Mixer::Frame input{};
  for (std::size_t index = 0; index < input.size(); ++index)
    input[index] = std::sin(static_cast<float>(index + 1));
  mixer.SetCoupling(0.f);
  Check(mixer.Process(input) == input, "zero coupling is exactly identity");
  for (const float coupling : {.1f, .37f, .8f, 1.f}) {
    mixer.SetCoupling(coupling);
    CheckNear(Energy(mixer.Process(input)), Energy(input), 3.e-6,
              "Givens coupling preserves frame energy");
  }
}

void TestDecayFilter() {
  using tfdsp::percussion::ThreeBandDecayFilter;
  ThreeBandDecayFilter filter;
  filter.Prepare(48000.f, 500.f, 5000.f);
  const tfdsp::percussion::ThreeBandDecayTimes equal{1.f, 1.f, 1.f};
  filter.SetDecayTimes(.01f, equal);
  const float gain = ThreeBandDecayFilter::GainForT60(.01f, 1.f);
  double largestError = 0.0;
  for (std::size_t sample = 0; sample < 10000; ++sample) {
    const float input = percussion_test::Sine(sample, 2173.f, 48000.f);
    largestError = std::max(largestError,
                            std::abs(static_cast<double>(
                                filter.Process(input) - gain * input)));
  }
  CheckNear(largestError, 0.0, 2.e-6,
            "equal three-band gains reconstruct scaled input exactly");
  CheckNear(ThreeBandDecayFilter::GainForT60(.25f, 1.f),
            std::pow(10.0, -.75), 1.e-7,
            "T60 maps to the expected per-path amplitude gain");
  Check(ThreeBandDecayFilter::GainForT60(
            .25f, std::numeric_limits<float>::infinity()) == 1.f,
        "infinite T60 is an explicit lossless traversal");
  filter.Reset();
  Check(filter.Process(std::numeric_limits<float>::denorm_min()) == 0.f,
        "decay filter flushes subnormal state to exact silence");

  filter.Reset();
  filter.SetDecayTimes(.01f, equal);
  const tfdsp::percussion::PassiveConstraintGains constraint{.8f, .5f, .5f,
                                                             .5f};
  largestError = 0.0;
  for (std::size_t sample = 0; sample < 10000; ++sample) {
    const float input = percussion_test::Sine(sample, 2173.f, 48000.f);
    largestError = std::max(largestError,
        std::abs(static_cast<double>(filter.Process(input, constraint) -
                                    .4f * gain * input)));
  }
  CheckNear(largestError, 0.0, 2.e-6,
            "equal constraint gains add exact passive loop attenuation");
}

void TestResonatorProjectionAndSubmix() {
  using Network = tfdsp::percussion::CoupledResonatorNetwork<3>;
  Network::Parameters parameters{};
  for (auto &line : parameters) {
    line.delaySamples = 12.f;
    line.decay = {10.f, 10.f, 10.f};
    line.inputGain = 1.f;
    line.outputGain = 1.f;
  }
  Network network;
  network.Prepare(48000.f, 64.f, parameters, 400.f, 5000.f);
  const Network::Frame excitation{1.f, 0.f, .25f};
  Network::Output output{};
  for (std::size_t sample = 0; sample <= 12; ++sample)
    output = network.ProcessProjected(sample == 0 ? 1.f : 0.f, excitation);
  CheckNear(output.lines[0], 1.0, 1.e-7,
            "projection excites the selected resonator line");
  CheckNear(output.lines[1], 0.0, 1.e-7,
            "zero projection leaves a resonator line unexcited");
  CheckNear(output.lines[2], .25, 1.e-7,
            "projection preserves continuous location weights");

  using Submix = tfdsp::percussion::ResonatorSubmix<3, 2>;
  Submix submix;
  Submix::Weights weights{{{{1.f, 1.f, 0.f}}, {{0.f, 0.f, 2.f}}}};
  submix.SetWeights(weights);
  const Submix::BusFrame buses = submix.Process(output.lines);
  CheckNear(buses[0], 1.0, 1.e-7,
            "resonator submix exposes an ordinary line group");
  CheckNear(buses[1], .5, 1.e-7,
            "resonator submix permits calibrated group weights");
}

void TestIndependentResonator() {
  using Network = tfdsp::percussion::CoupledResonatorNetwork<3>;
  Network::Parameters parameters{};
  for (auto &line : parameters) {
    line.inputGain = 0.f;
    line.outputGain = 0.f;
    line.delaySamples = 20.f;
    line.decay = {10.f, 10.f, 10.f};
  }
  parameters[0].inputGain = 1.f;
  parameters[0].outputGain = 1.f;
  Network network;
  network.Prepare(48000.f, 128.f, parameters, 400.f, 5000.f);
  network.SetCoupling(0.f);

  float firstReturn = 0.f;
  float secondReturn = 0.f;
  bool wetOnly = true;
  for (std::size_t sample = 0; sample <= 40; ++sample) {
    const float output = network.Process(sample == 0 ? 1.f : 0.f);
    if (sample < 20)
      wetOnly = wetOnly && output == 0.f;
    if (sample == 20)
      firstReturn = output;
    if (sample == 40)
      secondReturn = output;
  }
  const float expectedGain =
      tfdsp::percussion::ThreeBandDecayFilter::GainForT60(20.f / 48000.f, 10.f);
  Check(wetOnly, "resonator network has no dry leakage");
  CheckNear(firstReturn, 1.0, 1.e-7,
            "integer resonator returns at its configured delay");
  CheckNear(secondReturn, expectedGain, 2.e-6,
            "resonator recurrence follows its T60 gain");

  const float sanitized = network.Process(std::numeric_limits<float>::infinity());
  Check(std::isfinite(sanitized), "resonator sanitizes non-finite drive");
}

void TestStaticRetuneClearsEveryState() {
  using Network = tfdsp::percussion::CoupledResonatorNetwork<2>;
  Network::Parameters parameters{};
  for (std::size_t line = 0; line < parameters.size(); ++line) {
    parameters[line].delaySamples = 12.f + static_cast<float>(line);
    parameters[line].decay = {4.f, 3.f, 2.f};
  }
  Network network;
  network.Prepare(48000.f, 64.f, parameters, 400.f, 5000.f);
  network.SetCoupling(.7f);
  for (std::size_t sample = 0; sample < 256; ++sample)
    (void)network.Process(sample == 0 ? 1.f : 0.f);

  parameters[0].delaySamples = 19.25f;
  parameters[1].delaySamples = 23.5f;
  network.SetStaticParameters(parameters);
  double residual = 0.0;
  for (std::size_t sample = 0; sample < 256; ++sample)
    residual += std::abs(network.Process(0.f));
  CheckNear(residual, 0.0, 0.0,
            "static retuning clears delay and loss-filter state");
}

} // namespace

int main() {
  TestOrthogonalMixer();
  TestDecayFilter();
  TestResonatorProjectionAndSubmix();
  TestIndependentResonator();
  TestStaticRetuneClearsEveryState();
  if (percussion_test::failures == 0)
    std::cout << "All percussion resonator tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
