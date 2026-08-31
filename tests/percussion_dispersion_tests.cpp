#include "percussion_test_support.hpp"

#include "tfdsp/percussion/dispersion_loop.hpp"
#include "tfdsp/percussion/modulated_fractional_delay.hpp"
#include "tfdsp/percussion/self_phase_delay.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

void TestZeroDriveSelfPhaseReference() {
  tfdsp::percussion::SelfPhaseDelay nonlinear;
  tfdsp::percussion::ModulatedFractionalDelay reference;
  nonlinear.Prepare(48000.f, 64.f);
  reference.Prepare(64.f);
  tfdsp::percussion::SelfPhaseDelayParameters parameters;
  parameters.centreDelaySamples = 11.375f;
  parameters.maximumExcursionSamples = 3.f;
  parameters.drive = 0.f;
  nonlinear.SetParameters(parameters);
  double largestError = 0.0;
  for (std::size_t sample = 0; sample < 50000; ++sample) {
    const float input = percussion_test::Sine(sample, 4321.f, 48000.f);
    const float expected = reference.Process(input, 11.375f);
    const float actual = nonlinear.Process(input);
    largestError = std::max(largestError,
                            std::abs(static_cast<double>(actual - expected)));
  }
  CheckNear(largestError, 0.0, 1.e-7,
            "zero self-phase drive is the exact constant-delay reference");
}

tfdsp::percussion::DispersionLoopParameters LinearLoopParameters() {
  tfdsp::percussion::DispersionLoopParameters parameters;
  parameters.baseDelaySamples = 8.f;
  parameters.slowDelaySamples = 6.f;
  parameters.slowDepthSamples = 0.f;
  parameters.firstAllpassDelaySamples = 6.f;
  parameters.firstAllpassGain = 0.f;
  parameters.secondAllpassDelaySamples = 6.f;
  parameters.secondAllpassGain = 0.f;
  parameters.selfPhase.centreDelaySamples = 6.f;
  parameters.selfPhase.maximumExcursionSamples = 0.f;
  parameters.selfPhase.drive = 0.f;
  parameters.decay = {10.f, 10.f, 10.f};
  parameters.feedbackGain = .8f;
  return parameters;
}

void TestExplicitFeedbackCausality() {
  auto parameters = LinearLoopParameters();
  tfdsp::percussion::DispersionLoop loop;
  loop.Prepare(48000.f, 64.f, parameters);
  float first = 0.f;
  float second = 0.f;
  bool silentBeforePropagation = true;
  for (std::size_t sample = 0; sample <= 64; ++sample) {
    const float output = loop.Process(sample == 0 ? 1.f : 0.f);
    if (sample < 32)
      silentBeforePropagation = silentBeforePropagation && output == 0.f;
    if (sample == 32)
      first = output;
    if (sample == 64)
      second = output;
  }
  const float loss = tfdsp::percussion::ThreeBandDecayFilter::GainForT60(
      32.f / 48000.f, 10.f);
  Check(silentBeforePropagation,
        "dispersion loop has no direct or hidden early leakage");
  CheckNear(first, 1.0, 2.e-6,
            "dispersion impulse follows the declared serial propagation delay");
  CheckNear(second, .8 * loss, 3.e-6,
            "dispersion recurrence adds no hidden sample at the feedback sum");
}

void TestNonlinearStress() {
  tfdsp::percussion::DispersionLoop loop;
  tfdsp::percussion::DispersionLoopParameters parameters;
  parameters.feedbackGain = .995f;
  parameters.selfPhase.drive = 20.f;
  parameters.selfPhase.maximumExcursionSamples = 4.f;
  loop.Prepare(192000.f, 256.f, parameters);
  bool finite = true;
  float peak = 0.f;
  for (std::size_t sample = 0; sample < 1000000; ++sample) {
    const float drive = sample < 20000 ?
        percussion_test::Sine(sample, 17000.f, 192000.f) : 0.f;
    const float output = loop.Process(drive);
    finite = finite && std::isfinite(output);
    peak = std::max(peak, std::abs(output));
  }
  Check(finite, "nonlinear dispersion remains finite under long high feedback");
  Check(peak < 20.f, "nonlinear dispersion remains contractive without limiting");
  Check(std::isfinite(loop.Process(std::numeric_limits<float>::infinity())),
        "dispersion loop sanitizes non-finite drive");
}

} // namespace

int main() {
  TestZeroDriveSelfPhaseReference();
  TestExplicitFeedbackCausality();
  TestNonlinearStress();
  if (percussion_test::failures == 0)
    std::cout << "All percussion dispersion tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
