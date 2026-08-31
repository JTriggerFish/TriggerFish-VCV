#include "percussion_test_support.hpp"

#include "tfdsp/percussion/correlated_fm_burst.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

tfdsp::percussion::CorrelatedFmBurstParameters Parameters(
    const std::uint32_t seed = 1) {
  using Curve = tfdsp::percussion::TrajectoryCurve;
  tfdsp::percussion::CorrelatedFmBurstParameters parameters;
  parameters.amplitude.segments[0] = {1.f, .001f, Curve::Linear};
  parameters.amplitude.segments[1] = {.7f, .009f, Curve::Geometric};
  parameters.amplitude.segments[2] = {.001f, .04f, Curve::Geometric};
  parameters.amplitude.segments[3] = {0.f, .001f, Curve::Linear};
  parameters.amplitude.segmentCount = 4;
  parameters.carrierFrequencyHz.initialValue = 220.f;
  parameters.carrierFrequencyHz.segments[0] = {65.f, .05f, Curve::Geometric};
  parameters.carrierFrequencyHz.segmentCount = 1;
  parameters.frequencyDeviationHz.initialValue = 2400.f;
  parameters.frequencyDeviationHz.segments[0] = {20.f, .012f, Curve::Geometric};
  parameters.frequencyDeviationHz.segments[1] = {0.f, .001f, Curve::Linear};
  parameters.frequencyDeviationHz.segmentCount = 2;
  parameters.irregularCutoffHz = 7000.f;
  parameters.seed = seed;
  return parameters;
}

template <typename Burst>
std::vector<float> Render(const std::uint32_t seed, const bool periodic = false) {
  auto parameters = Parameters(seed);
  if (periodic) {
    parameters.carrierFrequencyHz.initialValue = 5000.f;
    parameters.carrierFrequencyHz.segmentCount = 0;
    parameters.frequencyDeviationHz.initialValue = 1800.f;
    parameters.frequencyDeviationHz.segmentCount = 0;
    parameters.periodicMix = 1.f;
    parameters.periodicModulatorHz = 1300.f;
  }
  Burst burst;
  burst.Prepare(48000.f);
  burst.Trigger(parameters);
  std::vector<float> output;
  while (burst.Active())
    output.push_back(burst.Process());
  return output;
}

void TestDeterministicFiniteBurst() {
  const auto first = Render<tfdsp::percussion::CorrelatedFmBurst>(4321);
  const auto repeated = Render<tfdsp::percussion::CorrelatedFmBurst>(4321);
  const auto different = Render<tfdsp::percussion::CorrelatedFmBurst>(4322);
  Check(first == repeated, "FM burst repeats exactly for a fixed seed");
  Check(first.size() == different.size() && first != different,
        "FM seed changes detail without changing event duration");
  double energy = 0.0;
  bool finite = true;
  for (const float sample : first) {
    energy += sample * sample;
    finite = finite && std::isfinite(sample);
  }
  Check(first.size() >= 2400 && finite && energy > 1.e-4,
        "FM burst renders a finite audible envelope and decimator tail");
}

template <std::size_t BinCount>
std::array<double, BinCount> Magnitude(const std::vector<float> &signal) {
  constexpr double TwoPi = 6.283185307179586476925286766559;
  std::array<double, BinCount> magnitude{};
  for (std::size_t bin = 0; bin < BinCount; ++bin) {
    double real = 0.0;
    double imaginary = 0.0;
    for (std::size_t sample = 0; sample < signal.size(); ++sample) {
      const double window = .5 - .5 * std::cos(
          TwoPi * sample / static_cast<double>(signal.size() - 1));
      const double phase = TwoPi * bin * sample / (2.0 * BinCount);
      real += window * signal[sample] * std::cos(phase);
      imaginary -= window * signal[sample] * std::sin(phase);
    }
    magnitude[bin] = std::hypot(real, imaginary);
  }
  return magnitude;
}

void TestOversamplingReference() {
  constexpr std::size_t BinCount = 512;
  const auto production = Render<tfdsp::percussion::CorrelatedFmBurst>(77, true);
  const auto reference =
      Render<tfdsp::percussion::CorrelatedFmBurstReference4x>(77, true);
  Check(production.size() == reference.size(),
        "FM oversampling references have equal host durations");
  const auto productionMagnitude = Magnitude<BinCount>(production);
  const auto referenceMagnitude = Magnitude<BinCount>(reference);
  const double referencePeak = *std::max_element(
      referenceMagnitude.begin(), referenceMagnitude.end());
  const double productionPeak = *std::max_element(
      productionMagnitude.begin(), productionMagnitude.end());
  double squaredDbError = 0.0;
  std::size_t comparedBins = 0;
  for (std::size_t bin = 0; bin < BinCount; ++bin) {
    if (referenceMagnitude[bin] < 1.e-3 * referencePeak)
      continue;
    const double referenceLevel = referenceMagnitude[bin] / referencePeak;
    const double productionLevel = productionMagnitude[bin] / productionPeak;
    const double errorDb = 20.0 * std::log10(
        std::max(productionLevel, 1.e-8) / referenceLevel);
    squaredDbError += errorDb * errorDb;
    ++comparedBins;
  }
  Check(comparedBins > 8, "FM comparison covers audible sidebands");
  CheckNear(std::sqrt(squaredDbError / comparedBins), 0.0, 1.0,
            "2x FM spectrum converges to its 4x offline reference");
}

void TestSupportedSampleRates() {
  for (const float sampleRate :
       {44100.f, 48000.f, 88200.f, 96000.f, 192000.f}) {
    tfdsp::percussion::CorrelatedFmBurst burst;
    burst.Prepare(sampleRate);
    burst.Trigger(Parameters(12));
    std::size_t samples = 0;
    double energy = 0.0;
    bool finite = true;
    while (burst.Active() && samples < static_cast<std::size_t>(sampleRate)) {
      const float output = burst.Process();
      energy += output * output;
      finite = finite && std::isfinite(output);
      ++samples;
    }
    const double duration = samples / sampleRate;
    Check(finite && energy > 1.e-5,
          "FM burst remains finite and audible at every supported rate");
    Check(duration > .050 && duration < .054,
          "FM envelope duration is sample-rate equivalent");
  }
}

} // namespace

int main() {
  TestDeterministicFiniteBurst();
  TestOversamplingReference();
  TestSupportedSampleRates();
  if (percussion_test::failures == 0)
    std::cout << "All percussion FM tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
