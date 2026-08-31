#include "percussion_test_support.hpp"

#include "tfdsp/percussion/biquad.hpp"
#include "tfdsp/percussion/biquad_design.hpp"
#include "tfdsp/percussion/radiation_filter.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

template <typename Filter>
double ToneGain(Filter &filter, const float frequencyHz, const float sampleRate) {
  double inputEnergy = 0.0;
  double outputEnergy = 0.0;
  for (std::size_t sample = 0; sample < 65536; ++sample) {
    const float input = percussion_test::Sine(sample, frequencyHz, sampleRate);
    const float output = filter.Process(input);
    if (sample < 4096)
      continue;
    inputEnergy += input * input;
    outputEnergy += output * output;
  }
  return std::sqrt(outputEnergy / inputEnergy);
}

void TestBiquadResponses() {
  for (const float sampleRate : {44100.f, 48000.f, 96000.f, 192000.f}) {
    tfdsp::percussion::Biquad peak;
    peak.SetCoefficients(tfdsp::percussion::biquad_design::Peaking(
        5000.f, 1.2f, 9.f, sampleRate));
    CheckNear(ToneGain(peak, 5000.f, sampleRate),
              std::pow(10.0, 9.0 / 20.0), .006,
              "peaking biquad reaches its requested centre gain");

    tfdsp::percussion::Biquad lowpass;
    lowpass.SetCoefficients(tfdsp::percussion::biquad_design::Lowpass(
        6000.f, .70710678f, sampleRate));
    Check(ToneGain(lowpass, 1000.f, sampleRate) > .99,
          "radiation low-pass preserves its lower passband");
    Check(ToneGain(lowpass, .4f * sampleRate, sampleRate) < .2,
          "radiation low-pass rejects near-Nyquist content");
  }
}

void TestRadiationChain() {
  tfdsp::percussion::RadiationFilter filter;
  tfdsp::percussion::RadiationFilterParameters parameters;
  parameters.lowCutHz = 300.f;
  parameters.colourFrequencyHz = 4500.f;
  parameters.colourGainDb = 6.f;
  parameters.highCutHz = 12000.f;
  filter.Prepare(48000.f, parameters);
  const double low = ToneGain(filter, 60.f, 48000.f);
  filter.Reset();
  const double colour = ToneGain(filter, 4500.f, 48000.f);
  filter.Reset();
  const double high = ToneGain(filter, 20000.f, 48000.f);
  Check(colour > 3.0 * low && colour > 3.0 * high,
        "radiation chain separates body colour from spectral extremes");
  Check(std::isfinite(filter.Process(std::numeric_limits<float>::infinity())),
        "radiation chain sanitizes non-finite input");
}

} // namespace

int main() {
  TestBiquadResponses();
  TestRadiationChain();
  if (percussion_test::failures == 0)
    std::cout << "All percussion radiation tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
