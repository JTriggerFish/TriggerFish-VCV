#include "percussion_test_support.hpp"

#include "tfdsp/percussion/contact_exciter.hpp"
#include "tfdsp/percussion/enveloped_noise_burst.hpp"
#include "tfdsp/percussion/finite_force_pulse.hpp"
#include "tfdsp/percussion/micro_contact_burst.hpp"
#include "tfdsp/percussion/spectral_tilt_filter.hpp"
#include "tfdsp/percussion/tonal_contact_chirp.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

template <typename Generator>
double RenderEnergy(Generator &generator) {
  double energy = 0.0;
  while (generator.Active()) {
    const float output = generator.Process();
    energy += output * output;
  }
  return energy;
}

void TestFiniteForcePulse() {
  tfdsp::percussion::FiniteForcePulse pulse;
  pulse.Prepare(48000.f);
  pulse.Trigger(.001f, .75f);
  std::size_t samples = 0;
  float peak = 0.f;
  bool nonnegative = true;
  while (pulse.Active()) {
    const float output = pulse.Process();
    peak = std::max(peak, output);
    nonnegative = nonnegative && output >= 0.f;
    ++samples;
  }
  Check(samples == 48, "half-sine pulse has its requested duration");
  CheckNear(peak, .75, .001, "half-sine pulse is peak normalized");
  Check(nonnegative, "half-sine contact applies compression in one direction");
  Check(pulse.Process() == 0.f, "half-sine pulse terminates exactly");

  pulse.Trigger(.001f, .25f);
  const double softEnergy = RenderEnergy(pulse);
  pulse.Trigger(.001f, 1.f);
  CheckNear(RenderEnergy(pulse) / softEnergy, 16.0, 1.e-4,
            "force-pulse energy follows squared hit strength");
}

void TestTonalContactChirp() {
  tfdsp::percussion::TonalContactChirp chirp;
  chirp.Prepare(48000.f);
  tfdsp::percussion::TonalContactChirpParameters parameters;
  parameters.durationSeconds = .006f;
  parameters.startFrequencyHz = 9000.f;
  parameters.endFrequencyHz = 1800.f;
  chirp.Trigger(parameters);
  std::size_t samples = 0;
  float first = 1.f;
  float last = 1.f;
  bool finite = true;
  while (chirp.Active()) {
    const float output = chirp.Process();
    if (samples == 0)
      first = output;
    last = output;
    finite = finite && std::isfinite(output);
    ++samples;
  }
  Check(samples == 288, "tonal contact chirp has its requested duration");
  CheckNear(first, 0.0, 1.e-7, "tonal contact chirp starts at zero");
  CheckNear(last, 0.0, 1.e-6, "tonal contact chirp ends at zero");
  Check(finite, "tonal contact chirp remains finite");

  parameters.startFrequencyHz = std::numeric_limits<float>::infinity();
  parameters.amplitude = std::numeric_limits<float>::quiet_NaN();
  chirp.Trigger(parameters);
  Check(std::isfinite(chirp.Process()),
        "tonal contact chirp sanitizes non-finite parameters");
}

double RenderMicroContacts(const std::uint32_t seed, const float amplitude,
                           float *capture = nullptr) {
  tfdsp::percussion::MicroContactBurst burst;
  burst.Prepare(48000.f);
  tfdsp::percussion::MicroContactBurstParameters parameters;
  parameters.durationSeconds = .05f;
  parameters.densityHz = 9000.f;
  parameters.amplitude = amplitude;
  parameters.seed = seed;
  burst.Trigger(parameters);
  double energy = 0.0;
  std::size_t sample = 0;
  while (burst.Active()) {
    const float output = burst.Process();
    if (capture)
      capture[sample] = output;
    energy += output * output;
    ++sample;
  }
  return energy;
}

void TestMicroContacts() {
  constexpr std::size_t Samples = 2400;
  float first[Samples]{};
  float repeated[Samples]{};
  float different[Samples]{};
  const double energy = RenderMicroContacts(1234, 1.f, first);
  RenderMicroContacts(1234, 1.f, repeated);
  RenderMicroContacts(1235, 1.f, different);
  bool identical = true;
  bool seedChangesDetail = false;
  bool bounded = true;
  for (std::size_t sample = 0; sample < Samples; ++sample) {
    identical = identical && first[sample] == repeated[sample];
    seedChangesDetail = seedChangesDetail || first[sample] != different[sample];
    bounded = bounded && std::abs(first[sample]) <= 1.f;
  }
  Check(identical, "micro-contact detail is repeatable for a fixed seed");
  Check(seedChangesDetail, "micro-contact seed changes only stochastic detail");
  Check(bounded, "micro-contact cluster has bounded amplitude");
  Check(energy > 1.e-5, "dense micro-contact cluster produces energy");
  CheckNear(RenderMicroContacts(1234, .25f) / energy, .0625, 1.e-6,
            "micro-contact energy follows squared hit strength");
}

void TestEnvelopedWhiteNoise() {
  tfdsp::percussion::EnvelopedNoiseBurst burst;
  burst.Prepare(48000.f);
  tfdsp::percussion::EnvelopedNoiseBurstParameters parameters;
  parameters.attackSeconds = 0.f;
  parameters.holdSeconds = 1.f;
  parameters.decaySeconds = 0.f;
  parameters.seed = 4321;
  burst.Trigger(parameters);
  double mean = 0.0;
  double energy = 0.0;
  double lagOne = 0.0;
  float previous = 0.f;
  std::size_t samples = 0;
  while (burst.Active()) {
    const float output = burst.Process();
    mean += output;
    energy += output * output;
    if (samples > 0)
      lagOne += output * previous;
    previous = output;
    ++samples;
  }
  Check(samples == 48000, "white-noise burst has exact envelope duration");
  Check(std::abs(mean / static_cast<double>(samples)) < .01,
        "white-noise burst has negligible DC bias");
  Check(std::abs(lagOne / energy) < .015,
        "white-noise burst has negligible one-sample correlation");

  parameters.attackSeconds = .001f;
  parameters.holdSeconds = .002f;
  parameters.decaySeconds = .003f;
  burst.Trigger(parameters);
  samples = 0;
  while (burst.Active()) {
    Check(std::abs(burst.Process()) <= 1.f,
          "white-noise envelope remains bounded");
    ++samples;
  }
  Check(samples == 288, "shaped white-noise burst sums attack, hold, and decay");
  Check(burst.Process() == 0.f, "white-noise burst terminates exactly");
}

void TestNoiseTilt() {
  tfdsp::percussion::SpectralTiltFilter neutral;
  neutral.Prepare(48000.f);
  neutral.SetTilt(0.f, 3000.f);
  double largestError = 0.0;
  for (std::size_t sample = 0; sample < 10000; ++sample) {
    const float input = percussion_test::Sine(sample, 7123.f, 48000.f);
    largestError = std::max(largestError,
                            std::abs(static_cast<double>(neutral.Process(input) - input)));
  }
  CheckNear(largestError, 0.0, 2.e-7,
            "zero noise tilt is an exact complementary reconstruction");

  tfdsp::percussion::SpectralTiltFilter dark;
  tfdsp::percussion::SpectralTiltFilter bright;
  dark.Prepare(48000.f);
  bright.Prepare(48000.f);
  dark.SetTilt(-12.f, 3000.f);
  bright.SetTilt(12.f, 3000.f);
  double darkHighEnergy = 0.0;
  double brightHighEnergy = 0.0;
  for (std::size_t sample = 0; sample < 20000; ++sample) {
    const float input = sample % 2 == 0 ? 1.f : -1.f;
    const float darkOutput = dark.Process(input);
    const float brightOutput = bright.Process(input);
    if (sample > 100)
      darkHighEnergy += darkOutput * darkOutput;
    if (sample > 100)
      brightHighEnergy += brightOutput * brightOutput;
  }
  Check(brightHighEnergy > 4.0 * darkHighEnergy,
        "positive noise tilt raises high-band energy relative to negative tilt");

  neutral.Reset();
  Check(neutral.Process(std::numeric_limits<float>::denorm_min()) == 0.f &&
            neutral.Process(0.f) == 0.f,
        "contact tilt filter flushes subnormal input and state");
}

void TestContactRouting() {
  tfdsp::percussion::ContactExciter exciter;
  exciter.Prepare(48000.f);
  tfdsp::percussion::ContactExciterParameters parameters;
  parameters.pulseDurationSeconds = .001f;
  parameters.pulseAmplitude = .5f;
  parameters.chirp.amplitude = 0.f;
  parameters.noise.amplitude = 0.f;
  parameters.microContacts.amplitude = 0.f;
  parameters.routing = {};
  parameters.routing.pulseDirect = 2.f;
  parameters.routing.pulseBody = 3.f;
  exciter.Trigger(parameters);
  bool separated = true;
  std::size_t activeSamples = 0;
  while (exciter.Active()) {
    const auto output = exciter.Process();
    separated = separated &&
        std::abs(output.bodyDrive - 1.5f * output.directRadiation) < 1.e-6f;
    ++activeSamples;
  }
  Check(separated, "contact exciter keeps direct and body routing explicit");
  Check(activeSamples >= 48,
        "contact exciter remains active until its longest primitive finishes");
  const auto silence = exciter.Process();
  Check(silence.directRadiation == 0.f && silence.bodyDrive == 0.f,
        "contact exciter returns to exact silence");
}

void TestContactSampleRatesAndBounds() {
  for (const float sampleRate : {44100.f, 48000.f, 88200.f, 96000.f, 192000.f}) {
    tfdsp::percussion::FiniteForcePulse pulse;
    pulse.Prepare(sampleRate);
    pulse.Trigger(.001f, std::numeric_limits<float>::max());
    std::size_t samples = 0;
    float peak = 0.f;
    while (pulse.Active()) {
      peak = std::max(peak, pulse.Process());
      ++samples;
    }
    Check(samples == static_cast<std::size_t>(std::lround(.001f * sampleRate)),
          "contact duration is equivalent at every supported sample rate");
    Check(std::isfinite(peak) && peak <= 16.f,
          "contact amplitude remains bounded for every finite control value");
    pulse.Trigger(.1f, 0.f);
    Check(!pulse.Active(), "zero-amplitude contact does no idle processing");
    pulse.Trigger(.1f, std::numeric_limits<float>::denorm_min());
    Check(!pulse.Active(),
          "subnormal contact strength is treated as exact silence");
  }

  tfdsp::percussion::ContactExciter exciter;
  exciter.Prepare(48000.f);
  tfdsp::percussion::ContactExciterParameters silent;
  silent.pulseAmplitude = 0.f;
  silent.chirp.amplitude = 0.f;
  silent.noise.amplitude = 0.f;
  silent.microContacts.amplitude = 0.f;
  exciter.Trigger(silent);
  Check(!exciter.Active(),
        "a zero-amplitude composite contact is immediately inactive");
}

} // namespace

int main() {
  TestFiniteForcePulse();
  TestTonalContactChirp();
  TestMicroContacts();
  TestEnvelopedWhiteNoise();
  TestNoiseTilt();
  TestContactRouting();
  TestContactSampleRatesAndBounds();
  if (percussion_test::failures == 0)
    std::cout << "All percussion contact tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
