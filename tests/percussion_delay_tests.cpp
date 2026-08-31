#include "percussion_test_support.hpp"

#include "tfdsp/percussion/modulated_fractional_delay.hpp"
#include "tfdsp/percussion/schroeder_allpass.hpp"
#include "tfdsp/percussion/static_fractional_delay.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

template <typename Processor>
double ToneGain(Processor &processor, const float frequencyHz,
                const float sampleRate) {
  double inputEnergy = 0.0;
  double outputEnergy = 0.0;
  for (std::size_t sample = 0; sample < 65536; ++sample) {
    const float input = percussion_test::Sine(sample, frequencyHz, sampleRate);
    const float output = processor.Process(input);
    if (sample < 4096)
      continue;
    inputEnergy += input * input;
    outputEnergy += output * output;
  }
  return std::sqrt(outputEnergy / inputEnergy);
}

void TestStaticDelay() {
  tfdsp::percussion::StaticFractionalDelay delay;
  delay.Prepare(128.f, 23.375f);
  const auto gain = ToneGain(delay, 7300.f, 48000.f);
  CheckNear(gain, 1.0, 2.e-4,
            "Thiran static delay preserves sinusoid magnitude");

  delay.Reset();
  double error = 0.0;
  double referenceEnergy = 0.0;
  for (std::size_t sample = 0; sample < 48000; ++sample) {
    const float input = percussion_test::Sine(sample, 127.f, 48000.f);
    const float output = delay.Process(input);
    if (sample < 2048)
      continue;
    const double reference = std::sin(
        6.283185307179586 * 127.0 *
        (static_cast<double>(sample) - 23.375) / 48000.0);
    error += (output - reference) * (output - reference);
    referenceEnergy += reference * reference;
  }
  CheckNear(std::sqrt(error / referenceEnergy), 0.0, 2.e-4,
            "Thiran static delay has the requested low-frequency phase delay");

  delay.SetDelaySamples(7.f);
  Check(delay.Process(0.f) == 0.f,
        "changing a static delay clears its previous state");

  delay.SetDelaySamples(7.5f);
  bool exactSilence = true;
  for (std::size_t sample = 0; sample < 32; ++sample)
    exactSilence = exactSilence &&
        delay.Process(sample == 0 ? std::numeric_limits<float>::denorm_min()
                                  : 0.f) == 0.f;
  Check(exactSilence,
        "fractional-delay storage and Thiran state flush subnormals");
}

void TestMovingDelay() {
  tfdsp::percussion::ModulatedFractionalDelay delay;
  delay.Prepare(128.f);
  double error = 0.0;
  double referenceEnergy = 0.0;
  float previous = 0.f;
  float largestStep = 0.f;
  for (std::size_t sample = 0; sample < 96000; ++sample) {
    const float time = static_cast<float>(sample) / 48000.f;
    const float tap = 31.f + 2.2f * std::sin(6.283185307179586f * .37f * time);
    const float input = percussion_test::Sine(sample, 911.f, 48000.f);
    const float output = delay.Process(input, tap);
    if (sample >= 4096) {
      const double phase = 6.283185307179586 * 911.0 *
                           (static_cast<double>(sample) - tap) / 48000.0;
      const double reference = std::sin(phase);
      error += (output - reference) * (output - reference);
      referenceEnergy += reference * reference;
      largestStep = std::max(largestStep, std::abs(output - previous));
    }
    previous = output;
  }
  CheckNear(std::sqrt(error / referenceEnergy), 0.0, 8.e-4,
            "moving sinc delay follows a continuously modulated analytic tone");
  Check(largestStep < .14f,
        "moving sinc delay crosses integer boundaries without spikes");

  const float sanitized = delay.Process(
      std::numeric_limits<float>::quiet_NaN(),
      std::numeric_limits<float>::infinity());
  Check(std::isfinite(sanitized), "moving delay sanitizes non-finite input");
}

void TestSchroederAllpass() {
  tfdsp::percussion::SchroederAllpass allpass;
  allpass.Prepare(64.f, 17.625f, .73f);
  const auto gain = ToneGain(allpass, 5300.f, 48000.f);
  CheckNear(gain, 1.0, 3.e-3,
            "fractional Schroeder allpass preserves sinusoid magnitude");

  allpass.Reset();
  double energy = 0.0;
  bool finite = true;
  for (std::size_t sample = 0; sample < 20000; ++sample) {
    const float output = allpass.Process(sample == 0 ? 1.f : 0.f);
    energy += output * output;
    finite = finite && std::isfinite(output);
  }
  Check(finite, "fractional Schroeder allpass remains finite");
  CheckNear(energy, 1.0, 2.e-3,
            "fractional Schroeder allpass preserves impulse energy");
}

void TestDelayAtEverySupportedRate() {
  for (const float sampleRate : {44100.f, 48000.f, 88200.f, 96000.f, 192000.f}) {
    const float scale = sampleRate / 48000.f;
    tfdsp::percussion::StaticFractionalDelay delay;
    delay.Prepare(512.f, 23.375f * scale);
    const auto gain = ToneGain(delay, .12f * sampleRate, sampleRate);
    CheckNear(gain, 1.0, 3.e-4,
              "static fractional delay preserves gain at supported rates");

    tfdsp::percussion::ModulatedFractionalDelay moving;
    moving.Prepare(512.f);
    bool finite = true;
    for (std::size_t sample = 0; sample < 8192; ++sample) {
      const float tap = scale * (31.f + 2.2f * std::sin(
          6.283185307179586f * .37f * sample / sampleRate));
      finite = finite && std::isfinite(moving.Process(
          percussion_test::Sine(sample, .08f * sampleRate, sampleRate), tap));
    }
    Check(finite, "moving delay remains finite at every supported sample rate");
  }
}

} // namespace

int main() {
  TestStaticDelay();
  TestMovingDelay();
  TestSchroederAllpass();
  TestDelayAtEverySupportedRate();
  if (percussion_test::failures == 0)
    std::cout << "All percussion delay tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
