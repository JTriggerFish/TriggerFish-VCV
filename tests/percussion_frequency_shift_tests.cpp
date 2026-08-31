#include "percussion_test_support.hpp"

#include "tfdsp/percussion/frequency_shifter.hpp"
#include "tfdsp/percussion/quadrature_oscillator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

double ToneMagnitude(const float *signal, const std::size_t count,
                     const float frequencyHz, const float sampleRate) {
  constexpr double TwoPi = 6.283185307179586476925286766559;
  double real = 0.0;
  double imaginary = 0.0;
  for (std::size_t sample = 0; sample < count; ++sample) {
    const double phase = TwoPi * frequencyHz *
                         static_cast<double>(sample) / sampleRate;
    real += signal[sample] * std::cos(phase);
    imaginary -= signal[sample] * std::sin(phase);
  }
  return 2.0 * std::hypot(real, imaginary) / static_cast<double>(count);
}

void CheckShift(const float sampleRate, const float inputHz,
                const float shiftHz) {
  constexpr std::size_t Warmup = 4096;
  constexpr std::size_t Analysis = 65536;
  static float output[Analysis];
  tfdsp::percussion::FrequencyShifter shifter;
  shifter.Prepare(sampleRate);
  shifter.SetShiftHz(shiftHz);
  for (std::size_t sample = 0; sample < Warmup + Analysis; ++sample) {
    const float input = percussion_test::Sine(sample, inputHz, sampleRate);
    const float shifted = shifter.Process(input);
    if (sample >= Warmup)
      output[sample - Warmup] = shifted;
  }

  const double wanted = ToneMagnitude(output, Analysis, inputHz + shiftHz,
                                      sampleRate);
  const double image = ToneMagnitude(output, Analysis, inputHz - shiftHz,
                                     sampleRate);
  CheckNear(wanted, 1.0, .006,
            "SSB frequency shifter preserves wanted tone level");
  Check(image < .002, "SSB frequency shifter rejects the opposite sideband");
}

void TestSignedShiftsAtSupportedRates() {
  for (const float sampleRate : {44100.f, 48000.f, 88200.f, 96000.f, 192000.f}) {
    CheckShift(sampleRate, .12f * sampleRate, .017f * sampleRate);
    CheckShift(sampleRate, .12f * sampleRate, -.019f * sampleRate);
  }
}

void TestZeroShiftIsExactDelay() {
  tfdsp::percussion::FrequencyShifter shifter;
  shifter.Prepare(48000.f);
  float history[tfdsp::percussion::FrequencyShifter::LatencySamples + 1]{};
  std::size_t index = 0;
  double largestError = 0.0;
  for (std::size_t sample = 0; sample < 20000; ++sample) {
    const float input = percussion_test::Sine(sample, 7331.f, 48000.f);
    history[index] = input;
    const std::size_t delayed =
        (index + 1) % std::size(history);
    const float output = shifter.Process(input);
    largestError = std::max(largestError,
                            std::abs(static_cast<double>(output - history[delayed])));
    index = delayed;
  }
  CheckNear(largestError, 0.0, 1.e-7,
            "zero frequency shift is an exact linear-phase delay");
}

void TestPhaseContinuousAutomation() {
  tfdsp::percussion::FrequencyShifter shifter;
  shifter.Prepare(48000.f);
  float previous = 0.f;
  float largestStep = 0.f;
  bool finite = true;
  for (std::size_t sample = 0; sample < 96000; ++sample) {
    const float shift = -900.f + 1800.f * static_cast<float>(sample) / 95999.f;
    shifter.SetShiftHz(shift);
    const float input = percussion_test::Sine(sample, 3100.f, 48000.f);
    const float output = shifter.Process(input);
    if (sample > 4096)
      largestStep = std::max(largestStep, std::abs(output - previous));
    previous = output;
    finite = finite && std::isfinite(output);
  }
  Check(finite, "frequency shifter remains finite through a signed sweep");
  Check(largestStep < .65f,
        "frequency shifter crosses zero without a discontinuity");

  shifter.SetShiftHz(std::numeric_limits<float>::quiet_NaN());
  Check(std::isfinite(shifter.Process(std::numeric_limits<float>::infinity())),
        "frequency shifter sanitizes non-finite controls and input");
}

void TestOscillatorNorm() {
  tfdsp::percussion::QuadratureOscillator oscillator;
  oscillator.Prepare(192000.f);
  oscillator.SetFrequencyHz(-37123.f);
  double largestError = 0.0;
  for (std::size_t sample = 0; sample < 1000000; ++sample) {
    const auto value = oscillator.Process();
    largestError = std::max(largestError,
        std::abs(value.cosine * value.cosine + value.sine * value.sine - 1.0));
  }
  CheckNear(largestError, 0.0, 8.e-6,
            "recursive quadrature oscillator retains unit magnitude");
}

} // namespace

int main() {
  TestSignedShiftsAtSupportedRates();
  TestZeroShiftIsExactDelay();
  TestPhaseContinuousAutomation();
  TestOscillatorNorm();
  if (percussion_test::failures == 0)
    std::cout << "All percussion frequency-shift tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
