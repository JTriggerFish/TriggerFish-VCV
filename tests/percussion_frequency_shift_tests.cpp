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

  shifter.Reset();
  bool exactSilence = true;
  for (std::size_t sample = 0; sample < 1024; ++sample)
    exactSilence = exactSilence &&
        shifter.Process(sample == 0 ? std::numeric_limits<float>::denorm_min()
                                    : 0.f) == 0.f;
  Check(exactSilence,
        "frequency-shifter filters and Hilbert history flush subnormals");
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

double ShiftedToneLevel(const float inputHz, const float shiftHz,
                        const float measuredHz) {
  constexpr std::size_t Warmup = 16384;
  constexpr std::size_t Count = 65536;
  static float output[Count];
  tfdsp::percussion::FrequencyShifter shifter;
  shifter.Prepare(48000.f);
  shifter.SetShiftHz(shiftHz);
  for (std::size_t sample = 0; sample < Warmup + Count; ++sample) {
    const float value = shifter.Process(
        percussion_test::Sine(sample, inputHz, 48000.f));
    if (sample >= Warmup)
      output[sample - Warmup] = value;
  }
  return ToneMagnitude(output, Count, measuredHz, 48000.f);
}

void TestTranslationBoundariesAreRejected() {
  const double upperFold = ShiftedToneLevel(23000.f, 4000.f, 21000.f);
  const double lowerFold = ShiftedToneLevel(1000.f, -3000.f, 2000.f);
  Check(upperFold < .12,
        "source limiting rejects content translated beyond Nyquist");
  Check(lowerFold < .12,
        "source limiting rejects content translated below DC");
}

void TestWideAutomationRemainsBounded() {
  tfdsp::percussion::FrequencyShifter shifter;
  shifter.Prepare(48000.f);
  float peak = 0.f;
  bool finite = true;
  for (std::size_t sample = 0; sample < 192000; ++sample) {
    const float phase = static_cast<float>(sample % 48000) / 47999.f;
    const float triangle = phase < .5f ? 4.f * phase - 1.f
                                      : 3.f - 4.f * phase;
    shifter.SetShiftHz(8000.f * triangle);
    const float output = shifter.Process(
        percussion_test::Sine(sample, 10000.f, 48000.f));
    if (sample > 4096)
      peak = std::max(peak, std::abs(output));
    finite = finite && std::isfinite(output);
  }
  Check(finite, "frequency-shift filter interpolation remains finite");
  Check(peak < 2.f,
        "wide frequency-shift automation has no coefficient-update spikes");
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
  TestTranslationBoundariesAreRejected();
  TestWideAutomationRemainsBounded();
  TestOscillatorNorm();
  if (percussion_test::failures == 0)
    std::cout << "All percussion frequency-shift tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
