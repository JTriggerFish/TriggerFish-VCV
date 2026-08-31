#include "percussion_test_support.hpp"

#include "tfdsp/percussion/dispersion_loop.hpp"
#include "tfdsp/percussion/modulated_fractional_delay.hpp"
#include "tfdsp/percussion/self_phase_delay.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

tfdsp::percussion::DispersionLoopParameters LinearLoopParameters();

template <typename Stage>
double ToneGain(Stage &stage, const float frequency, const float sampleRate) {
  constexpr std::size_t Warmup = 4096;
  constexpr std::size_t Count = 32768;
  double real = 0.0;
  double imaginary = 0.0;
  for (std::size_t sample = 0; sample < Warmup + Count; ++sample) {
    const float input = percussion_test::Sine(sample, frequency, sampleRate);
    const float output = stage.Process(input);
    if (sample < Warmup)
      continue;
    const double phase = 6.283185307179586 * frequency * sample / sampleRate;
    real += output * std::cos(phase);
    imaginary -= output * std::sin(phase);
  }
  return 2.0 * std::hypot(real, imaginary) / Count;
}

void TestZeroDriveSelfPhaseReference() {
  tfdsp::percussion::SelfPhaseDelayParameters parameters;
  parameters.centreDelaySamples = 11.375f;
  parameters.maximumExcursionSamples = 3.f;
  parameters.drive = 0.f;
  for (const float sampleRate : {44100.f, 48000.f, 88200.f, 96000.f, 192000.f}) {
    tfdsp::percussion::SelfPhaseDelay production;
    tfdsp::percussion::SelfPhaseDelayReference4x reference;
    production.Prepare(sampleRate, 64.f);
    reference.Prepare(sampleRate, 64.f);
    production.SetParameters(parameters);
    reference.SetParameters(parameters);
    for (const float frequency : {.04f * sampleRate, .12f * sampleRate,
                                  .24f * sampleRate, .38f * sampleRate}) {
      const double actual = ToneGain(production, frequency, sampleRate);
      const double expected = ToneGain(reference, frequency, sampleRate);
      CheckNear(actual, expected, .012,
                "2x zero-drive response follows the 4x reference");
      production.Reset();
      reference.Reset();
    }
  }
}

void TestSelfPhaseReportsItsEffectiveDelay() {
  tfdsp::percussion::SelfPhaseDelay stage;
  stage.Prepare(48000.f, 64.f);
  tfdsp::percussion::SelfPhaseDelayParameters parameters;
  parameters.centreDelaySamples = 200.f;
  stage.SetParameters(parameters);
  CheckNear(stage.CentreDelaySamples(), 64.0, 1.e-6,
            "self-phase reports its clamped maximum delay");
  parameters.centreDelaySamples =
      std::numeric_limits<float>::quiet_NaN();
  stage.SetParameters(parameters);
  CheckNear(stage.CentreDelaySamples(), 12.0, 1.e-6,
            "self-phase reports its host-rate fallback delay");
  stage.Reset();
  bool exactSilence = true;
  for (std::size_t sample = 0; sample < 256; ++sample)
    exactSilence = exactSilence &&
        stage.Process(sample == 0 ? std::numeric_limits<float>::denorm_min()
                                  : 0.f) == 0.f;
  Check(exactSilence,
        "self-phase filter, envelope and delay flush subnormals");

  auto loopParameters = LinearLoopParameters();
  loopParameters.selfPhase.centreDelaySamples = 200.f;
  tfdsp::percussion::DispersionLoop loop;
  loop.Prepare(48000.f, 64.f, loopParameters);
  CheckNear(loop.MinimumPropagationSamples(), 90.0, 1.e-6,
            "dispersion loss uses the effective self-phase path length");
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
  auto noFeedbackParameters = parameters;
  noFeedbackParameters.feedbackGain = 0.f;
  tfdsp::percussion::DispersionLoop loop;
  tfdsp::percussion::DispersionLoop noFeedback;
  loop.Prepare(48000.f, 64.f, parameters);
  noFeedback.Prepare(48000.f, 64.f, noFeedbackParameters);
  const auto minimum = static_cast<std::size_t>(
      std::lround(loop.MinimumPropagationSamples()));
  const auto nominal = static_cast<std::size_t>(
      std::lround(loop.NominalPropagationSamples()));
  std::size_t firstNonzero = 256;
  std::size_t feedbackNonzero = 256;
  std::size_t firstPeak = 0;
  float peak = 0.f;
  for (std::size_t sample = 0; sample < 256; ++sample) {
    const float input = sample == 0 ? 1.f : 0.f;
    const float output = loop.Process(input);
    const float direct = noFeedback.Process(input);
    if (std::abs(direct) > 1.e-9f && firstNonzero == 256)
      firstNonzero = sample;
    if (std::abs(output - direct) > 1.e-9f && feedbackNonzero == 256)
      feedbackNonzero = sample;
    if (std::abs(direct) > peak) {
      peak = std::abs(direct);
      firstPeak = sample;
    }
  }
  Check(firstNonzero == minimum,
        "dispersion output begins at its declared causal propagation delay");
  Check(feedbackNonzero == 2 * minimum,
        "dispersion feedback sum adds no hidden recurrence sample");
  Check(std::abs(static_cast<long long>(firstPeak) -
                 static_cast<long long>(nominal)) <= 2,
        "dispersion declares the oversampled stage's nominal latency");
}

template <typename Stage>
std::vector<float> RenderNonlinear(Stage &stage, const float sampleRate) {
  constexpr std::size_t Warmup = 4096;
  constexpr std::size_t Count = 32768;
  std::vector<float> output(Count);
  for (std::size_t sample = 0; sample < Warmup + Count; ++sample) {
    const float input = .55f * percussion_test::Sine(sample, 6717.f, sampleRate) +
                        .35f * percussion_test::Sine(sample, 11231.f, sampleRate);
    const float value = stage.Process(input);
    if (sample >= Warmup)
      output[sample - Warmup] = value;
  }
  return output;
}

std::vector<double> SparseSpectrum(const std::vector<float> &signal) {
  std::vector<double> spectrum;
  for (std::size_t bin = 128; bin < signal.size() / 2; bin += 128) {
    double real = 0.0;
    double imaginary = 0.0;
    for (std::size_t sample = 0; sample < signal.size(); ++sample) {
      const double window = .5 - .5 * std::cos(
          6.283185307179586 * sample / (signal.size() - 1));
      const double phase = 6.283185307179586 * bin * sample / signal.size();
      real += window * signal[sample] * std::cos(phase);
      imaginary -= window * signal[sample] * std::sin(phase);
    }
    spectrum.push_back(std::hypot(real, imaginary));
  }
  return spectrum;
}

double LevelMatchedSpectrumError(const std::vector<double> &actual,
                                 const std::vector<double> &reference) {
  double cross = 0.0;
  double referenceEnergy = 0.0;
  double actualEnergy = 0.0;
  for (std::size_t index = 0; index < actual.size(); ++index) {
    cross += actual[index] * reference[index];
    referenceEnergy += reference[index] * reference[index];
    actualEnergy += actual[index] * actual[index];
  }
  const double gain = cross / std::max(1.e-20, actualEnergy);
  double error = 0.0;
  for (std::size_t index = 0; index < actual.size(); ++index) {
    const double residual = gain * actual[index] - reference[index];
    error += residual * residual;
  }
  return std::sqrt(error / std::max(1.e-20, referenceEnergy));
}

void TestTwoTimesNonlinearityAgainstFourTimesReference() {
  constexpr float SampleRate = 48000.f;
  tfdsp::percussion::SelfPhaseDelayParameters parameters;
  parameters.centreDelaySamples = 14.f;
  parameters.maximumExcursionSamples = 4.f;
  parameters.drive = 14.f;
  parameters.toneHz = 14000.f;
  parameters.normalization = .65f;
  tfdsp::percussion::SelfPhaseDelayReference1x oneTimes;
  tfdsp::percussion::SelfPhaseDelay twoTimes;
  tfdsp::percussion::SelfPhaseDelayReference4x fourTimes;
  oneTimes.Prepare(SampleRate, 64.f);
  twoTimes.Prepare(SampleRate, 64.f);
  fourTimes.Prepare(SampleRate, 64.f);
  oneTimes.SetParameters(parameters);
  twoTimes.SetParameters(parameters);
  fourTimes.SetParameters(parameters);
  const auto reference = SparseSpectrum(RenderNonlinear(fourTimes, SampleRate));
  const double oneTimesError = LevelMatchedSpectrumError(
      SparseSpectrum(RenderNonlinear(oneTimes, SampleRate)), reference);
  const double twoTimesError = LevelMatchedSpectrumError(
      SparseSpectrum(RenderNonlinear(twoTimes, SampleRate)), reference);
  if (!(twoTimesError < oneTimesError && twoTimesError < .25))
    std::cerr << "self-phase spectral error 1x/2x: " << oneTimesError << "/"
              << twoTimesError << '\n';
  Check(twoTimesError < oneTimesError,
        "2x self-phase spectrum is closer to 4x than the host-rate core");
  Check(twoTimesError < .25,
        "2x self-phase spectrum remains close to the 4x offline reference");
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
  double earlyTailEnergy = 0.0;
  double lateTailEnergy = 0.0;
  for (std::size_t sample = 0; sample < 1000000; ++sample) {
    const float drive = sample < 20000 ?
        percussion_test::Sine(sample, 17000.f, 192000.f) : 0.f;
    const float output = loop.Process(drive);
    finite = finite && std::isfinite(output);
    peak = std::max(peak, std::abs(output));
    if (sample >= 40000 && sample < 140000)
      earlyTailEnergy += output * output;
    if (sample >= 900000)
      lateTailEnergy += output * output;
  }
  Check(finite, "nonlinear dispersion remains finite under long high feedback");
  Check(peak < 32.f, "nonlinear dispersion retains explicit internal headroom");
  Check(lateTailEnergy < .01 * earlyTailEnergy,
        "nonlinear dispersion contracts after excitation without limiting");
  Check(std::isfinite(loop.Process(std::numeric_limits<float>::infinity())),
        "dispersion loop sanitizes non-finite drive");
}

} // namespace

int main() {
  TestZeroDriveSelfPhaseReference();
  TestSelfPhaseReportsItsEffectiveDelay();
  TestExplicitFeedbackCausality();
  TestTwoTimesNonlinearityAgainstFourTimesReference();
  TestNonlinearStress();
  if (percussion_test::failures == 0)
    std::cout << "All percussion dispersion tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
