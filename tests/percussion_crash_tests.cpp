#include "percussion_test_support.hpp"

#include "tfdsp/percussion/crash_cymbal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using percussion_test::Check;

namespace {

std::vector<float> Render(const float strength, const float location,
                          const float hardness, const std::uint32_t seed,
                          const float seconds = 2.f,
                          const float sampleRate = 48000.f) {
  using namespace tfdsp::percussion;
  CrashCymbal cymbal;
  cymbal.Prepare(sampleRate, DefaultCrashCymbalParameters(sampleRate));
  cymbal.Trigger({strength, location, hardness, seed});
  std::vector<float> output(
      static_cast<std::size_t>(seconds * sampleRate));
  for (float &sample : output)
    sample = cymbal.Process();
  return output;
}

double Energy(const std::vector<float> &audio) {
  double energy = 0.0;
  for (const float sample : audio)
    energy += static_cast<double>(sample) * sample;
  return energy;
}

double Difference(const std::vector<float> &first,
                  const std::vector<float> &second) {
  double energy = 0.0;
  for (std::size_t sample = 0; sample < first.size(); ++sample) {
    const double delta = first[sample] - second[sample];
    energy += delta * delta;
  }
  return energy;
}

void TestDeterministicAndResponsive() {
  const auto first = Render(.8f, 1.f, .65f, 1234);
  const auto repeated = Render(.8f, 1.f, .65f, 1234);
  const auto differentSeed = Render(.8f, 1.f, .65f, 1235);
  Check(first == repeated, "crash rendering is deterministic for one seed");
  Check(Difference(first, differentSeed) > 1.e-7,
        "crash contact variation responds to its seed");

  const auto bell = Render(.8f, 0.f, .65f, 1234);
  const auto bow = Render(.8f, .5f, .65f, 1234);
  Check(Difference(bell, bow) > .01 * Energy(bow),
        "crash location changes the body projection audibly");

  const auto soft = Render(.8f, 1.f, .15f, 1234);
  const auto hard = Render(.8f, 1.f, .95f, 1234);
  Check(Difference(soft, hard) > .01 * Energy(hard),
        "crash hardness changes the contact audibly");
}

void TestVelocityEnergy() {
  const std::array<float, 4> strengths{.2f, .45f, .7f, 1.f};
  double previous = 0.0;
  for (const float strength : strengths) {
    const double current = Energy(Render(strength, 1.f, .65f, 81));
    if (!(current > previous))
      std::cerr << "crash velocity energy " << strength << ": " << current
                << " after " << previous << '\n';
    Check(current > previous,
          "crash energy increases across the velocity sweep");
    previous = current;
  }
}

void TestSparseModesArePlacedDirectly() {
  using namespace tfdsp::percussion;
  CrashCymbalFitParameters fit;
  fit.sparseFrequencyHz[0] = 731.f;
  fit.sparseFrequencyHz[1] = 1193.f;
  fit.bodyDecaySeconds.fill(4.25f);
  fit.sparseDecayRatio[0] = 1.f;
  const auto parameters = DefaultCrashCymbalParameters(48000.f, fit);
  Check(std::abs(parameters.sparseModes[0].frequencyHz - 731.f) < 1.e-5f &&
            std::abs(parameters.sparseModes[1].frequencyHz - 1193.f) < 1.e-5f,
        "crash sparse frequencies are independently and directly placed");
  Check(std::abs(parameters.sparseModes[0].decaySeconds - 4.25f) < 1.e-5f,
        "crash sparse decay follows the shared body T60 curve");
}

void TestImplementFamiliesAreDistinct() {
  using namespace tfdsp::percussion;
  constexpr float sampleRate = 48000.f;
  const auto render = [](const float implement, const float spread = .2f) {
    CrashCymbal cymbal;
    cymbal.Prepare(sampleRate, DefaultCrashCymbalParameters(sampleRate));
    cymbal.Trigger({.8f, .75f, .65f, 91, implement, spread});
    std::vector<float> result(12000);
    for (float &sample : result) sample = cymbal.Process();
    return result;
  };
  const auto brush = render(0.f);
  const auto mallet = render(.5f);
  const auto stick = render(1.f);
  const auto brushSweep = render(0.f, 1.f);
  Check(Difference(brush, mallet) > 1.e-5 * Energy(mallet),
        "brush and mallet contacts are distinct");
  Check(Difference(mallet, stick) > 1.e-5 * Energy(stick),
        "mallet and stick contacts are distinct");
  Check(Difference(brush, brushSweep) > .01 * Energy(brushSweep),
        "brush contact spread changes a tap into a sustained gesture");
}

void TestDefaultBodyCoversTheMeasuredLowRegion() {
  using namespace tfdsp::percussion;
  const auto parameters = DefaultCrashCymbalParameters(48000.f);
  const auto lowSparseModes = std::count_if(
      parameters.sparseModes.begin(), parameters.sparseModes.end(),
      [](const auto &mode) { return mode.frequencyHz < 500.f; });
  Check(lowSparseModes >= 1,
        "default crash retains a resolved plate ridge below 500 Hz");
  Check(parameters.denseModes.front().frequencyHz < 200.f,
        "default crash dense wash extends below the old 700 Hz gap");
  Check(parameters.observation[0].radiation.lowCutHz <= 50.f &&
            parameters.observation[1].radiation.lowCutHz <= 50.f &&
            parameters.observation[2].radiation.lowCutHz <= 50.f,
        "default crash observation does not remove its low plate body");
}

void TestContactCalibrationMacrosAreAudible() {
  using namespace tfdsp::percussion;
  constexpr float sampleRate = 48000.f;
  CrashCymbalFitParameters firstFit;
  CrashCymbalFitParameters secondFit;
  secondFit.contactNoiseGain = 2.f;
  secondFit.contactChirpFrequencyScale = 1.5f;
  CrashCymbal first;
  CrashCymbal second;
  first.Prepare(sampleRate, DefaultCrashCymbalParameters(sampleRate, firstFit));
  second.Prepare(sampleRate, DefaultCrashCymbalParameters(sampleRate, secondFit));
  first.Trigger({.8f, 1.f, .65f, 1234});
  second.Trigger({.8f, 1.f, .65f, 1234});
  double difference = 0.0;
  for (int sample = 0; sample < 720; ++sample) {
    const double delta = first.Process() - second.Process();
    difference += delta * delta;
  }
  Check(difference > 1.e-8,
        "crash contact calibration macros change the initial transient");
}

void TestMuteIsPassive() {
  using namespace tfdsp::percussion;
  constexpr float sampleRate = 48000.f;
  CrashCymbal natural;
  CrashCymbal muted;
  const auto parameters = DefaultCrashCymbalParameters(sampleRate);
  natural.Prepare(sampleRate, parameters);
  muted.Prepare(sampleRate, parameters);
  natural.Trigger({.9f, 1.f, .65f, 9});
  muted.Trigger({.9f, 1.f, .65f, 9});
  double naturalTail = 0.0;
  double mutedTail = 0.0;
  for (std::size_t sample = 0; sample < 2 * 48000; ++sample) {
    if (sample == 4800)
      muted.SetMute(1.f);
    const float first = natural.Process();
    const float second = muted.Process();
    if (sample > 12000) {
      naturalTail += static_cast<double>(first) * first;
      mutedTail += static_cast<double>(second) * second;
    }
  }
  Check(mutedTail < .25 * naturalTail,
        "crash mute removes stored energy instead of sustaining it");

  CrashCymbal silent;
  silent.Prepare(sampleRate, parameters);
  silent.SetMute(1.f);
  double silentEnergy = 0.0;
  for (std::size_t sample = 0; sample < 48000; ++sample) {
    const float output = silent.Process();
    silentEnergy += static_cast<double>(output) * output;
  }
  Check(silentEnergy == 0.0, "changing crash mute does not inject energy");
}

void TestFiniteAtSupportedRates() {
  for (const float sampleRate : {44100.f, 48000.f, 96000.f, 192000.f}) {
    const auto audio = Render(1.f, .73f, 1.f, 0xffffffffu, .5f, sampleRate);
    Check(std::all_of(audio.begin(), audio.end(), [](const float sample) {
      return std::isfinite(sample) && std::abs(sample) < 100.f;
    }), "crash remains finite and bounded across sample rates");
  }
}

void TestAnalysisFrameMatchesOutput() {
  using namespace tfdsp::percussion;
  CrashCymbal framed;
  CrashCymbal plain;
  const auto parameters = DefaultCrashCymbalParameters(48000.f);
  framed.Prepare(48000.f, parameters);
  plain.Prepare(48000.f, parameters);
  framed.Trigger({.8f, 1.f, .65f, 71});
  plain.Trigger({.8f, 1.f, .65f, 71});
  for (int sample = 0; sample < 4096; ++sample) {
    const auto frame = framed.ProcessFrame();
    Check(frame.output == plain.Process(),
          "crash analysis taps do not change production output");
    Check(std::isfinite(frame.directContact) &&
              std::isfinite(frame.dispersion) &&
              std::isfinite(frame.sparseModes) &&
              std::isfinite(frame.denseResidual),
          "crash analysis taps remain finite");
  }
}

void TestBodyBranchesCanBeAblatedIndependently() {
  using namespace tfdsp::percussion;
  CrashCymbalFitParameters fit;
  fit.sparseGain = 0.f;
  CrashCymbal noSparse;
  noSparse.Prepare(48000.f, DefaultCrashCymbalParameters(48000.f, fit));
  noSparse.Trigger({.8f, 1.f, .65f, 73});
  bool sparseSilent = true;
  for (int sample = 0; sample < 4096; ++sample)
    sparseSilent = sparseSilent && noSparse.ProcessFrame().sparseModes == 0.f;
  Check(sparseSilent, "zero sparse presentation gain ablates sparse DSP state");

  fit.sparseGain = .35f;
  fit.denseGain = 0.f;
  CrashCymbal noDense;
  noDense.Prepare(48000.f, DefaultCrashCymbalParameters(48000.f, fit));
  noDense.Trigger({.8f, 1.f, .65f, 73});
  bool denseSilent = true;
  for (int sample = 0; sample < 4096; ++sample)
    denseSilent = denseSilent && noDense.ProcessFrame().denseResidual == 0.f;
  Check(denseSilent, "zero dense presentation gain ablates dense DSP state");
}

void TestMaximumNonlinearBloomRemainsBounded() {
  using namespace tfdsp::percussion;
  CrashCymbalFitParameters fit;
  fit.dispersionFeedback = .995f;
  fit.dispersionDrive = 8.f;
  fit.dispersionExcursionSamples = 16.f;
  CrashCymbal cymbal;
  cymbal.Prepare(48000.f, DefaultCrashCymbalParameters(48000.f, fit));
  cymbal.Trigger({1.f, 1.f, 1.f, 91});
  for (int sample = 0; sample < 96000; ++sample) {
    const float output = cymbal.Process();
    Check(std::isfinite(output) && std::abs(output) < 100.f,
          "maximum crash self-phase settings remain bounded");
  }
}

void TestDenseVelocityLossIsPassive() {
  using namespace tfdsp::percussion;
  CrashCymbalFitParameters naturalFit;
  naturalFit.directGain = 0.f;
  naturalFit.sparseGain = 0.f;
  naturalFit.denseGain = 1.f;
  auto dampedFit = naturalFit;
  dampedFit.denseVelocityLossNepersPerSecond = 3.f;
  CrashCymbal natural;
  CrashCymbal damped;
  natural.Prepare(48000.f, DefaultCrashCymbalParameters(48000.f, naturalFit));
  damped.Prepare(48000.f, DefaultCrashCymbalParameters(48000.f, dampedFit));
  natural.Trigger({1.f, 1.f, .65f, 93});
  damped.Trigger({1.f, 1.f, .65f, 93});
  double naturalTail = 0.0;
  double dampedTail = 0.0;
  for (int sample = 0; sample < 96000; ++sample) {
    const float first = natural.ProcessFrame().denseResidual;
    const float second = damped.ProcessFrame().denseResidual;
    if (sample > 12000) {
      naturalTail += static_cast<double>(first) * first;
      dampedTail += static_cast<double>(second) * second;
    }
  }
  Check(dampedTail < naturalTail,
        "velocity-dependent dense loss only removes residual energy");
}

void TestRepeatedHitsAccumulateBodyEnergy() {
  using namespace tfdsp::percussion;
  constexpr int interval = 6000;
  constexpr int hitCount = 16;
  CrashCymbalFitParameters fit;
  fit.dispersionDrive = 6.f;
  fit.dispersionExcursionSamples = 12.f;
  fit.dispersionFeedback = .998f;
  fit.dispersionLowDecaySeconds = 2.f;
  fit.dispersionMiddleDecaySeconds = 1.5f;
  fit.dispersionHighDecaySeconds = .8f;
  CrashCymbal cymbal;
  cymbal.Prepare(48000.f, DefaultCrashCymbalParameters(48000.f, fit));
  std::array<double, hitCount> intervalEnergy{};
  for (int sample = 0; sample < hitCount * interval; ++sample) {
    const int hit = sample / interval;
    if (sample % interval == 0)
      cymbal.Trigger({.45f, 1.f, .65f, static_cast<std::uint32_t>(100 + hit)});
    const auto frame = cymbal.ProcessFrame();
    intervalEnergy[hit] +=
        static_cast<double>(frame.sparseModes) * frame.sparseModes +
        static_cast<double>(frame.denseResidual) * frame.denseResidual;
  }
  const double initial = intervalEnergy[0];
  const double earlyPeak = *std::max_element(
      intervalEnergy.begin() + 1, intervalEnergy.begin() + 5);
  const double lateMean = .25 * (
      intervalEnergy[hitCount - 4] + intervalEnergy[hitCount - 3] +
      intervalEnergy[hitCount - 2] + intervalEnergy[hitCount - 1]);
  if (!(earlyPeak > 1.5 * initial && lateMean > .5 * initial))
  {
    std::cerr << "crash repeated-hit initial/peak/late energy: " << initial
              << "/" << earlyPeak << "/" << lateMean << " intervals:";
    for (const double energy : intervalEnergy)
      std::cerr << ' ' << energy;
    std::cerr << '\n';
  }
  Check(earlyPeak > 1.5 * initial && lateMean > .5 * initial,
        "repeated crash hits build an initial swell and retain late energy");
}

} // namespace

int main() {
  TestDeterministicAndResponsive();
  TestVelocityEnergy();
  TestSparseModesArePlacedDirectly();
  TestImplementFamiliesAreDistinct();
  TestDefaultBodyCoversTheMeasuredLowRegion();
  TestContactCalibrationMacrosAreAudible();
  TestMuteIsPassive();
  TestFiniteAtSupportedRates();
  TestAnalysisFrameMatchesOutput();
  TestBodyBranchesCanBeAblatedIndependently();
  TestMaximumNonlinearBloomRemainsBounded();
  TestDenseVelocityLossIsPassive();
  TestRepeatedHitsAccumulateBodyEnergy();
  if (percussion_test::failures == 0)
    std::cout << "All percussion crash tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
