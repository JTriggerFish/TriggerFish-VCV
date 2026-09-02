#include "percussion_test_support.hpp"

#include "tfdsp/percussion/crash_cymbal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <utility>
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

double NormalizedDifferenceEnergy(const std::vector<float> &audio) {
  double signal = 0.0;
  double difference = 0.0;
  float previous = 0.f;
  for (const float sample : audio) {
    signal += static_cast<double>(sample) * sample;
    const double delta = sample - previous;
    difference += delta * delta;
    previous = sample;
  }
  return difference / std::max(signal, 1.e-30);
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

void TestVelocityChangesCymbalRegime() {
  using namespace tfdsp::percussion;
  constexpr float sampleRate = 48000.f;
  const auto quiet = Render(.25f, 1.f, .65f, 81, .75f);
  const auto loud = Render(1.f, 1.f, .65f, 81, .75f);
  const double quietBrightness = NormalizedDifferenceEnergy(quiet);
  const double loudBrightness = NormalizedDifferenceEnergy(loud);
  if (!(loudBrightness > 1.05 * quietBrightness)) {
    std::cerr << "crash velocity brightness quiet/loud: "
              << quietBrightness << '/' << loudBrightness << '\n';
  }
  Check(loudBrightness > 1.05 * quietBrightness,
        "strong crash strikes excite a brighter spectrum");

  struct BranchEnergy { double direct{}; double bloom{}; double dense{}; };
  const auto measure = [](const float strength) {
    CrashCymbal cymbal;
    cymbal.Prepare(sampleRate, DefaultCrashCymbalParameters(sampleRate));
    cymbal.Trigger({strength, 1.f, .65f, 81});
    BranchEnergy result;
    for (int sample = 0; sample < 24000; ++sample) {
      const auto frame = cymbal.ProcessFrame();
      result.direct += static_cast<double>(frame.directContact) *
          frame.directContact;
      result.bloom += static_cast<double>(frame.dispersion) * frame.dispersion;
      result.dense += static_cast<double>(frame.denseResidual) *
          frame.denseResidual;
    }
    return result;
  };
  const auto quietBranches = measure(.25f);
  const auto loudBranches = measure(1.f);
  const double directGrowth = loudBranches.direct / quietBranches.direct;
  const double bloomGrowth = loudBranches.bloom / quietBranches.bloom;
  const double denseGrowth = loudBranches.dense / quietBranches.dense;
  if (!(denseGrowth > 2.0 * bloomGrowth && bloomGrowth > 1.0)) {
    std::cerr << "crash velocity direct/bloom/dense growth: "
              << directGrowth << '/' << bloomGrowth << '/' << denseGrowth
              << '\n';
  }
  Check(denseGrowth > 2.0 * bloomGrowth && bloomGrowth > 1.0,
        "strong strikes increase bloom and brighten the modal wash");
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

void TestBodyDecayCurveUsesActiveErbKnots() {
  using namespace tfdsp::percussion;
  CrashCymbalFitParameters fit;
  fit.bodyDecayActive.fill(false);
  fit.bodyDecayActive.front() = true;
  fit.bodyDecayActive.back() = true;
  fit.bodyDecaySeconds.fill(1.f);
  fit.bodyDecayFrequencyHz[1] = 1000.f;
  fit.bodyDecaySeconds[1] = 8.f;
  fit.bodyDecayActive[1] = true;
  fit.sparseFrequencyHz[0] = 1000.f;
  fit.sparseDecayRatio[0] = 1.f;
  const auto withKnot = DefaultCrashCymbalParameters(48000.f, fit);
  Check(std::abs(withKnot.sparseModes[0].decaySeconds - 8.f) < 1.e-4f,
        "active body T60 knots are sampled by modal preparation");

  fit.bodyDecayActive[1] = false;
  const auto withoutKnot = DefaultCrashCymbalParameters(48000.f, fit);
  Check(std::abs(withoutKnot.sparseModes[0].decaySeconds - 1.f) < 1.e-4f,
        "inactive body T60 knots do not affect modal decay");
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
  const double brushToStickEnergy = Energy(brush) /
      std::max(Energy(stick), 1.e-30);
  Check(Difference(brush, mallet) > 1.e-5 * Energy(mallet),
        "brush and mallet contacts are distinct");
  Check(Difference(mallet, stick) > 1.e-5 * Energy(stick),
        "mallet and stick contacts are distinct");
  Check(Difference(brush, brushSweep) > .01 * Energy(brushSweep),
        "brush contact spread changes a tap into a sustained gesture");
  // The generic core preset is deliberately not implement-level matched; the
  // workbench preset owns that perceptual calibration and tests it separately.
  if (!(brushToStickEnergy > .01 && brushToStickEnergy < 4.0))
    std::cerr << "crash brush/stick energy ratio: "
              << brushToStickEnergy << '\n';
  Check(brushToStickEnergy > .01 && brushToStickEnergy < 4.0,
        "generic brush output remains audible without level matching");

  const auto contactShape = [](const float implement) {
    CrashCymbal cymbal;
    cymbal.Prepare(sampleRate, DefaultCrashCymbalParameters(sampleRate));
    cymbal.Trigger({.8f, .75f, .65f, 91, implement, .2f});
    std::array<double, 2> energy{};
    for (int sample = 0; sample < 12000; ++sample) {
      const double contact = cymbal.ProcessFrame().directContact;
      energy[sample < 96 ? 0 : 1] += contact * contact;
    }
    return energy;
  };
  const auto brushContact = contactShape(0.f);
  const auto stickContact = contactShape(1.f);
  const double brushTailRatio = brushContact[1] /
      std::max(brushContact[0], 1.e-30);
  const double stickTailRatio = stickContact[1] /
      std::max(stickContact[0], 1.e-30);
  Check(brushTailRatio > 4.0 * stickTailRatio,
        "brush contact is distributed instead of becoming a stick transient");
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

void TestDenseModeCountUsesNestedExtensionBank() {
  using namespace tfdsp::percussion;
  const auto active = [](const auto &modes) {
    return std::count_if(modes.begin(), modes.end(), [](const auto &mode) {
      return mode.inputGain != 0.f && mode.outputGain != 0.f;
    });
  };
  CrashCymbalFitParameters fit;
  fit.denseModeDensity = .5f;
  auto parameters = DefaultCrashCymbalParameters(48000.f, fit);
  Check(active(parameters.denseModes) == 1024 &&
            active(parameters.denseExtensionModes) == 0,
        "half-density crash uses exactly 1024 primary wash modes");
  fit.denseModeDensity = 1.f;
  parameters = DefaultCrashCymbalParameters(48000.f, fit);
  Check(active(parameters.denseModes) == 2048 &&
            active(parameters.denseExtensionModes) == 0,
        "factory crash retains its 2048-mode primary wash");
  fit.denseModeDensity = 1.5f;
  parameters = DefaultCrashCymbalParameters(48000.f, fit);
  Check(active(parameters.denseModes) == 2048 &&
            active(parameters.denseExtensionModes) == 1024,
        "higher crash density adds modes in the extension bank");
  fit.denseModeDensity = 2.f;
  parameters = DefaultCrashCymbalParameters(48000.f, fit);
  Check(active(parameters.denseModes) == 2048 &&
            active(parameters.denseExtensionModes) == 2048,
        "maximum crash density exposes exactly 4096 wash modes");
}

void TestDenseColourShapesExcitationEnergy() {
  using namespace tfdsp::percussion;
  const auto parameters = DefaultCrashCymbalParameters(48000.f);
  double squaredInputGain = 0.0;
  bool unityObservation = true;
  for (const auto &mode : parameters.denseModes) {
    squaredInputGain += static_cast<double>(mode.inputGain) * mode.inputGain;
    unityObservation = unityObservation &&
        (mode.outputGain == 0.f || mode.outputGain == 1.f);
  }
  Check(std::abs(squaredInputGain - .0025) < 1.e-6,
        "dense colour is a normalized modal excitation-energy curve");
  Check(unityObservation,
        "dense modal colour is not hidden in per-mode observation gains");
}

void TestUnifiedFieldExpandsAnchorsWithoutChangingDriveEnergy() {
  using namespace tfdsp::percussion;
  const auto measure = [](const CrashModalField::Parameters &modes) {
    std::pair<std::size_t, double> result{};
    for (const auto &mode : modes) {
      if (mode.inputGain == 0.f || mode.outputGain == 0.f)
        continue;
      ++result.first;
      result.second += static_cast<double>(mode.inputGain) * mode.inputGain;
    }
    return result;
  };
  CrashCymbalFitParameters coherentFit;
  coherentFit.fieldTurbulence = 0.f;
  auto diffuseFit = coherentFit;
  diffuseFit.fieldTurbulence = 1.f;
  const auto coherent = DefaultCrashCymbalParameters(48000.f, coherentFit);
  const auto diffuse = DefaultCrashCymbalParameters(48000.f, diffuseFit);
  const auto coherentMeasure = measure(coherent.modalField);
  const auto diffuseMeasure = measure(diffuse.modalField);
  Check(coherentMeasure.first == CrashSparseModeCount,
        "zero turbulence leaves exactly one coherent mode per anchor");
  Check(diffuseMeasure.first == CrashModalFieldModeCount,
        "maximum turbulence activates every modal packet member");
  Check(std::abs(coherentMeasure.second - diffuseMeasure.second) < 1.e-7,
        "modal packet expansion preserves normalized drive energy");

  auto selectiveFit = diffuseFit;
  selectiveFit.fieldTurbulenceScale[0] = 0.f;
  const auto selective = DefaultCrashCymbalParameters(48000.f, selectiveFit);
  const auto selectiveMeasure = measure(selective.modalField);
  Check(selectiveMeasure.first == CrashModalFieldModeCount -
            CrashPacketModeCount + 1,
        "a clean anchor retains its centre and disables only its satellites");
  Check(std::abs(selectiveMeasure.second - diffuseMeasure.second) < 1.e-7,
        "per-anchor turbulence changes coherence without changing drive energy");
}

void TestUnifiedFieldAcceptsConstructiveAnchorEditing() {
  using namespace tfdsp::percussion;
  CrashCymbalFitParameters fit;
  fit.sparseFrequencyHz[0] = 7300.f;
  fit.sparseFrequencyHz[1] = 20.f;
  fit.sparseAmplitude[0] = .75f;
  fit.sparseAmplitude[1] = .25f;
  fit.fieldTurbulenceScale.fill(0.f);
  fit.fieldTurbulenceScale[1] = 1.f;
  const auto parameters = DefaultCrashCymbalParameters(48000.f, fit);
  const auto activeInPacket = [&](const std::uint16_t packet) {
    return std::count_if(parameters.modalField.begin(),
                         parameters.modalField.end(),
                         [packet](const auto &mode) {
                           return mode.packet == packet && mode.inputGain != 0.f;
                         });
  };
  Check(activeInPacket(0) == CrashPacketModeCount && activeInPacket(1) == 1,
        "frequency sorting keeps each edited anchor's turbulence attached");

  fit.sparseAmplitude.fill(0.f);
  const auto cleared = DefaultCrashCymbalParameters(48000.f, fit);
  Check(std::none_of(cleared.modalField.begin(), cleared.modalField.end(),
                     [](const auto &mode) { return mode.inputGain != 0.f; }),
        "zero anchor energy silences the unified modal field");
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
              std::isfinite(frame.denseModes) &&
              std::isfinite(frame.turbulentResidual) &&
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
  bool turbulenceAudible = false;
  for (int sample = 0; sample < 4096; ++sample)
  {
    const auto frame = noDense.ProcessFrame();
    denseSilent = denseSilent && frame.denseModes == 0.f;
    turbulenceAudible = turbulenceAudible || frame.turbulentResidual != 0.f;
  }
  Check(denseSilent, "zero dense presentation gain ablates dense modal state");
  Check(turbulenceAudible,
        "turbulence remains independent of dense modal presentation gain");
}

void TestResolvedModesReceiveBloom() {
  using namespace tfdsp::percussion;
  constexpr float sampleRate = 48000.f;
  CrashCymbalFitParameters directFit;
  directFit.directGain = 0.f;
  directFit.sparseGain = 1.f;
  directFit.denseGain = 0.f;
  directFit.sparseBloomGain = 0.f;
  auto bloomFit = directFit;
  bloomFit.sparseBloomGain = 1.f;
  const auto render = [](const CrashCymbalFitParameters &fit) {
    CrashCymbal cymbal;
    cymbal.Prepare(sampleRate, DefaultCrashCymbalParameters(sampleRate, fit));
    cymbal.Trigger({.85f, 1.f, .65f, 79});
    std::vector<float> result(24000);
    for (float &sample : result)
      sample = cymbal.ProcessFrame().sparseModes;
    return result;
  };
  const auto direct = render(directFit);
  const auto bloomed = render(bloomFit);
  Check(Difference(direct, bloomed) > .01 * Energy(bloomed),
        "resolved-only cymbal receives nonlinear bloom excitation");
  Check(CrashCymbalFitParameters{}.sparseBloomGain > 0.f,
        "resolved bloom feed is enabled by default");
}

void TestBloomBodyLevelIsIndependentOfLoopCharacter() {
  using namespace tfdsp::percussion;
  constexpr float sampleRate = 48000.f;
  CrashCymbalFitParameters silentFit;
  silentFit.unifiedBodyEnabled = false;
  silentFit.directGain = 0.f;
  silentFit.sparseGain = 0.f;
  silentFit.denseGain = 1.f;
  silentFit.bodyBypassGain = 0.f;
  silentFit.turbulenceGain.fill(0.f);
  silentFit.bloomBodyGain = 0.f;
  const auto render = [](const CrashCymbalFitParameters &fit) {
    CrashCymbal cymbal;
    cymbal.Prepare(sampleRate, DefaultCrashCymbalParameters(sampleRate, fit));
    cymbal.Trigger({.9f, .8f, .7f, 71});
    double denseEnergy = 0.;
    double bloomEnergy = 0.;
    for (int sample = 0; sample < 24000; ++sample) {
      const auto frame = cymbal.ProcessFrame();
      denseEnergy += static_cast<double>(frame.denseModes) * frame.denseModes;
      bloomEnergy += static_cast<double>(frame.dispersion) * frame.dispersion;
    }
    return std::pair{denseEnergy, bloomEnergy};
  };
  const auto silent = render(silentFit);
  auto audibleFit = silentFit;
  audibleFit.bloomBodyGain = 1.f;
  const auto audible = render(audibleFit);
  Check(silent.first == 0. && audible.first > 0.,
        "bloom body level can completely remove dispersed body excitation");
  Check(silent.second > 0. && audible.second == silent.second,
        "bloom body level does not change dispersion-loop character");
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

void TestZeroStrengthTriggerIsANoOp() {
  using namespace tfdsp::percussion;
  const auto parameters = DefaultCrashCymbalParameters(48000.f);
  CrashCymbal control;
  CrashCymbal probed;
  control.Prepare(48000.f, parameters);
  probed.Prepare(48000.f, parameters);
  control.Trigger({.9f, .8f, .65f, 93});
  probed.Trigger({.9f, .8f, .65f, 93});
  for (int sample = 0; sample < 4096; ++sample) {
    Check(control.Process() == probed.Process(),
          "equal crash states remain sample-identical");
  }
  probed.Trigger({0.f, 0.f, 0.f, 0xffffffffu, 0.f, 1.f});
  for (int sample = 0; sample < 8192; ++sample)
    Check(control.Process() == probed.Process(),
          "zero-strength crash trigger cannot mutate a ringing body");
}

void TestRestrikeAddsWithoutRecolouringStoredEnergy() {
  using namespace tfdsp::percussion;
  constexpr int Onset = 4096;
  constexpr int Tail = 8192;
  CrashCymbalFitParameters fit;
  fit.dispersionDrive = 0.f;
  fit.dispersionExcursionSamples = 0.f;
  fit.turbulenceGain.fill(0.f);
  const auto parameters = DefaultCrashCymbalParameters(48000.f, fit);
  CrashCymbal combined;
  CrashCymbal original;
  CrashCymbal added;
  combined.Prepare(48000.f, parameters);
  original.Prepare(48000.f, parameters);
  added.Prepare(48000.f, parameters);
  combined.Trigger({.95f, .1f, .8f, 101});
  original.Trigger({.95f, .1f, .8f, 101});
  for (int sample = 0; sample < Onset; ++sample) {
    combined.Process();
    original.Process();
    added.Process();
  }
  combined.Trigger({.25f, .95f, .3f, 102});
  added.Trigger({.25f, .95f, .3f, 102});
  double error = 0.0;
  double energy = 0.0;
  for (int sample = 0; sample < Tail; ++sample) {
    const double sum = original.Process() + added.Process();
    const double actual = combined.Process();
    const double delta = actual - sum;
    error += delta * delta;
    energy += actual * actual;
  }
  if (!(error < 1.e-7 * std::max(energy, 1.e-30)))
    std::cerr << "crash linear restrike error/energy: " << error << '/'
              << energy << '\n';
  Check(error < 1.e-7 * std::max(energy, 1.e-30),
        "a new strike adds to, rather than recolours, stored linear body energy");
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
  TestVelocityChangesCymbalRegime();
  TestSparseModesArePlacedDirectly();
  TestBodyDecayCurveUsesActiveErbKnots();
  TestImplementFamiliesAreDistinct();
  TestDefaultBodyCoversTheMeasuredLowRegion();
  TestDenseModeCountUsesNestedExtensionBank();
  TestDenseColourShapesExcitationEnergy();
  TestUnifiedFieldExpandsAnchorsWithoutChangingDriveEnergy();
  TestUnifiedFieldAcceptsConstructiveAnchorEditing();
  TestContactCalibrationMacrosAreAudible();
  TestMuteIsPassive();
  TestFiniteAtSupportedRates();
  TestAnalysisFrameMatchesOutput();
  TestBodyBranchesCanBeAblatedIndependently();
  TestResolvedModesReceiveBloom();
  TestBloomBodyLevelIsIndependentOfLoopCharacter();
  TestMaximumNonlinearBloomRemainsBounded();
  TestZeroStrengthTriggerIsANoOp();
  TestRestrikeAddsWithoutRecolouringStoredEnergy();
  TestRepeatedHitsAccumulateBodyEnergy();
  if (percussion_test::failures == 0)
    std::cout << "All percussion crash tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
