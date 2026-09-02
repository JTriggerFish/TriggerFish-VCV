#include "percussion_test_support.hpp"

#include "tfdsp/percussion/modal_bank.hpp"
#include "tfdsp/percussion/modal_constraint.hpp"
#include "tfdsp/percussion/statistical_modal_cloud.hpp"
#include "tfdsp/percussion/turbulent_residual.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

void TestOneModeMatchesAnalyticRecurrence() {
  using Bank = tfdsp::percussion::ModalBank<1>;
  constexpr float sampleRate = 48000.f;
  constexpr float frequency = 1000.f;
  constexpr float decay = 1.f;
  Bank::Parameters parameters{{{frequency, decay, 1.f, 1.f, 0.f}}};
  Bank bank;
  bank.Prepare(sampleRate, parameters, 500.f, 5000.f);
  constexpr double TwoPi = 6.28318530717958647692;
  const double radius = std::exp(std::log(.001) / (decay * sampleRate));
  double maximumError = 0.0;
  double oneSecond = 0.0;
  for (std::size_t sample = 0; sample <= 48000; ++sample) {
    const double actual = bank.Process(sample == 0 ? 1.f : 0.f);
    const double expected = std::pow(radius, static_cast<double>(sample)) *
        std::cos(TwoPi * frequency * static_cast<double>(sample) / sampleRate);
    maximumError = std::max(maximumError, std::abs(actual - expected));
    if (sample == 48000)
      oneSecond = actual;
  }
  CheckNear(maximumError, 0.0, 2.e-4,
            "modal recurrence matches its declared analytic response");
  CheckNear(oneSecond, .001, 2.e-6,
            "modal recurrence reaches -60 dB at its declared T60");
}

void TestIndependentProjections() {
  using Bank = tfdsp::percussion::ModalBank<3>;
  Bank::Parameters parameters{{
      {500.f, 2.f, 1.f, 1.f, 0.f},
      {1300.f, 2.f, 1.f, 1.f, 0.f},
      {2700.f, 2.f, 1.f, 1.f, 0.f},
  }};
  Bank bank;
  bank.Prepare(48000.f, parameters, 700.f, 2000.f);
  const Bank::Projection excitation{0.f, .25f, 0.f};
  const Bank::Projection observation{1.f, 2.f, 1.f};
  CheckNear(bank.ProcessProjected(1.f, excitation, observation),
            .5, 1.e-7,
            "modal excitation and observation projections remain independent");
}

void TestPreparedProjectionIsSanitizedOnce() {
  using Bank = tfdsp::percussion::ModalBank<3>;
  Bank::Parameters parameters{{
      {500.f, 2.f, 1.f, 1.f, 0.f},
      {1300.f, 2.f, 1.f, 1.f, 0.f},
      {2700.f, 2.f, 1.f, 1.f, 0.f},
  }};
  Bank prepared;
  Bank reference;
  prepared.Prepare(48000.f, parameters, 700.f, 2000.f);
  reference.Prepare(48000.f, parameters, 700.f, 2000.f);
  const Bank::Projection unsafe{
      std::numeric_limits<float>::quiet_NaN(),
      std::numeric_limits<float>::infinity(), .25f};
  const Bank::Projection sanitized{0.f, 0.f, .25f};
  prepared.SetExcitationProjection(unsafe);
  for (int sample = 0; sample < 128; ++sample) {
    const float input = sample == 0 ? 1.f : 0.f;
    Check(prepared.ProcessExcited(input) ==
              reference.ProcessExcited(input, sanitized),
          "prepared modal projection retains per-sample safety semantics");
  }
}

void TestModalConstraintReferenceGain() {
  using namespace tfdsp::percussion;
  ModalConstraintController controller;
  controller.Prepare(48000.f, .001f, 0.f, 0.f);
  controller.SetTarget({.5f, .25f, .125f, .0625f});
  double broadband = 1.0;
  double low = 1.0;
  for (int sample = 0; sample < 48; ++sample) {
    const auto gain = controller.Process();
    broadband *= gain.broadband;
    low *= gain.low;
  }
  CheckNear(broadband, .5, 2.e-6,
            "modal constraint realizes broadband traversal attenuation");
  CheckNear(low, .25, 2.e-6,
            "modal constraint realizes band traversal attenuation");
}

void TestStatisticalCloudIsDeterministicAndNormalized() {
  using namespace tfdsp::percussion;
  constexpr std::size_t Count = 256;
  StatisticalModalCloudParameters parameters;
  parameters.outputGain = .7f;
  const auto first = MakeStatisticalModalCloud<Count>(48000.f, parameters);
  const auto repeated = MakeStatisticalModalCloud<Count>(48000.f, parameters);
  parameters.seed += 1;
  const auto different = MakeStatisticalModalCloud<Count>(48000.f, parameters);
  double squaredGain = 0.0;
  bool ordered = true;
  bool changed = false;
  for (std::size_t mode = 0; mode < Count; ++mode) {
    Check(first[mode].frequencyHz >= 650.f &&
              first[mode].frequencyHz <= 18000.f,
          "statistical modal frequencies stay inside the requested range");
    if (mode > 0)
      ordered = ordered &&
          first[mode].frequencyHz > first[mode - 1].frequencyHz;
    squaredGain += first[mode].outputGain * first[mode].outputGain;
    changed = changed ||
        first[mode].frequencyHz != different[mode].frequencyHz;
    Check(first[mode].frequencyHz == repeated[mode].frequencyHz &&
              first[mode].decaySeconds == repeated[mode].decaySeconds &&
              first[mode].outputGain == repeated[mode].outputGain &&
              first[mode].inputPhaseRadians ==
                  repeated[mode].inputPhaseRadians,
          "statistical modal cloud repeats exactly for one object seed");
  }
  Check(ordered, "statistical modal frequencies remain strictly ordered");
  Check(changed, "statistical modal object seed changes nuisance frequencies");
  CheckNear(squaredGain, .49, 2.e-6,
            "statistical modal output weights have count-independent energy");
}

void TestDenseCloudRemainsFinite() {
  using Bank = tfdsp::percussion::ModalBank<256>;
  const auto parameters =
      tfdsp::percussion::MakeStatisticalModalCloud<256>(48000.f, {});
  Bank bank;
  bank.Prepare(48000.f, parameters, 700.f, 6500.f);
  for (std::size_t sample = 0; sample < 480000; ++sample) {
    const float output = bank.Process(sample % 6000 == 0 ? .2f : 0.f);
    Check(std::isfinite(output) && std::abs(output) < 100.f,
          "dense modal cloud stays finite under repeated excitation");
  }
}

void TestDenseCloudGainEnvelopeShapesBroadBands() {
  using namespace tfdsp::percussion;
  constexpr std::size_t Count = 256;
  StatisticalModalCloudParameters flat;
  flat.tiltDbPerOctave = 0.f;
  flat.gainSpreadDb = 0.f;
  auto shaped = flat;
  for (std::size_t point = 0; point < shaped.gainEnvelopeDb.size(); ++point)
    shaped.gainEnvelopeDb[point] = -12.f + 24.f * static_cast<float>(point) /
        static_cast<float>(shaped.gainEnvelopeDb.size() - 1);
  const auto flatModes = MakeStatisticalModalCloud<Count>(48000.f, flat);
  const auto shapedModes = MakeStatisticalModalCloud<Count>(48000.f, shaped);
  double flatLow = 0.0;
  double flatHigh = 0.0;
  double shapedLow = 0.0;
  double shapedHigh = 0.0;
  for (std::size_t mode = 0; mode < Count / 4; ++mode) {
    flatLow += std::abs(flatModes[mode].outputGain);
    shapedLow += std::abs(shapedModes[mode].outputGain);
  }
  for (std::size_t mode = 3 * Count / 4; mode < Count; ++mode) {
    flatHigh += std::abs(flatModes[mode].outputGain);
    shapedHigh += std::abs(shapedModes[mode].outputGain);
  }
  Check(shapedHigh / shapedLow > 4.0 * flatHigh / flatLow,
        "modal-cloud gain envelope changes broad spectral balance");
}

void TestDenseCloudDecayEnvelopeShapesMiddleBand() {
  using namespace tfdsp::percussion;
  constexpr std::size_t Count = 256;
  StatisticalModalCloudParameters flat;
  flat.decaySpreadOctaves = 0.f;
  flat.lowDecaySeconds = 1.f;
  flat.highDecaySeconds = 1.f;
  auto shaped = flat;
  shaped.decayEnvelopeOctaves = {0.f, 0.f, 1.f, 1.f, 0.f, 0.f};

  const auto flatModes = MakeStatisticalModalCloud<Count>(48000.f, flat);
  const auto shapedModes = MakeStatisticalModalCloud<Count>(48000.f, shaped);
  Check(shapedModes[Count / 2].decaySeconds >
            1.8f * flatModes[Count / 2].decaySeconds,
        "modal-cloud decay envelope changes middle-band loss");
  CheckNear(shapedModes.front().decaySeconds,
            flatModes.front().decaySeconds, .02,
            "modal-cloud decay envelope preserves the low endpoint");
  CheckNear(shapedModes.back().decaySeconds,
            flatModes.back().decaySeconds, .02,
            "modal-cloud decay envelope preserves the high endpoint");
}

void TestDenseCloudDensityPreservesPlacementAndLevel() {
  using namespace tfdsp::percussion;
  constexpr std::size_t Count = 256;
  StatisticalModalCloudParameters dense;
  dense.modeDensity = 1.f;
  auto sparse = dense;
  sparse.modeDensity = .25f;
  const auto denseModes = MakeStatisticalModalCloud<Count>(48000.f, dense);
  const auto sparseModes = MakeStatisticalModalCloud<Count>(48000.f, sparse);
  std::size_t active = 0;
  double denseEnergy = 0.0;
  double sparseEnergy = 0.0;
  for (std::size_t mode = 0; mode < Count; ++mode) {
    Check(denseModes[mode].frequencyHz == sparseModes[mode].frequencyHz,
          "modal density does not relocate the cloud");
    active += sparseModes[mode].outputGain != 0.f;
    denseEnergy += denseModes[mode].outputGain * denseModes[mode].outputGain;
    sparseEnergy += sparseModes[mode].outputGain * sparseModes[mode].outputGain;
  }
  Check(active > 40 && active < 90,
        "modal density activates the expected stable subset");
  CheckNear(sparseEnergy, denseEnergy, 1.e-4,
            "modal density preserves normalized cloud level");
}

void TestTurbulentResidualStoresAndPassivelyLosesEnergy() {
  using namespace tfdsp::percussion;
  TurbulentResidualParameters parameters;
  parameters.gain = {1.f, 1.f, 1.f};
  parameters.decay = {.2f, .2f, .2f};
  TurbulentResidual residual;
  residual.Prepare(48000.f, parameters);
  Check(residual.Process(0.f) == 0.f,
        "unexcited turbulent residual is exactly silent");
  for (int sample = 0; sample < 256; ++sample)
    residual.Process(.5f);
  const float charged = residual.StoredEnergy();
  for (int sample = 0; sample < 256; ++sample)
    residual.Process(0.f, {.5f, 1.f, 1.f, 1.f});
  Check(charged > 0.f && residual.StoredEnergy() < charged,
        "turbulent residual stores drive energy and constraints only remove it");
}

void TestTurbulentResidualIsRepeatable() {
  using namespace tfdsp::percussion;
  TurbulentResidualParameters parameters;
  parameters.gain = {.2f, .2f, .2f};
  TurbulentResidual first;
  TurbulentResidual second;
  first.Prepare(48000.f, parameters);
  second.Prepare(48000.f, parameters);
  for (int sample = 0; sample < 4096; ++sample) {
    const float drive = sample < 512 ? .25f : 0.f;
    Check(first.Process(drive) == second.Process(drive),
          "turbulent residual repeats exactly for one seed");
  }
}

} // namespace

int main() {
  TestOneModeMatchesAnalyticRecurrence();
  TestIndependentProjections();
  TestPreparedProjectionIsSanitizedOnce();
  TestModalConstraintReferenceGain();
  TestStatisticalCloudIsDeterministicAndNormalized();
  TestDenseCloudRemainsFinite();
  TestDenseCloudGainEnvelopeShapesBroadBands();
  TestDenseCloudDecayEnvelopeShapesMiddleBand();
  TestDenseCloudDensityPreservesPlacementAndLevel();
  TestTurbulentResidualStoresAndPassivelyLosesEnergy();
  TestTurbulentResidualIsRepeatable();
  if (percussion_test::failures == 0)
    std::cout << "All percussion modal tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
