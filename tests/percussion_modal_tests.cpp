#include "percussion_test_support.hpp"

#include "tfdsp/percussion/modal_bank.hpp"
#include "tfdsp/percussion/modal_constraint.hpp"
#include "tfdsp/percussion/modal_energy_cascade.hpp"
#include "tfdsp/percussion/modal_packet_allocator.hpp"
#include "tfdsp/percussion/statistical_modal_cloud.hpp"
#include "tfdsp/percussion/stochastic_modal_field.hpp"
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

void TestTwoExcitationsShareOneStoredBody() {
  using Bank = tfdsp::percussion::ModalBank<3>;
  Bank::Parameters parameters{{
      {500.f, 2.f, 1.f, 1.f, 0.f},
      {1300.f, 2.f, 1.f, 1.f, 0.f},
      {2700.f, 2.f, 1.f, 1.f, 0.f},
  }};
  Bank bank;
  bank.Prepare(48000.f, parameters, 700.f, 2000.f);
  bank.SetExcitationProjection({0.f, .25f, 0.f});
  bank.SetSecondaryExcitationProjection({0.f, 0.f, .5f});
  CheckNear(bank.ProcessExcitedPair(1.f, 1.f), .75, 1.e-7,
            "two modal forces add through independent projections");
  Check(std::abs(bank.ProcessExcitedPair(0.f, 0.f)) > 0.f,
        "paired excitation writes one persistent modal state");
}

void TestEnergyNormalizedProjectionOnlyRedistributesDrive() {
  using Field = tfdsp::percussion::StochasticModalField<2>;
  Field::Parameters parameters{{
      {500.f, 2.f, .6f, 1.f, 0.f, 0.f, 0},
      {2000.f, 2.f, .8f, 1.f, 0.f, 0.f, 1},
  }};
  const auto storedAfterStrike = [&](const Field::Projection projection) {
    Field field;
    field.Prepare(48000.f, parameters, {}, 700.f, 1500.f);
    field.SetEnergyNormalizedExcitationProjection(projection);
    field.ProcessExcitedPair(1.f, 0.f);
    return field.StoredEnergy();
  };
  CheckNear(storedAfterStrike({1.f, 0.f}), 1.0, 2.e-7,
            "low-only projection preserves prepared drive energy");
  CheckNear(storedAfterStrike({0.f, .25f}), 1.0, 2.e-7,
            "high-only projection preserves prepared drive energy");
  CheckNear(storedAfterStrike({0.f, 0.f}), 0.0, 0.0,
            "an empty normalized projection remains silent");
}

void TestEnergyNormalizedProjectionSurvivesStaticRebuild() {
  using Field = tfdsp::percussion::StochasticModalField<3>;
  Field::Parameters parameters{{
      {500.f, 2.f, 0.f, 1.f, 0.f, 0.f, 0},
      {1000.f, 2.f, .6f, 1.f, 0.f, 0.f, 1},
      {2000.f, 2.f, .8f, 1.f, 0.f, 0.f, 2},
  }};
  Field field;
  field.Prepare(48000.f, parameters, {}, 700.f, 1500.f);
  field.SetEnergyNormalizedExcitationProjection({0.f, 1.f, 0.f});
  field.SetStaticParameters(parameters);
  field.ProcessExcitedPair(1.f, 0.f);
  CheckNear(field.StoredEnergy(), 1.0, 2.e-7,
            "normalized source-index projection survives a static rebuild");
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
  auto medium = dense;
  medium.modeDensity = .5f;
  const auto mediumModes = MakeStatisticalModalCloud<Count>(48000.f, medium);
  std::size_t active = 0;
  double denseEnergy = 0.0;
  double sparseEnergy = 0.0;
  for (std::size_t mode = 0; mode < Count; ++mode) {
    Check(denseModes[mode].frequencyHz == sparseModes[mode].frequencyHz,
          "modal density does not relocate the cloud");
    active += sparseModes[mode].outputGain != 0.f;
    Check(sparseModes[mode].outputGain == 0.f ||
              mediumModes[mode].outputGain != 0.f,
          "raising modal density only adds to the active subset");
    denseEnergy += denseModes[mode].outputGain * denseModes[mode].outputGain;
    sparseEnergy += sparseModes[mode].outputGain * sparseModes[mode].outputGain;
  }
  Check(active == Count / 4,
        "modal density activates the exact requested subset");
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

void TestStochasticPhaseBroadeningPreservesEnergy() {
  using Field = tfdsp::percussion::StochasticModalField<1>;
  Field::Parameters coherent{{
      {1200.f, 2.f, 1.f, 1.f, 0.f, 0.f, 0},
  }};
  auto diffused = coherent;
  diffused[0].phaseBandwidthHz = 1800.f;
  Field first;
  Field second;
  first.Prepare(48000.f, coherent, {}, 700.f, 6500.f);
  second.Prepare(48000.f, diffused, {}, 700.f, 6500.f);
  double outputDifference = 0.0;
  for (int sample = 0; sample < 4096; ++sample) {
    const float input = sample == 0 ? 1.f : 0.f;
    const double coherentOutput = first.ProcessExcitedPair(input, 0.f);
    const double diffusedOutput = second.ProcessExcitedPair(input, 0.f);
    const double difference = coherentOutput - diffusedOutput;
    outputDifference += difference * difference;
  }
  Check(outputDifference > 1.0,
        "phase bandwidth audibly decorrelates a modal ridge");
  CheckNear(second.StoredEnergy() / first.StoredEnergy(), 1.0, 1.5e-4,
            "phase broadening changes coherence without changing energy");
}

void TestLocalModalExchangeIsPassive() {
  using Field = tfdsp::percussion::StochasticModalField<4>;
  Field::Parameters parameters{{
      {700.f, 3.f, 1.f, 1.f, 0.f, 0.f, 0},
      {900.f, 3.f, .7f, 1.f, .4f, 0.f, 0},
      {1300.f, 3.f, .5f, 1.f, -.3f, 0.f, 1},
      {1700.f, 3.f, .3f, 1.f, .8f, 0.f, 1},
  }};
  Field independent;
  Field coupled;
  independent.Prepare(48000.f, parameters, {}, 600.f, 1500.f);
  coupled.Prepare(48000.f, parameters, {.02f, 17}, 600.f, 1500.f);
  double outputDifference = 0.0;
  double maximumSpontaneousGrowth = 0.0;
  double priorCoupledEnergy = 0.0;
  for (int sample = 0; sample < 4096; ++sample) {
    const float input = sample == 0 ? 1.f : 0.f;
    const double independentOutput = independent.ProcessExcitedPair(input, 0.f);
    const double coupledOutput = coupled.ProcessExcitedPair(input, 0.f);
    const double difference = independentOutput - coupledOutput;
    outputDifference += difference * difference;
    const double energy = coupled.StoredEnergy();
    if (sample > 0)
      maximumSpontaneousGrowth = std::max(
          maximumSpontaneousGrowth, energy - priorCoupledEnergy);
    priorCoupledEnergy = energy;
  }
  Check(outputDifference > 1.0,
        "local exchange moves energy within and between adjacent packets");
  CheckNear(coupled.StoredEnergy() / independent.StoredEnergy(), 1.0, 2.e-4,
            "local Givens exchange preserves stored modal energy");
  Check(maximumSpontaneousGrowth < 2.e-6,
        "local exchange never creates unforced modal energy");
}

void TestModalCascadeIsPassiveAndDoesNotRetransferArrivals() {
  using Cascade = tfdsp::percussion::ModalEnergyCascade<3>;
  const std::array<float, 3> frequencies{500.f, 2000.f, 8000.f};
  const std::array<float, 3> gains{1.f, 1.f, 1.f};
  const std::array<std::uint16_t, 3> packets{0, 1, 2};
  std::array<float, 3> real{1.f, 0.f, 0.f};
  std::array<float, 3> imaginary{};
  Cascade cascade;
  cascade.Prepare(100.f, frequencies, gains, packets, 3,
                  {12.f, 0.f, 1.f, 91});
  const auto energy = [&] {
    double result = 0.0;
    for (std::size_t index = 0; index < real.size(); ++index)
      result += static_cast<double>(real[index]) * real[index] +
          static_cast<double>(imaginary[index]) * imaginary[index];
    return result;
  };
  const double initialEnergy = energy();
  cascade.Process(real, imaginary);
  CheckNear(energy(), initialEnergy, 2.e-7,
            "modal cascade redistributes rather than creates energy");
  Check(real[1] != 0.f && real[2] == 0.f && imaginary[2] == 0.f,
        "modal cascade does not retransfer new arrivals in one sample");
  cascade.Process(real, imaginary);
  Check(real[2] != 0.f || imaginary[2] != 0.f,
        "modal cascade reaches successive bands progressively");
  CheckNear(energy(), initialEnergy, 3.e-7,
            "successive modal transfers remain energy preserving");
}

void TestModalCascadeIsIndependentOfPaintedAnchorDensity() {
  using Cascade = tfdsp::percussion::ModalEnergyCascade<17>;
  const auto highestEnergy = [](const std::size_t count) {
    std::array<float, 17> frequencies{};
    std::array<float, 17> gains{};
    std::array<std::uint16_t, 17> packets{};
    std::array<float, 17> real{};
    std::array<float, 17> imaginary{};
    for (std::size_t index = 0; index < count; ++index) {
      frequencies[index] = 500.f * std::exp2(
          4.f * static_cast<float>(index) /
          static_cast<float>(count - 1));
      gains[index] = 1.f / std::sqrt(static_cast<float>(count));
      packets[index] = static_cast<std::uint16_t>(index);
    }
    real[0] = 1.f;
    Cascade cascade;
    cascade.Prepare(1000.f, frequencies, gains, packets, count,
                    {4.f, 1.f, 0.f, 83});
    for (int sample = 0; sample < 600; ++sample)
      cascade.Process(real, imaginary);
    return real[count - 1] * real[count - 1] +
        imaginary[count - 1] * imaginary[count - 1];
  };
  CheckNear(highestEnergy(9), highestEnergy(17), 2.e-5,
            "intermediate painted anchors do not slow upward travel");
}

void TestModalCascadePreservesPacketTurbulenceWeights() {
  using Cascade = tfdsp::percussion::ModalEnergyCascade<3>;
  const std::array<float, 3> frequencies{500.f, 1000.f, 1100.f};
  const std::array<float, 3> gains{
      1.f, std::sqrt(.75f), std::sqrt(.25f)};
  const std::array<std::uint16_t, 3> packets{0, 1, 1};
  std::array<float, 3> real{1.f, 0.f, 0.f};
  std::array<float, 3> imaginary{};
  Cascade cascade;
  cascade.Prepare(1000.f, frequencies, gains, packets, 3,
                  {4.f, 0.f, 0.f, 89});
  cascade.Process(real, imaginary);
  const float coreEnergy = real[1] * real[1] + imaginary[1] * imaginary[1];
  const float satelliteEnergy =
      real[2] * real[2] + imaginary[2] * imaginary[2];
  CheckNear(coreEnergy / satelliteEnergy, 3.0, 2.e-5,
            "cascade arrivals retain the painted packet's centre-to-turbulence balance");
}

void TestModalFieldRejectsSplitPacketRuns() {
  using namespace tfdsp::percussion;
  using Field = StochasticModalField<3>;
  Field::Parameters parameters{{
      {500.f, 1.f, 1.f, 1.f, 0.f, 0.f, 7, 1.f},
      {1000.f, 1.f, 1.f, 1.f, 0.f, 0.f, 8, 1.f},
      {2000.f, 1.f, 1.f, 1.f, 0.f, 0.f, 7, 1.f},
  }};
  bool rejected = false;
  try {
    (void)PrepareStochasticModalField(
        48000.f, parameters, {}, 700.f, 6500.f);
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  Check(rejected, "modal field rejects non-contiguous packet membership");
}

void TestModalCascadeRateIncludesVerySlowTravel() {
  using Cascade = tfdsp::percussion::ModalEnergyCascade<2>;
  constexpr std::array<float, 2> frequencies{500.f, 1000.f};
  constexpr std::array<float, 2> gains{1.f, 1.f};
  constexpr std::array<std::uint16_t, 2> packets{0, 1};
  const auto upperEnergy = [&](const float rate) {
    std::array<float, 2> real{1.f, 0.f};
    std::array<float, 2> imaginary{};
    Cascade cascade;
    cascade.Prepare(1000.f, frequencies, gains, packets, 2,
                    {rate, 0.f, 0.f, 71});
    for (int sample = 0; sample < 100; ++sample)
      cascade.Process(real, imaginary);
    return real[1] * real[1] + imaginary[1] * imaginary[1];
  };
  const float stopped = upperEnergy(0.f);
  const float verySlow = upperEnergy(.1f);
  const float medium = upperEnergy(1.f);
  const float fast = upperEnergy(4.f);
  Check(stopped == 0.f && verySlow > 0.f && verySlow < .02f &&
            verySlow < medium && medium < fast,
        "modal cascade rate spans off, very slow, medium, and fast travel");
}

void TestModalCascadeEnergyAccelerationOnlyAccelerates() {
  using Cascade = tfdsp::percussion::ModalEnergyCascade<2>;
  constexpr std::array<float, 2> frequencies{500.f, 1000.f};
  constexpr float normalizedGain = .7071067811865475f;
  constexpr std::array<float, 2> gains{normalizedGain, normalizedGain};
  constexpr std::array<std::uint16_t, 2> packets{0, 1};
  const auto upperEnergy = [&](const float strength,
                               const float energyAcceleration) {
    std::array<float, 2> real{strength, 0.f};
    std::array<float, 2> imaginary{};
    Cascade cascade;
    cascade.Prepare(1000.f, frequencies, gains, packets, 2,
                    {2.f, energyAcceleration, 0.f, 97});
    for (int sample = 0; sample < 100; ++sample)
      cascade.Process(real, imaginary);
    return real[1] * real[1] + imaginary[1] * imaginary[1];
  };
  const float baseline = upperEnergy(.001f, 0.f);
  const float quietDependent = upperEnergy(.001f, 1.f);
  const float loudDependent = upperEnergy(1.f, 1.f);
  Check(quietDependent >= .999f * baseline,
        "energy acceleration never reduces the declared baseline rate");
  Check(loudDependent > 1.25f * upperEnergy(1.f, 0.f),
        "stored strike energy accelerates the field-wide cascade");
}

void TestModalPacketAllocationUsesOneSharedBoundedPool() {
  using namespace tfdsp::percussion;
  std::array<ModalPacketRequest, 32> requests{};
  for (std::size_t index = 0; index < 12; ++index)
    requests[index] = {5.f + static_cast<float>(index), 6.f, true};
  const auto centresOnly = AllocateModalPackets(requests, 512, 0.f);
  const auto dense = AllocateModalPackets(requests, 512, 1.f);
  Check(centresOnly.activeHandleCount == 12 &&
            centresOnly.stateCount == 12,
        "zero satellite density retains only painted centre handles");
  Check(dense.activeHandleCount == 12 && dense.stateCount <= 512 &&
            dense.stateCount > centresOnly.stateCount,
        "satellites share one bounded modal-state pool");
  for (std::size_t index = 12; index < requests.size(); ++index)
    Check(dense.sidebandPairs[index] == 0,
          "inactive handles never consume sideband states");
}

} // namespace

int main() {
  TestOneModeMatchesAnalyticRecurrence();
  TestIndependentProjections();
  TestPreparedProjectionIsSanitizedOnce();
  TestTwoExcitationsShareOneStoredBody();
  TestEnergyNormalizedProjectionOnlyRedistributesDrive();
  TestEnergyNormalizedProjectionSurvivesStaticRebuild();
  TestModalConstraintReferenceGain();
  TestStatisticalCloudIsDeterministicAndNormalized();
  TestDenseCloudRemainsFinite();
  TestDenseCloudGainEnvelopeShapesBroadBands();
  TestDenseCloudDecayEnvelopeShapesMiddleBand();
  TestDenseCloudDensityPreservesPlacementAndLevel();
  TestTurbulentResidualStoresAndPassivelyLosesEnergy();
  TestTurbulentResidualIsRepeatable();
  TestStochasticPhaseBroadeningPreservesEnergy();
  TestLocalModalExchangeIsPassive();
  TestModalCascadeIsPassiveAndDoesNotRetransferArrivals();
  TestModalCascadeIsIndependentOfPaintedAnchorDensity();
  TestModalCascadePreservesPacketTurbulenceWeights();
  TestModalFieldRejectsSplitPacketRuns();
  TestModalCascadeRateIncludesVerySlowTravel();
  TestModalCascadeEnergyAccelerationOnlyAccelerates();
  TestModalPacketAllocationUsesOneSharedBoundedPool();
  if (percussion_test::failures == 0)
    std::cout << "All percussion modal tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
