#include "percussion_test_support.hpp"
#include "tfdsp/percussion/crash_cymbal_parameters.hpp"
#include <array>
#include <cmath>

using namespace tfdsp::percussion;
using percussion_test::Check;
using percussion_test::CheckNear;

namespace {
using Field = StochasticModalField<12>;

Field::PreparedParameters PacketField(bool broad) {
  Field::Parameters modes{};
  constexpr std::array<float, 3> centres{100.f, 400.f, 1600.f};
  constexpr std::array<float, 4> narrowFrequency{.98f, .99f, 1.01f, 1.02f};
  constexpr std::array<float, 4> broadFrequency{.4f, .7f, 1.4f, 2.f};
  constexpr std::array<float, 4> narrowEnergy{.1f, .2f, .3f, .4f};
  constexpr std::array<float, 4> broadEnergy{.4f, .2f, .1f, .3f};
  for (std::size_t packet = 0; packet < 3; ++packet) {
    for (std::size_t member = 0; member < 4; ++member) {
      auto &mode = modes[packet*4+member];
      mode.frequencyHz = centres[packet] *
          (broad ? broadFrequency[member] : narrowFrequency[member]);
      mode.inputGain = std::sqrt(broad ? broadEnergy[member] : narrowEnergy[member]);
      mode.packet = static_cast<std::uint16_t>(packet);
      mode.transportFrequencyHz = centres[packet];
    }
  }
  return PrepareStochasticModalField(48000, modes, {}, 700, 6500);
}

void SameTransportDespiteChangingSidebands(bool diffusion) {
  const auto narrow = PacketField(false);
  const auto broad = PacketField(true);
  ModalEnergyCascade<12> first, second;
  const ModalEnergyCascadeParameters controls{3.3f, 1.f, .7f, 19, diffusion};
  first.Prepare(48000, narrow.transportFrequencyHz, narrow.inputGain,
                narrow.packet, 12, controls);
  second.Prepare(48000, broad.transportFrequencyHz, broad.inputGain,
                 broad.packet, 12, controls);
  std::array<float, 12> realA{}, imagA{}, realB{}, imagB{};
  for (std::size_t member = 0; member < 4; ++member) {
    realA[member] = narrow.inputGain[member];
    realB[member] = broad.inputGain[member];
  }
  for (int sample = 0; sample < 4800; ++sample) {
    first.Process(realA, imagA);
    second.Process(realB, imagB);
    if (sample % 480 != 0) continue;
    double total = 0;
    for (std::size_t packet = 0; packet < 3; ++packet) {
      double a = 0, b = 0;
      for (std::size_t member = packet*4; member < packet*4+4; ++member) {
        a += realA[member]*realA[member] + imagA[member]*imagA[member];
        b += realB[member]*realB[member] + imagB[member]*imagB[member];
      }
      CheckNear(a, b, 2.e-5, "sideband layout does not change packet transfer");
      total += a;
    }
    CheckNear(total, 1, 2.e-4, "transport remains passive");
  }
}

void CymbalUsesPaintedCentres() {
  CrashCymbalFitParameters fit;
  fit.fieldRelaxedTurbulence = true;
  fit.fieldTurbulence = 3;
  fit.fieldPacketSpreadErb = 12;
  fit.fieldTurbulenceSlopePerOctave = 1;
  fit.sparseTune = .75f;
  const auto result = DefaultCrashCymbalParameters(48000, fit);
  for (const auto &mode : result.modalField) {
    if (mode.inputGain == 0) continue;
    CheckNear(mode.transportFrequencyHz,
        fit.sparseFrequencyHz[mode.packet] * fit.sparseTune, .01,
        "even boundary-clamped sidebands retain the painted routing coordinate");
  }
  Field::Parameters modes{};
  const auto generic = PrepareStochasticModalField(48000, modes, {}, 700, 6500);
  Check(generic.transportFrequencyHz[0] == generic.frequencyHz[0],
        "generic modes without packet metadata retain their frequency coordinate");
}

void ExcitationMatchesDiffusionCells() {
  CrashCymbalFitParameters fit;
  fit.sparseAmplitude.fill(0);
  fit.sparseFrequencyHz[0] = 100;
  fit.sparseFrequencyHz[1] = 400;
  fit.sparseFrequencyHz[2] = 800;
  for (int i = 0; i < 3; ++i) fit.sparseAmplitude[i] = 1;
  fit.fieldTurbulence = 0;
  fit.bodyTiltDbPerOctave = 0;
  fit.bloomSpectralDiffusion = true;
  auto energy = [](const CrashCymbalParameters &p) {
    std::array<double, 4> result{};
    for (const auto &mode : p.modalField)
      if (mode.packet < 4) result[mode.packet] += mode.inputGain*mode.inputGain;
    return result;
  };
  const auto original = energy(DefaultCrashCymbalParameters(48000, fit));
  CheckNear(original[0], .25, 1.e-6, "flat excitation integrates the low cell");
  CheckNear(original[1], .35, 1.e-6, "flat excitation integrates the middle cell");
  CheckNear(original[2], .4, 1.e-6, "flat excitation integrates the upper cell");
  fit.sparseFrequencyHz[3] = 400;
  fit.sparseAmplitude[3] = 1;
  const auto duplicate = energy(DefaultCrashCymbalParameters(48000, fit));
  CheckNear(duplicate[0], original[0], 1.e-6, "duplicate does not steal low energy");
  CheckNear(duplicate[1]+duplicate[2], original[1], 1.e-6,
      "coincident handles split one excitation cell");
  CheckNear(duplicate[3], original[2], 1.e-6, "duplicate does not steal high energy");
  fit.bodyTiltDbPerOctave = -72;
  fit.bodyExcitationCentreHz = 40;
  fit.sparseFrequencyHz.fill(14000);
  const auto dark = energy(DefaultCrashCymbalParameters(48000, fit));
  CheckNear(dark[0]+dark[1]+dark[2]+dark[3], 1, 1.e-6,
      "extreme dark shelf does not become a hidden body attenuator");
}

void QuietPacketsRetainCoordinatesAndEnergy() {
  const std::array<float, 3> frequencies{100, 400, 1600};
  const std::array<std::uint16_t, 3> packets{0, 1, 2};
  const ModalEnergyCascadeParameters controls{3.3f, 0.f, 0.f, 19, true};
  std::array<float, 3> reference{}, referenceImaginary{};
  for (const bool quietExcitation : {false, true}) {
    for (const float amplitude : {1.f, 1.e-12f}) {
      ModalEnergyCascade<3> cascade;
      const std::array<float, 3> gains{1, quietExcitation ? 1.e-20f : 1.f,
                                        quietExcitation ? 1.e-20f : 1.f};
      cascade.Prepare(48000, frequencies, gains, packets, 3, controls);
      std::array<float, 3> real{amplitude, 0, 0}, imaginary{};
      for (int sample = 0; sample < 480; ++sample)
        cascade.Process(real, imaginary);
      if (!quietExcitation && amplitude == 1.f) {
        reference = real;
        referenceImaginary = imaginary;
      }
      double total = 0;
      for (std::size_t i = 0; i < 3; ++i) {
        const double energy = std::pow(double(real[i]) / amplitude, 2) +
                              std::pow(double(imaginary[i]) / amplitude, 2);
        const double expected = double(reference[i])*reference[i] +
            double(referenceImaginary[i])*referenceImaginary[i];
        CheckNear(energy, expected, 2.e-5,
            "quiet excitation and state preserve linear diffusion coordinates");
        total += energy;
      }
      CheckNear(total, 1, 2.e-5, "quiet packet diffusion conserves energy");
    }
  }
}
}

void SubBassModesKeepTheirPitch() {
  CrashCymbalFitParameters fit;
  fit.sparseAmplitude.fill(0);
  fit.sparseAmplitude[0] = 1;
  fit.fieldTurbulence = 0;
  fit.bodyDecayActive.fill(false);
  for (const float rate : {44100.f, 48000.f, 96000.f}) {
    for (const float hz : {1.f, 10.f, 16.3516f, 27.5f}) {
      fit.sparseFrequencyHz[0] = hz;
      const auto parameters = DefaultCrashCymbalParameters(rate, fit);
      CheckNear(parameters.modalField[0].frequencyHz, hz, 1.e-4,
          "sub-bass centre is not silently raised to 20 Hz");
      CheckNear(parameters.modalField[0].decaySeconds, fit.bodyDecaySeconds[0], 1.e-5,
          "existing low damping extends flat below its endpoint");
      const auto prepared = PrepareStochasticModalField(rate, parameters.modalField,
          StochasticModalFieldControls{}, 700.f, 6500.f);
      const double actual = std::atan2(prepared.sineCentre[0], prepared.cosineCentre[0])
          * rate / (2 * 3.14159265358979323846);
      CheckNear(actual, hz, 1.e-3, "prepared rotation retains accurate low pitch");
    }
  }
}

int main() {
  SubBassModesKeepTheirPitch();
  SameTransportDespiteChangingSidebands(false);
  SameTransportDespiteChangingSidebands(true);
  CymbalUsesPaintedCentres();
  ExcitationMatchesDiffusionCells();
  QuietPacketsRetainCoordinatesAndEnergy();
  return percussion_test::failures ? 1 : 0;
}
