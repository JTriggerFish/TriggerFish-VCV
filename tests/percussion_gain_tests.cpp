#include "percussion_test_support.hpp"
#include "tfdsp/percussion/contact_exciter.hpp"
#include "tfdsp/percussion/crash_cymbal.hpp"
#include "tfdsp/percussion/membrane_resonator.hpp"

#include <complex>

using namespace tfdsp::percussion;
using percussion_test::Check;
using percussion_test::CheckNear;

namespace {

ContactExciterParameters PulseOnly() {
  ContactExciterParameters result;
  result.chirp.amplitude = result.noise.amplitude =
      result.microContacts.amplitude = 0.f;
  return result;
}

void TestImpulseUnits() {
  for (const float rate : {44100.f, 48000.f, 96000.f, 192000.f}) {
    for (const float duration : {.0002f, .001f, .008f, .1f}) {
      ContactExciter contact;
      contact.Prepare(rate);
      auto parameters = PulseOnly();
      parameters.pulseDurationSeconds = duration;
      parameters.pulseAmplitude = 4.f;
      contact.Trigger(parameters);
      double impulse = 0;
      while (contact.Active()) impulse += contact.Process().bodyDrive;
      CheckNear(impulse, 4, .001, "pulse impulse is independent of duration/rate");
    }
  }
}

// A normalized, undamped test oscillator measures excitation in one fixed
// physical band. Comparing raw noise sum-of-squares would miss the old bug.
double NoiseResponse(const float rate, const bool micro) {
  double sum = 0;
  for (unsigned seed = 1; seed <= 128; ++seed) {
    ContactExciter contact;
    contact.Prepare(rate);
    auto p = PulseOnly();
    p.pulseAmplitude = 0;
    if (micro) {
      p.microContacts.amplitude = 1;
      p.microContacts.seed = seed;
      p.microContacts.durationSeconds = .05f;
    } else {
      p.noise.amplitude = 1;
      p.noise.seed = seed;
      p.noise.tiltDb = 0;
    }
    contact.Trigger(p);
    std::complex<double> state{};
    const auto rotation = std::polar(1.0, 6.283185307179586 * 300 / rate);
    while (contact.Active())
      state = rotation * state + static_cast<double>(contact.Process().bodyDrive);
    sum += std::norm(state);
  }
  return sum / 128;
}

void TestNoiseUnits() {
  for (const bool micro : {false, true}) {
    const double low = NoiseResponse(48000, micro);
    const double high = NoiseResponse(96000, micro);
    if (!(high / low > .65 && high / low < 1.45))
      std::cerr << "noise modal variance 48/96k, micro=" << micro << ": "
                << low << '/' << high << '\n';
    Check(high / low > .65 && high / low < 1.45,
          "stochastic excitation has sample-rate-consistent band energy");
  }
}

CrashCymbalFitParameters SingleMode(const float drive) {
  CrashCymbalFitParameters fit;
  fit.sparseAmplitude.fill(0);
  fit.sparseAmplitude[0] = 1;
  fit.sparseFrequencyHz[0] = 80;
  fit.fieldTurbulence = 0;
  fit.bodyExcitationGain = drive;
  fit.fieldGain = fit.outputGain = 1;
  fit.directGain = 0;
  fit.bodyRadiationEnabled = false;
  fit.contactChirpGain = fit.contactNoiseGain = fit.contactMicroGain = 0;
  return fit;
}

double SingleModeEnergy(const float rate, const float drive) {
  CrashCymbal voice;
  voice.Prepare(rate, DefaultCrashCymbalParameters(rate, SingleMode(drive)));
  voice.Trigger({.8f, .5f, .65f, 1, 1.f, .2f});
  double total = 0;
  double error = 0;
  float y1 = 0, y2 = 0;
  const auto parameters = DefaultCrashCymbalParameters(rate, SingleMode(drive));
  const double radius = std::exp(std::log(.001) /
      (parameters.modalField[0].decaySeconds * rate));
  const double coefficient = 2 * radius * std::cos(6.283185307179586 * 80 / rate);
  for (int sample = 0; sample < static_cast<int>(.2f * rate); ++sample) {
    const float y = voice.Process();
    if (sample > .02f * rate) {
      total += y*y;
      const double residual = y - coefficient*y1 + radius*radius*y2;
      error += residual*residual;
    }
    y2 = y1;
    y1 = y;
  }
  Check(error < 1.e-9 * total, "single mode remains a clean damped sine at 4x");
  return total / rate;
}

void TestSingleMode() {
  const double nominal = SingleModeEnergy(48000, 1);
  CheckNear(SingleModeEnergy(48000, 4) / nominal, 16, .001,
            "4x drive means 16x isolated linear energy, without compression");
  CheckNear(SingleModeEnergy(96000, 1) / nominal, 1, .03,
            "same strike has same modal response at 48/96k");
}

void TestDensityAndObservation() {
  auto fit = SingleMode(1);
  fit.fieldTurbulence = 1;
  for (const float density : {0.f, .1f, .5f, 1.f}) {
    fit.fieldSatelliteDensity = density;
    for (const float bar : {.25f, .5f, 1.f}) {
      fit.sparseAmplitude[0] = bar;
      const auto prepared = DefaultCrashCymbalParameters(48000, fit);
      double expectedPower = 0, inputNorm = 0;
      for (const auto &mode : prepared.modalField) {
        inputNorm += mode.inputGain * mode.inputGain;
        expectedPower += std::pow(mode.inputGain * mode.outputGain, 2);
      }
      CheckNear(inputNorm, 1, 1.e-6, "density and paint do not change input norm");
      CheckNear(expectedPower, bar*bar, 1.e-6,
                "packet refinement preserves incoherent power and literal bar gain");
    }
  }
}

void TestNoHiddenEnergyCeiling() {
  MembraneResonator<1> small, large;
  MembraneResonator<1>::Parameters modes{};
  small.Prepare(48000, modes);
  large.Prepare(48000, modes);
  small.Process({128}, 1);
  large.Process({256}, 1);
  Check(large.StoredEnergy() > 1000, "valid input is not clipped at an energy ceiling");
  float previous = large.StoredEnergy();
  for (int sample = 0; sample < 12000; ++sample) {
    const float low = small.Process({0}, 1);
    const float high = large.Process({0}, 1);
    CheckNear(high, 2*low, 1.e-5, "high-energy membrane remains linear");
    Check(large.StoredEnergy() <= previous,
          "unforced high-energy state dissipates without an energy cap");
    previous = large.StoredEnergy();
  }
}

} // namespace

int main() {
  TestImpulseUnits();
  TestNoiseUnits();
  TestSingleMode();
  TestDensityAndObservation();
  TestNoHiddenEnergyCeiling();
  return percussion_test::failures == 0 ? 0 : 1;
}
