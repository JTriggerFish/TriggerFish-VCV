#include "percussion_test_support.hpp"
#include "tfdsp/percussion/turbulence_profile.hpp"
#include "tfdsp/percussion/crash_cymbal.hpp"
#include <cmath>
#include <limits>

using namespace tfdsp::percussion;
using percussion_test::Check;
using percussion_test::CheckNear;

int main() {
  const auto at = [](float f, float local, bool relaxed) {
    return EvaluateTurbulence(f, 1.f, 1.f, 1000.f, local, relaxed);
  };
  CheckNear(at(2000, 1, true).intensity, 2, 1.e-6, "relaxed grows past one");
  CheckNear(at(2000, .1f, true).intensity, .2, 1.e-6, "local scales before energy mapping");
  CheckNear(at(2000, .1f, false).intensity, .1, 1.e-6, "classic remains capped");
  CheckNear(at(1000, 1, true).diffuseEnergy, .9, 1.e-6, "90 percent energy at one");
  Check(at(4000, 0, true).intensity == 0, "local zero stays coherent");
  Check(EvaluateTurbulence(15000, 0, 1, 40, 2, true).intensity == 0,
        "global zero disables turbulence at every frequency");
  float prior = -1.f;
  for (int i = 0; i <= 1000; ++i) {
    const auto r = EvaluateTurbulence(1000, i / 250.f, 0, 1000, 1, true);
    Check(r.diffuseEnergy >= prior && r.diffuseEnergy <= 1.f,
          "diffuse energy is bounded and monotonic");
    CheckNear((1.f-r.diffuseEnergy)+r.diffuseEnergy, 1.f, 1.e-7,
              "core plus diffuse energy is normalized");
    prior = r.diffuseEnergy;
  }
  for (float density : {0.f, .5f, 1.f}) {
    for (float level : {0.f, .7f, 1.f, 2.f, 4.f}) {
      CrashCymbalFitParameters fit;
      fit.fieldRelaxedTurbulence = true;
      fit.fieldTurbulence = level;
      fit.fieldTurbulenceSlopePerOctave = 1;
      fit.fieldTurbulenceCentreHz = 1000;
      fit.fieldSatelliteDensity = density;
      fit.bloomRateOctavesPerSecond = 0;
      const auto parameters = DefaultCrashCymbalParameters(48000, fit);
      double inputEnergy = 0;
      for (const auto &m : parameters.modalField) {
        Check(std::isfinite(m.inputGain) && std::isfinite(m.phaseBandwidthHz),
              "prepared modes remain finite across expanded range");
        inputEnergy += static_cast<double>(m.inputGain) * m.inputGain;
      }
      CheckNear(inputEnergy, 1, 3.e-6, "allocation preserves normalized drive energy");
      CrashModalField field;
      field.Prepare(48000, parameters.modalField, parameters.modalFieldControls, 700, 6500);
      field.ProcessExcitedPair(1, 0);
      double energy = field.StoredEnergy();
      for (int i=0;i<4800;++i) {
        const auto y = field.ProcessExcitedPair(0, 0);
        Check(std::isfinite(y), "relaxed field output is finite");
        const auto next = field.StoredEnergy();
        Check(next <= energy * 1.000002, "no energy generation in unforced tail");
        energy = next;
      }
    }
  }
  Check(std::isfinite(EvaluateTurbulence(std::numeric_limits<float>::infinity(),
      NAN, NAN, NAN, NAN, true).intensity), "invalid inputs are sanitized");
  return percussion_test::failures ? 1 : 0;
}
