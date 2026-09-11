#include "percussion_test_support.hpp"
#include "tfdsp/percussion/turbulence_profile.hpp"
#include "tfdsp/percussion/crash_cymbal.hpp"
#include <cmath>
#include <limits>

using namespace tfdsp::percussion;
using percussion_test::Check;
using percussion_test::CheckNear;

int main() {
  {
    CrashCymbalFitParameters fit;
    fit.sparseAmplitude.fill(0.f);
    fit.sparseAmplitude[0] = 1.f;
    fit.sparseFrequencyHz[0] = 8000.f;
    fit.fieldDistribution = ModalPacketDistribution::PlateCloud;
    fit.fieldTurbulence = 1.f;
    fit.fieldTurbulenceSlopePerOctave = 0.f;
    fit.fieldPacketSpreadErb = 4.f;
    fit.fieldSatelliteDensity = 1.f;
    const auto p = DefaultCrashCymbalParameters(48000.f, fit);
    double energy = 0;
    std::size_t active = 0;
    for (const auto &m : p.modalField) {
      if (m.inputGain == 0.f) continue;
      ++active;
      energy += double(m.inputGain)*m.inputGain;
      Check(m.packet == 0 && m.transportFrequencyHz == 8000.f,
            "one cloud handle owns all its side modes and transport coordinate");
      Check(m.frequencyHz > 4000.f && m.frequencyHz < 12000.f,
            "plate cloud keeps its wide, bounded support");
    }
    Check(active >= 510, "one handle can use the complete modal pool");
    CheckNear(energy, 1., 3.e-6, "dense single cloud preserves normalized excitation");
    const auto f = [](std::size_t pair) {
      return PacketSideFrequency(8000.f, 4.f, pair, 1.f,
          ModalPacketDistribution::PlateCloud, 2.f, 1.f, 23040.f);
    };
    CheckNear(f(1)-8000.f, .5f*(f(0)-8000.f), .002,
              "plate cloud positions are linear in Hz, not logarithmic");
  }
  for (std::size_t count : {1u, 2u, 3u, 15u, 16u}) {
    for (float depth : {0.f, .1f, .3f, 1.f}) {
      double energy = 0;
      for (std::size_t pair = 0; pair < count; ++pair) {
        const float scale = DoubletWeightScale(pair, count, depth);
        energy += scale * scale;
      }
      CheckNear(energy, double(count), 2.e-6, "doublet depth preserves pair energy including odd counts");
    }
  }
  {
    CrashCymbalFitParameters fit;
    const auto baseline = DefaultCrashCymbalParameters(48000, fit);
    fit.fieldPhaseTilt = -1.f;
    const auto tilted = DefaultCrashCymbalParameters(48000, fit);
    for (std::size_t i = 0; i < baseline.modalField.size(); ++i) {
      const auto &a = baseline.modalField[i], &b = tilted.modalField[i];
      Check(a.frequencyHz == b.frequencyHz && a.inputGain == b.inputGain &&
            a.outputGain == b.outputGain && a.decaySeconds == b.decaySeconds,
            "blur colour leaves placement, damping and energy unchanged");
      if (a.frequencyHz > 0 && a.inputGain != 0)
        CheckNear(b.phaseBandwidthHz, a.phaseBandwidthHz * 1000.f / a.frequencyHz,
                  std::max(.001, double(b.phaseBandwidthHz) * 1.e-5),
                  "blur colour follows explicit 1 kHz pivot");
    }
  }
  for (const auto layout : {ModalPacketDistribution::Scattered,
       ModalPacketDistribution::Even, ModalPacketDistribution::Doublets,
       ModalPacketDistribution::PlateCloud}) {
    for (const float centre : {1.f, 1000.f, 20000.f}) {
      for (std::size_t pair = 0; pair < 100; ++pair) {
        for (const float side : {-1.f, 1.f}) {
          const float f = PacketSideFrequency(centre, 4, pair, side, layout, 6, .96f, 23040);
          Check(f >= 1 && f <= 23040 && std::isfinite(f), "packet support bounded");
        }
      }
    }
  }
  const float first = PacketSideFrequency(5000, 2, 0, 1,
      ModalPacketDistribution::Doublets, 6, 1, 23040);
  const float second = PacketSideFrequency(5000, 2, 1, 1,
      ModalPacketDistribution::Doublets, 6, 1, 23040);
  CheckNear(second - first, 6, .002, "doublet spacing is in Hz, not relative pitch");
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
