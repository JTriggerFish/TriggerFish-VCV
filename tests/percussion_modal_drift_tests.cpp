#include "percussion_test_support.hpp"
#include "tfdsp/percussion/smooth_modal_drift.hpp"
#include "tfdsp/percussion/stochastic_modal_field.hpp"
#include <array>
#include <cmath>
#include <limits>

using percussion_test::Check;
using percussion_test::CheckNear;
using tfdsp::percussion::SmoothModalDrift;

void TestTrajectory() {
  SmoothModalDrift<2> drift;
  drift.Prepare(48000.f, {10000.f, 10000.f}, {1.f, 0.f}, 5.f, 8.f, 7);
  std::array<float, 48000> first{};
  float maximumStep = 0.f, prior = 0.f;
  for (std::size_t i = 0; i < first.size(); ++i) {
    first[i] = drift.NextAngle(0);
    Check(std::abs(first[i]) <= 6.283186f * 500.f / 48000.f,
          "frequency drift is bounded by stated depth");
    Check(drift.NextAngle(1) == 0.f, "clean modes stay fixed");
    if (i) maximumStep = std::max(maximumStep, std::abs(first[i]-prior));
    prior = first[i];
  }
  Check(maximumStep < .00005f, "smooth across random target boundaries");
  drift.Reset();
  for (std::size_t i = 0; i < first.size(); ++i) {
    Check(first[i] == drift.NextAngle(0), "reset reproduces drift trajectory");
    drift.NextAngle(1);
  }
}

void TestEnergyAndPreparedReplay() {
  using Field = tfdsp::percussion::StochasticModalField<2>;
  Field::Parameters modes{{{700.f, 3.f, 1.f, 1.f, 0.f, 0.f},
                                {8000.f, 3.f, .5f, 1.f, 0.f, 0.f}}};
  for (auto &mode : modes) mode.exchangeAmount = 0.f;
  tfdsp::percussion::StochasticModalFieldControls controls{};
  controls.driftDepthHz = 10.f;
  controls.driftKnotsPerSecond = 40.f;
  Field field, replay;
  field.Prepare(48000.f, modes, controls, 500.f, 5000.f);
  replay.LoadPrepared(tfdsp::percussion::PrepareStochasticModalField(
      48000.f, modes, controls, 500.f, 5000.f));
  double priorEnergy = 1.25;
  float peakDifference = 0.f;
  Field staticField;
  staticField.Prepare(48000.f, modes, {}, 500.f, 5000.f);
  for (int i = 0; i < 48000; ++i) {
    const float y = field.ProcessExcitedPair(i == 0 ? 1.f : 0.f, 0.f);
    Check(y == replay.ProcessExcitedPair(i == 0 ? 1.f : 0.f, 0.f),
          "prepared renderer reproduces drift");
    peakDifference = std::max(peakDifference, std::abs(y -
        staticField.ProcessExcitedPair(i == 0 ? 1.f : 0.f, 0.f)));
    const double energy = field.StoredEnergy();
    Check(energy <= priorEnergy * 1.000001, "drift cannot inject tail energy");
    priorEnergy = energy;
  }
  Check(peakDifference > .1f, "nonzero drift actually affects audio");
  CheckNear(field.StoredEnergy(), 1.25 * std::pow(.001, 2./3.), .0002,
            "drift retains the declared T60");
}

void TestAbsoluteHz() {
  SmoothModalDrift<3> drift;
  drift.PrepareHz(48000.f, {125.f, 8000.f, .5f}, 1.f, .5f, 9);
  std::array<double, 2> squared{};
  constexpr float hzPerAngle = 48000.f / 6.28318530718f;
  for (int sample=0; sample<480000; ++sample) {
    for (std::size_t i=0; i<3; ++i) {
      const float hz = drift.NextAngle(i)*hzPerAngle;
      Check(std::abs(hz) <= (i==2 ? .5f : 1.f), "Hz wander and boundary remain bounded");
      if (i<2) squared[i] += hz*hz;
    }
  }
  Check(squared[0]>1000 && squared[1]>1000, "both bass and treble wander");
  Check(squared[1]/squared[0]<4, "treble does not inherit fractional-frequency broadening");
}

int main() {
  TestTrajectory();
  TestAbsoluteHz();
  TestEnergyAndPreparedReplay();
  SmoothModalDrift<1> drift;
  drift.Prepare(48000.f, {std::numeric_limits<float>::infinity()}, {1.f},
                0.f, 8.f, 1);
  Check(!drift.Enabled() && drift.NextAngle(0) == 0.f,
        "zero depth and malformed frequency are safe");
  return percussion_test::failures ? 1 : 0;
}
