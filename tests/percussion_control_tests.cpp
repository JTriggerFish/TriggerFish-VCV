#include "percussion_test_support.hpp"

#include "tfdsp/percussion/breakpoint_trajectory.hpp"
#include "tfdsp/percussion/passive_constraint.hpp"

#include <array>
#include <cmath>
#include <limits>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

void TestLinearAndGeometricSegments() {
  using Trajectory = tfdsp::percussion::BreakpointTrajectory<4>;
  using Curve = tfdsp::percussion::TrajectoryCurve;
  Trajectory trajectory;
  trajectory.Prepare(1000.f);
  Trajectory::Segments segments{};
  segments[0] = {1.f, .010f, Curve::Linear};
  segments[1] = {4.f, .002f, Curve::Geometric};
  trajectory.Start(0.f, segments, 2);
  for (int sample = 1; sample <= 10; ++sample)
    CheckNear(trajectory.Process(), .1 * sample, 2.e-6,
              "linear trajectory reaches every expected sample");
  CheckNear(trajectory.Process(), 2.0, 2.e-6,
            "geometric trajectory follows a constant ratio");
  CheckNear(trajectory.Process(), 4.0, 2.e-6,
            "geometric trajectory reaches its endpoint exactly");
  Check(!trajectory.Active(), "trajectory becomes inactive after its last segment");
}

void TestTrajectoryRetriggerAndSanitization() {
  using Trajectory = tfdsp::percussion::BreakpointTrajectory<3>;
  Trajectory trajectory;
  trajectory.Prepare(1000.f);
  Trajectory::Segments segments{};
  segments[0] = {1.f, .010f};
  trajectory.Start(0.f, segments, 1);
  for (int sample = 0; sample < 5; ++sample)
    (void)trajectory.Process();
  CheckNear(trajectory.Value(), .5, 1.e-6,
            "trajectory exposes its current retrigger value");
  segments[0] = {0.f, .005f};
  trajectory.RetriggerFromCurrent(segments, 1);
  CheckNear(trajectory.Process(), .4, 1.e-6,
            "retrigger starts continuously from the current value");

  segments[0] = {std::numeric_limits<float>::quiet_NaN(), 0.f};
  trajectory.Start(std::numeric_limits<float>::infinity(), segments, 1);
  Check(trajectory.Value() == 0.f && !trajectory.Active(),
        "trajectory sanitizes non-finite and zero-duration input");

  segments[0] = {-std::numeric_limits<float>::max(), .002f};
  trajectory.Start(std::numeric_limits<float>::max(), segments, 1);
  Check(std::isfinite(trajectory.Process()) &&
            std::isfinite(trajectory.Process()),
        "trajectory remains finite between opposite extreme endpoints");
}

void TestPassiveConstraintSmoothing() {
  using namespace tfdsp::percussion;
  DynamicLossController controller;
  controller.Prepare(1000.f, 0.f, .1f);
  controller.SetTarget({.8f, .9f, .7f, .2f});
  const PassiveConstraintGains attacked = controller.Process();
  CheckNear(attacked.broadband, .8, 1.e-7,
            "constraint attack can apply within one sample");
  CheckNear(attacked.high, .2, 1.e-7,
            "constraint attack applies every band together");

  controller.SetTarget({1.f, 1.f, 1.f, 1.f});
  const PassiveConstraintGains firstRelease = controller.Process();
  Check(firstRelease.broadband > attacked.broadband &&
            firstRelease.broadband < 1.f,
        "constraint release is smooth and monotonic");
  for (int sample = 0; sample < 1000; ++sample)
    (void)controller.Process();
  CheckNear(controller.Current().high, 1.0, 1.e-4,
            "constraint release converges to the neutral state");

  controller.SetTarget({2.f, .2f, .9f,
                        std::numeric_limits<float>::quiet_NaN()});
  const PassiveConstraintGains bounded = controller.Process();
  Check(bounded.broadband <= 1.f && bounded.low <= 1.f &&
            bounded.middle <= bounded.low && bounded.high <= bounded.middle,
        "constraint gains remain passive and spectrally ordered");
}

} // namespace

int main() {
  TestLinearAndGeometricSegments();
  TestTrajectoryRetriggerAndSanitization();
  TestPassiveConstraintSmoothing();
  if (percussion_test::failures == 0)
    std::cout << "All percussion control tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
