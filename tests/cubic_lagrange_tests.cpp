#include "percussion_test_support.hpp"

#include "tfdsp/cubic_fractional_delay.hpp"
#include "tfdsp/cubic_lagrange_interpolator.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

namespace {

using percussion_test::Check;
using percussion_test::CheckNear;

void TestCubicPolynomialIsExact() {
  for (float delay = 2.f; delay <= 30.f; delay += .125f) {
    const float output = tfdsp::ReadCubicLagrange(
        delay, 64, [](const std::size_t distance) {
          const float time = -static_cast<float>(distance);
          return .25f * time * time * time - .5f * time * time + 2.f * time;
        });
    const float expected =
        -.25f * delay * delay * delay - .5f * delay * delay - 2.f * delay;
    CheckNear(output, expected, 2.e-3,
              "cubic Lagrange exactly reconstructs cubic data");
  }
}

void TestIntegerBoundaryContinuity() {
  constexpr float Omega = .21f;
  for (std::size_t integer = 3; integer < 30; ++integer) {
    const auto read = [=](const std::size_t distance) {
      return std::sin(-Omega * static_cast<float>(distance));
    };
    const float left = tfdsp::ReadCubicLagrange(integer - 1.e-4f, 64, read);
    const float right = tfdsp::ReadCubicLagrange(integer + 1.e-4f, 64, read);
    Check(std::abs(left - right) < 5.e-5f,
          "cubic Lagrange is continuous across integer delays");
  }
}

void TestInvalidInputsAreSafe() {
  const float nan = std::numeric_limits<float>::quiet_NaN();
  const float fallback = tfdsp::ReadCubicLagrange(
      nan, 64, [](const std::size_t distance) {
        return static_cast<float>(distance);
      });
  CheckNear(fallback, 2.0, 1.e-7,
            "non-finite cubic delay uses the minimum causal tap");
  Check(tfdsp::ReadCubicLagrange(
            2.f, 0, [](const std::size_t) { return 1.f; }) == 0.f,
        "undersized cubic storage cannot underflow its capacity");

  tfdsp::CubicFractionalDelay delay;
  delay.Prepare(16);
  delay.Push(1.f);
  Check(std::isfinite(delay.Read(nan)),
        "fractional-delay wrapper sanitizes a non-finite read position");

  tfdsp::LiveFractionalDelay live;
  live.Prepare(16);
  Check(live.Read(nan, nan) == 0.f,
        "live fractional delay sanitizes position and current sample");
}

} // namespace

int main() {
  TestCubicPolynomialIsExact();
  TestIntegerBoundaryContinuity();
  TestInvalidInputsAreSafe();
  if (percussion_test::failures == 0)
    std::cout << "All cubic Lagrange tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
