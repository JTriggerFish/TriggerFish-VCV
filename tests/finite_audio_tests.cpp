#include "tfdsp/finite_audio.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace {

void Check(const bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

template <typename Sample> void TestClassification() {
  const Sample normal = Sample{.25};
  Check(tfdsp::FiniteNormalOrZero(normal) == normal,
        "positive normal value is unchanged");
  Check(tfdsp::FiniteNormalOrZero(-normal) == -normal,
        "negative normal value is unchanged");
  Check(tfdsp::FiniteNormalOrZero(std::numeric_limits<Sample>::min()) ==
            std::numeric_limits<Sample>::min(),
        "minimum normal value is unchanged");
  Check(tfdsp::FiniteNormalOrZero(std::numeric_limits<Sample>::denorm_min()) ==
            Sample{},
        "positive subnormal is flushed");
  Check(tfdsp::FiniteNormalOrZero(
            -std::numeric_limits<Sample>::denorm_min()) == Sample{},
        "negative subnormal is flushed");
  Check(tfdsp::FiniteNormalOrZero(std::numeric_limits<Sample>::infinity()) ==
            Sample{},
        "positive infinity is rejected");
  Check(tfdsp::FiniteNormalOrZero(
            -std::numeric_limits<Sample>::infinity()) == Sample{},
        "negative infinity is rejected");
  Check(tfdsp::FiniteNormalOrZero(
            std::numeric_limits<Sample>::quiet_NaN()) == Sample{},
        "NaN is rejected");
  Check(!std::signbit(tfdsp::FiniteNormalOrZero(-Sample{})),
        "negative zero becomes positive zero");
}

} // namespace

int main() {
  TestClassification<float>();
  TestClassification<double>();
  std::cout << "Finite-audio tests passed\n";
}
