#include "percussion_test_support.hpp"
#include "tfdsp/percussion/modal_spectral_diffusion.hpp"
#include <numeric>

using percussion_test::Check;
using percussion_test::CheckNear;
using tfdsp::percussion::ModalSpectralDiffusion;

namespace {
template <std::size_t N> double Sum(const std::array<float, N> &x) {
  return std::accumulate(x.begin(), x.end(), 0.0);
}

void ConservationAndEquilibrium() {
  ModalSpectralDiffusion<5> diffusion;
  const std::array<float, 5> frequencies{40, 170, 360, 800, 950};
  const std::array<float, 5> weights{.2f, .2f, .2f, .2f, .2f};
  diffusion.Prepare(frequencies, weights, 5, 2000);
  std::array<float, 5> energy{}, result{};
  for (std::size_t i = 0; i < 5; ++i) energy[i] = float(diffusion.CellWidth(i));
  diffusion.Process(energy, result, 1, 1000, 1);
  for (std::size_t i = 0; i < 5; ++i)
    CheckNear(result[i], energy[i], 1.e-7, "flat density is equilibrium on unequal cells");
  for (double scale : {1.e-9, 1.0, 1.e6}) {
    energy = {float(scale), 0, float(.2*scale), 0, 0};
    for (double step : {0., 1.e-6, .01, 1.e6}) {
      diffusion.Process(energy, result, 1, step, 1);
      CheckNear(Sum(result) / Sum(energy), 1, 2.e-7, "stiff step conserves energy");
      for (float x : result) Check(std::isfinite(x) && x >= 0, "positive finite energy");
    }
  }
  energy.fill(0);
  diffusion.Process(energy, result, 1, 10, 1);
  Check(Sum(result) == 0, "silence cannot self-excite");
}

void DirectionAndStrength() {
  ModalSpectralDiffusion<2> diffusion;
  diffusion.Prepare({250, 750}, {.5f, .5f}, 2, 2000);
  std::array<float, 2> low{1, 0}, high{0, 1}, loud{4, 0}, a{}, b{}, c{};
  diffusion.Process(low, a, 1, .001, 1);
  diffusion.Process(high, b, 1, .001, 1);
  diffusion.Process(loud, c, 1, .001, 1);
  Check(a[1] > 0 && b[0] > 0, "gradient permits both transfer directions");
  Check(a[0] > .9f, "initial transfer does not abruptly remove the low body");
  Check(c[1]/4 > a[1]*10, "quadratic conductivity strengthens with strike energy");
  diffusion.Process(low, a, 1, 0, 1);
  Check(a == low, "zero strength bypass is exact");
  for (double exponent : {.25, .5, .75}) {
    std::array<float, 2> quiet{.01f, 0}, quietNext{};
    diffusion.Process(low, a, 1, .001, exponent);
    diffusion.Process(quiet, quietNext, 1, .001, exponent);
    Check(quietNext[1]/.01 < a[1]*.2,
        "intermediate nonlinearities have no constant low-energy leak");
    CheckNear(Sum(quietNext), .01, 1.e-8, "nonlinear exponent remains passive");
    Check(quietNext[0] >= 0 && quietNext[1] >= 0, "fractional exponent is positive");
  }
}

void CoincidentHandles() {
  ModalSpectralDiffusion<3> split;
  ModalSpectralDiffusion<2> original;
  split.Prepare({250, 250, 750}, {.2f, .3f, .5f}, 3, 2000);
  original.Prepare({250, 750}, {.5f, .5f}, 2, 2000);
  std::array<float, 3> a{};
  std::array<float, 2> b{};
  for (const std::array<float, 3> start : {
      std::array<float, 3>{.4f, .6f, 0}, {0, 0, 1}}) {
    split.Process(start, a, 1, .1, 1);
    original.Process({start[0]+start[1], start[2]}, b, 1, .1, 1);
    CheckNear(a[0]+a[1], b[0], 1.e-7, "duplicate handles do not change transfer");
    CheckNear(a[2], b[1], 1.e-7, "duplicate handles preserve destination energy");
    CheckNear(a[0]/a[1], 2./3., 1.e-6, "silent duplicate cell uses input weights");
  }
}

template <std::size_t N> std::array<double, 4> RefinedSpectrum() {
  ModalSpectralDiffusion<N> diffusion;
  std::array<float, N> frequency{}, weight{}, energy{}, next{};
  for (std::size_t i = 0; i < N; ++i) {
    const double x = (double(i)+.5)/N;
    frequency[i] = float(x*1000);
    weight[i] = 1.f/N;
    energy[i] = float((.1+std::exp(-6*x))/N);
  }
  diffusion.Prepare(frequency, weight, N, 2000);
  for (int step = 0; step < 100; ++step) {
    diffusion.Process(energy, next, 1, .0001, 1);
    energy = next;
  }
  std::array<double, 4> bands{};
  for (std::size_t i = 0; i < N; ++i) bands[4*i/N] += energy[i];
  return bands;
}

void GridConvergence() {
  const auto coarse = RefinedSpectrum<16>();
  const auto fine = RefinedSpectrum<32>();
  const auto reference = RefinedSpectrum<128>();
  double coarseError = 0, fineError = 0;
  for (std::size_t i = 0; i < 4; ++i) {
    coarseError += std::abs(coarse[i]-reference[i]);
    fineError += std::abs(fine[i]-reference[i]);
  }
  Check(fineError < .4*coarseError, "frequency-cell refinement converges, not serial slowdown");
}

std::array<float, 3> TimedSpectrum(const int steps) {
  ModalSpectralDiffusion<3> diffusion;
  diffusion.Prepare({100, 400, 800}, {.3f, .3f, .4f}, 3, 2000);
  std::array<float, 3> energy{.25f, .1f, 0}, next{};
  for (int i=0; i<steps; ++i) {
    diffusion.Process(energy, next, 1, .1/steps, 1);
    energy = next;
  }
  return energy;
}

void TimeConvergence() {
  const auto coarse = TimedSpectrum(100);
  const auto fine = TimedSpectrum(200);
  const auto reference = TimedSpectrum(1600);
  double a = 0, b = 0;
  for (std::size_t i=0; i<3; ++i) {
    a += std::abs(coarse[i]-reference[i]);
    b += std::abs(fine[i]-reference[i]);
  }
  Check(b < .65*a, "semi-implicit time refinement converges");
}

void SampleRateIndependentGeometry() {
  ModalSpectralDiffusion<3> low, high;
  low.Prepare({100, 400, 800}, {.3f, .3f, .4f}, 3, 48000);
  high.Prepare({100, 400, 800}, {.3f, .3f, .4f}, 3, 96000);
  std::array<float, 3> a{}, b{};
  low.Process({1, .2f, 0}, a, 1, .001, 1);
  high.Process({1, .2f, 0}, b, 1, .001, 1);
  Check(a == b, "higher Nyquist must not widen the last represented frequency cell");
}
}

int main() {
  ConservationAndEquilibrium();
  DirectionAndStrength();
  CoincidentHandles();
  GridConvergence();
  TimeConvergence();
  SampleRateIndependentGeometry();
  return percussion_test::failures ? 1 : 0;
}
