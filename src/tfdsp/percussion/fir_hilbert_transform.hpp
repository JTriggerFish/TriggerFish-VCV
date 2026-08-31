#pragma once

#include <array>
#include <cmath>
#include <cstddef>

namespace tfdsp::percussion {

struct AnalyticSample {
  float real{};
  float quadrature{};
};

// Linear-phase Blackman-windowed Hilbert transformer. Antisymmetry and the
// zero-valued even taps reduce the 255-tap reference design to 64 multiplies.
class FirHilbertTransform {
public:
  static constexpr std::size_t TapCount = 255;
  static constexpr std::size_t LatencySamples = (TapCount - 1) / 2;

  FirHilbertTransform() { MakeCoefficients(); }

  void Reset() noexcept {
    samples_.fill(0.f);
    writeIndex_ = 0;
  }

  AnalyticSample Process(float input) noexcept {
    if (!std::isfinite(input))
      input = 0.f;
    samples_[writeIndex_] = input;
    if (++writeIndex_ == TapCount)
      writeIndex_ = 0;

    float quadrature = 0.f;
    for (std::size_t pair = 0; pair < PairCount; ++pair) {
      const std::size_t offset = 2 * pair + 1;
      quadrature += coefficients_[pair] *
                    (Read(LatencySamples + offset) -
                     Read(LatencySamples - offset));
    }
    return {Read(LatencySamples), quadrature};
  }

private:
  static constexpr std::size_t PairCount = (LatencySamples + 1) / 2;
  static constexpr double Pi = 3.1415926535897932384626433832795;

  void MakeCoefficients() noexcept {
    constexpr double denominator = TapCount - 1;
    for (std::size_t pair = 0; pair < PairCount; ++pair) {
      const std::size_t offset = 2 * pair + 1;
      const double tap = static_cast<double>(LatencySamples + offset);
      const double window = .42 - .5 * std::cos(2.0 * Pi * tap / denominator) +
                            .08 * std::cos(4.0 * Pi * tap / denominator);
      coefficients_[pair] =
          static_cast<float>(window * 2.0 /
                             (Pi * static_cast<double>(offset)));
    }
  }

  float Read(const std::size_t delay) const noexcept {
    const std::size_t latest = writeIndex_ == 0 ? TapCount - 1 : writeIndex_ - 1;
    const std::size_t wrapped = delay % TapCount;
    const std::size_t index = latest >= wrapped
                                  ? latest - wrapped
                                  : latest + TapCount - wrapped;
    return samples_[index];
  }

  std::array<float, TapCount> samples_{};
  std::array<float, PairCount> coefficients_{};
  std::size_t writeIndex_{};
};

} // namespace tfdsp::percussion
