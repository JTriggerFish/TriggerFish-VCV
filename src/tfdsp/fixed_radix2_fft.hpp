#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

namespace tfdsp {

template <typename Sample, std::size_t Size> class FixedRadix2Fft {
public:
  static_assert(Size >= 2 && (Size & (Size - 1)) == 0,
                "FFT size must be a power of two");
  using Spectrum = std::array<std::complex<Sample>, Size>;

  void Prepare() noexcept {
    constexpr double Pi = 3.141592653589793238462643383279502884;
    for (std::size_t index = 0; index < roots_.size(); ++index) {
      const double angle = -2.0 * Pi * static_cast<double>(index) / Size;
      roots_[index] = {static_cast<Sample>(std::cos(angle)),
                       static_cast<Sample>(std::sin(angle))};
    }
    for (std::size_t index = 0; index < Size; ++index) {
      std::size_t value = index;
      std::size_t reversed = 0;
      for (std::size_t bit = 1; bit < Size; bit <<= 1) {
        reversed = (reversed << 1) | (value & 1);
        value >>= 1;
      }
      bitReversed_[index] = reversed;
    }
  }

  void Transform(Spectrum &values, const bool inverse) const noexcept {
    for (std::size_t index = 0; index < Size; ++index)
      if (index < bitReversed_[index])
        std::swap(values[index], values[bitReversed_[index]]);
    for (std::size_t length = 2; length <= Size; length <<= 1) {
      const std::size_t half = length / 2;
      const std::size_t rootStep = Size / length;
      for (std::size_t start = 0; start < Size; start += length)
        for (std::size_t offset = 0; offset < half; ++offset) {
          const auto root = inverse ? std::conj(roots_[offset * rootStep])
                                    : roots_[offset * rootStep];
          const auto even = values[start + offset];
          const auto odd = values[start + offset + half] * root;
          values[start + offset] = even + odd;
          values[start + offset + half] = even - odd;
        }
    }
    if (inverse)
      for (auto &value : values)
        value /= static_cast<Sample>(Size);
  }

private:
  std::array<std::complex<Sample>, Size / 2> roots_{};
  std::array<std::size_t, Size> bitReversed_{};
};

} // namespace tfdsp
