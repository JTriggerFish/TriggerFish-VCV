#pragma once

#include "finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace tfdsp {

// Fixed-width delay bank with interleaved histories. Tap gathering is scalar
// on baseline SSE/WebAssembly; interpolation arithmetic and writes are batched.
template <std::size_t LineCount> class CubicFractionalDelayBank {
public:
  using Frame = std::array<float, LineCount>;

  void Prepare(const std::size_t capacity) {
    if (capacity < 8)
      throw std::invalid_argument("fractional delay capacity is too small");
    capacity_ = capacity;
    buffer_.assign(capacity * LineCount, 0.f);
    writeIndex_ = 0;
  }

  void Reset() noexcept {
    std::fill(buffer_.begin(), buffer_.end(), 0.f);
    writeIndex_ = 0;
  }

  Frame Read(const Frame &delaySamples) const noexcept {
    if (buffer_.empty())
      return {};
    std::array<Frame, 4> taps{};
    Frame fraction{};
    for (std::size_t line = 0; line < LineCount; ++line)
      Gather(line, delaySamples[line], fraction[line], taps);

    Frame output{};
    for (std::size_t line = 0; line < LineCount; ++line) {
      const float mu = fraction[line];
      const std::array<float, 4> coefficients{
          -mu * (mu - 1.f) * (mu - 2.f) / 6.f,
          (mu + 1.f) * (mu - 1.f) * (mu - 2.f) / 2.f,
          -(mu + 1.f) * mu * (mu - 2.f) / 2.f,
          (mu + 1.f) * mu * (mu - 1.f) / 6.f};
      float value = 0.f;
      for (std::size_t tap = 0; tap < taps.size(); ++tap)
        value += coefficients[tap] * taps[tap][line];
      output[line] = value;
    }
    return output;
  }

  void Push(const Frame &input) noexcept {
    if (buffer_.empty())
      return;
    const std::size_t base = writeIndex_ * LineCount;
    for (std::size_t line = 0; line < LineCount; ++line)
      buffer_[base + line] = FiniteNormalOrZero(input[line]);
    if (++writeIndex_ == capacity_)
      writeIndex_ = 0;
  }

private:
  float ReadInteger(const std::size_t line,
                    const std::size_t distance) const noexcept {
    const std::size_t index = writeIndex_ >= distance
                                  ? writeIndex_ - distance
                                  : writeIndex_ + capacity_ - distance;
    return buffer_[index * LineCount + line];
  }

  void Gather(const std::size_t line, float delaySamples, float &fraction,
              std::array<Frame, 4> &taps) const noexcept {
    if (!std::isfinite(delaySamples))
      delaySamples = 2.f;
    delaySamples = std::clamp(delaySamples, 2.f,
                              static_cast<float>(capacity_ - 3));
    const auto integer = static_cast<std::size_t>(std::floor(delaySamples));
    fraction = delaySamples - static_cast<float>(integer);
    taps[0][line] = ReadInteger(line, integer - 1);
    taps[1][line] = ReadInteger(line, integer);
    taps[2][line] = ReadInteger(line, integer + 1);
    taps[3][line] = ReadInteger(line, integer + 2);
  }

  std::vector<float> buffer_{};
  std::size_t capacity_{};
  std::size_t writeIndex_{};
};

} // namespace tfdsp
