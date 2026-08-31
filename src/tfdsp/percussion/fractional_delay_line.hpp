#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace tfdsp::percussion {

// Circular delay storage with an efficient high-quality read for moving taps.
// Prepare() allocates and warms the shared interpolation table; audio processing
// performs no allocation or transcendental operations.
class FractionalDelayLine {
public:
  static constexpr std::size_t InterpolationTaps = 12;
  static constexpr std::size_t InterpolationPhases = 2048;
  static constexpr float MinimumSincDelaySamples = 6.f;

  void Prepare(const float maximumDelaySamples) {
    if (!std::isfinite(maximumDelaySamples) ||
        maximumDelaySamples < MinimumSincDelaySamples)
      throw std::invalid_argument("fractional delay capacity is too small");
    maximumDelaySamples_ = maximumDelaySamples;
    const auto samples = static_cast<std::size_t>(
        std::ceil(maximumDelaySamples)) + InterpolationTaps + 2;
    buffer_.assign(samples, 0.f);
    (void)CoefficientTable();
    Reset();
  }

  void Reset() noexcept {
    std::fill(buffer_.begin(), buffer_.end(), 0.f);
    writeIndex_ = 0;
  }

  void Push(const float input) noexcept {
    if (buffer_.empty())
      return;
    buffer_[writeIndex_] = tfdsp::FiniteNormalOrZero(input);
    if (++writeIndex_ == buffer_.size())
      writeIndex_ = 0;
  }

  float ReadInteger(const std::size_t delaySamples) const noexcept {
    if (buffer_.empty())
      return 0.f;
    const auto bounded = std::clamp<std::size_t>(
        delaySamples, 1, static_cast<std::size_t>(maximumDelaySamples_));
    return ReadUnchecked(bounded);
  }

  float ReadSinc(float delaySamples) const noexcept {
    if (buffer_.empty())
      return 0.f;
    delaySamples = std::clamp(delaySamples, MinimumSincDelaySamples,
                              maximumDelaySamples_);
    const auto integer = static_cast<std::size_t>(std::floor(delaySamples));
    const float fraction = delaySamples - static_cast<float>(integer);
    if (fraction < 1.e-7f)
      return ReadInteger(integer);
    const float tablePosition = fraction * InterpolationPhases;
    const auto phase = static_cast<std::size_t>(tablePosition);
    const float blend = tablePosition - static_cast<float>(phase);
    const auto &first = CoefficientTable()[phase];
    const auto &second = CoefficientTable()[phase + 1];

    float output = 0.f;
    for (std::size_t tap = 0; tap < InterpolationTaps; ++tap) {
      const auto delay = integer + tap - 5;
      const float coefficient = first[tap] + blend * (second[tap] - first[tap]);
      output += coefficient * ReadUnchecked(delay);
    }
    return output;
  }

  float MaximumDelaySamples() const noexcept { return maximumDelaySamples_; }

private:
  using Coefficients = std::array<float, InterpolationTaps>;
  using Table = std::array<Coefficients, InterpolationPhases + 1>;
  static constexpr double Pi = 3.1415926535897932384626433832795;

  static double Sinc(const double value) noexcept {
    if (std::abs(value) < 1.e-12)
      return 1.0;
    const double radians = Pi * value;
    return std::sin(radians) / radians;
  }

  static Table MakeCoefficientTable() {
    Table table{};
    for (std::size_t phase = 0; phase <= InterpolationPhases; ++phase) {
      const double fraction = static_cast<double>(phase) / InterpolationPhases;
      double sum = 0.0;
      for (std::size_t tap = 0; tap < InterpolationTaps; ++tap) {
        const double offset = static_cast<double>(tap) - 5.0 - fraction;
        const double coefficient = Sinc(offset) * Sinc(offset / 6.0);
        table[phase][tap] = static_cast<float>(coefficient);
        sum += coefficient;
      }
      for (auto &coefficient : table[phase])
        coefficient = static_cast<float>(coefficient / sum);
    }
    return table;
  }

  static const Table &CoefficientTable() {
    static const Table table = MakeCoefficientTable();
    return table;
  }

  float ReadUnchecked(const std::size_t delaySamples) const noexcept {
    const auto wrappedDelay = delaySamples % buffer_.size();
    const auto index = writeIndex_ >= wrappedDelay
                           ? writeIndex_ - wrappedDelay
                           : writeIndex_ + buffer_.size() - wrappedDelay;
    return buffer_[index];
  }

  std::vector<float> buffer_{};
  std::size_t writeIndex_{};
  float maximumDelaySamples_{MinimumSincDelaySamples};
};

} // namespace tfdsp::percussion
