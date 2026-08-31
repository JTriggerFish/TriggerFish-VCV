#pragma once

#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace tfdsp::percussion {

enum class TrajectoryCurve { Linear, Geometric };

struct TrajectorySegment {
  float targetValue{};
  float durationSeconds{};
  TrajectoryCurve curve{TrajectoryCurve::Linear};
};

// Fixed-capacity breakpoint trajectory for fitted audio-rate controls. Segment
// setup is off the sample loop; processing is allocation-free and retriggering
// may begin either from an explicit value or the currently rendered value.
template <std::size_t MaximumSegments> class BreakpointTrajectory {
public:
  using Segments = std::array<TrajectorySegment, MaximumSegments>;

  void Prepare(const float sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("trajectory sample rate must be positive");
    sampleRate_ = sampleRate;
    Reset();
  }

  void Reset(const float value = 0.f) noexcept {
    value_ = tfdsp::FiniteNormalOrZero(value);
    segmentStart_ = value_;
    segmentIndex_ = segmentSample_ = segmentSamples_ = segmentCount_ = 0;
    linearStep_ = 0.f;
    geometricRatio_ = 1.f;
    geometric_ = false;
  }

  void Start(const float initialValue, const Segments &segments,
             const std::size_t segmentCount) noexcept {
    value_ = tfdsp::FiniteNormalOrZero(initialValue);
    Load(segments, segmentCount);
  }

  void RetriggerFromCurrent(const Segments &segments,
                            const std::size_t segmentCount) noexcept {
    Load(segments, segmentCount);
  }

  float Process() noexcept {
    if (!Active())
      return value_;
    ++segmentSample_;
    if (geometric_)
      value_ *= geometricRatio_;
    else
      value_ = static_cast<float>(
          static_cast<double>(segmentStart_) +
          static_cast<double>(linearStep_) * segmentSample_);
    if (segmentSample_ >= segmentSamples_) {
      value_ = segments_[segmentIndex_].targetValue;
      ++segmentIndex_;
      BeginSegment();
    }
    value_ = tfdsp::FiniteNormalOrZero(value_);
    return value_;
  }

  float Value() const noexcept { return value_; }
  bool Active() const noexcept { return segmentIndex_ < segmentCount_; }

private:
  static TrajectorySegment Sanitize(TrajectorySegment segment) noexcept {
    segment.targetValue = tfdsp::FiniteNormalOrZero(segment.targetValue);
    segment.durationSeconds = std::clamp(
        std::isfinite(segment.durationSeconds) ? segment.durationSeconds : 0.f,
        0.f, 60.f);
    return segment;
  }

  void Load(const Segments &segments, const std::size_t count) noexcept {
    segmentCount_ = std::min(count, MaximumSegments);
    for (std::size_t index = 0; index < segmentCount_; ++index)
      segments_[index] = Sanitize(segments[index]);
    segmentIndex_ = 0;
    BeginSegment();
  }

  void BeginSegment() noexcept {
    while (segmentIndex_ < segmentCount_ &&
           segments_[segmentIndex_].durationSeconds <= 0.f) {
      value_ = segments_[segmentIndex_].targetValue;
      ++segmentIndex_;
    }
    if (!Active())
      return;
    segmentStart_ = value_;
    segmentSample_ = 0;
    segmentSamples_ = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::lround(
               segments_[segmentIndex_].durationSeconds * sampleRate_)));
    const float target = segments_[segmentIndex_].targetValue;
    const double step =
        (static_cast<double>(target) - segmentStart_) / segmentSamples_;
    linearStep_ = static_cast<float>(std::clamp(
        step, -static_cast<double>(std::numeric_limits<float>::max()),
        static_cast<double>(std::numeric_limits<float>::max())));
    geometric_ = segments_[segmentIndex_].curve ==
                     TrajectoryCurve::Geometric &&
                 segmentStart_ != 0.f && target != 0.f &&
                 std::signbit(segmentStart_) == std::signbit(target);
    geometricRatio_ = 1.f;
    if (geometric_) {
      const double ratio = std::exp(
          (std::log(std::abs(static_cast<double>(target))) -
           std::log(std::abs(static_cast<double>(segmentStart_)))) /
          segmentSamples_);
      geometric_ = std::isfinite(ratio) &&
                   ratio >= std::numeric_limits<float>::min() &&
                   ratio <= std::numeric_limits<float>::max();
      if (geometric_)
        geometricRatio_ = static_cast<float>(ratio);
    }
  }

  Segments segments_{};
  float sampleRate_{48000.f};
  float value_{};
  float segmentStart_{};
  float linearStep_{};
  float geometricRatio_{1.f};
  std::size_t segmentIndex_{};
  std::size_t segmentSample_{};
  std::size_t segmentSamples_{};
  std::size_t segmentCount_{};
  bool geometric_{};
};

} // namespace tfdsp::percussion
