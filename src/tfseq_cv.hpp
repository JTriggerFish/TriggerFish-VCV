#pragma once

#include "tfseq.hpp"

#include <algorithm>
#include <cmath>

namespace tfseq {

// CV knots run in score time independently of note density. The current
// interpolation segment is cached, so only knot crossings scan the lane.
class CvLanePlayer {
public:
  void setEvent(const RuntimeEvent &event, std::size_t lane) noexcept {
    lane_ = lane;
    sequence_ = event.cvSequence;
    interpolation_ = event.cvInterpolation[lane];
    power_ = event.cvPower[lane];
    originBeat_ = event.cvOriginBeat;
    cycle_ = event.cvCycle;
    seed_ = event.cvSeed;
    if (sequence_ &&
        sequence_->alignment[static_cast<std::size_t>(CvCursorLane(lane))] ==
            LaneAlignment::Free &&
        !sequence_->cv[lane].empty()) {
      refresh(event.beat);
    } else {
      sequence_ = nullptr;
      value_ = from_ = target_ = event.cvValue[lane];
      beginBeat_ = endBeat_ = event.beat;
    }
  }

  // The last value remains usable while waiting for the new program's event.
  void detach() noexcept {
    sequence_ = nullptr;
    from_ = target_ = value_;
    beginBeat_ = endBeat_ = 0.0;
  }

  float process(double beat) noexcept {
    if (sequence_ && (beat >= endBeat_ || beat < beginBeat_))
      refresh(beat);
    if (interpolation_ == CvInterpolation::Step || endBeat_ <= beginBeat_) {
      value_ = target_;
    } else {
      double amount = (beat - beginBeat_) / (endBeat_ - beginBeat_);
      amount = std::clamp(amount, 0.0, 1.0);
      if (interpolation_ == CvInterpolation::Smooth)
        amount = amount * amount * (3.0 - 2.0 * amount);
      else if (interpolation_ == CvInterpolation::Power)
        amount = std::pow(amount, power_);
      value_ = static_cast<float>(
          from_ + (static_cast<double>(target_) - from_) * amount);
    }
    return value_;
  }

private:
  void refresh(double beat) noexcept;
  const Sequence *sequence_ = nullptr;
  std::size_t lane_ = 0;
  double originBeat_ = 0.0;
  std::uint64_t cycle_ = 0;
  std::uint64_t seed_ = 1;
  float value_ = 0.f;
  float from_ = 0.f;
  float target_ = 0.f;
  double beginBeat_ = 0.0;
  double endBeat_ = 0.0;
  float power_ = 1.f;
  CvInterpolation interpolation_ = CvInterpolation::Step;
};

} // namespace tfseq
