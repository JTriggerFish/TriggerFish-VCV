// Frozen scalar implementation for optimisation-equivalence tests.
#pragma once

#include "tfdsp/percussion/deterministic_random.hpp"
#include "tfdsp/percussion/bounded_modal_motion.hpp"
#include "tfdsp/finite_audio.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace tfdsp::percussion {

// Bounded phase displacement, not an accumulated random frequency walk.
// Delta(displacement) rotates the stored complex state without changing its
// norm. At modest depth the coherent carrier survives alongside moving wings.
template <std::size_t Count> class ReferenceBoundedModalMotion {
public:
  void Prepare(float rate, const std::array<float, Count> &frequencies,
               const std::array<std::uint16_t, Count> &packets,
               std::size_t active, ModalMotionControls controls,
               std::uint32_t seed) {
    if (!std::isfinite(rate) || rate < 1000.f)
      throw std::invalid_argument("modal motion requires sample rate >= 1 kHz");
    active_ = std::min(active, Count);
    depth_ = std::clamp(tfdsp::FiniteNormalOrZero(controls.depthRadians), 0.f, 3.f);
    const float speed = std::clamp(
        tfdsp::FiniteNormalOrZero(controls.knotsPerSecond), .1f, 200.f);
    sharing_ = std::clamp(tfdsp::FiniteNormalOrZero(controls.packetSharing), 0.f, 1.f);
    increment_ = speed / rate;
    seed_ = seed;
    packetCount_ = 0;
    for (std::size_t i = 0; i < active_; ++i) {
      if (i == 0 || packets[i] != packets[i-1]) ++packetCount_;
      packetIndex_[i] = packetCount_-1;
      const float f = std::clamp(tfdsp::FiniteNormalOrZero(frequencies[i]), 1.f, .49f*rate);
      // Quintic derivative <= 1.875, target span <= 2, speed jitter <= 1.25.
      // Bound instantaneous deviation near DC/Nyquist; not a strict AA filter.
      depthByMode_[i] = std::min(depth_, 6.28318530718f *
          std::min(f-1.f, .49f*rate-f) / (4.6875f*speed));
    }
    Reset();
  }

  bool Enabled() const noexcept { return depth_ > 0.f; }

  void Reset() noexcept {
    random_.Seed(seed_);
    for (std::size_t i = 0; i < active_+packetCount_; ++i) {
      auto &t = trajectory_[i];
      t = {random_.Uniform(), random_.Bipolar(), random_.Bipolar(),
           increment_ * (.75f+.5f*random_.Uniform())};
      value_[i] = Value(t);
    }
    for (std::size_t i = 0; i < active_; ++i) previous_[i] = Displacement(i);
  }

  void BeginSample() noexcept {
    for (std::size_t i = 0; i < active_+packetCount_; ++i) {
      auto &t = trajectory_[i];
      t.phase += t.step;
      if (t.phase >= 1.f) {
        t.phase -= 1.f;
        t.from = t.to;
        t.to = random_.Bipolar();
        t.step = increment_ * (.75f+.5f*random_.Uniform());
      }
      value_[i] = Value(t);
    }
  }

  float NextAngle(std::size_t mode) noexcept {
    const float phase = Displacement(mode);
    const float delta = phase-previous_[mode];
    previous_[mode] = phase;
    return delta;
  }

  void RotateCoefficients(std::size_t mode, float &cosine, float &sine) noexcept {
    const float angle = NextAngle(mode), a2 = angle*angle;
    // Cayley rotation, unit norm algebraically; bounded small-angle tan series.
    const float tangent = std::abs(angle) > .3f ? std::tan(.5f*angle) :
        .5f*angle*(1.f+a2*(1.f/12.f+a2/120.f));
    const float inverse = 1.f/(1.f+tangent*tangent);
    const float c = (1.f-tangent*tangent)*inverse, s = 2.f*tangent*inverse;
    const float prior = cosine;
    cosine = c*prior-s*sine;
    sine = s*prior+c*sine;
  }

private:
  struct Trajectory { float phase{}, from{}, to{}, step{}; };
  static float Value(const Trajectory &t) noexcept {
    const float x = t.phase;
    return t.from + x*x*x*(10.f+x*(-15.f+6.f*x))*(t.to-t.from);
  }
  float Displacement(std::size_t mode) const noexcept {
    return depthByMode_[mode] * ((1.f-sharing_)*value_[mode] +
        sharing_*value_[active_+packetIndex_[mode]]);
  }
  std::array<Trajectory, 2*Count> trajectory_{};
  std::array<float, 2*Count> value_{};
  std::array<float, Count> previous_{}, depthByMode_{};
  std::array<std::size_t, Count> packetIndex_{};
  DeterministicRandom random_{};
  std::uint32_t seed_{};
  std::size_t active_{}, packetCount_{};
  float increment_{}, depth_{}, sharing_{};
};

} // namespace tfdsp::percussion
