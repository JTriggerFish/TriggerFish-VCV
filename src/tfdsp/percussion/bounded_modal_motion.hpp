#pragma once

#include "deterministic_random.hpp"
#include "tfdsp/finite_audio.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace tfdsp::percussion {

struct ModalMotionControls {
  float depthRadians{};
  float knotsPerSecond{40.f};
  float packetSharing{.5f};
};

// Bounded phase displacement, not an accumulated random frequency walk.
// Delta(displacement) rotates the stored complex state without changing its
// norm. At modest depth the coherent carrier survives alongside moving wings.
template <std::size_t Count> class BoundedModalMotion {
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
    // Same derivative bound as the frequency guard below; leave rounding
    // margin below the existing .3-radian scalar tan() threshold.
    smallAngles_ = 4.6875f * speed * depth_ / rate < .29f;
    seed_ = seed;
    packetCount_ = 0;
    for (std::size_t i = 0; i < active_; ++i) {
      if (i == 0 || packets[i] != packets[i-1]) packetBegin_[packetCount_++] = i;
      packetIndex_[i] = packetCount_-1;
      const float f = std::clamp(tfdsp::FiniteNormalOrZero(frequencies[i]), 1.f, .49f*rate);
      // Quintic derivative <= 1.875, target span <= 2, speed jitter <= 1.25.
      // Bound instantaneous deviation near DC/Nyquist; not a strict AA filter.
      depthByMode_[i] = std::min(depth_, 6.28318530718f *
          std::min(f-1.f, .49f*rate-f) / (4.6875f*speed));
    }
    packetBegin_[packetCount_] = active_;
    Reset();
  }

  bool Enabled() const noexcept { return depth_ > 0.f; }

  void Reset() noexcept {
    random_.Seed(seed_);
    for (std::size_t i = 0; i < active_+packetCount_; ++i) {
      phase_[i] = random_.Uniform();
      from_[i] = random_.Bipolar();
      to_[i] = random_.Bipolar();
      step_[i] = increment_ * (.75f+.5f*random_.Uniform());
      value_[i] = Value(i);
    }
    for (std::size_t i = 0; i < active_; ++i) previous_[i] = Displacement(i);
  }

  void BeginSample() noexcept {
    const auto count = active_ + packetCount_;
    // Separate the regular SIMD work from infrequent random knot transitions.
    for (std::size_t i = 0; i < count; ++i) phase_[i] += step_[i];
    for (std::size_t i = 0; i < count; ++i) {
      if (phase_[i] >= 1.f) {
        phase_[i] -= 1.f;
        from_[i] = to_[i];
        to_[i] = random_.Bipolar();
        step_[i] = increment_ * (.75f+.5f*random_.Uniform());
      }
    }
    for (std::size_t i = 0; i < count; ++i) value_[i] = Value(i);
  }

  // Hoist the scalar packet lookup and rare large-angle path out of the modal
  // propagation loop. Contiguous rotations leave that loop SIMD-friendly.
  void PrepareRotations() noexcept {
    for (std::size_t packet = 0; packet < packetCount_; ++packet) {
      const float shared = sharing_ * value_[active_+packet];
      for (std::size_t i = packetBegin_[packet]; i < packetBegin_[packet+1]; ++i) {
        const float phase = depthByMode_[i] * ((1.f-sharing_)*value_[i]+shared);
        angle_[i] = phase-previous_[i];
        previous_[i] = phase;
      }
    }
    if (smallAngles_) MakeRotations<true>();
    else MakeRotations<false>();
  }

  void RotatePrepared(std::size_t mode, float &cosine, float &sine) const noexcept {
    const float prior = cosine;
    cosine = rotationCos_[mode]*prior-rotationSin_[mode]*sine;
    sine = rotationSin_[mode]*prior+rotationCos_[mode]*sine;
  }

  float NextAngle(std::size_t mode) noexcept {
    const float phase = Displacement(mode);
    const float delta = phase-previous_[mode];
    previous_[mode] = phase;
    return delta;
  }

private:
  // Unchanged Cayley rotation: unit norm algebraically, not a cheaper
  // first-order approximation that would inject energy over a long tail.
  template <bool Small> void MakeRotations() noexcept {
    for (std::size_t i = 0; i < active_; ++i) {
      const float a = angle_[i], a2 = a*a;
      float tangent = .5f*a*(1.f+a2*(1.f/12.f+a2/120.f));
      if constexpr (!Small) {
        if (std::abs(a) > .3f) tangent = std::tan(.5f*a);
      }
      const float inverse = 1.f/(1.f+tangent*tangent);
      rotationCos_[i] = (1.f-tangent*tangent)*inverse;
      rotationSin_[i] = 2.f*tangent*inverse;
    }
  }
  float Value(std::size_t i) const noexcept {
    const float x = phase_[i];
    return from_[i] + x*x*x*(10.f+x*(-15.f+6.f*x))*(to_[i]-from_[i]);
  }
  float Displacement(std::size_t mode) const noexcept {
    return depthByMode_[mode] * ((1.f-sharing_)*value_[mode] +
        sharing_*value_[active_+packetIndex_[mode]]);
  }
  std::array<float, 2*Count> phase_{}, from_{}, to_{}, step_{};
  std::array<float, 2*Count> value_{};
  std::array<float, Count> previous_{}, depthByMode_{};
  std::array<float, Count> angle_{}, rotationCos_{}, rotationSin_{};
  std::array<std::size_t, Count> packetIndex_{};
  std::array<std::size_t, Count+1> packetBegin_{};
  DeterministicRandom random_{};
  std::uint32_t seed_{};
  std::size_t active_{}, packetCount_{};
  float increment_{}, depth_{}, sharing_{};
  bool smallAngles_{};
};

} // namespace tfdsp::percussion
