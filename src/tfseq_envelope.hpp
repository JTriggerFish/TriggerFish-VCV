#pragma once

#include "tfseq.hpp"

#include <algorithm>
#include <cmath>

namespace tfseq {

inline bool CvEnvelopeTriggers(EventKind kind) noexcept {
  return kind == EventKind::Attack;
}

inline float CvEnvelopePeak(const CvEnvelopeSpec &spec, float velocity,
                            bool accented) noexcept {
  float peak = spec.depth;
  if (spec.followVelocity)
    peak *= velocity;
  if (accented)
    peak *= spec.accentMultiplier;
  return peak;
}

inline float CvEnvelopeOutput(float base, float envelope,
                              const CvEnvelopeSpec &spec) noexcept {
  if (!spec.enabled)
    return base;
  return spec.composition == CvEnvelopeComposition::Add ? base + envelope
                                                        : envelope;
}

// A sample-by-sample, retrigger-safe CV envelope. Segment endpoints are exact;
// curve changes affect their shape while preserving their duration.
class CvEnvelopeEngine {
public:
  void reset() noexcept {
    stage_ = Stage::Idle;
    value_ = 0.f;
    start_ = 0.f;
    target_ = 0.f;
    peak_ = 0.f;
    progress_ = 0.0;
    lastGate_ = false;
    modeKnown_ = false;
  }

  float process(bool gateHigh, bool trigger, float peak,
                const CvEnvelopeSpec &spec, double deltaSeconds,
                double deltaBeats) noexcept {
    if (!spec.enabled) {
      reset();
      return 0.f;
    }
    if (!modeKnown_ || spec.mode != mode_) {
      reset();
      mode_ = spec.mode;
      modeKnown_ = true;
    }

    const bool gateRose = gateHigh && !lastGate_;
    const bool gateFell = !gateHigh && lastGate_;
    if (spec.mode == CvEnvelopeMode::Ad) {
      if (trigger)
        beginAttack(peak, spec);
    } else if (spec.mode == CvEnvelopeMode::Ar) {
      if (gateRose)
        beginAttack(peak, spec);
      else if (gateFell)
        begin(Stage::Release, 0.f);
    } else {
      if ((trigger && gateHigh) || gateRose)
        beginAttack(peak, spec);
      else if (gateFell)
        begin(Stage::Release, 0.f);
    }
    lastGate_ = gateHigh;

    double remainingSeconds = std::max(0.0, deltaSeconds);
    double remainingBeats = std::max(0.0, deltaBeats);
    for (int transition = 0; transition < 4; ++transition) {
      const auto duration = stageDuration(spec);
      if (stage_ == Stage::Idle || stage_ == Stage::Hold ||
          stage_ == Stage::Sustain)
        break;
      const double delta = duration.unit == CvEnvelopeTimeUnit::Seconds
                               ? remainingSeconds
                               : remainingBeats;
      if (duration.value > 0.0) {
        const double required = duration.value * (1.0 - progress_);
        const double consumedFraction =
            delta > 0.0 ? std::min(1.0, required / delta) : 1.0;
        progress_ = std::min(1.0, progress_ + delta / duration.value);
        const float shaped = shape(static_cast<float>(progress_), spec.curve);
        value_ = start_ + (target_ - start_) * shaped;
        if (progress_ >= 1.0) {
          remainingSeconds *= 1.0 - consumedFraction;
          remainingBeats *= 1.0 - consumedFraction;
        }
      } else {
        progress_ = 1.0;
        value_ = target_;
      }
      if (progress_ < 1.0)
        break;
      value_ = target_;
      advance(spec);
    }
    return value_;
  }

  float value() const noexcept { return value_; }

private:
  enum class Stage { Idle, Attack, Decay, Hold, Sustain, Release };

  static float shape(float phase, float curve) noexcept {
    phase = std::clamp(phase, 0.f, 1.f);
    curve = std::clamp(curve, -1.f, 1.f);
    const float amount = curve <= 0.f ? 3.f * (curve + 1.f) : 3.f + 5.f * curve;
    if (amount < 1e-5f)
      return phase;
    return std::expm1(-amount * phase) / std::expm1(-amount);
  }

  void begin(Stage stage, float target) noexcept {
    stage_ = stage;
    start_ = value_;
    target_ = target;
    progress_ = 0.0;
  }

  void beginAttack(float peak, const CvEnvelopeSpec &spec) noexcept {
    peak_ = peak;
    begin(Stage::Attack, peak_);
    if (spec.attack.value == 0.0)
      value_ = peak_;
  }

  CvEnvelopeTime stageDuration(const CvEnvelopeSpec &spec) const noexcept {
    switch (stage_) {
    case Stage::Attack:
      return spec.attack;
    case Stage::Decay:
      return spec.decay;
    case Stage::Release:
      return spec.release;
    default:
      return {};
    }
  }

  void advance(const CvEnvelopeSpec &spec) noexcept {
    switch (stage_) {
    case Stage::Attack:
      if (spec.mode == CvEnvelopeMode::Ad)
        begin(Stage::Decay, 0.f);
      else if (spec.mode == CvEnvelopeMode::Ar)
        begin(Stage::Hold, peak_);
      else
        begin(Stage::Decay, peak_ * spec.sustain);
      break;
    case Stage::Decay:
      if (spec.mode == CvEnvelopeMode::Ad) {
        stage_ = Stage::Idle;
        value_ = 0.f;
      } else {
        stage_ = Stage::Sustain;
        value_ = peak_ * spec.sustain;
      }
      break;
    case Stage::Release:
      stage_ = Stage::Idle;
      value_ = 0.f;
      break;
    default:
      break;
    }
  }

  Stage stage_ = Stage::Idle;
  CvEnvelopeMode mode_ = CvEnvelopeMode::Ad;
  float value_ = 0.f;
  float start_ = 0.f;
  float target_ = 0.f;
  float peak_ = 0.f;
  double progress_ = 0.0;
  bool lastGate_ = false;
  bool modeKnown_ = false;
};

} // namespace tfseq
