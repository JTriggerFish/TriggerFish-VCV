#pragma once

#include "biquad.hpp"

#include <algorithm>
#include <cmath>

namespace tfdsp::percussion::biquad_design {

inline float BoundFrequency(float frequencyHz, const float sampleRate) noexcept {
  if (!std::isfinite(frequencyHz))
    frequencyHz = 1000.f;
  return std::clamp(frequencyHz, 1.f, .49f * sampleRate);
}

inline BiquadCoefficients Normalize(const float b0, const float b1,
                                    const float b2, const float a0,
                                    const float a1, const float a2) noexcept {
  const float inverseA0 = 1.f / a0;
  return {b0 * inverseA0, b1 * inverseA0, b2 * inverseA0,
          a1 * inverseA0, a2 * inverseA0};
}

inline BiquadCoefficients Lowpass(float frequencyHz, float q,
                                  const float sampleRate) noexcept {
  frequencyHz = BoundFrequency(frequencyHz, sampleRate);
  q = std::clamp(std::isfinite(q) ? q : .707f, .1f, 20.f);
  const float omega = 6.283185307179586f * frequencyHz / sampleRate;
  const float cosine = std::cos(omega);
  const float alpha = std::sin(omega) / (2.f * q);
  return Normalize(.5f * (1.f - cosine), 1.f - cosine,
                   .5f * (1.f - cosine), 1.f + alpha,
                   -2.f * cosine, 1.f - alpha);
}

inline BiquadCoefficients Highpass(float frequencyHz, float q,
                                   const float sampleRate) noexcept {
  frequencyHz = BoundFrequency(frequencyHz, sampleRate);
  q = std::clamp(std::isfinite(q) ? q : .707f, .1f, 20.f);
  const float omega = 6.283185307179586f * frequencyHz / sampleRate;
  const float cosine = std::cos(omega);
  const float alpha = std::sin(omega) / (2.f * q);
  return Normalize(.5f * (1.f + cosine), -(1.f + cosine),
                   .5f * (1.f + cosine), 1.f + alpha,
                   -2.f * cosine, 1.f - alpha);
}

inline BiquadCoefficients Peaking(float frequencyHz, float q, float gainDb,
                                  const float sampleRate) noexcept {
  frequencyHz = BoundFrequency(frequencyHz, sampleRate);
  q = std::clamp(std::isfinite(q) ? q : 1.f, .1f, 20.f);
  gainDb = std::clamp(std::isfinite(gainDb) ? gainDb : 0.f, -36.f, 36.f);
  const float amplitude = std::pow(10.f, gainDb / 40.f);
  const float omega = 6.283185307179586f * frequencyHz / sampleRate;
  const float cosine = std::cos(omega);
  const float alpha = std::sin(omega) / (2.f * q);
  return Normalize(1.f + alpha * amplitude, -2.f * cosine,
                   1.f - alpha * amplitude, 1.f + alpha / amplitude,
                   -2.f * cosine, 1.f - alpha / amplitude);
}

} // namespace tfdsp::percussion::biquad_design
