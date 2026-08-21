#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace tfdsp {

// A deliberately diffuse octave shifter for use inside a reverb feedback
// path. Eight uniformly staggered grains avoid the exposed two-head crossfade
// of a conventional delay-line shifter. Their equal durations preserve that
// spacing (and therefore a constant overlap sum); each restart gets only a
// small deterministic random look-back and pitch offset, turning a fixed comb
// of modulation sidebands into a quiet, decorrelated skirt without allowing
// grain phases to bunch together and expose clicks.
//
// An eighth-order Butterworth low-pass precedes the grain buffer. Its cutoff is
// derived from the +12-semitone ratio so out-of-band input cannot repeatedly
// alias back into the audible range as it circulates through the reverb.
class WindowedPitchShifter {
public:
  static constexpr std::size_t GrainCount = 8;

  void Prepare(const double sampleRate, const float windowSeconds = 0.120f,
               const float initialPhase = 0.f,
               const std::uint32_t randomSeed = 0x6d2b79f5u) {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument("pitch shifter sample rate must be positive");
    if (!std::isfinite(windowSeconds) || windowSeconds <= 0.f)
      throw std::invalid_argument("pitch shifter window must be positive");

    sampleRate_ = sampleRate;
    windowSamples_ =
        std::max(64.f, windowSeconds * static_cast<float>(sampleRate));
    initialPhase_ = WrapPhase(initialPhase);
    initialSeed_ = randomSeed != 0 ? randomSeed : 0x6d2b79f5u;
    lookbackJitterSamples_ = 0.005f * static_cast<float>(sampleRate_);
    safetySamples_ = 0.010f * static_cast<float>(sampleRate_);
    const float maximumLookback =
        1.10f * windowSamples_ + safetySamples_ + lookbackJitterSamples_;
    buffer_.assign(static_cast<std::size_t>(std::ceil(maximumLookback)) + 16,
                   0.f);

    // For an octave-up ratio, input above Nyquist / 2 cannot have a legal
    // output. A conservative transition band gives the eighth-order filter
    // room to reject it instead of allowing metallic aliases into the loop.
    constexpr float PitchRatio = 2.f;
    const float cutoff =
        std::min(0.45f * static_cast<float>(sampleRate_),
                 0.72f * 0.5f * static_cast<float>(sampleRate_) / PitchRatio);
    constexpr std::array<float, 4> ButterworthQ{0.5097955791f, 0.6013448869f,
                                                0.8999762231f, 2.5629154477f};
    for (std::size_t stage = 0; stage < antiAlias_.size(); ++stage)
      antiAlias_[stage].SetLowpass(sampleRate_, cutoff, ButterworthQ[stage]);
    Reset();
  }

  void Reset() noexcept {
    std::fill(buffer_.begin(), buffer_.end(), 0.f);
    for (auto &filter : antiAlias_)
      filter.Reset();
    writeIndex_ = 0;
    randomState_ = initialSeed_;
    for (std::size_t index = 0; index < grains_.size(); ++index) {
      auto &grain = grains_[index];
      grain.phase =
          WrapPhase(initialPhase_ + static_cast<float>(index) / GrainCount);
      StartGrain(grain, grain.phase);
    }
  }

  float Process(const float input) noexcept {
    return ProcessInternal(input, true);
  }

  // Keep the input history, anti-alias filters, grain phases and deterministic
  // random sequence current without paying for interpolated grain reads and
  // window evaluation. This makes a zero-gain shimmer branch cheap while
  // allowing it to become audible immediately when its gain is raised.
  void Advance(const float input) noexcept { ProcessInternal(input, false); }

  double SampleRate() const noexcept { return sampleRate_; }
  float WindowSamples() const noexcept { return windowSamples_; }

private:
  float ProcessInternal(const float input, const bool render) noexcept {
    if (buffer_.empty())
      return 0.f;

    float filtered = std::isfinite(input) ? input : 0.f;
    for (auto &filter : antiAlias_)
      filtered = filter.Process(filtered);
    buffer_[writeIndex_] = filtered;

    float output = 0.f;
    float windowSum = 0.f;
    for (auto &grain : grains_) {
      if (render) {
        const float window = Hann(grain.phase);
        output += window * Read(grain.delaySamples);
        windowSum += window;
      }

      grain.delaySamples -= grain.pitchRatio - 1.f;
      const float nextPhase = grain.phase + 1.f / grain.durationSamples;
      if (nextPhase >= 1.f)
        StartGrain(grain, 0.f);
      else
        grain.phase = nextPhase;
    }

    if (++writeIndex_ == buffer_.size())
      writeIndex_ = 0;
    if (!render)
      return 0.f;
    output /= std::max(windowSum, 1.e-6f);
    return std::isfinite(output) ? output : 0.f;
  }
  static constexpr float Pi = 3.14159265358979323846f;
  static constexpr float MinimumDelaySamples = 4.f;
  static constexpr float PitchWanderCents = 3.f;

  struct Biquad {
    float b0{1.f};
    float b1{};
    float b2{};
    float a1{};
    float a2{};
    float z1{};
    float z2{};

    void SetLowpass(const double sampleRate, const float cutoff,
                    const float q) noexcept {
      const float omega = 2.f * Pi * cutoff / static_cast<float>(sampleRate);
      const float cosine = std::cos(omega);
      const float sine = std::sin(omega);
      const float alpha = sine / (2.f * q);
      const float inverseA0 = 1.f / (1.f + alpha);
      b0 = 0.5f * (1.f - cosine) * inverseA0;
      b1 = (1.f - cosine) * inverseA0;
      b2 = b0;
      a1 = -2.f * cosine * inverseA0;
      a2 = (1.f - alpha) * inverseA0;
    }

    float Process(const float input) noexcept {
      const float output = b0 * input + z1;
      z1 = b1 * input - a1 * output + z2;
      z2 = b2 * input - a2 * output;
      return output;
    }

    void Reset() noexcept { z1 = z2 = 0.f; }
  };

  struct Grain {
    float phase{};
    float delaySamples{};
    float pitchRatio{2.f};
    float durationSamples{5'760.f};
  };

  static float WrapPhase(float phase) noexcept {
    phase -= std::floor(phase);
    return phase < 0.f ? phase + 1.f : phase;
  }

  static float Hann(const float phase) noexcept {
    const float sine = std::sin(Pi * phase);
    return sine * sine;
  }

  std::uint32_t RandomBits() noexcept {
    randomState_ ^= randomState_ << 13;
    randomState_ ^= randomState_ >> 17;
    randomState_ ^= randomState_ << 5;
    return randomState_;
  }

  float BipolarRandom() noexcept {
    constexpr float Scale = 1.f / static_cast<float>(0x00ffffffu);
    return 2.f * static_cast<float>(RandomBits() & 0x00ffffffu) * Scale - 1.f;
  }

  void StartGrain(Grain &grain, const float phase) noexcept {
    const float cents = PitchWanderCents * BipolarRandom();
    grain.pitchRatio = 2.f * std::exp2(cents / 1200.f);
    grain.durationSamples = windowSamples_;
    grain.phase = phase;
    const float jitter = lookbackJitterSamples_ * BipolarRandom();
    const float elapsed = phase * grain.durationSamples;
    grain.delaySamples = grain.durationSamples * (grain.pitchRatio - 1.f) +
                         safetySamples_ + jitter -
                         elapsed * (grain.pitchRatio - 1.f);
  }

  float Read(float delaySamples) const noexcept {
    delaySamples = std::clamp(delaySamples, MinimumDelaySamples,
                              static_cast<float>(buffer_.size() - 3));
    const auto integer = static_cast<std::size_t>(std::floor(delaySamples));
    const float mu = delaySamples - static_cast<float>(integer);
    const float c0 = -mu * (mu - 1.f) * (mu - 2.f) / 6.f;
    const float c1 = (mu + 1.f) * (mu - 1.f) * (mu - 2.f) / 2.f;
    const float c2 = -(mu + 1.f) * mu * (mu - 2.f) / 2.f;
    const float c3 = (mu + 1.f) * mu * (mu - 1.f) / 6.f;
    const auto sample = [this, integer](const int offset) {
      const auto distance =
          static_cast<std::size_t>(static_cast<int>(integer) + offset);
      return buffer_[(writeIndex_ + buffer_.size() - distance) %
                     buffer_.size()];
    };
    return c0 * sample(-1) + c1 * sample(0) + c2 * sample(1) + c3 * sample(2);
  }

  double sampleRate_{48'000.0};
  float windowSamples_{5'760.f};
  float initialPhase_{};
  float lookbackJitterSamples_{240.f};
  float safetySamples_{480.f};
  std::uint32_t initialSeed_{0x6d2b79f5u};
  std::uint32_t randomState_{initialSeed_};
  std::array<Biquad, 4> antiAlias_{};
  std::array<Grain, GrainCount> grains_{};
  std::vector<float> buffer_{};
  std::size_t writeIndex_{};
};

} // namespace tfdsp
