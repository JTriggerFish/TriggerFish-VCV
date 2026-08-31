#pragma once

#include "deterministic_random.hpp"
#include "spectral_tilt_filter.hpp"
#include "stochastic_event_scheduler.hpp"
#include "tfdsp/finite_audio.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace tfdsp::percussion {

struct MicroContactProcessParameters {
  float densityHz{4000.f};
  float contactDurationSeconds{.0007f};
  float amplitude{1.f};
  float tiltDb{};
  float tiltPivotHz{4000.f};
  float attackSeconds{.001f};
  float releaseSeconds{.015f};
  float densityNormalization{1.f};
  std::uint32_t seed{1};
};

// A renewal-event stream whose events raise smooth overlapping contact
// windows over one correlated noise field. Low rates expose discrete contacts;
// dense rates fuse into a continuous scrape without hard grain boundaries.
class MicroContactProcess {
public:
  void Prepare(const float sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate < 1.f)
      throw std::invalid_argument("micro-contact process rate must be positive");
    sampleRate_ = sampleRate;
    scheduler_.Prepare(sampleRate);
    tilt_.Prepare(sampleRate);
    Reset();
  }

  void Reset() noexcept {
    scheduler_.Reset(1);
    random_.Seed(1);
    tilt_.Reset();
    for (auto &voice : voices_)
      voice = {};
    gate_ = targetGate_ = 0.f;
    clusterSamplesRemaining_ = 0;
    streaming_ = false;
  }

  void StartStream(const MicroContactProcessParameters &parameters) noexcept {
    Configure(parameters);
    streaming_ = true;
    targetGate_ = 1.f;
  }

  void StopStream() noexcept {
    streaming_ = false;
    targetGate_ = 0.f;
  }

  void SetDensityHz(const float densityHz) noexcept {
    densityHz_ = std::clamp(
        std::isfinite(densityHz) ? densityHz : 0.f, 0.f,
        .45f * sampleRate_);
    scheduler_.SetRateHz(densityHz_);
    UpdateDensityScale();
  }

  void SetAmplitude(const float amplitude) noexcept {
    amplitude_ = std::clamp(tfdsp::FiniteNormalOrZero(amplitude), 0.f, 16.f);
  }

  void SetTilt(const float tiltDb, const float pivotHz) noexcept {
    tilt_.SetTilt(tiltDb, pivotHz);
  }

  void TriggerCluster(const MicroContactProcessParameters &parameters,
                      float durationSeconds) noexcept {
    Configure(parameters);
    durationSeconds = std::clamp(
        std::isfinite(durationSeconds) ? durationSeconds : 0.f, 0.f, 10.f);
    clusterSamplesRemaining_ = static_cast<std::size_t>(
        std::lround(durationSeconds * sampleRate_));
    targetGate_ = clusterSamplesRemaining_ > 0 ? 1.f : 0.f;
  }

  float Process() noexcept {
    if (!Active())
      return 0.f;
    const bool scheduling = streaming_ || clusterSamplesRemaining_ > 0;
    if (scheduling && scheduler_.Process())
      TriggerContact();
    if (clusterSamplesRemaining_ > 0 && --clusterSamplesRemaining_ == 0)
      targetGate_ = streaming_ ? 1.f : 0.f;
    const float coefficient = targetGate_ > gate_ ? attackCoefficient_
                                                   : releaseCoefficient_;
    gate_ += coefficient * (targetGate_ - gate_);
    gate_ = tfdsp::FiniteNormalOrZero(gate_);
    const float envelope = ContactEnvelope();
    const float noise = tilt_.Process(random_.Bipolar());
    return tfdsp::FiniteNormalOrZero(
        amplitude_ * densityScale_ * gate_ * envelope * noise);
  }

  bool Active() const noexcept {
    if (streaming_ || clusterSamplesRemaining_ > 0 || gate_ > 1.e-5f)
      return true;
    for (const auto &voice : voices_)
      if (voice.remaining > 0)
        return true;
    return false;
  }

private:
  static constexpr std::size_t VoiceCount = 12;

  struct Voice {
    float sine{};
    float cosine{1.f};
    float rotationSine{};
    float rotationCosine{1.f};
    float amplitude{};
    std::size_t remaining{};
  };

  static float Coefficient(const float seconds, const float sampleRate) noexcept {
    const float bounded = std::clamp(
        std::isfinite(seconds) ? seconds : 0.f, 0.f, 10.f);
    return bounded == 0.f ? 1.f : -std::expm1(-1.f / (bounded * sampleRate));
  }

  void Configure(const MicroContactProcessParameters &parameters) noexcept {
    Reset();
    scheduler_.Reset(parameters.seed);
    random_.Seed(parameters.seed ^ 0x9e3779b9u);
    densityNormalization_ = std::clamp(
        std::isfinite(parameters.densityNormalization)
            ? parameters.densityNormalization : 1.f,
        0.f, 1.f);
    SetDensityHz(parameters.densityHz);
    contactSamples_ = std::max<std::size_t>(3,
        static_cast<std::size_t>(std::lround(std::clamp(
            std::isfinite(parameters.contactDurationSeconds)
                ? parameters.contactDurationSeconds : 0.f,
            1.f / sampleRate_, .1f) * sampleRate_)));
    SetAmplitude(parameters.amplitude);
    SetTilt(parameters.tiltDb, parameters.tiltPivotHz);
    attackCoefficient_ = Coefficient(parameters.attackSeconds, sampleRate_);
    releaseCoefficient_ = Coefficient(parameters.releaseSeconds, sampleRate_);
    UpdateDensityScale();
  }

  void UpdateDensityScale() noexcept {
    const float occupancy = std::max(
        1.f, densityHz_ * static_cast<float>(contactSamples_) / sampleRate_);
    densityScale_ = std::pow(occupancy, -.5f * densityNormalization_);
  }

  void TriggerContact() noexcept {
    Voice *selected = &voices_[0];
    for (auto &voice : voices_) {
      if (voice.remaining == 0) {
        selected = &voice;
        break;
      }
      if (voice.remaining < selected->remaining)
        selected = &voice;
    }
    constexpr double Pi = 3.1415926535897932384626433832795;
    const double step = Pi / static_cast<double>(contactSamples_ - 1);
    selected->rotationSine = static_cast<float>(std::sin(step));
    selected->rotationCosine = static_cast<float>(std::cos(step));
    selected->sine = 0.f;
    selected->cosine = 1.f;
    selected->amplitude = .5f + .5f * random_.Uniform();
    selected->remaining = contactSamples_;
  }

  float ContactEnvelope() noexcept {
    float envelope = 0.f;
    for (auto &voice : voices_) {
      if (voice.remaining == 0)
        continue;
      envelope += voice.amplitude * voice.sine;
      const float nextSine = voice.sine * voice.rotationCosine +
                             voice.cosine * voice.rotationSine;
      voice.cosine = voice.cosine * voice.rotationCosine -
                     voice.sine * voice.rotationSine;
      voice.sine = nextSine;
      --voice.remaining;
    }
    return tfdsp::FiniteNormalOrZero(envelope);
  }

  std::array<Voice, VoiceCount> voices_{};
  StochasticEventScheduler scheduler_{};
  DeterministicRandom random_{};
  SpectralTiltFilter tilt_{};
  float sampleRate_{48000.f};
  float amplitude_{1.f};
  float densityHz_{};
  float densityScale_{1.f};
  float densityNormalization_{1.f};
  float gate_{};
  float targetGate_{};
  float attackCoefficient_{1.f};
  float releaseCoefficient_{1.f};
  std::size_t contactSamples_{1};
  std::size_t clusterSamplesRemaining_{};
  bool streaming_{};
};

} // namespace tfdsp::percussion
