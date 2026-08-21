#pragma once

#include "early_reflections_worker.hpp"
#include "late_reverb.hpp"
#include "reverb_defaults.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

namespace tfdsp {

struct RoomReverbControls {
  float space{reverb_defaults::Space};
  float aspect{reverb_defaults::Aspect};
  // Horizontal centring avoids a default ear bias; depth asymmetry avoids the
  // largest set of coincident image-source paths at the exact room centre.
  std::array<float, 3> listener{reverb_defaults::Listener};
  float preDelay{reverb_defaults::PreDelay};
  float decay{reverb_defaults::Decay};
  float damping{reverb_defaults::Damping};
  float diffusion{reverb_defaults::Diffusion};
  float modulation{reverb_defaults::Modulation};
  float shimmer{reverb_defaults::Shimmer};
  float width{reverb_defaults::Width};
  float earlyLevelDb{reverb_defaults::EarlyLevelDb};
  float tailLevelDb{reverb_defaults::TailLevelDb};
  float lowCut{reverb_defaults::LowCut};
  float highCut{reverb_defaults::HighCut};
};

class RoomReverb {
public:
  static constexpr std::size_t MaximumSources = EarlyReflectionMaximumSources;
  static constexpr std::size_t PartitionSize = 128;
  static constexpr double MaximumPreDelaySeconds = 0.250;
  using InputFrame = std::array<float, MaximumSources>;
  using SourcePosition = std::array<float, 3>;
  using SourcePositions = std::array<SourcePosition, MaximumSources>;
  using StereoFrame = std::array<float, 2>;
  using Convolver =
      EarlyReflectionConvolver<float, PartitionSize, MaximumSources>;

private:
  static constexpr float EarlyCalibrationGain = 1.f;
  static constexpr float PositionBalanceLimitDb = 9.f;
  static constexpr float PreDelayTransitionSeconds = 0.050f;
  static constexpr std::size_t ControlUpdateInterval = 64;
  class FractionalDelay {
  public:
    void Prepare(const std::size_t capacity) {
      if (capacity < 8)
        throw std::invalid_argument("room pre-delay capacity is too small");
      buffer_.assign(capacity, 0.f);
      writeIndex_ = 0;
    }

    void Reset() noexcept {
      std::fill(buffer_.begin(), buffer_.end(), 0.f);
      writeIndex_ = 0;
    }

    float Read(float delaySamples) const noexcept {
      if (buffer_.empty())
        return 0.f;
      delaySamples =
          std::clamp(delaySamples, 2.f, static_cast<float>(buffer_.size() - 3));
      const auto integer = static_cast<std::size_t>(std::floor(delaySamples));
      const float mu = delaySamples - static_cast<float>(integer);
      const std::array<float, 4> coefficients{
          -mu * (mu - 1.f) * (mu - 2.f) / 6.f,
          (mu + 1.f) * (mu - 1.f) * (mu - 2.f) / 2.f,
          -(mu + 1.f) * mu * (mu - 2.f) / 2.f,
          (mu + 1.f) * mu * (mu - 1.f) / 6.f,
      };
      constexpr std::array<int, 4> offsets{-1, 0, 1, 2};
      float value = 0.f;
      for (std::size_t tap = 0; tap < coefficients.size(); ++tap) {
        const auto distance =
            static_cast<std::size_t>(static_cast<int>(integer) + offsets[tap]);
        const auto index =
            (writeIndex_ + buffer_.size() - distance) % buffer_.size();
        value += coefficients[tap] * buffer_[index];
      }
      return value;
    }

    void Push(const float value) noexcept {
      buffer_[writeIndex_] = std::isfinite(value) ? value : 0.f;
      if (++writeIndex_ == buffer_.size())
        writeIndex_ = 0;
    }

  private:
    std::vector<float> buffer_{};
    std::size_t writeIndex_{};
  };

  struct SubmittedScene {
    float space{};
    float aspect{};
    std::array<float, 3> listener{};
    float damping{};
    SourcePositions positions{};
    std::size_t sourceCount{};
  };

  struct AutomaticBalance {
    InputFrame tailSend{};
    std::size_t sourceCount{};
  };

  // Fixed-capacity single-producer/single-consumer mailbox. The worker owns a
  // slot only while it is Writing and the audio thread owns it only while it
  // is Reading, so the POD payload is never accessed concurrently. This keeps
  // allocation, reference-count destruction, and library locks off the audio
  // thread while still allowing the newest completed scene to supersede older
  // ones.
  class AutomaticBalanceMailbox {
  public:
    bool Publish(const AutomaticBalance &value) noexcept {
      for (auto &slot : slots_) {
        auto expected = SlotState::Free;
        if (!slot.state.compare_exchange_strong(
                expected, SlotState::Writing, std::memory_order_acquire,
                std::memory_order_relaxed))
          continue;
        slot.value = value;
        slot.sequence = nextSequence_.fetch_add(1, std::memory_order_relaxed);
        slot.state.store(SlotState::Ready, std::memory_order_release);
        return true;
      }
      return false;
    }

    bool ConsumeLatest(AutomaticBalance &value,
                       std::uint64_t &sequence) noexcept {
      bool found = false;
      for (auto &slot : slots_) {
        auto expected = SlotState::Ready;
        if (!slot.state.compare_exchange_strong(
                expected, SlotState::Reading, std::memory_order_acquire,
                std::memory_order_relaxed))
          continue;
        if (!found || slot.sequence > sequence) {
          value = slot.value;
          sequence = slot.sequence;
          found = true;
        }
        slot.state.store(SlotState::Free, std::memory_order_release);
      }
      return found;
    }

    // The worker must be stopped before this is called.
    void Reset() noexcept {
      for (auto &slot : slots_)
        slot.state.store(SlotState::Free, std::memory_order_relaxed);
      nextSequence_.store(1, std::memory_order_relaxed);
    }

  private:
    enum class SlotState : std::uint8_t { Free, Writing, Ready, Reading };
    struct Slot {
      std::atomic<SlotState> state{SlotState::Free};
      AutomaticBalance value{};
      std::uint64_t sequence{};
    };
    std::array<Slot, 4> slots_{};
    std::atomic<std::uint64_t> nextSequence_{1};
  };

  double sampleRate_{48'000.0};
  EarlyReflectionConfig earlyConfig_{};
  Convolver convolver_{};
  LateReverb late_{};
  std::array<FractionalDelay, MaximumSources> preDelays_{};
  std::array<float, 2> highCutState_{};
  std::array<float, 2> lowCutState_{};
  std::unique_ptr<EarlyReflectionWorker> worker_{};
  SubmittedScene submitted_{};
  bool hasSubmittedScene_{};
  std::size_t requestCountdown_{};
  std::size_t requestIntervalSamples_{2'400};
  AutomaticBalanceMailbox balanceMailbox_{};
  std::uint64_t appliedBalanceSequence_{};
  InputFrame currentTailSend_{};
  InputFrame fromTailSend_{};
  InputFrame targetTailSend_{};
  float balanceTransitionPhase_{1.f};
  float balanceTransitionIncrement_{1.f};
  std::size_t balancePollCountdown_{};
  float wetSizeGain_{1.f};
  float targetWetSizeGain_{1.f};
  bool wetSizeGainInitialized_{};
  InputFrame currentEarlyPositionGain_{};
  InputFrame currentTailPositionGain_{};
  InputFrame targetEarlyPositionGain_{};
  InputFrame targetTailPositionGain_{};
  bool positionGainInitialized_{};
  float positionGainSmoothingCoefficient_{1.f};
  EarlyReflectionRoom controlRoom_{};
  float earlyGain_{1.f};
  float tailGain_{1.f};
  float highCutAlpha_{1.f};
  float lowCutAlpha_{1.f};
  std::size_t controlCountdown_{};
  float currentPreDelaySamples_{};
  float fromPreDelaySamples_{};
  float toPreDelaySamples_{};
  float pendingPreDelaySamples_{};
  float preDelayTransitionPhase_{1.f};
  float preDelayTransitionIncrement_{1.f};
  bool preDelayInitialized_{};

  static float ClampControl(const float value) noexcept {
    return std::clamp(std::isfinite(value) ? value : 0.f, 0.f, 1.f);
  }

  static float Smoothstep(const float value) noexcept {
    const float limited = ClampControl(value);
    return limited * limited * (3.f - 2.f * limited);
  }

  static float DecibelsToGain(const float decibels) noexcept {
    if (!std::isfinite(decibels) || decibels <= -60.f)
      return 0.f;
    return std::pow(10.f, std::clamp(decibels, -60.f, 6.f) / 20.f);
  }

  static float UnclampedDecibelsToGain(const float decibels) noexcept {
    return std::pow(10.f, decibels / 20.f);
  }

  static double ReferenceSourceListenerDistanceMetres() noexcept {
    double squared = 0.0;
    for (std::size_t axis = 0; axis < 3; ++axis) {
      const double difference =
          static_cast<double>(reverb_defaults::Source[axis] -
                              reverb_defaults::Listener[axis]) *
          static_cast<double>(reverb_defaults::RoomDimensionsMetres[axis]);
      squared += difference * difference;
    }
    return std::sqrt(squared);
  }

  static std::array<float, 2>
  PositionBalanceGains(const EarlyReflectionRoom &room,
                       const SourcePosition &normalizedSource) noexcept {
    double squared = 0.0;
    for (std::size_t axis = 0; axis < 3; ++axis) {
      const double source =
          static_cast<double>(ClampControl(normalizedSource[axis])) *
          room.dimensionsMetres[axis];
      const double difference = source - room.listenerPositionMetres[axis];
      squared += difference * difference;
    }
    // In a diffuse-field model late energy varies much less with receiver
    // distance than the direct and discrete-reflection fields. Since this
    // module deliberately leaves the dry path untouched, encode that distance
    // cue as a bounded, equal-and-opposite ER/tail balance around the factory
    // placement. One metre of ratio change follows the usual 20 log10 law;
    // the clamp avoids singular behaviour when the two markers coincide.
    const double distance = std::max(0.25, std::sqrt(squared));
    const double reference =
        std::max(0.25, ReferenceSourceListenerDistanceMetres());
    const float balanceDb = static_cast<float>(std::clamp(
        20.0 * std::log10(distance / reference),
        -static_cast<double>(PositionBalanceLimitDb),
        static_cast<double>(PositionBalanceLimitDb)));
    return {UnclampedDecibelsToGain(-0.5f * balanceDb),
            UnclampedDecibelsToGain(0.5f * balanceDb)};
  }

  void UpdatePreDelayTarget(const float targetSamples) noexcept {
    if (!preDelayInitialized_) {
      currentPreDelaySamples_ = fromPreDelaySamples_ = toPreDelaySamples_ =
          pendingPreDelaySamples_ = targetSamples;
      preDelayTransitionPhase_ = 1.f;
      preDelayInitialized_ = true;
      return;
    }
    pendingPreDelaySamples_ = targetSamples;
    if (preDelayTransitionPhase_ < 1.f ||
        std::abs(pendingPreDelaySamples_ - currentPreDelaySamples_) < 0.25f)
      return;
    fromPreDelaySamples_ = currentPreDelaySamples_;
    toPreDelaySamples_ = pendingPreDelaySamples_;
    preDelayTransitionPhase_ = 0.f;
  }

  float ReadPreDelay(const std::size_t source, const float liveSample,
                     const float delaySamples) const noexcept {
    return delaySamples < 0.5f ? liveSample
                              : preDelays_[source].Read(delaySamples);
  }

  void AdvancePreDelayTransition() noexcept {
    if (preDelayTransitionPhase_ >= 1.f)
      return;
    preDelayTransitionPhase_ =
        std::min(1.f, preDelayTransitionPhase_ + preDelayTransitionIncrement_);
    if (preDelayTransitionPhase_ < 1.f)
      return;
    currentPreDelaySamples_ = toPreDelaySamples_;
    if (std::abs(pendingPreDelaySamples_ - currentPreDelaySamples_) >= 0.25f) {
      fromPreDelaySamples_ = currentPreDelaySamples_;
      toPreDelaySamples_ = pendingPreDelaySamples_;
      preDelayTransitionPhase_ = 0.f;
    }
  }

  static float WetSizeCalibration(const EarlyReflectionRoom &room) noexcept {
    const double volume = room.dimensionsMetres[0] * room.dimensionsMetres[1] *
                          room.dimensionsMetres[2];
    constexpr double ReferenceVolume =
        static_cast<double>(reverb_defaults::RoomDimensionsMetres[0]) *
        static_cast<double>(reverb_defaults::RoomDimensionsMetres[1]) *
        static_cast<double>(reverb_defaults::RoomDimensionsMetres[2]);
    // Image-source pressure falls approximately with inverse path length. A
    // cube-root volume ratio is the characteristic linear scale, and keeps a
    // Size sweep from also becoming an unintended wet/dry-distance sweep.
    return static_cast<float>(std::clamp(
        std::cbrt(std::max(volume, 1.e-9) / ReferenceVolume), 0.30, 3.0));
  }

  static AutomaticBalance
  MakeAutomaticBalance(const EarlyReflectionImpulseResponse &response) {
    AutomaticBalance balance;
    balance.tailSend.fill(1.f);
    balance.sourceCount = response.sourceCount;
    // The absolute reference belongs to the exported late topology. Geometry
    // changes only the pre-loop source send; it never changes feedback poles.
    constexpr double epsilon = 1.e-12;
    for (std::size_t source = 0; source < response.sourceCount; ++source) {
      const auto &handoff = response.sourceHandoffs[source];
      double weightedPower = 0.0;
      constexpr std::array<double, EarlyReflectionBandCount> weights{
          0.16, 0.29, 0.34, 0.21};
      for (std::size_t band = 0; band < EarlyReflectionBandCount; ++band)
        weightedPower += weights[band] * handoff.bandHandoffPower[band];
      if (!(weightedPower > epsilon))
        weightedPower = handoff.broadbandHandoffPower;
      const double matched =
          std::sqrt(std::max(weightedPower, epsilon) /
                    std::max<double>(late_reverb_coefficients::UnitHandoffPower,
                                     epsilon));
      balance.tailSend[source] =
          static_cast<float>(std::clamp(matched, 0.35, 3.0));
    }
    return balance;
  }

  void UpdateAutomaticBalance() noexcept {
    if (balancePollCountdown_ == 0) {
      AutomaticBalance latest;
      auto sequence = appliedBalanceSequence_;
      if (balanceMailbox_.ConsumeLatest(latest, sequence) &&
          sequence != appliedBalanceSequence_) {
        appliedBalanceSequence_ = sequence;
        fromTailSend_ = currentTailSend_;
        targetTailSend_ = latest.tailSend;
        balanceTransitionPhase_ = 0.f;
      }
      balancePollCountdown_ = 64;
    }
    --balancePollCountdown_;
    if (balanceTransitionPhase_ >= 1.f)
      return;
    const float fade = balanceTransitionPhase_ * balanceTransitionPhase_ *
                       (3.f - 2.f * balanceTransitionPhase_);
    for (std::size_t source = 0; source < MaximumSources; ++source)
      currentTailSend_[source] =
          fromTailSend_[source] +
          fade * (targetTailSend_[source] - fromTailSend_[source]);
    balanceTransitionPhase_ =
        std::min(1.f, balanceTransitionPhase_ + balanceTransitionIncrement_);
    if (balanceTransitionPhase_ >= 1.f)
      currentTailSend_ = targetTailSend_;
  }

  static StereoFrame ApplyWidth(const StereoFrame &stereo,
                                const float control) noexcept {
    const float width = 1.5f * Smoothstep(control);
    constexpr float inverseRootTwo = 0.7071067811865475f;
    const float mid = (stereo[0] + stereo[1]) * inverseRootTwo;
    const float side = width * (stereo[0] - stereo[1]) * inverseRootTwo;
    return {(mid + side) * inverseRootTwo, (mid - side) * inverseRootTwo};
  }

  static bool Different(const float left, const float right,
                        const float threshold = 0.005f) noexcept {
    return std::abs(left - right) >= threshold;
  }

  bool SceneChanged(const RoomReverbControls &controls,
                    const SourcePositions &positions,
                    const std::size_t sourceCount) const noexcept {
    if (!hasSubmittedScene_ || sourceCount != submitted_.sourceCount ||
        Different(controls.space, submitted_.space) ||
        Different(controls.aspect, submitted_.aspect) ||
        Different(controls.damping, submitted_.damping))
      return true;
    for (std::size_t axis = 0; axis < 3; ++axis)
      if (Different(controls.listener[axis], submitted_.listener[axis]))
        return true;
    for (std::size_t source = 0; source < sourceCount; ++source)
      for (std::size_t axis = 0; axis < 3; ++axis)
        if (Different(positions[source][axis],
                      submitted_.positions[source][axis]))
          return true;
    return false;
  }

  void UpdateControlTargets(const RoomReverbControls &controls,
                            const SourcePositions &positions,
                            const std::size_t sourceCount) noexcept {
    controlRoom_ = MakeRoom(controls);
    for (std::size_t source = 0; source < MaximumSources; ++source) {
      const auto target = source < sourceCount
                              ? PositionBalanceGains(controlRoom_,
                                                     positions[source])
                              : std::array<float, 2>{1.f, 1.f};
      targetEarlyPositionGain_[source] = target[0];
      targetTailPositionGain_[source] = target[1];
      if (!positionGainInitialized_) {
        currentEarlyPositionGain_[source] = target[0];
        currentTailPositionGain_[source] = target[1];
      }
    }
    positionGainInitialized_ = true;

    targetWetSizeGain_ = WetSizeCalibration(controlRoom_);
    if (!wetSizeGainInitialized_) {
      wetSizeGain_ = targetWetSizeGain_;
      wetSizeGainInitialized_ = true;
    }
    earlyGain_ =
        EarlyCalibrationGain * DecibelsToGain(controls.earlyLevelDb);
    tailGain_ = DecibelsToGain(controls.tailLevelDb);

    const float highControl = ClampControl(controls.highCut);
    const float lowControl = ClampControl(controls.lowCut);
    const float highCutHz =
        std::min(0.45f * static_cast<float>(sampleRate_),
                 1'000.f * std::exp(highControl * std::log(20.f)));
    const float lowCutHz = 20.f * std::exp(lowControl * std::log(50.f));
    highCutAlpha_ =
        1.f - std::exp(-2.f * static_cast<float>(EarlyReflectionPi) *
                       highCutHz / static_cast<float>(sampleRate_));
    lowCutAlpha_ =
        1.f - std::exp(-2.f * static_cast<float>(EarlyReflectionPi) *
                       lowCutHz / static_cast<float>(sampleRate_));
  }

  void SubmitScene(const RoomReverbControls &controls,
                   const SourcePositions &positions,
                   const std::size_t sourceCount) noexcept {
    if (!worker_ || sourceCount == 0 ||
        !SceneChanged(controls, positions, sourceCount))
      return;
    EarlyReflectionBuildRequest request;
    request.config = earlyConfig_;
    request.room = MakeRoom(controls);
    request.sourceCount = sourceCount;
    for (std::size_t source = 0; source < sourceCount; ++source) {
      std::array<double, 3> normalized{};
      for (std::size_t axis = 0; axis < 3; ++axis)
        normalized[axis] = ClampControl(positions[source][axis]);
      request.sources[source] =
          MakeEarlyReflectionSource(request.room, normalized);
    }
    request.materials =
        MakeEarlyReflectionMaterials(ClampControl(controls.damping));
    request.convolutionLatencySamples = Convolver::LatencySamples;
    request.transitionSeconds = 0.100;
    if (worker_->Submit(request) == 0)
      return;

    submitted_.space = controls.space;
    submitted_.aspect = controls.aspect;
    submitted_.listener = controls.listener;
    submitted_.damping = controls.damping;
    submitted_.positions = positions;
    submitted_.sourceCount = sourceCount;
    hasSubmittedScene_ = true;
  }

  StereoFrame FilterWet(StereoFrame wet) noexcept {
    for (std::size_t channel = 0; channel < wet.size(); ++channel) {
      highCutState_[channel] +=
          highCutAlpha_ * (wet[channel] - highCutState_[channel]);
      wet[channel] = highCutState_[channel];
      lowCutState_[channel] +=
          lowCutAlpha_ * (wet[channel] - lowCutState_[channel]);
      wet[channel] -= lowCutState_[channel];
    }
    return wet;
  }

public:
  RoomReverb() { SetSampleRate(sampleRate_); }

  ~RoomReverb() {
    if (worker_)
      worker_->Stop();
  }

  RoomReverb(const RoomReverb &) = delete;
  RoomReverb &operator=(const RoomReverb &) = delete;

  static EarlyReflectionRoom MakeRoom(const RoomReverbControls &controls) {
    auto room = MakeEarlyReflectionRoom(ClampControl(controls.space));
    const double aspect =
        std::exp((2.0 * ClampControl(controls.aspect) - 1.0) * std::log(1.8));
    const double rootAspect = std::sqrt(aspect);
    room.dimensionsMetres[0] *= rootAspect;
    room.dimensionsMetres[1] /= rootAspect;
    for (std::size_t axis = 0; axis < 3; ++axis) {
      const double fraction =
          std::clamp(static_cast<double>(ClampControl(controls.listener[axis])),
                     0.02, 0.98);
      room.listenerPositionMetres[axis] =
          fraction * room.dimensionsMetres[axis];
    }
    room.Validate();
    return room;
  }

  void SetSampleRate(const double sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument("room reverb sample rate must be positive");
    if (worker_) {
      worker_->Stop();
      worker_.reset();
    }
    balanceMailbox_.Reset();
    sampleRate_ = sampleRate;
    earlyConfig_ = {};
    earlyConfig_.sampleRate = sampleRate_;
    late_.SetSampleRate(sampleRate_);

    EarlyReflectionRoom maximumRoom;
    maximumRoom.dimensionsMetres = {25.0, 34.0, 8.0};
    maximumRoom.listenerPositionMetres = {0.02, 0.02, 0.02};
    convolver_.Prepare(sampleRate_, MaximumEarlyReflectionImpulseSamples(
                                        earlyConfig_, maximumRoom));
    const auto preDelayCapacity = static_cast<std::size_t>(std::ceil(
                                      MaximumPreDelaySeconds * sampleRate_)) +
                                  8;
    for (auto &delay : preDelays_)
      delay.Prepare(preDelayCapacity);
    requestIntervalSamples_ = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::llround(sampleRate_ / 20.0)));
    balanceTransitionIncrement_ =
        1.f / std::max(1.f, 0.100f * static_cast<float>(sampleRate_));
    positionGainSmoothingCoefficient_ =
        1.f - std::exp(-1.f / (0.050f * static_cast<float>(sampleRate_)));
    preDelayTransitionIncrement_ =
        1.f / std::max(1.f, PreDelayTransitionSeconds *
                               static_cast<float>(sampleRate_));
    requestCountdown_ = 0;
    hasSubmittedScene_ = false;
    worker_ = std::make_unique<EarlyReflectionWorker>(
        20.0, [this](const EarlyReflectionImpulseResponse &response,
                     const double transitionSeconds) {
          if (!convolver_.PrepareAndQueueImpulseResponse(response,
                                                         transitionSeconds))
            return false;
          balanceMailbox_.Publish(MakeAutomaticBalance(response));
          return true;
        });
    Reset();
  }

  void Reset() noexcept {
    convolver_.Reset();
    late_.Reset();
    for (auto &delay : preDelays_)
      delay.Reset();
    highCutState_.fill(0.f);
    lowCutState_.fill(0.f);
    currentTailSend_.fill(1.f);
    fromTailSend_.fill(1.f);
    targetTailSend_.fill(1.f);
    balanceTransitionPhase_ = 1.f;
    balancePollCountdown_ = 0;
    appliedBalanceSequence_ = 0;
    wetSizeGain_ = 1.f;
    wetSizeGainInitialized_ = false;
    currentEarlyPositionGain_.fill(1.f);
    currentTailPositionGain_.fill(1.f);
    targetEarlyPositionGain_.fill(1.f);
    targetTailPositionGain_.fill(1.f);
    positionGainInitialized_ = false;
    targetWetSizeGain_ = 1.f;
    earlyGain_ = tailGain_ = 1.f;
    highCutAlpha_ = lowCutAlpha_ = 1.f;
    controlCountdown_ = 0;
    currentPreDelaySamples_ = fromPreDelaySamples_ = toPreDelaySamples_ =
        pendingPreDelaySamples_ = 0.f;
    preDelayTransitionPhase_ = 1.f;
    preDelayInitialized_ = false;
  }

  double SampleRate() const noexcept { return sampleRate_; }
  bool HasSubmittedScene() const noexcept { return hasSubmittedScene_; }

  bool WaitForLatestScene(const std::chrono::milliseconds timeout) {
    if (!worker_)
      return false;
    const auto result = worker_->WaitForLatestResult(timeout);
    return result && result->Succeeded() && result->publishedToConvolver;
  }

  StereoFrame Process(const InputFrame &input, const SourcePositions &positions,
                      const std::size_t sourceCount,
                      const RoomReverbControls &controls) noexcept {
    const std::size_t activeSources = std::min(sourceCount, MaximumSources);
    if (controlCountdown_ == 0) {
      UpdateControlTargets(controls, positions, activeSources);
      controlCountdown_ = ControlUpdateInterval;
    }
    --controlCountdown_;
    if (requestCountdown_ == 0) {
      SubmitScene(controls, positions, activeSources);
      requestCountdown_ = requestIntervalSamples_;
    }
    --requestCountdown_;

    const float preDelaySamples = MaximumPreDelaySeconds *
                                  ClampControl(controls.preDelay) *
                                  static_cast<float>(sampleRate_);
    UpdatePreDelayTarget(preDelaySamples);
    InputFrame delayed{};
    for (std::size_t source = 0; source < MaximumSources; ++source) {
      const float sample =
          source < activeSources && std::isfinite(input[source]) ? input[source]
                                                                 : 0.f;
      if (preDelayTransitionPhase_ >= 1.f) {
        delayed[source] =
            ReadPreDelay(source, sample, currentPreDelaySamples_);
      } else {
        const float fade = preDelayTransitionPhase_ * preDelayTransitionPhase_ *
                           (3.f - 2.f * preDelayTransitionPhase_);
        const float from =
            ReadPreDelay(source, sample, fromPreDelaySamples_);
        const float to = ReadPreDelay(source, sample, toPreDelaySamples_);
        delayed[source] = from + fade * (to - from);
      }
      preDelays_[source].Push(sample);
    }
    AdvancePreDelayTransition();

    InputFrame earlyInput{};
    InputFrame lateInputGain{};
    for (std::size_t source = 0; source < MaximumSources; ++source) {
      currentEarlyPositionGain_[source] +=
          positionGainSmoothingCoefficient_ *
          (targetEarlyPositionGain_[source] -
           currentEarlyPositionGain_[source]);
      currentTailPositionGain_[source] +=
          positionGainSmoothingCoefficient_ *
          (targetTailPositionGain_[source] -
           currentTailPositionGain_[source]);
      earlyInput[source] = delayed[source] * currentEarlyPositionGain_[source];
    }

    const auto early = convolver_.Process(earlyInput, activeSources);
    UpdateAutomaticBalance();
    for (std::size_t source = 0; source < MaximumSources; ++source)
      lateInputGain[source] =
          currentTailSend_[source] * currentTailPositionGain_[source];
    LateReverbControls lateControls;
    lateControls.decay = controls.decay;
    lateControls.damping = controls.damping;
    lateControls.diffusion = controls.diffusion;
    lateControls.modulation = controls.modulation;
    lateControls.shimmer = controls.shimmer;
    lateControls.listener = controls.listener;
    for (std::size_t axis = 0; axis < 3; ++axis)
      lateControls.roomDimensionsMetres[axis] =
          static_cast<float>(controlRoom_.dimensionsMetres[axis]);
    const auto tail = late_.Process(delayed, positions, activeSources,
                                    lateControls, lateInputGain);
    wetSizeGain_ += balanceTransitionIncrement_ *
                    (targetWetSizeGain_ - wetSizeGain_);
    const auto wet = ApplyWidth(
        {wetSizeGain_ * (earlyGain_ * early[0] + tailGain_ * tail[0]),
         wetSizeGain_ * (earlyGain_ * early[1] + tailGain_ * tail[1])},
        controls.width);
    return FilterWet(wet);
  }
};

} // namespace tfdsp
