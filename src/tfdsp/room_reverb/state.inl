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
  static constexpr float PreDelayTransitionSeconds = 0.050f;
  static constexpr std::size_t ControlUpdateInterval = 64;
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
    InputFrame handoffStartSeconds{};
    std::size_t sourceCount{};
    std::size_t sceneSequence{};
  };

  // Fixed-capacity single-producer/single-consumer mailbox. The worker owns a
  // slot only while it is Writing and the audio thread owns it only while it
  // is Reading, so the POD payload is never accessed concurrently. This keeps
  // allocation, reference-count destruction, and library locks off the audio
  // thread. Future-scene metadata remains queued until the convolver adopts
  // the FIR bank carrying the same scene sequence.
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
        slot.state.store(SlotState::Ready, std::memory_order_release);
        return true;
      }
      return false;
    }

    bool ConsumeForScene(AutomaticBalance &value,
                         const std::size_t sceneSequence) noexcept {
      bool found = false;
      for (auto &slot : slots_) {
        auto expected = SlotState::Ready;
        if (!slot.state.compare_exchange_strong(
                expected, SlotState::Reading, std::memory_order_acquire,
                std::memory_order_relaxed))
          continue;
        if (slot.value.sceneSequence > sceneSequence) {
          slot.state.store(SlotState::Ready, std::memory_order_release);
          continue;
        }
        if (slot.value.sceneSequence == sceneSequence) {
          value = slot.value;
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
    }

  private:
    enum class SlotState : std::uint8_t { Free, Writing, Ready, Reading };
    struct Slot {
      std::atomic<SlotState> state{SlotState::Free};
      AutomaticBalance value{};
    };
    std::array<Slot, 8> slots_{};
  };

  double sampleRate_{48'000.0};
  EarlyReflectionConfig earlyConfig_{};
  Convolver convolver_{};
  LateReverb late_{};
  std::array<LiveFractionalDelay, MaximumSources> preDelays_{};
  std::array<LiveFractionalDelay, MaximumSources> lateAlignmentDelays_{};
  std::array<float, 2> highCutState_{};
  std::array<float, 2> lowCutState_{};
  std::unique_ptr<EarlyReflectionWorker> worker_{};
  SubmittedScene submitted_{};
  bool hasSubmittedScene_{};
  std::size_t requestCountdown_{};
  std::size_t requestIntervalSamples_{2'400};
  AutomaticBalanceMailbox balanceMailbox_{};
  std::size_t appliedBalanceSequence_{};
  InputFrame currentTailSend_{};
  InputFrame fromTailSend_{};
  InputFrame targetTailSend_{};
  InputFrame handoffStartSeconds_{};
  InputFrame currentLateAlignmentSamples_{};
  InputFrame fromLateAlignmentSamples_{};
  InputFrame toLateAlignmentSamples_{};
  InputFrame pendingLateAlignmentSamples_{};
  InputFrame lateAlignmentTransitionPhase_{};
  std::array<bool, MaximumSources> lateAlignmentInitialized_{};
  float lateAlignmentTransitionIncrement_{1.f};
  float balanceTransitionPhase_{1.f};
  float balanceTransitionIncrement_{1.f};
  float wetSizeGain_{1.f};
  float targetWetSizeGain_{1.f};
  bool wetSizeGainInitialized_{};
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
  using DirectGains =
      std::array<StereoFrame, MaximumSources>;
  DirectGains directGains_{};
  DirectGains directGainSteps_{};
  bool directGainsInitialized_{};

