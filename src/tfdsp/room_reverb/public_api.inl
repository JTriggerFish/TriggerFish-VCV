public:
  RoomReverb() { SetSampleRate(sampleRate_); }

  ~RoomReverb() {
    if (worker_)
      worker_->Stop();
  }

  RoomReverb(const RoomReverb &) = delete;
  RoomReverb &operator=(const RoomReverb &) = delete;

  static EarlyReflectionRoom MakeRoom(const RoomReverbControls &controls) {
    auto room = MakeRoomFromSanitizedControls(controls);
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
    targetTailSend_.fill(0.f);
    handoffStartSeconds_.fill(0.f);
    appliedBalanceSequence_ = 0;
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
    const auto alignmentCapacity = static_cast<std::size_t>(std::ceil(
                                       earlyConfig_.responseTimeSeconds *
                                       sampleRate_)) +
                                   8;
    for (auto &delay : lateAlignmentDelays_)
      delay.Prepare(alignmentCapacity);
    requestIntervalSamples_ = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::llround(sampleRate_ / 20.0)));
    balanceTransitionIncrement_ =
        1.f / std::max(1.f, 0.100f * static_cast<float>(sampleRate_));
    lateAlignmentTransitionIncrement_ = balanceTransitionIncrement_;
    preDelayTransitionIncrement_ =
        1.f / std::max(1.f, PreDelayTransitionSeconds *
                               static_cast<float>(sampleRate_));
    requestCountdown_ = 0;
    hasSubmittedScene_ = false;
    worker_ = std::make_unique<EarlyReflectionWorker>(
        20.0, [this](const EarlyReflectionImpulseResponse &response,
                     const double transitionSeconds,
                     const std::size_t sceneSequence) {
          if (!convolver_.PrepareAndQueueImpulseResponse(response,
                                                         transitionSeconds,
                                                         sceneSequence))
            return false;
          auto balance = MakeAutomaticBalance(response);
          balance.sceneSequence = sceneSequence;
          return balanceMailbox_.Publish(balance);
        });
    Reset();
  }

  void Reset() noexcept {
    convolver_.Reset();
    late_.Reset();
    for (auto &delay : preDelays_)
      delay.Reset();
    for (auto &delay : lateAlignmentDelays_)
      delay.Reset();
    highCutState_.fill(0.f);
    lowCutState_.fill(0.f);
    // Keep the accepted scene's non-audio metadata paired with the active ER
    // bank. Reset clears delay/filter history, but it must not silently revert
    // handoff gain and onset until some future geometry build happens.
    currentTailSend_ = targetTailSend_;
    fromTailSend_ = targetTailSend_;
    currentLateAlignmentSamples_.fill(0.f);
    fromLateAlignmentSamples_.fill(0.f);
    toLateAlignmentSamples_.fill(0.f);
    pendingLateAlignmentSamples_.fill(0.f);
    lateAlignmentTransitionPhase_.fill(1.f);
    lateAlignmentInitialized_.fill(false);
    UpdateLateAlignmentTargets();
    balanceTransitionPhase_ = 1.f;
    wetSizeGain_ = 1.f;
    wetSizeGainInitialized_ = false;
    targetWetSizeGain_ = 1.f;
    earlyGain_ = tailGain_ = 1.f;
    highCutAlpha_ = lowCutAlpha_ = 1.f;
    controlCountdown_ = 0;
    requestCountdown_ = 0;
    currentPreDelaySamples_ = fromPreDelaySamples_ = toPreDelaySamples_ =
        pendingPreDelaySamples_ = 0.f;
    preDelayTransitionPhase_ = 1.f;
    preDelayInitialized_ = false;
    directGains_ = {};
    directGainSteps_ = {};
    directGainsInitialized_ = false;
  }

  double SampleRate() const noexcept { return sampleRate_; }
  void SetLateReverbFlavour(const LateReverbFlavour flavour) noexcept {
    late_.SetFlavour(flavour);
  }
  void SetLateReverbFlavourImmediate(
      const LateReverbFlavour flavour) noexcept {
    late_.SetFlavourImmediate(flavour);
  }
  LateReverbFlavour LateFlavour() const noexcept { return late_.Flavour(); }
  bool HasSubmittedScene() const noexcept { return hasSubmittedScene_; }

  bool WaitForLatestScene(const std::chrono::milliseconds timeout) {
    if (!worker_)
      return false;
    const auto result = worker_->WaitForLatestResult(timeout);
    return result && result->Succeeded() && result->publishedToConvolver;
  }

  RoomReverbFrame Process(const InputFrame &input,
                          const SourcePositions &positions,
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

    const StereoFrame direct = RenderDirect(input, activeSources);

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

    const auto early = convolver_.Process(delayed, activeSources);
    UpdateAutomaticBalance();
    float lateInput = 0.f;
    for (std::size_t source = 0; source < MaximumSources; ++source) {
      const float aligned = ProcessLateAlignment(source, delayed[source]);
      lateInput += aligned * currentTailSend_[source];
    }
    LateReverbControls lateControls;
    lateControls.decay = controls.decay;
    lateControls.damping = controls.damping;
    lateControls.diffusion = controls.diffusion;
    lateControls.modulation = controls.modulation;
    lateControls.shimmer = controls.shimmer;
    for (std::size_t axis = 0; axis < 3; ++axis)
      lateControls.roomDimensionsMetres[axis] =
          static_cast<float>(controlRoom_.dimensionsMetres[axis]);
    const auto tail = late_.Process(lateInput, lateControls);
    wetSizeGain_ += balanceTransitionIncrement_ *
                    (targetWetSizeGain_ - wetSizeGain_);
    const auto wet = ApplyWidth(
        {wetSizeGain_ * (earlyGain_ * early[0] + tailGain_ * tail[0]),
         wetSizeGain_ * (earlyGain_ * early[1] + tailGain_ * tail[1])},
        controls.width);
    return {direct, FilterWet(wet)};
  }
