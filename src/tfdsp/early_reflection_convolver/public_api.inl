public:
  void Prepare(const double sampleRate,
               const std::size_t maximumImpulseSamples) {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument(
          "convolver sample rate must be positive and finite");
    if (maximumImpulseSamples == 0)
      throw std::invalid_argument("convolver FIR capacity must be non-zero");
    sampleRate_ = sampleRate;
    maximumImpulseSamples_ = maximumImpulseSamples;
    maximumPartitions_ = std::max<std::size_t>(
        1, maximumImpulseSamples > HeadSize
               ? (maximumImpulseSamples - HeadSize + PartitionSize - 1) /
                     PartitionSize
               : 1);
    fft_.Prepare();
    for (auto &bank : banks_)
      for (auto &kernel : bank.spectra)
        kernel.assign(maximumPartitions_, Spectrum{});
    for (auto &history : inputHistory_)
      history.assign(maximumPartitions_, Spectrum{});
    for (auto &state : bankStates_)
      state.store(BankState::Free, std::memory_order_relaxed);
    readyBankIndex_.store(-1, std::memory_order_relaxed);
    currentBankIndex_ = InvalidBank;
    previousBankIndex_ = InvalidBank;
    activeSceneSequence_ = 0;
    renderedSceneSequence_ = 0;
    prepared_ = true;
    Reset();
  }

  void Reset() noexcept {
    for (auto &history : inputHistory_)
      for (auto &spectrum : history)
        spectrum.fill({});
    for (auto &block : previousInputBlock_)
      block.fill(Sample{});
    for (auto &history : headHistory_)
      history.fill(Sample{});
    for (auto &block : inputBlock_)
      block.fill(Sample{});
    for (auto &block : currentOutputBlock_)
      block.fill(Sample{});
    for (auto &block : previousOutputBlock_)
      block.fill(Sample{});
    blockPosition_ = 0;
    blockSourceCount_ = 0;
    historyWriteIndex_ = 0;
    headWriteIndex_ = 0;
    sourceFlushBlocks_.fill(0);
    transitionSamples_ = 0;
    transitionRemaining_ = 0;
    renderedSceneSequence_ = activeSceneSequence_;
  }

  void SetImpulseResponse(const EarlyReflectionImpulseResponse &response,
                          const std::size_t sceneSequence = 0) {
    if (!prepared_)
      throw std::logic_error("prepare the ER convolver before loading an FIR");
    readyBankIndex_.store(-1, std::memory_order_relaxed);
    for (auto &state : bankStates_)
      state.store(BankState::Free, std::memory_order_relaxed);
    currentBankIndex_ = 0;
    previousBankIndex_ = InvalidBank;
    bankStates_[currentBankIndex_].store(BankState::Preparing,
                                         std::memory_order_relaxed);
    LoadBank(response, currentBankIndex_);
    banks_[currentBankIndex_].sceneSequence = sceneSequence;
    activeSceneSequence_ = sceneSequence;
    bankStates_[currentBankIndex_].store(BankState::Active,
                                         std::memory_order_release);
    Reset();
  }

  bool
  PrepareAndQueueImpulseResponse(const EarlyReflectionImpulseResponse &response,
                                 const double transitionSeconds = 0.100,
                                 const std::size_t sceneSequence = 0) {
    if (!std::isfinite(transitionSeconds) || transitionSeconds < 0.0)
      throw std::invalid_argument(
          "ER transition time must be finite and non-negative");
    if (!prepared_)
      throw std::logic_error("prepare the ER convolver before loading an FIR");

    std::size_t bankIndex = InvalidBank;
    for (std::size_t candidate = 0; candidate < BankCount; ++candidate) {
      BankState expected = BankState::Free;
      if (bankStates_[candidate].compare_exchange_strong(
              expected, BankState::Preparing, std::memory_order_acq_rel,
              std::memory_order_acquire)) {
        bankIndex = candidate;
        break;
      }
    }
    if (bankIndex == InvalidBank)
      return false;

    try {
      LoadBank(response, bankIndex);
      banks_[bankIndex].sceneSequence = sceneSequence;
      banks_[bankIndex].transitionSamples =
          std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(
                                       transitionSeconds * sampleRate_)));
    } catch (...) {
      bankStates_[bankIndex].store(BankState::Free, std::memory_order_release);
      throw;
    }
    bankStates_[bankIndex].store(BankState::Ready, std::memory_order_release);
    const int replaced = readyBankIndex_.exchange(static_cast<int>(bankIndex),
                                                  std::memory_order_acq_rel);
    if (replaced >= 0) {
      BankState expected = BankState::Ready;
      bankStates_[static_cast<std::size_t>(replaced)].compare_exchange_strong(
          expected, BankState::Free, std::memory_order_acq_rel,
          std::memory_order_acquire);
    }
    return true;
  }

  // Scene identity used for the most recently returned audio sample.
  std::size_t RenderedSceneSequence() const noexcept {
    return renderedSceneSequence_;
  }

  void
  TransitionToImpulseResponse(const EarlyReflectionImpulseResponse &response,
                              const double transitionSeconds = 0.100,
                              const std::size_t sceneSequence = 0) {
    if (currentBankIndex_ == InvalidBank) {
      SetImpulseResponse(response, sceneSequence);
      return;
    }
    if (!PrepareAndQueueImpulseResponse(response, transitionSeconds,
                                        sceneSequence))
      throw std::runtime_error("no free ER convolution bank is available");
  }

  OutputFrame Process(const InputFrame &input,
                      const std::size_t sourceCount) noexcept {
    if (!prepared_)
      return {};
    renderedSceneSequence_ = activeSceneSequence_;
    const std::size_t activeSources = std::min(sourceCount, MaximumSources);
    for (std::size_t source = 0; source < MaximumSources; ++source) {
      const Sample sample =
          source < activeSources && std::isfinite(input[source])
              ? FiniteNormalOrZero(input[source])
              : Sample{};
      inputBlock_[source][blockPosition_] = sample;
      headHistory_[source][headWriteIndex_] = sample;
    }
    blockSourceCount_ = std::max(blockSourceCount_, activeSources);

    const OutputFrame currentHead =
        currentBankIndex_ != InvalidBank ? RenderHead(banks_[currentBankIndex_])
                                         : OutputFrame{};
    const OutputFrame previousHead =
        previousBankIndex_ != InvalidBank
            ? RenderHead(banks_[previousBankIndex_])
            : OutputFrame{};
    OutputFrame output{};
    for (std::size_t channel = 0; channel < EarlyReflectionOutputCount;
         ++channel) {
      output[channel] =
          currentHead[channel] + currentOutputBlock_[channel][blockPosition_];
      if (transitionRemaining_ > 0) {
        const std::size_t completed = transitionSamples_ - transitionRemaining_;
        const double linear =
            transitionSamples_ <= 1
                ? 1.0
                : static_cast<double>(completed) /
                      static_cast<double>(transitionSamples_ - 1);
        const Sample amount =
            static_cast<Sample>(EarlyReflectionSmoothstep(linear));
        const Sample previous = previousHead[channel] +
                                previousOutputBlock_[channel][blockPosition_];
        output[channel] = previous +
            amount * (output[channel] -
                      previous);
      }
      output[channel] = FiniteNormalOrZero(output[channel]);
    }

    if (transitionRemaining_ > 0)
      --transitionRemaining_;
    if (++headWriteIndex_ == HeadSize)
      headWriteIndex_ = 0;
    ++blockPosition_;
    if (blockPosition_ == PartitionSize) {
      blockPosition_ = 0;
      ProcessBlock();
    }
    return output;
  }

  bool IsPrepared() const noexcept { return prepared_; }
  bool IsTransitioning() const noexcept {
    return transitionRemaining_ > 0 ||
           readyBankIndex_.load(std::memory_order_acquire) >= 0;
  }
  std::size_t MaximumImpulseSamples() const noexcept {
    return maximumImpulseSamples_;
  }
