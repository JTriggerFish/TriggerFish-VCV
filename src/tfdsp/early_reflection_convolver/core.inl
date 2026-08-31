  static_assert(std::is_floating_point_v<Sample>,
                "EarlyReflectionConvolver requires floating-point samples");
  static_assert(PartitionSize >= 16 &&
                    (PartitionSize & (PartitionSize - 1)) == 0,
                "ER convolution partition size must be a power of two");
  static_assert(MaximumSources >= 1 &&
                    MaximumSources <= EarlyReflectionMaximumSources,
                "ER convolver source capacity must lie in [1, 8]");
  static constexpr std::size_t FftSize = 2 * PartitionSize;
  static constexpr std::size_t MaximumKernelCount =
      EarlyReflectionOutputCount * MaximumSources;
  static constexpr std::size_t BankCount = 4;
  static constexpr std::size_t InvalidBank = BankCount;

public:
  using InputFrame = std::array<Sample, MaximumSources>;
  using OutputFrame = std::array<Sample, EarlyReflectionOutputCount>;
  // The direct-form head is sample synchronous. The partitioned tail retains
  // an internal PartitionSize latency and stores coefficients shifted left by
  // that amount, so the complete convolver has zero implementation latency.
  static constexpr std::size_t LatencySamples = 0;
  static constexpr std::size_t HeadSize = PartitionSize;

private:
  using Fft = FixedRadix2Fft<Sample, FftSize>;
  using Spectrum = typename Fft::Spectrum;

  struct KernelBank {
    std::size_t sourceCount{};
    std::size_t partitionCount{};
    std::size_t transitionSamples{};
    std::size_t sceneSequence{};
    std::array<std::array<Sample, HeadSize>, MaximumKernelCount> head{};
    std::array<std::size_t, MaximumKernelCount> headTapCount{};
    std::array<std::vector<Spectrum>, MaximumKernelCount> spectra{};
  };

  enum class BankState : std::uint8_t { Free, Preparing, Ready, Active };

  double sampleRate_{};
  std::size_t maximumImpulseSamples_{};
  std::size_t maximumPartitions_{};
  std::size_t blockPosition_{};
  std::size_t blockSourceCount_{};
  std::size_t historyWriteIndex_{};
  std::size_t currentBankIndex_{InvalidBank};
  std::size_t previousBankIndex_{InvalidBank};
  std::size_t transitionSamples_{};
  std::size_t transitionRemaining_{};
  std::size_t activeSceneSequence_{};
  std::size_t renderedSceneSequence_{};
  bool prepared_{};
  Fft fft_{};
  std::array<KernelBank, BankCount> banks_{};
  std::array<std::atomic<BankState>, BankCount> bankStates_{};
  std::atomic<int> readyBankIndex_{-1};
  std::array<std::vector<Spectrum>, MaximumSources> inputHistory_{};
  std::array<std::array<Sample, HeadSize>, MaximumSources> headHistory_{};
  std::size_t headWriteIndex_{};
  std::array<std::size_t, MaximumSources> sourceFlushBlocks_{};
  std::array<std::array<Sample, PartitionSize>, MaximumSources>
      previousInputBlock_{};
  std::array<std::array<Sample, PartitionSize>, MaximumSources> inputBlock_{};
  std::array<std::array<Sample, PartitionSize>, EarlyReflectionOutputCount>
      currentOutputBlock_{};
  std::array<std::array<Sample, PartitionSize>, EarlyReflectionOutputCount>
      previousOutputBlock_{};

  void LoadBank(const EarlyReflectionImpulseResponse &response,
                const std::size_t bankIndex) {
    response.Validate();
    if (!prepared_)
      throw std::logic_error("prepare the ER convolver before loading an FIR");
    if (std::abs(response.sampleRate - sampleRate_) > 1.0e-9)
      throw std::invalid_argument(
          "ER FIR sample rate does not match the convolver");
    if (response.sourceCount > MaximumSources)
      throw std::invalid_argument(
          "ER FIR exceeds the convolver source capacity");
    if (response.Size() > maximumImpulseSamples_)
      throw std::invalid_argument(
          "ER FIR exceeds the prepared convolution capacity");
    auto &bank = banks_[bankIndex];
    bank.sourceCount = response.sourceCount;
    bank.partitionCount =
        response.Size() > HeadSize
            ? (response.Size() - HeadSize + PartitionSize - 1) / PartitionSize
            : 0;
    for (auto &head : bank.head)
      head.fill(Sample{});
    bank.headTapCount.fill(0);
    for (std::size_t output = 0; output < EarlyReflectionOutputCount; ++output)
      for (std::size_t source = 0; source < response.sourceCount; ++source) {
        const std::size_t kernelIndex = output * MaximumSources + source;
        const auto &impulse =
            response.kernels[response.KernelIndex(output, source)];
        for (std::size_t sample = 0;
             sample < std::min(HeadSize, impulse.size()); ++sample) {
          bank.head[kernelIndex][sample] =
              FiniteNormalOrZero(static_cast<Sample>(impulse[sample]));
          if (std::abs(impulse[sample]) > 1.0e-20)
            bank.headTapCount[kernelIndex] = sample + 1;
        }
        for (std::size_t partition = 0; partition < bank.partitionCount;
             ++partition) {
          auto &spectrum = bank.spectra[kernelIndex][partition];
          spectrum.fill({});
          for (std::size_t sample = 0; sample < PartitionSize; ++sample) {
            const std::size_t impulseIndex =
                HeadSize + partition * PartitionSize + sample;
            if (impulseIndex < impulse.size())
              spectrum[sample] = FiniteNormalOrZero(
                  static_cast<Sample>(impulse[impulseIndex]));
          }
          fft_.Transform(spectrum, false);
        }
      }
  }

  OutputFrame RenderHead(const KernelBank &bank) const noexcept {
    OutputFrame output{};
    for (std::size_t channel = 0; channel < EarlyReflectionOutputCount;
         ++channel)
      for (std::size_t source = 0; source < bank.sourceCount; ++source) {
        const auto &kernel = bank.head[channel * MaximumSources + source];
        const std::size_t tapCount =
            bank.headTapCount[channel * MaximumSources + source];
        for (std::size_t tap = 0; tap < tapCount; ++tap) {
          const std::size_t historyIndex =
              (headWriteIndex_ + HeadSize - tap) % HeadSize;
          output[channel] += kernel[tap] * headHistory_[source][historyIndex];
        }
      }
    return output;
  }

  void RenderBank(const KernelBank &bank,
                  std::array<std::array<Sample, PartitionSize>,
                             EarlyReflectionOutputCount> &output) noexcept {
    for (std::size_t channel = 0; channel < EarlyReflectionOutputCount;
         ++channel) {
      Spectrum accumulated{};
      for (std::size_t source = 0; source < bank.sourceCount; ++source) {
        if (sourceFlushBlocks_[source] == 0)
          continue;
        const std::size_t kernel = channel * MaximumSources + source;
        for (std::size_t partition = 0; partition < bank.partitionCount;
             ++partition) {
          const std::size_t historyIndex =
              (historyWriteIndex_ + maximumPartitions_ - partition) %
              maximumPartitions_;
          const auto &input = inputHistory_[source][historyIndex];
          const auto &filter = bank.spectra[kernel][partition];
          for (std::size_t bin = 0; bin < FftSize; ++bin)
            accumulated[bin] += input[bin] * filter[bin];
        }
      }
      fft_.Transform(accumulated, true);
      for (std::size_t sample = 0; sample < PartitionSize; ++sample)
        output[channel][sample] = FiniteNormalOrZero(
            accumulated[PartitionSize + sample].real());
    }
  }

  void ReleasePreviousBank() noexcept {
    if (previousBankIndex_ == InvalidBank)
      return;
    bankStates_[previousBankIndex_].store(BankState::Free,
                                          std::memory_order_release);
    previousBankIndex_ = InvalidBank;
  }

  void AdoptLatestReadyBank() noexcept {
    if (previousBankIndex_ != InvalidBank)
      return;
    const int candidate =
        readyBankIndex_.exchange(-1, std::memory_order_acq_rel);
    if (candidate < 0)
      return;
    const auto candidateIndex = static_cast<std::size_t>(candidate);
    BankState expected = BankState::Ready;
    if (!bankStates_[candidateIndex].compare_exchange_strong(
            expected, BankState::Active, std::memory_order_acq_rel,
            std::memory_order_acquire))
      return;

    if (currentBankIndex_ == InvalidBank) {
      currentBankIndex_ = candidateIndex;
      activeSceneSequence_ = banks_[candidateIndex].sceneSequence;
      return;
    }
    previousBankIndex_ = currentBankIndex_;
    currentBankIndex_ = candidateIndex;
    activeSceneSequence_ = banks_[candidateIndex].sceneSequence;
    transitionSamples_ = banks_[candidateIndex].transitionSamples;
    transitionRemaining_ = transitionSamples_;
  }

  void ProcessBlock() noexcept {
    for (std::size_t source = 0; source < MaximumSources; ++source) {
      const bool hasCurrentInput = source < blockSourceCount_;
      if (hasCurrentInput)
        sourceFlushBlocks_[source] = maximumPartitions_ + 1;
      if (!hasCurrentInput && sourceFlushBlocks_[source] == 0) {
        previousInputBlock_[source].fill(Sample{});
        continue;
      }
      Spectrum spectrum{};
      for (std::size_t sample = 0; sample < PartitionSize; ++sample) {
        spectrum[sample] = previousInputBlock_[source][sample];
        spectrum[PartitionSize + sample] = inputBlock_[source][sample];
      }
      fft_.Transform(spectrum, false);
      inputHistory_[source][historyWriteIndex_] = spectrum;
      previousInputBlock_[source] = inputBlock_[source];
      if (!hasCurrentInput)
        --sourceFlushBlocks_[source];
    }
    blockSourceCount_ = 0;
    if (transitionRemaining_ == 0)
      ReleasePreviousBank();
    AdoptLatestReadyBank();
    if (currentBankIndex_ != InvalidBank)
      RenderBank(banks_[currentBankIndex_], currentOutputBlock_);
    else
      for (auto &output : currentOutputBlock_)
        output.fill(Sample{});
    if (previousBankIndex_ != InvalidBank && transitionRemaining_ > 0)
      RenderBank(banks_[previousBankIndex_], previousOutputBlock_);

    ++historyWriteIndex_;
    if (historyWriteIndex_ == maximumPartitions_)
      historyWriteIndex_ = 0;
  }

