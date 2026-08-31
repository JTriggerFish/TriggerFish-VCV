public:
  LateReverb() { SetSampleRate(sampleRate_); }

  void SetSampleRate(const double sampleRate) {
    if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
      throw std::invalid_argument("late reverb sample rate must be positive");
    sampleRate_ = sampleRate;
    const auto &baseRatios = LateMainRatios(LateReverbFlavour::Base);
    const auto &optimizedRatios = LateMainRatios(LateReverbFlavour::Optimized);
    const float maximumRatio =
        std::max(*std::max_element(baseRatios.begin(), baseRatios.end()),
                 *std::max_element(optimizedRatios.begin(),
                                   optimizedRatios.end()));
    const auto mainCapacity = static_cast<std::size_t>(std::ceil(
        (0.078 * maximumRatio + MaximumFdnModulationSeconds) * sampleRate_ +
        8.0));
    for (std::size_t line = 0; line < FdnLineCount; ++line) {
      mainDelays_[line].Prepare(mainCapacity);
      mainDecayFilters_[line].Prepare(sampleRate_);
    }
    feedbackMatrix_.Prepare(sampleRate_);
    for (std::size_t line = 0; line < FdnLineCount; ++line)
      modulators_[line].Prepare(
          sampleRate_, late_reverb_coefficients::ModulationRateHz[line],
          0x51f15e5du + static_cast<std::uint32_t>(line) * 0x9e3779b9u);
    constexpr std::array<std::uint32_t,
                         late_reverb_coefficients::ShimmerBusCount>
        ShimmerSeeds{0x8f3a21d7u, 0xc4b82e39u, 0x71d934abu, 0xe29c6f15u};
    for (std::size_t bus = 0; bus < shimmerShifters_.size(); ++bus)
      shimmerShifters_[bus].Prepare(
          sampleRate_, 0.120f,
          static_cast<float>(bus) / static_cast<float>(shimmerShifters_.size()),
          ShimmerSeeds[bus]);
    shimmerHighpassAlpha_ =
        1.f - std::exp(-2.f * Pi * 250.f / static_cast<float>(sampleRate_));
    const float shimmerLowpassCutoff =
        std::min(6'500.f, 0.20f * static_cast<float>(sampleRate_));
    shimmerLowpassAlpha_ = 1.f -
                           std::exp(-2.f * Pi * shimmerLowpassCutoff /
                                    static_cast<float>(sampleRate_));
    mainDelayIncrement_ =
        1.f / std::max(1.f, MainDelayTransitionSeconds *
                                static_cast<float>(sampleRate_));
    mainRatioTransitionIncrement_ = mainDelayIncrement_;
    Reset();
  }

  void Reset() noexcept {
    for (auto &delay : mainDelays_)
      delay.Reset();
    feedbackMatrix_.Reset();
    for (auto &modulator : modulators_)
      modulator.Reset();
    for (auto &filter : mainDecayFilters_)
      filter.Reset();
    for (auto &shifter : shimmerShifters_)
      shifter.Reset();
    shimmerHighpassState_.fill(0.f);
    for (auto &stage : shimmerLowpassState_)
      stage.fill(0.f);
    currentMeanDelaySamples_ = fromMeanDelaySamples_ = toMeanDelaySamples_ =
        pendingMeanDelaySamples_ = static_cast<float>(0.010 * sampleRate_);
    mainDelayPhase_ = 1.f;
    mainDelayInitialized_ = false;
    mainRatioFrom_ = mainRatioTarget_ = LateMainRatios(flavour_);
    mainRatioTransitionPhase_ = 1.f;
    controlCountdown_ = 0;
    controlRoomScale_ = 1.f;
    controlModulationDepth_ = 0.f;
    controlModulationAmount_ = 0.f;
    controlDecaySeconds_ = controlHighT60_ = controlLowT60_ = 1.f;
    controlDiffusion_ = 1.f;
    controlShimmerAmount_ = 0.f;
  }

  double SampleRate() const noexcept { return sampleRate_; }

  void SetFlavour(const LateReverbFlavour flavour) noexcept {
    if (flavour == flavour_)
      return;
    if (mainRatioTransitionPhase_ >= 0.5f)
      mainRatioFrom_ = mainRatioTarget_;
    mainRatioTarget_ = LateMainRatios(flavour);
    mainRatioTransitionPhase_ = 0.f;
    flavour_ = flavour;
    feedbackMatrix_.SetFlavour(flavour);
  }

  // Select a coefficient set without a transition. This is intended for
  // restoring a saved flavour before the first audio sample; live changes use
  // SetFlavour() so the stored tail is crossfaded instead of reinterpreted.
  void SetFlavourImmediate(const LateReverbFlavour flavour) noexcept {
    flavour_ = flavour;
    mainRatioFrom_ = mainRatioTarget_ = LateMainRatios(flavour);
    mainRatioTransitionPhase_ = 1.f;
    feedbackMatrix_.SetFlavourImmediate(flavour);
  }

  LateReverbFlavour Flavour() const noexcept { return flavour_; }

  float IntrinsicFirstOutputSeconds(
      const std::array<float, 3> &roomDimensionsMetres) const noexcept {
    const auto &ratios = LateMainRatios(flavour_);
    return MeanFreeTimeSeconds(roomDimensionsMetres) *
           *std::min_element(ratios.begin(), ratios.end());
  }

  // Static wall-domain probe retained for differentiable-reference parity.
  // Production stereo rendering uses the fixed line-domain vectors below.
  WallFrame ProcessWallFrame(const WallFrame &inputWalls,
                             const LateReverbControls &controls) noexcept {
    PollControlState(controls);
    return ProcessWallLoop(inputWalls);
  }

  StereoFrame Process(const float input,
                      const LateReverbControls &controls) noexcept {
    PollControlState(controls);
    LineFrame injection{};
    const float safeInput = FiniteNormalOrZero(input);
    for (std::size_t line = 0; line < FdnLineCount; ++line)
      injection[line] =
          safeInput * late_reverb_coefficients::FixedInputVector[line];
    const auto delayed = ProcessFdnLoop(injection);
    StereoFrame stereo{};
    for (std::size_t channel = 0; channel < stereo.size(); ++channel)
      for (std::size_t line = 0; line < FdnLineCount; ++line)
        stereo[channel] +=
            late_reverb_coefficients::FixedStereoOutput[channel][line] *
            delayed[line];
    for (auto &value : stereo)
      value *= late_reverb_coefficients::TailOutputCalibration;

    if (!std::isfinite(stereo[0]) || !std::isfinite(stereo[1])) {
      Reset();
      return {};
    }
    return stereo;
  }
