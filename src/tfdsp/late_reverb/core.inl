public:
  static constexpr std::size_t FdnLineCount =
      late_reverb_coefficients::LineCount;
  static constexpr std::size_t WallCount = late_reverb_coefficients::WallCount;
  using StereoFrame = std::array<float, 2>;
  using LineFrame = std::array<float, FdnLineCount>;
  using WallFrame = std::array<float, WallCount>;

private:
  static constexpr float Pi = 3.14159265358979323846f;
  static constexpr float SpeedOfSound = 343.f;
  static constexpr float MaximumFdnModulationSeconds = 0.000500f;
  // Shimmer adds a bounded octave return to four orthogonal coordinates inside
  // the main velvet loop without replacing the unshifted room feedback. The
  // remaining twelve coordinates are unaltered.
  static constexpr float MaximumShimmerLoopGain = 0.85f;
  static constexpr float MainDelayTransitionSeconds = 0.050f;
  static constexpr std::size_t ControlUpdateInterval = 64;

  double sampleRate_{48'000.0};
  CubicFractionalDelayBank<FdnLineCount> mainDelays_{};
  VelvetFeedbackMatrix feedbackMatrix_{};
  std::array<SmoothRandomModulator, FdnLineCount> modulators_{};
  MultibandDecayFilterBank<FdnLineCount> mainDecayFilters_{};
  std::array<WindowedPitchShifter, late_reverb_coefficients::ShimmerBusCount>
      shimmerShifters_{};
  std::array<float, late_reverb_coefficients::ShimmerBusCount>
      shimmerHighpassState_{};
  std::array<std::array<float, late_reverb_coefficients::ShimmerBusCount>, 2>
      shimmerLowpassState_{};
  float shimmerHighpassAlpha_{};
  float shimmerLowpassAlpha_{};
  float currentMeanDelaySamples_{480.f};
  float fromMeanDelaySamples_{480.f};
  float toMeanDelaySamples_{480.f};
  float pendingMeanDelaySamples_{480.f};
  float mainDelayPhase_{1.f};
  float mainDelayIncrement_{1.f};
  bool mainDelayInitialized_{};
  LateMainDelayRatios mainRatioFrom_ =
      LateMainRatios(DefaultLateReverbFlavour);
  LateMainDelayRatios mainRatioTarget_ =
      LateMainRatios(DefaultLateReverbFlavour);
  float mainRatioTransitionPhase_{1.f};
  float mainRatioTransitionIncrement_{1.f};
  LateReverbFlavour flavour_{DefaultLateReverbFlavour};
  std::size_t controlCountdown_{};
  float controlRoomScale_{1.f};
  float controlModulationDepth_{};
  float controlModulationAmount_{};
  float controlDecaySeconds_{1.f};
  float controlHighT60_{1.f};
  float controlLowT60_{1.f};
  float controlDiffusion_{1.f};
  float controlShimmerAmount_{};

  static float ClampControl(const float value) noexcept {
    return std::clamp(std::isfinite(value) ? value : 0.f, 0.f, 1.f);
  }
  static float Smoothstep(const float value) noexcept {
    const float x = ClampControl(value);
    return x * x * (3.f - 2.f * x);
  }

  void UpdateControlState(const LateReverbControls &controls) noexcept {
    const float meanFreeTime =
        MeanFreeTimeSeconds(controls.roomDimensionsMetres);
    UpdateMeanDelay(meanFreeTime * static_cast<float>(sampleRate_));
    controlRoomScale_ =
        meanFreeTime /
        MeanFreeTimeSeconds(reverb_defaults::RoomDimensionsMetres);
    const float modulation = ClampControl(controls.modulation);
    controlModulationAmount_ = modulation * modulation;
    controlModulationDepth_ = controlModulationAmount_ *
                              MaximumFdnModulationSeconds *
                              static_cast<float>(sampleRate_);
    const float damping = Smoothstep(controls.damping);
    controlDecaySeconds_ =
        0.25f * std::exp(ClampControl(controls.decay) * std::log(32.f));
    controlHighT60_ = controlDecaySeconds_ * std::pow(0.20f, damping);
    controlLowT60_ = controlDecaySeconds_ * std::pow(0.72f, damping);
    controlDiffusion_ = ClampControl(controls.diffusion);
    controlShimmerAmount_ = Smoothstep(controls.shimmer);
  }

  void PollControlState(const LateReverbControls &controls) noexcept {
    if (controlCountdown_ == 0) {
      UpdateControlState(controls);
      controlCountdown_ = ControlUpdateInterval;
    }
    --controlCountdown_;
  }

  void UpdateMeanDelay(const float targetSamples) noexcept {
    if (!mainDelayInitialized_) {
      currentMeanDelaySamples_ = fromMeanDelaySamples_ = toMeanDelaySamples_ =
          pendingMeanDelaySamples_ = targetSamples;
      mainDelayPhase_ = 1.f;
      mainDelayInitialized_ = true;
      return;
    }
    pendingMeanDelaySamples_ = targetSamples;
    if (mainDelayPhase_ < 1.f ||
        std::abs(targetSamples - currentMeanDelaySamples_) < 0.25f)
      return;
    fromMeanDelaySamples_ = currentMeanDelaySamples_;
    toMeanDelaySamples_ = targetSamples;
    mainDelayPhase_ = 0.f;
  }

  LineFrame ReadMainDelaysForRatios(
      const LateMainDelayRatios &ratios,
      const LineFrame &modulation) const noexcept {
    LineFrame delaySamples{};
    for (std::size_t line = 0; line < FdnLineCount; ++line)
      delaySamples[line] =
          ratios[line] * currentMeanDelaySamples_ + modulation[line];
    if (mainDelayPhase_ >= 1.f)
      return mainDelays_.Read(delaySamples);
    LineFrame fromSamples{};
    LineFrame toSamples{};
    for (std::size_t line = 0; line < FdnLineCount; ++line) {
      fromSamples[line] =
          ratios[line] * fromMeanDelaySamples_ + modulation[line];
      toSamples[line] = ratios[line] * toMeanDelaySamples_ + modulation[line];
    }
    const float fade =
        mainDelayPhase_ * mainDelayPhase_ * (3.f - 2.f * mainDelayPhase_);
    const auto from = mainDelays_.Read(fromSamples);
    const auto to = mainDelays_.Read(toSamples);
    LineFrame output{};
    for (std::size_t line = 0; line < FdnLineCount; ++line)
      output[line] = from[line] + fade * (to[line] - from[line]);
    return output;
  }

  LineFrame ReadMainDelays(const LineFrame &modulation) const noexcept {
    if (mainRatioTransitionPhase_ >= 1.f)
      return ReadMainDelaysForRatios(mainRatioTarget_, modulation);
    const float phase = mainRatioTransitionPhase_;
    const float fade = phase * phase * (3.f - 2.f * phase);
    const auto from = ReadMainDelaysForRatios(mainRatioFrom_, modulation);
    const auto to = ReadMainDelaysForRatios(mainRatioTarget_, modulation);
    LineFrame output{};
    for (std::size_t line = 0; line < FdnLineCount; ++line)
      output[line] = from[line] + fade * (to[line] - from[line]);
    return output;
  }

  static float
  MeanFreeTimeSeconds(const std::array<float, 3> &dimensions) noexcept {
    const float length = std::clamp(dimensions[0], 1.f, 40.f);
    const float width = std::clamp(dimensions[1], 1.f, 40.f);
    const float height = std::clamp(dimensions[2], 1.f, 40.f);
    const float volume = length * width * height;
    const float surface =
        2.f * (length * width + length * height + width * height);
    return std::clamp(4.f * volume / (surface * SpeedOfSound), 0.003f, 0.078f);
  }

  static std::array<float, FdnLineCount>
  ProjectWalls(const std::array<float, WallCount> &walls) noexcept {
    std::array<float, FdnLineCount> frame{};
    for (std::size_t line = 0; line < FdnLineCount; ++line)
      for (std::size_t wall = 0; wall < WallCount; ++wall)
        frame[line] += late_reverb_coefficients::InverseRootLineCount *
                       late_reverb_coefficients::WallProjection[line][wall] *
                       walls[wall];
    return frame;
  }

  LineFrame ProcessFdnLoop(const LineFrame &injection) noexcept {
    LineFrame modulation{};
    for (std::size_t line = 0; line < FdnLineCount; ++line)
      modulation[line] = controlModulationDepth_ * modulators_[line].Next();
    const LineFrame delayed = ReadMainDelays(modulation);

    const float attenuationMeanSamples =
        mainDelayPhase_ < 1.f
            ? std::max(fromMeanDelaySamples_, toMeanDelaySamples_)
            : currentMeanDelaySamples_;
    std::array<float, FdnLineCount> pathSeconds{};
    for (std::size_t line = 0; line < FdnLineCount; ++line) {
      const float attenuationRatio =
          mainRatioTransitionPhase_ < 1.f
              ? std::max(mainRatioFrom_[line], mainRatioTarget_[line])
              : mainRatioTarget_[line];
      pathSeconds[line] = attenuationRatio * attenuationMeanSamples /
                          static_cast<float>(sampleRate_);
    }
    auto attenuated = mainDecayFilters_.Process(
        delayed, pathSeconds, controlLowT60_, controlDecaySeconds_,
        controlHighT60_);

    // The octave branch belongs inside the production velvet loop. Projecting
    // a bounded return onto four orthonormal coordinates and then running the
    // complete 16-line VFM avoids a direct grainy output layer and a separate
    // recursive tank. At zero gain the shifters advance their filters, input
    // history and grain state but omit window evaluation and interpolated
    // reads, so automation remains immediate without paying the full render
    // cost for an inaudible branch.
    std::array<float, late_reverb_coefficients::ShimmerBusCount> shimmerBus{};
    std::array<float, late_reverb_coefficients::ShimmerBusCount> shiftedBus{};
    const bool renderShimmer = controlShimmerAmount_ > 1.e-6f;
    for (std::size_t bus = 0; bus < shimmerBus.size(); ++bus) {
      for (std::size_t line = 0; line < FdnLineCount; ++line)
        shimmerBus[bus] +=
            late_reverb_coefficients::ShimmerProjection[bus][line] *
            attenuated[line];
      shimmerHighpassState_[bus] +=
          shimmerHighpassAlpha_ *
          (shimmerBus[bus] - shimmerHighpassState_[bus]);
      shimmerHighpassState_[bus] =
          FiniteNormalOrZero(shimmerHighpassState_[bus]);
      const float highpassed = shimmerBus[bus] - shimmerHighpassState_[bus];
      float shifted = 0.f;
      if (renderShimmer)
        shifted = shimmerShifters_[bus].Process(highpassed);
      else
        shimmerShifters_[bus].Advance(highpassed);
      for (auto &stage : shimmerLowpassState_) {
        stage[bus] += shimmerLowpassAlpha_ * (shifted - stage[bus]);
        stage[bus] = FiniteNormalOrZero(stage[bus]);
        shifted = stage[bus];
      }
      shiftedBus[bus] = FiniteNormalOrZero(shifted);
    }

    const float shimmerGain =
        MaximumShimmerLoopGain * controlShimmerAmount_;
    if (shimmerGain > 1.e-6f) {
      for (std::size_t line = 0; line < FdnLineCount; ++line)
        for (std::size_t bus = 0; bus < shimmerBus.size(); ++bus)
          attenuated[line] +=
              late_reverb_coefficients::ShimmerProjection[bus][line] *
              (shimmerGain * shiftedBus[bus]);
    }

    auto feedback = feedbackMatrix_.Process(
        attenuated, controlDiffusion_, controlDecaySeconds_, controlHighT60_,
        controlLowT60_, controlRoomScale_, controlModulationAmount_);

    for (std::size_t line = 0; line < FdnLineCount; ++line)
      feedback[line] += .42f * injection[line];
    mainDelays_.Push(feedback);

    if (mainDelayPhase_ < 1.f) {
      mainDelayPhase_ = std::min(1.f, mainDelayPhase_ + mainDelayIncrement_);
      if (mainDelayPhase_ >= 1.f) {
        currentMeanDelaySamples_ = toMeanDelaySamples_;
        if (std::abs(pendingMeanDelaySamples_ - currentMeanDelaySamples_) >=
            0.25f) {
          fromMeanDelaySamples_ = currentMeanDelaySamples_;
          toMeanDelaySamples_ = pendingMeanDelaySamples_;
          mainDelayPhase_ = 0.f;
        }
      }
    }
    if (mainRatioTransitionPhase_ < 1.f) {
      mainRatioTransitionPhase_ =
          std::min(1.f, mainRatioTransitionPhase_ +
                            mainRatioTransitionIncrement_);
      if (mainRatioTransitionPhase_ >= 1.f)
        mainRatioFrom_ = mainRatioTarget_;
    }
    for (const auto value : delayed)
      if (!std::isfinite(value)) {
        Reset();
        return {};
      }
    return delayed;
  }

  WallFrame ProcessWallLoop(const WallFrame &inputWalls) noexcept {
    const auto delayed = ProcessFdnLoop(ProjectWalls(inputWalls));
    WallFrame outputWalls{};
    for (std::size_t wall = 0; wall < WallCount; ++wall)
      for (std::size_t line = 0; line < FdnLineCount; ++line)
        outputWalls[wall] +=
            late_reverb_coefficients::InverseRootLineCount *
            late_reverb_coefficients::WallProjection[line][wall] *
            delayed[line];
    return outputWalls;
  }
