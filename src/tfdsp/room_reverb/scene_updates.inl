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
    controlRoom_ = MakeRoomFromSanitizedControls(controls);
    UpdateDirectGains(positions, sourceCount);
    UpdateLateAlignmentTargets();
    targetWetSizeGain_ = WetSizeCalibration(controlRoom_);
    if (!wetSizeGainInitialized_) {
      wetSizeGain_ = targetWetSizeGain_;
      wetSizeGainInitialized_ = true;
    }
    const float balance = ClampControl(controls.balance);
    earlyGain_ =
        EarlyCalibrationGain * std::min(1.f, 2.f * (1.f - balance));
    tailGain_ = std::min(1.f, 2.f * balance);

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
    request.room = MakeRoomFromSanitizedControls(controls);
    request.sourceCount = sourceCount;
    for (std::size_t source = 0; source < sourceCount; ++source) {
      std::array<double, 3> normalized{};
      for (std::size_t axis = 0; axis < 3; ++axis)
        normalized[axis] = ClampControl(positions[source][axis]);
      request.sources[source] =
          detail::MakeEarlyReflectionSourceFromValidatedInputs(request.room,
                                                               normalized);
    }
    request.materials = detail::MakeEarlyReflectionMaterialsFromValidatedInput(
        ClampControl(controls.damping));
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
      highCutState_[channel] = FiniteNormalOrZero(highCutState_[channel]);
      wet[channel] = highCutState_[channel];
      lowCutState_[channel] +=
          lowCutAlpha_ * (wet[channel] - lowCutState_[channel]);
      lowCutState_[channel] = FiniteNormalOrZero(lowCutState_[channel]);
      wet[channel] -= lowCutState_[channel];
      wet[channel] = FiniteNormalOrZero(wet[channel]);
    }
    return wet;
  }

