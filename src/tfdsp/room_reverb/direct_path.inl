  static float ClampControl(const float value) noexcept {
    return std::clamp(std::isfinite(value) ? value : 0.f, 0.f, 1.f);
  }

  static EarlyReflectionRoom
  MakeRoomFromSanitizedControls(const RoomReverbControls &controls) noexcept {
    auto room = detail::MakeEarlyReflectionRoomFromValidatedInputs(
        ClampControl(controls.space), EarlyReflectionRoomFamily::Room);
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
    return room;
  }

  static float Smoothstep(const float value) noexcept {
    const float limited = ClampControl(value);
    return limited * limited * (3.f - 2.f * limited);
  }

  static DirectGains CalculateDirectGains(const SourcePositions &positions,
                                          const std::size_t sourceCount,
                                          const EarlyReflectionRoom &room) noexcept {
    DirectGains gains{};
    const std::size_t active = std::min(sourceCount, MaximumSources);
    const double roomScale = std::cbrt(
        room.dimensionsMetres[0] * room.dimensionsMetres[1] *
        room.dimensionsMetres[2]);
    // Factory source/listener distance divided by the cube-root factory room
    // volume. Derive it from the authoritative scene so changing the factory
    // room cannot silently invalidate positioned-direct unity calibration.
    static const double referenceRho = [] {
      double distanceSquared = 0.0;
      double volume = 1.0;
      for (std::size_t axis = 0; axis < 3; ++axis) {
        const double dimension =
            static_cast<double>(reverb_defaults::RoomDimensionsMetres[axis]);
        const double difference =
            (static_cast<double>(reverb_defaults::Source[axis]) -
             static_cast<double>(reverb_defaults::Listener[axis])) *
            dimension;
        distanceSquared += difference * difference;
        volume *= dimension;
      }
      return std::sqrt(distanceSquared) / std::cbrt(volume);
    }();
    for (std::size_t source = 0; source < active; ++source) {
      std::array<double, 3> difference{};
      for (std::size_t axis = 0; axis < difference.size(); ++axis) {
        constexpr double sourceMargin = 0.001;
        const double normalizedSource = std::clamp(
            static_cast<double>(ClampControl(positions[source][axis])),
            sourceMargin, 1.0 - sourceMargin);
        const double sourceMetres =
            normalizedSource * room.dimensionsMetres[axis];
        difference[axis] = sourceMetres - room.listenerPositionMetres[axis];
      }
      const double horizontal = std::hypot(difference[0], difference[1]);
      const double lateral = horizontal > 1.0e-12
                                 ? std::clamp(difference[0] / horizontal,
                                              -1.0, 1.0)
                                 : 0.0;
      const double pan = 0.5 * (lateral + 1.0);
      const double distance =
          std::sqrt(difference[0] * difference[0] +
                    difference[1] * difference[1] +
                    difference[2] * difference[2]);
      const double rho = distance / std::max(roomScale, 1.0e-12);
      const float distanceGain = static_cast<float>(std::clamp(
          referenceRho / std::max(rho, 0.05), 0.25, 2.0));
      gains[source][0] = distanceGain * static_cast<float>(
                                           std::cos(0.5 * EarlyReflectionPi * pan));
      gains[source][1] = distanceGain * static_cast<float>(
                                           std::sin(0.5 * EarlyReflectionPi * pan));
    }
    return gains;
  }

  void UpdateDirectGains(const SourcePositions &positions,
                         const std::size_t sourceCount) noexcept {
    const auto targets =
        CalculateDirectGains(positions, sourceCount, controlRoom_);
    if (!directGainsInitialized_) {
      directGains_ = targets;
      directGainSteps_ = {};
      directGainsInitialized_ = true;
      return;
    }
    for (std::size_t source = 0; source < MaximumSources; ++source)
      for (std::size_t channel = 0; channel < 2; ++channel)
        directGainSteps_[source][channel] =
            (targets[source][channel] - directGains_[source][channel]) /
            static_cast<float>(ControlUpdateInterval);
  }

  StereoFrame RenderDirect(const InputFrame &input,
                           const std::size_t sourceCount) noexcept {
    StereoFrame direct{};
    const std::size_t active = std::min(sourceCount, MaximumSources);
    for (std::size_t source = 0; source < active; ++source) {
      const float sample =
          std::isfinite(input[source]) ? input[source] : 0.f;
      direct[0] += sample * directGains_[source][0];
      direct[1] += sample * directGains_[source][1];
    }
    for (std::size_t source = 0; source < MaximumSources; ++source)
      for (std::size_t channel = 0; channel < 2; ++channel)
        directGains_[source][channel] += directGainSteps_[source][channel];
    return direct;
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
    return preDelays_[source].Read(delaySamples, liveSample);
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

