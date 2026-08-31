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
    // A source not present in this scene has no accepted handoff yet. Keep its
    // late send closed until a matching FIR/metadata scene is adopted.
    balance.tailSend.fill(0.f);
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
      balance.handoffStartSeconds[source] =
          static_cast<float>(handoff.finalStartSeconds);
    }
    return balance;
  }

  void SetLateAlignmentTarget(const std::size_t source,
                              const float samples) noexcept {
    const float target = std::max(0.f, samples);
    if (!lateAlignmentInitialized_[source]) {
      currentLateAlignmentSamples_[source] =
          fromLateAlignmentSamples_[source] =
              toLateAlignmentSamples_[source] =
                  pendingLateAlignmentSamples_[source] = target;
      lateAlignmentTransitionPhase_[source] = 1.f;
      lateAlignmentInitialized_[source] = true;
      return;
    }
    pendingLateAlignmentSamples_[source] = target;
    if (lateAlignmentTransitionPhase_[source] < 1.f ||
        std::abs(target - currentLateAlignmentSamples_[source]) < 0.25f)
      return;
    fromLateAlignmentSamples_[source] = currentLateAlignmentSamples_[source];
    toLateAlignmentSamples_[source] = target;
    lateAlignmentTransitionPhase_[source] = 0.f;
  }

  void UpdateLateAlignmentTargets() noexcept {
    std::array<float, 3> dimensions{};
    for (std::size_t axis = 0; axis < dimensions.size(); ++axis)
      dimensions[axis] =
          static_cast<float>(controlRoom_.dimensionsMetres[axis]);
    const float intrinsic = late_.IntrinsicFirstOutputSeconds(dimensions);
    for (std::size_t source = 0; source < MaximumSources; ++source) {
      if (handoffStartSeconds_[source] <= 0.f &&
          !lateAlignmentInitialized_[source])
        continue;
      SetLateAlignmentTarget(
          source, std::max(0.f, handoffStartSeconds_[source] - intrinsic) *
                      static_cast<float>(sampleRate_));
    }
  }

  float ProcessLateAlignment(const std::size_t source,
                             const float input) noexcept {
    const auto read = [&](const float delay) {
      return lateAlignmentDelays_[source].Read(delay, input);
    };
    float output = read(currentLateAlignmentSamples_[source]);
    auto &phase = lateAlignmentTransitionPhase_[source];
    if (phase < 1.f) {
      const float fade = phase * phase * (3.f - 2.f * phase);
      output = read(fromLateAlignmentSamples_[source]) +
               fade * (read(toLateAlignmentSamples_[source]) -
                       read(fromLateAlignmentSamples_[source]));
      phase = std::min(1.f, phase + lateAlignmentTransitionIncrement_);
      if (phase >= 1.f) {
        currentLateAlignmentSamples_[source] =
            toLateAlignmentSamples_[source];
        if (std::abs(pendingLateAlignmentSamples_[source] -
                     currentLateAlignmentSamples_[source]) >= 0.25f) {
          fromLateAlignmentSamples_[source] =
              currentLateAlignmentSamples_[source];
          toLateAlignmentSamples_[source] =
              pendingLateAlignmentSamples_[source];
          phase = 0.f;
        }
      }
    }
    lateAlignmentDelays_[source].Push(input);
    return output;
  }

  void UpdateAutomaticBalance() noexcept {
    AutomaticBalance latest;
    const auto activeScene = convolver_.RenderedSceneSequence();
    if (activeScene > appliedBalanceSequence_ &&
        balanceMailbox_.ConsumeForScene(latest, activeScene)) {
      const bool firstAcceptedScene = appliedBalanceSequence_ == 0;
      appliedBalanceSequence_ = activeScene;
      targetTailSend_ = latest.tailSend;
      handoffStartSeconds_ = latest.handoffStartSeconds;
      UpdateLateAlignmentTargets();
      if (firstAcceptedScene) {
        currentTailSend_ = fromTailSend_ = targetTailSend_;
        balanceTransitionPhase_ = 1.f;
      } else {
        fromTailSend_ = currentTailSend_;
        balanceTransitionPhase_ = 0.f;
      }
    }
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

