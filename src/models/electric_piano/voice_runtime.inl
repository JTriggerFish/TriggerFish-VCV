	double Step(double pitchVolts, double gateVolts, double velocity,
		bool sustain, const ElectricPianoControls& controls,
		const ElectricPianoModulation& modulation = {}, bool retrigger = false)
	{
		if (!std::isfinite(pitchVolts))
			pitchVolts = 0.0;
		if (!std::isfinite(gateVolts))
			gateVolts = 0.0;
		velocity = std::clamp(std::isfinite(velocity) ? velocity : 0.8,
			0.0, 1.0);
		const double exponentialPitch = std::clamp(
			std::isfinite(modulation.exponentialPitch) ?
				modulation.exponentialPitch : 0.0, -4.0, 4.0);
		const double linearFrequencyRatio = std::clamp(
			std::isfinite(modulation.linearFrequencyRatio) ?
				modulation.linearFrequencyRatio : 0.0, -4.0, 4.0);
		const double phaseModulation = std::clamp(
			std::isfinite(modulation.phaseRadians) ?
				modulation.phaseRadians : 0.0, -2.0 * TwoPi, 2.0 * TwoPi);

		const bool gate = gateVolts >= 1.0;
		bool newSilentStrike = false;
		if (gate && (!lastGate_ || retrigger))
		{
			const bool wasAudible = IsAudible();
			if (!wasAudible)
			{
				newSilentStrike = true;
				// A dormant voice is not stepped, so its smoothed controls otherwise
				// remain frozen at the preceding note's values. Re-prime them from the
				// current panel state before this fresh attack; in particular, do not
				// spend the next note's transient gliding old pickup/coupling settings
				// or the slow Decay smoother toward values chosen during silence.
				controlsInitialized_ = false;
				// A silent key has no phase history to preserve. The oversampled
				// modulation controls are primed below, so an already-present DC PM
				// voltage becomes initial phase rather than a one-sample impulse.
				modulationPhasorReal_.fill(1.0);
				modulationPhasorImaginary_.fill(0.0);
				modulationPhasorNormalizationCountdown_ = 4096;
				modulationPathActive_ = false;
				modulationControlsPrimed_ = false;
				exponentialFmPathActive_ = false;
				linearFmPathActive_ = false;
				phaseModulationPathActive_ = false;
				exponentialFmInterpolatorPrimed_ = false;
				linearFmInterpolatorPrimed_ = false;
				phaseModulationInterpolatorPrimed_ = false;
				modulationCrossfade_ = 0.0;
				verticalPositionInterpolator_->Reset();
				verticalVelocityInterpolator_->Reset();
				horizontalPositionInterpolator_->Reset();
				horizontalVelocityInterpolator_->Reset();
				previousHostExponentialPitch_ = exponentialPitch;
				previousHostLinearFrequencyRatio_ = linearFrequencyRatio;
				previousHostPhaseModulation_ = phaseModulation;
			}
			Strike(pitchVolts, velocity, controls);
		}
		else if (!gate && lastGate_)
		{
			keyHeld_ = false;
			keyReleaseNoise_ = 0.035 + 0.075 * latchedVelocity_;
			TriggerKeyReleaseMechanics();
			if (!sustain)
			{
				damperNoise_ = 0.030 + 0.090 * latchedVelocity_;
				TriggerDamperMechanics();
			}
		}
		if (lastSustain_ && !sustain && !keyHeld_)
		{
			damperNoise_ = std::max(damperNoise_,
				0.030 + 0.090 * latchedVelocity_);
			TriggerDamperMechanics();
		}
		lastGate_ = gate;
		lastSustain_ = sustain;
		if (gate)
			keyHeld_ = true;

		const auto activeControls = SmoothActiveControls(controls);
		const double boundedPitch = std::clamp(pitchVolts, -6.0, 6.0);
		bool pitchChanged = false;
		if (boundedPitch != cachedPitch_ || coefficientsDirty_)
		{
			currentFundamental_ = ElectricPianoReferenceFrequency *
				std::exp2(boundedPitch);
			cachedPitch_ = boundedPitch;
			pitchChanged = true;
		}
		const bool coefficientUpdateTick =
			modeCoefficientUpdateCountdown_ <= 0;
		if (coefficientUpdateTick)
			modeCoefficientUpdateCountdown_ = modeCoefficientUpdatePeriod_ - 1;
		else
			--modeCoefficientUpdateCountdown_;
		const bool damped = !keyHeld_ && !sustain;
		RefreshModeCoefficients(currentFundamental_, damped, activeControls,
			coefficientUpdateTick || pitchChanged);
		RefreshTimbreCoefficients(activeControls,
			coefficientUpdateTick || pitchChanged);

		if (exponentialPitch != 0.0)
			exponentialFmPathActive_ = true;
		if (linearFrequencyRatio != 0.0)
			linearFmPathActive_ = true;
		if (phaseModulation != 0.0)
			phaseModulationPathActive_ = true;
		modulationPathActive_ = exponentialFmPathActive_ ||
			linearFmPathActive_ || phaseModulationPathActive_;
		if (newSilentStrike && modulationPathActive_)
			modulationCrossfade_ = 1.0;

		if (contactActive_)
			AdvanceCoupledHammerAndModes();
		else
			AdvanceFreeModes();

		Eigen::Array<double, PickupOversamplingFactor, 1> exponentialPitchFrames;
		Eigen::Array<double, PickupOversamplingFactor, 1> linearRatioFrames;
		Eigen::Array<double, PickupOversamplingFactor, 1> phaseModulationFrames;
		if (modulationPathActive_)
		{
			exponentialPitchFrames.setZero();
			linearRatioFrames.setZero();
			phaseModulationFrames.setZero();
			if (!modulationControlsPrimed_)
			{
				previousOversampledPhaseModulation_ =
					previousHostPhaseModulation_;
				for (std::size_t mode = 0; mode < modes_.size(); ++mode)
				{
					const double initialPhase = ModulationModeWeights[mode] *
						previousOversampledPhaseModulation_;
					modulationPhasorReal_[mode] = std::cos(initialPhase);
					modulationPhasorImaginary_[mode] = std::sin(initialPhase);
				}
				modulationControlsPrimed_ = true;
			}
			if (exponentialFmPathActive_)
			{
				if (!exponentialFmInterpolatorPrimed_)
				{
					exponentialFmInterpolator_->PrimeUpsample(
						previousHostExponentialPitch_);
					exponentialFmInterpolatorPrimed_ = true;
				}
				exponentialPitchFrames = exponentialFmInterpolator_->Upsample(
					exponentialPitch);
			}
			if (linearFmPathActive_)
			{
				if (!linearFmInterpolatorPrimed_)
				{
					linearFmInterpolator_->PrimeUpsample(
						previousHostLinearFrequencyRatio_);
					linearFmInterpolatorPrimed_ = true;
				}
				linearRatioFrames = linearFmInterpolator_->Upsample(
					linearFrequencyRatio);
			}
			if (phaseModulationPathActive_)
			{
				if (!phaseModulationInterpolatorPrimed_)
				{
					phaseModulationInterpolator_->PrimeUpsample(
						previousHostPhaseModulation_);
					phaseModulationInterpolatorPrimed_ = true;
				}
				phaseModulationFrames = phaseModulationInterpolator_->Upsample(
					phaseModulation);
			}
		}
		previousHostExponentialPitch_ = exponentialPitch;
		previousHostLinearFrequencyRatio_ = linearFrequencyRatio;
		previousHostPhaseModulation_ = phaseModulation;

		energy_ = 0.0;
		for (std::size_t index = 0; index < modes_.size(); ++index)
		{
			energy_ += modes_[index].real * modes_[index].real +
				modes_[index].imaginary * modes_[index].imaginary;
		}

		const double white = WhiteNoise();
		mechanicsLowPass_ += mechanicsLowPassCoefficient_ *
			(white - mechanicsLowPass_);
		mechanicsHighPass_ += mechanicsHighPassCoefficient_ *
			(white - mechanicsHighPass_);
		const double middleNoise = mechanicsHighPass_ - mechanicsLowPass_;
		const double hammerColour = (1.0 - hammerNoiseColour_) *
			mechanicsLowPass_ + hammerNoiseColour_ * middleNoise;
		// Key mechanisms are impacts, not continuously running noise sources.
		// These paired, inharmonic resonators restart at a defined phase for
		// every event, keeping chords tight. A little filtered noise prevents
		// the transients from sounding like tiny tuned oscillators.
		const double hammerCarrier = hammerNoise_ > 1.0e-8 ?
			0.26 * AdvanceMechanicalResonator(mechanicalResonators_[0]) +
			0.13 * AdvanceMechanicalResonator(mechanicalResonators_[1]) +
			0.06 * hammerColour : 0.0;
		const double releaseCarrier = keyReleaseNoise_ > 1.0e-8 ?
			0.22 * AdvanceMechanicalResonator(mechanicalResonators_[2]) +
			0.11 * AdvanceMechanicalResonator(mechanicalResonators_[3]) +
			0.05 * mechanicsLowPass_ : 0.0;
		const double damperCarrier = damperNoise_ > 1.0e-8 ?
			0.20 * AdvanceMechanicalResonator(mechanicalResonators_[4]) +
			0.13 * AdvanceMechanicalResonator(mechanicalResonators_[5]) +
			0.07 * middleNoise : 0.0;
		const double mechanicalSignal = mechanicsLevel_ *
			(0.74 * hammerNoise_ * hammerCarrier +
				0.20 * keyReleaseNoise_ * releaseCarrier +
				0.28 * damperNoise_ * damperCarrier);

		RenderPickup(exponentialPitchFrames, linearRatioFrames,
			phaseModulationFrames, mechanicalSignal);
		hammerNoise_ *= hammerNoiseDecay_;
		keyReleaseNoise_ *= keyReleaseNoiseDecay_;
		damperNoise_ *= damperNoiseDecay_;

		if (energy_ < 1.0e-14 && hammerNoise_ < 1.0e-8 &&
			keyReleaseNoise_ < 1.0e-8 && damperNoise_ < 1.0e-8)
		{
			for (auto& mode : modes_)
				mode = {};
			hammerNoise_ = keyReleaseNoise_ = damperNoise_ = 0.0;
		}

		const double result = pickupOutput_;
		return std::isfinite(result) ? result : 0.0;
	}

	double Energy() const { return energy_; }
	double Activity() const
	{
		return energy_ + hammerNoise_ * hammerNoise_ +
			keyReleaseNoise_ * keyReleaseNoise_ +
			damperNoise_ * damperNoise_ +
			pickupOutput_ * pickupOutput_;
	}
	bool IsAudible() const { return Activity() > 1.0e-12; }
	bool GateHigh() const { return lastGate_; }
	double NotePitch() const { return notePitch_; }
	bool ContactActive() const { return contactActive_; }
	double ContactAge() const { return contactAge_; }
	double ModeFrequencyRatio(std::size_t index) const
	{
		return index < modeRatio_.size() ? modeRatio_[index] : 0.0;
	}
	double ModeAmplitudeLifetimeSeconds(std::size_t index) const
	{
		if (index >= modeRadius_.size() || modeRadius_[index] <= 0.0 ||
			modeRadius_[index] >= 1.0)
			return 0.0;
		return -1.0 / (sampleRate_ * std::log(modeRadius_[index]));
	}
	double ModeDisplacementAmplitude(std::size_t index) const
	{
		return index < modes_.size() ?
			std::hypot(modes_[index].real, modes_[index].imaginary) : 0.0;
	}
	double ModePickupDisplacementAmplitude(std::size_t index) const
	{
		return index < modes_.size() ? ModeDisplacementAmplitude(index) *
			std::abs(modeOutputWeight_[index]) : 0.0;
	}
	bool ModeRendered(std::size_t index) const
	{
		return index < modeActive_.size() && modeActive_[index];
	}
	bool ModeParticipatesInContact(std::size_t index) const
	{
		return index < contactModeActive_.size() && contactModeActive_[index];
	}
	double ContactModeProjection(std::size_t index) const
	{
		return index < contactModeShape_.size() ? contactModeShape_[index] : 0.0;
	}
	double StrikePosition() const { return currentStrikePosition_; }
	double ContactWidthMetres() const { return contactWidthMetres_; }
	static std::array<double, 2> MagneticPickupGradient(double vertical,
		double horizontal, double gap)
	{
		return MagneticFluxGradient(vertical, horizontal, gap);
	}
	static std::array<double, 2> DirectMagneticPickupGradient(
		double vertical, double horizontal, double gap)
	{
		constexpr double HorizontalFieldScale = 0.62;
		const double radial = std::sqrt(gap * gap +
			HorizontalFieldScale * horizontal * horizontal + 0.020);
		const auto field = EvaluateCalibratedPoleGradient(vertical, radial);
		return {field.vertical, field.radial * HorizontalFieldScale * horizontal /
			std::max(1.0e-9, radial)};
	}

