private:

	void RenderPickup(
		const Eigen::Array<double, PickupOversamplingFactor, 1>& exponentialPitchFrames,
		const Eigen::Array<double, PickupOversamplingFactor, 1>& linearRatioFrames,
		const Eigen::Array<double, PickupOversamplingFactor, 1>& phaseModulationFrames,
		const double mechanicalSignal)
	{
		// Preserve the established, cheaper aggregate interpolation path whenever
		// FM/PM has never been active on this note. If modulation arrives on a held
		// note, keep it alive briefly and crossfade to the per-mode 4x readout; this
		// avoids both a timing discontinuity and permanent baseline CPU overhead.
		const bool legacyPickupNeeded = !modulationPathActive_ ||
			modulationCrossfade_ < 1.0;
		Eigen::Array<double, PickupOversamplingFactor, 1> legacyVerticalPositions;
		Eigen::Array<double, PickupOversamplingFactor, 1> legacyVerticalVelocities;
		Eigen::Array<double, PickupOversamplingFactor, 1> legacyHorizontalPositions;
		Eigen::Array<double, PickupOversamplingFactor, 1> legacyHorizontalVelocities;
		if (legacyPickupNeeded)
		{
			double verticalPosition = 0.0;
			double verticalVelocity = 0.0;
			double horizontalPosition = 0.0;
			double horizontalVelocity = 0.0;
			for (std::size_t mode = 0; mode < modes_.size(); ++mode)
			{
				if (!modeActive_[mode])
					continue;
				const double verticalWeight = modeOutputWeight_[mode] *
					modeBandlimitGain_[mode];
				const double horizontalWeight = modeHorizontalWeight_[mode] *
					modeBandlimitGain_[mode];
				const double velocityScale = TineDisplacementScale *
					pickupExcursionScale_ *
					modeAngularFrequency_[mode];
				verticalPosition += verticalWeight * modes_[mode].real;
				verticalVelocity -= verticalWeight * velocityScale *
					modes_[mode].imaginary;
				horizontalPosition += horizontalWeight * modes_[mode].real;
				horizontalVelocity -= horizontalWeight * velocityScale *
					modes_[mode].imaginary;
			}
			legacyVerticalPositions = verticalPositionInterpolator_->Upsample(
				TineDisplacementScale * pickupExcursionScale_ * verticalPosition);
			legacyVerticalVelocities = verticalVelocityInterpolator_->Upsample(
				verticalVelocity);
			legacyHorizontalPositions = horizontalPositionInterpolator_->Upsample(
				TineDisplacementScale * pickupExcursionScale_ * horizontalPosition);
			legacyHorizontalVelocities = horizontalVelocityInterpolator_->Upsample(
				horizontalVelocity);
		}
		Eigen::Array<double, PickupOversamplingFactor, 1> pickupValues;
		const double pickupRate = sampleRate_ *
			static_cast<double>(PickupOversamplingFactor);
		for (int frame = 0; frame < PickupOversamplingFactor; ++frame)
		{
			if (!modulationPathActive_)
			{
				const auto legacyGradient = MagneticFluxGradient(
					pickupVerticalOffset_ + legacyVerticalPositions(frame),
					pickupHorizontalOffset_ + legacyHorizontalPositions(frame),
					pickupGap_);
				const double legacyEmf = legacyGradient[0] *
					legacyVerticalVelocities(frame) + legacyGradient[1] *
					legacyHorizontalVelocities(frame);
				pickupValues(frame) = std::clamp(
					pickupVoltageScale_ * legacyEmf, -12.0, 12.0) +
					mechanicalSignal;
				continue;
			}
			double verticalPosition = 0.0;
			double verticalVelocity = 0.0;
			double horizontalPosition = 0.0;
			double horizontalVelocity = 0.0;
			const double frameExponentialPitch = std::clamp(
				exponentialPitchFrames(frame), -4.0, 4.0);
			const double frameLinearRatio = std::clamp(linearRatioFrames(frame),
				-4.0, 4.0);
			const double framePhaseModulation = std::clamp(
				phaseModulationFrames(frame), -2.0 * TwoPi, 2.0 * TwoPi);
			const double bodyExponentialRatio =
				frameExponentialPitch == 0.0 ? 1.0 :
				Exp2Taylor5(static_cast<float>(frameExponentialPitch));
			const double phaseModulationFrequency = modulationPathActive_ ?
				(framePhaseModulation - previousOversampledPhaseModulation_) *
					pickupRate / TwoPi : 0.0;
			for (std::size_t mode = 0;
				mode < ElectricPianoToneBarSubMode; ++mode)
			{
				if (!modeActive_[mode])
					continue;
				const double physicalFrequency = currentFundamental_ *
					modeRatio_[mode];
				double renderedFrequency = physicalFrequency;
				double renderedReal = pickupModeFrames_[frame][mode].real;
				double renderedImaginary =
					pickupModeFrames_[frame][mode].imaginary;
				double renderedBandlimit = modeBandlimitGain_[mode];
				if (modulationPathActive_ && ModulationModeWeights[mode] != 0.0)
				{
					const double modulationWeight = ModulationModeWeights[mode];
					const double weightedExponentialPitch = std::clamp(
						modulationWeight * frameExponentialPitch, -4.0, 4.0);
					const double exponentialRatio = modulationWeight == 1.0 ?
						bodyExponentialRatio :
						(weightedExponentialPitch == 0.0 ? 1.0 :
							Exp2Taylor5(static_cast<float>(weightedExponentialPitch)));
					const double frequencyOffset = physicalFrequency *
						(exponentialRatio - 1.0 +
							modulationWeight * frameLinearRatio);
					const double phaseIncrement = TwoPi * frequencyOffset / pickupRate +
						modulationWeight * (framePhaseModulation -
							previousOversampledPhaseModulation_);
					const auto rotation = ModulationPhaseRotation(phaseIncrement);
					const double oldPhasorReal = modulationPhasorReal_[mode];
					const double oldPhasorImaginary =
						modulationPhasorImaginary_[mode];
					modulationPhasorReal_[mode] = rotation[0] * oldPhasorReal -
						rotation[1] * oldPhasorImaginary;
					modulationPhasorImaginary_[mode] = rotation[1] * oldPhasorReal +
						rotation[0] * oldPhasorImaginary;
					const double phaseCosine = modulationPhasorReal_[mode];
					const double phaseSine = modulationPhasorImaginary_[mode];
					const double physicalReal = renderedReal;
					const double physicalImaginary = renderedImaginary;
					renderedReal = phaseCosine * physicalReal -
						phaseSine * physicalImaginary;
					renderedImaginary = phaseSine * physicalReal +
						phaseCosine * physicalImaginary;
					renderedFrequency = physicalFrequency + frequencyOffset +
						modulationWeight * phaseModulationFrequency;
					renderedBandlimit = std::min(renderedBandlimit,
						ElectricPianoModeBandlimitGain(
							std::abs(renderedFrequency), sampleRate_));
				}
				const double verticalWeight = modeOutputWeight_[mode] *
					renderedBandlimit;
				const double horizontalWeight = modeHorizontalWeight_[mode] *
					renderedBandlimit;
				const double velocityScale = TineDisplacementScale *
					pickupExcursionScale_ * TwoPi * renderedFrequency;
				verticalPosition += verticalWeight * renderedReal;
				verticalVelocity -= verticalWeight * velocityScale *
					renderedImaginary;
				horizontalPosition += horizontalWeight * renderedReal;
				horizontalVelocity -= horizontalWeight * velocityScale *
					renderedImaginary;
			}
			// The separately measured tone-bar submode remains physical under the
			// creative FM/PM inputs. Keep it out of the hot ten-coordinate modulation
			// loop: besides avoiding meaningless phasor work, this preserves the
			// compiler's established vectorization for polyphonic audio-rate FM.
			constexpr std::size_t subMode = ElectricPianoToneBarSubMode;
			if (modeActive_[subMode])
			{
				const double renderedBandlimit = modeBandlimitGain_[subMode];
				const double verticalWeight = modeOutputWeight_[subMode] *
					renderedBandlimit;
				const double horizontalWeight = modeHorizontalWeight_[subMode] *
					renderedBandlimit;
				const double renderedFrequency = currentFundamental_ *
					modeRatio_[subMode];
				const double velocityScale = TineDisplacementScale *
					pickupExcursionScale_ * TwoPi * renderedFrequency;
				const double renderedReal =
					pickupModeFrames_[frame][subMode].real;
				const double renderedImaginary =
					pickupModeFrames_[frame][subMode].imaginary;
				verticalPosition += verticalWeight * renderedReal;
				verticalVelocity -= verticalWeight * velocityScale *
					renderedImaginary;
				horizontalPosition += horizontalWeight * renderedReal;
				horizontalVelocity -= horizontalWeight * velocityScale *
					renderedImaginary;
			}
			if (modulationPathActive_)
				previousOversampledPhaseModulation_ = framePhaseModulation;
			double modulatedPickup = 0.0;
			if (modulationPathActive_)
			{
				const auto gradient = MagneticFluxGradient(
					pickupVerticalOffset_ + TineDisplacementScale *
						pickupExcursionScale_ * verticalPosition,
					pickupHorizontalOffset_ + TineDisplacementScale *
						pickupExcursionScale_ * horizontalPosition,
					pickupGap_);
				const double emf = gradient[0] * verticalVelocity +
					gradient[1] * horizontalVelocity;
				modulatedPickup = pickupVoltageScale_ * emf;
			}
			double selectedPickup = modulatedPickup;
			if (legacyPickupNeeded)
			{
				const auto legacyGradient = MagneticFluxGradient(
					pickupVerticalOffset_ + legacyVerticalPositions(frame),
					pickupHorizontalOffset_ + legacyHorizontalPositions(frame),
					pickupGap_);
				const double legacyEmf = legacyGradient[0] *
					legacyVerticalVelocities(frame) + legacyGradient[1] *
					legacyHorizontalVelocities(frame);
				const double legacyPickup = pickupVoltageScale_ * legacyEmf;
				selectedPickup = modulationPathActive_ ?
					legacyPickup + modulationCrossfade_ *
						(modulatedPickup - legacyPickup) : legacyPickup;
			}
			pickupValues(frame) = std::clamp(selectedPickup, -12.0, 12.0) +
				mechanicalSignal;
		}
		if (modulationPathActive_ && modulationCrossfade_ < 1.0)
			modulationCrossfade_ = std::min(1.0, modulationCrossfade_ +
				1.0 / std::max(1.0, 0.006 * sampleRate_));
		if (modulationPathActive_ &&
			--modulationPhasorNormalizationCountdown_ <= 0)
		{
			for (std::size_t mode = 0; mode < modes_.size(); ++mode)
			{
				const double magnitude = std::hypot(modulationPhasorReal_[mode],
					modulationPhasorImaginary_[mode]);
				if (magnitude > 0.0)
				{
					modulationPhasorReal_[mode] /= magnitude;
					modulationPhasorImaginary_[mode] /= magnitude;
				}
			}
			modulationPhasorNormalizationCountdown_ = 4096;
		}
		// The seventh-order 4x decimator is the pickup reconstruction filter.
		// A former host-rate 16.5 kHz one-pole was not part of the pickup or
		// Peterson circuit and added 0.6--1.9 dB of avoidable upper-band loss.
		// Keep the decimated signal directly; alias rejection is verified at the
		// actual 4x boundary by the register-wide spectral tests.
		pickupOutput_ = pickupDecimator_->Downsample(pickupValues);
	}

	struct Mode
	{
		double real{};
		double imaginary{};
	};

	struct MechanicalResonator
	{
		Mode state{};
		double realCoefficient = 1.0;
		double imaginaryCoefficient{};
	};

