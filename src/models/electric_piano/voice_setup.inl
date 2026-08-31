public:
	ElectricPianoVoice()
		: verticalPositionInterpolator_(CreateX4Resampler_Cheby7()),
		  verticalVelocityInterpolator_(CreateX4Resampler_Cheby7()),
		  horizontalPositionInterpolator_(CreateX4Resampler_Cheby7()),
		  horizontalVelocityInterpolator_(CreateX4Resampler_Cheby7()),
		  exponentialFmInterpolator_(CreateX4Resampler_Cheby7()),
		  linearFmInterpolator_(CreateX4Resampler_Cheby7()),
		  phaseModulationInterpolator_(CreateX4Resampler_Cheby7()),
		  pickupDecimator_(CreateX4Resampler_Cheby7())
	{
	}

	void SetNoiseSeed(std::uint32_t seed)
	{
		noiseState_ = seed != 0 ? seed : 0x6d2b79f5u;
		modeCoefficientUpdateCountdown_ = static_cast<int>(noiseState_ %
			static_cast<std::uint32_t>(modeCoefficientUpdatePeriod_));
	}

	void SetSampleRate(double sampleRate)
	{
		if (!std::isfinite(sampleRate))
			throw std::invalid_argument("electric-piano sample rate must be finite");
		sampleRate_ = std::clamp(sampleRate, 8000.0, 768000.0);
		pickupOutput_ = 0.0;
		hammerNoiseDecay_ = std::exp(-1.0 / (0.0045 * sampleRate_));
		keyReleaseNoiseDecay_ = std::exp(-1.0 / (0.012 * sampleRate_));
		damperNoiseDecay_ = std::exp(-1.0 / (0.018 * sampleRate_));
		mechanicsLowPassCoefficient_ = 1.0 - std::exp(
			-TwoPi * 420.0 / sampleRate_);
		mechanicsHighPassCoefficient_ = 1.0 - std::exp(
			-TwoPi * 3100.0 / sampleRate_);
		controlSmoothingCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.006 * sampleRate_));
		// A decay change alters the radii of every resonant mode. Treating it
		// like a fast tone control creates audible amplitude modulation on held
		// notes, so it gets its own deliberately slow, still-live transition.
		decaySmoothingCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.180 * sampleRate_));
		modeCoefficientUpdatePeriod_ = std::max(1,
			static_cast<int>(std::round(sampleRate_ / 1000.0)));
		modeCoefficientUpdateCountdown_ = static_cast<int>(noiseState_ %
			static_cast<std::uint32_t>(modeCoefficientUpdatePeriod_));
		for (auto& resonator : mechanicalResonators_)
			resonator.state = {};
		coefficientsDirty_ = true;
		cachedTone_ = -1.0;
		exponentialFmInterpolator_->Reset();
		linearFmInterpolator_->Reset();
		phaseModulationInterpolator_->Reset();
		verticalPositionInterpolator_->Reset();
		verticalVelocityInterpolator_->Reset();
		horizontalPositionInterpolator_->Reset();
		horizontalVelocityInterpolator_->Reset();
		exponentialFmInterpolatorPrimed_ = false;
		linearFmInterpolatorPrimed_ = false;
		phaseModulationInterpolatorPrimed_ = false;
		pickupDecimator_->Reset();
	}

	void Reset()
	{
		for (auto& mode : modes_)
			mode = {};
		lastGate_ = false;
		lastSustain_ = false;
		keyHeld_ = false;
		latchedVelocity_ = 0.0;
		keyPosition_ = 0.5;
		pickupExcursionScale_ = 1.0;
		pickupOutput_ = 0.0;
		hammerNoise_ = 0.0;
		keyReleaseNoise_ = 0.0;
		damperNoise_ = 0.0;
		mechanicsLowPass_ = 0.0;
		mechanicsHighPass_ = 0.0;
		for (auto& resonator : mechanicalResonators_)
			resonator = {};
		hammerPosition_ = 0.0;
		hammerVelocity_ = 0.0;
		hammerIncomingVelocity_ = MinimumHammerContactVelocity;
		hammerInverseMass_ = 1.0 /
			ElectricPianoPublishedMechanicalData::HammerMassKg;
		contactStiffness_ = 0.0;
		contactExponent_ =
			ElectricPianoPublishedMechanicalData::ContactExponent;
		contactLoss_ = 0.0;
		contactAge_ = 0.0;
		contactActive_ = false;
		contactEngaged_ = false;
		contactModeShape_.fill(0.0);
		modeActive_.fill(false);
		contactModeActive_.fill(false);
		modeInverseMass_.fill(1.0);
		coupledToneBarRatio_.fill(0.0);
		tineModalMass_ = 1.0e-4;
		coupledForkInitialized_ = false;
		energy_ = 0.0;
		modulationPhasorReal_.fill(1.0);
		modulationPhasorImaginary_.fill(0.0);
		modulationPhasorNormalizationCountdown_ = 4096;
		previousOversampledPhaseModulation_ = 0.0;
		previousHostExponentialPitch_ = 0.0;
		previousHostLinearFrequencyRatio_ = 0.0;
		previousHostPhaseModulation_ = 0.0;
		modulationPathActive_ = false;
		modulationControlsPrimed_ = false;
		exponentialFmPathActive_ = false;
		linearFmPathActive_ = false;
		phaseModulationPathActive_ = false;
		exponentialFmInterpolatorPrimed_ = false;
		linearFmInterpolatorPrimed_ = false;
		phaseModulationInterpolatorPrimed_ = false;
		modulationCrossfade_ = 0.0;
		for (auto& frame : pickupModeFrames_)
			frame.fill({});
		controlsInitialized_ = false;
		coefficientsDirty_ = true;
		timbreDirty_ = true;
		exponentialFmInterpolator_->Reset();
		linearFmInterpolator_->Reset();
		phaseModulationInterpolator_->Reset();
		verticalPositionInterpolator_->Reset();
		verticalVelocityInterpolator_->Reset();
		horizontalPositionInterpolator_->Reset();
		horizontalVelocityInterpolator_->Reset();
		pickupDecimator_->Reset();
	}

