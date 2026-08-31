	static std::array<double, 2> ModulationPhaseRotation(double angle)
	{
		// Normal musical FM moves only a small fraction of a cycle per 4x frame.
		// This seventh/sixth-order local rotation is effectively exact in that
		// range and avoids two transcendental calls per mode. Extreme modulation
		// retains the full library implementation rather than range-reducing a
		// polynomial outside its error bound.
		if (std::abs(angle) < 0.35)
		{
			const double squared = angle * angle;
			const double sine = angle * (1.0 + squared * (-1.0 / 6.0 +
				squared * (1.0 / 120.0 - squared / 5040.0)));
			const double cosine = 1.0 + squared * (-1.0 / 2.0 + squared *
				(1.0 / 24.0 - squared / 720.0));
			return {cosine, sine};
		}
		return {std::cos(angle), std::sin(angle)};
	}

	static constexpr double TwoPi = 6.283185307179586476925286766559;
	static constexpr double FieldFalloff = 1.30;
	// Modal coordinates and velocities are now metres and metres/second. The
	// pickup field below is expressed in millimetres. PickupExcursion is the
	// separately documented fit needed because no measured flux map is available.
	static constexpr double TineDisplacementScale = 1000.0 *
		ElectricPianoMechanicalTrim::PickupExcursion;
	// The published SI stiffness is roughly five orders above the former
	// model-space fit. Sixty-four symplectic substeps keep collision spectra
	// sample-rate invariant; this path runs only for the few milliseconds of
	// actual contact and does not affect sustained-voice CPU.
	static constexpr int ContactOversamplingFactor = 64;
	static constexpr double MinimumHammerContactVelocity = 0.012;
	static constexpr double ReferencePickupGap = 1.6;
	static constexpr double ReferencePickupVertical = 0.34 + 0.22 * 0.55 * 0.55 +
		0.020 * ((60.0 - 28.0) / 72.0);
	static constexpr double ReferencePickupHorizontal = 0.10 + 0.035 *
		((60.0 - 28.0) / 72.0);
	// The first three coordinates form the audible coupled tine/tone-bar body.
	// Broadly FM-modulating the intentionally inharmonic, short-lived attack
	// modes produces a sideband cloud rather than the coherent overtones players
	// expect from FM. Keep those impact coordinates physical and let the nonlinear
	// magnetic pickup derive the upper modulation spectrum from the body motion.
	static constexpr std::array<double, ElectricPianoModeCount>
		ModulationModeWeights{
			1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	static constexpr std::array<double, 8> CantileverWaveNumbers{
		1.875104, 4.694091, 7.854757, 10.995541, 14.137168, 17.278760,
		20.420352, 23.561945};
	static constexpr double IntrinsicForkDecaySeconds = 6.8;
	static constexpr double TransverseModeDecaySeconds = 1.9;

	std::array<Mode, ElectricPianoModeCount> modes_{};
	std::array<double, ElectricPianoModeCount> modeRealCoefficient_{};
	std::array<double, ElectricPianoModeCount> modeImaginaryCoefficient_{};
	std::array<double, ElectricPianoModeCount> modeRadius_{};
	std::array<double, ElectricPianoModeCount> modePickupRealCoefficient_{};
	std::array<double, ElectricPianoModeCount> modePickupImaginaryCoefficient_{};
	std::array<double, ElectricPianoModeCount> modePickupRadius_{};
	std::array<double, ElectricPianoModeCount> modeSubstepRealCoefficient_{};
	std::array<double, ElectricPianoModeCount> modeSubstepImaginaryCoefficient_{};
	std::array<double, ElectricPianoModeCount> modeSubstepRadius_{};
	std::array<double, ElectricPianoModeCount> modeAngularFrequency_{};
	std::array<double, ElectricPianoModeCount> modeOutputWeight_{};
	std::array<double, ElectricPianoModeCount> modeHorizontalWeight_{};
	std::array<double, ElectricPianoModeCount> modeRatio_{};
	std::array<double, ElectricPianoModeCount> modeBandlimitGain_{};
	std::array<bool, ElectricPianoModeCount> modeActive_{};
	std::array<bool, ElectricPianoModeCount> contactModeActive_{};
	std::array<double, ElectricPianoModeCount> contactModeShape_{};
	std::array<double, ElectricPianoModeCount> modeInverseMass_{};
	std::array<double, 2> coupledToneBarRatio_{};
	std::array<double, 2> coupledSupportLossFactor_{};
	double sampleRate_ = 48000.0;
	double latchedVelocity_{};
	double keyPosition_ = 0.5;
	double currentStrikePosition_ = 0.22;
	double contactWidthMetres_ = 0.009;
	double tineModalMass_ = 1.0e-4;
	double notePitch_{};
	double keyPickupSensitivity_ = 1.0;
	double pickupExcursionScale_ = 1.0;
	double pickupOutput_{};
	double pickupGap_ = 1.0;
	double pickupVerticalOffset_ = ReferencePickupVertical;
	double pickupHorizontalOffset_ = ReferencePickupHorizontal;
	double inverseReferenceGradient_ = 1.0 / ReferenceGradientMagnitude();
	double pickupAlignmentGain_ = 1.0;
	double pickupVoltageScale_{};
	double mechanicsLevel_{};
	double hammerNoiseColour_ = 0.5;
	double hammerNoiseDecay_ = 0.9979;
	double keyReleaseNoiseDecay_ = 0.9986;
	double damperNoiseDecay_ = 0.9992;
	double hammerNoise_{};
	double keyReleaseNoise_{};
	double damperNoise_{};
	double mechanicsLowPass_{};
	double mechanicsHighPass_{};
	double mechanicsLowPassCoefficient_ = 0.054;
	double mechanicsHighPassCoefficient_ = 0.334;
	std::array<MechanicalResonator, 6> mechanicalResonators_{};
	double hammerPosition_{};
	double hammerVelocity_{};
	double hammerIncomingVelocity_ = 0.004;
	double hammerInverseMass_ =
		1.0 / ElectricPianoPublishedMechanicalData::HammerMassKg;
	double contactStiffness_{};
	double contactExponent_ =
		ElectricPianoPublishedMechanicalData::ContactExponent;
	double contactLoss_{};
	double contactAge_{};
	double controlSmoothingCoefficient_ = 0.0035;
	double decaySmoothingCoefficient_ = 0.00012;
	int modeCoefficientUpdatePeriod_ = 48;
	int modeCoefficientUpdateCountdown_{};
	ElectricPianoControls smoothedControls_{};
	double energy_{};
	std::array<std::array<Mode, ElectricPianoModeCount>, PickupOversamplingFactor>
		pickupModeFrames_{};
	std::array<double, ElectricPianoModeCount> modulationPhasorReal_{};
	std::array<double, ElectricPianoModeCount> modulationPhasorImaginary_{};
	double previousOversampledPhaseModulation_{};
	double previousHostExponentialPitch_{};
	double previousHostLinearFrequencyRatio_{};
	double previousHostPhaseModulation_{};
	double modulationCrossfade_{};
	int modulationPhasorNormalizationCountdown_ = 4096;
	double cachedFundamental_{};
	double currentFundamental_{};
	double cachedPitch_ = -1000.0;
	double cachedDecay_{};
	double cachedBell_{};
	double cachedCoupling_ = -1.0;
	double cachedRelease_{};
	double cachedBodyWeight_ = -1.0;
	double cachedBellWeight_ = -1.0;
	double cachedCouplingWeight_ = -1.0;
	double cachedProximity_ = -1.0;
	double cachedTone_ = -1.0;
	double cachedMechanics_ = -1.0;
	std::uint32_t noiseState_ = 0x6d2b79f5u;
	bool cachedDamped_{};
	bool coefficientsDirty_ = true;
	bool timbreDirty_ = true;
	bool controlsInitialized_{};
	bool coupledForkInitialized_{};
	bool contactActive_{};
	bool contactEngaged_{};
	bool lastGate_{};
	bool lastSustain_{};
	bool keyHeld_{};
	bool modulationPathActive_{};
	bool modulationControlsPrimed_{};
	bool exponentialFmPathActive_{};
	bool linearFmPathActive_{};
	bool phaseModulationPathActive_{};
	bool exponentialFmInterpolatorPrimed_{};
	bool linearFmInterpolatorPrimed_{};
	bool phaseModulationInterpolatorPrimed_{};
	std::unique_ptr<X4Resampler_Order7> verticalPositionInterpolator_;
	std::unique_ptr<X4Resampler_Order7> verticalVelocityInterpolator_;
	std::unique_ptr<X4Resampler_Order7> horizontalPositionInterpolator_;
	std::unique_ptr<X4Resampler_Order7> horizontalVelocityInterpolator_;
	std::unique_ptr<X4Resampler_Order7> exponentialFmInterpolator_;
	std::unique_ptr<X4Resampler_Order7> linearFmInterpolator_;
	std::unique_ptr<X4Resampler_Order7> phaseModulationInterpolator_;
	std::unique_ptr<X4Resampler_Order7> pickupDecimator_;
