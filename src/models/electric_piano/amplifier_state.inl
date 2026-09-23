private:
	struct PowerChannel
	{
		PetersonPowerAmplifier circuit;
		double supplyCurrent{};
	};

	double ProcessPowerModule(double input, PowerChannel& channel,
		double supplyRail)
	{
		const auto result = channel.circuit.Step(input, supplyRail);
		channel.supplyCurrent = result.positiveRailCurrent +
			result.negativeRailCurrent;
		return result.voltage;
	}

	// With the module's 5x pickup-sum feed, 0.13 V/model-unit makes a calibrated
	// hard single note present about 0.42 Vpk to the preamp at maximum Drive.
	// After nominal preamp-gain normalization, Figure 11-9 retains the existing
	// voltage calibration while real preamp compression reduces its overload
	// demand. The default remains in both circuits' feedback-linear range.
	static constexpr double InputVoltsPerModelUnit = 0.13;
	static constexpr double PowerClosedLoopGain = 57.0; // 1 + 5.6k / 100
	// Small-signal gain of Q1/Q2 plus the centred Figure 11-8 nonlinear tone
	// feedback stage at 250 Hz. This is a unit calibration only: there is no
	// normalization between the two physical stages, so interstage overload and
	// gain compression remain in the circuit path.
	static constexpr double NominalFullPreampGain = 6.3230;
	static constexpr double NominalSupplyRail = 35.0;
	static constexpr double MinimumSupplyRail = 27.0;
	static constexpr double SupplyResistance = 0.72;
	static constexpr double NominalGain = 1.12;

	double sampleRate_ = 48000.0;
	PetersonPreAmplifier preamp_{};
	PetersonTonePreAmplifier tonePreamp_{};
	std::array<PowerChannel, 2> powerChannels_{};
	double supplyRail_ = NominalSupplyRail;
	double vibratoPhase_{};
	double vibratoLamp_{};
	double cachedDrive_ = -1.0;
	double cachedBass_ = -1.0;
	double cachedTreble_ = -1.0;
	double cachedVibratoSpeed_ = -1.0;
	double smoothedDrive_{};
	double driveSmoothingCoefficient_ = 0.0035;
	double controlSmoothingCoefficient_ = 0.002;
	double vibratoLampCoefficient_ = 0.002;
	double supplyDischargeCoefficient_ = 0.0069;
	double supplyRechargeCoefficient_ = 0.00038;
	double driveGain_ = 1.0;
	double makeupGain_ = NominalGain;
	double bassPotPosition_ = 0.5;
	double treblePotPosition_ = 0.5;
	double vibratoFrequency_ = 3.0;
	double smoothedOutputVolume_ = 0.50;
	double smoothedBass_ = 0.50;
	double smoothedTreble_ = 0.50;
	double smoothedVibratoSpeed_ = 0.32;
	double smoothedVibratoIntensity_{};
	bool controlsInitialized_{};
	bool rightPowerDormant_{true};
	std::unique_ptr<X2Resampler_Order7> inputInterpolator_;
	std::unique_ptr<X2Resampler_Order7> leftOutputDecimator_;
	std::unique_ptr<X2Resampler_Order7> rightOutputDecimator_;
