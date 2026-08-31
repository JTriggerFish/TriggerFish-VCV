public:
	ElectricPianoAmplifier()
		: inputInterpolator_(CreateX2Resampler_Chebychev7()),
		  leftOutputDecimator_(CreateX2Resampler_Chebychev7()),
		  rightOutputDecimator_(CreateX2Resampler_Chebychev7())
	{
	}

	void SetSampleRate(double sampleRate)
	{
		if (!std::isfinite(sampleRate))
			throw std::invalid_argument("electric-piano amplifier sample rate must be finite");
		sampleRate_ = std::clamp(sampleRate, 8000.0, 768000.0);
		const double oversampledRate = sampleRate_ *
			X2Resampler_Order7::ResamplingFactor;
		// Figure 11-10 specifies 3000 uF reservoirs. Effective winding,
		// rectifier and 0.5 Ohm emitter-path resistance determine the fast
		// droop; recharge is slower than discharge between mains peaks.
		supplyDischargeCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.0030 * sampleRate_));
		supplyRechargeCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.055 * sampleRate_));
		for (auto& channel : powerChannels_)
			channel.circuit.SetSampleRate(oversampledRate);
		preamp_.SetSampleRate(oversampledRate);
		tonePreamp_.SetSampleRate(oversampledRate);
		driveSmoothingCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.008 * sampleRate_));
		controlSmoothingCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.012 * sampleRate_));
		vibratoLampCoefficient_ = 1.0 - std::exp(
			-1.0 / (0.010 * sampleRate_));
		Reset();
	}

	void Reset()
	{
		for (auto& channel : powerChannels_)
		{
			channel.circuit.Reset(NominalSupplyRail);
			channel.supplyCurrent = 0.0;
		}
		preamp_.Reset();
		tonePreamp_.Reset();
		supplyRail_ = NominalSupplyRail;
		vibratoPhase_ = 0.0;
		vibratoLamp_ = 0.0;
		cachedDrive_ = -1.0;
		cachedBass_ = -1.0;
		cachedTreble_ = -1.0;
		cachedVibratoSpeed_ = -1.0;
		controlsInitialized_ = false;
		rightPowerDormant_ = true;
		inputInterpolator_->Reset();
		leftOutputDecimator_->Reset();
		rightOutputDecimator_->Reset();
	}

	// Small-signal impedance of the electrical speaker equivalent used by the
	// power-stage solve. This is deliberately public so tests and future fitting
	// tools can validate the load independently of the nonlinear amplifier.
	// It is not an acoustic cabinet response.
	static std::complex<double> ElectricalLoadImpedance(double frequency)
	{
		return PetersonPowerAmplifier::ElectricalLoadImpedance(frequency);
	}

	static double DriveGainForPosition(double drive)
	{
		drive = std::clamp(std::isfinite(drive) ? drive : 0.32, 0.0, 1.0);
		// This is an explicit extension, not a historical Peterson control. The
		// eased-quadratic dB law keeps the default and middle of the travel in the
		// clean headroom region, then approaches compression progressively. Its
		// reciprocal is
		// applied at Figure 11-8's real interstage volume node, before Figure 11-9:
		// Drive can exercise the preamp without turning the power amp into the
		// compensator or making polyphonic peaks hit its rails. No additional
		// clipper is introduced.
		// With the complete tone-feedback circuit present, the former smoothstep
		// law already delivered 20 dB by 70% travel. Simultaneous hard notes then
		// crossed both preamp headroom boundaries together. A 24 dB eased span
		// retains a deliberately overdriven end stop while keeping moderate chord
		// settings out of that abrupt shared overload region. Blending 40% of a
		// quadratic with 60% smoothstep avoids moving all audible change to the last
		// few degrees of travel.
		const double driveShape = drive * drive * (2.2 - 1.2 * drive);
		const double driveDb = 24.0 * driveShape;
		return std::pow(10.0, driveDb / 20.0);
	}

	static double BassPotPositionForControl(double control)
	{
		control = std::clamp(std::isfinite(control) ? control : 0.5, 0.0, 1.0);
		// The historical 100 k linear pot is highly bunched electrically: most
		// of both cut and boost occurs in the last few degrees at either end.
		// Preserve its centre and end stops, but use a reverse-quarter-power UI
		// law so equal panel travel produces a much more even dB trajectory.
		const double displacement = 2.0 * std::abs(control - 0.5);
		const double mappedDisplacement = std::sqrt(std::sqrt(displacement));
		return 0.5 + std::copysign(0.5 * mappedDisplacement,
			control - 0.5);
	}

	static double TreblePotPositionForControl(double control)
	{
		control = std::clamp(std::isfinite(control) ? control : 0.5, 0.0, 1.0);
		const double displacement = 2.0 * std::abs(control - 0.5);
		const double mappedDisplacement = std::pow(displacement, 0.70);
		return 0.5 + std::copysign(0.5 * mappedDisplacement,
			control - 0.5);
	}

	// Explicit circuit/Rack-domain boundary used by calibration tests. The
	// caller supplies the amplifier's model-domain input (the module currently
	// feeds five times the summed pickup signal).
	static double ModelInputToCircuitVolts(double modelInput, double drive)
	{
		return InputVoltsPerModelUnit * DriveGainForPosition(drive) * modelInput;
	}

	static double CircuitPowerOutputToModel(double circuitVolts)
	{
		return NominalGain * circuitVolts /
			(PowerClosedLoopGain * InputVoltsPerModelUnit);
	}

	double SupplyRailVoltage() const { return supplyRail_; }

	std::array<double, 2> Step(double input,
		const ElectricPianoControls& controls)
	{
		if (!std::isfinite(input))
			input = 0.0;
		auto control = [](double value, double fallback)
		{
			return std::clamp(std::isfinite(value) ? value : fallback, 0.0, 1.0);
		};
		const double drive = control(controls.drive, 0.32);
		const double outputVolume = control(controls.outputVolume, 0.50);
		const double bass = control(controls.amplifierBass, 0.50);
		const double treble = control(controls.amplifierTreble, 0.50);
		const double vibratoSpeed = control(controls.vibratoSpeed, 0.32);
		const double vibratoIntensity = control(controls.vibratoIntensity, 0.0);
		if (!controlsInitialized_)
		{
			smoothedDrive_ = drive;
			smoothedOutputVolume_ = outputVolume;
			smoothedBass_ = bass;
			smoothedTreble_ = treble;
			smoothedVibratoSpeed_ = vibratoSpeed;
			smoothedVibratoIntensity_ = vibratoIntensity;
			controlsInitialized_ = true;
		}
		smoothedDrive_ += driveSmoothingCoefficient_ *
			(drive - smoothedDrive_);
		smoothedOutputVolume_ += controlSmoothingCoefficient_ *
			(outputVolume - smoothedOutputVolume_);
		smoothedBass_ += controlSmoothingCoefficient_ *
			(bass - smoothedBass_);
		smoothedTreble_ += controlSmoothingCoefficient_ *
			(treble - smoothedTreble_);
		smoothedVibratoSpeed_ += controlSmoothingCoefficient_ *
			(vibratoSpeed - smoothedVibratoSpeed_);
		smoothedVibratoIntensity_ += controlSmoothingCoefficient_ *
			(vibratoIntensity - smoothedVibratoIntensity_);
		const bool stereoPowerActive = smoothedVibratoIntensity_ > 1.0e-7;
		if (stereoPowerActive && rightPowerDormant_)
		{
			// Both channels have received identical drive while vibrato was off.
			// Materialize the second circuit state only when stereo routing starts.
			powerChannels_[1] = powerChannels_[0];
			rightPowerDormant_ = false;
		}
		else if (!stereoPowerActive)
		{
			rightPowerDormant_ = true;
		}
		if (smoothedDrive_ != cachedDrive_)
		{
			// The historical volume control did not create an independent clipper.
			// This extended Drive control raises the signal presented to the preamp;
			// its small-signal gain is compensated at the real volume node before the
			// power modules. The eased law is deliberately shallow around the clean
			// default and steeper through the deliberate upper-range transition.
			driveGain_ = DriveGainForPosition(smoothedDrive_);
			makeupGain_ = NominalGain;
			cachedDrive_ = smoothedDrive_;
		}
		if (smoothedVibratoSpeed_ != cachedVibratoSpeed_)
		{
			vibratoFrequency_ = 1.5 * std::pow(8.0, smoothedVibratoSpeed_);
			cachedVibratoSpeed_ = smoothedVibratoSpeed_;
		}
		if (smoothedBass_ != cachedBass_)
		{
			bassPotPosition_ = BassPotPositionForControl(smoothedBass_);
			cachedBass_ = smoothedBass_;
		}
		if (smoothedTreble_ != cachedTreble_)
		{
			treblePotPosition_ = TreblePotPositionForControl(smoothedTreble_);
			cachedTreble_ = smoothedTreble_;
		}
		// Interpolate before Figure 11-8. The complete preamp -> tone -> Figure
		// 11-9 path runs at 2x; only controls and the slow supply remain at host
		// rate. The pickup/voice model has its own independent 4x path.
		const auto inputs = inputInterpolator_->Upsample(input);
		Eigen::Array<double, X2Resampler_Order7::ResamplingFactor, 1> leftOutputs;
		Eigen::Array<double, X2Resampler_Order7::ResamplingFactor, 1> rightOutputs;
		const double compensatedPreampDenominator =
			driveGain_ * NominalFullPreampGain;
		const double pan = smoothedVibratoIntensity_ * vibratoLamp_;
		const double leftDrive = stereoPowerActive ?
			std::sqrt(std::max(0.0, 1.0 + pan)) : 1.0;
		const double rightDrive = stereoPowerActive ?
			std::sqrt(std::max(0.0, 1.0 - pan)) : 1.0;
		const double modelScale = makeupGain_ /
			(PowerClosedLoopGain * InputVoltsPerModelUnit);
		for (int index = 0; index < X2Resampler_Order7::ResamplingFactor; ++index)
		{
			// The piano/pickup signal is expressed in model units. Convert it to
			// the physical voltage presented to the +25 V Peterson input transistor.
			// Figure 11-8 places that stage before its tone/feedback and volume
			// sections. Extended Drive is compensated at the pre-power volume node.
			const double physicalInput = InputVoltsPerModelUnit * inputs(index);
			// Drive raises the voltage presented to the Figure 11-8 Q1/Q2 circuit.
			// Its reciprocal is applied at Figure 11-8's volume node below, so the
			// control changes preamp headroom rather than becoming a volume control.
			// Q1/Q2 feed the real active feedback network. Do not normalize between
			// those stages: their physical interstage voltage is what determines when
			// the Darlington feedback pair and Q5 leave their small-signal region.
			const double inputStageOutput = preamp_.Step(
				driveGain_ * physicalInput).voltage;
			const double toneOutput = tonePreamp_.Step(inputStageOutput,
				bassPotPosition_, treblePotPosition_).voltage;
			// Figure 11-8 places its 100 kOhm volume control before the channel
			// output buffers and Figure 11-9. Use that physical node for Drive's
			// reciprocal gain. Putting the compensation after Figure 11-9 preserved
			// the output level but still drove chord peaks into the power rails.
			const double preampOutput = toneOutput /
				compensatedPreampDenominator;

			// In Figure 11-8 the optical network routes the preamp signal to two
			// outputs, each of which feeds its own Figure 11-9 power module. Keep
			// those power stages independent so their crossover and overload history
			// follow the channel currents rather than a post-distortion pan law.
			const double leftVoltage = ProcessPowerModule(preampOutput * leftDrive,
				powerChannels_[0], supplyRail_);
			double rightVoltage = leftVoltage;
			if (stereoPowerActive)
			{
				rightVoltage = ProcessPowerModule(preampOutput * rightDrive,
					powerChannels_[1], supplyRail_);
			}
			leftOutputs(index) = modelScale * leftVoltage;
			rightOutputs(index) = modelScale * rightVoltage;
		}
		const double left = leftOutputDecimator_->Downsample(leftOutputs);
		const double right = rightOutputDecimator_->Downsample(rightOutputs);
		const double totalCurrent = stereoPowerActive ?
			powerChannels_[0].supplyCurrent + powerChannels_[1].supplyCurrent :
			2.0 * powerChannels_[0].supplyCurrent;
		const double loadedRail = std::clamp(NominalSupplyRail -
			SupplyResistance * totalCurrent, MinimumSupplyRail,
			NominalSupplyRail);
		const double supplyCoefficient = loadedRail < supplyRail_ ?
			supplyDischargeCoefficient_ : supplyRechargeCoefficient_;
		supplyRail_ += supplyCoefficient * (loadedRail - supplyRail_);

		vibratoPhase_ += vibratoFrequency_ / sampleRate_;
		if (vibratoPhase_ >= 1.0)
			vibratoPhase_ -= 1.0;
		// Keep the physical lamp/LDR state running while Intensity is down. This
		// costs one slow nonlinear evaluation but makes enabling Vibrato resume the
		// continuously evolving oscillator instead of restarting its lamp envelope.
		const double sine = tfdsp::SinTwoPi(vibratoPhase_);
		const double lampTarget = tfdsp::TanhPade76(1.6 * sine) /
			tfdsp::TanhPade76(1.6);
		vibratoLamp_ += vibratoLampCoefficient_ *
			(lampTarget - vibratoLamp_);
		const double outputGain = 2.0 * smoothedOutputVolume_;
		return {
			outputGain * (std::isfinite(left) ? left : 0.0),
			outputGain * (std::isfinite(right) ? right : 0.0)};
	}

