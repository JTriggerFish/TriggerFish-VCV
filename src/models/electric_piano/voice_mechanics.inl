	static double Clamp01(double value)
	{
		return std::clamp(std::isfinite(value) ? value : 0.0, 0.0, 1.0);
	}

	void TriggerMechanicalResonator(MechanicalResonator& resonator,
		double frequency, double amplitude = 1.0)
	{
		const double boundedFrequency = std::clamp(frequency, 20.0,
			0.43 * sampleRate_);
		const double angle = TwoPi * boundedFrequency / sampleRate_;
		resonator.realCoefficient = std::cos(angle);
		resonator.imaginaryCoefficient = std::sin(angle);
		// A cosine-phase start puts the physical impact at the event sample;
		// it cannot wait for an arbitrary random carrier crossing.
		resonator.state.real = amplitude;
		resonator.state.imaginary = 0.0;
	}

	static double AdvanceMechanicalResonator(
		MechanicalResonator& resonator)
	{
		const double output = resonator.state.real;
		const double oldReal = resonator.state.real;
		const double oldImaginary = resonator.state.imaginary;
		resonator.state.real = resonator.realCoefficient * oldReal -
			resonator.imaginaryCoefficient * oldImaginary;
		resonator.state.imaginary = resonator.imaginaryCoefficient * oldReal +
			resonator.realCoefficient * oldImaginary;
		return output;
	}

	void TriggerKeyReleaseMechanics()
	{
		const double baseFrequency = 310.0 + 450.0 * keyPosition_;
		TriggerMechanicalResonator(mechanicalResonators_[2], baseFrequency);
		TriggerMechanicalResonator(mechanicalResonators_[3],
			1.86 * baseFrequency);
	}

	void TriggerDamperMechanics()
	{
		const double baseFrequency = 950.0 + 850.0 * keyPosition_;
		TriggerMechanicalResonator(mechanicalResonators_[4], baseFrequency);
		TriggerMechanicalResonator(mechanicalResonators_[5],
			1.52 * baseFrequency);
	}

	ElectricPianoControls SmoothActiveControls(
		const ElectricPianoControls& controls)
	{
		if (!controlsInitialized_)
		{
			// Only these continuously active controls are consumed below. Sanitize
			// them while priming so a non-finite API caller cannot poison the
			// smoother's persistent state before Clamp01() sees the target.
			smoothedControls_.body = Clamp01(controls.body);
			smoothedControls_.bell = Clamp01(controls.bell);
			smoothedControls_.coupling = Clamp01(controls.coupling);
			smoothedControls_.tone = Clamp01(controls.tone);
			smoothedControls_.proximity = Clamp01(controls.proximity);
			smoothedControls_.decay = Clamp01(controls.decay);
			smoothedControls_.release = Clamp01(controls.release);
			smoothedControls_.mechanics = Clamp01(controls.mechanics);
			controlsInitialized_ = true;
		}
		auto smooth = [&](double& current, double target, double coefficient)
		{
			target = Clamp01(target);
			current += coefficient * (target - current);
		};
		smooth(smoothedControls_.body, controls.body,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.bell, controls.bell,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.coupling, controls.coupling,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.tone, controls.tone,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.proximity, controls.proximity,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.decay, controls.decay,
			decaySmoothingCoefficient_);
		smooth(smoothedControls_.release, controls.release,
			controlSmoothingCoefficient_);
		smooth(smoothedControls_.mechanics, controls.mechanics,
			controlSmoothingCoefficient_);
		return smoothedControls_;
	}

	void Strike(double pitchVolts, double velocity,
		const ElectricPianoControls& controls)
	{
		const auto key = MakeElectricPianoKeyProfile(pitchVolts);
		notePitch_ = std::clamp(pitchVolts, -6.0, 6.0);
		keyPosition_ = key.keyboardPosition;
		keyPickupSensitivity_ = key.pickupSensitivity;
		const double bassExcursionPosition = std::clamp(keyPosition_ /
			((60.0 - 28.0) / 72.0), 0.0, 1.0);
		const double bassExcursionSmooth = bassExcursionPosition *
			bassExcursionPosition * (3.0 - 2.0 * bassExcursionPosition);
		// The published gap range tops out at 1/8 inch, yet the reduced uniform-
		// beam coordinates over-predict the excursion of the longest, tapered
		// tines. Keep this visible as a bass-only calibration awaiting measured
		// tine-tip trajectories; it reaches unity at C4 and never boosts motion.
		pickupExcursionScale_ = 0.55 + 0.45 * bassExcursionSmooth;
		coefficientsDirty_ = true;
		timbreDirty_ = true;
		const double velocityCurve = Clamp01(controls.velocityCurve);
		const double dynamics = Clamp01(controls.dynamics);
		const double hammer = Clamp01(controls.hammer);
		const double strikePosition = Clamp01(controls.strikePosition);
		const double gamma = std::exp2(1.0 - 2.0 * velocityCurve);
		const double curvedVelocity = velocity > 0.0 ?
			std::pow(velocity, gamma) : 0.0;
		latchedVelocity_ = curvedVelocity;
		const double compressedVelocity = curvedVelocity /
			(0.28 + 0.72 * curvedVelocity);
		const double dynamicAmplitude = (1.0 - dynamics) *
			compressedVelocity + dynamics * curvedVelocity;

		// The hammer is an independent finite mass. MIDI velocity sets its
		// incoming speed; no post-contact amplitude or brightness envelope is
		// prescribed. The key-dependent term represents the graduated neoprene
		// tips used across the real keyboard, while HAMMER moves around that
		// baseline by changing contact material properties only.
		const double factoryDurometer =
			ElectricPianoFactoryHammerTipDurometer(keyPosition_);
		const double factoryHardness = (factoryDurometer - 30.0) / 70.0;
		const double contactHardness = std::clamp(factoryHardness +
			0.70 * (hammer - 0.52), 0.0, 1.0);
		hammerInverseMass_ = 1.0 /
			ElectricPianoPublishedMechanicalData::HammerMassKg;
		hammerPosition_ = 0.0;
		const double requestedHammerVelocity =
			ElectricPianoPublishedMechanicalData::
				MaximumHammerVelocityMetresPerSecond *
			ElectricPianoMechanicalTrim::HammerVelocity *
			std::pow(dynamicAmplitude,
				ElectricPianoMechanicalTrim::HammerVelocityCurveExponent);
		// Strikes below this physical speed are many orders below the audible
		// output range but can otherwise spend seconds in an almost static
		// oversampled collision when driven by tiny positive CV residue. Classify
		// them as zero-energy before contact begins; never truncate an active
		// collision based on elapsed time.
		const bool audibleHammerStrike =
			requestedHammerVelocity >= MinimumHammerContactVelocity;
		hammerVelocity_ = audibleHammerStrike ? requestedHammerVelocity : 0.0;
		hammerIncomingVelocity_ = std::max(MinimumHammerContactVelocity,
			requestedHammerVelocity);
		// The service manual documents five graduated tip grades (30, 50, 70,
		// 90 durometer, then wrapped extra-hard) but not their force curves.  The
		// published k is therefore the middle-C/default anchor.  The modest zoned
		// slope and the symmetric panel span are explicit listening trims around
		// that anchor, not substitutes for the sourced SI stiffness.
		const double graduatedTipOctaves = 1.5 * std::log2(
			factoryDurometer / 50.0);
		const double panelHardnessOctaves = 5.0 * (hammer - 0.52);
		contactStiffness_ =
			ElectricPianoPublishedMechanicalData::
				ContactStiffnessNewtonPerMetrePower *
			ElectricPianoMechanicalTrim::ContactStiffness *
			std::exp2(graduatedTipOctaves + panelHardnessOctaves);
		contactExponent_ =
			ElectricPianoPublishedMechanicalData::ContactExponent;
		// Sonderbo reports lambda separately from k. Scale lambda with the same
		// tip-grade factor so lambda/k (seconds/metre) remains the measured 0.6
		// while Hammer changes compliance rather than inventing restitution.
		contactLoss_ =
			ElectricPianoPublishedMechanicalData::ContactDampingWeight /
			ElectricPianoPublishedMechanicalData::
				ContactStiffnessNewtonPerMetrePower;
		contactAge_ = 0.0;
		contactActive_ = audibleHammerStrike;
		contactEngaged_ = false;
		currentStrikePosition_ = StrikePositionFromControl(keyPosition_,
			strikePosition);
		// Falaize's distributed-force reference uses a 15 mm hammer width. The
		// effective loaded strip is smaller than the complete tip, so retain the
		// calibrated 6--12 mm hard-to-soft interval. Deriving this longitudinal
		// width from Hunt--Crossley indentation assumed a spherical 4 mm tip that
		// is neither present in the source nor representative of the block-shaped
		// Rhodes neoprene tip; it narrowed the footprint to sub-5.5 mm values and
		// changed the previously validated bass modal balance.
		contactWidthMetres_ = 0.006 + 0.006 * (1.0 - contactHardness);
		const double fundamentalProjection = FiniteContactModeProjection(0,
			keyPosition_, currentStrikePosition_, contactWidthMetres_);
		contactModeShape_[0] = fundamentalProjection;
		contactModeShape_[1] = fundamentalProjection;
		contactModeShape_[2] = -0.055 * fundamentalProjection;
		tineModalMass_ = key.tineModalMassKg;
		modeInverseMass_[2] = 1.0 / (key.tineModalMassKg * 1.65);
		for (std::size_t index = ElectricPianoAttackModeBegin;
			index < ElectricPianoAttackModeEnd; ++index)
		{
			// Generalized force projection is a property of contact geometry and the
			// resonator mode, never of MIDI velocity.  Integrate each mode over a
			// finite neoprene patch; contact duration then supplies velocity-dependent
			// bandwidth through the Hunt-Crossley force itself.
			contactModeShape_[index] = FiniteContactModeProjection(index,
				keyPosition_, currentStrikePosition_, contactWidthMetres_) *
				AttackModeEnergyCalibration(index, keyPosition_);
			modeInverseMass_[index] = 1.0 /
				(key.tineModalMassKg *
					AttackModeModalMassMultiplier(index, keyPosition_));
		}
		// The tone bar has a separately observed sub-fundamental mode in addition
		// to its strong played fundamental. In normal-coordinate form the tine
		// participation projects both the hammer force and pickup observation;
		// applying it reciprocally avoids an arbitrary output-only oscillator.
		const double subModeParticipation =
			ToneBarSubModeTineParticipation(keyPosition_) *
			ElectricPianoToneBarSubModeCouplingScale(controls.coupling,
				keyPosition_);
		contactModeShape_[ElectricPianoToneBarSubMode] =
			subModeParticipation * fundamentalProjection;
		modeInverseMass_[ElectricPianoToneBarSubMode] =
			1.0 / key.tineModalMassKg;

		// Mechanics remains a separate sound-design layer. Tip hardness has only
		// a restrained influence on its colour so it cannot duplicate the modal
		// attack as a second, disconnected click.
		hammerNoiseColour_ = 0.24 + 0.22 * contactHardness;
		const double hammerMechanicalFrequency = 680.0 +
			1250.0 * keyPosition_ + 620.0 * contactHardness;
		TriggerMechanicalResonator(mechanicalResonators_[0],
			hammerMechanicalFrequency);
		TriggerMechanicalResonator(mechanicalResonators_[1],
			2.17 * hammerMechanicalFrequency);
		hammerNoise_ = dynamicAmplitude * (0.085 + 0.075 * curvedVelocity) *
			(0.82 + 0.18 * keyPosition_);
		keyReleaseNoise_ = 0.0;
		damperNoise_ = 0.0;
		keyHeld_ = true;
		energy_ = 1.0;
	}

