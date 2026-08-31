	void RefreshModeCoefficients(double fundamental, bool damped,
		const ElectricPianoControls& controls, bool controlTick)
	{
		const double decay = Clamp01(controls.decay);
		const double bell = Clamp01(controls.bell);
		const double coupling = Clamp01(controls.coupling);
		const double release = Clamp01(controls.release);
		const bool stateChange = coefficientsDirty_ ||
			fundamental != cachedFundamental_ || damped != cachedDamped_;
		if (!stateChange && !controlTick)
			return;
		const bool couplingChanged = coefficientsDirty_ ||
			coupling != cachedCoupling_;
		const bool frequencyChanged = coefficientsDirty_ ||
			fundamental != cachedFundamental_ || couplingChanged;
		const bool decayChanged = coefficientsDirty_ ||
			damped != cachedDamped_ || decay != cachedDecay_ ||
			bell != cachedBell_ || couplingChanged || release != cachedRelease_;
		if (!frequencyChanged && !decayChanged)
			return;

		const bool transformCoupledState = coupledForkInitialized_ &&
			frequencyChanged &&
			modeAngularFrequency_[0] > 0.0 && modeAngularFrequency_[1] > 0.0;
		double oldTinePosition = 0.0;
		double oldToneBarPosition = 0.0;
		double oldTineVelocity = 0.0;
		double oldToneBarVelocity = 0.0;
		if (transformCoupledState)
		{
			for (std::size_t mode = 0; mode < 2; ++mode)
			{
				const double modalVelocity = -modeAngularFrequency_[mode] *
					modes_[mode].imaginary;
				oldTinePosition += modes_[mode].real;
				oldToneBarPosition += coupledToneBarRatio_[mode] *
					modes_[mode].real;
				oldTineVelocity += modalVelocity;
				oldToneBarVelocity += coupledToneBarRatio_[mode] * modalVelocity;
			}
		}

		const double keyboardDecayScale = 1.18 - 0.50 * keyPosition_;
		if (couplingChanged)
		{
			// Coupled-fork eigenvectors and support reaction depend on Coupling and
			// key geometry, but not instantaneous pitch. Cache them so future FM only
			// updates oscillator angles and the state-preserving velocity transform.
			const auto coupledFork = MakeElectricPianoCoupledForkProfile(coupling,
				keyPosition_);
			for (std::size_t index = 0; index < 2; ++index)
			{
				modeRatio_[index] = coupledFork.frequencyRatios[index];
				coupledToneBarRatio_[index] =
					coupledFork.toneBarDisplacementRatios[index];
				coupledSupportLossFactor_[index] =
					coupledFork.supportReactionLossFactors[index];
				modeInverseMass_[index] =
					coupledFork.inverseModalMassRatios[index] /
					std::max(1.0e-9, tineModalMass_);
			}
			modeRatio_[2] = 1.0018 + 0.0010 * (keyPosition_ - 0.5);
			for (std::size_t index = ElectricPianoAttackModeBegin;
				index < ElectricPianoAttackModeEnd; ++index)
				modeRatio_[index] = KeyboardModeRatio(index, keyPosition_);
			modeRatio_[ElectricPianoToneBarSubMode] =
				ToneBarSubModeFrequencyRatio(keyPosition_);
		}

		if (decayChanged)
		{
			const double releaseSeconds =
				ElectricPianoDamperReleaseSeconds(release);
			const double decayScale = std::pow(4.0, 2.0 * (decay - 0.5));
			// Published modal lifetimes are the default calibration. Bell may extend
			// or shorten the attack modes deliberately, but its factory 0.52 position
			// must not silently add the former 1.76% hidden correction.
			const double bellDecay = 1.0 + 0.38 * (bell - 0.52);
			// Material loss is small. Most fundamental loss occurs when an
			// unbalanced fork drives its resilient mounting block and harp. The
			// support coordinate has been analytically reduced into the reaction
			// factors above; anti-phase tine/bar motion cancels that force and raises
			// Q, while an isolated tine rings briefly. The keyboard term represents
			// increasing mounting loss toward the short, stiff treble assemblies.
			const double intrinsicForkLifetime = IntrinsicForkDecaySeconds *
				keyboardDecayScale;
			const double supportLossRate = 0.55 + 4.45 * keyPosition_;
			for (std::size_t index = 0; index < modes_.size(); ++index)
			{
				double lifetime;
				if (index < 2)
				{
					const double lossRate = 1.0 / intrinsicForkLifetime +
						supportLossRate * coupledSupportLossFactor_[index];
					lifetime = decayScale / std::max(1.0e-9, lossRate);
				}
				else if (index == ElectricPianoToneBarSubMode)
					lifetime = ToneBarSubModeDecaySeconds(keyPosition_) *
						decayScale;
				else
					lifetime = (index == 2 ? TransverseModeDecaySeconds *
						keyboardDecayScale :
						AttackModeDecaySeconds(index, keyPosition_)) *
						decayScale;
				if (index >= ElectricPianoAttackModeBegin &&
					index < ElectricPianoAttackModeEnd)
					lifetime *= bellDecay;
				if (damped)
					lifetime = std::min(lifetime, releaseSeconds *
						(index >= ElectricPianoAttackModeBegin &&
							index < ElectricPianoAttackModeEnd ? 0.55 : 1.0));
				modeRadius_[index] = std::exp(-1.0 /
					(std::max(0.002, lifetime) * sampleRate_));
				modePickupRadius_[index] = std::exp(-1.0 /
					(std::max(0.002, lifetime) * sampleRate_ *
						static_cast<double>(PickupOversamplingFactor)));
				modeSubstepRadius_[index] = std::exp(-1.0 /
					(std::max(0.002, lifetime) * sampleRate_ *
						static_cast<double>(ContactOversamplingFactor)));
			}
		}

		for (std::size_t index = 0; index < modes_.size(); ++index)
		{
			const double frequency = fundamental * modeRatio_[index];
			const double normalizedFrequency = frequency / sampleRate_;
			modeBandlimitGain_[index] = ElectricPianoModeBandlimitGain(
				frequency, sampleRate_);
			modeActive_[index] = normalizedFrequency < 0.49 &&
				modeBandlimitGain_[index] > 0.0;
			// A mode need not be audible to load the hammer. The contact solve runs
			// at 64x host rate, so retain the already modelled ultrasonic modes in
			// the point mobility during collision while keeping them out of pickup
			// reconstruction. The old single activity flag removed most higher
			// modes from upper-key collisions and concentrated their energy in the
			// last audible coordinate. This split changes no bass mode that was
			// already below the pickup bandwidth.
			if (contactActive_)
				contactModeActive_[index] = frequency < 0.49 * sampleRate_ *
					static_cast<double>(ContactOversamplingFactor);
			// Contact-only coefficients are needed only during the short collision.
			// Avoid paying their trigonometric update cost on every pitch-modulated
			// free-decay sample after the hammer has separated.
			const bool collisionModeNeeded = contactActive_ &&
				contactModeActive_[index];
			if (!modeActive_[index] && !collisionModeNeeded)
			{
				modes_[index] = {};
				modeRealCoefficient_[index] = 0.0;
				modeImaginaryCoefficient_[index] = 0.0;
				modePickupRealCoefficient_[index] = 0.0;
				modePickupImaginaryCoefficient_[index] = 0.0;
				modeSubstepRealCoefficient_[index] = 0.0;
				modeSubstepImaginaryCoefficient_[index] = 0.0;
				modeAngularFrequency_[index] = 0.0;
				continue;
			}
			const double angle = TwoPi * frequency / sampleRate_;
			if (modeActive_[index])
			{
				modeRealCoefficient_[index] = modeRadius_[index] * std::cos(angle);
				modeImaginaryCoefficient_[index] = modeRadius_[index] * std::sin(angle);
				const double pickupAngle = angle /
					static_cast<double>(PickupOversamplingFactor);
				modePickupRealCoefficient_[index] = modePickupRadius_[index] *
					std::cos(pickupAngle);
				modePickupImaginaryCoefficient_[index] = modePickupRadius_[index] *
					std::sin(pickupAngle);
			}
			else
			{
				modeRealCoefficient_[index] = 0.0;
				modeImaginaryCoefficient_[index] = 0.0;
				modePickupRealCoefficient_[index] = 0.0;
				modePickupImaginaryCoefficient_[index] = 0.0;
				if (!contactActive_)
					modes_[index] = {};
			}
			const double substepAngle = angle /
				static_cast<double>(ContactOversamplingFactor);
			modeSubstepRealCoefficient_[index] = modeSubstepRadius_[index] *
				std::cos(substepAngle);
			modeSubstepImaginaryCoefficient_[index] = modeSubstepRadius_[index] *
				std::sin(substepAngle);
			modeAngularFrequency_[index] = TwoPi * frequency;
		}

		if (transformCoupledState && modeActive_[0] && modeActive_[1])
		{
			const double ratioDifference = coupledToneBarRatio_[0] -
				coupledToneBarRatio_[1];
			if (std::abs(ratioDifference) > 1.0e-9)
			{
				modes_[0].real = (oldToneBarPosition -
					coupledToneBarRatio_[1] * oldTinePosition) /
					ratioDifference;
				modes_[1].real = oldTinePosition - modes_[0].real;
				const double modeZeroVelocity = (oldToneBarVelocity -
					coupledToneBarRatio_[1] * oldTineVelocity) /
					ratioDifference;
				const double modeOneVelocity = oldTineVelocity - modeZeroVelocity;
				modes_[0].imaginary = -modeZeroVelocity /
					modeAngularFrequency_[0];
				modes_[1].imaginary = -modeOneVelocity /
					modeAngularFrequency_[1];
			}
		}
		coupledForkInitialized_ = true;
		cachedFundamental_ = fundamental;
		cachedDecay_ = decay;
		cachedBell_ = bell;
		cachedCoupling_ = coupling;
		cachedRelease_ = release;
		cachedDamped_ = damped;
		coefficientsDirty_ = false;
	}

