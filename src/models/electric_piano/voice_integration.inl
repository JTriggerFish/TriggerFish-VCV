	void AdvanceFreeModes()
	{
		if (!modulationPathActive_)
		{
			for (std::size_t index = 0;
				index < ElectricPianoToneBarSubMode; ++index)
			{
				if (!modeActive_[index])
					continue;
				const double oldReal = modes_[index].real;
				const double oldImaginary = modes_[index].imaginary;
				modes_[index].real = modeRealCoefficient_[index] * oldReal -
					modeImaginaryCoefficient_[index] * oldImaginary;
				modes_[index].imaginary = modeImaginaryCoefficient_[index] *
					oldReal + modeRealCoefficient_[index] * oldImaginary;
			}
			constexpr std::size_t subMode = ElectricPianoToneBarSubMode;
			if (modeActive_[subMode])
			{
				const double oldReal = modes_[subMode].real;
				const double oldImaginary = modes_[subMode].imaginary;
				modes_[subMode].real = modeRealCoefficient_[subMode] * oldReal -
					modeImaginaryCoefficient_[subMode] * oldImaginary;
				modes_[subMode].imaginary = modeImaginaryCoefficient_[subMode] *
					oldReal + modeRealCoefficient_[subMode] * oldImaginary;
			}
			return;
		}
		// The pickup nonlinearity and creative FM/PM are evaluated at 4x. Advance
		// the physical resonators at that same cadence so their modal coordinates
		// do not have to be reconstructed from a host-rate aggregate signal.
		for (int frame = 0; frame < PickupOversamplingFactor; ++frame)
		{
			for (std::size_t index = 0;
				index < ElectricPianoToneBarSubMode; ++index)
			{
				if (!modeActive_[index])
				{
					pickupModeFrames_[frame][index] = {};
					continue;
				}
				const double oldReal = modes_[index].real;
				const double oldImaginary = modes_[index].imaginary;
				modes_[index].real = modePickupRealCoefficient_[index] * oldReal -
					modePickupImaginaryCoefficient_[index] * oldImaginary;
				modes_[index].imaginary = modePickupImaginaryCoefficient_[index] *
					oldReal + modePickupRealCoefficient_[index] * oldImaginary;
				pickupModeFrames_[frame][index] = modes_[index];
			}
			constexpr std::size_t subMode = ElectricPianoToneBarSubMode;
			if (!modeActive_[subMode])
				pickupModeFrames_[frame][subMode] = {};
			else
			{
				const double oldReal = modes_[subMode].real;
				const double oldImaginary = modes_[subMode].imaginary;
				modes_[subMode].real = modePickupRealCoefficient_[subMode] * oldReal -
					modePickupImaginaryCoefficient_[subMode] * oldImaginary;
				modes_[subMode].imaginary = modePickupImaginaryCoefficient_[subMode] *
					oldReal + modePickupRealCoefficient_[subMode] * oldImaginary;
				pickupModeFrames_[frame][subMode] = modes_[subMode];
			}
		}
	}

	void AdvanceCoupledHammerAndModes()
	{
		// A free hammer has no preload that could sustain contact, so a valid
		// collision must end through physical separation. Do not impose an elapsed-
		// time cutoff: soft bass strikes legitimately remain compressed longest.
		// Only non-finite state is treated as a numerical failure below.
		const double timeStep = 1.0 /
			(sampleRate_ * static_cast<double>(ContactOversamplingFactor));
		auto contactPointState = [&]()
		{
			double position = 0.0;
			double velocity = 0.0;
			for (std::size_t index = 0; index < modes_.size(); ++index)
			{
				if (!contactModeActive_[index])
					continue;
				position -= contactModeShape_[index] * modes_[index].real;
				velocity += contactModeShape_[index] *
					modeAngularFrequency_[index] * modes_[index].imaginary;
			}
			return std::array<double, 2>{position, velocity};
		};
		if (modulationPathActive_)
			for (auto& frame : pickupModeFrames_)
				frame.fill({});
		for (int substep = 0; substep < ContactOversamplingFactor; ++substep)
		{
			if (contactActive_)
			{
				const auto tine = contactPointState();
				const double compression = hammerPosition_ - tine[0];
				const double relativeVelocity = hammerVelocity_ - tine[1];
				if (!std::isfinite(compression) ||
					!std::isfinite(relativeVelocity))
				{
					for (auto& mode : modes_)
						mode = {};
					hammerPosition_ = 0.0;
					hammerVelocity_ = 0.0;
					contactActive_ = false;
					contactEngaged_ = false;
				}
				else
				{
					double force = 0.0;
					if (compression > 0.0)
					{
						contactEngaged_ = true;
						// Hunt--Crossley damping may reach zero on a sufficiently fast
						// unloading stroke, but a contact cannot pull the hammer back
						// toward the tine.  Do not cap the measured loading branch.
						const double hysteresis =
							1.0 + contactLoss_ * relativeVelocity;
						if (hysteresis <= 0.0 && relativeVelocity < 0.0)
						{
							// The dissipative fit has crossed the point at which its
							// algebraic force would become tensile. A real neoprene tip
							// cannot pull the receding hammer, so this is separation,
							// not a zero-force state held at positive penetration.
							hammerPosition_ = tine[0];
							contactActive_ = false;
							contactEngaged_ = false;
						}
						else
							force = contactStiffness_ *
								std::pow(compression, contactExponent_) * hysteresis;
					}
					if (!std::isfinite(force))
					{
						for (auto& mode : modes_)
							mode = {};
						hammerPosition_ = 0.0;
						hammerVelocity_ = 0.0;
						contactActive_ = false;
						contactEngaged_ = false;
					}
					else
					{
						const double impulse = force * timeStep;
						hammerVelocity_ -= impulse * hammerInverseMass_;
						for (std::size_t index = 0; index < modes_.size(); ++index)
						{
							if (!contactModeActive_[index] ||
								!(modeAngularFrequency_[index] > 0.0))
								continue;
							modes_[index].imaginary += contactModeShape_[index] *
								modeInverseMass_[index] * impulse /
								modeAngularFrequency_[index];
						}
						hammerPosition_ += timeStep * hammerVelocity_;
					}
				}
			}

			for (std::size_t index = 0; index < modes_.size(); ++index)
			{
				if (!contactModeActive_[index])
					continue;
				const double oldReal = modes_[index].real;
				const double oldImaginary = modes_[index].imaginary;
				modes_[index].real = modeSubstepRealCoefficient_[index] *
					oldReal - modeSubstepImaginaryCoefficient_[index] *
					oldImaginary;
				modes_[index].imaginary = modeSubstepImaginaryCoefficient_[index] *
					oldReal + modeSubstepRealCoefficient_[index] * oldImaginary;
			}
			if (modulationPathActive_ && (substep + 1) %
				(ContactOversamplingFactor / PickupOversamplingFactor) == 0)
			{
				const int frame = (substep + 1) /
					(ContactOversamplingFactor / PickupOversamplingFactor) - 1;
				for (std::size_t index = 0; index < modes_.size(); ++index)
					pickupModeFrames_[frame][index] = modeActive_[index] ?
						modes_[index] : Mode{};
			}

			if (contactActive_)
			{
				contactAge_ += timeStep;
				const auto tine = contactPointState();
				const double compression = hammerPosition_ - tine[0];
				const double relativeVelocity = hammerVelocity_ - tine[1];
				if (!std::isfinite(hammerPosition_) ||
					!std::isfinite(hammerVelocity_) ||
					!std::isfinite(compression) ||
					!std::isfinite(relativeVelocity))
				{
					for (auto& mode : modes_)
						mode = {};
					hammerPosition_ = 0.0;
					hammerVelocity_ = 0.0;
					contactActive_ = false;
					contactEngaged_ = false;
				}
				else if (contactEngaged_ && compression <= 0.0 &&
					relativeVelocity < 0.0)
				{
					contactActive_ = false;
					contactEngaged_ = false;
				}
			}
		}
		// Contact-only coordinates cannot be advanced at host rate once the
		// hammer separates. Their remaining ultrasonic energy is outside the
		// pickup model and is intentionally discarded here, after it has affected
		// the complete physical collision rather than before it.
		if (!contactActive_)
			for (std::size_t index = 0; index < modes_.size(); ++index)
				if (!modeActive_[index])
					modes_[index] = {};
	}

