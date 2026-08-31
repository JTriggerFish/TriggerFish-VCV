	static double AttackModeEnergyCalibration(std::size_t index,
		double keyboardPosition)
	{
		// The measured F1/F3 decay constants describe modal residues whose
		// magnitudes were fitted at the same time. Substituting only those decay
		// constants changes a mode's total observed energy in direct proportion to
		// its lifetime. Preserve the existing reduced model's broadband impulse
		// calibration by scaling the contact projection by sqrt(tau_beam/tau_fit).
		// This is an explicit residue normalization; it does not alter the loaded
		// eigenfunction, generalized mass, frequency, or measured decay constant.
		const double fittedLifetime = AttackModeDecaySeconds(index,
			keyboardPosition);
		const double beamLifetime = BeamAttackModeDecaySeconds(index,
			keyboardPosition);
		return std::sqrt(std::clamp(beamLifetime /
			std::max(1.0e-9, fittedLifetime), 0.05, 4.0));
	}

	static double ToneBarSubModeFrequencyRatio(double keyboardPosition)
	{
		// Gabrielli/Cantarini SLDV measurements locate the tone-bar submode at
		// 0.83*f0 for F1 and 0.58*f0 for F3. The openly available measurements
		// stop there. Above F3, aligned direct-harp recordings show the strongest
		// early submode settling near 0.4--0.5*f0 rather than rising toward f0.
		constexpr double F1Position = (29.0 - 28.0) / 72.0;
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		if (keyboardPosition <= F1Position)
			return 0.83;
		if (keyboardPosition <= F3Position)
		{
			const double linear = (keyboardPosition - F1Position) /
				(F3Position - F1Position);
			const double smooth = linear * linear * (3.0 - 2.0 * linear);
			return std::exp(std::log(0.83) + smooth *
				(std::log(0.58) - std::log(0.83)));
		}
		const double upperPosition = (keyboardPosition - F3Position) /
			(1.0 - F3Position);
		return 0.42 + 0.16 * std::exp(-5.0 * upperPosition);
	}

	static double ToneBarSubModeDecaySeconds(double keyboardPosition)
	{
		// Table 8.6 gives amplitude slopes of -9.1 dB/s at F1 and
		// -138 dB/s at F3. Convert with tau=20*log10(e)/|slope|. The
		// upper continuation is constrained by direct-harp spectrograms, where
		// the line falls by roughly 25--40 dB between 100 and 250 ms.
		constexpr double DbPerNeper = 8.685889638065037;
		constexpr double F1Lifetime = DbPerNeper / 9.1;
		constexpr double F3Lifetime = DbPerNeper / 138.0;
		constexpr double F1Position = (29.0 - 28.0) / 72.0;
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		if (keyboardPosition <= F1Position)
			return F1Lifetime;
		if (keyboardPosition <= F3Position)
		{
			const double linear = (keyboardPosition - F1Position) /
				(F3Position - F1Position);
			const double smooth = linear * linear * (3.0 - 2.0 * linear);
			return std::exp(std::log(F1Lifetime) + smooth *
				(std::log(F3Lifetime) - std::log(F1Lifetime)));
		}
		const double upperPosition = (keyboardPosition - F3Position) /
			(1.0 - F3Position);
		return F3Lifetime * std::pow(0.035 / F3Lifetime, upperPosition);
	}

	static double ToneBarSubModeTineParticipation(double keyboardPosition)
	{
		// A normal mode is observed and forced through the same tine component.
		// The two SLDV anchors imply a weak F1 residue (-28 dB relative to the
		// played mode) and a much stronger F3 residue (-4.5 dB). Above F3 the
		// direct-harp corpus shows a rapid reduction; by the upper register the
		// submode and its pickup sidebands are 40--60 dB below the fundamental.
		// These participation values are an initial reciprocal modal fit; unlike
		// the former observation gain they enter twice, at force and pickup.
		constexpr double F1Position = (29.0 - 28.0) / 72.0;
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		constexpr std::array<double, 8> KeyPositions{
			F1Position, F3Position, (55.0 - 28.0) / 72.0,
			(59.0 - 28.0) / 72.0, (65.0 - 28.0) / 72.0,
			(76.0 - 28.0) / 72.0, (88.0 - 28.0) / 72.0, 1.0};
		constexpr std::array<double, 8> Participations{
			0.177, 0.56, 0.40, 0.25, 0.15, 0.060, 0.045, 0.040};
		if (keyboardPosition <= KeyPositions.front())
			return Participations.front();
		for (std::size_t upper = 1; upper < KeyPositions.size(); ++upper)
		{
			if (keyboardPosition > KeyPositions[upper])
				continue;
			const double linear = (keyboardPosition - KeyPositions[upper - 1]) /
				(KeyPositions[upper] - KeyPositions[upper - 1]);
			const double smooth = linear * linear * (3.0 - 2.0 * linear);
			return std::exp(std::log(Participations[upper - 1]) + smooth *
				(std::log(Participations[upper]) -
					std::log(Participations[upper - 1])));
		}
		return Participations.back();
	}

	static double StrikePositionFromControl(double keyboardPosition,
		double strikeControl)
	{
		// The service procedure does not describe a straight line: it asks the
		// technician to set C4, F3 and C3 independently, then accept the intervening
		// keys when they fall close to maximum power without a thunk. A linear
		// 0.40L--0.20L approximation put C4 and much of the upper register close to
		// higher-mode nodes. A full-keyboard render sweep found the same optimum the
		// manual describes: a fairly flat bass line, a pronounced bend around the
		// C3--C4 checkpoints, and a gentler approach toward the treble clamp.
		//
		// Positions remain explicit listening trims because no published factory
		// jig dimensions identify the contact point on every tine. Keeping the six
		// checkpoints here makes a later measured-harp calibration local and avoids
		// hiding the correction in per-mode gains.
		constexpr std::array<double, 6> KeyboardCheckpoints{
			0.0, (48.0 - 28.0) / 72.0, (53.0 - 28.0) / 72.0,
			(60.0 - 28.0) / 72.0, (84.0 - 28.0) / 72.0, 1.0};
		constexpr std::array<double, 6> FactoryStrikePositions{
			0.38, 0.29, 0.22, 0.205, 0.16, 0.14};
		const double position = std::clamp(keyboardPosition, 0.0, 1.0);
		std::size_t lower = 0;
		while (lower + 2 < KeyboardCheckpoints.size() &&
			position > KeyboardCheckpoints[lower + 1])
			++lower;
		const double interval = KeyboardCheckpoints[lower + 1] -
			KeyboardCheckpoints[lower];
		const double linearFraction = std::clamp(
			(position - KeyboardCheckpoints[lower]) / interval, 0.0, 1.0);
		const double smoothFraction = linearFraction * linearFraction *
			(3.0 - 2.0 * linearFraction);
		const double factoryPosition = FactoryStrikePositions[lower] +
			smoothFraction * (FactoryStrikePositions[lower + 1] -
				FactoryStrikePositions[lower]);
		// The panel now represents a signed physical displacement of the harp's
		// striking line about the factory point. A full-range movement tapers from
		// 6 mm in the bass to 1 mm in the treble, consistent with the much shorter
		// upper tines. Unlike interpolation toward 0.04L/0.96L, this is smooth at
		// centre and cannot fling adjacent high-register settings across several
		// modal nodes. Centre remains bit-for-bit the calibrated factory position.
		const double tineLength =
			ElectricPianoPublishedMechanicalData::LongestTineMetres * std::pow(
				ElectricPianoPublishedMechanicalData::ShortestTineMetres /
					ElectricPianoPublishedMechanicalData::LongestTineMetres,
				position);
		const double maximumOffsetMetres = 0.006 + position * (0.001 - 0.006);
		const double signedControl = 2.0 *
			(std::clamp(strikeControl, 0.0, 1.0) - 0.5);
		return std::clamp(factoryPosition + signedControl *
			maximumOffsetMetres / tineLength, 0.04, 0.96);
	}

	static double UniformCantileverModeShape(std::size_t waveNumberIndex,
		double position)
	{
		const double beta = CantileverWaveNumbers[waveNumberIndex];
		const double sigma = (std::cosh(beta) + std::cos(beta)) /
			(std::sinh(beta) + std::sin(beta));
		auto shape = [&](double normalizedPosition)
		{
			const double argument = beta * normalizedPosition;
			return std::cosh(argument) - std::cos(argument) - sigma *
				(std::sinh(argument) - std::sin(argument));
		};
		return std::clamp(shape(std::clamp(position, 0.0, 1.0)) / shape(1.0),
			-1.5, 1.5);
	}

	static double CantileverModeShape(std::size_t index, double position)
	{
		const std::size_t waveNumberIndex = index < 3 ? 0 : index - 2;
		return UniformCantileverModeShape(waveNumberIndex, position);
	}

	static double FiniteContactModeProjection(std::size_t index,
		double keyboardPosition, double centre, double contactWidthMetres)
	{
		// Integrate the calibrated physical contact strip over the modal shape. The
		// symmetric binomial weights approximate the peaked pressure of a compliant
		// contact; this runs once per strike, not in the contact substeps.
		const double tineLength =
			ElectricPianoPublishedMechanicalData::LongestTineMetres * std::pow(
				ElectricPianoPublishedMechanicalData::ShortestTineMetres /
					ElectricPianoPublishedMechanicalData::LongestTineMetres,
				std::clamp(keyboardPosition, 0.0, 1.0));
		const double halfWidth = std::clamp(0.5 * contactWidthMetres /
			std::max(1.0e-6, tineLength), 0.001, 0.14);
		constexpr std::array<double, 5> Offsets{-1.0, -0.5, 0.0, 0.5, 1.0};
		constexpr std::array<double, 5> Weights{1.0, 4.0, 6.0, 4.0, 1.0};
		const bool springLoadedMode = index >= 3 && index <= 5;
		const auto loadedReduction = springLoadedMode ?
			SpringLoadedModeReduction(index, keyboardPosition) :
			std::array<double, 9>{};
		auto modeShape = [&](double position)
		{
			if (!springLoadedMode)
				return CantileverModeShape(index, position);
			double shape = 0.0;
			for (std::size_t basis = 0; basis < 8; ++basis)
				shape += loadedReduction[basis] *
					UniformCantileverModeShape(basis, position);
			return std::clamp(shape, -1.5, 1.5);
		};
		double projection = 0.0;
		for (std::size_t sample = 0; sample < Offsets.size(); ++sample)
			projection += Weights[sample] * modeShape(
				centre + halfWidth * Offsets[sample]);
		return projection / 16.0;
	}

	struct MagneticFieldSample
	{
		double vertical{};
		double radial{};
	};

