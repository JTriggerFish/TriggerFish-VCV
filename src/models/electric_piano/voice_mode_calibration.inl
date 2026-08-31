	static double KeyboardModeRatio(std::size_t index,
		double keyboardPosition)
	{
		if (index < 3)
			return 1.0;
		const std::size_t waveNumberIndex = index - 2;
		const double uniformRatio = std::pow(
			CantileverWaveNumbers[waveNumberIndex] /
			CantileverWaveNumbers[0], 2.0);
		auto applyShearCorrection = [&](double ratio)
		{
			// The measured F1/F3 ratios already contain the real tine's shear and
			// rotary-inertia effects. Above the final measurement, correct the
			// Euler--Bernoulli extrapolation with the first-order Timoshenko/Rayleigh
			// factor cited by Pfeifle and Sonderbo. It is normalized by the
			// fundamental so pitch is unchanged; only the increasingly sensitive
			// upper-mode ratios of short tines move. Blend in only beyond F3 so the
			// measured anchors and the successful bass remain exact.
			constexpr double PoissonRatio = 0.29;
			constexpr double CircularShearCoefficient =
				6.0 * (1.0 + PoissonRatio) / (7.0 + 6.0 * PoissonRatio);
			const double youngsToShearRatio = 2.0 * (1.0 + PoissonRatio);
			const double tineLength =
				ElectricPianoPublishedMechanicalData::LongestTineMetres * std::pow(
					ElectricPianoPublishedMechanicalData::ShortestTineMetres /
						ElectricPianoPublishedMechanicalData::LongestTineMetres,
					std::clamp(keyboardPosition, 0.0, 1.0));
			const double radiusToLength =
				ElectricPianoPublishedMechanicalData::TineRadiusMetres /
				std::max(1.0e-9, tineLength);
			const double correctionCoefficient = 0.25 * radiusToLength *
				radiusToLength * (1.0 + youngsToShearRatio /
					CircularShearCoefficient);
			const double correction = std::sqrt(
				(1.0 + correctionCoefficient * CantileverWaveNumbers[0] *
					CantileverWaveNumbers[0]) /
				(1.0 + correctionCoefficient *
					CantileverWaveNumbers[waveNumberIndex] *
					CantileverWaveNumbers[waveNumberIndex]));
			constexpr double F3Position = (53.0 - 28.0) / 72.0;
			constexpr double A4Position = (69.0 - 28.0) / 72.0;
			const double blendLinear = std::clamp((keyboardPosition - F3Position) /
				(A4Position - F3Position), 0.0, 1.0);
			const double blend = blendLinear * blendLinear *
				(3.0 - 2.0 * blendLinear);
			return ratio * std::exp(blend * std::log(correction));
		};
		// Gabrielli/Cantarini report the first measured overtones of F1 at
		// 7.2/20.6 and those of F3 at 7.4/20.7/38.7. A Rayleigh--Ritz cantilever with
		// one sliding point mass fits both anchors with tuning-mass fractions
		// 0.146/0.177 and positions 0.859L/0.840L. Extrapolating those physical
		// parameters (mass fraction versus tine length and spring position versus
		// key) yields the compact ratio polynomials below. This replaces the former
		// return to an unloaded A4 beam even though every real key retains a tuning
		// spring. The fit is explicitly provisional above F3 until the article's
		// full-key data or new SLDV measurements are available. The very lowest key
		// retains the calibrated 7.11/20.25 ratios.
		constexpr double F1Position = (29.0 - 28.0) / 72.0;
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		auto smoothInterpolate = [](double position, double lowerPosition,
			double upperPosition, double lowerValue, double upperValue)
		{
			const double linear = std::clamp((position - lowerPosition) /
				(upperPosition - lowerPosition), 0.0, 1.0);
			const double smooth = linear * linear * (3.0 - 2.0 * linear);
			return lowerValue + smooth * (upperValue - lowerValue);
		};
		if (waveNumberIndex >= 1 && waveNumberIndex <= 3)
		{
			const std::size_t measuredIndex = waveNumberIndex - 1;
			constexpr std::array<double, 3> LowestKeyRatios{
				7.11, 20.25, 34.386};
			constexpr std::array<double, 3> F1Ratios{7.2, 20.6, 34.386};
			constexpr std::array<double, 3> F3Ratios{7.4, 20.7, 38.7};
			if (keyboardPosition < F1Position)
				return applyShearCorrection(smoothInterpolate(keyboardPosition, 0.0,
					F1Position, LowestKeyRatios[measuredIndex],
					F1Ratios[measuredIndex]));
			if (keyboardPosition < F3Position)
				return applyShearCorrection(smoothInterpolate(keyboardPosition,
					F1Position, F3Position, F1Ratios[measuredIndex],
					F3Ratios[measuredIndex]));
			const double pointMassPosition = std::clamp(
				(keyboardPosition - F3Position) / (1.0 - F3Position), 0.0, 1.0);
			const double squared = pointMassPosition * pointMassPosition;
			constexpr std::array<double, 3> LinearPointMassCoefficients{
				0.3887, 0.2667, -1.1750};
			constexpr std::array<double, 3> QuadraticPointMassCoefficients{
				0.0022, -0.3675, 0.6190};
			return applyShearCorrection(F3Ratios[measuredIndex] +
				LinearPointMassCoefficients[measuredIndex] * pointMassPosition +
				QuadraticPointMassCoefficients[measuredIndex] * squared);
		}
		return applyShearCorrection(uniformRatio);
	}

	static std::array<double, 9> SpringLoadedModeReduction(std::size_t index,
		double keyboardPosition)
	{
		// Offline Rayleigh--Ritz reduction of a uniform cantilever plus the sliding
		// tuning mass fitted at F1/F3. Components 0--7 expand a tip-normalized loaded
		// eigenfunction in the first eight tip-normalized cantilever shapes; component
		// 8 is generalized modal mass divided by the quarter-beam reference mass.
		// Fifth-order keyed polynomials reproduce the 73-key eigensolve with maximum
		// coefficient error below 6.5e-6 and mass-ratio error below 2.3e-5. Keeping
		// this solve out of the audio path also preserves cheap continuous tuning.
		constexpr std::array<std::array<std::array<double, 6>, 9>, 3>
			ReductionPolynomials{{
			{{
				{{-0.0956129, 0.0456194, 0.00777974, -0.000116554,
					-9.46565e-05, -3.09906e-06}},
				{{1.10324, -0.0455643, -0.0104353, -9.22499e-05,
					0.000148787, 1.19656e-05}},
				{{-0.00398375, -0.00128955, 0.00183961, 0.000294257,
					-1.9185e-05, -8.84302e-06}},
				{{-0.00220983, 0.000352183, 0.000614798, 2.54259e-05,
					-1.9119e-05, -2.65702e-06}},
				{{-0.000963916, 0.000374115, 0.000208078, -2.96839e-05,
					-1.20937e-05, -3.05088e-07}},
				{{-0.000370971, 0.000262402, 4.75427e-05, -3.62169e-05,
					-5.54762e-06, 7.79864e-07}},
				{{-0.000102739, 0.000160785, -1.67049e-05, -2.82371e-05,
					-6.24785e-07, 1.14076e-06}},
				{{8.68671e-06, 8.48909e-05, -3.78117e-05, -1.67413e-05,
					2.43924e-06, 1.01795e-06}},
				{{1.24536, -0.130942, -0.0154338, 0.00239892,
					0.00042329, -1.50138e-05}}
			}},
			{{
				{{0.0856358, 0.0718697, 0.00239745, -0.000800856,
					0.000103003, 2.45192e-05}},
				{{0.0264046, 0.00725223, -0.0122908, -0.00115264,
					9.91097e-05, -1.35701e-05}},
				{{0.856621, -0.10777, 0.0105401, 0.00509314,
					-0.000285909, -0.000179625}},
				{{0.0204186, 0.0217557, 0.00242873, -0.00193329,
					-0.00012666, 0.000105396}},
				{{0.00748901, 0.00635952, -0.000622121, -0.00081385,
					2.26418e-05, 4.0268e-05}},
				{{0.00274673, 0.00146424, -0.00106307, -0.000356234,
					6.33963e-05, 1.89727e-05}},
				{{0.000746877, -0.000241582, -0.00086532, -9.23743e-05,
					6.8684e-05, 6.35504e-06}},
				{{-6.27144e-05, -0.000689937, -0.000524997, 5.61063e-05,
					5.57334e-05, -2.31667e-06}},
				{{0.758331, -0.147027, 0.0427399, 0.00513243,
					-0.00150902, -5.12784e-05}}
			}},
			{{
				{{0.179255, 0.069498, -0.000834428, -9.1712e-05,
					-0.000210303, -0.000145463}},
				{{0.0500679, -0.00904389, -0.0123546, -0.000331001,
					5.78134e-07, 7.40167e-07}},
				{{-0.0772817, -0.0930271, -0.0253526, 0.000664447,
					0.000802763, 0.000230014}},
				{{0.750814, 0.0193002, 0.0563247, -0.000129465,
					-0.00115623, 0.000269954}},
				{{0.0702597, 0.022177, -0.00752142, -0.000515268,
					0.000267363, -0.000236803}},
				{{0.021748, 0.000528111, -0.00541273, -0.000221897,
					0.000127871, -7.32499e-05}},
				{{0.00559402, -0.00461574, -0.00332659, 0.000191674,
					0.000100441, -2.94211e-05}},
				{{-0.000457813, -0.00481666, -0.00152236, 0.000433223,
					6.75157e-05, -1.57717e-05}},
				{{0.680285, 0.112127, 0.101452, 0.00504365,
					0.00165225, -7.54704e-05}}
			}}
		}};
		std::array<double, 9> reduction{};
		if (index < 3 || index > 5)
		{
			reduction[8] = 1.0;
			return reduction;
		}
		const double coordinate = 2.0 * std::clamp(keyboardPosition, 0.0, 1.0) - 1.0;
		const auto& polynomials = ReductionPolynomials[index - 3];
		for (std::size_t component = 0; component < reduction.size(); ++component)
		{
			double value = polynomials[component].back();
			for (std::size_t term = polynomials[component].size() - 1;
				term-- > 0;)
				value = value * coordinate + polynomials[component][term];
			reduction[component] = value;
		}
		reduction[8] = std::clamp(reduction[8], 0.5, 1.5);
		return reduction;
	}

	static double AttackModeModalMassMultiplier(std::size_t index,
		double keyboardPosition)
	{
		// Modes 3--5 use the spring-loaded Rayleigh--Ritz eigenvectors fitted to
		// the measured F1/F3 ratios. Their generalized masses therefore differ from
		// the quarter-beam mass used by KeyProfile. Higher, unconstrained modes keep
		// the uniform-beam normalization until measurements justify loading them.
		if (index < 3 || index > 5)
			return 1.0;
		return SpringLoadedModeReduction(index, keyboardPosition)[8];
	}

	static double BeamAttackModeDecaySeconds(std::size_t index,
		double keyboardPosition)
	{
		const double tineLength =
			ElectricPianoPublishedMechanicalData::LongestTineMetres * std::pow(
				ElectricPianoPublishedMechanicalData::ShortestTineMetres /
				ElectricPianoPublishedMechanicalData::LongestTineMetres,
				std::clamp(keyboardPosition, 0.0, 1.0));
		const std::size_t waveNumberIndex = index - 2;
		const double spatialWaveNumber =
			CantileverWaveNumbers[waveNumberIndex] / tineLength;
		const double lossRate =
			ElectricPianoPublishedMechanicalData::FrequencyIndependentLossPerSecond +
			ElectricPianoPublishedMechanicalData::
				FrequencyDependentLossSquareMetresPerSecond *
			spatialWaveNumber * spatialWaveNumber;
		return 1.0 / std::max(1.0e-9, lossRate);
	}

	static double AttackModeDecaySeconds(std::size_t index,
		double keyboardPosition)
	{
		const double beamDecay = BeamAttackModeDecaySeconds(index,
			keyboardPosition);
		const std::size_t waveNumberIndex = index - 2;
		if (waveNumberIndex < 1 || waveNumberIndex > 3)
			return beamDecay;

		// Table 8.6 of Cantarini/Gabrielli gives amplitude-decay slopes in
		// dB/s. Convert with tau = 20 log10(e)/|slope| and use the resulting
		// measured lifetime as a multiplier on Sonderbo's distributed-loss law.
		// Retaining the beam law between anchors preserves register scaling; the
		// multiplier admits the observed non-monotonic mode-dependent losses that a
		// single sigma0 + sigma1*k^2 curve cannot reproduce. The F3 multiplier must
		// not remain frozen above the last measured damping anchor: doing that made
		// modes four and five roughly eleven and five times too long throughout the
		// treble. In the absence of treble damping measurements, relax the measured
		// anomaly independently by A4 to the sourced distributed-loss law. Do not tie
		// this provisional loss coordinate to the point-mass frequency extrapolation:
		// a tuning spring can shift frequency without preserving F3's anomalous loss.
		constexpr double DbPerNeper = 8.685889638065037;
		constexpr double F1Position = (29.0 - 28.0) / 72.0;
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		constexpr std::array<double, 3> F1DecaySeconds{
			DbPerNeper / 21.1, DbPerNeper / 67.7, 0.0};
		constexpr std::array<double, 3> F3DecaySeconds{
			DbPerNeper / 294.0, DbPerNeper / 37.0, DbPerNeper / 161.0};
		const std::size_t measuredIndex = waveNumberIndex - 1;
		const double f1Multiplier = F1DecaySeconds[measuredIndex] > 0.0 ?
			F1DecaySeconds[measuredIndex] /
				BeamAttackModeDecaySeconds(index, F1Position) : 1.0;
		const double f3Multiplier = F3DecaySeconds[measuredIndex] /
			BeamAttackModeDecaySeconds(index, F3Position);
		if (keyboardPosition > F3Position)
		{
			constexpr double A4Position = (69.0 - 28.0) / 72.0;
			const double relaxationLinear = std::clamp(
				(keyboardPosition - F3Position) / (A4Position - F3Position),
				0.0, 1.0);
			const double relaxation = relaxationLinear * relaxationLinear *
				(3.0 - 2.0 * relaxationLinear);
			const double measuredLossPerturbation = 1.0 - relaxation;
			return beamDecay * std::exp(measuredLossPerturbation *
				std::log(std::max(1.0e-12, f3Multiplier)));
		}
		const double linear = std::clamp((keyboardPosition - F1Position) /
			(F3Position - F1Position), 0.0, 1.0);
		const double smooth = linear * linear * (3.0 - 2.0 * linear);
		const double logMultiplier = std::log(std::max(1.0e-12, f1Multiplier)) +
			smooth * (std::log(std::max(1.0e-12, f3Multiplier)) -
				std::log(std::max(1.0e-12, f1Multiplier)));
		return beamDecay * std::exp(logMultiplier);
	}

