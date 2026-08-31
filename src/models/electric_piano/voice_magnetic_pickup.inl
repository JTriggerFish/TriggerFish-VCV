	static MagneticFieldSample EvaluateCalibratedPoleGradient(double vertical,
		double radialDistance)
	{
		// Pfeifle derives a three-part integral over an idealised pole
		// cross-section, but publishes neither production dimensions nor a flux scan.
		// Falaize's sourced circular-pole reduction was also evaluated here; as a
		// keyboard-wide replacement for this asymmetric pole it over-produced bass H2
		// by 5--14x and reduced rather than restored the treble harmonic ladder. Keep
		// the measured-listening-calibrated asymmetric edge reduction until a real 2D
		// pickup scan exists. It remains one scalar flux construction, differentiated
		// consistently for arbitrary two-dimensional tine motion.
		const double radialFalloff = PowNegative1p3(radialDistance);
		MagneticFieldSample result;
		auto accumulateEdge = [&](double edgePosition, double weight,
			double edgeRadius)
		{
			constexpr double GapBroadening = 0.180;
			const double width = edgeRadius + GapBroadening * radialDistance;
			const double inverseWidth = 1.0 / std::max(1.0e-9, width);
			const double displacement = vertical - edgePosition;
			const double argument = displacement * inverseWidth;
			const double transition = TanhPade76(argument);
			const double transitionDerivative = 1.0 - transition * transition;
			const double amplitude = weight * radialFalloff;
			result.vertical += amplitude * transitionDerivative * inverseWidth;
			const double amplitudeDerivative = -FieldFalloff * amplitude /
				std::max(1.0e-9, radialDistance);
			const double argumentDerivative = -displacement * GapBroadening *
				inverseWidth * inverseWidth;
			result.radial += amplitudeDerivative * transition + amplitude *
				transitionDerivative * argumentDerivative;
		};
		accumulateEdge(0.34, 1.0, 0.030);
		accumulateEdge(-0.46, 0.20, 0.095);
		return result;
	}

	struct MagneticFieldTable
	{
		// The primary edge is only 0.03 mm wide at zero clearance. A dense
		// process-wide table keeps interpolation below the direct-field error bound
		// without changing the four-lookups-per-sample runtime cost.
		static constexpr std::size_t VerticalSamples = 241;
		static constexpr std::size_t RadialSamples = 135;
		static constexpr double VerticalMinimum = -1.2;
		static constexpr double VerticalMaximum = 1.8;
		static constexpr double RadialMinimum = 0.45;
		static constexpr double RadialMaximum = 3.8;
		std::array<MagneticFieldSample,
			VerticalSamples * RadialSamples> samples{};

		MagneticFieldTable()
		{
			for (std::size_t verticalIndex = 0;
				verticalIndex < VerticalSamples; ++verticalIndex)
			{
				const double vertical = VerticalMinimum +
					(VerticalMaximum - VerticalMinimum) *
						static_cast<double>(verticalIndex) /
						static_cast<double>(VerticalSamples - 1);
				for (std::size_t radialIndex = 0;
					radialIndex < RadialSamples; ++radialIndex)
				{
					const double radial = RadialMinimum +
						(RadialMaximum - RadialMinimum) *
							static_cast<double>(radialIndex) /
							static_cast<double>(RadialSamples - 1);
					samples[verticalIndex * RadialSamples + radialIndex] =
						EvaluateCalibratedPoleGradient(vertical, radial);
				}
			}
		}

		MagneticFieldSample Interpolate(double vertical, double radial) const
		{
			const double verticalCoordinate = std::clamp(
				(vertical - VerticalMinimum) * (VerticalSamples - 1) /
					(VerticalMaximum - VerticalMinimum), 0.0,
				static_cast<double>(VerticalSamples - 1));
			const double radialCoordinate = std::clamp(
				(radial - RadialMinimum) * (RadialSamples - 1) /
					(RadialMaximum - RadialMinimum), 0.0,
				static_cast<double>(RadialSamples - 1));
			const std::size_t vertical0 = std::min(
				static_cast<std::size_t>(verticalCoordinate), VerticalSamples - 2);
			const std::size_t radial0 = std::min(
				static_cast<std::size_t>(radialCoordinate), RadialSamples - 2);
			const double verticalFraction = verticalCoordinate -
				static_cast<double>(vertical0);
			const double radialFraction = radialCoordinate -
				static_cast<double>(radial0);
			auto at = [&](std::size_t v, std::size_t r) -> const MagneticFieldSample&
			{
				return samples[v * RadialSamples + r];
			};
			auto blend = [&](double MagneticFieldSample::*component)
			{
				const double lower = at(vertical0, radial0).*component +
					radialFraction * (at(vertical0, radial0 + 1).*component -
						at(vertical0, radial0).*component);
				const double upper = at(vertical0 + 1, radial0).*component +
					radialFraction * (at(vertical0 + 1, radial0 + 1).*component -
						at(vertical0 + 1, radial0).*component);
				return lower + verticalFraction * (upper - lower);
			};
			return {blend(&MagneticFieldSample::vertical),
				blend(&MagneticFieldSample::radial)};
		}
	};

	static std::array<double, 2> MagneticFluxGradient(double vertical,
		double horizontal, double gap)
	{
		// The table is two-dimensional in alignment and radial clearance. Horizontal
		// tine motion changes that clearance, so the second component follows by the
		// chain rule. Construction occurs once per process; every voice then pays
		// only a bilinear lookup at the 4x pickup rate.
		constexpr double HorizontalFieldScale = 0.62;
		const double radialDistance = std::sqrt(gap * gap +
			HorizontalFieldScale * horizontal * horizontal + 0.020);
		static const MagneticFieldTable Table;
		const auto field = Table.Interpolate(vertical, radialDistance);
		return {field.vertical, field.radial * HorizontalFieldScale * horizontal /
			std::max(1.0e-9, radialDistance)};
	}

	static double MagneticTrajectoryRms(double verticalOffset,
		double horizontalOffset, double gap, double nominalExcursion)
	{
		// The normalization trajectory is dimensionless in time, so its angular
		// frequency cancels from the ratio. A small quadrature component represents
		// the measured two-polarization orbit without claiming a per-key ellipse.
		// Average low, nominal and high excursions geometrically: one arbitrary
		// MIDI velocity then cannot dominate the pickup calibration.
		constexpr std::array<double, 3> ExcursionScales{0.45, 1.0, 1.8};
		constexpr std::array<double, 12> PhaseCosines{
			0.9659258263, 0.7071067812, 0.2588190451, -0.2588190451,
			-0.7071067812, -0.9659258263, -0.9659258263, -0.7071067812,
			-0.2588190451, 0.2588190451, 0.7071067812, 0.9659258263};
		constexpr std::array<double, 12> PhaseSines{
			0.2588190451, 0.7071067812, 0.9659258263, 0.9659258263,
			0.7071067812, 0.2588190451, -0.2588190451, -0.7071067812,
			-0.9659258263, -0.9659258263, -0.7071067812, -0.2588190451};
		double logarithmicRms = 0.0;
		for (double excursionScale : ExcursionScales)
		{
			const double verticalAmplitude = nominalExcursion * excursionScale;
			const double horizontalAmplitude = 0.12 * verticalAmplitude;
			double energy = 0.0;
			for (std::size_t sample = 0; sample < PhaseCosines.size(); ++sample)
			{
				const double cosine = PhaseCosines[sample];
				const double sine = PhaseSines[sample];
				const auto gradient = MagneticFluxGradient(
					verticalOffset + verticalAmplitude * cosine,
					horizontalOffset + horizontalAmplitude * sine, gap);
				const double emf = gradient[0] * (-verticalAmplitude * sine) +
					gradient[1] * (horizontalAmplitude * cosine);
				energy += emf * emf;
			}
			const double rms = std::sqrt(energy /
				static_cast<double>(PhaseCosines.size()));
			// The explicit branch both handles a hypothetical NaN and makes the
			// strictly-positive logarithm precondition visible to static analysis.
			const double positiveRms = rms > 1.0e-12 ? rms : 1.0e-12;
			logarithmicRms += std::log(positiveRms);
		}
		return std::exp(logarithmicRms /
			static_cast<double>(ExcursionScales.size()));
	}

	static double ReferenceGradientMagnitude()
	{
		const auto gradient = MagneticFluxGradient(ReferencePickupVertical,
			ReferencePickupHorizontal, ReferencePickupGap);
		return std::sqrt(gradient[0] * gradient[0] +
			gradient[1] * gradient[1]);
	}

	double WhiteNoise()
	{
		noiseState_ ^= noiseState_ << 13;
		noiseState_ ^= noiseState_ >> 17;
		noiseState_ ^= noiseState_ << 5;
		return 2.0 * (static_cast<double>(noiseState_) /
			static_cast<double>(UINT32_MAX)) - 1.0;
	}

