	void RefreshTimbreCoefficients(const ElectricPianoControls& controls,
		bool controlTick)
	{
		// The nonlinear trajectory normalizer below is deliberately control-rate.
		// All inputs have already been smoothed at audio rate, and 1 kHz coefficient
		// refresh is far above knob/CV bandwidth; evaluating dozens of field points
		// for every active voice on every host sample would make this calibration a
		// significant steady CPU regression while a control is moving.
		if (!timbreDirty_ && !controlTick)
			return;
		const double body = Clamp01(controls.body);
		const double bell = Clamp01(controls.bell);
		const double coupling = Clamp01(controls.coupling);
		const double proximity = Clamp01(controls.proximity);
		const double tone = Clamp01(controls.tone);
		const double mechanics = Clamp01(controls.mechanics);
		if (!timbreDirty_ && body == cachedBodyWeight_ &&
			bell == cachedBellWeight_ &&
			coupling == cachedCouplingWeight_ &&
			proximity == cachedProximity_ &&
			tone == cachedTone_ &&
			mechanics == cachedMechanics_)
			return;

		// Both coupled normal coordinates are normalized to the same unit tine
		// displacement. The pickup must therefore observe them with identical
		// gain; changing their weights independently would break the underlying
		// two-body system and turn Coupling back into a hidden level control.
		modeOutputWeight_[0] = 0.58 + 0.52 * body;
		modeOutputWeight_[1] = modeOutputWeight_[0];
		modeOutputWeight_[2] = 0.06 + 0.12 * body;
		// Bell is a deliberate balance control around the calibrated physical
		// residue, not a correction for the default model. The upper keyboard has
		// progressively fewer attack coordinates below Nyquist, so normalize panel
		// sensitivity from a 24 dB bass/middle span toward 36 dB at C6 rather than
		// leaving the same knob effectively dead there. The midpoint remains unity.
		constexpr double F3Position = (53.0 - 28.0) / 72.0;
		constexpr double MiddleCPosition = (60.0 - 28.0) / 72.0;
		constexpr double C6Position = (84.0 - 28.0) / 72.0;
		const double upperBellLinear = std::clamp(
			(keyPosition_ - F3Position) / (C6Position - F3Position), 0.0, 1.0);
		const double upperBellSmooth = upperBellLinear * upperBellLinear *
			(3.0 - 2.0 * upperBellLinear);
		const double bellOctaveSpan = 4.0 + 2.0 * upperBellSmooth;
		const double bellGain = std::exp2(bellOctaveSpan * (bell - 0.52));
		for (std::size_t index = ElectricPianoAttackModeBegin;
			index < ElectricPianoAttackModeEnd; ++index)
			modeOutputWeight_[index] =
				ElectricPianoMechanicalTrim::AttackModeOutput *
				bellGain;

		// Imperfect transverse coupling gives the tine tip a shallow elliptical
		// orbit. The near-fundamental third mode carries most of the horizontal
		// motion; attack modes alternate orientation as observed on a real fork.
		modeHorizontalWeight_[0] = 0.010;
		modeHorizontalWeight_[1] = -0.012;
		modeHorizontalWeight_[2] = 0.42 + 0.20 * body;
		modeOutputWeight_[ElectricPianoToneBarSubMode] =
			modeOutputWeight_[1] * ToneBarSubModeTineParticipation(keyPosition_) *
			ElectricPianoToneBarSubModeCouplingScale(coupling,
				keyPosition_);
		for (std::size_t index = ElectricPianoAttackModeBegin;
			index < ElectricPianoAttackModeEnd; ++index)
			modeHorizontalWeight_[index] = (index % 2 == 0 ? -1.0 : 1.0) *
			0.097 * bellGain;
		modeHorizontalWeight_[ElectricPianoToneBarSubMode] = 0.0;

		// The service adjustment called "timbre" is the tine's vertical alignment
		// to the pickup pole. Moving toward the pole centre suppresses the linear
		// fundamental relative to curvature-generated harmonics. PROXIMITY controls
		// the independent front-to-back distance and therefore magnetic curvature.
		// Each real key is voiced by moving its pickup. Trajectory normalization
		// below makes this a dynamic-timbre calibration, not a register-dependent
		// gain contour.
		// The service range is expressed directly in millimetres: ordinary
		// voicing is 1/16--1/8 inch, while post-1972 middle/upper pickups can be
		// brought as close as 0.020 inch. The panel midpoint is the keyed neutral
		// setup below; each half is logarithmic so the close end retains useful
		// resolution.
		constexpr double MinimumServiceGapMillimetres = 0.020 * 25.4;
		constexpr double DefaultServiceGapMillimetres = 1.6;
		constexpr double MaximumServiceGapMillimetres = 0.125 * 25.4;
		// A serviced harp is individually voiced. Long bass tines need the ordinary
		// 1/16-inch clearance, while the service manual explicitly permits 0.020 inch
		// in the middle/upper range to retain dynamic response from their much smaller
		// travel. Preserve the successful bass and bark calibration exactly through
		// middle C, then graduate the neutral per-key screw setting toward 0.52 mm at
		// the top, still just above the documented 0.020-inch service limit. The
		// former 0.60 mm endpoint left the highest keys with little pickup-curvature
		// sparkle even after the mechanical topology was corrected. This remains a
		// physical gap change, not a register EQ or output crossfade.
		constexpr double TrebleNeutralGapMillimetres = 0.52;
		const double trebleGapPosition = std::clamp(
			(keyPosition_ - MiddleCPosition) / (1.0 - MiddleCPosition), 0.0, 1.0);
		const double trebleGapSmooth = trebleGapPosition * trebleGapPosition *
			(3.0 - 2.0 * trebleGapPosition);
		const double neutralKeyGap = DefaultServiceGapMillimetres * std::pow(
			TrebleNeutralGapMillimetres / DefaultServiceGapMillimetres,
			trebleGapSmooth);
		// Only the closest endpoint is additionally zoned: bass keys retain the
		// ordinary 1/16-inch service limit.
		const double closeGapTransition = std::clamp(keyPosition_ /
			((52.0 - 28.0) / 72.0), 0.0, 1.0);
		const double closeGapSmooth = closeGapTransition * closeGapTransition *
			(3.0 - 2.0 * closeGapTransition);
		const double minimumKeyGap = 0.0625 * 25.4 + closeGapSmooth *
			(MinimumServiceGapMillimetres - 0.0625 * 25.4);
		pickupGap_ = proximity < 0.48 ?
			MaximumServiceGapMillimetres * std::pow(
				neutralKeyGap / MaximumServiceGapMillimetres,
				proximity / 0.48) :
			neutralKeyGap * std::pow(
				minimumKeyGap / neutralKeyGap,
				(proximity - 0.48) / 0.52);
		// At minimum Tone the tine rests close to the primary edge's maximum
		// gradient. Raising Tone moves it onto one flank, trading odd-harmonic
		// symmetry for the stronger even/sideband mixture used when voicing a
		// Rhodes for bite.
		// In the treble, most inharmonic attack coordinates lie above pickup
		// bandwidth and raw modal Bell gain cannot create a useful control. Real
		// Rhodes voicing also changes bell character by moving the pickup vertically
		// relative to the tine. Add that second physical route only above middle C,
		// where it smoothly takes over from the mechanical residues. The signed
		// offset is exactly zero at Bell's calibrated 0.52 default. Trajectory and
		// alignment normalization below remove the associated broadband level shift,
		// leaving curvature-generated harmonic balance rather than another VCA.
		const double pickupBellLinear = std::clamp(
			(keyPosition_ - MiddleCPosition) /
				(C6Position - MiddleCPosition), 0.0, 1.0);
		const double pickupBellBlend = pickupBellLinear * pickupBellLinear *
			(3.0 - 2.0 * pickupBellLinear);
		constexpr double PickupBellAlignmentSpanMillimetres = 0.09;
		const double neutralPickupVerticalOffset = 0.34 + 0.22 * tone * tone +
			0.020 * keyPosition_;
		const double minimumPickupVerticalOffset = 0.34 + 0.020 * keyPosition_;
		const double maximumPickupVerticalOffset = 0.56 + 0.020 * keyPosition_;
		pickupVerticalOffset_ = std::clamp(neutralPickupVerticalOffset +
			PickupBellAlignmentSpanMillimetres * pickupBellBlend * (bell - 0.52),
			minimumPickupVerticalOffset, maximumPickupVerticalOffset);
		pickupHorizontalOffset_ = 0.10 + 0.035 * keyPosition_;
		const auto alignmentGradient = MagneticFluxGradient(
			neutralPickupVerticalOffset, pickupHorizontalOffset_,
			ReferencePickupGap);
		const double alignmentMagnitude = std::sqrt(
			alignmentGradient[0] * alignmentGradient[0] +
			alignmentGradient[1] * alignmentGradient[1]);
		// Small-signal normalization alone cannot calibrate a nonlinear pickup: it
		// made the same panel movement cut bass keys while boosting upper keys by
		// roughly 6 dB. Normalize the RMS EMF of a representative keyed trajectory
		// instead. The nominal excursion follows the existing, documented bass
		// excursion trim; sampling three amplitudes makes the correction robust to
		// velocity without imposing an envelope on the actual signal. Curvature and
		// harmonic balance are retained because this is one scalar per pickup setup,
		// never an instantaneous waveshaper inversion.
		const double nominalExcursionMillimetres =
			(0.08 + 0.32 * std::pow(1.0 - keyPosition_, 1.5)) *
			pickupExcursionScale_;
		const double referenceTrajectory = MagneticTrajectoryRms(
			neutralPickupVerticalOffset, pickupHorizontalOffset_,
			ReferencePickupGap, nominalExcursionMillimetres);
		const double currentTrajectory = MagneticTrajectoryRms(
			pickupVerticalOffset_, pickupHorizontalOffset_, pickupGap_,
			nominalExcursionMillimetres);
		// Compare Bell's shifted observation against the same key's neutral pickup
		// placement, rather than allowing the shifted placement to establish a new
		// gain reference. This jointly normalizes Bell and Proximity across the
		// representative trajectories and prevents extreme Tone/Proximity settings
		// from turning their interaction into several decibels of hidden gain.
		const double trajectoryGain = referenceTrajectory /
			std::max(1.0e-9, currentTrajectory);
		// Preserve a restrained service-like sensitivity change: the close end may
		// become at most about 1 dB louder after trajectory normalization, instead
		// of acting as another output-level control.
		// The three-amplitude normalizer samples a generic ellipse; real keyed
		// trajectories retain a small register-dependent residual at close gaps.
		// This measured correction raises the under-normalized middle register
		// without reviving the former upper-key level jump.
		const double intentionalSensitivity = std::pow(
			neutralKeyGap / pickupGap_, 0.08);
		pickupAlignmentGain_ = std::clamp(ReferenceGradientMagnitude() /
			std::max(1.0e-6, alignmentMagnitude) * trajectoryGain *
			intentionalSensitivity, 0.04, 7.0);
		pickupVoltageScale_ = 0.130 *
			ElectricPianoMechanicalTrim::PickupOutput * keyPickupSensitivity_ *
			pickupAlignmentGain_ * inverseReferenceGradient_;
		// Preserve a wide sound-design range while making the calibrated default
		// a quiet mechanical contribution rather than a parallel percussion layer.
		mechanicsLevel_ = 0.060 * mechanics * (0.30 + 0.70 * mechanics);
		cachedBodyWeight_ = body;
		cachedBellWeight_ = bell;
		cachedCouplingWeight_ = coupling;
		cachedProximity_ = proximity;
		cachedTone_ = tone;
		cachedMechanics_ = mechanics;
		timbreDirty_ = false;
	}

