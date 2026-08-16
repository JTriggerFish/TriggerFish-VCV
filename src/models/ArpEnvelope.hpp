#pragma once

#include <algorithm>
#include <cmath>

namespace tfdsp
{

// Retriggerable envelope based on the ARP 2600's 4020 ADSR and board-4 AR
// circuits. At the default curve, the attack follows a capacitor charging
// toward 15 V and crossing the 10 V peak threshold. Decay and release
// use the approximately three-time-constant interval reported for the 4020.
// The curve control varies those normalized exponentials while preserving the
// selected segment duration and continuous output.
class ArpEnvelope
{
public:
	enum class Mode
	{
		Adsr,
		Ar,
	};
	enum class Stage
	{
		Idle,
		Attack,
		Decay,
		Sustain,
		Hold,
		Release,
	};

	void SetSampleRate(double sampleRate)
	{
		if (std::isfinite(sampleRate) && sampleRate > 0.0)
			_sampleRate = sampleRate;
	}

	void SetMode(Mode mode)
	{
		if (mode == _mode)
			return;
		_mode = mode;
		if (!_gateHigh)
			return;
		if (_mode == Mode::Ar)
		{
			if (_stage != Stage::Attack && _stage != Stage::Hold)
				BeginAttack(_lastCurve);
		}
		else if (_stage == Stage::Hold)
		{
			BeginStage(Stage::Decay);
		}
	}

	Mode GetMode() const { return _mode; }

	void Reset()
	{
		_value = 0.0;
		_phase = 0.0;
		_gateHigh = false;
		_triggerHigh = false;
		_stage = Stage::Idle;
	}

	double Step(double gateVolts, double triggerVolts, double attackSeconds,
		double decaySeconds, double sustain, double releaseSeconds,
		double curve = 0.0, bool autoGateTrigger = true)
	{
		gateVolts = std::isfinite(gateVolts) ? gateVolts : 0.0;
		triggerVolts = std::isfinite(triggerVolts) ? triggerVolts : 0.0;
		attackSeconds = SanitizeTime(attackSeconds);
		decaySeconds = SanitizeTime(decaySeconds);
		releaseSeconds = SanitizeTime(releaseSeconds);
		sustain = std::isfinite(sustain) ? std::clamp(sustain, 0.0, 1.0) : 0.0;
		curve = std::isfinite(curve) ? std::clamp(curve, -1.0, 1.0) : 0.0;
		_lastCurve = curve;

		const bool gateRising = !_gateHigh && gateVolts >= GateHighVolts;
		const bool gateFalling = _gateHigh && gateVolts <= GateLowVolts;
		const bool triggerRising = !_triggerHigh &&
			triggerVolts >= TriggerHighVolts;
		if (gateRising)
			_gateHigh = true;
		else if (gateFalling)
			_gateHigh = false;
		if (!_triggerHigh && triggerVolts >= TriggerHighVolts)
			_triggerHigh = true;
		else if (_triggerHigh && triggerVolts <= TriggerLowVolts)
			_triggerHigh = false;

		// A gate edge provides the practical gate-only trigger used by the 2600
		// keyboard interface. A patched trigger can retrigger ADSR while its gate
		// remains high. The separate AR circuit ignores trigger.
		if (_mode == Mode::Ar)
		{
			if (gateRising)
				BeginAttack(curve);
		}
		else if (_gateHigh && triggerRising)
		{
			BeginAttack(curve);
		}
		else if (gateRising)
		{
			if (autoGateTrigger)
				BeginAttack(curve);
			else
				BeginStage(Stage::Decay);
		}
		if (gateFalling)
			BeginStage(Stage::Release);

		switch (_stage)
		{
		case Stage::Idle:
			_value = 0.0;
			break;
		case Stage::Attack:
			AdvanceSegment(1.0, attackSeconds,
				_mode == Mode::Ar ? FallingCurve(curve) : AttackCurve(curve));
			if (_phase >= 1.0)
			{
				_value = 1.0;
				if (_mode == Mode::Ar)
					BeginStage(Stage::Hold);
				else
					BeginStage(Stage::Decay);
			}
			break;
		case Stage::Decay:
			AdvanceSegment(sustain, decaySeconds, FallingCurve(curve));
			if (_phase >= 1.0)
			{
				_value = sustain;
				BeginStage(Stage::Sustain);
			}
			break;
		case Stage::Sustain:
			_value = sustain;
			break;
		case Stage::Hold:
			_value = 1.0;
			break;
		case Stage::Release:
			AdvanceSegment(0.0, releaseSeconds, FallingCurve(curve));
			if (_phase >= 1.0)
			{
				_value = 0.0;
				BeginStage(Stage::Idle);
			}
			break;
		}

		if (!_gateHigh && _stage != Stage::Idle && _stage != Stage::Release)
			BeginStage(Stage::Release);
		return _value;
	}

	double Value() const { return _value; }
	bool GateHigh() const { return _gateHigh; }
	Stage GetStage() const { return _stage; }

	static constexpr double MinimumAttackSeconds = 0.0014;
	static constexpr double MaximumAttackSeconds = 5.0;
	static constexpr double MinimumDecaySeconds = 0.0064;
	static constexpr double MaximumDecaySeconds = 6.0;
	static constexpr double MinimumReleaseSeconds = 0.00052;
	static constexpr double MaximumReleaseSeconds = 6.0;
	static constexpr double HardwareAttackTarget = 1.5;

	static double AttackCurve(double curve)
	{
		const double hardwareMagnitude =
			-std::log(1.0 - 1.0 / HardwareAttackTarget);
		return -CurveMagnitude(hardwareMagnitude, curve, 0.1,
			6.907755278982137);
	}

	static double FallingCurve(double curve)
	{
		return -CurveMagnitude(2.995732273553991, curve, 0.25, 8.0);
	}

	static double NormalizedCurve(double phase, double coefficient)
	{
		phase = std::clamp(phase, 0.0, 1.0);
		if (std::abs(coefficient) < 1.0e-8)
			return phase;
		return std::expm1(coefficient * phase) / std::expm1(coefficient);
	}

private:
	static constexpr double GateHighVolts = 1.0;
	static constexpr double GateLowVolts = 0.1;
	static constexpr double TriggerHighVolts = 1.0;
	static constexpr double TriggerLowVolts = 0.1;

	double _sampleRate{48000.0};
	double _value{};
	double _phase{};
	double _lastCurve{};
	bool _gateHigh{};
	bool _triggerHigh{};
	Mode _mode{Mode::Adsr};
	Stage _stage{Stage::Idle};

	static double SanitizeTime(double seconds)
	{
		if (!std::isfinite(seconds))
			return 0.1;
		return std::clamp(seconds, 1.0e-5, 60.0);
	}

	static double CurveMagnitude(double hardware, double curve,
		double minimum, double maximum)
	{
		curve = std::clamp(curve, -1.0, 1.0);
		if (curve < 0.0)
			return hardware * std::pow(hardware / minimum, curve);
		return hardware * std::pow(maximum / hardware, curve);
	}

	static double InverseNormalizedCurve(double value, double coefficient)
	{
		value = std::clamp(value, 0.0, 1.0);
		if (std::abs(coefficient) < 1.0e-8)
			return value;
		return std::log1p(value * std::expm1(coefficient)) / coefficient;
	}

	void BeginStage(Stage stage)
	{
		_stage = stage;
		_phase = 0.0;
	}

	void BeginAttack(double curve)
	{
		_stage = Stage::Attack;
		const double shape = _mode == Mode::Ar ? FallingCurve(curve) :
			AttackCurve(curve);
		_phase = InverseNormalizedCurve(_value, shape);
	}

	void AdvanceSegment(double target, double durationSeconds,
		double curveCoefficient)
	{
		const double oldPhase = _phase;
		double nextPhase = std::min(1.0,
			oldPhase + 1.0 / (_sampleRate * durationSeconds));
		if (nextPhase >= 1.0 - 1.0e-12)
			nextPhase = 1.0;
		const double oldCurve = NormalizedCurve(oldPhase, curveCoefficient);
		const double nextCurve = NormalizedCurve(nextPhase, curveCoefficient);
		const double remaining = 1.0 - oldCurve;
		const double fraction = remaining > 1.0e-12 ?
			std::clamp((nextCurve - oldCurve) / remaining, 0.0, 1.0) : 1.0;
		_value += fraction * (target - _value);
		_phase = nextPhase;
	}
};

} // namespace tfdsp
