#pragma once

#include <algorithm>
#include <cmath>

namespace tfdsp
{

// Retriggerable RC-shaped envelope following the ARP 2600's 4020 ADSR and
// board-4 AR behavior. Segment times are specified to -60 dB from the target,
// which gives finite, predictable panel values while retaining exponential RC
// curves. Output is normalized; the Rack wrapper exposes the hardware's 10 V.
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
				_stage = Stage::Attack;
		}
		else if (_stage == Stage::Hold)
		{
			_stage = Stage::Decay;
		}
	}

	Mode GetMode() const { return _mode; }

	void Reset()
	{
		_value = 0.0;
		_gateHigh = false;
		_triggerHigh = false;
		_stage = Stage::Idle;
	}

	double Step(double gateVolts, double triggerVolts, double attackSeconds,
		double decaySeconds, double sustain, double releaseSeconds)
	{
		gateVolts = std::isfinite(gateVolts) ? gateVolts : 0.0;
		triggerVolts = std::isfinite(triggerVolts) ? triggerVolts : 0.0;
		attackSeconds = SanitizeTime(attackSeconds);
		decaySeconds = SanitizeTime(decaySeconds);
		releaseSeconds = SanitizeTime(releaseSeconds);
		sustain = std::isfinite(sustain) ? std::clamp(sustain, 0.0, 1.0) : 0.0;

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

		// The original 4020 requires gate and trigger, but derives a trigger
		// from an externally supplied gate edge. A separate trigger retriggers
		// attack only while gate is present. The AR circuit ignores trigger.
		if (gateRising || (_mode == Mode::Adsr && _gateHigh && triggerRising))
			_stage = Stage::Attack;
		if (gateFalling)
			_stage = Stage::Release;

		switch (_stage)
		{
		case Stage::Idle:
			_value = 0.0;
			break;
		case Stage::Attack:
			Approach(1.0, attackSeconds);
			if (_value >= 1.0 - CompletionError)
			{
				_value = 1.0;
				_stage = _mode == Mode::Ar ? Stage::Hold : Stage::Decay;
			}
			break;
		case Stage::Decay:
			Approach(sustain, decaySeconds);
			if (std::abs(_value - sustain) <= CompletionError)
			{
				_value = sustain;
				_stage = Stage::Sustain;
			}
			break;
		case Stage::Sustain:
			_value = sustain;
			break;
		case Stage::Hold:
			_value = 1.0;
			break;
		case Stage::Release:
			Approach(0.0, releaseSeconds);
			if (_value <= CompletionError)
			{
				_value = 0.0;
				_stage = Stage::Idle;
			}
			break;
		}

		if (!_gateHigh && _stage != Stage::Idle && _stage != Stage::Release)
			_stage = Stage::Release;
		return _value;
	}

	double Value() const { return _value; }
	bool GateHigh() const { return _gateHigh; }
	Stage GetStage() const { return _stage; }

	static constexpr double MinimumAttackSeconds = 0.0014;
	static constexpr double MaximumAttackSeconds = 1.5;
	static constexpr double MinimumDecaySeconds = 0.0064;
	static constexpr double MaximumDecaySeconds = 6.0;
	static constexpr double MinimumReleaseSeconds = 0.00052;
	static constexpr double MaximumReleaseSeconds = 6.0;

private:
	static constexpr double GateHighVolts = 1.0;
	static constexpr double GateLowVolts = 0.1;
	static constexpr double TriggerHighVolts = 1.0;
	static constexpr double TriggerLowVolts = 0.1;
	static constexpr double CompletionError = 1.0e-3;
	static constexpr double TimeToTau = 1.0 / 6.907755278982137;

	double _sampleRate{48000.0};
	double _value{};
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

	void Approach(double target, double durationSeconds)
	{
		const double tau = std::max(durationSeconds * TimeToTau, 1.0e-9);
		const double coefficient = -std::expm1(-1.0 / (_sampleRate * tau));
		_value += coefficient * (target - _value);
	}
};

} // namespace tfdsp
