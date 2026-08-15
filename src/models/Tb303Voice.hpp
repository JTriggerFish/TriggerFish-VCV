#pragma once

#include <algorithm>
#include <cmath>

#include "AnalogOutputStage.hpp"
#include "OtaVca.hpp"

namespace tfdsp
{

class Tb303AccentSweep
{
public:
	enum class Mode
	{
		Off,
		Fast,
		Normal,
		Slow,
	};

	void SetSampleRate(double sampleRate)
	{
		if (std::isfinite(sampleRate) && sampleRate > 0.0)
			_sampleRate = sampleRate;
	}

	void Reset()
	{
		_capacitor = 0.0;
	}

	double Step(double source, double resonance)
	{
		if (!std::isfinite(source))
			source = 0.0;
		if (!std::isfinite(resonance))
			resonance = 0.0;
		source = std::clamp(source, 0.0, 1.0);
		resonance = std::clamp(resonance, 0.0, 1.0);

		if (_mode == Mode::Off)
		{
			const double coefficient = -std::expm1(-1.0 /
				(_sampleRate * ReleaseSeconds));
			_capacitor += coefficient * (0.0 - _capacitor);
			return 0.0;
		}

		// At the anti-clockwise end, the 47k source resistor and 100k
		// summing resistor give 100/(47+100) of the accented MEG directly.
		constexpr double DirectGain = 100.0 / 147.0;

		// At the clockwise end, C13 (1 uF) charges through 147k while
		// loaded by the 100k summing resistor: gain 100/(147+100) and
		// tau=(147k||100k)*1uF=59.5 ms. When the diode opens, C13 drains
		// through the 100k load, giving the documented 100 ms release.
		constexpr double CapacitorGain = 100.0 / 247.0;
		constexpr double AttackSeconds = (147.0 * 100.0 / 247.0) * 1.0e-3;
		double target = CapacitorGain * source;
		double attackSeconds = AttackSeconds;
		if (_mode == Mode::Fast)
		{
			// The exact Version-1.x switch component values are not public.
			// A leaky peak detector reproduces the documented behavior: a
			// strong first accent and progressively smaller repeated accents.
			target = source;
			attackSeconds = 0.010;
		}
		else if (_mode == Mode::Slow)
		{
			// The manual describes Slow as a much longer response capable of
			// twice Normal's sweep. Keep this explicitly behavioral rather
			// than inventing undocumented Devil Fish component values.
			target = 2.0 * CapacitorGain * source;
			attackSeconds = 4.0 * AttackSeconds;
		}
		const double releaseSeconds = _mode == Mode::Slow ?
			4.0 * ReleaseSeconds : ReleaseSeconds;
		const double time = target > _capacitor ? attackSeconds : releaseSeconds;
		const double coefficient = -std::expm1(-1.0 / (_sampleRate * time));
		_capacitor += coefficient * (target - _capacitor);
		if (std::abs(_capacitor) < 1.0e-15 && target == 0.0)
			_capacitor = 0.0;

		const double clockwise = _mode == Mode::Fast ?
			DirectGain * std::max(source - _capacitor, 0.0) : _capacitor;
		return (1.0 - resonance) * DirectGain * source +
			resonance * clockwise;
	}

	void SetMode(Mode mode)
	{
		if (_mode != mode)
		{
			_mode = mode;
			Reset();
		}
	}

	Mode GetMode() const { return _mode; }

	double CapacitorState() const { return _capacitor; }

private:
	static constexpr double ReleaseSeconds = 0.100;
	double _sampleRate{48000.0};
	double _capacitor{};
	Mode _mode{Mode::Normal};
};

class Tb303Articulation
{
public:
	enum class Mode
	{
		Stock,
		DevilFish,
	};

	struct Output
	{
		double mainEnvelope{};
		double filterAccent{};
		double volumeEnvelope{};
		double vcaAccent{};
	};

	void SetSampleRate(double sampleRate)
	{
		if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
			return;
		_sampleRate = sampleRate;
		_accentSweep.SetSampleRate(sampleRate);
		_vcaAccentCoefficient = -std::expm1(-1.0 /
			(_sampleRate * VcaAccentTimeSeconds));
	}

	void SetMode(Mode mode)
	{
		_mode = mode;
	}

	Mode GetMode() const { return _mode; }

	void SetAccentSweepMode(Tb303AccentSweep::Mode mode)
	{
		_accentSweep.SetMode(mode);
	}

	Tb303AccentSweep::Mode GetAccentSweepMode() const
	{
		return _accentSweep.GetMode();
	}

	void Reset()
	{
		_mainEnvelope = 0.0;
		_volumeEnvelope = 0.0;
		_vcaAccent = 0.0;
		_attackStart = 0.0;
		_releaseStart = 0.0;
		_stageTime = 0.0;
		_gateHigh = false;
		_volumeStage = VolumeStage::Idle;
		_accentSweep.Reset();
	}

	Output Step(double gateVolts, double accentVolts, double resonance,
		double normalDecaySeconds, double accentDecaySeconds,
		double vcaDecayControl)
	{
		if (!std::isfinite(gateVolts))
			gateVolts = 0.0;
		if (!std::isfinite(accentVolts))
			accentVolts = 0.0;
		const double accent = std::clamp(accentVolts / 10.0, 0.0, 1.0);
		resonance = std::isfinite(resonance) ?
			std::clamp(resonance, 0.0, 1.0) : 0.0;
		normalDecaySeconds = std::isfinite(normalDecaySeconds) ?
			std::clamp(normalDecaySeconds, 0.030, 3.0) : 0.5;
		accentDecaySeconds = std::isfinite(accentDecaySeconds) ?
			std::clamp(accentDecaySeconds, 0.030, 3.0) : 0.2;

		const bool rising = !_gateHigh && gateVolts >= 1.0;
		const bool falling = _gateHigh && gateVolts <= 0.1;
		if (rising)
		{
			_gateHigh = true;
			_mainEnvelope = 1.0;
			_attackStart = _volumeEnvelope;
			_stageTime = 0.0;
			_volumeStage = VolumeStage::Delay;
		}
		else if (falling)
		{
			_gateHigh = false;
			_releaseStart = _volumeEnvelope;
			_stageTime = 0.0;
			_volumeStage = _mode == Mode::Stock ?
				VolumeStage::ReleaseHold : VolumeStage::ReleaseExponential;
		}

		const double main = _mainEnvelope;
		const double accentSource = main * accent;
		const double filterAccent = _accentSweep.Step(accentSource, resonance);
		_vcaAccent += _vcaAccentCoefficient * (accentSource - _vcaAccent);

		const double decaySeconds = std::exp(
			(1.0 - accent) * std::log(normalDecaySeconds) +
			accent * std::log(accentDecaySeconds));
		_mainEnvelope *= std::exp(-1.0 / (_sampleRate * decaySeconds));
		if (_mainEnvelope < 1.0e-15)
			_mainEnvelope = 0.0;

		StepVolumeEnvelope(vcaDecayControl);
		return {main, filterAccent, _volumeEnvelope, _vcaAccent};
	}

	double AccentMemory() const { return _accentSweep.CapacitorState(); }

	static void MapVcaDecay(double control, double& decaySeconds,
		double& sustain)
	{
		control = std::isfinite(control) ? std::clamp(control, 0.0, 1.0) : 0.5;
		if (control <= 0.5)
		{
			constexpr double MinimumDecay = 0.016;
			constexpr double StockDecay = 3.5;
			decaySeconds = MinimumDecay * std::pow(StockDecay / MinimumDecay,
				control / 0.5);
			sustain = 0.0;
		}
		else
		{
			decaySeconds = 3.5;
			sustain = (control - 0.5) / 0.5;
		}
	}

private:
	enum class VolumeStage
	{
		Idle,
		Delay,
		Attack,
		Decay,
		ReleaseHold,
		ReleaseLinear,
		ReleaseExponential,
	};

	static constexpr double StockDelaySeconds = 0.004;
	static constexpr double DevilFishDelaySeconds = 0.0005;
	static constexpr double AttackSeconds = 0.003;
	static constexpr double StockReleaseHoldSeconds = 0.008;
	static constexpr double StockReleaseLinearSeconds = 0.008;
	// 8 ms to -60 dB, expressed as an exponential time constant.
	static constexpr double DevilFishReleaseTimeConstant = 0.001158118618408672;
	static constexpr double VcaAccentTimeSeconds = 47000.0 * 33.0e-9;

	double _sampleRate{48000.0};
	double _mainEnvelope{};
	double _volumeEnvelope{};
	double _vcaAccent{};
	double _vcaAccentCoefficient{-std::expm1(-1.0 /
		(48000.0 * VcaAccentTimeSeconds))};
	double _attackStart{};
	double _releaseStart{};
	double _stageTime{};
	bool _gateHigh{};
	Mode _mode{Mode::Stock};
	VolumeStage _volumeStage{VolumeStage::Idle};
	Tb303AccentSweep _accentSweep{};

	void StepVolumeEnvelope(double vcaDecayControl)
	{
		const double delta = 1.0 / _sampleRate;
		_stageTime += delta;
		double decaySeconds = 3.5;
		double sustain = 0.0;
		MapVcaDecay(vcaDecayControl, decaySeconds, sustain);

		switch (_volumeStage)
		{
		case VolumeStage::Idle:
			_volumeEnvelope = 0.0;
			break;
		case VolumeStage::Delay:
		{
			const double delay = _mode == Mode::Stock ?
				StockDelaySeconds : DevilFishDelaySeconds;
			if (_stageTime >= delay)
			{
				_stageTime -= delay;
				_volumeStage = VolumeStage::Attack;
			}
			break;
		}
		case VolumeStage::Attack:
		{
			const double progress = std::clamp(_stageTime / AttackSeconds, 0.0, 1.0);
			_volumeEnvelope = _attackStart + (1.0 - _attackStart) * progress;
			if (progress >= 1.0)
			{
				_stageTime = 0.0;
				_volumeStage = VolumeStage::Decay;
			}
			break;
		}
		case VolumeStage::Decay:
		{
			const double coefficient = -std::expm1(-delta / decaySeconds);
			_volumeEnvelope += coefficient * (sustain - _volumeEnvelope);
			break;
		}
		case VolumeStage::ReleaseHold:
			if (_stageTime >= StockReleaseHoldSeconds)
			{
				_stageTime -= StockReleaseHoldSeconds;
				_releaseStart = _volumeEnvelope;
				_volumeStage = VolumeStage::ReleaseLinear;
			}
			break;
		case VolumeStage::ReleaseLinear:
		{
			const double progress = std::clamp(
				_stageTime / StockReleaseLinearSeconds, 0.0, 1.0);
			_volumeEnvelope = _releaseStart * (1.0 - progress);
			if (progress >= 1.0)
			{
				_volumeEnvelope = 0.0;
				_volumeStage = VolumeStage::Idle;
			}
			break;
		}
		case VolumeStage::ReleaseExponential:
			_volumeEnvelope *= std::exp(-delta / DevilFishReleaseTimeConstant);
			if (_volumeEnvelope < 1.0e-6)
			{
				_volumeEnvelope = 0.0;
				_volumeStage = VolumeStage::Idle;
			}
			break;
		}
		_volumeEnvelope = std::clamp(_volumeEnvelope, 0.0, 1.0);
	}
};

class Tb303Vca
{
public:
	void SetSampleRate(double sampleRate)
	{
		if (!std::isfinite(sampleRate) || sampleRate <= 0.0)
			return;
		_sampleRate = sampleRate;
		_outputCouplingCoefficient = -std::expm1(-2.0 * Pi *
			OutputCouplingCutoffHz / _sampleRate);
	}

	void Reset()
	{
		_outputLowPass = 0.0;
		_lastControlCurrent = 0.0;
		_core.Reset();
	}

	double Step(double audioRackVolts, double baseControl,
		double accentControl)
	{
		if (!std::isfinite(audioRackVolts) || !std::isfinite(baseControl) ||
			!std::isfinite(accentControl))
		{
			Reset();
			return 0.0;
		}
		baseControl = std::clamp(baseControl, 0.0, 1.0);
		accentControl = std::clamp(accentControl, 0.0, 1.0);
		_lastControlCurrent = BaseControlCurrentAmps * baseControl +
			AccentControlCurrentAmps * accentControl;

		const double differential = RackToDifferentialScale * audioRackVolts;
		const double outputCurrent = _core.ProcessCurrent(differential,
			_lastControlCurrent);
		// The 220k pin-6 load converts OTA current to the hardware voltage. The
		// output mirror and buffer progressively lose compliance before an ideal
		// current source could reach the supply rail. Split this physical node
		// from the Rack calibration so its overload can be modeled explicitly.
		const double physicalVolts = PhysicalOutputLoadOhms * outputCurrent;
		const double limitedPhysicalVolts = SoftPhysicalCompliance(physicalVolts);
		const double rackVolts = RackOutputCalibration * limitedPhysicalVolts;

		// C38 is the BA662 buffer's 1 uF output coupling capacitor. The
		// following 50k volume control gives a 3.18 Hz corner, below the 303's
		// already modeled filter coupling losses.
		_outputLowPass += _outputCouplingCoefficient *
			(rackVolts - _outputLowPass);
		return rackVolts - _outputLowPass;
	}

	double LastControlCurrent() const { return _lastControlCurrent; }
	const OtaVcaCore& Core() const { return _core; }

	static constexpr double RackToDifferentialScale = 0.001414213562373095;
	static constexpr double BaseControlCurrentAmps = 20.0e-6;
	static constexpr double AccentControlCurrentAmps = 20.0e-6;
	static constexpr double PhysicalOutputLoadOhms = 220000.0;
	static constexpr double RackOutputCalibration = 9.818181818181818;

private:
	static constexpr double Pi = 3.14159265358979323846;
	static constexpr double OutputCouplingCutoffHz = 1.0 /
		(2.0 * Pi * 50000.0 * 1.0e-6);

	double _sampleRate{48000.0};
	double _outputCouplingCoefficient{-std::expm1(-2.0 * Pi *
		OutputCouplingCutoffHz / 48000.0)};
	double _outputLowPass{};
	double _lastControlCurrent{};
	OtaVcaCore _core{};

	static double SoftPhysicalCompliance(double voltage)
	{
		// Rack calibration maps these physical voltages to the same 8 V knee
		// and 11 V asymptote used by the modeled output buffer. These values are
		// a conservative headroom calibration; original BA662 output-compliance
		// curves are not available over the required operating range.
		constexpr double knee = AnalogOutputStage::KneeVolts /
			RackOutputCalibration;
		constexpr double rail = AnalogOutputStage::RailVolts /
			RackOutputCalibration;
		constexpr double headroom = rail - knee;
		const double magnitude = std::abs(voltage);
		if (magnitude <= knee)
			return voltage;
		return std::copysign(knee + headroom * std::tanh(
			(magnitude - knee) / headroom), voltage);
	}
};

} // namespace tfdsp
