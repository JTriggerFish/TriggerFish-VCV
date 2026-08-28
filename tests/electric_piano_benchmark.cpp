#include "models/ElectricPiano.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace
{
	constexpr double SampleRate = 48000.0;
	constexpr double TwoPi = 6.2831853071795864769;
	constexpr int RenderSamples = 48000;
	constexpr int WarmupSamples = 12000;

	struct Measurement
	{
		double milliseconds{};
		double realtimeLoad{};
		double checksum{};
	};

	enum class ModulationMode
	{
		Fixed,
		LegacyPitch,
		Exponential,
		Linear,
		Phase
	};

	template <std::size_t VoiceCount>
	Measurement MeasureVoices(ModulationMode mode,
		bool includeAmplifier = false)
	{
		std::array<tfdsp::ElectricPianoVoice, VoiceCount> voices;
		tfdsp::ElectricPianoAmplifier amplifier;
		amplifier.SetSampleRate(SampleRate);
		std::array<double, VoiceCount> pitches{};
		tfdsp::ElectricPianoControls controls;
		for (std::size_t voice = 0; voice < VoiceCount; ++voice)
		{
			voices[voice].SetSampleRate(SampleRate);
			voices[voice].SetNoiseSeed(0x6d2b79f5u +
				static_cast<std::uint32_t>(voice) * 0x9e3779b9u);
			pitches[voice] = -2.0 + 4.0 * static_cast<double>(voice) /
				static_cast<double>(VoiceCount);
		}

		std::vector<double> modulation(RenderSamples);
		for (int sample = 0; sample < RenderSamples; ++sample)
			modulation[static_cast<std::size_t>(sample)] =
				mode != ModulationMode::Fixed ?
				0.01 * std::sin(TwoPi * 5.0 * sample / SampleRate) : 0.0;

		for (int sample = 0; sample < WarmupSamples; ++sample)
		{
			double pickup = 0.0;
			for (std::size_t voice = 0; voice < VoiceCount; ++voice)
				pickup += voices[voice].Step(pitches[voice], 10.0, 0.85, false,
					controls);
			if (includeAmplifier)
				amplifier.Step(5.0 * pickup, controls);
		}

		double checksum = 0.0;
		const auto begin = std::chrono::steady_clock::now();
		for (int sample = 0; sample < RenderSamples; ++sample)
		{
			const double offset = modulation[static_cast<std::size_t>(sample)];
			double pickup = 0.0;
			for (std::size_t voice = 0; voice < VoiceCount; ++voice)
			{
				tfdsp::ElectricPianoModulation voiceModulation;
				if (mode == ModulationMode::Exponential)
					voiceModulation.exponentialPitch = offset;
				else if (mode == ModulationMode::Linear)
					voiceModulation.linearFrequencyRatio = 10.0 * offset;
				else if (mode == ModulationMode::Phase)
					voiceModulation.phaseRadians = 50.0 * offset;
				const double pitch = pitches[voice] +
					(mode == ModulationMode::LegacyPitch ? offset : 0.0);
				pickup += voices[voice].Step(pitch, 10.0, 0.85, false,
					controls, voiceModulation);
			}
			if (includeAmplifier)
			{
				const auto output = amplifier.Step(5.0 * pickup, controls);
				checksum += output[0] + output[1];
			}
			else
				checksum += pickup;
		}
		const auto end = std::chrono::steady_clock::now();
		const double milliseconds = std::chrono::duration<double, std::milli>(
			end - begin).count();
		return {milliseconds, milliseconds / 1000.0, checksum};
	}

	template <typename Render>
	Measurement BestOf(Render&& render)
	{
		Measurement best{std::numeric_limits<double>::infinity(), 0.0, 0.0};
		for (int repeat = 0; repeat < 3; ++repeat)
		{
			const auto measurement = render();
			if (measurement.milliseconds < best.milliseconds)
				best = measurement;
		}
		return best;
	}

	void Print(const char* label, const Measurement& measurement)
	{
		std::cout << std::left << std::setw(27) << label << std::right <<
			std::fixed << std::setprecision(2) << std::setw(9) <<
			measurement.milliseconds << " ms  " << std::setw(7) <<
			100.0 * measurement.realtimeLoad << "% realtime  checksum " <<
			std::setprecision(9) << measurement.checksum << '\n';
	}
}

int main()
{
	Print("one voice, fixed pitch", BestOf([]
	{
		return MeasureVoices<1>(ModulationMode::Fixed);
	}));
	Print("16 voices, fixed pitch", BestOf([]
	{
		return MeasureVoices<16>(ModulationMode::Fixed);
	}));
	Print("16 voices + shared amp", BestOf([]
	{
		return MeasureVoices<16>(ModulationMode::Fixed, true);
	}));
	Print("16 voices, pitch modulated", BestOf([]
	{
		return MeasureVoices<16>(ModulationMode::LegacyPitch);
	}));
	Print("16 voices, EXP FM", BestOf([]
	{
		return MeasureVoices<16>(ModulationMode::Exponential);
	}));
	Print("16 voices, linear TZ FM", BestOf([]
	{
		return MeasureVoices<16>(ModulationMode::Linear);
	}));
	Print("16 voices, phase modulated", BestOf([]
	{
		return MeasureVoices<16>(ModulationMode::Phase);
	}));
}
