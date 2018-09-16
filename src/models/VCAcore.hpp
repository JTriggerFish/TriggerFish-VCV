#pragma once

#include <memory>
#include <cmath>
#include <functional>
#include <array>
#include <random>
#include "../dsp/filters.hpp"
#include "../dsp/noise.hpp"
#include "OTA1PoleIntegrator.hpp"
#include "Transistor1PoleIntegrator.hpp"

template<typename Float>
class TanhBlock
{
private:
	Float _x1{};
public:
	Float Process(Float x)
	{
		Float y = Tanh<Float>::Value(x, _x1);
		_x1 = x;
		return y;
	}
};

template<typename Oversampler, typename Model>
class VCACore
{
private:
	VCACore(const VCACore&) = delete;
	VCACore& operator=(const VCACore&) = delete;

	//Oversampling of audio and cv:--------------------------------------------------
	float _sampleRate{};
	static constexpr unsigned int ResamplingFactor{ Oversampler::ResamplingFactor };
	std::unique_ptr<Oversampler> _audioResampler;
	std::unique_ptr<Oversampler> _cvResampler;
	//-------------------------------------------------------------------------------

	//Models for audio and cv inputs:
	//[0] : audioPath
	//[1] : cvPath
	std::array<Model, 2> _models{};
	std::array<double, 2> _rolloffs{ {Model::DefaultRolloff,Model::DefaultRolloff} };
	std::array<double, 2> _g{};//Normalised and prewarped rolloffs

	dsp::PinkNoiseSource _noise{};
	double _noiseLevel{ 1.0e-10 };
	double _noiseStdDev{};

	double _cvScaling{3.0}; // For additional cv staturation. TODO trimmer for this ?
	double _powerSupplyVoltage{ 12.0 };
	TanhBlock<double> _outputStage{};


public:
	VCACore(std::function<Oversampler*()> resamplerCreator) :
		_audioResampler{ resamplerCreator() }, _cvResampler{ resamplerCreator() }
	{
	}
	void SetSampleRate(const float f0)
	{
		_sampleRate = f0 * ResamplingFactor;
		//g = w^~_c * T  = 2* tan( wc *T / 2 ) ( prewarping )
		//i.e g = 2 * tan( pi / 2 * f / (fo/2)) = tan(pi / 2 * fc)
		for(unsigned int i=0; i < 2; ++i)
		{
			_g[i] = _rolloffs[i] / (0.5*_sampleRate);
			_g[i] = 2 * std::tan(PI / 2 * _g[i]);
		}
		//Conserve the power spectral density independently of the sample rate
		_noiseStdDev = std::sqrt( _noiseLevel * _sampleRate / 2);
	}
	float Step(const float audio, const float cv, const float finalGain)
	{
		if (_sampleRate <= 0.f)
		{
			throw std::runtime_error("Sample rate invalid or not initialized");
		}
		const double noise = _noiseStdDev * _noise.Step();
		double input = noise + audio;

		std::array<double, ResamplingFactor> audioA = _audioResampler->Upsample(input);
		std::array<double, ResamplingFactor> cvA = _cvResampler->Upsample(double(_cvScaling*cv));

		Step(audioA, cvA, finalGain);

		return float(_audioResampler->Downsample(audioA));
	}
private:
	inline void Step(std::array<double, ResamplingFactor>& audio, const std::array<double, ResamplingFactor>& cv, const double finalGain)
	{
		for (unsigned int i = 0; i < ResamplingFactor; ++i)
		{
			//TODO: 
			//  kill signal when CV << noise floor ?

			std::array<double, 2> audioAndCv;
			audioAndCv[0] = audio[i];
			audioAndCv[1] = cv[i];

			Model::StepDual(_models, audioAndCv, _g);

			audio[i] = audioAndCv[0] * audioAndCv[1] / _cvScaling;
		}
		//Apply final gain and saturate to power supply voltage
		for (unsigned int i = 0; i < ResamplingFactor; ++i)
		{
			audio[i] = _powerSupplyVoltage * _outputStage.Process(finalGain * audio[i] / _powerSupplyVoltage);
		}
	}
};
template<typename T>
using VCA_OTACore = ::VCACore<T, OTA1PoleIntegrator>;
template<typename T>
using VCA_TransistorCore = ::VCACore<T, Transistor1PoleIntegrator>;
