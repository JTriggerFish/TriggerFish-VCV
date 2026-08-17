#pragma once

#include <Eigen/Dense>
#include <memory>
#include <array>
#include <cmath>
#include <functional>
#include <random>
#include "../tfdsp/filters.hpp"
#include "../tfdsp/noise.hpp"
#include "OTA1PoleIntegrator.hpp"
#include "tfdsp/rail.hpp"
#include "Transistor1PoleIntegrator.hpp"


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
	std::unique_ptr<Oversampler> _exponentialCvResampler;
	//-------------------------------------------------------------------------------

	//Models for audio and cv inputs:
	//[0] : audioPath
	//[1] : cvPath
	std::array<Model, 2> _models{};
	Eigen::Array<double, 2, 1> _rolloffs;
	Eigen::Array<double, 2, 1> _g;//Normalised and prewarped rolloffs

	tfdsp::PinkNoiseSource _noise{};
	double _noiseLevel{ 1.0e-10 };
	double _noiseStdDev{};

	double _cvScaling{3.0}; // For additional cv staturation. TODO trimmer for this ?
	double _powerSupplyVoltage{ 12.0 };
	TanhBlock<double, ResamplingFactor> _outputStage{};
	float _lastControl{};


public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	explicit VCACore(std::function<std::unique_ptr<Oversampler>()> resamplerCreator) :
		_audioResampler{ resamplerCreator() }, _cvResampler{ resamplerCreator() },
		_exponentialCvResampler{ resamplerCreator() }
	{
		_rolloffs << Model::DefaultRolloff, Model::DefaultRolloff;
	}
	void SetSampleRate(const float f0)
	{
		_sampleRate = f0 * ResamplingFactor;
		//g = w^~_c * T  = 2* tan( wc *T / 2 ) ( prewarping )
		//i.e g = 2 * tan( pi / 2 * f / (fo/2)) = tan(pi / 2 * fc)
		_g = 2.0 * Eigen::tan(tfdsp::PI / 2.0 * _rolloffs / (0.5 * _sampleRate));

		//Conserve the power spectral density independently of the sample rate
		_noiseStdDev = std::sqrt( _noiseLevel * _sampleRate / 2);
	}
	void Reset()
	{
		for (auto& model : _models)
			model.Reset();
		_audioResampler->Reset();
		_cvResampler->Reset();
		_exponentialCvResampler->Reset();
		_noise.Reset();
		_outputStage.Reset();
		_lastControl = 0.0f;
	}
	float Step(const float audio, const float cv, const float finalGain)
	{
		if (_sampleRate <= 0.f)
		{
			throw std::runtime_error("Sample rate invalid or not initialized");
		}
		if (!std::isfinite(audio) || !std::isfinite(cv) || !std::isfinite(finalGain))
		{
			Reset();
			return 0.0f;
		}
		const double noise = _noiseStdDev * _noise.Step();
		double input = noise + audio;

		auto audioA = _audioResampler->Upsample(input);
		auto cvA = _cvResampler->Upsample(double(_cvScaling*cv));

		Step(audioA, cvA, finalGain);
		const Eigen::Array<double, ResamplingFactor, 1> normalizedCv =
			cvA / _cvScaling;
		_lastControl = static_cast<float>(
			_cvResampler->Downsample(normalizedCv));

		const float output = static_cast<float>(
			tfdsp::RackOutputAdapter::ProcessPostDecimation(
				_audioResampler->Downsample(audioA)));
		if (!std::isfinite(output))
		{
			Reset();
			return 0.0f;
		}
		return output;
	}
	float LastControl() const { return _lastControl; }

	float StepControls(const float audio, const float linearCv,
		const float exponentialCv, const float exponentialBase,
		const float finalGain)
	{
		if (_sampleRate <= 0.f)
			throw std::runtime_error("Sample rate invalid or not initialized");
		if (!std::isfinite(audio) || !std::isfinite(linearCv) ||
			!std::isfinite(exponentialCv) || !std::isfinite(exponentialBase) ||
			!std::isfinite(finalGain))
		{
			Reset();
			return 0.0f;
		}

		const double noise = _noiseStdDev * _noise.Step();
		auto audioValues = _audioResampler->Upsample(noise + audio);
		auto linearValues = _cvResampler->Upsample(linearCv);
		const auto exponentialValues =
			_exponentialCvResampler->Upsample(exponentialCv);
		const double base = std::max<double>(exponentialBase, 1.0e-6);
		for (unsigned int i = 0; i < ResamplingFactor; ++i)
		{
			const double linear = std::clamp(linearValues(i), 0.0, 1.0);
			const double exponent = std::clamp(exponentialValues(i), 0.0, 1.0);
			const double shaped = std::abs(base - 1.0) < 1.0e-8 ? exponent :
				(std::pow(base, exponent) - 1.0) / (base - 1.0);
			linearValues(i) = _cvScaling *
				std::clamp(linear + shaped, 0.0, 1.0);
		}
		Step(audioValues, linearValues, finalGain);
		const Eigen::Array<double, ResamplingFactor, 1> normalizedCv =
			linearValues / _cvScaling;
		_lastControl = static_cast<float>(
			_cvResampler->Downsample(normalizedCv));

		const float output = static_cast<float>(
			tfdsp::RackOutputAdapter::ProcessPostDecimation(
				_audioResampler->Downsample(audioValues)));
		if (!std::isfinite(output))
		{
			Reset();
			return 0.0f;
		}
		return output;
	}
private:
	inline void Step(Eigen::Ref<Eigen::Array<double, ResamplingFactor, 1>> audio, const Eigen::Array<double, ResamplingFactor, 1>& cv, const double finalGain)
	{
		for (unsigned int i = 0; i < ResamplingFactor; ++i)
		{
			Eigen::Array<double, 2 ,1> audioAndCv;
			audioAndCv(0) = audio(i);
			audioAndCv(1) = cv(i);

			Model::StepDual(_models, audioAndCv, _g);

			audio(i) = audioAndCv(0) * audioAndCv(1) / _cvScaling;
		}
		//Apply final gain and saturate to power supply voltage
		audio = _powerSupplyVoltage *  _outputStage.Process(finalGain / _powerSupplyVoltage * audio);
		for (unsigned int i = 0; i < ResamplingFactor; ++i)
			audio(i) = tfdsp::RackOutputAdapter::ProcessOversampled(audio(i));
	}
};
template<typename T>
using VCA_OTACore = ::VCACore<T, OTA1PoleIntegrator>;
template<typename T>
using VCA_TransistorCore = ::VCACore<T, Transistor1PoleIntegrator>;
