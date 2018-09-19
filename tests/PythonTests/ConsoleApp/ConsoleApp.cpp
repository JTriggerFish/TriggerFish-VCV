// ConsoleApp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include "../../../src/dsp/filters.hpp"
#include "../../../src/dsp/sampleRate.hpp"
#include "../../../src/models/VCAcore.hpp"
#include "../../../src/models/VdpOscillator.hpp"

int main()
{

	auto vcaNoOversampling = std::unique_ptr<VCA_OTACore<dsp::DummyResampler>>{ new VCA_OTACore<dsp::DummyResampler>(dsp::CreateDummyResampler) };
	auto vca2xButterWorth = std::unique_ptr<VCA_OTACore<dsp::X2Resampler_Order5>>{ new VCA_OTACore<dsp::X2Resampler_Order5>(dsp::CreateX2Resampler_Butterworth5) };
	auto vca2xCheby7 = std::unique_ptr<VCA_OTACore<dsp::X2Resampler_Order7>>{ new VCA_OTACore<dsp::X2Resampler_Order7>(dsp::CreateX2Resampler_Chebychev7) };
	auto vca2xCheby9 = std::unique_ptr<VCA_OTACore<dsp::X2Resampler_Order9>>{ new VCA_OTACore<dsp::X2Resampler_Order9>(dsp::CreateX2Resampler_Chebychev9) };

	auto vca2x_Tr_Cheby7 = std::unique_ptr<VCA_TransistorCore<dsp::X2Resampler_Order7>>{ new VCA_TransistorCore<dsp::X2Resampler_Order7>(dsp::CreateX2Resampler_Chebychev7) };

	auto sampleRate = 48000.0;

	VdpOscillator<dsp::X2Resampler_Order7, 3> vdpo2x{ dsp::CreateX2Resampler_Chebychev7 };
	VdpOscillator<dsp::X4Resampler_Order7, 3> vdpo4x{ dsp::CreateX4Resampler_Cheby7 };

	vdpo2x.SetSampleRate(sampleRate);
	vdpo4x.SetSampleRate(sampleRate);

	auto T = 10.0;

	auto phase = 0.0;
	auto t = 0.0;

	auto f = 120.0;
	auto a = 5.0;

	vca2xCheby7->SetSampleRate(sampleRate);
	vca2x_Tr_Cheby7->SetSampleRate(sampleRate);

	std::vector<double> out( size_t(T*sampleRate) );

	size_t i = 0;
	auto mu = 3.0;
	auto w = 2 * PI * 1200;

	//auto w = 2 * PI * 3500;
	//auto wInc = 2 * PI * (4200-3500) / (T*sampleRate);
	while(t < T)
	{
		auto x = a * std::sin(2 * PI * phase);
		////auto y = vca2xCheby7->Step(x, 1.0);
		auto y = vca2x_Tr_Cheby7->Step(x, 1.0, 1.0);

		//auto y = vdpo2x.Step(0, mu, w);
		//auto y = vdpo4x.Step(0, mu, w);
		out[i] = y;
		//if(out[i] >= 12.0 || out[i] <= -12.0)
		//{
		//	auto breakHere = 0;
		//}

		i++;
		//w += wInc;
		t += 1.0 / sampleRate;
		phase += f / sampleRate;
		if (phase >= 1.0)
			phase -= 1.0;
	}
}
