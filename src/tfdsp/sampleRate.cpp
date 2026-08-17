#include "sampleRate.hpp"


namespace tfdsp {

std::unique_ptr<X2Resampler_Order5> CreateX2Resampler_Butterworth5()
{
	Eigen::Array<double, 1, 1> directCoeffs;
	Eigen::Array<double, 1, 1> delayedCoeffs;

	directCoeffs << 1.0 / (5 + 2 * std::sqrt(5));
	delayedCoeffs << 5 - 2 * std::sqrt(5);

	return std::make_unique<X2Resampler_Order5>(directCoeffs, delayedCoeffs);
}
std::unique_ptr<X2Resampler_Order7> CreateX2Resampler_Chebychev7()
{
	Eigen::Array<double, 2, 1> directCoeffs;
	Eigen::Array<double, 1, 1> delayedCoeffs;
	directCoeffs << 0.081430023176616115, 0.70977080010248506;
	delayedCoeffs << 0.31565984021666094;
	return std::make_unique<X2Resampler_Order7>(directCoeffs, delayedCoeffs);
}
std::unique_ptr<X2Resampler_Order9> CreateX2Resampler_Chebychev9()
{
	Eigen::Array<double, 2, 1> directCoeffs;
	Eigen::Array<double, 2, 1> delayedCoeffs;
	directCoeffs << 0.079866426236357438, 0.54532365107113168;
	delayedCoeffs << 0.28382934487410966, 0.83441189148073658;
	return std::make_unique<X2Resampler_Order9>(directCoeffs, delayedCoeffs);
}
std::unique_ptr<DummyResampler> CreateDummyResampler()
{
	return std::make_unique<DummyResampler>();
}
std::unique_ptr<X4Resampler_Order7> CreateX4Resampler_Cheby7()
{
	return std::make_unique<X4Resampler_Order7>(CreateX2Resampler_Chebychev7);
}

std::unique_ptr<X16Resampler_Order7> CreateX16Resampler_Cheby7()
{
	return std::make_unique<X16Resampler_Order7>(
		CreateX2Resampler_Chebychev7);
}

}
