#include "sampleRate.hpp"


namespace dsp {

X2Resampler_Order5* CreateX2Resampler_Butterworth5()
{
    return new PolyphaseIIR_X2Resampler<1,1> (
		{{1.0 / (5 + 2 * std::sqrt(5))}},
		{{5 - 2 * std::sqrt(5)}}
    );
}
X2Resampler_Order7* CreateX2Resampler_Chebychev7()
{
   return new PolyphaseIIR_X2Resampler<2,1> (
		{{0.081430023176616115, 0.70977080010248506 }},
		{{0.31565984021666094 }}
    );
}
X2Resampler_Order9* CreateX2Resampler_Chebychev9()
{
    return new PolyphaseIIR_X2Resampler<2,2> (
		{{ 0.079866426236357438, 0.54532365107113168}},
		{{0.28382934487410966, 0.83441189148073658 }}
    );
}
DummyResampler* CreateDummyResampler()
{
    return new DummyResampler();
}
X4Resampler_Order7* CreateX4Resampler_Cheby7()
{
	return new X4Resampler_Order7(CreateX2Resampler_Chebychev7);
}

}