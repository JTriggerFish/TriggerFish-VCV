#include "percussion_test_support.hpp"
#include "reference_bounded_modal_motion.hpp"
#include "tfdsp/percussion/bounded_modal_motion.hpp"
#include "tfdsp/percussion/stochastic_modal_field.hpp"
#include "tfdsp/fixed_radix2_fft.hpp"
#include <complex>
#include <vector>

using percussion_test::Check;
using percussion_test::CheckNear;
using namespace tfdsp::percussion;

void TestBoundAndSharing() {
  BoundedModalMotion<3> motion;
  motion.Prepare(48000.f, {3000.f, 3010.f, 6000.f}, {0,0,1}, 3, {.6f,80.f,1.f}, 7);
  std::array<double,3> integral{};
  std::vector<float> replay;
  bool distinct = false;
  for (int i=0;i<96000;++i) {
    motion.BeginSample();
    const float a=motion.NextAngle(0), b=motion.NextAngle(1), c=motion.NextAngle(2);
    Check(a==b,"shared packet phase preserves internal beat relationships");
    distinct |= a!=c;
    replay.push_back(a);
    integral[0]+=a; integral[1]+=b; integral[2]+=c;
    for (double x:integral) Check(std::abs(x)<=1.20001,"phase displacement does not random-walk");
  }
  Check(distinct,"different packets do not share a global LFO");
  motion.Reset();
  for(float a:replay) {
    motion.BeginSample();
    Check(a==motion.NextAngle(0),"motion reset is deterministic");
    motion.NextAngle(1); motion.NextAngle(2);
  }
}

void TestCarrierAndEnergy() {
  using Field=StochasticModalField<1>;
  for (float rate:{32000.f,44100.f,48000.f,96000.f}) {
    Field::Parameters mode{{{6000.f,3.f,1.f,1.f,0.f,0.f}}};
    StochasticModalFieldControls controls{};
    controls.motion={.6f,80.f,.5f};
    Field field, replay, blur;
    auto prepared=PrepareStochasticModalField(rate,mode,controls,500.f,5000.f);
    field.LoadPrepared(prepared); replay.LoadPrepared(prepared);
    mode[0].phaseBandwidthHz=20.f;
    blur.Prepare(rate,mode,{},500.f,5000.f);
    std::complex<double> carrier{}, blurred{};
    const int frames=static_cast<int>(rate);
    for(int i=0;i<frames;++i) {
      const float y=field.ProcessExcitedPair(i==0?1.f:0.f,0.f);
      Check(y==replay.ProcessExcitedPair(i==0?1.f:0.f,0.f),"prepared motion reproduces audio");
      const double expected=std::pow(.001,2.*i/(3.*rate));
      CheckNear(field.StoredEnergy(),expected,.004,"movement preserves declared damping");
      const auto phase=std::polar(1./std::sqrt(expected),-6.28318530718*6000*i/rate);
      carrier+=double(y)*phase;
      blurred+=double(blur.ProcessExcitedPair(i==0?1.f:0.f,0.f))*phase;
    }
    const double retained=2.*std::abs(carrier)/frames;
    Check(retained>.85,"bounded motion retains a strong coherent ridge");
    Check(retained>2.*std::abs(blurred)/frames+.2,"motion retains more carrier than random phase blur");
  }
}

void TestUpperSidebands() {
  // A complex carrier distinguishes sidebands crossing Nyquist from its real
  // signal's intentional negative-frequency image. Hann removes window edges.
  constexpr std::size_t size=32768;
  static tfdsp::FixedRadix2Fft<double,size> fft;
  static decltype(fft)::Spectrum spectrum;
  fft.Prepare();
  for(float rate:{32000.f,44100.f,48000.f,96000.f}) {
    BoundedModalMotion<1> motion;
    motion.Prepare(rate,{15000.f},{0},1,{3.f,200.f,.25f},17);
    double displacement=0;
    for(std::size_t i=0;i<size;++i) {
      motion.BeginSample(); displacement+=motion.NextAngle(0);
      const double window=.5-.5*std::cos(6.28318530718*i/(size-1));
      spectrum[i]=std::polar(window,6.28318530718*15000*i/rate+displacement);
    }
    fft.Transform(spectrum,false);
    double wanted=0,wrapped=0;
    for(std::size_t i=0;i<size;++i)
      (i<size/2?wanted:wrapped)+=std::norm(spectrum[i]);
    Check(wrapped/wanted<1.e-6,"15 kHz maximum movement has < -60 dB wrapped sideband energy");
  }
}

template <std::size_t Count> void TestPreparedRotationEquivalence(float rate) {
  std::array<float, Count> frequencies{};
  std::array<std::uint16_t, Count> packets{};
  for (std::size_t i=0; i<Count; ++i) {
    frequencies[i] = 1.f + .48f*rate*float(i+1)/float(Count+1);
    packets[i] = static_cast<std::uint16_t>(i/5);
  }
  BoundedModalMotion<Count> actual;
  ReferenceBoundedModalMotion<Count> expected;
  for (float sharing : {0.f, .5f, 1.f}) {
    actual.Prepare(rate,frequencies,packets,Count,{3.f,200.f,sharing},123);
    expected.Prepare(rate,frequencies,packets,Count,{3.f,200.f,sharing},123);
    for (int frame=0; frame<4000; ++frame) {
      actual.BeginSample(); expected.BeginSample();
      actual.PrepareRotations();
      for (std::size_t i=0; i<Count; ++i) {
        float ac=.7f, as=-.4f, ec=ac, es=as;
        actual.RotatePrepared(i,ac,as);
        expected.RotateCoefficients(i,ec,es);
        CheckNear(ac,ec,2.e-7,"prepared rotation matches scalar cosine");
        CheckNear(as,es,2.e-7,"prepared rotation matches scalar sine");
      }
    }
  }
}

int main() {
  TestBoundAndSharing();
  TestCarrierAndEnergy();
  TestUpperSidebands();
  TestPreparedRotationEquivalence<9>(1000.f);
  TestPreparedRotationEquivalence<17>(8000.f);
  TestPreparedRotationEquivalence<513>(48000.f);
  return percussion_test::failures?1:0;
}
