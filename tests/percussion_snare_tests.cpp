#include "percussion_test_support.hpp"

#include "tfdsp/percussion/snare_drum.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

using percussion_test::Check;

std::vector<float> Render(const std::uint32_t seed) {
  tfdsp::percussion::SnareDrum drum;
  drum.Prepare(48000.f, tfdsp::percussion::DefaultSnareDrumParameters());
  drum.Trigger({.8f, .42f, .68f, 1.f, .2f, seed});
  std::vector<float> result(48000);
  for (auto &sample : result) sample = drum.Process();
  return result;
}

std::vector<float> RenderStrength(const float strength) {
  tfdsp::percussion::SnareDrum drum;
  drum.Prepare(48000.f, tfdsp::percussion::DefaultSnareDrumParameters());
  drum.Trigger({strength, .3f, .65f, 1.f, .2f, 1575});
  std::vector<float> result(5760);
  for (auto &sample : result)
    sample = drum.Process();
  return result;
}

double Energy(const std::vector<float> &audio) {
  double result = 0.;
  for (const float sample : audio) result += sample * sample;
  return result;
}

void TestDeterminismAndSources() {
  const auto first = Render(17);
  Check(first == Render(17), "snare rendering is deterministic");
  Check(first != Render(18), "snare seed varies its stochastic sources");
  Check(std::all_of(first.begin(), first.end(), [](const float sample) {
    return std::isfinite(sample);
  }) && Energy(first) > 1.e-5, "snare output is finite and audible");

  tfdsp::percussion::SnareDrum drum;
  drum.Prepare(48000.f, tfdsp::percussion::DefaultSnareDrumParameters());
  drum.Trigger({.8f, .4f, .6f, 1.f, .2f, 29});
  double bodyEnergy = 0.;
  double wireEnergy = 0.;
  for (std::size_t sample = 0; sample < 24000; ++sample) {
    const auto frame = drum.ProcessFrame();
    bodyEnergy += frame.body * frame.body;
    wireEnergy += frame.wires * frame.wires;
  }
  Check(bodyEnergy > 1.e-5 && wireEnergy > 1.e-5,
        "one hit excites both membrane and body-driven wires");
}

void TestVelocityResponse() {
  const double soft = Energy(RenderStrength(.32f));
  const double medium = Energy(RenderStrength(.65f));
  const double hard = Energy(RenderStrength(.92f));
  Check(soft < .2 * medium && hard > .8 * medium,
        "snare velocity layers retain useful soft-to-medium dynamics");
}

void TestRoutingAndPreparedState() {
  using namespace tfdsp::percussion;
  auto parameters = DefaultSnareDrumParameters();
  parameters.routing.Set(static_cast<std::size_t>(
      SnareDrumRoute::BodyToWires), 0.f);
  SnareDrum dry;
  dry.Prepare(48000.f, parameters);
  dry.Trigger({.8f, .5f, .6f, 1.f, .2f, 7});
  double wireEnergy = 0.;
  for (std::size_t sample = 0; sample < 12000; ++sample) {
    const auto frame = dry.ProcessFrame();
    wireEnergy += frame.wires * frame.wires;
  }
  Check(wireEnergy < 1.e-20, "disabled body-to-wire route is silent");

  const auto prepared = PrepareSnareDrumParameters(
      48000.f, DefaultSnareDrumParameters());
  SnareDrum first;
  SnareDrum second;
  first.Prepare(prepared);
  second.Prepare(prepared);
  const MembraneDrumHit hit{.75f, .25f, .7f, 1.f, .15f, 41};
  first.Trigger(hit);
  second.Trigger(hit);
  bool identical = true;
  for (std::size_t sample = 0; sample < 24000; ++sample)
    identical = identical && first.Process() == second.Process();
  Check(identical, "prepared snare state is sample-identical");
}

void TestRatesAndRepeatedHits() {
  using namespace tfdsp::percussion;
  for (const float sampleRate : {44100.f, 48000.f, 96000.f, 192000.f}) {
    SnareDrum drum;
    drum.Prepare(sampleRate, DefaultSnareDrumParameters());
    float peak = 0.f;
    for (std::size_t sample = 0;
         sample < static_cast<std::size_t>(sampleRate); ++sample) {
      if (sample % static_cast<std::size_t>(.125f * sampleRate) == 0)
        drum.Trigger({.7f, .45f, .65f, 1.f, .2f,
                      static_cast<std::uint32_t>(sample + 1)});
      peak = std::max(peak, std::abs(drum.Process()));
    }
    Check(std::isfinite(peak) && peak > .01f && peak < 4.f,
          "repeated snare hits remain finite across sample rates");
  }
}

} // namespace

int main() {
  TestDeterminismAndSources();
  TestVelocityResponse();
  TestRoutingAndPreparedState();
  TestRatesAndRepeatedHits();
  return percussion_test::failures == 0 ? 0 : 1;
}
