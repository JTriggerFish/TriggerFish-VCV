#include "percussion_test_support.hpp"

#include "tfdsp/percussion/wire_rack.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

using percussion_test::Check;

double Energy(const std::vector<float> &audio) {
  double result = 0.;
  for (const float sample : audio) result += sample * sample;
  return result;
}

std::vector<float> Render(const float sampleRate, const bool fastMotion) {
  tfdsp::percussion::WireRack rack;
  rack.Prepare(sampleRate, {});
  std::vector<float> output(static_cast<std::size_t>(.5f * sampleRate));
  for (std::size_t sample = 0; sample < output.size(); ++sample) {
    const float time = static_cast<float>(sample) / sampleRate;
    const float drive = fastMotion
        ? (sample < static_cast<std::size_t>(.08f * sampleRate)
               ? .2f * std::sin(6.283185307179586f * 190.f * time) : 0.f)
        : std::min(.2f, .2f * time / .25f);
    output[sample] = rack.Process(drive);
  }
  return output;
}

void TestBodyDrivenResponse() {
  using namespace tfdsp::percussion;
  WireRack first;
  WireRack second;
  first.Prepare(48000.f, {});
  second.Prepare(48000.f, {});
  double energy = 0.;
  bool identical = true;
  for (std::size_t sample = 0; sample < 48000; ++sample) {
    const float drive = sample < 6000
        ? .12f * percussion_test::Sine(sample, 185.f, 48000.f) : 0.f;
    const float a = first.Process(drive);
    const float b = second.Process(drive);
    identical = identical && a == b;
    energy += a * a;
  }
  Check(identical, "wire rack is deterministic after reset");
  Check(energy > 1.e-4, "membrane motion excites the wire rack");
  Check(first.StoredEnergy() < 1.f,
        "wire modes decay after body motion stops");
}

void TestStaticAndSlowMotionArePassive() {
  using namespace tfdsp::percussion;
  WireRack rack;
  rack.Prepare(48000.f, {});
  double lateEnergy = 0.;
  for (std::size_t sample = 0; sample < 48000; ++sample) {
    const float output = rack.Process(.2f);
    if (sample > 40000) lateEnergy += output * output;
  }
  Check(lateEnergy < 1.e-8,
        "constant displacement does not sustain wire radiation");
  Check(Energy(Render(48000.f, true)) > 20. * Energy(Render(48000.f, false)),
        "fast membrane motion drives wires far more than a slow closure");
}

void TestOrdinaryMotionDoesNotReachSafetyCeiling() {
  using namespace tfdsp::percussion;
  const auto prepared = PrepareWireRackParameters(48000.f, {});
  WireRack rack;
  rack.Prepare(prepared);
  float maximumEnergy = 0.f;
  for (std::size_t sample = 0; sample < 24000; ++sample) {
    const float drive = sample < 6000
        ? .2f * percussion_test::Sine(sample, 190.f, 48000.f) : 0.f;
    rack.Process(drive);
    maximumEnergy = std::max(maximumEnergy, rack.StoredEnergy());
  }
  if (!(maximumEnergy < .25f * prepared.maximumModalEnergy))
    std::cerr << "wire ordinary-motion energy/capacity: " << maximumEnergy
              << '/' << prepared.maximumModalEnergy << '\n';
  Check(maximumEnergy < .25f * prepared.maximumModalEnergy,
        "ordinary wire motion stays clear of the safety ceiling");
}

void TestRatesAndBounds() {
  using namespace tfdsp::percussion;
  for (const float sampleRate : {44100.f, 48000.f, 96000.f, 192000.f}) {
    WireRack rack;
    const auto prepared = PrepareWireRackParameters(sampleRate, {});
    rack.Prepare(prepared);
    float peak = 0.f;
    float maximumEnergy = 0.f;
    for (std::size_t sample = 0;
         sample < static_cast<std::size_t>(sampleRate); ++sample) {
      const float drive = .8f * percussion_test::Sine(
          sample, 210.f, sampleRate);
      peak = std::max(peak, std::abs(rack.Process(drive)));
      maximumEnergy = std::max(maximumEnergy, rack.StoredEnergy());
    }
    if (!(std::isfinite(peak) && peak > .01f && peak < 2.f))
      std::cerr << "wire peak/energy at " << sampleRate << ": " << peak
                << '/' << maximumEnergy << '\n';
    Check(std::isfinite(peak) && peak > .01f && peak < 2.f,
          "wire rack remains finite and calibrated across sample rates");
    Check(maximumEnergy <= prepared.maximumModalEnergy + 1.e-4f,
          "wire rack respects its stored-energy capacity");
  }
}

} // namespace

int main() {
  TestBodyDrivenResponse();
  TestStaticAndSlowMotionArePassive();
  TestOrdinaryMotionDoesNotReachSafetyCeiling();
  TestRatesAndBounds();
  return percussion_test::failures == 0 ? 0 : 1;
}
