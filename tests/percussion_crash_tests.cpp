#include "percussion_test_support.hpp"

#include "tfdsp/percussion/crash_cymbal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using percussion_test::Check;

namespace {

std::vector<float> Render(const float strength, const float location,
                          const float hardness, const std::uint32_t seed,
                          const float seconds = 2.f,
                          const float sampleRate = 48000.f) {
  using namespace tfdsp::percussion;
  CrashCymbal cymbal;
  cymbal.Prepare(sampleRate, DefaultCrashCymbalParameters(sampleRate));
  cymbal.Trigger({strength, location, hardness, seed});
  std::vector<float> output(
      static_cast<std::size_t>(seconds * sampleRate));
  for (float &sample : output)
    sample = cymbal.Process();
  return output;
}

double Energy(const std::vector<float> &audio) {
  double energy = 0.0;
  for (const float sample : audio)
    energy += static_cast<double>(sample) * sample;
  return energy;
}

double Difference(const std::vector<float> &first,
                  const std::vector<float> &second) {
  double energy = 0.0;
  for (std::size_t sample = 0; sample < first.size(); ++sample) {
    const double delta = first[sample] - second[sample];
    energy += delta * delta;
  }
  return energy;
}

void TestDeterministicAndResponsive() {
  const auto first = Render(.8f, 1.f, .65f, 1234);
  const auto repeated = Render(.8f, 1.f, .65f, 1234);
  const auto differentSeed = Render(.8f, 1.f, .65f, 1235);
  Check(first == repeated, "crash rendering is deterministic for one seed");
  Check(Difference(first, differentSeed) > 1.e-7,
        "crash contact variation responds to its seed");

  const auto bell = Render(.8f, 0.f, .65f, 1234);
  const auto bow = Render(.8f, .5f, .65f, 1234);
  Check(Difference(bell, bow) > .01 * Energy(bow),
        "crash location changes the body projection audibly");

  const auto soft = Render(.8f, 1.f, .15f, 1234);
  const auto hard = Render(.8f, 1.f, .95f, 1234);
  Check(Difference(soft, hard) > .01 * Energy(hard),
        "crash hardness changes the contact audibly");
}

void TestVelocityEnergy() {
  const std::array<float, 4> strengths{.2f, .45f, .7f, 1.f};
  double previous = 0.0;
  for (const float strength : strengths) {
    const double current = Energy(Render(strength, 1.f, .65f, 81));
    Check(current > previous,
          "crash energy increases across the velocity sweep");
    previous = current;
  }
}

void TestMuteIsPassive() {
  using namespace tfdsp::percussion;
  constexpr float sampleRate = 48000.f;
  CrashCymbal natural;
  CrashCymbal muted;
  const auto parameters = DefaultCrashCymbalParameters(sampleRate);
  natural.Prepare(sampleRate, parameters);
  muted.Prepare(sampleRate, parameters);
  natural.Trigger({.9f, 1.f, .65f, 9});
  muted.Trigger({.9f, 1.f, .65f, 9});
  double naturalTail = 0.0;
  double mutedTail = 0.0;
  for (std::size_t sample = 0; sample < 2 * 48000; ++sample) {
    if (sample == 4800)
      muted.SetMute(1.f);
    const float first = natural.Process();
    const float second = muted.Process();
    if (sample > 12000) {
      naturalTail += static_cast<double>(first) * first;
      mutedTail += static_cast<double>(second) * second;
    }
  }
  Check(mutedTail < .25 * naturalTail,
        "crash mute removes stored energy instead of sustaining it");

  CrashCymbal silent;
  silent.Prepare(sampleRate, parameters);
  silent.SetMute(1.f);
  double silentEnergy = 0.0;
  for (std::size_t sample = 0; sample < 48000; ++sample) {
    const float output = silent.Process();
    silentEnergy += static_cast<double>(output) * output;
  }
  Check(silentEnergy == 0.0, "changing crash mute does not inject energy");
}

void TestFiniteAtSupportedRates() {
  for (const float sampleRate : {44100.f, 48000.f, 96000.f, 192000.f}) {
    const auto audio = Render(1.f, .73f, 1.f, 0xffffffffu, .5f, sampleRate);
    Check(std::all_of(audio.begin(), audio.end(), [](const float sample) {
      return std::isfinite(sample) && std::abs(sample) < 100.f;
    }), "crash remains finite and bounded across sample rates");
  }
}

} // namespace

int main() {
  TestDeterministicAndResponsive();
  TestVelocityEnergy();
  TestMuteIsPassive();
  TestFiniteAtSupportedRates();
  if (percussion_test::failures == 0)
    std::cout << "All percussion crash tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
