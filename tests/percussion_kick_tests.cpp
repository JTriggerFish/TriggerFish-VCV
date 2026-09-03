#include "percussion_test_support.hpp"

#include "tfdsp/percussion/compact_kick.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace {

using percussion_test::Check;

std::vector<float> Render(const float sampleRate, const float strength,
                          const float hardness, const std::uint32_t seed,
                          const std::array<float, 3> routes = {1.f, 1.f, 1.f}) {
  auto parameters = tfdsp::percussion::DefaultCompactKickParameters();
  for (std::size_t index = 0; index < routes.size(); ++index)
    parameters.routing.Set(index, routes[index]);
  tfdsp::percussion::CompactKick kick;
  kick.Prepare(sampleRate, parameters);
  kick.Trigger({strength, hardness, seed});
  std::vector<float> result(static_cast<std::size_t>(sampleRate));
  for (float &sample : result) sample = kick.Process();
  return result;
}

double Energy(const std::vector<float> &audio) {
  double result = 0.0;
  for (const float sample : audio) result += sample * sample;
  return result;
}

double Difference(const std::vector<float> &left,
                  const std::vector<float> &right) {
  double result = 0.0;
  for (std::size_t index = 0; index < left.size(); ++index) {
    const double difference = left[index] - right[index];
    result += difference * difference;
  }
  return result;
}

double HighFrequencyEnergy(const std::vector<float> &audio,
                           const float sampleRate) {
  const double coefficient = 1.0 - std::exp(-6.283185307179586 * 2500.0 /
                                             sampleRate);
  double low = 0.0;
  double result = 0.0;
  for (const float sample : audio) {
    low += coefficient * (sample - low);
    const double high = sample - low;
    result += high * high;
  }
  return result;
}

void TestDeterminismAndRates() {
  for (const float sampleRate : {44100.f, 48000.f, 96000.f, 192000.f}) {
    const auto first = Render(sampleRate, .8f, .5f, 17);
    const auto repeated = Render(sampleRate, .8f, .5f, 17);
    const auto variation = Render(sampleRate, .8f, .5f, 18);
    Check(first == repeated, "compact kick is deterministic for a fixed seed");
    Check(Difference(first, variation) > 1.e-7,
          "compact kick seed changes the noisy realization");
    Check(std::all_of(first.begin(), first.end(), [](const float sample) {
      return std::isfinite(sample);
    }), "compact kick remains finite at every supported sample rate");
    Check(Energy(first) > 1.e-5, "compact kick remains audible");
  }
}

void TestVelocityAndHardness() {
  const auto silent = Render(48000.f, 0.f, 1.f, 7);
  const auto whisper = Render(48000.f, .001f, 1.f, 7);
  const auto soft = Render(48000.f, .2f, .2f, 7);
  const auto strong = Render(48000.f, 1.f, .2f, 7);
  const auto hard = Render(48000.f, 1.f, 1.f, 7);
  Check(Energy(strong) > 4.0 * Energy(soft),
        "kick strength materially raises event energy");
  Check(Energy(silent) < 1.e-20,
        "zero kick strength is an exactly silent event");
  Check(Energy(whisper) < .01 * Energy(soft),
        "kick energy approaches silence continuously with strength");
  const double strongHigh = HighFrequencyEnergy(strong, 48000.f);
  const double hardHigh = HighFrequencyEnergy(hard, 48000.f);
  if (!(hardHigh > 1.1 * strongHigh))
    std::cerr << "kick hard/soft high-frequency energy ratio "
              << hardHigh / strongHigh << '\n';
  Check(hardHigh > 1.1 * strongHigh,
        "kick hardness materially raises attack brightness");
}

void TestVoiceExhaustion() {
  auto controls = tfdsp::percussion::CompactKickControls{};
  controls.bodyDecaySeconds = 2.f;
  const auto parameters =
      tfdsp::percussion::DefaultCompactKickParameters(controls);
  tfdsp::percussion::CompactKick baseline;
  tfdsp::percussion::CompactKick stolen;
  baseline.Prepare(48000.f, parameters);
  stolen.Prepare(48000.f, parameters);
  for (std::size_t hit = 0; hit < tfdsp::percussion::CompactKick::VoiceCount;
       ++hit) {
    baseline.Trigger({.8f, .5f, static_cast<std::uint32_t>(100 + hit)});
    stolen.Trigger({.8f, .5f, static_cast<std::uint32_t>(100 + hit)});
    for (std::size_t sample = 0; sample < 300; ++sample) {
      baseline.Process();
      stolen.Process();
    }
  }
  stolen.Trigger({.8f, .5f, 999});
  const float expected = baseline.Process();
  const float replaced = stolen.Process();
  Check(std::abs(replaced - expected) < .05f,
        "voice exhaustion preserves the stolen waveform at the boundary");
}

void TestRoutes() {
  const auto complete = Render(48000.f, .8f, .5f, 9);
  const auto primary = Render(48000.f, .8f, .5f, 9, {1.f, 0.f, 0.f});
  const auto secondary = Render(48000.f, .8f, .5f, 9, {0.f, 1.f, 0.f});
  const auto click = Render(48000.f, .8f, .5f, 9, {0.f, 0.f, 1.f});
  const auto silent = Render(48000.f, .8f, .5f, 9, {0.f, 0.f, 0.f});
  Check(Energy(primary) > 1.e-5 && Energy(secondary) > 1.e-7 &&
            Energy(click) > 1.e-8,
        "each compact-kick source route is independently audible");
  Check(Difference(primary, secondary) > 1.e-5 &&
            Difference(complete, primary) > 1.e-6,
        "compact-kick sources make distinct contributions");
  Check(Energy(silent) < 1.e-20,
        "disabling every kick source route produces silence");
}

void TestOverlappingHits() {
  const auto parameters = tfdsp::percussion::DefaultCompactKickParameters();
  tfdsp::percussion::CompactKick single;
  tfdsp::percussion::CompactKick repeated;
  single.Prepare(48000.f, parameters);
  repeated.Prepare(48000.f, parameters);
  single.Trigger({.8f, .5f, 31});
  repeated.Trigger({.8f, .5f, 31});
  double singleEnergy = 0.0;
  double repeatedEnergy = 0.0;
  for (std::size_t sample = 0; sample < 12000; ++sample) {
    if (sample == 2400) repeated.Trigger({.8f, .5f, 32});
    const float one = single.Process();
    const float two = repeated.Process();
    if (sample >= 2400) {
      singleEnergy += one * one;
      repeatedEnergy += two * two;
    }
  }
  Check(repeatedEnergy > 1.3 * singleEnergy,
        "retriggering adds a voice without resetting the active kick tail");
}

} // namespace

int main() {
  TestDeterminismAndRates();
  TestVelocityAndHardness();
  TestRoutes();
  TestOverlappingHits();
  TestVoiceExhaustion();
  if (percussion_test::failures == 0)
    std::cout << "All compact kick tests passed\n";
  return percussion_test::failures == 0 ? 0 : 1;
}
