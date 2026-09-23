#include "models/ElectricPiano.hpp"

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <stdexcept>

namespace {
std::atomic<std::size_t> trackedAllocations{};
thread_local bool trackAllocations{};
int failures{};
}

void *operator new(const std::size_t size) {
  if (trackAllocations)
    trackedAllocations.fetch_add(1, std::memory_order_relaxed);
  if (void *memory = std::malloc(size))
    return memory;
  throw std::bad_alloc{};
}

void *operator new[](const std::size_t size) { return ::operator new(size); }
void operator delete(void *memory) noexcept { std::free(memory); }
void operator delete[](void *memory) noexcept { ::operator delete(memory); }
void operator delete(void *memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void *memory, std::size_t) noexcept {
  ::operator delete(memory);
}

namespace {

void Check(const bool condition, const char *message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << '\n';
  }
}

template <typename Function>
void CheckRejectsNonFiniteRate(Function function, const char *message) {
  bool rejected = false;
  try {
    function();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  Check(rejected, message);
}

void TestSampleRateContracts() {
  tfdsp::ElectricPianoVoice voice;
  tfdsp::ElectricPianoAmplifier amplifier;
  CheckRejectsNonFiniteRate(
      [&] { voice.SetSampleRate(std::numeric_limits<double>::quiet_NaN()); },
      "electric-piano voice rejects a NaN sample rate");
  CheckRejectsNonFiniteRate(
      [&] { amplifier.SetSampleRate(std::numeric_limits<double>::infinity()); },
      "electric-piano amplifier rejects an infinite sample rate");
}

void TestAudioPathDoesNotAllocate() {
  tfdsp::ElectricPianoVoice voice;
  tfdsp::ElectricPianoAmplifier amplifier;
  tfdsp::ElectricPianoControls controls;
  voice.SetSampleRate(48000.0);
  amplifier.SetSampleRate(48000.0);
  for (std::size_t sample = 0; sample < 256; ++sample) {
    const double gate = sample < 128 ? 10.0 : 0.0;
    (void)amplifier.Step(voice.Step(0.0, gate, .75, false, controls), controls);
  }
  trackedAllocations.store(0, std::memory_order_relaxed);
  trackAllocations = true;
  bool finite = true;
  for (std::size_t sample = 0; sample < 4096; ++sample) {
    const double gate = sample < 2048 ? 10.0 : 0.0;
    const double voiceOutput = voice.Step(0.0, gate, .75, false, controls);
    const auto output = amplifier.Step(voiceOutput, controls);
    finite = finite && std::isfinite(output[0]) && std::isfinite(output[1]);
  }
  trackAllocations = false;
  Check(finite, "electric-piano audio path remains finite");
  Check(trackedAllocations.load(std::memory_order_relaxed) == 0,
        "electric-piano voice and amplifier allocate nothing while processing");
}

} // namespace

int main() {
  TestSampleRateContracts();
  TestAudioPathDoesNotAllocate();
  if (failures == 0)
    std::cout << "All electric-piano realtime tests passed\n";
  return failures == 0 ? 0 : 1;
}
