#include "tfdsp/late_reverb.hpp"
#include "tfdsp/room_reverb.hpp"

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <new>

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

void TestLateReverbDoesNotAllocate() {
  tfdsp::LateReverb reverb;
  tfdsp::LateReverbControls controls;
  reverb.SetSampleRate(48000.0);
  trackedAllocations.store(0, std::memory_order_relaxed);
  trackAllocations = true;
  bool finite = true;
  for (std::size_t sample = 0; sample < 8192; ++sample) {
    const float phase = static_cast<float>(sample) / 8191.f;
    controls.diffusion = phase;
    controls.modulation = 1.f - phase;
    controls.decay = .2f + .8f * phase;
    const auto output = reverb.Process(sample == 0 ? 1.f : 0.f, controls);
    finite = finite && std::isfinite(output[0]) && std::isfinite(output[1]);
  }
  trackAllocations = false;
  Check(finite, "late reverb remains finite under simultaneous automation");
  Check(trackedAllocations.load(std::memory_order_relaxed) == 0,
        "late reverb allocates nothing while processing");
}

void TestRoomReverbDoesNotAllocate() {
  tfdsp::RoomReverb reverb;
  tfdsp::RoomReverbControls controls;
  tfdsp::RoomReverb::SourcePositions positions{};
  reverb.SetSampleRate(48000.0);
  positions[0] = {.25, .4, .5};
  trackedAllocations.store(0, std::memory_order_relaxed);
  trackAllocations = true;
  bool finite = true;
  for (std::size_t sample = 0; sample < 8192; ++sample) {
    const float phase = static_cast<float>(sample) / 8191.f;
    controls.preDelay = phase;
    controls.diffusion = 1.f - phase;
    positions[0][0] = .2 + .6 * phase;
    tfdsp::RoomReverb::InputFrame input{};
    input[0] = sample == 0 ? 1.f : 0.f;
    const auto output = reverb.Process(input, positions, 1, controls);
    finite = finite && std::isfinite(output.direct[0]) &&
             std::isfinite(output.direct[1]) && std::isfinite(output.wet[0]) &&
             std::isfinite(output.wet[1]);
  }
  trackAllocations = false;
  Check(finite, "room reverb remains finite under simultaneous automation");
  Check(trackedAllocations.load(std::memory_order_relaxed) == 0,
        "room reverb allocates nothing while processing");
}

} // namespace

int main() {
  TestLateReverbDoesNotAllocate();
  TestRoomReverbDoesNotAllocate();
  if (failures == 0)
    std::cout << "All reverb realtime tests passed\n";
  return failures == 0 ? 0 : 1;
}
