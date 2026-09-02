#include "crash_api.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

int failures{};

void Check(const bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

std::vector<float> Render(const std::uint32_t handle, const std::uint32_t seed,
                          const float location, const std::size_t blockSize) {
  constexpr std::size_t FrameCount = 8192;
  std::vector<float> result(FrameCount);
  Check(tf_crash_reset(handle) == 1, "reset accepts a live handle");
  Check(tf_crash_trigger(handle, .8f, location, .65f, seed) == 1,
        "trigger accepts a live handle");
  for (std::size_t first = 0; first < result.size(); first += blockSize) {
    const auto count = std::min(blockSize, result.size() - first);
    Check(tf_crash_process(handle, result.data() + first,
                           static_cast<std::uint32_t>(count)) == 1,
          "block processing succeeds");
  }
  return result;
}

double Difference(const std::vector<float> &first,
                  const std::vector<float> &second) {
  double result = 0.0;
  for (std::size_t index = 0; index < first.size(); ++index) {
    const double delta = first[index] - second[index];
    result += delta * delta;
  }
  return result;
}

} // namespace

int main() {
  Check(tf_crash_api_version() == 1, "API version is explicit");
  Check(tf_crash_macro_count() == 10, "the compact macro surface is versioned");
  Check(tf_crash_macro_name(0) != nullptr &&
            tf_crash_macro_unit(0) != nullptr,
        "macro labels and units are available through the C ABI");
  Check(tf_crash_macro_name(tf_crash_macro_count()) == nullptr,
        "out-of-range macro metadata is rejected");
  Check(tf_crash_create(0.f) == 0, "invalid sample rates are rejected");
  const auto handle = tf_crash_create(48000.f);
  Check(handle != 0, "a renderer session can be created");

  const auto first = Render(handle, 17, .8f, 256);
  const auto repeated = Render(handle, 17, .8f, 256);
  const auto whole = Render(handle, 17, .8f, first.size());
  const auto variation = Render(handle, 18, .8f, 256);
  Check(first == repeated, "equal events render deterministically");
  Check(first == whole, "rendering is independent of host block size");
  Check(Difference(first, variation) > 1.e-8,
        "the event seed changes stochastic contact");
  Check(std::all_of(first.begin(), first.end(), [](const float sample) {
    return std::isfinite(sample);
  }), "rendered samples remain finite");

  const auto controlHandle = tf_crash_create(48000.f);
  Check(controlHandle != 0, "a second live renderer can be created");
  std::vector<float> liveFirst(4096), liveSecond(4096), controlSecond(4096);
  Check(tf_crash_reset(handle) && tf_crash_reset(controlHandle),
        "live restrike test resets both renderers");
  Check(tf_crash_trigger(handle, .8f, .8f, .65f, 31) &&
            tf_crash_trigger(controlHandle, .8f, .8f, .65f, 31),
        "live restrike test starts equal states");
  Check(tf_crash_process(handle, liveFirst.data(), liveFirst.size()) &&
            tf_crash_process(controlHandle, liveFirst.data(), liveFirst.size()),
        "live restrike prefix renders");
  Check(tf_crash_trigger(handle, .7f, .8f, .65f, 32),
        "a running cymbal accepts another strike");
  Check(tf_crash_process(handle, liveSecond.data(), liveSecond.size()) &&
            tf_crash_process(controlHandle, controlSecond.data(),
                             controlSecond.size()),
        "restruck and free tails render");
  Check(Difference(liveSecond, controlSecond) > 1.e-6,
        "a restrike adds energy without resetting the existing tail");
  tf_crash_destroy(controlHandle);

  const float quieterLevel = tf_crash_macro_default(0) - 6.f;
  Check(tf_crash_macro_set(handle, 0, quieterLevel) == 1,
        "a valid macro edit is accepted");
  Check(tf_crash_macro_commit(handle) == 1,
        "macro edits prepare a new renderer state explicitly");
  const auto quieter = Render(handle, 17, .8f, 256);
  const double quietEnergy =
      Difference(quieter, std::vector<float>(quieter.size()));
  const double baselineEnergy =
      Difference(first, std::vector<float>(first.size()));
  Check(quietEnergy / baselineEnergy > .24 &&
            quietEnergy / baselineEnergy < .26,
        "model level applies the declared decibel scaling");
  Check(tf_crash_macro_set(handle, tf_crash_macro_count(), 0.f) == 0,
        "invalid macro indices are rejected");

  tf_crash_destroy(handle);
  Check(tf_crash_process(handle, nullptr, 0) == 0,
        "destroyed handles cannot be reused");
  return failures == 0 ? 0 : 1;
}
