#include "percussion_api.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

int failures{};

void Check(const bool condition, const char *message) {
  if (condition)
    return;
  std::cerr << "FAIL: " << message << '\n';
  ++failures;
}

std::vector<float> Render(const std::uint32_t handle, const std::uint32_t seed,
                          const std::size_t blockSize) {
  std::vector<float> result(24000);
  Check(tf_percussion_reset(handle), "recipe reset succeeds");
  Check(tf_percussion_trigger(handle, .8f, .5f, .6f, 1.f, .2f, seed),
        "recipe trigger succeeds");
  for (std::size_t first = 0; first < result.size(); first += blockSize) {
    const auto count = std::min(blockSize, result.size() - first);
    Check(tf_percussion_process(handle, result.data() + first,
                                static_cast<std::uint32_t>(count)),
          "recipe block processing succeeds");
  }
  return result;
}

double Energy(const std::vector<float> &audio) {
  double result = 0.0;
  for (const float sample : audio)
    result += sample * sample;
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

std::vector<std::uint8_t> ExportPrepared(const std::uint32_t handle) {
  std::vector<std::uint8_t> result(tf_percussion_prepared_size(handle));
  Check(!result.empty() && tf_percussion_export_prepared(
                               handle, result.data(),
                               static_cast<std::uint32_t>(result.size())),
        "prepared recipe exports");
  return result;
}

void CheckPreparedRoundTrip(const std::uint32_t preparedHandle,
                            const std::uint32_t recipe,
                            const char *description) {
  const auto blob = ExportPrepared(preparedHandle);
  const auto restored = tf_percussion_create_unprepared(recipe, 48000.f);
  Check(restored != 0, "unprepared recipe session can be created");
  Check(tf_percussion_apply_prepared(restored, blob.data(),
                                     static_cast<std::uint32_t>(blob.size())),
        "prepared recipe installs");
  Check(Render(preparedHandle, 37, 128) == Render(restored, 37, 128),
        description);

  auto corrupt = blob;
  corrupt.front() ^= 0xff;
  Check(
      !tf_percussion_apply_prepared(restored, corrupt.data(),
                                    static_cast<std::uint32_t>(corrupt.size())),
      "corrupt prepared recipe is rejected");
  Check(!tf_percussion_apply_prepared(
            restored, blob.data(), static_cast<std::uint32_t>(blob.size() - 1)),
        "truncated prepared recipe is rejected");
  tf_percussion_destroy(restored);
}

std::uint32_t ParameterIndex(const std::uint32_t handle,
                             const std::string_view key) {
  for (std::uint32_t index = 0; index < tf_percussion_parameter_count(handle);
       ++index) {
    if (key == tf_percussion_parameter_key(handle, index))
      return index;
  }
  return UINT32_MAX;
}

} // namespace

int main() {
  Check(tf_percussion_api_version() == 1,
        "unreleased percussion API remains version one");
  Check(tf_percussion_recipe_count() == 4,
        "four compiled percussion recipes are registered");
  Check(
      std::string_view(tf_percussion_recipe_key(0)) == "metal.cymbal.v1" &&
          std::string_view(tf_percussion_recipe_key(1)) == "drum.kick-fm.v1" &&
          std::string_view(tf_percussion_recipe_key(2)) == "drum.membrane.v1" &&
          std::string_view(tf_percussion_recipe_key(3)) == "drum.snare.v1",
      "recipe keys are stable");
  Check(tf_percussion_recipe_key(4) == nullptr &&
            tf_percussion_create(4, 48000.f) == 0,
        "unknown recipes are rejected");

  const auto metallic = tf_percussion_create(0, 48000.f);
  Check(metallic != 0 && tf_percussion_parameter_count(metallic) == 126,
        "metallic recipe exposes only active unified-model parameters");
  Check(ParameterIndex(metallic, "body_tone_wash") == UINT32_MAX &&
            ParameterIndex(metallic, "dense_mode_density") == UINT32_MAX,
        "legacy no-op metallic controls are absent from the recipe API");
  CheckPreparedRoundTrip(metallic, 0,
                         "prepared metallic recipe is sample-identical");
  tf_percussion_destroy(metallic);

  const auto kick = tf_percussion_create(1, 48000.f);
  Check(kick != 0 && tf_percussion_recipe(kick) == 1,
        "compact kick session can be created");
  Check(tf_percussion_parameter_count(kick) == 15,
        "compact kick exposes its bounded control surface");
  Check(tf_percussion_route_count(kick) == 3,
        "compact kick exposes three source routes");
  const auto level = ParameterIndex(kick, "model_level_db");
  const auto pitch = ParameterIndex(kick, "fundamental_hz");
  Check(level < 15 && pitch < 15,
        "compact kick parameters have stable identifiers");
  Check(tf_percussion_parameter_scale(kick, pitch) == 1,
        "compact kick pitch declares logarithmic control scaling");

  const auto first = Render(kick, 17, 128);
  const auto repeated = Render(kick, 17, 128);
  const auto whole = Render(kick, 17, first.size());
  const auto variation = Render(kick, 18, 128);
  Check(first == repeated, "compact kick API is deterministic");
  Check(first == whole, "compact kick API is host-block independent");
  Check(Difference(first, variation) > 1.e-7,
        "compact kick API forwards event seeds");
  Check(std::all_of(first.begin(), first.end(),
                    [](const float sample) { return std::isfinite(sample); }) &&
            Energy(first) > 1.e-5,
        "compact kick API renders finite audible output");
  CheckPreparedRoundTrip(kick, 1, "prepared compact kick is sample-identical");

  for (std::uint32_t route = 0; route < 3; ++route)
    Check(tf_percussion_route_set(kick, route, 0.f),
          "compact kick route can be disabled");
  Check(tf_percussion_commit(kick), "compact kick routing commits");
  Check(Energy(Render(kick, 17, 128)) < 1.e-20,
        "disabling every compact kick source produces silence");
  for (std::uint32_t route = 0; route < 3; ++route)
    Check(tf_percussion_route_set(kick, route, 1.f),
          "compact kick route can be restored");
  Check(tf_percussion_parameter_set(
            kick, level, tf_percussion_parameter_default(kick, level) - 6.f) &&
            tf_percussion_commit(kick),
        "compact kick parameters commit through the shared API");
  const double quietRatio = Energy(Render(kick, 17, 128)) / Energy(first);
  Check(quietRatio > .24 && quietRatio < .26,
        "compact kick model level follows its decibel law");
  Check(!tf_percussion_parameter_set(kick, 99, 0.f) &&
            !tf_percussion_route_set(kick, 99, 1.f),
        "out-of-range compact kick edits are rejected");
  tf_percussion_destroy(kick);
  Check(!tf_percussion_process(kick, nullptr, 0),
        "destroyed recipe handles cannot be reused");

  const auto membrane = tf_percussion_create(2, 48000.f);
  Check(membrane != 0 && tf_percussion_recipe(membrane) == 2,
        "membrane session can be created");
  Check(tf_percussion_parameter_count(membrane) == 33 &&
            tf_percussion_route_count(membrane) == 5,
        "membrane exposes its bounded controls and routing");
  const auto membranePitch = ParameterIndex(membrane, "fundamental_hz");
  const auto eqMode = ParameterIndex(membrane, "equalizer_mode");
  const auto directVelocity =
      ParameterIndex(membrane, "direct_velocity_exponent");
  const auto bodyVelocity =
      ParameterIndex(membrane, "body_velocity_exponent");
  Check(membranePitch < 33 && eqMode < 33,
        "membrane parameters have stable identifiers");
  Check(tf_percussion_parameter_minimum(membrane, directVelocity) == 1.f &&
            tf_percussion_parameter_minimum(membrane, bodyVelocity) == 1.f &&
            tf_percussion_parameter_default(membrane, directVelocity) == 1.f &&
            tf_percussion_parameter_default(membrane, bodyVelocity) == 1.f,
        "membrane velocity controls cannot compress the input curve");
  Check(tf_percussion_parameter_scale(membrane, eqMode) == 3,
        "membrane output EQ declares a discrete choice control");
  Check(!tf_percussion_parameter_set(membrane, eqMode, 1.5f) &&
            tf_percussion_parameter_set(membrane, eqMode, 2.f) &&
            tf_percussion_parameter_set(membrane, eqMode, 1.f),
        "fractional choices are rejected and integral choices are accepted");
  const auto membraneFirst = Render(membrane, 71, 128);
  Check(Energy(membraneFirst) > 1.e-5 &&
            membraneFirst == Render(membrane, 71, 37),
        "membrane API is audible, deterministic, and block independent");
  CheckPreparedRoundTrip(membrane, 2,
                         "prepared membrane recipe is sample-identical");
  tf_percussion_destroy(membrane);

  const auto snare = tf_percussion_create(3, 48000.f);
  Check(snare != 0 && tf_percussion_recipe(snare) == 3,
        "snare session can be created");
  Check(tf_percussion_parameter_count(snare) == 50 &&
            tf_percussion_route_count(snare) == 7,
        "snare exposes membrane, wire, and routing controls");
  const auto wireLevel = ParameterIndex(snare, "wire_level");
  const auto wireDensity = ParameterIndex(snare, "wire_density");
  const auto snarePitch = ParameterIndex(snare, "fundamental_hz");
  Check(wireLevel < 50 && wireDensity < 50,
        "snare wire parameters have stable identifiers");
  Check(tf_percussion_parameter_default(snare, snarePitch) == 185.f &&
            tf_percussion_parameter_get(snare, snarePitch) == 185.f,
        "snare metadata and initialized state share fitted defaults");
  const auto snareFirst = Render(snare, 81, 128);
  Check(Energy(snareFirst) > 1.e-5 && snareFirst == Render(snare, 81, 37),
        "snare API is audible, deterministic, and block independent");
  CheckPreparedRoundTrip(snare, 3, "prepared snare is sample-identical");
  Check(tf_percussion_route_set(snare, 5, 0.f) && tf_percussion_commit(snare),
        "snare body-to-wire route can be disabled");
  const auto withoutWires = Render(snare, 81, 128);
  Check(Difference(snareFirst, withoutWires) > 1.e-6,
        "body-driven wire route changes snare output");
  tf_percussion_destroy(snare);

  std::array<std::uint32_t, 4> pool{};
  for (auto &pooled : pool)
    pooled = tf_percussion_create(0, 48000.f);
  Check(std::all_of(pool.begin(), pool.end(),
                    [](const auto pooled) { return pooled != 0; }) &&
            tf_percussion_create(0, 48000.f) == 0,
        "recipe sessions have a bounded four-instance pool");
  tf_percussion_destroy(pool[1]);
  pool[1] = tf_percussion_create(1, 48000.f);
  Check(pool[1] != 0 && tf_percussion_recipe(pool[1]) == 1,
        "destroyed recipe sessions are immediately reusable");
  for (const auto pooled : pool)
    tf_percussion_destroy(pooled);

  return failures == 0 ? 0 : 1;
}
