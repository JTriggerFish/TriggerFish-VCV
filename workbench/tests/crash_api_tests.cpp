#include "crash_api.hpp"
#include "crash_macros.hpp"
#include "percussion_api.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string_view>
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
                          const float location, const std::size_t blockSize,
                          const float implement = .75f,
                          const std::size_t frameCount = 8192,
                          const float hardness = .65f) {
  std::vector<float> result(frameCount);
  Check(tf_crash_reset(handle) == 1, "reset accepts a live handle");
  Check(tf_crash_trigger(handle, .8f, location, hardness, implement, .2f, seed) == 1,
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

double Energy(const std::vector<float> &audio) {
  double result = 0.0;
  for (const float sample : audio)
    result += static_cast<double>(sample) * sample;
  return result;
}

double HighFrequencyFraction(const std::vector<float> &audio) {
  constexpr double coefficient = 1.0 -
      std::exp(-6.283185307179586 * 5000.0 / 48000.0);
  double low = 0.0;
  double highEnergy = 0.0;
  double totalEnergy = 0.0;
  for (const float sample : audio) {
    low += coefficient * (sample - low);
    const double high = sample - low;
    highEnergy += high * high;
    totalEnergy += static_cast<double>(sample) * sample;
  }
  return highEnergy / std::max(totalEnergy, 1.e-30);
}

double CrestFactor(const std::vector<float> &audio) {
  const auto peak = *std::max_element(
      audio.begin(), audio.end(), [](const float first, const float second) {
        return std::abs(first) < std::abs(second);
      });
  const double rms = std::sqrt(Energy(audio) /
      std::max<std::size_t>(1, audio.size()));
  return std::abs(peak) / std::max(rms, 1.e-30);
}

bool SetRoutes(const std::uint32_t handle,
               const std::array<bool, 3> &enabled) {
  bool accepted = true;
  for (std::size_t index = 0; index < enabled.size(); ++index) {
    accepted &= tf_crash_route_enable(
                    handle, static_cast<std::uint32_t>(index), enabled[index]) == 1;
  }
  return accepted && tf_crash_macro_commit(handle) == 1;
}

} // namespace

int main() {
  Check(tf_crash_api_version() == 1, "unreleased API remains version one");
  Check(tf_crash_route_count() == 3,
        "the metallic recipe exposes three optional routes");
  Check(tf_crash_macro_count() == tfworkbench::CrashMacroCount,
        "the fitting surface contains only active unified parameters");
  Check(std::abs(tf_crash_macro_default(0) + 6.f) < 1.e-6f &&
            tf_crash_macro_maximum(0) == 0.f,
        "the crash workbench uses an attenuation-only -6 dB level with unity base gain");
  Check(tf_crash_macro_key(0) != nullptr &&
            tf_crash_macro_name(0) != nullptr &&
            tf_crash_macro_unit(0) != nullptr,
        "macro labels and units are available through the C ABI");
  Check(tf_crash_macro_scale(0) == 0 && tf_crash_macro_scale(2) == 1,
        "macro scales are available through the C ABI");
  Check(tf_crash_macro_name(tf_crash_macro_count()) == nullptr,
        "out-of-range macro metadata is rejected");
  auto painted = tfworkbench::DefaultCrashMacros();
  std::size_t firstFrequency = painted.size();
  std::size_t firstLevel = painted.size();
  std::size_t secondLevel = painted.size();
  std::size_t fieldTurbulence = painted.size();
  std::size_t bloomDiffusion = painted.size();
  std::size_t bloomRate = painted.size();
  std::size_t bloomEnergyAcceleration = painted.size();
  std::size_t bodyExcitation = painted.size();
  std::size_t firstDecayActive = painted.size();
  std::size_t firstModeTurbulence = painted.size();
  for (std::size_t index = 0; index < painted.size(); ++index) {
    const auto &key = tfworkbench::CrashMacroDescription(index).key;
    if (key == "resolved_frequency_0") firstFrequency = index;
    else if (key == "resolved_level_0") firstLevel = index;
    else if (key == "resolved_level_1") secondLevel = index;
    else if (key == "field_turbulence") fieldTurbulence = index;
    else if (key == "bloom_phase_diffusion") bloomDiffusion = index;
    else if (key == "bloom_rate") bloomRate = index;
    else if (key == "bloom_energy_acceleration")
      bloomEnergyAcceleration = index;
    else if (key == "body_excitation") bodyExcitation = index;
    else if (key == "body_decay_active_1") firstDecayActive = index;
    else if (key == "resolved_turbulence_0") firstModeTurbulence = index;
  }
  Check(firstFrequency < painted.size() && firstLevel < painted.size() &&
            secondLevel < painted.size() &&
            fieldTurbulence < painted.size() && bloomDiffusion < painted.size() &&
            bloomRate < painted.size() &&
            bloomEnergyAcceleration < painted.size() &&
            bodyExcitation < painted.size() &&
            firstDecayActive < painted.size() &&
            firstModeTurbulence < painted.size(),
        "modal-field macros have stable keys");
  Check(tfworkbench::CrashMacroDescription(firstFrequency).maximum == 15000.f,
        "the constructive modal editor stops at 15 kHz");
  painted[firstFrequency] = 177.f;
  painted[firstLevel] = 12.f;
  painted[secondLevel] = -12.f;
  painted[firstModeTurbulence] = 0.f;
  painted[bodyExcitation] = .2f;
  const auto baseFit = tfworkbench::MetallicWorkbenchBaseFit();
  Check(std::abs(baseFit.outputGain - 1.f) < 1.e-5f,
        "the neutral metallic base has no hidden output calibration gain");
  const auto paintedFit = tfworkbench::ApplyCrashMacros(baseFit, painted);
  const auto defaultFit = tfworkbench::ApplyCrashMacros(
      baseFit, tfworkbench::DefaultCrashMacros());
  Check(std::abs(paintedFit.sparseFrequencyHz[0] - 177.f) < 1.e-5f &&
            paintedFit.sparseAmplitude[0] > paintedFit.sparseAmplitude[1],
        "resolved editor directly places and levels resolved modes");
  Check(paintedFit.fieldTurbulenceScale[0] == 0.f,
        "each modal anchor owns a turbulence-response scaler");
  Check(std::abs(paintedFit.bodyExcitationGain - .2f) < 1.e-6f,
        "body excitation is an explicit modal-input parameter");
  Check(std::abs(defaultFit.bloomRateOctavesPerSecond - 2.f) < 1.e-5f &&
            defaultFit.contactDurationScale == .65f,
        "default macros resolve explicit contact and bloom values");
  auto noBloom = tfworkbench::DefaultCrashMacros();
  noBloom[bloomRate] = 0.f;
  const auto noBloomFit = tfworkbench::ApplyCrashMacros(baseFit, noBloom);
  Check(noBloomFit.bloomRateOctavesPerSecond == 0.f &&
            noBloomFit.bloomEnergyAcceleration ==
                defaultFit.bloomEnergyAcceleration,
        "bloom rate reaches zero without changing its energy response");
  Check(std::count(defaultFit.bodyDecayActive.begin(),
                   defaultFit.bodyDecayActive.end(), true) == 0,
        "the T60 envelope defaults to no optional interior knots");
  auto slowMacros = tfworkbench::DefaultCrashMacros();
  auto fastMacros = slowMacros;
  slowMacros[bloomRate] = 1.f;
  fastMacros[bloomRate] = 12.f;
  const auto slowFit = tfworkbench::ApplyCrashMacros(baseFit, slowMacros);
  const auto fastFit = tfworkbench::ApplyCrashMacros(baseFit, fastMacros);
  Check(slowFit.bloomRateOctavesPerSecond <
            fastFit.bloomRateOctavesPerSecond &&
            slowFit.bloomEnergyAcceleration == fastFit.bloomEnergyAcceleration &&
            slowFit.bloomPhaseDiffusion == fastFit.bloomPhaseDiffusion,
        "bloom rate is independent of energy response and phase diffusion");
  auto cleared = tfworkbench::DefaultCrashMacros();
  for (std::size_t mode = 0;
       mode < tfworkbench::ResolvedModePointCount; ++mode) {
    const auto index = firstLevel + mode;
    cleared[index] = tfworkbench::CrashMacroDescription(index).minimum;
  }
  const auto clearedFit = tfworkbench::ApplyCrashMacros(baseFit, cleared);
  Check(std::all_of(clearedFit.sparseAmplitude.begin(),
                    clearedFit.sparseAmplitude.end(),
                    [](const float amplitude) { return amplitude == 0.f; }),
        "clearing the modal editor deactivates every anchor exactly");
  Check(tf_crash_create(0.f) == 0, "invalid sample rates are rejected");
  const auto handle = tf_crash_create(48000.f);
  Check(handle != 0, "a renderer session can be created");
  Check(tf_percussion_parameter_count(handle) == tf_crash_macro_count(),
        "legacy and generic APIs expose the same complete control surface");
  for (std::uint32_t index = 0; index < tf_crash_macro_count(); ++index) {
    const auto *legacyKey = tf_crash_macro_key(index);
    const auto *genericKey = tf_percussion_parameter_key(handle, index);
    Check(legacyKey != nullptr && genericKey != nullptr &&
              std::string_view(legacyKey) == genericKey,
          "legacy and generic APIs expose controls in identical key order");
  }
  for (std::uint32_t index = 0; index < tf_crash_route_count(); ++index) {
    Check(tf_crash_route_enabled(handle, index) == 1,
          "metallic recipe routes default to enabled");
  }

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

  Check(SetRoutes(handle, {false, true, false}),
        "the direct-contact-only route is accepted");
  const auto directOnly = Render(handle, 17, .8f, 256);
  Check(Energy(directOnly) > 1.e-12,
        "the direct-contact-only route remains audible");
  Check(Difference(first, directOnly) > 1.e-8,
        "the direct route is distinct from the complete recipe");

  Check(SetRoutes(handle, {true, false, true}),
        "the modal-body-only route is accepted");
  const auto modalBody = Render(handle, 17, .8f, 256);
  Check(Energy(modalBody) > 1.e-12,
        "the modal body path remains audible");

  Check(SetRoutes(handle, {true, false, false}),
        "a silent observation ablation is accepted by the low-level API");
  const auto silent = Render(handle, 17, .8f, 256);
  Check(Energy(silent) < 1.e-20,
        "disabling both observation feeds produces silence");
  Check(tf_crash_route_enable(handle, tf_crash_route_count(), 1) == 0,
        "invalid route indices are rejected");
  Check(tf_crash_route_enable(handle, 0, 2) == 0,
        "non-boolean route states are rejected");
  Check(SetRoutes(handle, {true, true, true}),
        "factory routing can be restored after ablation");

  const auto brush = Render(handle, 17, .8f, 256, 0.f, 48000);
  const auto softBrush = Render(handle, 17, .8f, 256, 0.f, 48000, 0.f);
  const auto stiffBrush = Render(handle, 17, .8f, 256, 0.f, 48000, 1.f);
  const auto middleBrush = Render(handle, 17, .8f, 256, 0.f, 48000, .5f);
  const auto adjacentBrush = Render(handle, 17, .8f, 256, 0.f, 48000, .51f);
  const auto stick = Render(handle, 17, .8f, 256, 1.f, 48000);
  const double brushToStick = Energy(brush) / std::max(Energy(stick), 1.e-30);
  const double softToStick = Energy(softBrush) / std::max(Energy(stick), 1.e-30);
  const double stiffToStick = Energy(stiffBrush) / std::max(Energy(stick), 1.e-30);
  const double softBrightness = HighFrequencyFraction(softBrush);
  const double stiffBrightness = HighFrequencyFraction(stiffBrush);
  if (!(brushToStick > .25 && brushToStick < 4 &&
        softToStick > .25 && softToStick < 4 &&
        stiffToStick > .25 && stiffToStick < 4))
    std::cerr << "workbench brush/stick energy ratios: " << softToStick
              << ", " << brushToStick << ", " << stiffToStick << '\n';
  Check(brushToStick > .25 && brushToStick < 4 &&
            softToStick > .25 && softToStick < 4 &&
            stiffToStick > .25 && stiffToStick < 4,
        "brush/stick integrated levels stay within 6 dB without implement makeup");
  Check(stiffBrightness > 1.35 * softBrightness,
        "bristle stiffness audibly increases high-frequency articulation");
  Check(CrestFactor(stiffBrush) < 12.,
        "stiff bristle contacts remain a fused brush texture");
  Check(Difference(middleBrush, adjacentBrush) <
            .1 * Difference(softBrush, stiffBrush),
        "small bristle-stiffness moves preserve the contact realization");

  Check(tf_crash_macro_set(handle, fieldTurbulence, 0.f) &&
            tf_crash_macro_commit(handle),
        "the unified field can collapse to coherent anchors");
  const auto coherent = Render(handle, 17, .8f, 256);
  Check(tf_crash_macro_set(handle, fieldTurbulence, 1.f) &&
            tf_crash_macro_commit(handle),
        "the unified field accepts maximum turbulence");
  const auto diffuse = Render(handle, 17, .8f, 256);
  const double coherentEnergy =
      Difference(coherent, std::vector<float>(coherent.size()));
  Check(Difference(coherent, diffuse) > 1.e-4 * coherentEnergy,
        "unified turbulence materially changes the body response");
  Check(tf_crash_macro_set(handle, bloomDiffusion, 0.f) &&
            tf_crash_macro_commit(handle),
        "bloom transfer phase diffusion can be disabled");
  const auto focusedBloom = Render(handle, 17, .8f, 256);
  Check(tf_crash_macro_set(handle, bloomDiffusion, 1.f) &&
            tf_crash_macro_commit(handle),
        "bloom transfer can restore full phase diffusion");
  const auto diffuseBloom = Render(handle, 17, .8f, 256);
  const double bloomEnergy =
      Difference(focusedBloom, std::vector<float>(focusedBloom.size()));
  Check(Difference(focusedBloom, diffuseBloom) > 1.e-4 * bloomEnergy,
        "bloom diffusion materially changes the body excitation");
  for (std::size_t index = 0; index < tf_crash_macro_count(); ++index)
    Check(tf_crash_macro_set(handle, index, tf_crash_macro_default(index)),
          "factory controls can be restored after ablation");
  Check(tf_crash_macro_commit(handle),
        "restored factory controls prepare successfully");

  const auto controlHandle = tf_crash_create(48000.f);
  Check(controlHandle != 0, "a second live renderer can be created");
  std::vector<float> liveFirst(4096), liveSecond(4096), controlSecond(4096);
  Check(tf_crash_reset(handle) && tf_crash_reset(controlHandle),
        "live restrike test resets both renderers");
  Check(tf_crash_trigger(handle, .8f, .8f, .65f, .75f, .2f, 31) &&
            tf_crash_trigger(controlHandle, .8f, .8f, .65f, .75f, .2f, 31),
        "live restrike test starts equal states");
  Check(tf_crash_process(handle, liveFirst.data(), liveFirst.size()) &&
            tf_crash_process(controlHandle, liveFirst.data(), liveFirst.size()),
        "live restrike prefix renders");
  Check(tf_crash_trigger(handle, .7f, .8f, .65f, .75f, .2f, 32),
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
