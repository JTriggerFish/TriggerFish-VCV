#include "tfdsp/reverb_output.hpp"
#include "tfdsp/room_reverb.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double SampleRate = 48'000.0;

void Check(const bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}

void TestRoomGeometryDerivesHeightFromSizeAndUsesAspectAndListener() {
  tfdsp::RoomReverbControls controls;
  const auto neutral = tfdsp::RoomReverb::MakeRoom(controls);
  controls.aspect = 1.f;
  controls.listener = {0.2f, 0.3f, 0.4f};
  const auto shaped = tfdsp::RoomReverb::MakeRoom(controls);
  Check(shaped.dimensionsMetres[0] > neutral.dimensionsMetres[0] &&
            shaped.dimensionsMetres[1] < neutral.dimensionsMetres[1],
        "aspect must reshape the footprint while preserving its scale");
  controls.space = 1.f;
  const auto large = tfdsp::RoomReverb::MakeRoom(controls);
  Check(large.dimensionsMetres[2] > neutral.dimensionsMetres[2],
        "Size must derive a monotonically taller ceiling");
  for (std::size_t axis = 0; axis < 3; ++axis)
    Check(std::abs(shaped.listenerPositionMetres[axis] /
                       shaped.dimensionsMetres[axis] -
                   controls.listener[axis]) < 1.e-6,
          "listener controls must map to normalized room coordinates");
}

void TestRealtimeSceneConstructionSanitizesNonFiniteControls() {
  tfdsp::RoomReverb reverb;
  reverb.SetSampleRate(SampleRate);
  tfdsp::RoomReverbControls controls;
  controls.space = std::numeric_limits<float>::quiet_NaN();
  controls.aspect = std::numeric_limits<float>::infinity();
  controls.damping = -std::numeric_limits<float>::infinity();
  controls.listener = {std::numeric_limits<float>::quiet_NaN(),
                       std::numeric_limits<float>::infinity(),
                       -std::numeric_limits<float>::infinity()};
  tfdsp::RoomReverb::SourcePositions positions{};
  positions[0] = {std::numeric_limits<float>::quiet_NaN(),
                  std::numeric_limits<float>::infinity(),
                  -std::numeric_limits<float>::infinity()};
  tfdsp::RoomReverb::InputFrame input{};
  input[0] = 1.f;
  const auto output = reverb.Process(input, positions, 1, controls);
  Check(std::isfinite(output.direct[0]) && std::isfinite(output.direct[1]) &&
            std::isfinite(output.wet[0]) && std::isfinite(output.wet[1]),
        "the noexcept realtime scene path must sanitize non-finite controls");
}

struct RenderResult {
  std::vector<tfdsp::RoomReverb::StereoFrame> response{};
  double firstHalfEnergy{};
  double secondHalfEnergy{};
  float peak{};
  std::size_t firstNonzero{};
};

RenderResult Render(tfdsp::RoomReverbControls controls,
                    const std::size_t samples = 48'000,
                    const tfdsp::RoomReverb::SourcePosition &sourcePosition =
                        tfdsp::reverb_defaults::Source) {
  tfdsp::RoomReverb reverb;
  reverb.SetSampleRate(SampleRate);
  tfdsp::RoomReverb::SourcePositions positions{};
  positions[0] = sourcePosition;
  tfdsp::RoomReverb::InputFrame silence{};
  reverb.Process(silence, positions, 1, controls);
  Check(reverb.WaitForLatestScene(std::chrono::seconds(3)),
        "background ER scene must build and publish");
  // Adoption happens on an audio block boundary; allow the matching FIR and
  // handoff metadata to become the accepted scene before clearing histories.
  for (std::size_t sample = 0; sample < 256; ++sample)
    reverb.Process(silence, positions, 1, controls);
  reverb.Reset();

  RenderResult result;
  result.response.resize(samples);
  result.firstNonzero = samples;
  for (std::size_t sample = 0; sample < samples; ++sample) {
    tfdsp::RoomReverb::InputFrame input{};
    input[0] = sample == 0 ? 1.f : 0.f;
    result.response[sample] =
        reverb.Process(input, positions, 1, controls).wet;
    for (const float value : result.response[sample]) {
      Check(std::isfinite(value), "combined room response must remain finite");
      result.peak = std::max(result.peak, std::abs(value));
      if (std::abs(value) > 1.e-8f)
        result.firstNonzero = std::min(result.firstNonzero, sample);
      const double energy = static_cast<double>(value) * value;
      if (sample < samples / 2)
        result.firstHalfEnergy += energy;
      else
        result.secondHalfEnergy += energy;
    }
  }
  return result;
}

double WindowEnergy(const RenderResult &render, const std::size_t first = 0,
                    std::size_t last = static_cast<std::size_t>(-1)) {
  double energy = 0.0;
  last = std::min(last, render.response.size());
  for (std::size_t sample = first; sample < last; ++sample)
    for (const float value : render.response[sample])
      energy += static_cast<double>(value) * value;
  return energy;
}

double DifferenceEnergy(const RenderResult &left, const RenderResult &right) {
  const std::size_t samples =
      std::min(left.response.size(), right.response.size());
  double energy = 0.0;
  for (std::size_t sample = 0; sample < samples; ++sample)
    for (std::size_t channel = 0; channel < 2; ++channel) {
      const double difference =
          left.response[sample][channel] - right.response[sample][channel];
      energy += difference * difference;
    }
  return energy;
}

void CheckResponseChanged(const RenderResult &left, const RenderResult &right,
                          const std::string &control,
                          const double minimumRelativeEnergy = 1.e-7) {
  const double reference =
      std::max({WindowEnergy(left), WindowEnergy(right), 1.e-20});
  const double relative = DifferenceEnergy(left, right) / reference;
  Check(relative > minimumRelativeEnergy,
        control +
            " must audibly change the impulse response; relative "
            "difference=" +
            std::to_string(relative));
}

double BandEnergy(const RenderResult &render, const double lowerFrequency,
                  const double upperFrequency, const std::size_t first = 0) {
  std::array<double, 2> lowState{};
  std::array<double, 2> highState{};
  const double lowAlpha = 1.0 - std::exp(-2.0 * 3.14159265358979323846 *
                                         lowerFrequency / SampleRate);
  const double highAlpha = 1.0 - std::exp(-2.0 * 3.14159265358979323846 *
                                          upperFrequency / SampleRate);
  double energy = 0.0;
  for (std::size_t sample = 0; sample < render.response.size(); ++sample)
    for (std::size_t channel = 0; channel < 2; ++channel) {
      const double value = render.response[sample][channel];
      lowState[channel] += lowAlpha * (value - lowState[channel]);
      highState[channel] += highAlpha * (value - highState[channel]);
      if (sample >= first) {
        const double band = highState[channel] - lowState[channel];
        energy += band * band;
      }
    }
  return energy;
}

double SideToMidEnergy(const RenderResult &render) {
  double mid = 0.0;
  double side = 0.0;
  for (const auto &frame : render.response) {
    const double m = frame[0] + frame[1];
    const double s = frame[0] - frame[1];
    mid += m * m;
    side += s * s;
  }
  return side / std::max(mid, 1.e-20);
}

void TestCanonicalFactoryDefaults() {
  tfdsp::RoomReverb processor;
  Check(processor.LateFlavour() == tfdsp::LateReverbFlavour::Optimized,
        "the complete room must default to the Optimized late-tail flavour");
  processor.SetLateReverbFlavour(tfdsp::LateReverbFlavour::Base);
  Check(processor.LateFlavour() == tfdsp::LateReverbFlavour::Base,
        "the complete room must forward the Base late-tail selection");
  processor.SetLateReverbFlavourImmediate(
      tfdsp::LateReverbFlavour::Optimized);
  Check(processor.LateFlavour() == tfdsp::LateReverbFlavour::Optimized,
        "the complete room must forward an immediate restored flavour");

  const tfdsp::RoomReverbControls room;
  Check(room.space == tfdsp::reverb_defaults::Space &&
            room.aspect == tfdsp::reverb_defaults::Aspect &&
            room.listener == tfdsp::reverb_defaults::Listener &&
            room.preDelay == tfdsp::reverb_defaults::PreDelay &&
            room.decay == tfdsp::reverb_defaults::Decay &&
            room.damping == tfdsp::reverb_defaults::Damping &&
            room.diffusion == tfdsp::reverb_defaults::Diffusion &&
            room.modulation == tfdsp::reverb_defaults::Modulation &&
            room.shimmer == tfdsp::reverb_defaults::Shimmer &&
            room.width == tfdsp::reverb_defaults::Width &&
            room.balance == tfdsp::reverb_defaults::Balance &&
            room.lowCut == tfdsp::reverb_defaults::LowCut &&
            room.highCut == tfdsp::reverb_defaults::HighCut,
        "RoomReverbControls must use the canonical module baseline");

  const tfdsp::LateReverbControls late;
  Check(
      late.decay == room.decay && late.damping == room.damping &&
          late.diffusion == room.diffusion &&
          late.modulation == room.modulation && late.shimmer == room.shimmer &&
          late.roomDimensionsMetres ==
              tfdsp::reverb_defaults::RoomDimensionsMetres,
      "standalone late-reverb defaults must match the complete room baseline");

  const auto mapped = tfdsp::RoomReverb::MakeRoom(room);
  for (std::size_t axis = 0; axis < 3; ++axis)
    Check(std::abs(mapped.dimensionsMetres[axis] -
                   tfdsp::reverb_defaults::RoomDimensionsMetres[axis]) < 1.e-5,
          "canonical late-room dimensions must match MakeRoom");
}

void TestProgressiveDefaultSourceSpread() {
  using tfdsp::reverb_defaults::ProgressiveSourceX;
  Check(ProgressiveSourceX(0, 1) == 0.5f,
        "one default source must be horizontally centred");

  float previousHalfSpan = 0.f;
  for (std::size_t sourceCount = 2; sourceCount <= 8; ++sourceCount) {
    const float left = ProgressiveSourceX(0, sourceCount);
    const float right = ProgressiveSourceX(sourceCount - 1, sourceCount);
    Check(std::abs(left + right - 1.f) < 1.e-6f,
          "default polyphonic source spread must remain symmetric");
    Check(left < 0.5f && right > 0.5f,
          "multiple default sources must straddle the room centre");
    Check(0.5f - left > previousHalfSpan,
          "each added source must progressively widen the default scene");
    previousHalfSpan = 0.5f - left;
    for (std::size_t source = 1; source < sourceCount; ++source)
      Check(ProgressiveSourceX(source - 1, sourceCount) <
                ProgressiveSourceX(source, sourceCount),
            "default source positions must increase from left to right");
  }
}

void TestAcousticPresetDefinitions() {
  using tfdsp::reverb_defaults::MediumHall;
  using tfdsp::reverb_defaults::SmallRoom;
  using tfdsp::reverb_defaults::Superlush;
  const auto makeRoom = [](const tfdsp::reverb_defaults::ReverbPreset &preset) {
    tfdsp::RoomReverbControls controls;
    controls.space = preset.space;
    controls.aspect = preset.aspect;
    controls.listener = preset.listener;
    return tfdsp::RoomReverb::MakeRoom(controls);
  };
  const auto volume = [](const tfdsp::EarlyReflectionRoom &room) {
    return room.dimensionsMetres[0] * room.dimensionsMetres[1] *
           room.dimensionsMetres[2];
  };
  const auto small = makeRoom(SmallRoom);
  const auto medium = makeRoom(MediumHall);
  const auto lush = makeRoom(Superlush);
  Check(volume(small) < volume(medium) && volume(medium) < volume(lush),
        "Small Room, Medium Hall, and Superlush must have increasing volume");
  Check(SmallRoom.decay < MediumHall.decay &&
            MediumHall.decay < Superlush.decay,
        "preset decay times must increase with their intended acoustic scale");
  Check(SmallRoom.diffusion < MediumHall.diffusion &&
            MediumHall.diffusion < Superlush.diffusion,
        "preset diffusion must progress from room ambience to maximum lushness");
  Check(SmallRoom.modulation < MediumHall.modulation &&
            MediumHall.modulation < Superlush.modulation &&
            Superlush.modulation == 1.f &&
            std::abs(Superlush.shimmer - 0.45f) < 1.e-7f,
        "Superlush must use maximum modulation and a strong 45% shimmer while "
        "room presets remain restrained");
  Check(SmallRoom.mix < MediumHall.mix && MediumHall.mix < Superlush.mix &&
            std::abs(Superlush.mix - 0.40f) < 1.e-7f,
        "Superlush must add more wet texture while remaining dry-forward");
  Check(SmallRoom.balance < MediumHall.balance && MediumHall.balance == 0.5f &&
            Superlush.balance > MediumHall.balance,
        "Small Room must lean toward early reflections, Medium Hall must use "
        "the inferred centre, and Superlush must lean toward the late tail");
  const auto lowCutHertz = [](const float value) {
    return 20.f * std::exp(value * std::log(50.f));
  };
  const auto highCutHertz = [](const float value) {
    return 1'000.f * std::exp(value * std::log(20.f));
  };
  Check(std::abs(lowCutHertz(SmallRoom.lowCut) - 40.f) < 0.01f &&
            std::abs(highCutHertz(SmallRoom.highCut) - 9'000.f) < 1.f &&
            std::abs(lowCutHertz(Superlush.lowCut) - 120.f) < 0.01f &&
            std::abs(highCutHertz(Superlush.highCut) - 8'000.f) < 1.f,
        "room and texture preset filters must retain their documented cutoffs");

  for (const auto &[name, preset] :
       {std::pair<const char *, tfdsp::reverb_defaults::ReverbPreset>{
            "Medium Hall", MediumHall},
        {"Small Room", SmallRoom},
        {"Superlush", Superlush}}) {
    tfdsp::RoomReverbControls controls;
    controls.space = preset.space;
    controls.aspect = preset.aspect;
    controls.listener = preset.listener;
    controls.preDelay = preset.preDelay;
    controls.decay = preset.decay;
    controls.damping = preset.damping;
    controls.diffusion = preset.diffusion;
    controls.modulation = preset.modulation;
    controls.shimmer = preset.shimmer;
    controls.width = preset.width;
    controls.balance = preset.balance;
    controls.lowCut = preset.lowCut;
    controls.highCut = preset.highCut;
    const auto response = Render(controls, 96'000, preset.source);
    Check(response.peak > 1.e-5f && response.peak < 4.f &&
              response.firstHalfEnergy > 0.0 && response.secondHalfEnergy > 0.0,
          std::string(name) +
              " preset must produce a useful, bounded, persistent wet response");
  }
}

void TestCombinedEarlyAndLateResponse() {
  const auto combined = Render({});
  Check(combined.firstNonzero < 8'000,
        "the room response must produce an audible onset");
  Check(combined.peak > 1.e-4f && combined.peak < 4.f,
        "the default room response must have useful bounded level");
  Check(combined.firstHalfEnergy > 0.0 && combined.secondHalfEnergy > 0.0,
        "the room response must contain both an onset and a persistent tail");
  bool stereo = false;
  for (const auto &frame : combined.response)
    stereo |= std::abs(frame[0] - frame[1]) > 1.e-7f;
  Check(stereo, "positioned early and late fields must decode to stereo");

  tfdsp::RoomReverbControls monoControls;
  monoControls.width = 0.f;
  const auto mono = Render(monoControls, 12'000);
  for (const auto &frame : mono.response)
    Check(std::abs(frame[0] - frame[1]) < 1.e-6f,
          "zero Width must collapse the complete wet room response to mono");
}

void TestEarlyTailIsolationAndPreDelay() {
  tfdsp::RoomReverbControls earlyOnly;
  earlyOnly.balance = 0.f;
  const auto early = Render(earlyOnly, 24'000);
  tfdsp::RoomReverbControls tailOnly;
  tailOnly.balance = 1.f;
  const auto tail = Render(tailOnly, 24'000);
  tfdsp::RoomReverbControls inferred;
  inferred.balance = 0.5f;
  const auto centre = Render(inferred, 24'000);
  Check(early.firstHalfEnergy > 1.e-6,
        "the left balance endpoint must expose the geometric response");
  Check(tail.firstHalfEnergy > 1.e-6,
        "the right balance endpoint must expose the FDN response");
  for (std::size_t sample = 0; sample < centre.response.size(); ++sample)
    for (std::size_t channel = 0; channel < 2; ++channel)
      Check(std::abs(centre.response[sample][channel] -
                     early.response[sample][channel] -
                     tail.response[sample][channel]) < 2.e-5f,
            "the centre balance position must preserve the inferred sum");

  tfdsp::RoomReverbControls delayed = tailOnly;
  delayed.preDelay = 0.5f;
  const auto withPreDelay = Render(delayed, 24'000);
  Check(withPreDelay.firstNonzero > tail.firstNonzero + 5'000,
        "shared pre-delay must postpone the wet response");

  std::cout << "Room response calibration: early/tail energy 0-20ms="
            << WindowEnergy(early, 0, 960) << "/" << WindowEnergy(tail, 0, 960)
            << " 20-80ms=" << WindowEnergy(early, 960, 3'840) << "/"
            << WindowEnergy(tail, 960, 3'840)
            << " 80-200ms=" << WindowEnergy(early, 3'840, 9'600) << "/"
            << WindowEnergy(tail, 3'840, 9'600)
            << " 200-500ms=" << WindowEnergy(early, 9'600, 24'000) << "/"
            << WindowEnergy(tail, 9'600, 24'000) << " peaks=" << early.peak
            << "/" << tail.peak << '\n';
}

void TestGeometryControlsChangeTheResponseIndependently() {
  constexpr std::size_t samples = 18'000;
  const tfdsp::RoomReverbControls baselineControls;
  const auto baseline = Render(baselineControls, samples);

  std::array<RenderResult, 2> sizeResponses;
  for (std::size_t index = 0; index < sizeResponses.size(); ++index) {
    auto controls = baselineControls;
    controls.space = index == 0 ? 0.f : 1.f;
    sizeResponses[index] = Render(controls, samples);
    CheckResponseChanged(baseline, sizeResponses[index], "Size");
  }
  Check(sizeResponses[0].firstNonzero < baseline.firstNonzero &&
            baseline.firstNonzero < sizeResponses[1].firstNonzero,
        "Size must monotonically postpone the first geometric arrival");

  for (const float aspect : {0.f, 1.f}) {
    auto controls = baselineControls;
    controls.aspect = aspect;
    CheckResponseChanged(baseline, Render(controls, samples), "Aspect");
    const auto room = tfdsp::RoomReverb::MakeRoom(controls);
    const auto reference = tfdsp::RoomReverb::MakeRoom(baselineControls);
    Check(std::abs(room.dimensionsMetres[0] * room.dimensionsMetres[1] -
                   reference.dimensionsMetres[0] *
                       reference.dimensionsMetres[1]) < 1.e-5,
          "Aspect must preserve floor area while changing the IR");
  }

}

void TestSizeSweepScalesAllDimensionsWithoutBecomingALevelControl() {
  constexpr std::size_t samples = 72'000;
  constexpr std::array<float, 5> sizes{0.f, 0.25f, 0.5f, 0.75f, 1.f};
  std::array<double, sizes.size()> energy{};
  std::array<tfdsp::EarlyReflectionRoom, sizes.size()> rooms{};
  for (std::size_t index = 0; index < sizes.size(); ++index) {
    tfdsp::RoomReverbControls controls;
    controls.space = sizes[index];
    rooms[index] = tfdsp::RoomReverb::MakeRoom(controls);
    energy[index] = WindowEnergy(Render(controls, samples));
    if (index > 0)
      for (std::size_t axis = 0; axis < 3; ++axis)
        Check(rooms[index].dimensionsMetres[axis] >
                  rooms[index - 1].dimensionsMetres[axis],
              "Size must monotonically increase length, width, and height");
  }

  const double minimum = *std::min_element(energy.begin(), energy.end());
  const double maximum = *std::max_element(energy.begin(), energy.end());
  const double spreadDb = 10.0 * std::log10(maximum / std::max(minimum, 1.e-20));
  Check(spreadDb < 6.0,
        "Size compensation must prevent room scale from becoming an "
        "unintended wet-level control; spread=" +
            std::to_string(spreadDb) + " dB");
}

void TestSourceAndListenerPositionSweepsChangeTheResponse() {
  constexpr std::size_t samples = 18'000;
  const tfdsp::RoomReverbControls baselineControls;
  const auto baseline = Render(baselineControls, samples);
  for (const std::size_t axis : {std::size_t{0}, std::size_t{1}}) {
    for (const float coordinate : {0.15f, 0.85f}) {
      auto source = tfdsp::reverb_defaults::Source;
      source[axis] = coordinate;
      CheckResponseChanged(baseline, Render(baselineControls, samples, source),
                           axis == 0 ? "Source X" : "Source Y");

      auto controls = baselineControls;
      controls.listener[axis] = coordinate;
      CheckResponseChanged(baseline, Render(controls, samples),
                           axis == 0 ? "Listener X" : "Listener Y");
    }
  }
}

std::array<float, 2>
RenderDirectAt(const tfdsp::RoomReverbControls &controls,
               const tfdsp::RoomReverb::SourcePosition &position) {
  tfdsp::RoomReverb reverb;
  reverb.SetSampleRate(SampleRate);
  tfdsp::RoomReverb::InputFrame input{};
  input[0] = 1.f;
  tfdsp::RoomReverb::SourcePositions positions{};
  positions[0] = position;
  return reverb.Process(input, positions, 1, controls).direct;
}

void TestPositionedDirectSound() {
  tfdsp::RoomReverbControls controls;
  const auto left = RenderDirectAt(controls, {0.15f, 0.35f, 0.45f});
  const auto right = RenderDirectAt(controls, {0.85f, 0.35f, 0.45f});
  Check(std::abs(left[0]) > std::abs(left[1]) &&
            std::abs(right[1]) > std::abs(right[0]),
        "a source left of the listener favours the left speaker and a source "
        "right of the listener favours the right speaker on its input sample");

  tfdsp::RoomReverb polyphonic;
  polyphonic.SetSampleRate(SampleRate);
  tfdsp::RoomReverb::InputFrame polyInput{};
  polyInput[0] = polyInput[1] = 1.f;
  tfdsp::RoomReverb::SourcePositions polyPositions{};
  polyPositions[0] = {0.15f, 0.35f, 0.45f};
  polyPositions[1] = {0.85f, 0.35f, 0.45f};
  const auto poly =
      polyphonic.Process(polyInput, polyPositions, 2, controls).direct;
  Check(std::abs(poly[0] - (left[0] + right[0])) < 1.e-6f &&
            std::abs(poly[1] - (left[1] + right[1])) < 1.e-6f,
        "polyphonic sources are independently positioned before summation");

  const auto near = RenderDirectAt(controls, {0.5f, 0.66f, 0.45f});
  const auto factory =
      RenderDirectAt(controls, tfdsp::reverb_defaults::Source);
  const auto far = RenderDirectAt(controls, {0.05f, 0.05f, 0.45f});
  const auto power = [](const std::array<float, 2> &frame) {
    return frame[0] * frame[0] + frame[1] * frame[1];
  };
  Check(power(near) > power(factory) && power(factory) > power(far),
        "bounded positioned-direct gain decreases with normalized distance");
  Check(std::abs(std::sqrt(power(factory)) - 1.f) < 1.e-5f,
        "factory positioned-direct distance has unity gain");
  Check(std::sqrt(power(near)) <= 2.f + 1.e-6f &&
            std::sqrt(power(far)) >= 0.25f - 1.e-6f,
        "positioned-direct distance gain stays within +6/-12 dB bounds");

  controls.preDelay = 1.f;
  const auto factoryWithMaximumPreDelay =
      RenderDirectAt(controls, tfdsp::reverb_defaults::Source);
  Check(std::abs(factoryWithMaximumPreDelay[0] - factory[0]) < 1.e-7f &&
            std::abs(factoryWithMaximumPreDelay[1] - factory[1]) < 1.e-7f,
        "wet pre-delay must not delay or otherwise alter the direct sample");

}

void TestResetPreservesAcceptedHandoffMetadata() {
  tfdsp::RoomReverb reverb;
  reverb.SetSampleRate(SampleRate);
  tfdsp::RoomReverbControls controls;
  controls.balance = 1.f;
  tfdsp::RoomReverb::SourcePositions positions{};
  positions[0] = {0.18f, 0.22f, 0.45f};
  tfdsp::RoomReverb::InputFrame silence{};
  reverb.Process(silence, positions, 1, controls);
  Check(reverb.WaitForLatestScene(std::chrono::seconds(3)),
        "reset metadata test must receive its prepared scene");

  // Let the queued FIR become active and the handoff gain finish its normal
  // scene-change crossfade before comparing two clean resets.
  for (std::size_t sample = 0; sample < 6'000; ++sample)
    reverb.Process(silence, positions, 1, controls);

  const auto renderAfterReset = [&] {
    reverb.Reset();
    std::vector<tfdsp::RoomReverb::StereoFrame> response(24'000);
    for (std::size_t sample = 0; sample < response.size(); ++sample) {
      tfdsp::RoomReverb::InputFrame input{};
      input[0] = sample == 0 ? 1.f : 0.f;
      response[sample] = reverb.Process(input, positions, 1, controls).wet;
    }
    return response;
  };

  const auto first = renderAfterReset();
  const auto second = renderAfterReset();
  for (std::size_t sample = 0; sample < first.size(); ++sample)
    for (std::size_t channel = 0; channel < 2; ++channel)
      Check(std::abs(first[sample][channel] - second[sample][channel]) < 1.e-7f,
            "Reset must retain the accepted handoff gain and alignment");
}

void TestPreDelayAutomationCrossfadesWithoutAnEnergyHole() {
  tfdsp::RoomReverb reverb;
  reverb.SetSampleRate(SampleRate);
  tfdsp::RoomReverbControls controls;
  controls.balance = 1.f;
  controls.decay = 0.8f;
  tfdsp::RoomReverb::SourcePositions positions{};
  positions[0] = tfdsp::reverb_defaults::Source;
  tfdsp::RoomReverb::InputFrame silence{};
  reverb.Process(silence, positions, 1, controls);
  Check(reverb.WaitForLatestScene(std::chrono::seconds(3)),
        "pre-delay automation test scene must publish");
  reverb.Reset();

  constexpr std::size_t changeSample = 14'400;
  constexpr std::size_t window = 480;
  constexpr std::size_t sampleCount = 24'000;
  std::array<double, 10> before{};
  std::array<double, 10> after{};
  std::uint32_t noise = 0x243f6a88u;
  for (std::size_t sample = 0; sample < sampleCount; ++sample) {
    noise ^= noise << 13;
    noise ^= noise >> 17;
    noise ^= noise << 5;
    tfdsp::RoomReverb::InputFrame input{};
    input[0] = 0.20f *
               (2.f * static_cast<float>(noise & 0x00ffffffu) /
                    static_cast<float>(0x00ffffffu) -
                1.f);
    if (sample >= changeSample)
      controls.preDelay = 1.f;
    const auto output = reverb.Process(input, positions, 1, controls).wet;
    const double energy = static_cast<double>(output[0]) * output[0] +
                          static_cast<double>(output[1]) * output[1];
    if (sample >= changeSample - before.size() * window &&
        sample < changeSample)
      before[(sample - (changeSample - before.size() * window)) / window] +=
          energy;
    if (sample >= changeSample &&
        sample < changeSample + after.size() * window)
      after[(sample - changeSample) / window] += energy;
  }
  double beforeMean = 0.0;
  for (const double value : before)
    beforeMean += value;
  beforeMean /= before.size();
  const double minimumAfter = *std::min_element(after.begin(), after.end());
  Check(minimumAfter > 0.03 * beforeMean,
        "moving wet pre-delay must crossfade live read heads instead of "
        "abandoning buffer history; minimum-window ratio=" +
            std::to_string(minimumAfter / std::max(beforeMean, 1.e-20)));
}

void TestPreDelaySweepTranslatesTheWetImpulse() {
  std::array<RenderResult, 3> responses;
  const std::array<float, 3> values{
      0.f,
      1.f / static_cast<float>(tfdsp::RoomReverb::MaximumPreDelaySeconds *
                               SampleRate),
      0.25f};
  for (std::size_t index = 0; index < values.size(); ++index) {
    tfdsp::RoomReverbControls controls;
    controls.preDelay = values[index];
    responses[index] = Render(controls, 30'000);
  }
  for (std::size_t index = 1; index < responses.size(); ++index) {
    const auto measured =
        static_cast<long long>(responses[index].firstNonzero) -
        static_cast<long long>(responses[0].firstNonzero);
    const auto expected = static_cast<long long>(
        std::llround(values[index] * tfdsp::RoomReverb::MaximumPreDelaySeconds *
                     SampleRate));
    const auto tolerance = index == 1 ? 0ll : 2ll;
    Check(std::abs(measured - expected) <= tolerance,
          "Pre-delay must translate the wet IR by the requested samples; "
          "measured=" +
              std::to_string(measured) +
              ", expected=" + std::to_string(expected));
  }
}

void TestDecayAndDampingSweepsHaveExpectedSpectralEffects() {
  std::array<double, 3> decayTailEnergy{};
  const std::array<float, 3> decayValues{0.25f, tfdsp::reverb_defaults::Decay,
                                         0.85f};
  for (std::size_t index = 0; index < decayValues.size(); ++index) {
    tfdsp::RoomReverbControls controls;
    controls.balance = 1.f;
    controls.decay = decayValues[index];
    const auto response = Render(controls, 96'000);
    decayTailEnergy[index] = WindowEnergy(response, 24'000);
  }
  Check(decayTailEnergy[0] < decayTailEnergy[1] &&
            decayTailEnergy[1] < decayTailEnergy[2],
        "Decay must monotonically increase late impulse-response energy");

  std::array<double, 3> highToMid{};
  const std::array<float, 3> dampingValues{0.f, tfdsp::reverb_defaults::Damping,
                                           1.f};
  for (std::size_t index = 0; index < dampingValues.size(); ++index) {
    tfdsp::RoomReverbControls controls;
    controls.balance = 1.f;
    controls.damping = dampingValues[index];
    const auto response = Render(controls, 96'000);
    highToMid[index] =
        BandEnergy(response, 4'000.0, 12'000.0, 4'800) /
        std::max(BandEnergy(response, 400.0, 2'000.0, 4'800), 1.e-20);
  }
  Check(highToMid[0] > highToMid[1] && highToMid[1] > highToMid[2],
        "Damping must monotonically shorten the high-frequency tail relative "
        "to the midband");
}

void TestDiffusionModulationAndShimmerSweeps() {
  constexpr std::size_t samples = 48'000;
  tfdsp::RoomReverbControls referenceControls;
  referenceControls.balance = 1.f;
  const auto baseline = Render(referenceControls, samples);

  for (const float diffusion : {0.f, 0.5f, 1.f}) {
    if (diffusion == tfdsp::reverb_defaults::Diffusion)
      continue;
    auto controls = referenceControls;
    controls.diffusion = diffusion;
    CheckResponseChanged(baseline, Render(controls, samples), "Diffusion");
  }

  std::array<RenderResult, 3> modulationResponses;
  const std::array<float, 3> modulationValues{
      0.f, tfdsp::reverb_defaults::Modulation, 1.f};
  for (std::size_t index = 0; index < modulationValues.size(); ++index) {
    auto controls = referenceControls;
    controls.modulation = modulationValues[index];
    modulationResponses[index] = Render(controls, samples);
  }
  CheckResponseChanged(modulationResponses[0], modulationResponses[1],
                       "Default modulation", 1.e-12);
  CheckResponseChanged(modulationResponses[1], modulationResponses[2],
                       "Extended modulation");

  std::array<RenderResult, 3> shimmerResponses;
  const std::array<float, 3> shimmerValues{0.f, 0.5f, 1.f};
  for (std::size_t index = 0; index < shimmerValues.size(); ++index) {
    auto controls = referenceControls;
    controls.shimmer = shimmerValues[index];
    shimmerResponses[index] = Render(controls, samples);
  }
  CheckResponseChanged(shimmerResponses[0], shimmerResponses[1], "Shimmer");
  CheckResponseChanged(shimmerResponses[1], shimmerResponses[2], "Shimmer");
  for (std::size_t index = 1; index < shimmerResponses.size(); ++index)
    Check(WindowEnergy(shimmerResponses[index], 6'000) >=
              0.5 * WindowEnergy(shimmerResponses[0], 6'000),
          "Shimmer must preserve the underlying late tail at every tested "
          "setting");
}

void TestWidthLevelAndFilterSweeps() {
  std::array<double, 3> sideRatios{};
  const std::array<float, 3> widths{0.f, tfdsp::reverb_defaults::Width, 1.f};
  for (std::size_t index = 0; index < widths.size(); ++index) {
    tfdsp::RoomReverbControls controls;
    controls.width = widths[index];
    sideRatios[index] = SideToMidEnergy(Render(controls, 24'000));
  }
  Check(sideRatios[0] < 1.e-12 && sideRatios[0] < sideRatios[1] &&
            sideRatios[1] < sideRatios[2],
        "Width must monotonically increase wet side energy from mono through "
        "native to 150%");

  std::array<double, 3> lowBandRatios{};
  const std::array<float, 3> lowCuts{0.f, 0.5f, 1.f};
  for (std::size_t index = 0; index < lowCuts.size(); ++index) {
    tfdsp::RoomReverbControls controls;
    controls.lowCut = lowCuts[index];
    const auto response = Render(controls, 48'000);
    lowBandRatios[index] =
        BandEnergy(response, 25.0, 120.0) /
        std::max(BandEnergy(response, 300.0, 2'000.0), 1.e-20);
  }
  Check(lowBandRatios[0] > lowBandRatios[1] &&
            lowBandRatios[1] > lowBandRatios[2],
        "Low cut must monotonically attenuate bass relative to the midband");

  std::array<double, 3> highBandRatios{};
  const std::array<float, 3> highCuts{0.f, tfdsp::reverb_defaults::HighCut,
                                      1.f};
  for (std::size_t index = 0; index < highCuts.size(); ++index) {
    tfdsp::RoomReverbControls controls;
    controls.highCut = highCuts[index];
    const auto response = Render(controls, 48'000);
    highBandRatios[index] =
        BandEnergy(response, 8'000.0, 18'000.0) /
        std::max(BandEnergy(response, 500.0, 4'000.0), 1.e-20);
  }
  Check(highBandRatios[0] < highBandRatios[1] &&
            highBandRatios[1] < highBandRatios[2],
        "High cut must monotonically restore treble as its cutoff rises");
}

void TestMixAndOutputLevelLaws() {
  constexpr std::array<float, 2> direct{{2.f, -0.5f}};
  constexpr std::array<float, 2> wet{{3.f, -1.f}};
  const auto dryOnly = tfdsp::MixReverbOutput(direct, wet, 0.f, 0.f);
  const auto wetOnly = tfdsp::MixReverbOutput(direct, wet, 1.f, 0.f);
  Check(std::abs(dryOnly[0] - direct[0]) < 1.e-6f &&
            std::abs(dryOnly[1] - direct[1]) < 1.e-6f,
        "Mix at zero must be exactly dry");
  Check(std::abs(wetOnly[0] - wet[0]) < 1.e-6f &&
            std::abs(wetOnly[1] - wet[1]) < 1.e-6f,
        "Mix at maximum must be exactly wet");

  const auto baseline =
      tfdsp::MixReverbOutput(direct, wet, tfdsp::reverb_defaults::Mix, 0.f);
  const auto lower =
      tfdsp::MixReverbOutput(direct, wet, tfdsp::reverb_defaults::Mix, -6.f);
  const auto higher =
      tfdsp::MixReverbOutput(direct, wet, tfdsp::reverb_defaults::Mix, 6.f);
  const float sixDb = std::pow(10.f, 6.f / 20.f);
  for (std::size_t channel = 0; channel < 2; ++channel) {
    Check(std::abs(lower[channel] * sixDb - baseline[channel]) < 1.e-5f,
          "Output level -6 dB must apply the exact gain law");
    Check(std::abs(higher[channel] / sixDb - baseline[channel]) < 1.e-5f,
          "Output level +6 dB must apply the exact gain law");
  }
}

void TestRoomReverbAtEverySupportedRate() {
  for (const double sampleRate : {44'100.0, 48'000.0, 88'200.0,
                                  96'000.0, 192'000.0}) {
    tfdsp::RoomReverb reverb;
    reverb.SetSampleRate(sampleRate);
    tfdsp::RoomReverbControls controls;
    tfdsp::RoomReverb::SourcePositions positions{};
    positions[0] = tfdsp::reverb_defaults::Source;
    bool finite = true;
    const auto count = static_cast<std::size_t>(sampleRate * .02);
    for (std::size_t sample = 0; sample < count; ++sample) {
      tfdsp::RoomReverb::InputFrame input{};
      input[0] = sample == 0 ? 1.f : 0.f;
      const auto output = reverb.Process(input, positions, 1, controls);
      finite = finite && std::isfinite(output.direct[0]) &&
               std::isfinite(output.direct[1]) &&
               std::isfinite(output.wet[0]) && std::isfinite(output.wet[1]);
    }
    Check(finite, "room reverb remains finite at every supported sample rate");
  }
}

void TestRoomReverbFlushesSubnormals() {
  tfdsp::RoomReverb reverb;
  tfdsp::RoomReverbControls controls;
  tfdsp::RoomReverb::SourcePositions positions{};
  positions[0] = tfdsp::reverb_defaults::Source;
  bool exactSilence = true;
  for (std::size_t sample = 0; sample < 8192; ++sample) {
    tfdsp::RoomReverb::InputFrame input{};
    if (sample == 0)
      input[0] = std::numeric_limits<float>::denorm_min();
    const auto output = reverb.Process(input, positions, 1, controls);
    exactSilence = exactSilence && output.direct[0] == 0.f &&
                   output.direct[1] == 0.f && output.wet[0] == 0.f &&
                   output.wet[1] == 0.f;
  }
  Check(exactSilence,
        "complete room-reverb direct, ER and tail paths flush subnormals");
}

} // namespace

int main() {
  TestCanonicalFactoryDefaults();
  TestProgressiveDefaultSourceSpread();
  TestAcousticPresetDefinitions();
  TestRoomGeometryDerivesHeightFromSizeAndUsesAspectAndListener();
  TestRealtimeSceneConstructionSanitizesNonFiniteControls();
  TestCombinedEarlyAndLateResponse();
  TestEarlyTailIsolationAndPreDelay();
  TestGeometryControlsChangeTheResponseIndependently();
  TestSizeSweepScalesAllDimensionsWithoutBecomingALevelControl();
  TestSourceAndListenerPositionSweepsChangeTheResponse();
  TestPositionedDirectSound();
  TestResetPreservesAcceptedHandoffMetadata();
  TestPreDelayAutomationCrossfadesWithoutAnEnergyHole();
  TestPreDelaySweepTranslatesTheWetImpulse();
  TestDecayAndDampingSweepsHaveExpectedSpectralEffects();
  TestDiffusionModulationAndShimmerSweeps();
  TestWidthLevelAndFilterSweeps();
  TestMixAndOutputLevelLaws();
  TestRoomReverbAtEverySupportedRate();
  TestRoomReverbFlushesSubnormals();
  std::cout << "Room reverb tests passed\n";
}
