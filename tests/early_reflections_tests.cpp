#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <new>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "tfdsp/early_reflections.hpp"
#include "tfdsp/early_reflections_worker.hpp"
#include "tfdsp/scene_pack.hpp"

namespace {
std::atomic<std::size_t> trackedAllocations{};
thread_local bool trackAllocations{};
} // namespace

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
int failures = 0;

void Check(const bool condition, const std::string &message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << '\n';
  }
}

void CheckNear(const double actual, const double expected,
               const double tolerance, const std::string &message) {
  Check(std::abs(actual - expected) <= tolerance,
        message + " (actual " + std::to_string(actual) + ", expected " +
            std::to_string(expected) + ")");
}

template <typename Function>
void CheckThrows(Function &&function, const std::string &message) {
  try {
    function();
  } catch (const std::exception &) {
    return;
  }
  Check(false, message);
}

tfdsp::EarlyReflectionRoom TestRoom() {
  return {{6.0, 8.0, 3.0}, {3.0, 5.0, 1.35}};
}

tfdsp::EarlyReflectionMaterials UniformMaterials(const double gain = 0.8) {
  tfdsp::EarlyReflectionMaterials materials;
  for (auto &surface : materials.reflectionAmplitudes)
    surface.fill(gain);
  materials.airAbsorptionDbPerMetre.fill(0.0);
  return materials;
}

tfdsp::EarlyReflectionConfig TestConfig(const double responseSeconds = 0.08) {
  tfdsp::EarlyReflectionConfig config;
  config.sampleRate = 48000.0;
  config.responseTimeSeconds = responseSeconds;
  config.materialFilterTailSeconds = 0.005;
  return config;
}

const tfdsp::EarlyReflectionPath *
FindPath(const std::vector<tfdsp::EarlyReflectionPath> &paths,
         const std::size_t source, const std::array<int, 3> index) {
  const auto found = std::find_if(
      paths.begin(), paths.end(), [&](const tfdsp::EarlyReflectionPath &path) {
        return path.sourceIndex == source && path.imageIndex == index;
      });
  return found == paths.end() ? nullptr : &*found;
}

void TestDefaultsAndLimits() {
  Check(tfdsp::EarlyReflectionMaximumSources == 8,
        "the public ER source limit is eight");
  Check(tfdsp::EarlyReflectionImageCountThroughOrder(1) == 6,
        "order-one image count is six");
  Check(tfdsp::EarlyReflectionImageCountThroughOrder(2) == 24,
        "order-two cumulative image count is 24");
  Check(tfdsp::EarlyReflectionImageCountThroughOrder(4) == 128,
        "order-four cumulative image count is 128");

  const auto room = TestRoom();
  const auto sources = tfdsp::MakeDefaultEarlyReflectionSources(room, 1, 0.5);
  Check(sources.size() == 1, "the default scene contains one source");
  CheckNear(sources[0].positionMetres[0], 0.5 * room.dimensionsMetres[0],
            1.0e-12, "the single default source is horizontally centred");
  CheckThrows([&] { tfdsp::MakeDefaultEarlyReflectionSources(room, 9, 0.5); },
              "nine sources are rejected");
  Check(tfdsp::MakeDefaultEarlyReflectionSources(room, 8, 0.5).size() == 8,
        "eight sources are accepted");

  auto undersizedRoom = room;
  undersizedRoom.dimensionsMetres[0] = 0.999;
  CheckThrows([&] { undersizedRoom.Validate(); },
              "every room dimension must be at least one metre");
}

void TestSubSamplePathUsesCausalInterpolation() {
  auto config = TestConfig(0.021);
  config.sampleRate = 50.0;
  config.crossoverHz = {1.0, 5.0, 10.0};
  config.minimumHandoffSeconds = 0.020;
  config.handoffOverlapSeconds = 0.001;
  config.analysisWindowSeconds = 0.010;
  config.analysisHopSeconds = 0.001;
  config.analysisStableSeconds = 0.001;
  config.materialFilterTailSeconds = 0.020;
  const auto room = TestRoom();
  const std::vector<tfdsp::EarlyReflectionSource> sources{
      tfdsp::MakeEarlyReflectionSource(room, {0.25, 0.30, 0.40})};
  const auto materials = UniformMaterials();
  const auto paths =
      tfdsp::EnumerateEarlyReflectionPaths(config, room, sources, materials);
  Check(!paths.empty() &&
            paths.front().propagationSeconds * config.sampleRate < 1.0,
        "the regression scene contains a valid sub-sample reflection");
  const auto response = tfdsp::BuildEarlyReflectionImpulseResponse(
      config, room, sources, materials);
  Check(response.Size() > 0, "a sub-sample reflection builds a causal FIR");
  for (const auto &kernel : response.kernels)
    for (const double sample : kernel)
      Check(std::isfinite(sample),
            "causal sub-sample interpolation remains finite");
}

void TestImageSourceGeometryAndMaterials() {
  const auto config = TestConfig();
  const auto room = TestRoom();
  const std::vector<tfdsp::EarlyReflectionSource> sources{
      tfdsp::MakeEarlyReflectionSource(room, {0.25, 0.30, 0.40})};
  const auto materials = UniformMaterials(0.8);
  const auto paths =
      tfdsp::EnumerateEarlyReflectionPaths(config, room, sources, materials);
  const auto mixing = tfdsp::PredictEarlyReflectionMixing(config, room);
  Check(!paths.empty(), "the image-source model produces reflected paths");
  Check(FindPath(paths, 0, {0, 0, 0}) == nullptr,
        "the direct path is excluded from the ER response");

  const auto *lowerX = FindPath(paths, 0, {-1, 0, 0});
  const auto *upperX = FindPath(paths, 0, {1, 0, 0});
  const auto *lowerY = FindPath(paths, 0, {0, -1, 0});
  const auto *lowerZ = FindPath(paths, 0, {0, 0, -1});
  const auto *upperZ = FindPath(paths, 0, {0, 0, 1});
  Check(lowerX != nullptr && upperX != nullptr,
        "both first-order x-wall images are present");
  if (lowerX != nullptr) {
    CheckNear(lowerX->imagePositionMetres[0], -sources[0].positionMetres[0],
              1.0e-12, "the lower-wall image coordinate is analytic");
    Check(lowerX->surfaceCounts[2] == 1 && lowerX->reflectionOrder == 1,
          "the x-wall bounce is counted once");
    const double expectedGain = 0.8 / lowerX->distanceMetres;
    for (const double gain : lowerX->bandGains)
      CheckNear(gain, expectedGain, 1.0e-12,
                "uniform wall material gives the analytic path gain");
    CheckNear(lowerX->outputGains[0] * lowerX->outputGains[0] +
                  lowerX->outputGains[1] * lowerX->outputGains[1],
              1.0, 1.0e-12, "stereo panning is equal-power");
  }
  if (upperX != nullptr)
    CheckNear(upperX->imagePositionMetres[0],
              2.0 * room.dimensionsMetres[0] - sources[0].positionMetres[0],
              1.0e-12, "the upper-wall image coordinate is analytic");
  Check(lowerY != nullptr && lowerY->surfaceCounts[3] == 1,
        "front/rear-wall bounces use their material group");
  Check(lowerZ != nullptr && lowerZ->surfaceCounts[0] == 1,
        "floor bounces use the floor material group");
  Check(upperZ != nullptr && upperZ->surfaceCounts[1] == 1,
        "ceiling bounces use the ceiling material group");

  for (const auto &path : paths) {
    Check(path.excessDelaySeconds <= mixing.generationHorizonSeconds + 1.0e-12,
          "every enumerated path lies inside its direct-relative time window");
    Check(path.reflectionOrder ==
              static_cast<std::size_t>(std::abs(path.imageIndex[0]) +
                                       std::abs(path.imageIndex[1]) +
                                       std::abs(path.imageIndex[2])),
          "reflection order matches the image lattice index");
  }

  const auto brightPaths = tfdsp::EnumerateEarlyReflectionPaths(
      config, room, sources, tfdsp::MakeEarlyReflectionMaterials(0.0));
  const auto darkPaths = tfdsp::EnumerateEarlyReflectionPaths(
      config, room, sources, tfdsp::MakeEarlyReflectionMaterials(1.0));
  Check(brightPaths.size() == darkPaths.size(),
        "damping changes path spectra while preserving geometry");
  if (!brightPaths.empty() && brightPaths.size() == darkPaths.size()) {
    const auto &bright = brightPaths.front().bandGains;
    const auto &dark = darkPaths.front().bandGains;
    Check(dark[3] < bright[3],
          "the damped material reduces high-frequency reflection gain");
    Check(dark[3] / dark[0] < bright[3] / bright[0],
          "damping steepens the reflected high-to-low spectral slope");
  }
}

void TestAdaptiveMixingPredictionAndRelativeTime() {
  auto config = TestConfig(0.15);
  const auto room = tfdsp::MakeEarlyReflectionRoom(0.45);
  const auto neutral = tfdsp::PredictEarlyReflectionMixing(config, room, 0.5);
  const double length = room.dimensionsMetres[0];
  const double width = room.dimensionsMetres[1];
  const double height = room.dimensionsMetres[2];
  const double expectedVolume = length * width * height;
  const double expectedSurface =
      2.0 * (length * width + length * height + width * height);
  CheckNear(neutral.roomVolumeCubicMetres, expectedVolume, 1.0e-12,
            "mixing prediction uses the analytic room volume");
  CheckNear(neutral.roomSurfaceSquareMetres, expectedSurface, 1.0e-12,
            "mixing prediction uses the analytic room surface area");
  CheckNear(neutral.averageSeconds,
            0.001 * (20.0 * expectedVolume / expectedSurface + 12.0) *
                neutral.diffusionMultiplier,
            1.0e-12, "T50 follows the V/S perceptual predictor");
  const auto low = tfdsp::PredictEarlyReflectionMixing(config, room, 0.0);
  const auto high = tfdsp::PredictEarlyReflectionMixing(config, room, 1.0);
  Check(low.averageSeconds > neutral.averageSeconds &&
            neutral.averageSeconds > high.averageSeconds,
        "higher Diffusion produces an earlier mixing prediction");
  Check(neutral.generationHorizonSeconds >= neutral.conservativeSeconds,
        "generation includes the analysis safety region");

  const std::vector<tfdsp::EarlyReflectionSource> sources{
      tfdsp::MakeEarlyReflectionSource(room, {0.5, 0.60, 0.45}),
      tfdsp::MakeEarlyReflectionSource(room, {0.2, 0.10, 0.70})};
  const auto paths = tfdsp::EnumerateEarlyReflectionPaths(
      config, room, sources, UniformMaterials(), 0.5);
  std::array<bool, 2> sawPath{};
  for (const auto &path : paths) {
    const double direct =
        tfdsp::EarlyReflectionDirectDistance(room, sources[path.sourceIndex]);
    CheckNear(path.excessDelaySeconds,
              (path.distanceMetres - direct) / config.speedOfSound, 1.0e-12,
              "path horizon is measured relative to that source's direct path");
    Check(path.excessDelaySeconds <= neutral.generationHorizonSeconds + 1.0e-12,
          "direct-relative path lies inside the adaptive horizon");
    sawPath[path.sourceIndex] = true;
  }
  Check(sawPath[0] && sawPath[1],
        "near and far sources both receive adaptive ER paths");
}

void TestEchoDensityStatisticsAndRefinement() {
  constexpr std::size_t count = 4096;
  std::vector<double> left(count);
  std::vector<double> right(count);
  std::mt19937 generator(0x5eed);
  std::normal_distribution<double> normal;
  for (std::size_t sample = 0; sample < count; ++sample) {
    left[sample] = normal(generator);
    right[sample] = normal(generator);
  }
  const auto diffuse =
      tfdsp::CalculateEarlyReflectionDensity(left, right, 0, count);
  Check(diffuse.valid, "Gaussian late energy produces valid density metrics");
  Check(std::abs(diffuse.normalizedEchoDensity - 1.0) < 0.08,
        "Gaussian late energy has normalized echo density near one");
  Check(std::abs(diffuse.excessKurtosis) < 0.15,
        "Gaussian late energy has excess kurtosis near zero");

  std::fill(left.begin(), left.end(), 0.0);
  std::fill(right.begin(), right.end(), 0.0);
  left[count / 2] = 1.0;
  const auto sparse =
      tfdsp::CalculateEarlyReflectionDensity(left, right, 0, count);
  Check(sparse.valid && sparse.normalizedEchoDensity < 0.1 &&
            sparse.excessKurtosis > 100.0,
        "a sparse echo is rejected as a diffuse field");

  tfdsp::EarlyReflectionImpulseResponse response;
  response.sampleRate = 48000.0;
  response.sourceCount = 1;
  response.kernels = {std::vector<double>(count), std::vector<double>(count)};
  generator.seed(1234);
  for (std::size_t sample = 0; sample < count; ++sample) {
    response.kernels[0][sample] = normal(generator);
    response.kernels[1][sample] = normal(generator);
  }
  response.mixingPrediction.generationHorizonSeconds = 0.075;
  response.sourceHandoffs.resize(1);
  response.sourceHandoffs[0].directPropagationSeconds = 0.005;
  response.sourceHandoffs[0].predictedStartSeconds = 0.030;
  response.sourceHandoffs[0].predictedEndSeconds = 0.050;
  response.sourceHandoffs[0].finalStartSeconds = 0.030;
  response.sourceHandoffs[0].finalEndSeconds = 0.050;
  auto analysisConfig = TestConfig(0.08);
  tfdsp::RefineEarlyReflectionHandoffs(analysisConfig, response);
  Check(response.sourceHandoffs[0].detectedFromResponse,
        "stable Gaussian response statistics refine the mixing time");
  Check(response.sourceHandoffs[0].finalEndSeconds >=
            response.sourceHandoffs[0].finalStartSeconds +
                analysisConfig.handoffOverlapSeconds - 1.0e-12,
        "refined handoff retains the minimum overlap");
}

void TestAllPathsAndMultipleSources() {
  const auto config = TestConfig(0.15);
  const auto room = tfdsp::MakeEarlyReflectionRoom(0.45);
  const auto materials = tfdsp::MakeEarlyReflectionMaterials(0.5);
  const auto one = tfdsp::MakeDefaultEarlyReflectionSources(room, 1, 0.5);
  const auto two = std::vector<tfdsp::EarlyReflectionSource>{one[0], one[0]};
  const auto onePaths =
      tfdsp::EnumerateEarlyReflectionPaths(config, room, one, materials);
  const auto twoPaths =
      tfdsp::EnumerateEarlyReflectionPaths(config, room, two, materials);
  Check(onePaths.size() > 128,
        "the time-bounded generator is not capped at order four or 128 paths");
  Check(twoPaths.size() == 2 * onePaths.size(),
        "each independent source gets a complete image-source response");

  const auto response = tfdsp::BuildEarlyReflectionImpulseResponse(
      config, room, two, materials, 128);
  Check(response.sourceCount == 2 && response.kernels.size() == 4,
        "two mono sources produce a 2x2 source-to-output FIR matrix");
  for (std::size_t output = 0; output < 2; ++output) {
    const auto &first = response.kernels[response.KernelIndex(output, 0)];
    const auto &second = response.kernels[response.KernelIndex(output, 1)];
    Check(first == second,
          "coincident sources generate identical independent FIR kernels");
  }
}

void TestGeneratedFirMatchesAnalyticTapTrain() {
  const auto config = TestConfig(0.045);
  const auto room = TestRoom();
  const std::vector<tfdsp::EarlyReflectionSource> sources{
      tfdsp::MakeEarlyReflectionSource(room, {0.32, 0.27, 0.43})};
  const auto materials = UniformMaterials(0.73);
  constexpr std::size_t convolutionLatency = 128;
  const auto paths =
      tfdsp::EnumerateEarlyReflectionPaths(config, room, sources, materials);
  const auto response = tfdsp::BuildEarlyReflectionImpulseResponse(
      config, room, sources, materials, convolutionLatency);
  Check(response.appliedLatencyCompensationSamples == convolutionLatency,
        "normal ER timing fully absorbs the convolution latency");
  CheckNear(response.ResidualLatencySeconds(), 0.0, 0.0,
            "fully compensated convolution has no residual timing offset");

  std::vector<std::vector<double>> expected(
      2, std::vector<double>(response.Size()));
  for (const auto &path : paths) {
    const double delay =
        path.propagationSeconds * config.sampleRate -
        static_cast<double>(response.appliedLatencyCompensationSamples);
    const auto integerDelay = static_cast<std::size_t>(std::floor(delay));
    const auto coefficients =
        tfdsp::EarlyReflectionLagrange4(
            delay - static_cast<double>(integerDelay));
    constexpr std::array<int, 4> offsets{-1, 0, 1, 2};
    for (std::size_t output = 0; output < 2; ++output) {
      const double gain = path.outputGains[output] * path.bandGains[0];
      for (std::size_t tap = 0; tap < 4; ++tap) {
        const auto index = static_cast<std::size_t>(
            static_cast<std::ptrdiff_t>(integerDelay) + offsets[tap]);
        if (index < expected[output].size())
          expected[output][index] += gain * coefficients[tap];
      }
    }
  }
  const auto &handoff = response.sourceHandoffs[0];
  for (auto &kernel : expected)
    for (std::size_t sample = 0; sample < kernel.size(); ++sample) {
      const double relativeSeconds =
          (static_cast<double>(sample) +
           static_cast<double>(response.appliedLatencyCompensationSamples)) /
              config.sampleRate -
          handoff.directPropagationSeconds;
      if (relativeSeconds >= handoff.finalEndSeconds)
        kernel[sample] = 0.0;
      else if (relativeSeconds > handoff.finalStartSeconds) {
        const double position =
            (relativeSeconds - handoff.finalStartSeconds) /
            (handoff.finalEndSeconds - handoff.finalStartSeconds);
        kernel[sample] *= std::cos(0.5 * tfdsp::EarlyReflectionPi *
                                   tfdsp::EarlyReflectionSmoothstep(position));
      }
    }
  for (std::size_t output = 0; output < 2; ++output)
    for (std::size_t sample = 0; sample < response.Size(); ++sample)
      CheckNear(
          response.kernels[response.KernelIndex(output, 0)][sample],
          expected[output][sample], 2.0e-12,
          "equal band gains reconstruct the exact fractional-delay tap train");
}

double MagnitudeAt(const std::vector<double> &impulse, const double frequency,
                   const double sampleRate) {
  double real = 0.0;
  double imaginary = 0.0;
  for (std::size_t sample = 0; sample < impulse.size(); ++sample) {
    const double phase = -2.0 * tfdsp::EarlyReflectionPi * frequency *
                         static_cast<double>(sample) / sampleRate;
    real += impulse[sample] * std::cos(phase);
    imaginary += impulse[sample] * std::sin(phase);
  }
  return std::hypot(real, imaginary);
}

void TestUnequalMaterialBandsProduceSpectralTilt() {
  const auto config = TestConfig(0.060);
  const auto room = TestRoom();
  const std::vector<tfdsp::EarlyReflectionSource> sources{
      tfdsp::MakeEarlyReflectionSource(room, {0.25, 0.30, 0.40})};
  auto materials = UniformMaterials(1.0);
  materials.reflectionAmplitudes[2] = {1.0, 0.72, 0.30, 0.05};
  const auto paths =
      tfdsp::EnumerateEarlyReflectionPaths(config, room, sources, materials);
  const auto *sideWall = FindPath(paths, 0, {-1, 0, 0});
  Check(sideWall != nullptr,
        "the spectral material test has a first-order side-wall path");
  if (sideWall == nullptr)
    return;
  const std::vector<tfdsp::EarlyReflectionPath> onePath{*sideWall};
  const auto response = tfdsp::BuildEarlyReflectionImpulseResponse(
      config, room, sources, materials, 128, 0.5, &onePath);
  const auto &left = response.kernels[response.KernelIndex(0, 0)];
  const auto &right = response.kernels[response.KernelIndex(1, 0)];
  const double low = std::hypot(MagnitudeAt(left, 100.0, config.sampleRate),
                                MagnitudeAt(right, 100.0, config.sampleRate));
  const double high = std::hypot(MagnitudeAt(left, 10000.0, config.sampleRate),
                                 MagnitudeAt(right, 10000.0,
                                             config.sampleRate));
  Check(low > 0.0 && high < 0.4 * low,
        "unequal material bands produce the requested high-frequency damping");
}

tfdsp::EarlyReflectionImpulseResponse
SyntheticResponse(const std::size_t sources, const std::size_t size,
                  const std::vector<std::vector<double>> &kernels) {
  tfdsp::EarlyReflectionImpulseResponse response;
  response.sampleRate = 48000.0;
  response.sourceCount = sources;
  response.imagePathCount = 1;
  response.analysisPathCount = 1;
  response.requestedConvolutionLatencySamples = 16;
  response.appliedLatencyCompensationSamples = 16;
  response.mixingPrediction.generationHorizonSeconds = 0.001;
  response.sourceHandoffs.resize(sources);
  for (auto &handoff : response.sourceHandoffs) {
    handoff.finalStartSeconds = 0.0005;
    handoff.finalEndSeconds = 0.001;
    handoff.imagePathCount = 1;
    handoff.analysisPathCount = 1;
  }
  response.kernels = kernels;
  for (auto &kernel : response.kernels)
    kernel.resize(size);
  response.Validate();
  return response;
}

void TestPartitionedConvolutionAndSuperposition() {
  constexpr std::size_t length = 91;
  std::vector<std::vector<double>> kernels(4, std::vector<double>(length, 0.0));
  kernels[0][0] = 0.7;
  kernels[0][17] = -0.2;
  kernels[0][63] = 0.11;
  kernels[1][3] = -0.4;
  kernels[1][31] = 0.25;
  kernels[2][1] = 0.3;
  kernels[2][48] = -0.17;
  kernels[3][0] = -0.6;
  kernels[3][20] = 0.09;
  const auto response = SyntheticResponse(2, length, kernels);

  tfdsp::EarlyReflectionConvolver<double, 16, 8> convolver;
  convolver.Prepare(48000.0, 128);
  convolver.SetImpulseResponse(response);
  constexpr std::size_t inputLength = 160;
  constexpr std::size_t renderLength = inputLength + length + 32;
  std::array<std::vector<double>, 2> input{std::vector<double>(inputLength),
                                           std::vector<double>(inputLength)};
  for (std::size_t sample = 0; sample < inputLength; ++sample) {
    input[0][sample] = std::sin(0.071 * static_cast<double>(sample));
    input[1][sample] = 0.4 * std::cos(0.113 * static_cast<double>(sample));
  }
  std::array<std::vector<double>, 2> expected{
      std::vector<double>(renderLength), std::vector<double>(renderLength)};
  for (std::size_t output = 0; output < 2; ++output)
    for (std::size_t source = 0; source < 2; ++source)
      for (std::size_t sample = 0; sample < inputLength; ++sample)
        for (std::size_t tap = 0; tap < length; ++tap)
          if (sample + tap + 16 < renderLength)
            expected[output][sample + tap + 16] +=
                input[source][sample] * kernels[output * 2 + source][tap];

  for (std::size_t sample = 0; sample < renderLength; ++sample) {
    tfdsp::EarlyReflectionConvolver<double, 16, 8>::InputFrame frame{};
    if (sample < inputLength) {
      frame[0] = input[0][sample];
      frame[1] = input[1][sample];
    }
    const auto output = convolver.Process(frame, 2);
    CheckNear(
        output[0], expected[0][sample], 2.0e-10,
        "partitioned left output matches direct multi-source convolution");
    CheckNear(
        output[1], expected[1][sample], 2.0e-10,
        "partitioned right output matches direct multi-source convolution");
  }
}

void TestPhysicalTimingWithCompensatedLatency() {
  const auto config = TestConfig(0.04);
  const auto room = TestRoom();
  const std::vector<tfdsp::EarlyReflectionSource> sources{
      tfdsp::MakeEarlyReflectionSource(room, {0.4, 0.3, 0.4})};
  const auto materials = UniformMaterials(0.8);
  const auto paths =
      tfdsp::EnumerateEarlyReflectionPaths(config, room, sources, materials);
  const auto response = tfdsp::BuildEarlyReflectionImpulseResponse(
      config, room, sources, materials, 16);
  tfdsp::EarlyReflectionConvolver<double, 16, 8> convolver;
  convolver.Prepare(config.sampleRate,
                    tfdsp::MaximumEarlyReflectionImpulseSamples(config, room));
  convolver.SetImpulseResponse(response);

  const std::size_t renderLength = response.Size() + 32;
  std::array<std::vector<double>, 2> rendered{
      std::vector<double>(renderLength), std::vector<double>(renderLength)};
  for (std::size_t sample = 0; sample < renderLength; ++sample) {
    tfdsp::EarlyReflectionConvolver<double, 16, 8>::InputFrame input{};
    input[0] = sample == 0 ? 1.0 : 0.0;
    const auto output = convolver.Process(input, 1);
    rendered[0][sample] = output[0];
    rendered[1][sample] = output[1];
  }
  for (std::size_t output = 0; output < 2; ++output)
    for (std::size_t tap = 0; tap < response.Size(); ++tap)
      CheckNear(
          rendered[output][tap + 16],
          response.kernels[response.KernelIndex(output, 0)][tap], 2.0e-10,
          "the runtime convolution latency equals the compensated FIR shift");

  const double firstPhysicalArrival =
      paths.front().propagationSeconds * config.sampleRate;
  const double firstRenderedTap =
      firstPhysicalArrival -
      static_cast<double>(response.appliedLatencyCompensationSamples) + 16.0;
  CheckNear(firstRenderedTap, firstPhysicalArrival, 1.0e-12,
            "physical propagation is preserved after latency compensation");
}

void TestRealtimeSafetyAndValidation() {
  std::vector<std::vector<double>> firstKernels(2,
                                                std::vector<double>(48, 0.0));
  firstKernels[0][0] = 1.0;
  firstKernels[1][0] = 1.0;
  auto secondKernels = firstKernels;
  secondKernels[0][7] = 0.5;
  secondKernels[1][11] = -0.25;
  const auto first = SyntheticResponse(1, 48, firstKernels);
  const auto second = SyntheticResponse(1, 48, secondKernels);
  tfdsp::EarlyReflectionConvolver<float, 16, 8> convolver;
  convolver.Prepare(48000.0, 64);
  convolver.SetImpulseResponse(first);
  convolver.TransitionToImpulseResponse(second, 0.001);

  trackedAllocations.store(0, std::memory_order_relaxed);
  trackAllocations = true;
  float checksum = 0.0f;
  for (std::size_t sample = 0; sample < 512; ++sample) {
    tfdsp::EarlyReflectionConvolver<float, 16, 8>::InputFrame input{};
    input[0] = sample == 0 ? 1.0f : 0.01f;
    const auto output = convolver.Process(input, 1);
    checksum += output[0] + output[1];
  }
  trackAllocations = false;
  Check(trackedAllocations.load(std::memory_order_relaxed) == 0,
        "audio-thread processing and FIR crossfading allocate no memory");
  Check(std::isfinite(checksum), "realtime output remains finite");
  Check(!convolver.IsTransitioning(), "the scheduled FIR transition completes");

  auto wrongLatency = first;
  wrongLatency.requestedConvolutionLatencySamples = 128;
  CheckThrows([&] { convolver.SetImpulseResponse(wrongLatency); },
              "an FIR generated for the wrong convolution latency is rejected");
  CheckThrows(
      [&] {
        auto tooLong = first;
        for (auto &kernel : tooLong.kernels)
          kernel.resize(65);
        convolver.SetImpulseResponse(tooLong);
      },
      "an FIR larger than the prepared capacity is rejected");
}

void TestNumericalImpulseResponseTransition() {
  constexpr std::size_t partition = 16;
  std::vector<std::vector<double>> oldKernels(
      2, std::vector<double>(partition, 0.0));
  std::vector<std::vector<double>> newKernels = oldKernels;
  oldKernels[0][0] = oldKernels[1][0] = 1.0;
  newKernels[0][0] = newKernels[1][0] = 3.0;
  const auto oldResponse = SyntheticResponse(1, partition, oldKernels);
  const auto newResponse = SyntheticResponse(1, partition, newKernels);
  tfdsp::EarlyReflectionConvolver<double, partition, 8> convolver;
  convolver.Prepare(48000.0, partition);
  convolver.SetImpulseResponse(oldResponse);
  for (std::size_t sample = 0; sample < 4 * partition; ++sample) {
    tfdsp::EarlyReflectionConvolver<double, partition, 8>::InputFrame input{};
    input[0] = 1.0;
    convolver.Process(input, 1);
  }
  constexpr std::size_t transitionSamples = 32;
  convolver.TransitionToImpulseResponse(
      newResponse, static_cast<double>(transitionSamples) / 48000.0);
  for (std::size_t sample = 0; sample < partition + transitionSamples + 8;
       ++sample) {
    tfdsp::EarlyReflectionConvolver<double, partition, 8>::InputFrame input{};
    input[0] = 1.0;
    const auto output = convolver.Process(input, 1);
    double expected = 1.0;
    if (sample >= partition) {
      const std::size_t transitionSample = sample - partition;
      if (transitionSample >= transitionSamples)
        expected = 3.0;
      else
        expected = 1.0 + 2.0 * tfdsp::EarlyReflectionSmoothstep(
                                   static_cast<double>(transitionSample) /
                                   static_cast<double>(transitionSamples - 1));
    }
    CheckNear(output[0], expected, 2.0e-10,
              "the FIR transition follows its expected smooth waveform");
    CheckNear(output[1], expected, 2.0e-10,
              "both FIR outputs use the same transition waveform");
  }
}

void TestInactiveSourceTailAndHistoryFlush() {
  constexpr std::size_t partition = 16;
  constexpr std::size_t impulseSize = 64;
  std::vector<std::vector<double>> kernels(
      4, std::vector<double>(impulseSize, 0.0));
  kernels[1][47] = 0.75;
  kernels[3][47] = -0.5;
  const auto response = SyntheticResponse(2, impulseSize, kernels);
  tfdsp::EarlyReflectionConvolver<double, partition, 8> convolver;
  convolver.Prepare(48000.0, impulseSize);
  convolver.SetImpulseResponse(response);

  constexpr std::size_t firstRenderLength = 160;
  for (std::size_t sample = 0; sample < firstRenderLength; ++sample) {
    tfdsp::EarlyReflectionConvolver<double, partition, 8>::InputFrame input{};
    input[1] = sample == 0 ? 1.0 : 0.0;
    const auto output = convolver.Process(input, sample == 0 ? 2 : 1);
    const double expectedLeft = sample == partition + 47 ? 0.75 : 0.0;
    const double expectedRight = sample == partition + 47 ? -0.5 : 0.0;
    CheckNear(output[0], expectedLeft, 2.0e-10,
              "a removed source retains its remaining left FIR tail");
    CheckNear(output[1], expectedRight, 2.0e-10,
              "a removed source retains its remaining right FIR tail");
  }

  for (std::size_t sample = 0; sample < 2 * impulseSize; ++sample) {
    tfdsp::EarlyReflectionConvolver<double, partition, 8>::InputFrame input{};
    const auto output = convolver.Process(input, 2);
    CheckNear(output[0], 0.0, 2.0e-10,
              "reactivating a flushed source has no stale left history");
    CheckNear(output[1], 0.0, 2.0e-10,
              "reactivating a flushed source has no stale right history");
  }
}

void TestScenePackCompactionAndChaining() {
  tfdsp::ScenePackInput first;
  first.localConnected = {true, false, true, false};
  first.local[0] = {1.0f, 2.0f, 3.0f, 4.0f};
  first.local[2] = {3.0f, 6.0f, 7.0f, 8.0f};
  const auto firstOutput = tfdsp::PackSceneSources(first);
  Check(firstOutput.sourceCount == 2 && firstOutput.sources[0].audio == 1.0f &&
            firstOutput.sources[1].audio == 3.0f,
        "Scene Pack compacts connected local lanes in lane order");

  tfdsp::ScenePackInput chained;
  chained.bus = firstOutput.sources;
  chained.busCount = firstOutput.sourceCount;
  chained.localConnected.fill(true);
  for (std::size_t lane = 0; lane < chained.local.size(); ++lane)
    chained.local[lane] = {static_cast<float>(10 + lane), 5.0f, 5.0f, 5.0f};
  const auto chainedOutput = tfdsp::PackSceneSources(chained);
  Check(chainedOutput.sourceCount == 6 &&
            chainedOutput.sources[0].audio == 1.0f &&
            chainedOutput.sources[2].audio == 10.0f,
        "a chained Scene Pack preserves bus channels before local lanes");

  chained.busCount = 7;
  for (std::size_t source = 0; source < 7; ++source)
    chained.bus[source].audio = static_cast<float>(source);
  const auto limited = tfdsp::PackSceneSources(chained);
  Check(limited.sourceCount == 8 && limited.sources[7].audio == 10.0f,
        "Scene Pack appends only enough local lanes to reach eight sources");
}

void TestRateLimitedCoalescingWorkerAndPreparedBankHandoff() {
  using Convolver = tfdsp::EarlyReflectionConvolver<float, 16, 8>;
  auto request = tfdsp::EarlyReflectionBuildRequest{};
  request.config = TestConfig(0.08);
  request.room = TestRoom();
  request.SetSources(
      tfdsp::MakeDefaultEarlyReflectionSources(request.room, 1, 0.5));
  request.materials = UniformMaterials();
  request.convolutionLatencySamples = Convolver::LatencySamples;
  request.transitionSeconds = 0.001;

  Convolver convolver;
  convolver.Prepare(request.config.sampleRate,
                    tfdsp::MaximumEarlyReflectionImpulseSamples(request.config,
                                                                request.room));
  {
    tfdsp::EarlyReflectionWorker worker(
        20.0, [&](const tfdsp::EarlyReflectionImpulseResponse &response,
                  const double transitionSeconds) {
          return convolver.PrepareAndQueueImpulseResponse(response,
                                                          transitionSeconds);
        });
    const auto sequence = worker.Submit(request);
    auto result = worker.WaitForLatestResult(std::chrono::milliseconds(1000));
    Check(result.has_value() && result->Succeeded() &&
              result->sequence == sequence && result->publishedToConvolver,
          "worker builds the FIR and partition spectra away from audio "
          "processing");
    Check(convolver.IsTransitioning(),
          "prepared worker bank waits for an audio block boundary");
  }

  trackedAllocations.store(0, std::memory_order_relaxed);
  trackAllocations = true;
  float checksum = 0.0f;
  for (std::size_t sample = 0; sample < 4096; ++sample) {
    Convolver::InputFrame input{};
    input[0] = sample == 0 ? 1.0f : 0.0f;
    const auto output = convolver.Process(input, 1);
    checksum += std::abs(output[0]) + std::abs(output[1]);
  }
  trackAllocations = false;
  Check(trackedAllocations.load(std::memory_order_relaxed) == 0,
        "adopting and rendering a worker-prepared bank allocates no memory");
  Check(checksum > 0.0f && std::isfinite(checksum),
        "worker-prepared FIR reaches the audio renderer");

  tfdsp::EarlyReflectionWorker coalescingWorker(5.0);
  const auto firstSequence = coalescingWorker.Submit(request);
  auto first =
      coalescingWorker.WaitForLatestResult(std::chrono::milliseconds(1000));
  Check(first.has_value() && first->sequence == firstSequence &&
            first->Succeeded(),
        "rate-limited worker completes its first request");
  trackedAllocations.store(0, std::memory_order_relaxed);
  trackAllocations = true;
  std::size_t nonblockingSequence = 0;
  for (std::size_t update = 0; update < 32; ++update) {
    request.diffusion = static_cast<double>(update % 2);
    nonblockingSequence = coalescingWorker.Submit(request);
  }
  trackAllocations = false;
  Check(nonblockingSequence != 0 &&
            trackedAllocations.load(std::memory_order_relaxed) == 0,
        "audio-thread request publication is nonallocating and retains capacity");
  const auto submittedAt = std::chrono::steady_clock::now();
  request.diffusion = 0.0;
  coalescingWorker.Submit(request);
  request.diffusion = 0.5;
  coalescingWorker.Submit(request);
  request.diffusion = 1.0;
  const auto latestSequence = coalescingWorker.Submit(request);
  std::optional<tfdsp::EarlyReflectionBuildResult> latest;
  const auto deadline = submittedAt + std::chrono::seconds(2);
  while (std::chrono::steady_clock::now() < deadline) {
    auto candidate =
        coalescingWorker.WaitForLatestResult(std::chrono::milliseconds(250));
    if (candidate && candidate->sequence == latestSequence) {
      latest = std::move(candidate);
      break;
    }
  }
  const double elapsed = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - submittedAt)
                             .count();
  Check(latest.has_value() && latest->Succeeded(),
        "coalescing worker ultimately builds the latest submitted scene");
  Check(latest &&
            std::abs(latest->response.mixingPrediction.diffusionMultiplier -
                     0.75) < 1.0e-12,
        "superseded scene requests do not replace the latest result");
  Check(elapsed >= 0.14,
        "worker start rate respects its configured five-update-per-second cap");

  request.materials = tfdsp::MakeEarlyReflectionMaterials(1.0);
  const auto dampingSequence = coalescingWorker.Submit(request);
  auto dampingResult =
      coalescingWorker.WaitForLatestResult(std::chrono::milliseconds(1000));
  Check(dampingResult && dampingResult->Succeeded() &&
            dampingResult->sequence == dampingSequence &&
            dampingResult->geometryReused,
        "material-only updates reuse cached image-source geometry");

  request.materials.reflectionAmplitudes[0][0] = 1.01;
  const auto invalidMaterialSequence = coalescingWorker.Submit(request);
  auto invalidMaterialResult =
      coalescingWorker.WaitForLatestResult(std::chrono::milliseconds(1000));
  Check(invalidMaterialResult && !invalidMaterialResult->Succeeded() &&
            invalidMaterialResult->sequence == invalidMaterialSequence &&
            invalidMaterialResult->geometryReused,
        "cached-geometry builds still validate changed materials");
}
} // namespace

int main() {
  TestDefaultsAndLimits();
  TestSubSamplePathUsesCausalInterpolation();
  TestImageSourceGeometryAndMaterials();
  TestAdaptiveMixingPredictionAndRelativeTime();
  TestEchoDensityStatisticsAndRefinement();
  TestAllPathsAndMultipleSources();
  TestGeneratedFirMatchesAnalyticTapTrain();
  TestUnequalMaterialBandsProduceSpectralTilt();
  TestPartitionedConvolutionAndSuperposition();
  TestPhysicalTimingWithCompensatedLatency();
  TestRealtimeSafetyAndValidation();
  TestNumericalImpulseResponseTransition();
  TestInactiveSourceTailAndHistoryFlush();
  TestScenePackCompactionAndChaining();
  TestRateLimitedCoalescingWorkerAndPreparedBankHandoff();
  if (failures != 0) {
    std::cerr << failures << " early-reflection test(s) failed\n";
    return 1;
  }
  std::cout << "All early-reflection tests passed\n";
  return 0;
}
