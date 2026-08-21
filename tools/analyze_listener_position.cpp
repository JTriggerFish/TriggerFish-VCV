#include "tfdsp/room_reverb.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>
#include <tuple>
#include <vector>

namespace {

struct Candidate {
  std::array<double, 3> position{};
  double score{};
  double worst{};
};

double SourceScore(const std::vector<tfdsp::EarlyReflectionPath> &paths,
                   const std::size_t sourceIndex) {
  constexpr double binSeconds = 0.00025;
  constexpr std::size_t bins = 600;
  std::array<double, bins> energyByTime{};
  std::array<double, 2> stereoEnergy{};
  double total = 0.0;
  for (const auto &path : paths) {
    if (path.sourceIndex != sourceIndex)
      continue;
    double bandEnergy = 0.0;
    for (const double gain : path.bandGains)
      bandEnergy += gain * gain;
    bandEnergy /= path.bandGains.size();
    const auto bin = std::min(
        bins - 1,
        static_cast<std::size_t>(path.excessDelaySeconds / binSeconds));
    energyByTime[bin] += bandEnergy;
    for (std::size_t channel = 0; channel < stereoEnergy.size(); ++channel)
      stereoEnergy[channel] +=
          bandEnergy * path.outputGains[channel] * path.outputGains[channel];
    total += bandEnergy;
  }
  if (total <= 0.0)
    return std::numeric_limits<double>::infinity();
  double concentration = 0.0;
  double peak = 0.0;
  for (const double energy : energyByTime) {
    const double normalized = energy / total;
    concentration += normalized * normalized;
    peak = std::max(peak, normalized);
  }
  const double balance =
      std::abs(stereoEnergy[0] - stereoEnergy[1]) / total;
  return peak + 3.0 * concentration + 0.10 * balance;
}

Candidate Evaluate(const std::array<double, 3> &listener) {
  constexpr std::array<std::array<float, 2>, 3> roomControls{{
      {{0.30f, 0.35f}},
      {{0.50f, 0.50f}},
      {{0.72f, 0.67f}},
  }};
  constexpr std::array<std::array<double, 3>, 3> sourcePositions{{
      {{0.50, 0.35, 0.42}},
      {{0.30, 0.28, 0.38}},
      {{0.72, 0.44, 0.52}},
  }};
  tfdsp::EarlyReflectionConfig config;
  const auto materials = tfdsp::MakeEarlyReflectionMaterials(0.18);
  double sum = 0.0;
  double worst = 0.0;
  std::size_t cases = 0;
  for (const auto &roomControl : roomControls) {
    tfdsp::RoomReverbControls controls;
    controls.space = roomControl[0];
    controls.aspect = roomControl[1];
    controls.listener = {static_cast<float>(listener[0]),
                         static_cast<float>(listener[1]),
                         static_cast<float>(listener[2])};
    const auto room = tfdsp::RoomReverb::MakeRoom(controls);
    std::vector<tfdsp::EarlyReflectionSource> sources;
    for (const auto &position : sourcePositions)
      sources.push_back(tfdsp::MakeEarlyReflectionSource(room, position));
    const auto paths = tfdsp::EnumerateEarlyReflectionPaths(
        config, room, sources, materials, 0.75);
    for (std::size_t source = 0; source < sources.size(); ++source) {
      const double score = SourceScore(paths, source);
      sum += score;
      worst = std::max(worst, score);
      ++cases;
    }
  }
  const double mean = sum / static_cast<double>(cases);
  return {listener, mean + 0.35 * worst, worst};
}

} // namespace

int main(const int argc, const char *const *argv) {
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "previous " << Evaluate({0.50, 0.68, 0.45}).score << '\n';
  std::cout << "default  " << Evaluate({0.5, 0.682, 0.45}).score << '\n';
  std::cout << "centre   " << Evaluate({0.50, 0.50, 0.45}).score << '\n';
  if (argc > 1 && std::string_view(argv[1]) == "--summary")
    return 0;

  std::vector<Candidate> candidates;
  // Keep at least 15% of each room dimension between the listener and a
  // boundary.  Positions closer than this can make the direct sound and one
  // wall reflection nearly coincident even when the global collision score is
  // otherwise attractive.
  for (double x = 0.15; x <= 0.8501; x += 0.025)
    for (double y = 0.15; y <= 0.8501; y += 0.025)
      for (double z = 0.20; z <= 0.8001; z += 0.025)
        candidates.push_back(Evaluate({x, y, z}));
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate &left, const Candidate &right) {
              return left.score < right.score;
            });
  for (std::size_t index = 0; index < 12; ++index) {
    const auto &candidate = candidates[index];
    std::cout << index + 1 << "  " << candidate.position[0] << "  "
              << candidate.position[1] << "  " << candidate.position[2]
              << "  score=" << candidate.score
              << " worst=" << candidate.worst << '\n';
  }
}
