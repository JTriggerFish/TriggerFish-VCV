// Developer-only inspection of the actual C++ packet preparation.
// stdin: sample rate, then whitespace-separated parameter-key/value pairs.
#include "crash_macros.hpp"
#include "tfdsp/percussion/crash_cymbal_parameters.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main() {
  try {
    float rate{};
    if (!(std::cin >> rate) || !std::isfinite(rate) || rate < 8000.f)
      throw std::invalid_argument("expected sample rate >= 8000");
    auto values = tfworkbench::DefaultCrashMacros();
    std::string key;
    while (std::cin >> key) {
      float value{};
      if (!(std::cin >> value) || !std::isfinite(value))
        throw std::invalid_argument("invalid value for " + key);
      bool found = false;
      for (std::size_t i = 0; i < values.size(); ++i) {
        const auto &d = tfworkbench::CrashMacroDescription(i);
        if (d.key != key) continue;
        if (value < d.minimum || value > d.maximum)
          throw std::invalid_argument("parameter outside UI range: " + key);
        values[i] = value;
        found = true;
        break;
      }
      if (!found) throw std::invalid_argument("unknown parameter: " + key);
    }
    const auto fit = tfworkbench::ApplyCrashMacros(
        tfworkbench::MetallicWorkbenchBaseFit(), values);
    const auto prepared = tfdsp::percussion::DefaultCrashCymbalParameters(rate, fit);
    std::cout << std::setprecision(10) << "[";
    bool first = true;
    for (const auto &m : prepared.modalField) {
      if (m.inputGain == 0.f || m.outputGain == 0.f) continue;
      if (!first) std::cout << ",";
      first = false;
      std::cout << "{\"frequency\":" << m.frequencyHz
        << ",\"packet\":" << m.packet
        << ",\"centre\":" << m.transportFrequencyHz
        << ",\"input\":" << m.inputGain
        << ",\"output\":" << m.outputGain
        << ",\"phase\":" << m.inputPhaseRadians
        << ",\"blur_hz\":" << m.phaseBandwidthHz
        << ",\"t60\":" << m.decaySeconds << "}";
    }
    std::cout << "]\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
