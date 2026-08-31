#pragma once

#include "early_reflection_scene.hpp"
#include "finite_audio.hpp"

namespace tfdsp {
namespace early_reflection_detail {
class TwoPoleLowPass {
  double coefficient_{};
  double state1_{};
  double state2_{};

public:
  void Prepare(const double cutoffHz, const double sampleRate) noexcept {
    coefficient_ =
        1.0 - std::exp(-2.0 * EarlyReflectionPi * cutoffHz / sampleRate);
    state1_ = 0.0;
    state2_ = 0.0;
  }

  double Process(const double input) noexcept {
    const double safeInput = FiniteNormalOrZero(input);
    state1_ += coefficient_ * (safeInput - state1_);
    state2_ += coefficient_ * (state1_ - state2_);
    state1_ = FiniteNormalOrZero(state1_);
    state2_ = FiniteNormalOrZero(state2_);
    return state2_;
  }
};

class ComplementaryBandFilter {
  std::array<TwoPoleLowPass, EarlyReflectionBandCount - 1> lowPasses_{};
  std::size_t band_{};

public:
  void Prepare(const EarlyReflectionConfig &config,
               const std::size_t band) noexcept {
    band_ = band;
    for (std::size_t crossover = 0; crossover < lowPasses_.size(); ++crossover)
      lowPasses_[crossover].Prepare(config.crossoverHz[crossover],
                                    config.sampleRate);
  }

  double Process(const double input) noexcept {
    double residual = input;
    for (std::size_t crossover = 0; crossover < lowPasses_.size();
         ++crossover) {
      const double low = lowPasses_[crossover].Process(residual);
      if (band_ == crossover)
        return low;
      residual -= low;
    }
    return residual;
  }
};
} // namespace early_reflection_detail

} // namespace tfdsp
