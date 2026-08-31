#pragma once

#include "fir_hilbert_transform.hpp"
#include "quadrature_oscillator.hpp"

namespace tfdsp::percussion {

// Single-sideband frequency shifter. Positive and negative offsets share the
// same phase-continuous oscillator; zero offset is the exactly delayed input.
class FrequencyShifter {
public:
  static constexpr std::size_t LatencySamples =
      FirHilbertTransform::LatencySamples;

  void Prepare(const float sampleRate) {
    oscillator_.Prepare(sampleRate);
    Reset();
  }

  void Reset() noexcept {
    hilbert_.Reset();
    oscillator_.Reset();
  }

  void SetShiftHz(const float shiftHz) noexcept {
    oscillator_.SetFrequencyHz(shiftHz);
  }

  float Process(const float input) noexcept {
    const auto analytic = hilbert_.Process(input);
    const auto rotation = oscillator_.Process();
    return analytic.real * rotation.cosine -
           analytic.quadrature * rotation.sine;
  }

private:
  FirHilbertTransform hilbert_{};
  QuadratureOscillator oscillator_{};
};

} // namespace tfdsp::percussion
