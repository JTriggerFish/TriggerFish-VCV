#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cmath>

namespace tfdsp::percussion {

// Finite-volume energy diffusion on sorted packet centres, not oscillator
// phases. Coordinates are f/1000 Hz; energy is measured in unit-strike units.
// See docs/TfPercussion-spectral-diffusion.md for the equation and limitations.
// No allocation, iteration to convergence, output normalization or energy sink.
template <std::size_t Capacity> class ModalSpectralDiffusion {
public:
  template <typename Reference = float>
  void Prepare(const std::array<float, Capacity> &frequencyHz,
               const std::array<Reference, Capacity> &referenceEnergy,
               const std::size_t count, const float sampleRate) noexcept {
    packets_ = std::min(count, Capacity);
    bins_ = 0;
    reference_.fill(0.0);
    for (std::size_t i = 0; i < packets_; ++i) {
      const double x = std::clamp(double(frequencyHz[i]), 0.0,
                                  double(sampleRate) * .5) / 1000.0;
      if (!bins_ || x != centre_[bins_ - 1]) centre_[bins_++] = x;
      const auto bin = bins_ - 1;
      bin_[i] = bin;
      reference_[bin] += referenceEnergy[i];
      packetReference_[i] = referenceEnergy[i];
    }
    const double upper = bins_ > 1 ? centre_[bins_-1] +
        .5 * (centre_[bins_-1]-centre_[bins_-2]) : double(sampleRate) / 2000.0;
    for (std::size_t i = 0; i < bins_; ++i) {
      const double left = i ? .5 * (centre_[i-1] + centre_[i]) : 0.0;
      // Extend the last represented cell by half its neighbour spacing. A
      // Nyquist-wide empty cell would change diffusion when only Fs changes.
      const double right = i+1 < bins_ ? .5 * (centre_[i] + centre_[i+1])
          : std::min(upper, double(sampleRate) / 2000.0);
      width_[i] = right - left;
      if (i+1 < bins_)
        geometry_[i] = right / (centre_[i+1] - centre_[i]);
    }
  }

  // A semi-implicit step: freeze nonnegative conductivities at the old state
  // and solve the resulting tridiagonal M-matrix. Closed boundary fluxes make
  // the step conservative; its inverse is nonnegative even for stiff steps.
  void Process(const std::array<float, Capacity> &energy,
               std::array<float, Capacity> &result,
               const double referenceEnergy, const double step,
               const double nonlinearity) noexcept {
    std::copy_n(energy.begin(), packets_, result.begin());
    if (bins_ < 2 || !(step > 0) || !(referenceEnergy > 0)) return;
    std::fill_n(original_.begin(), bins_, 0.0);
    for (std::size_t i = 0; i < packets_; ++i)
      original_[bin_[i]] += double(energy[i]) / referenceEnergy;
    for (std::size_t i = 0; i < bins_; ++i)
      density_[i] = original_[i] / width_[i];
    BuildConductance(step, std::clamp(nonlinearity, 0.0, 1.0));
    Solve();
    for (std::size_t i = 0; i < packets_; ++i) {
      const auto bin = bin_[i];
      // Coincident handles share one cell. Existing energy proportions are
      // retained; an initially silent cell uses declared excitation weights.
      const double share = original_[bin] > 0
          ? (double(energy[i]) / referenceEnergy) / original_[bin]
          : ReferenceShare(i);
      result[i] = float(density_[bin] * width_[bin] * referenceEnergy * share);
    }
  }

  double CellWidth(const std::size_t packet) const noexcept {
    return width_[bin_[packet]];
  }

  // Quadrature mass for sampling a continuous excitation density on this
  // grid. Coincident handles split one cell rather than duplicating its mass.
  double ExcitationCellWeight(const std::size_t packet) const noexcept {
    const auto bin = bin_[packet];
    return width_[bin] * ReferenceShare(packet);
  }

private:
  // Normalize actual weights, not a floored denominator: a very dark strike
  // spectrum must not prevent a silent upper cell receiving stored energy.
  double ReferenceShare(const std::size_t packet) const noexcept {
    const double total = reference_[bin_[packet]];
    return total > 0 ? packetReference_[packet] / total : 0.0;
  }

  void BuildConductance(const double step, const double nonlinearity) noexcept {
    for (std::size_t i = 0; i+1 < bins_; ++i) {
      const double a = density_[i], b = density_[i+1];
      // Divided difference of rho^3/3: consistent rho^2 conductivity, including
      // interfaces with zero density. At n=0 this becomes linear diffusion.
      const double quadratic = (a*a + a*b + b*b) / 3.0;
      // Interpolate the energy exponent, not a constant-conductivity floor.
      // n=0 is linear diffusion; n=1 is the quadratic plate closure. Every
      // n>0 vanishes with energy, so quiet strikes and late tails do not keep
      // leaking into a wash at an unrelated, energy-independent rate.
      const double conductivity = nonlinearity == 0 ? 1.0 :
          nonlinearity == 1 ? quadratic : std::pow(quadratic, nonlinearity);
      conductance_[i] = step * geometry_[i] * conductivity;
    }
    conductance_[bins_ - 1] = 0.0;
  }

  void Solve() noexcept {
    // Cancellation-free Thomas elimination: keep the positive diagonal
    // remainder separately from its right-edge conductance.
    double remainder = width_[0];
    diagonal_[0] = remainder + conductance_[0];
    rhs_[0] = original_[0];
    for (std::size_t i = 1; i < bins_; ++i) {
      const double left = conductance_[i-1];
      const double ratio = left / (remainder + left);
      remainder = width_[i] + ratio * remainder;
      diagonal_[i] = remainder + conductance_[i];
      rhs_[i] = original_[i] + ratio * rhs_[i-1];
    }
    density_[bins_-1] = rhs_[bins_-1] / diagonal_[bins_-1];
    for (std::size_t i = bins_-1; i-- > 0;)
      density_[i] = (rhs_[i] + conductance_[i] * density_[i+1]) / diagonal_[i];
  }

  std::array<std::size_t, Capacity> bin_{};
  std::array<double, Capacity> centre_{}, width_{}, geometry_{}, reference_{};
  std::array<double, Capacity> packetReference_{}, original_{}, density_{};
  std::array<double, Capacity> conductance_{}, diagonal_{}, rhs_{};
  std::size_t packets_{}, bins_{};
};

} // namespace tfdsp::percussion
