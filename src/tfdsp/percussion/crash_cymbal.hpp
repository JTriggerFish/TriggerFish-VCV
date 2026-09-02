#pragma once

#include "crash_cymbal_parameters.hpp"
#include "modal_constraint.hpp"
#include "passive_constraint.hpp"
#include "tfdsp/finite_audio.hpp"

#include <array>
#include <cstdint>

namespace tfdsp::percussion {

struct CrashCymbalHit {
  float strength{.8f};
  float location{1.f};
  float hardness{.65f};
  std::uint32_t seed{1};
  float implement{.75f};
  // Contact gesture length. A stick/mallet remains an impact; with a brush,
  // this extends one trigger into a correlated stream of bristle contacts.
  float contactSpread{.2f};
};

struct CrashCymbalFrame {
  float directContact{};
  float dispersion{};
  float sparseModes{};
  float denseModes{};
  float turbulentResidual{};
  float denseResidual{};
  float output{};
};

// Mono crash body. Stereo is an observation/presentation concern and never
// duplicates the contact, dispersion, or modal state owned here.
class CrashCymbal {
public:
  void Prepare(float sampleRate, const CrashCymbalParameters &parameters);
  void Reset() noexcept;
  void Trigger(const CrashCymbalHit &hit) noexcept;
  CrashCymbalFrame ProcessFrame() noexcept;
  float Process() noexcept;
  void SetMute(float amount) noexcept;

  float MinimumBodyDelaySamples() const noexcept;

private:
  ContactExciterParameters ContactParameters(
      const CrashCymbalHit &hit) const noexcept;
  void SetExcitationProjection(float location, float strength) noexcept;

  ContactExciter contact_{};
  DispersionLoop dispersion_{};
  TurbulentResidual turbulence_{};
  CrashSparseModes sparseModes_{};
  CrashDenseModes denseModes_{};
  CrashDenseModes denseExtensionModes_{};
  ObservationModel<4> observation_{};
  DynamicLossController delayConstraint_{};
  ModalConstraintController modalConstraint_{};
  CrashCymbalParameters parameters_{};
  CrashSparseModes::Projection sparseProjection_{};
  CrashDenseModes::Projection denseProjection_{};
  CrashDenseModes::Projection denseExtensionProjection_{};
  float bodyDriveScale_{1.f};
  float sampleRate_{48000.f};
  float sparseBloomGain_{};
  float denseBypassGain_{};
  bool sparseEnabled_{true};
  bool denseEnabled_{true};
  bool denseExtensionEnabled_{};
  bool turbulenceEnabled_{};
};

} // namespace tfdsp::percussion
