#pragma once

#include "crash_cymbal_parameters.hpp"
#include "modal_constraint.hpp"
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
  float bloomTransferEnergy{};
  float modalBody{};
  float output{};
};

// Mono crash body. Stereo is an observation/presentation concern and never
// duplicates the contact, dispersion, or modal state owned here.
class CrashCymbal {
public:
  void Prepare(float sampleRate, const CrashCymbalParameters &parameters);
  void Prepare(const CrashCymbalPreparedParameters &parameters);
  void Reset() noexcept;
  void Trigger(const CrashCymbalHit &hit) noexcept;
  CrashCymbalFrame ProcessFrame() noexcept;
  float Process() noexcept;
  void SetMute(float amount) noexcept;

  double StoredBodyEnergy() const noexcept;
  float StoredBodyEnergyCentroidHz() const noexcept;

private:
  void PrepareComponents(float sampleRate,
                         const CrashCymbalParameters &parameters);
  ContactExciterParameters ContactParameters(
      const CrashCymbalHit &hit) const noexcept;
  void SetExcitationProjection(float location, float strength) noexcept;

  ContactExciter contact_{};
  CrashModalField modalField_{};
  ObservationModel<2> observation_{};
  ModalConstraintController modalConstraint_{};
  CrashCymbalParameters parameters_{};
  CrashModalField::Projection fieldProjection_{};
  float bodyExcitationGain_{.05f};
  float bodyDriveScale_{1.f};
  MetallicPlateRouting routing_{};
  float sampleRate_{48000.f};
  bool hasProcessedSinceReset_{};
};

} // namespace tfdsp::percussion
