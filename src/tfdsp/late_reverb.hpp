#pragma once

#include "cubic_fractional_delay_bank.hpp"
#include "finite_audio.hpp"

#include "late_reverb_coefficient_sets.hpp"
#include "multiband_decay_filter_bank.hpp"
#include "reverb_defaults.hpp"
#include "smooth_random_modulator.hpp"
#include "velvet_feedback_matrix.hpp"
#include "windowed_pitch_shifter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace tfdsp {

struct LateReverbControls {
  float decay{reverb_defaults::Decay};
  float damping{reverb_defaults::Damping};
  float diffusion{reverb_defaults::Diffusion};
  float modulation{reverb_defaults::Modulation};
  float shimmer{reverb_defaults::Shimmer};
  std::array<float, 3> roomDimensionsMetres{
      reverb_defaults::RoomDimensionsMetres};
};

class LateReverb {
#include "late_reverb/core.inl"
#include "late_reverb/public_api.inl"
};

} // namespace tfdsp
