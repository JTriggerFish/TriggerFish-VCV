#pragma once

#include "cubic_fractional_delay.hpp"
#include "finite_audio.hpp"

#include "early_reflections_worker.hpp"
#include "late_reverb.hpp"
#include "reverb_defaults.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

namespace tfdsp {

struct RoomReverbFrame {
  std::array<float, 2> direct{};
  std::array<float, 2> wet{};
};

struct RoomReverbControls {
  float space{reverb_defaults::Space};
  float aspect{reverb_defaults::Aspect};
  // Horizontal centring avoids a default ear bias; depth asymmetry avoids the
  // largest set of coincident image-source paths at the exact room centre.
  std::array<float, 3> listener{reverb_defaults::Listener};
  float preDelay{reverb_defaults::PreDelay};
  float decay{reverb_defaults::Decay};
  float damping{reverb_defaults::Damping};
  float diffusion{reverb_defaults::Diffusion};
  float modulation{reverb_defaults::Modulation};
  float shimmer{reverb_defaults::Shimmer};
  float width{reverb_defaults::Width};
  float balance{reverb_defaults::Balance};
  float lowCut{reverb_defaults::LowCut};
  float highCut{reverb_defaults::HighCut};
};

class RoomReverb {
#include "room_reverb/state.inl"
#include "room_reverb/direct_path.inl"
#include "room_reverb/automatic_balance.inl"
#include "room_reverb/scene_updates.inl"
#include "room_reverb/public_api.inl"
};

} // namespace tfdsp
