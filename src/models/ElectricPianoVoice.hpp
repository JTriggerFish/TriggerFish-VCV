#pragma once

#include "ElectricPianoParameters.hpp"

namespace tfdsp {

class ElectricPianoVoice
{
	// This constant must precede the responsibility fragments because it appears
	// in several of their fixed-size Eigen and std::array types.
	static constexpr int PickupOversamplingFactor =
		X4Resampler_Order7::ResamplingFactor;

#include "electric_piano/voice_setup.inl"
#include "electric_piano/voice_runtime.inl"
#include "electric_piano/voice_pickup.inl"
#include "electric_piano/voice_mechanics.inl"
#include "electric_piano/voice_timbre.inl"
#include "electric_piano/voice_modes.inl"
#include "electric_piano/voice_integration.inl"
#include "electric_piano/voice_mode_calibration.inl"
#include "electric_piano/voice_contact_calibration.inl"
#include "electric_piano/voice_magnetic_pickup.inl"
#include "electric_piano/voice_state.inl"
};

} // namespace tfdsp
