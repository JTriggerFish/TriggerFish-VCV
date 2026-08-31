#pragma once

#include "ElectricPianoParameters.hpp"

namespace tfdsp {

// Figure 11-8's transistor preamplifier and Figure 11-9's power stages are
// solved as nonlinear circuits; only the lamp/LDR routing is reduced.
class ElectricPianoAmplifier
{
#include "electric_piano/amplifier_public.inl"
#include "electric_piano/amplifier_state.inl"
};

} // namespace tfdsp
