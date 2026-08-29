#pragma once

#include "tfseq.hpp"

#include <array>
#include <cstddef>

namespace tfseq {

struct VoicingResult {
  std::array<int, MaximumPolyphony> semitones{};
  std::size_t count = 0;
};

// Realize one semantic chord at an integer-semitone root. The routine is
// deterministic, allocation-free, and intended to run only at chord attacks.
// previous may be empty for the first chord in a phrase.
VoicingResult RealizeChordVoicing(
    const ChordValue &chord, VoicingStyle style, int rootSemitone,
    const std::array<int, MaximumPolyphony> &previous,
    std::size_t previousCount) noexcept;

// Conservative event-workspace sizing for a named chord under a recipe.
std::size_t MaximumVoicingCount(const ChordValue &chord,
                                VoicingStyle style) noexcept;

} // namespace tfseq
