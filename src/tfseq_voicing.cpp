#include "tfseq_voicing.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <tuple>

namespace tfseq {
namespace {

struct ToneSet {
  std::array<ChordTone, MaximumPolyphony> tones{};
  std::size_t count = 0;
};

struct VoicingRecipe {
  std::size_t nominalVoices;
  bool omitRoot;
  int registerOffset;
  int maximumSpan;
  std::array<int, 5> omissionOrder;
  std::array<int, 5> additionOrder;
};

constexpr VoicingRecipe BasicRecipe{
    0, false, 6, 24, {5, 1, 9, 11, 13}, {13, 11, 9, 5, 1}};
constexpr VoicingRecipe Rootless3Recipe{
    3, true, -2, 18, {1, 5, 9, 11, 13}, {13, 11, 9, 5, 1}};
constexpr VoicingRecipe Rootless4Recipe{
    4, true, -2, 18, {1, 5, 9, 11, 13}, {13, 11, 9, 5, 1}};

// A reduced `alt` voicing keeps its altered ninth before its altered fifth;
// at four voices both survive beside the guide tones. Ordinary extended
// chords instead protect the highest explicitly written colour.
constexpr std::array<int, 5> AlteredOmissionOrder{1, 11, 13, 5, 9};
constexpr std::array<int, 5> AlteredAdditionOrder{9, 11, 13, 5, 1};

const VoicingRecipe &RecipeFor(const VoicingStyle style) noexcept {
  if (style == VoicingStyle::Rootless3Notes)
    return Rootless3Recipe;
  if (style == VoicingStyle::Rootless4Notes)
    return Rootless4Recipe;
  return BasicRecipe;
}

struct AlteredTensionPair {
  ChordTone ninth;
  ChordTone fifth;
};

constexpr std::array<AlteredTensionPair, 4> AlteredTensions{{
    {{9, -1, 13}, {11, 1, 18}},
    {{9, 1, 15}, {11, 1, 18}},
    {{9, -1, 13}, {13, -1, 20}},
    {{9, 1, 15}, {13, -1, 20}},
}};

struct CandidateScore {
  std::int64_t countDifference = std::numeric_limits<std::int64_t>::max();
  std::int64_t totalMotion = std::numeric_limits<std::int64_t>::max();
  std::int64_t largestLeap = std::numeric_limits<std::int64_t>::max();
  std::int64_t movingVoices = std::numeric_limits<std::int64_t>::max();
  std::int64_t registerDistance = std::numeric_limits<std::int64_t>::max();
  int recipeRank = std::numeric_limits<int>::max();

  auto tuple() const noexcept {
    return std::tie(countDifference, totalMotion, largestLeap, movingVoices,
                    registerDistance, recipeRank);
  }
};

int PositiveMod12(const int value) noexcept {
  const int remainder = value % 12;
  return remainder < 0 ? remainder + 12 : remainder;
}

bool HasDegree(const ToneSet &set, const int degree) noexcept {
  for (std::size_t index = 0; index < set.count; ++index)
    if (set.tones[index].degree == degree)
      return true;
  return false;
}

void AddTone(ToneSet &set, const ChordTone tone) noexcept {
  if (set.count >= set.tones.size())
    return;
  const int pitchClass = PositiveMod12(tone.semitones);
  for (std::size_t index = 0; index < set.count; ++index)
    if (PositiveMod12(set.tones[index].semitones) == pitchClass)
      return;
  set.tones[set.count++] = tone;
}

ToneSet FormulaTones(const ChordValue &chord) noexcept {
  ToneSet result;
  for (const auto &tone : chord.tones)
    AddTone(result, tone);
  return result;
}

// `alt` denotes alternatives, not a simultaneous b9/#9 cluster. Each branch
// contains one altered ninth and one altered fifth. Explicit factor overrides
// bypass these automatic branches and may deliberately request any legal set.
ToneSet AlteredBranch(const ChordValue &chord, const int branch) noexcept {
  auto result = FormulaTones(chord);
  const auto &tensions = AlteredTensions[static_cast<std::size_t>(branch)];
  AddTone(result, tensions.ninth);
  AddTone(result, tensions.fifth);
  return result;
}

bool IsGuideTone(const ChordTone &tone, const ToneSet &pool) noexcept {
  if (tone.degree == 3)
    return true;
  // A suspension takes the third's harmonic role, while a sixth is the
  // second guide tone of a rootless 6/9 voicing. Treating either as a guide
  // keeps the named colour when the root and optional fifth are omitted.
  if (!HasDegree(pool, 3) && (tone.degree == 2 || tone.degree == 4))
    return true;
  if (tone.degree == 6 && !HasDegree(pool, 7))
    return true;
  return tone.degree == 7 && HasDegree(pool, 7);
}

int DegreeRank(const ChordTone &tone, const std::array<int, 5> &order,
               const bool protectAlteredFifth) noexcept {
  if (protectAlteredFifth && tone.degree == 5 && tone.accidental != 0)
    return 100;
  const auto found = std::find(order.begin(), order.end(), tone.degree);
  return found == order.end()
             ? 99
             : static_cast<int>(std::distance(order.begin(), found));
}

void EraseTone(ToneSet &set, const std::size_t erased) noexcept {
  for (std::size_t index = erased; index + 1 < set.count; ++index)
    set.tones[index] = set.tones[index + 1];
  if (set.count > 0)
    --set.count;
}

ToneSet ApplyRecipe(const ToneSet &pool, const ChordValue &chord,
                    const VoicingStyle style,
                    const std::size_t previousCount) noexcept {
  if (chord.hasFactorOverride)
    return pool;

  const auto &recipe = RecipeFor(style);
  const bool hasSeventh = HasDegree(pool, 7);
  const bool hasSixthGuides = HasDegree(pool, 3) && HasDegree(pool, 6);
  // Seventh chords use third/seventh (or suspension/seventh) guides. Sixth
  // chords use third/sixth guides. Plain triads still use the basic fallback
  // instead of degenerating into an ambiguous rootless dyad.
  const bool rootless = recipe.omitRoot && (hasSeventh || hasSixthGuides);
  const std::size_t nominal =
      recipe.nominalVoices == 0 ? pool.count : recipe.nominalVoices;
  const std::size_t minimum = std::min<std::size_t>(3, pool.count);
  const std::size_t maximum =
      rootless ? std::min(nominal, pool.count) : pool.count;
  std::size_t target = std::clamp(nominal, minimum, maximum);
  if (previousCount >= minimum && previousCount <= maximum)
    target = previousCount;

  ToneSet selected;
  if (rootless) {
    for (std::size_t index = 0; index < pool.count; ++index)
      if (IsGuideTone(pool.tones[index], pool))
        AddTone(selected, pool.tones[index]);
    std::array<std::size_t, MaximumPolyphony> order{};
    for (std::size_t index = 0; index < pool.count; ++index)
      order[index] = index;
    std::stable_sort(
        order.begin(), order.begin() + pool.count,
        [&](const std::size_t left, const std::size_t right) {
          const auto &additionOrder =
              chord.altered ? AlteredAdditionOrder : recipe.additionOrder;
          return DegreeRank(pool.tones[left], additionOrder, false) <
                 DegreeRank(pool.tones[right], additionOrder, false);
        });
    for (std::size_t position = 0;
         position < pool.count && selected.count < target; ++position) {
      const auto &tone = pool.tones[order[position]];
      if (tone.degree != 1)
        AddTone(selected, tone);
    }
  } else {
    selected = pool;
    while (selected.count > target) {
      std::size_t best = selected.count;
      int bestRank = std::numeric_limits<int>::max();
      for (std::size_t index = 0; index < selected.count; ++index) {
        if (IsGuideTone(selected.tones[index], selected))
          continue;
        const auto &omissionOrder =
            chord.altered ? AlteredOmissionOrder : recipe.omissionOrder;
        const int rank = DegreeRank(selected.tones[index], omissionOrder, true);
        if (rank < bestRank) {
          bestRank = rank;
          best = index;
        }
      }
      if (best == selected.count)
        break;
      EraseTone(selected, best);
    }
  }
  return selected;
}

CandidateScore
ScoreCandidate(const std::array<int, MaximumPolyphony> &candidate,
               const std::size_t candidateCount,
               const std::array<int, MaximumPolyphony> &previous,
               const std::size_t previousCount, const int registerTarget,
               const int recipeRank) noexcept {
  CandidateScore score;
  score.countDifference = std::abs(static_cast<int>(candidateCount) -
                                   static_cast<int>(previousCount));
  score.totalMotion = 0;
  score.largestLeap = 0;
  score.movingVoices = 0;
  if (previousCount > 0 && candidateCount > 0) {
    std::array<int, MaximumPolyphony> orderedPrevious = previous;
    std::sort(orderedPrevious.begin(), orderedPrevious.begin() + previousCount);
    const std::size_t paired = std::min(candidateCount, previousCount);
    for (std::size_t index = 0; index < paired; ++index) {
      const std::size_t candidateIndex =
          paired == 1 ? candidateCount / 2
                      : index * (candidateCount - 1) / (paired - 1);
      const std::size_t previousIndex =
          paired == 1 ? previousCount / 2
                      : index * (previousCount - 1) / (paired - 1);
      const std::int64_t movement = std::llabs(
          static_cast<long long>(candidate[candidateIndex]) -
          static_cast<long long>(orderedPrevious[previousIndex]));
      score.totalMotion += movement;
      score.largestLeap = std::max(score.largestLeap, movement);
      if (movement != 0)
        ++score.movingVoices;
    }
    score.totalMotion += score.countDifference * 12;
  }
  const int low = candidateCount ? candidate[0] : registerTarget;
  const int high =
      candidateCount ? candidate[candidateCount - 1] : registerTarget;
  score.registerDistance =
      std::llabs(static_cast<long long>(low) + static_cast<long long>(high) -
                 2LL * static_cast<long long>(registerTarget));
  score.recipeRank = recipeRank;
  return score;
}

void ConsiderCandidate(const std::array<int, MaximumPolyphony> &candidate,
                       const std::size_t count,
                       const std::array<int, MaximumPolyphony> &previous,
                       const std::size_t previousCount,
                       const int registerTarget, const int maximumSpan,
                       const int recipeRank, CandidateScore &bestScore,
                       VoicingResult &best) noexcept {
  if (count == 0 || candidate[count - 1] - candidate[0] > maximumSpan)
    return;
  const auto score = ScoreCandidate(candidate, count, previous, previousCount,
                                    registerTarget, recipeRank);
  if (best.count == 0 || score.tuple() < bestScore.tuple()) {
    bestScore = score;
    best.count = count;
    std::copy_n(candidate.begin(), count, best.semitones.begin());
  }
}

VoicingResult RealizeToneSet(const ToneSet &selected, const VoicingStyle style,
                             const int rootSemitone,
                             const std::array<int, MaximumPolyphony> &previous,
                             const std::size_t previousCount) noexcept {
  VoicingResult result;
  if (selected.count == 0)
    return result;

  CandidateScore bestScore;
  const auto &recipe = RecipeFor(style);
  const int registerTarget = rootSemitone + recipe.registerOffset;

  // On the first basic chord, honour the promised obvious root-position stack.
  // Thereafter inversions participate in the same deterministic policy.
  if (previousCount == 0 && style == VoicingStyle::Basic) {
    result.count = selected.count;
    for (std::size_t index = 0; index < selected.count; ++index)
      result.semitones[index] = rootSemitone + selected.tones[index].semitones;
    std::sort(result.semitones.begin(),
              result.semitones.begin() + result.count);
    return result;
  }

  struct OrderedTone {
    int pitchClass = 0;
    int degree = 0;
  };
  std::array<OrderedTone, MaximumPolyphony> orderedTones{};
  for (std::size_t index = 0; index < selected.count; ++index) {
    orderedTones[index] = {PositiveMod12(selected.tones[index].semitones),
                           selected.tones[index].degree};
  }
  std::sort(orderedTones.begin(), orderedTones.begin() + selected.count,
            [](const OrderedTone &left, const OrderedTone &right) {
              return left.pitchClass < right.pitchClass;
            });

  const bool useRootlessForms =
      style != VoicingStyle::Basic && HasDegree(selected, 3) &&
      HasDegree(selected, 7) && !HasDegree(selected, 1);

  int rank = 0;
  if (useRootlessForms) {
    auto considerForm = [&](const std::array<int, MaximumPolyphony> &degrees,
                            const std::size_t degreeCount) {
      std::array<bool, MaximumPolyphony> used{};
      std::array<int, MaximumPolyphony> orderedPitchClasses{};
      std::size_t orderedCount = 0;
      auto appendDegree = [&](const int degree) {
        for (std::size_t tone = 0; tone < selected.count; ++tone) {
          if (!used[tone] && selected.tones[tone].degree == degree) {
            used[tone] = true;
            orderedPitchClasses[orderedCount++] =
                PositiveMod12(selected.tones[tone].semitones);
            return;
          }
        }
      };
      for (std::size_t index = 0; index < degreeCount; ++index)
        appendDegree(degrees[index]);
      for (std::size_t tone = 0; tone < selected.count; ++tone) {
        if (!used[tone])
          orderedPitchClasses[orderedCount++] =
              PositiveMod12(selected.tones[tone].semitones);
      }
      if (orderedCount != selected.count)
        return;
      std::array<int, MaximumPolyphony> relative{};
      relative[0] = orderedPitchClasses[0];
      for (std::size_t voice = 1; voice < orderedCount; ++voice) {
        int pitch = orderedPitchClasses[voice];
        while (pitch <= relative[voice - 1])
          pitch += 12;
        relative[voice] = pitch;
      }
      for (int octave = -3; octave <= 3; ++octave) {
        std::array<int, MaximumPolyphony> candidate{};
        for (std::size_t voice = 0; voice < orderedCount; ++voice)
          candidate[voice] = rootSemitone + relative[voice] + octave * 12;
        ConsiderCandidate(candidate, orderedCount, previous, previousCount,
                          registerTarget, recipe.maximumSpan, rank++, bestScore,
                          result);
      }
    };

    const int colour = HasDegree(selected, 13)   ? 13
                       : HasDegree(selected, 5)  ? 5
                       : HasDegree(selected, 11) ? 11
                                                 : 0;
    const int tension = HasDegree(selected, 9)                    ? 9
                        : HasDegree(selected, 11) && colour != 11 ? 11
                                                                  : 0;
    std::array<int, MaximumPolyphony> aForm{3, colour, 7, tension};
    std::array<int, MaximumPolyphony> bForm{7, tension, 3, colour};
    considerForm(aForm, 4);
    considerForm(bForm, 4);
    if (result.count != 0)
      return result;
  }

  for (std::size_t inversion = 0; inversion < selected.count; ++inversion) {
    // Bill Evans-style A/B forms begin on the third or seventh. Restricting
    // the candidate family here prevents a low fifth or ninth from winning a
    // purely geometric score while keeping explicit non-guide-tone formulas
    // usable as ordinary close inversions.
    std::array<int, MaximumPolyphony> relative{};
    int prior = orderedTones[inversion].pitchClass;
    relative[0] = prior;
    for (std::size_t voice = 1; voice < selected.count; ++voice) {
      int pitchClass =
          orderedTones[(inversion + voice) % selected.count].pitchClass;
      while (pitchClass <= prior)
        pitchClass += 12;
      relative[voice] = pitchClass;
      prior = pitchClass;
    }
    for (int octave = -3; octave <= 3; ++octave) {
      std::array<int, MaximumPolyphony> candidate{};
      for (std::size_t voice = 0; voice < selected.count; ++voice)
        candidate[voice] = rootSemitone + relative[voice] + octave * 12;
      ConsiderCandidate(candidate, selected.count, previous, previousCount,
                        registerTarget, recipe.maximumSpan, rank++, bestScore,
                        result);
    }
  }
  return result;
}

} // namespace

VoicingResult
RealizeChordVoicing(const ChordValue &chord, const VoicingStyle style,
                    const int rootSemitone,
                    const std::array<int, MaximumPolyphony> &previous,
                    const std::size_t previousCount) noexcept {
  VoicingResult best;
  CandidateScore bestScore;
  const int branches = chord.altered && !chord.hasFactorOverride ? 4 : 1;
  for (int branch = 0; branch < branches; ++branch) {
    const auto pool =
        branches == 1 ? FormulaTones(chord) : AlteredBranch(chord, branch);
    const auto selected = ApplyRecipe(pool, chord, style, previousCount);
    const auto candidate =
        RealizeToneSet(selected, style, rootSemitone, previous, previousCount);
    if (candidate.count == 0)
      continue;
    const auto score = ScoreCandidate(
        candidate.semitones, candidate.count, previous, previousCount,
        rootSemitone + RecipeFor(style).registerOffset, branch);
    if (best.count == 0 || score.tuple() < bestScore.tuple()) {
      best = candidate;
      bestScore = score;
    }
  }
  return best;
}

std::size_t MaximumVoicingCount(const ChordValue &chord,
                                const VoicingStyle style) noexcept {
  if (chord.hasFactorOverride)
    return chord.tones.size();
  const auto &recipe = RecipeFor(style);
  const auto available = chord.altered ? 5 : chord.tones.size();
  return recipe.nominalVoices == 0
             ? available
             : std::min<std::size_t>(recipe.nominalVoices, available);
}

} // namespace tfseq
