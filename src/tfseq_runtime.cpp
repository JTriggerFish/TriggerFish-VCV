#include "tfseq.hpp"
#include "tfseq_voicing.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace tfseq {
namespace {

bool TransformEnabled(const Transform &transform, std::uint64_t cycle,
                      std::uint64_t seed, std::uint64_t salt) noexcept {
  (void)salt;
  if (transform.condition == TransformCondition::Always)
    return true;
  if (transform.condition == TransformCondition::Every)
    return transform.period > 0 && (cycle + 1) % transform.period == 0;
  // Source identity, rather than the destination lane's salt, makes a
  // sequence-level conditional transform take one coherent decision across
  // notes, articulation, and control lanes. Lane-local transforms naturally
  // have distinct source spans.
  const auto sourceKey = static_cast<std::uint64_t>(static_cast<std::uint32_t>(
                             transform.span.begin + 1))
                             << 32 |
                         static_cast<std::uint32_t>(transform.span.end + 1);
  std::uint64_t hash = seed ^ (cycle + 0x9e3779b97f4a7c15ULL) ^ sourceKey;
  hash ^= hash >> 30;
  hash *= 0xbf58476d1ce4e5b9ULL;
  hash ^= hash >> 27;
  hash *= 0x94d049bb133111ebULL;
  hash ^= hash >> 31;
  const double unit = static_cast<double>(hash >> 11) / 9007199254740992.0;
  return unit < transform.probability;
}

std::size_t TransformIndex(std::size_t index, std::size_t size,
                           const std::vector<Transform> &transforms,
                           std::uint64_t cycle, std::uint64_t seed,
                           std::uint64_t salt) noexcept {
  if (size == 0)
    return 0;
  std::int64_t mapped = static_cast<std::int64_t>(index % size);
  for (std::size_t transformIndex = 0; transformIndex < transforms.size();
       ++transformIndex) {
    const auto &transform = transforms[transformIndex];
    if (!TransformEnabled(transform, cycle, seed,
                          salt + transformIndex * 0x10001ULL))
      continue;
    if (transform.kind == TransformKind::Reverse)
      mapped = static_cast<std::int64_t>(size) - 1 - mapped;
    else if (transform.kind == TransformKind::Rotate) {
      mapped += transform.integer;
      mapped %= static_cast<std::int64_t>(size);
      if (mapped < 0)
        mapped += static_cast<std::int64_t>(size);
    }
  }
  return static_cast<std::size_t>(mapped);
}

std::int64_t SaturatingAdd(std::int64_t left, std::int64_t right) noexcept {
  if (right > 0 && left > std::numeric_limits<std::int64_t>::max() - right)
    return std::numeric_limits<std::int64_t>::max();
  if (right < 0 && left < std::numeric_limits<std::int64_t>::min() - right)
    return std::numeric_limits<std::int64_t>::min();
  return left + right;
}

std::int64_t Transpose(const std::vector<Transform> &transforms,
                       TransformKind kind, std::uint64_t cycle,
                       std::uint64_t seed, std::uint64_t salt) noexcept {
  std::int64_t result = 0;
  for (std::size_t index = 0; index < transforms.size(); ++index) {
    const auto &transform = transforms[index];
    if (transform.kind == kind &&
        TransformEnabled(transform, cycle, seed, salt + index * 0x10001ULL))
      result = SaturatingAdd(result, transform.integer);
  }
  return result;
}

double TimeScale(const std::vector<Transform> &transforms, std::uint64_t cycle,
                 std::uint64_t seed) noexcept {
  double result = 1.0;
  for (std::size_t index = 0; index < transforms.size(); ++index) {
    const auto &transform = transforms[index];
    if (!TransformEnabled(transform, cycle, seed, 0x9000 + index))
      continue;
    if (transform.kind == TransformKind::Fast)
      result /= transform.number;
    else if (transform.kind == TransformKind::Slow)
      result *= transform.number;
  }
  return std::isfinite(result) && result > 0.0 ? result : 1.0;
}

double LaneRate(const std::vector<Transform> &transforms, std::uint64_t cycle,
                std::uint64_t seed) noexcept {
  double result = 1.0;
  for (std::size_t index = 0; index < transforms.size(); ++index) {
    const auto &transform = transforms[index];
    if (transform.kind == TransformKind::Rate &&
        TransformEnabled(transform, cycle, seed, 0x9050 + index))
      result *= transform.number;
  }
  return std::isfinite(result) && result > 0.0 ? result : 1.0;
}

struct SwingSettings {
  double ratio = 0.5;
  double subdivisionBeats = 0.0;
};

SwingSettings Swing(const std::vector<Transform> &transforms,
                    std::uint64_t cycle, std::uint64_t seed) noexcept {
  SwingSettings result;
  for (std::size_t index = 0; index < transforms.size(); ++index) {
    const auto &transform = transforms[index];
    if (transform.kind != TransformKind::Swing ||
        !TransformEnabled(transform, cycle, seed, 0x9100 + index))
      continue;
    result.ratio = transform.number;
    result.subdivisionBeats = transform.swingSubdivisionBeats;
  }
  return result;
}

double TimingRandom(std::uint64_t seed, std::uint64_t cycle,
                    const Transform &transform,
                    std::size_t eventIndex) noexcept {
  std::uint64_t hash = seed ^ (cycle * 0x9e3779b97f4a7c15ULL) ^
                       (eventIndex * 0xbf58476d1ce4e5b9ULL) ^
                       static_cast<std::uint64_t>(transform.span.begin + 1);
  hash ^= hash >> 30;
  hash *= 0xbf58476d1ce4e5b9ULL;
  hash ^= hash >> 27;
  hash *= 0x94d049bb133111ebULL;
  hash ^= hash >> 31;
  return static_cast<double>(hash >> 11) / 9007199254740992.0;
}

void ApplyTimingOffsets(RuntimeEvent &event,
                        const std::vector<Transform> &transforms,
                        std::uint64_t cycle, std::uint64_t seed,
                        std::size_t eventIndex) noexcept {
  for (std::size_t index = 0; index < transforms.size(); ++index) {
    const auto &transform = transforms[index];
    if ((transform.kind != TransformKind::Early &&
         transform.kind != TransformKind::Late) ||
        !TransformEnabled(transform, cycle, seed, 0x9200 + index))
      continue;
    double amount = transform.number;
    if (transform.randomAmount)
      amount *= TimingRandom(seed, cycle, transform, eventIndex);
    const double direction =
        transform.kind == TransformKind::Early ? -1.0 : 1.0;
    if (transform.timeUnit == TimeUnit::Milliseconds)
      event.timingOffsetMilliseconds += direction * amount;
    else
      event.beat += direction * amount;
  }
}

std::int64_t ScaleDegreeInterval(const std::vector<int> &intervals,
                                 int degree) noexcept {
  if (degree == 0 || intervals.empty())
    return 0;
  const int direction = degree < 0 ? -1 : 1;
  const std::int64_t zeroBased = std::abs(degree) - 1LL;
  const auto scaleSize = static_cast<std::int64_t>(intervals.size());
  const auto octaves = zeroBased / scaleSize;
  const auto index = static_cast<std::size_t>(zeroBased % scaleSize);
  return direction * (12 * octaves + intervals[index]);
}

std::uint64_t MixRandom(std::uint64_t value) noexcept {
  value ^= value >> 30;
  value *= 0xbf58476d1ce4e5b9ULL;
  value ^= value >> 27;
  value *= 0x94d049bb133111ebULL;
  value ^= value >> 31;
  return value;
}

double RandomUnit(std::uint64_t seed, std::uint64_t key, std::uint64_t identity,
                  std::uint64_t salt, std::uint64_t draw = 0) noexcept {
  const auto bits = MixRandom(seed ^ MixRandom(key + 0x9e3779b97f4a7c15ULL) ^
                              MixRandom(identity + 0xa0761d6478bd642fULL) ^
                              salt ^ MixRandom(draw + 0xd1b54a32d192ed03ULL));
  // The half-step keeps both logarithm inputs away from zero and one.
  return (static_cast<double>(bits >> 11) + 0.5) / 9007199254740992.0;
}

double RandomNormal(std::uint64_t seed, std::uint64_t key,
                    std::uint64_t identity, std::uint64_t salt) noexcept {
  constexpr double TwoPi = 6.283185307179586476925286766559;
  const double radius =
      std::sqrt(-2.0 * std::log(RandomUnit(seed, key, identity, salt, 0)));
  return std::clamp(
      radius * std::cos(TwoPi * RandomUnit(seed, key, identity, salt, 1)),
      -NormalDeviationLimit, NormalDeviationLimit);
}

double RandomInclusiveInteger(double low, double high, std::uint64_t seed,
                              std::uint64_t key, std::uint64_t identity,
                              std::uint64_t salt) noexcept {
  const long double first = low;
  const long double count = static_cast<long double>(high) - first + 1.0L;
  const long double offset = std::floor(
      static_cast<long double>(RandomUnit(seed, key, identity, salt)) * count);
  return static_cast<double>(
      std::min(first + offset, static_cast<long double>(high)));
}

double SampleScalarItem(const ScalarItem &item, std::uint64_t seed,
                        std::uint64_t key, std::uint64_t salt) noexcept {
  double value = item.value;
  if (item.randomDistribution == ScalarItem::RandomDistribution::Uniform) {
    value = item.randomInteger
                ? RandomInclusiveInteger(item.randomFirst, item.randomSecond,
                                         seed, key, item.randomIdentity, salt)
                : item.randomFirst +
                      (item.randomSecond - item.randomFirst) *
                          RandomUnit(seed, key, item.randomIdentity, salt);
  } else if (item.randomDistribution ==
             ScalarItem::RandomDistribution::Normal) {
    value =
        item.randomFirst +
        item.randomSecond * RandomNormal(seed, key, item.randomIdentity, salt);
    if (item.randomInteger)
      value = std::round(value);
  }
  return std::clamp(value, item.randomMinimum, item.randomMaximum);
}

PitchValue SampleRandomPitch(const PitchItem &item, const Sequence &sequence,
                             std::uint64_t seed, std::uint64_t key,
                             std::uint64_t salt) noexcept {
  double sampled = 1.0;
  if (item.randomDefaultRange) {
    sampled = RandomInclusiveInteger(
        1.0, static_cast<double>(sequence.scale.intervals.size()), seed, key,
        item.randomIdentity, salt);
  } else if (item.randomDistribution ==
             PitchItem::RandomDistribution::Uniform) {
    sampled = RandomInclusiveInteger(item.randomFirst, item.randomSecond, seed,
                                     key, item.randomIdentity, salt);
  } else {
    sampled = std::round(
        item.randomFirst +
        item.randomSecond * RandomNormal(seed, key, item.randomIdentity, salt));
  }
  sampled = std::clamp(sampled,
                       static_cast<double>(std::numeric_limits<int>::min() + 1),
                       static_cast<double>(std::numeric_limits<int>::max()));
  const auto integer = static_cast<int>(sampled);
  PitchValue pitch;
  pitch.span = item.span;
  if (item.randomDomain == PitchItem::RandomDomain::ScaleDegree) {
    pitch.degree = integer;
    return pitch;
  }

  const std::int64_t absolute =
      static_cast<std::int64_t>(sequence.scale.tonicSemitone) + integer;
  std::int64_t octave = absolute / 12;
  std::int64_t pitchClass = absolute % 12;
  if (pitchClass < 0) {
    pitchClass += 12;
    --octave;
  }
  pitch.absolute = true;
  pitch.pitchClass = static_cast<int>(pitchClass);
  pitch.octaveOffset = static_cast<int>(
      std::clamp<std::int64_t>(octave, std::numeric_limits<int>::min(),
                               std::numeric_limits<int>::max()));
  return pitch;
}

std::int64_t ModulateDegrees(const std::vector<Transform> &transforms,
                             const Scale &scale, std::uint64_t cycle,
                             std::uint64_t seed, std::uint64_t salt) noexcept {
  std::int64_t result = 0;
  for (std::size_t index = 0; index < transforms.size(); ++index) {
    const auto &transform = transforms[index];
    if (transform.kind == TransformKind::ModulateDegree &&
        TransformEnabled(transform, cycle, seed, salt + index * 0x10001ULL))
      result = SaturatingAdd(
          result, ScaleDegreeInterval(transform.modulationIntervals.empty()
                                          ? scale.intervals
                                          : transform.modulationIntervals,
                                      transform.integer));
  }
  return result;
}

std::int64_t TransposeToKey(const std::vector<Transform> &transforms,
                            const Sequence &sequence, std::uint64_t cycle,
                            std::uint64_t seed, std::uint64_t salt) noexcept {
  std::int64_t result = 0;
  for (std::size_t index = 0; index < transforms.size(); ++index) {
    const auto &transform = transforms[index];
    if (transform.kind != TransformKind::TransposeKey ||
        !TransformEnabled(transform, cycle, seed, salt + index * 0x10001ULL))
      continue;
    int interval = (transform.integer - sequence.keyPitchClass + 12) % 12;
    // A target key denotes the nearest conventional transposition. The exact
    // tritone is kept upward; octave transforms remain available explicitly.
    if (interval > 6)
      interval -= 12;
    result = SaturatingAdd(result, interval);
  }
  return result;
}

float PitchVolts(const PitchValue &pitch, const Sequence &sequence,
                 int sequenceOctave, std::int64_t degreeShift,
                 std::int64_t degreeTransposeSemitones,
                 std::int64_t semitoneTranspose,
                 std::int64_t octaveTranspose) noexcept {
  const long double chromaticTranspose =
      static_cast<long double>(degreeTransposeSemitones) +
      static_cast<long double>(semitoneTranspose) +
      12.0L * static_cast<long double>(octaveTranspose);
  const int baseOctave = pitch.hasOctave ? pitch.octave : sequenceOctave;
  if (pitch.absolute) {
    const long double midi =
        (static_cast<long double>(baseOctave) + pitch.octaveOffset + 1.0L) *
            12.0L +
        pitch.pitchClass + chromaticTranspose;
    return static_cast<float>((midi - 60.0L) / 12.0L);
  }
  const int scaleSize = static_cast<int>(sequence.scale.intervals.size());
  if (scaleSize == 0)
    return 0.f;
  const auto zeroBased =
      SaturatingAdd(static_cast<std::int64_t>(pitch.degree) - 1, degreeShift);
  auto octave = zeroBased / scaleSize;
  auto index = zeroBased % scaleSize;
  if (index < 0) {
    index += scaleSize;
    --octave;
  }
  const long double tonicMidi =
      (static_cast<long double>(baseOctave) + pitch.octaveOffset + 1.0L) *
          12.0L +
      sequence.scale.tonicSemitone;
  const long double midi =
      tonicMidi + 12.0L * static_cast<long double>(octave) +
      sequence.scale.intervals[static_cast<std::size_t>(index)] +
      pitch.accidental + chromaticTranspose;
  return static_cast<float>((midi - 60.0L) / 12.0L);
}

int BoundedSemitone(const float pitchVolts) noexcept {
  constexpr int Margin = 64;
  const double semitones = static_cast<double>(pitchVolts) * 12.0;
  if (!std::isfinite(semitones))
    return std::signbit(semitones) ? std::numeric_limits<int>::min() + Margin
                                   : std::numeric_limits<int>::max() - Margin;
  const double bounded = std::clamp(
      semitones, static_cast<double>(std::numeric_limits<int>::min() + Margin),
      static_cast<double>(std::numeric_limits<int>::max() - Margin));
  return static_cast<int>(std::lround(bounded));
}

template <typename T>
const T &Cycled(const std::vector<T> &items, std::uint64_t cursor,
                const std::vector<Transform> &transforms, std::uint64_t cycle,
                std::uint64_t seed, std::uint64_t salt) noexcept {
  const auto index =
      TransformIndex(static_cast<std::size_t>(cursor), items.size(), transforms,
                     cycle, seed, salt);
  return items[index];
}

struct PitchTimelineSelection {
  const PitchTimelineStep *step = nullptr;
  std::size_t logicalIndex = 0;
};

PitchTimelineSelection
SelectPitchTimeline(const Sequence &sequence, double phraseBeat,
                    const std::vector<Transform> &transforms,
                    std::uint64_t cycle, std::uint64_t seed) noexcept {
  PitchTimelineSelection result;
  if (sequence.pitchTimeline.empty())
    return result;
  phraseBeat = std::clamp(phraseBeat, 0.0,
                          std::max(0.0, sequence.pitchTimelineBeats - 1e-12));
  double cursor = 0.0;
  for (std::size_t logical = 0; logical < sequence.pitchTimeline.size();
       ++logical) {
    const auto mapped = TransformIndex(logical, sequence.pitchTimeline.size(),
                                       transforms, cycle, seed, 0x3f00);
    const auto &step = sequence.pitchTimeline[mapped];
    const double span = sequence.subdivisionBeats * step.durationMultiplier;
    if (phraseBeat < cursor + span ||
        logical + 1 == sequence.pitchTimeline.size()) {
      result.step = &step;
      result.logicalIndex = logical;
      return result;
    }
    cursor += span;
  }
  return result;
}

bool AlignedItemIndex(std::size_t itemCount, LaneAlignment alignment,
                      std::size_t alignmentSplit, std::uint64_t structuralCell,
                      std::uint64_t structuralCellCount,
                      std::size_t &index) noexcept {
  if (alignment == LaneAlignment::Free || itemCount == 0 ||
      structuralCellCount == 0 || structuralCell >= structuralCellCount)
    return false;
  if (alignment == LaneAlignment::Left) {
    if (structuralCell >= itemCount)
      return false;
    index = static_cast<std::size_t>(structuralCell);
    return true;
  }
  if (alignment == LaneAlignment::Edges) {
    const auto leftCount = std::min(alignmentSplit, itemCount);
    const auto rightCount = itemCount - leftCount;
    const bool hasLeft = structuralCell < leftCount;
    const auto visibleRight =
        std::min<std::uint64_t>(rightCount, structuralCellCount);
    const auto rightStartCell = structuralCellCount - visibleRight;
    const bool hasRight = rightCount > 0 && structuralCell >= rightStartCell;

    // Presence omissions can make the realized Notes pass shorter than the
    // two explicit edge groups together. Preserve the values nearest each
    // outer edge; where the groups overlap, the closest edge owns the cell.
    // A central tie goes to the left group for a stable deterministic result.
    const bool chooseLeft =
        hasLeft && (!hasRight ||
                    structuralCell <= structuralCellCount - 1 - structuralCell);
    if (chooseLeft) {
      index = static_cast<std::size_t>(structuralCell);
      return true;
    }
    if (hasRight) {
      const auto omittedRight = rightCount - visibleRight;
      index = leftCount + static_cast<std::size_t>(
                              omittedRight + structuralCell - rightStartCell);
      return true;
    }
    return false;
  }
  if (structuralCellCount >= itemCount) {
    const auto padding = structuralCellCount - itemCount;
    if (structuralCell < padding)
      return false;
    index = static_cast<std::size_t>(structuralCell - padding);
    return true;
  }
  index = itemCount - static_cast<std::size_t>(structuralCellCount) +
          static_cast<std::size_t>(structuralCell);
  return true;
}

double Scalar(const std::vector<ScalarItem> &items, std::uint64_t &cursor,
              const std::vector<Transform> &transforms, std::uint64_t cycle,
              std::uint64_t seed, double fallback, SourceSpan &span,
              std::uint64_t salt, LaneAlignment alignment = LaneAlignment::Free,
              std::uint64_t structuralCell = 0,
              std::uint64_t structuralCellCount = 0,
              std::size_t alignmentSplit = 0, SourceSpan alignmentSpan = {},
              bool *milliseconds = nullptr, double scoreBeat = 0.0,
              bool *noteValue = nullptr) noexcept {
  if (milliseconds)
    *milliseconds = false;
  if (noteValue)
    *noteValue = false;
  if (items.empty())
    return fallback;
  const bool phaseRate = std::any_of(
      transforms.begin(), transforms.end(), [](const Transform &transform) {
        return transform.kind == TransformKind::Rate;
      });
  std::uint64_t position = 0;
  std::uint64_t randomKey = 0;
  if (alignment != LaneAlignment::Free) {
    const auto mappedCell =
        TransformIndex(static_cast<std::size_t>(structuralCell),
                       static_cast<std::size_t>(structuralCellCount),
                       transforms, cycle, seed, salt);
    std::size_t alignedIndex = 0;
    if (!AlignedItemIndex(items.size(), alignment, alignmentSplit, mappedCell,
                          structuralCellCount, alignedIndex)) {
      span = alignmentSpan;
      return fallback;
    }
    randomKey = structuralCell;
    const auto &item = items[alignedIndex];
    span = item.span;
    if (milliseconds)
      *milliseconds = item.isMilliseconds;
    if (noteValue)
      *noteValue = item.isNoteValue;
    randomKey ^= MixRandom(cycle);
    return item.isDefault ? fallback
                          : SampleScalarItem(item, seed, randomKey, salt);
  } else if (phaseRate) {
    const double whole = std::max(
        0.0, std::floor(scoreBeat * LaneRate(transforms, cycle, seed)));
    position = static_cast<std::uint64_t>(
        std::fmod(whole, static_cast<double>(items.size())));
    randomKey = static_cast<std::uint64_t>(whole);
  } else {
    randomKey = cursor;
    position = cursor++;
  }
  randomKey ^= MixRandom(cycle);
  const auto &item = Cycled(items, position, transforms, cycle, seed, salt);
  span = item.span;
  if (milliseconds)
    *milliseconds = item.isMilliseconds;
  if (noteValue)
    *noteValue = item.isNoteValue;
  return item.isDefault ? fallback
                        : SampleScalarItem(item, seed, randomKey, salt);
}

struct CvSample {
  float value = 0.f;
  float target = 0.f;
  double targetBeat = 0.0;
  SourceSpan span;
};

CvSample SampleCv(const std::vector<ScalarItem> &items,
                  const std::vector<Transform> &transforms,
                  LaneAlignment alignment, std::uint64_t structuralCell,
                  std::uint64_t structuralCellCount, std::size_t alignmentSplit,
                  SourceSpan alignmentSpan, double scoreBeat, double cycleBeat,
                  double cellSpan, CvInterpolation interpolation, double power,
                  std::uint64_t cycle, std::uint64_t seed,
                  std::uint64_t salt) noexcept {
  CvSample sample;
  sample.targetBeat = scoreBeat + cellSpan;
  if (items.empty())
    return sample;
  if (alignment != LaneAlignment::Free) {
    sample.span = alignmentSpan;
    if (structuralCellCount == 0)
      return sample;
    const auto virtualCount = structuralCellCount;
    const auto base = structuralCell % virtualCount;
    if (alignment == LaneAlignment::Edges) {
      const auto mappedCell =
          TransformIndex(static_cast<std::size_t>(base),
                         static_cast<std::size_t>(virtualCount), transforms,
                         cycle, seed, salt);
      std::size_t itemIndex = 0;
      if (!AlignedItemIndex(items.size(), alignment, alignmentSplit, mappedCell,
                            virtualCount, itemIndex))
        return sample;
      const auto &item = items[itemIndex];
      sample.span = item.span;
      if (item.isDefault)
        return sample;
      const auto value = static_cast<float>(SampleScalarItem(
          item, seed, structuralCell ^ MixRandom(cycle), salt));
      sample.value = value;
      sample.target = value;
      return sample;
    }
    const ScalarItem *previous = nullptr;
    std::uint64_t previousCell = base;
    for (std::uint64_t distance = 0; distance < virtualCount; ++distance) {
      const auto cell = (base + virtualCount - distance) % virtualCount;
      const auto mappedCell =
          TransformIndex(static_cast<std::size_t>(cell),
                         static_cast<std::size_t>(virtualCount), transforms,
                         cycle, seed, salt);
      std::size_t itemIndex = 0;
      if (!AlignedItemIndex(items.size(), alignment, alignmentSplit, mappedCell,
                            virtualCount, itemIndex))
        continue;
      const auto &candidate = items[itemIndex];
      if (distance == 0)
        sample.span = candidate.span;
      if (!candidate.isDefault) {
        previous = &candidate;
        previousCell = cell;
        break;
      }
    }
    if (!previous)
      return sample;
    const auto key = previousCell ^ MixRandom(cycle);
    sample.value =
        static_cast<float>(SampleScalarItem(*previous, seed, key, salt));
    sample.target = sample.value;
    return sample;
  }
  const double rate = LaneRate(transforms, cycle, seed);
  const double phase = std::max(0.0, cycleBeat * rate);
  const auto base = static_cast<std::size_t>(
      std::fmod(std::floor(phase), static_cast<double>(items.size())));
  auto at = [&](std::size_t distance, bool backward) -> const ScalarItem & {
    const auto size = items.size();
    const auto origin = base;
    const auto shift = distance % size;
    const auto position =
        backward ? (origin >= shift ? origin - shift : size - (shift - origin))
                 : (shift < size - origin ? origin + shift
                                          : shift - (size - origin));
    return Cycled(items, position, transforms, cycle, seed, salt);
  };
  sample.span = at(0, false).span;

  std::size_t previousDistance = items.size();
  for (std::size_t distance = 0; distance < items.size(); ++distance) {
    if (!at(distance, true).isDefault) {
      previousDistance = distance;
      break;
    }
  }
  if (previousDistance == items.size())
    return sample;
  const auto &previous = at(previousDistance, true);
  auto randomKey = [cycle](double knot) {
    const auto bounded = std::clamp(
        knot, static_cast<double>(std::numeric_limits<std::int64_t>::min()),
        static_cast<double>(std::numeric_limits<std::int64_t>::max()));
    return static_cast<std::uint64_t>(static_cast<std::int64_t>(bounded)) ^
           MixRandom(cycle);
  };
  const double previousPhase =
      std::floor(phase) - static_cast<double>(previousDistance);
  const double previousValue =
      SampleScalarItem(previous, seed, randomKey(previousPhase), salt);
  sample.value = static_cast<float>(previousValue);
  sample.target = sample.value;
  if (interpolation == CvInterpolation::Step)
    return sample;

  std::size_t nextDistance = items.size();
  for (std::size_t distance = 1; distance <= items.size(); ++distance) {
    if (!at(distance, false).isDefault) {
      nextDistance = distance;
      break;
    }
  }
  if (nextDistance > items.size())
    return sample;
  const auto &next = at(nextDistance, false);
  const double nextPhase =
      std::floor(phase) + static_cast<double>(nextDistance);
  const double nextValue =
      SampleScalarItem(next, seed, randomKey(nextPhase), salt);
  sample.target = static_cast<float>(nextValue);
  double amount = (phase - previousPhase) / (nextPhase - previousPhase);
  amount = std::clamp(amount, 0.0, 1.0);
  if (interpolation == CvInterpolation::Smooth)
    amount = amount * amount * (3.0 - 2.0 * amount);
  else if (interpolation == CvInterpolation::Power)
    amount = std::pow(amount, power);
  sample.value =
      static_cast<float>(previousValue + (nextValue - previousValue) * amount);
  sample.targetBeat = scoreBeat + (nextPhase - phase) / rate;
  return sample;
}

ArticulationKind EffectiveArticulation(const ArticulationAtom &atom,
                                       std::uint64_t articulationCursor,
                                       std::size_t atomIndex,
                                       std::uint64_t seed) noexcept {
  auto passes = [&](float probability, std::uint64_t group,
                    std::uint64_t level) {
    std::uint64_t probabilityHash =
        seed ^ (articulationCursor * 0x9e3779b9ULL) ^
        ((group != 0 ? group : atomIndex + 1) * 0x85ebca6bULL) ^
        (level * 0xc2b2ae35ULL);
    probabilityHash ^= probabilityHash >> 33;
    const float random = static_cast<float>(probabilityHash & 0xffffffU) /
                         static_cast<float>(0x1000000U);
    return random < probability;
  };
  if (!passes(atom.probability, atom.probabilityGroup, 1))
    return ArticulationKind::Rest;
  for (const auto &gate : atom.enclosingProbabilityGates) {
    // A prepared group id is unique per expanded group instance. Keeping the
    // hash level fixed makes one outer group decision coherent even when some
    // of its members also belong to nested groups.
    if (!passes(gate.probability, gate.group, 2))
      return ArticulationKind::Rest;
  }
  return atom.kind;
}

bool StepIsPresent(const ArticulationStep &step, std::uint64_t cycle,
                   std::uint64_t seed) noexcept {
  return step.presenceProbability >= 1.f ||
         RandomUnit(seed, cycle, step.presenceIdentity, 0xb170U) <
             step.presenceProbability;
}

std::uint64_t SurvivingStructuralCellCount(const Sequence &sequence,
                                           std::uint64_t cycle,
                                           std::uint64_t seed) noexcept {
  std::uint64_t result = 0;
  for (const auto &step : sequence.articulation) {
    if (StepIsPresent(step, cycle, seed))
      result += static_cast<std::uint64_t>(step.cellCount);
  }
  return result;
}

std::size_t PreparedVoiceCount(const PitchItem &item,
                               const VoicingStyle style) noexcept {
  if (item.randomDomain != PitchItem::RandomDomain::None)
    return 1;
  std::size_t result = 0;
  for (const auto &choice : item.values) {
    const auto harmonicVoices =
        choice.meaning == ChordValue::Meaning::JazzSymbol ||
                choice.meaning == ChordValue::Meaning::RomanSymbol
            ? MaximumVoicingCount(choice, style)
            : choice.voices.size();
    const auto voices = harmonicVoices + (choice.hasBass ? 1u : 0u);
    if (voices == 0 || (result != 0 && result != voices))
      return 0;
    result = voices;
  }
  return result;
}

} // namespace

double SchedulingLookaheadBeats(const CompiledProgram &program,
                                bool periodKnown, double periodSamples,
                                double sampleRateHz) noexcept {
  double millisecondsAsBeats = 0.0;
  if (periodKnown && periodSamples > 0.0 && sampleRateHz > 0.0) {
    millisecondsAsBeats =
        program.maximumEarlyMilliseconds * 0.001 * sampleRateHz / periodSamples;
  }
  const double lookahead =
      1.0 + program.maximumEarlyBeats + millisecondsAsBeats;
  return std::isfinite(lookahead) && lookahead >= 1.0 ? lookahead : 1.0;
}

void Runtime::setProgram(const CompiledProgram *program) noexcept {
  program_ = program;
  states_ = program ? program->stateWorkspace.data() : nullptr;
  stateCount_ = program ? program->stateWorkspace.size() : 0;
  reset();
}

bool Runtime::canPreserveCurrentPhase(
    const CompiledProgram *program) const noexcept {
  if (!program_ || !program)
    return false;
  const auto &previousSemantic = program_->semantic();
  const auto &nextSemantic = program->semantic();
  if (arrangementPart_ >= previousSemantic.arrangement.size())
    return false;
  const auto previousSequenceIndex =
      previousSemantic.arrangement[arrangementPart_].sequence;
  if (previousSequenceIndex >= previousSemantic.sequences.size())
    return false;
  const auto &previousSequence =
      previousSemantic.sequences[previousSequenceIndex];
  const auto sameSequence = [&](const SemanticProgram &semantic,
                                const ArrangementPart &part) {
    if (part.sequence >= semantic.sequences.size())
      return false;
    const auto &candidate = semantic.sequences[part.sequence];
    return candidate.stableId == previousSequence.stableId &&
           candidate.name == previousSequence.name;
  };
  std::size_t occurrence = 0;
  for (std::size_t part = 0; part < arrangementPart_; ++part)
    if (sameSequence(previousSemantic, previousSemantic.arrangement[part]))
      ++occurrence;
  std::size_t nextOccurrence = 0;
  for (const auto &part : nextSemantic.arrangement) {
    if (!sameSequence(nextSemantic, part))
      continue;
    if (nextOccurrence++ != occurrence)
      continue;
    const auto &nextSequence = nextSemantic.sequences[part.sequence];
    return nextSequence.separateRhythm == previousSequence.separateRhythm;
  }
  return false;
}

Runtime::ReplacementResult
Runtime::replaceProgram(const CompiledProgram *program, double beat) noexcept {
  const bool activePhaseCompatible = canPreserveCurrentPhase(program);
  const auto previousArrangementPart = arrangementPart_;
  const auto previousArrangementCycle = arrangementCycle_;
  const auto previousPartStartBeat = partStartBeat_;
  const auto *previousStates = states_;
  const auto previousStateCount = stateCount_;
  const auto *previousProgram = program_;
  auto *nextStates = program ? program->stateWorkspace.data() : nullptr;
  const auto nextStateCount = program ? program->stateWorkspace.size() : 0;
  if (nextStates)
    std::fill_n(nextStates, nextStateCount, SequencePlaybackState{});
  if (program && previousProgram && nextStates && previousStates) {
    const auto &nextSequences = program->semantic().sequences;
    const auto &previousSequences = previousProgram->semantic().sequences;
    const auto &nextOrder = program->stateTransferOrder;
    const auto &previousOrder = previousProgram->stateTransferOrder;
    std::size_t nextPosition = 0;
    std::size_t previousPosition = 0;
    while (nextPosition < nextOrder.size() &&
           previousPosition < previousOrder.size()) {
      const auto nextIndex = nextOrder[nextPosition];
      const auto previousIndex = previousOrder[previousPosition];
      if (nextIndex >= nextStateCount || previousIndex >= previousStateCount ||
          nextIndex >= nextSequences.size() ||
          previousIndex >= previousSequences.size())
        break;
      const auto &nextSequence = nextSequences[nextIndex];
      const auto &previousSequence = previousSequences[previousIndex];
      const bool nextBeforePrevious =
          nextSequence.stableId < previousSequence.stableId ||
          (nextSequence.stableId == previousSequence.stableId &&
           nextSequence.name < previousSequence.name);
      const bool previousBeforeNext =
          previousSequence.stableId < nextSequence.stableId ||
          (previousSequence.stableId == nextSequence.stableId &&
           previousSequence.name < nextSequence.name);
      if (nextBeforePrevious) {
        ++nextPosition;
      } else if (previousBeforeNext) {
        ++previousPosition;
      } else {
        if (nextSequence.separateRhythm == previousSequence.separateRhythm) {
          nextStates[nextIndex] = previousStates[previousIndex];
          // Presence decisions and the realized aligned-lane length belong to
          // the newly compiled structure. Recompute the cached count on its
          // first event after activation.
          nextStates[nextIndex].structuralCellCount = 0;
          if (nextSequence.voicing != previousSequence.voicing ||
              nextSequence.hasKey != previousSequence.hasKey ||
              nextSequence.keyPitchClass != previousSequence.keyPitchClass) {
            nextStates[nextIndex].previousVoicingCount = 0;
          }
        }
        ++nextPosition;
        ++previousPosition;
      }
    }
  }

  program_ = program;
  states_ = nextStates;
  stateCount_ = nextStateCount;
  const auto arrangementSize =
      program_ ? program_->semantic().arrangement.size() : 0;
  if (arrangementSize == 0) {
    arrangementPart_ = 0;
    arrangementCycle_ = 0;
    partStartBeat_ = beat;
    return ReplacementResult::RestartedCurrentTerm;
  }

  // Arrangement indices are not identities: inserting or reordering a term
  // before the current one must not jump playback to an unrelated section.
  // Match the active named sequence and its occurrence among duplicate terms.
  bool transferredPart = false;
  if (previousProgram) {
    const auto &previousSemantic = previousProgram->semantic();
    const auto &previousArrangement = previousSemantic.arrangement;
    if (previousArrangementPart < previousArrangement.size()) {
      const auto previousSequenceIndex =
          previousArrangement[previousArrangementPart].sequence;
      if (previousSequenceIndex < previousSemantic.sequences.size()) {
        const auto &previousSequence =
            previousSemantic.sequences[previousSequenceIndex];
        const auto sameSequence = [&](const SemanticProgram &semantic,
                                      const ArrangementPart &part) {
          if (part.sequence >= semantic.sequences.size())
            return false;
          const auto &candidate = semantic.sequences[part.sequence];
          return candidate.stableId == previousSequence.stableId &&
                 candidate.name == previousSequence.name;
        };

        std::size_t previousOccurrence = 0;
        for (std::size_t part = 0; part < previousArrangementPart; ++part) {
          if (sameSequence(previousSemantic, previousArrangement[part]))
            ++previousOccurrence;
        }

        const auto &nextSemantic = program_->semantic();
        std::size_t nextOccurrence = 0;
        for (std::size_t part = 0; part < arrangementSize; ++part) {
          if (!sameSequence(nextSemantic, nextSemantic.arrangement[part]))
            continue;
          if (nextOccurrence++ != previousOccurrence)
            continue;
          arrangementPart_ = part;
          transferredPart = true;
          break;
        }
      }
    }
  }

  if (!transferredPart) {
    // The active term was removed. Keep the nearest structural location, but
    // start that replacement term cleanly at the activation beat.
    arrangementPart_ = std::min(previousArrangementPart, arrangementSize - 1);
    arrangementCycle_ = 0;
    partStartBeat_ = beat;
    return ReplacementResult::RestartedCurrentTerm;
  }

  const int cycles = program_->semantic().arrangement[arrangementPart_].cycles;
  arrangementCycle_ = cycles < 0
                          ? std::max(0, previousArrangementCycle)
                          : std::clamp(previousArrangementCycle, 0, cycles - 1);
  partStartBeat_ = previousPartStartBeat;
  if (!activePhaseCompatible) {
    const auto activeSequence =
        program_->semantic().arrangement[arrangementPart_].sequence;
    if (activeSequence < stateCount_)
      states_[activeSequence] = {};
    partStartBeat_ = beat;
    return ReplacementResult::RestartedCurrentTerm;
  }
  return ReplacementResult::PreservedPhase;
}

void Runtime::reset() noexcept {
  if (states_)
    std::fill_n(states_, stateCount_, SequencePlaybackState{});
  arrangementPart_ = 0;
  arrangementCycle_ = 0;
  partStartBeat_ = 0.0;
}

void Runtime::snapshot(Runtime &destination) const noexcept {
  destination = *this;
  if (!program_) {
    destination.states_ = nullptr;
    destination.stateCount_ = 0;
    return;
  }
  auto &snapshot = program_->snapshotStateWorkspace;
  const auto count = std::min(stateCount_, snapshot.size());
  if (states_ && count > 0)
    std::copy_n(states_, count, snapshot.data());
  destination.states_ = snapshot.data();
  destination.stateCount_ = count;
}

const Sequence *Runtime::currentSequence(SourceSpan &partSpan,
                                         std::uint64_t &cycle) noexcept {
  if (!program_ || program_->semantic().arrangement.empty())
    return nullptr;
  const auto &semantic = program_->semantic();
  if (arrangementPart_ >= semantic.arrangement.size())
    arrangementPart_ = 0;
  const auto &part = semantic.arrangement[arrangementPart_];
  if (part.sequence >= semantic.sequences.size())
    return nullptr;
  const auto &sequence = semantic.sequences[part.sequence];
  cycle =
      part.sequence < stateCount_ ? states_[part.sequence].completedCycles : 0;
  partSpan = part.span.valid() ? part.span : sequence.nameSpan;
  return &sequence;
}

StepEvents Runtime::next(double beat) noexcept {
  StepEvents output;
  if (!program_)
    return output;
  output.events = program_->stepWorkspace.data();
  output.capacity = program_->stepWorkspace.size();
  const auto &semantic = program_->semantic();
  const auto seed = semantic.seed;
  SourceSpan partSpan;
  std::uint64_t cycle = 0;
  const Sequence *sequence = nullptr;
  SequencePlaybackState *statePointer = nullptr;
  const ArticulationStep *step = nullptr;
  std::size_t sequenceIndex = 0;

  // Presence probability removes a top-level step before Duration or any
  // event-indexed lane is sampled. Walk across zero-time omissions here so the
  // caller always receives a positive-duration step.
  for (;;) {
    sequence = currentSequence(partSpan, cycle);
    if (!sequence || sequence->articulation.empty())
      return output;
    sequenceIndex =
        static_cast<std::size_t>(sequence - semantic.sequences.data());
    if (sequenceIndex >= stateCount_)
      return output;
    statePointer = &states_[sequenceIndex];
    auto &candidateState = *statePointer;
    if (sequence->separateRhythm && candidateState.pitchTimelineBeat >=
                                        sequence->pitchTimelineBeats - 1e-9) {
      const auto &part = semantic.arrangement[arrangementPart_];
      bool continuesSameSequence = false;
      if (part.cycles < 0 || arrangementCycle_ + 1 < part.cycles) {
        continuesSameSequence = true;
      } else {
        const auto nextPart =
            (arrangementPart_ + 1) % semantic.arrangement.size();
        continuesSameSequence =
            semantic.arrangement[nextPart].sequence == sequenceIndex;
      }
      const auto completedCycles = candidateState.completedCycles + 1;
      if (continuesSameSequence) {
        candidateState.completedCycles = completedCycles;
        candidateState.freeLaneBeat += std::max(0.0, beat - partStartBeat_);
        candidateState.articulation = 0;
        candidateState.structuralCell = 0;
        candidateState.structuralCellCount = 0;
        candidateState.pitchTimelineBeat = 0.0;
        candidateState.activePitchTimelineStep =
            std::numeric_limits<std::size_t>::max();
        candidateState.pitchEntryConsumed = false;
      } else {
        const auto lastBaseDuration = candidateState.lastBaseDuration;
        candidateState = {};
        candidateState.completedCycles = completedCycles;
        candidateState.lastBaseDuration = lastBaseDuration;
      }
      partStartBeat_ = beat;
      if (part.cycles >= 0) {
        ++arrangementCycle_;
        if (arrangementCycle_ >= part.cycles) {
          arrangementCycle_ = 0;
          arrangementPart_ =
              (arrangementPart_ + 1) % semantic.arrangement.size();
        }
      }
      continue;
    }
    const auto articulationLane =
        sequence->separateRhythm ? CursorLane::Rhythm : CursorLane::Notes;
    step = &Cycled(
        sequence->articulation, candidateState.articulation++,
        sequence->transforms[static_cast<std::size_t>(articulationLane)], cycle,
        seed, 0x3000);
    const bool completesPass =
        candidateState.articulation % sequence->articulation.size() == 0;
    if (StepIsPresent(*step, cycle, seed)) {
      if (candidateState.structuralCellCount == 0) {
        candidateState.structuralCellCount =
            SurvivingStructuralCellCount(*sequence, cycle, seed);
      }
      break;
    }
    if (sequence->separateRhythm)
      continue;
    if (!completesPass)
      continue;

    const auto &part = semantic.arrangement[arrangementPart_];
    bool continuesSameSequence = false;
    if (part.cycles < 0 || arrangementCycle_ + 1 < part.cycles) {
      continuesSameSequence = true;
    } else {
      const auto nextPart =
          (arrangementPart_ + 1) % semantic.arrangement.size();
      continuesSameSequence =
          semantic.arrangement[nextPart].sequence == sequenceIndex;
    }
    const auto completedCycles = candidateState.completedCycles + 1;
    if (continuesSameSequence) {
      candidateState.completedCycles = completedCycles;
      candidateState.freeLaneBeat += std::max(0.0, beat - partStartBeat_);
      candidateState.structuralCell = 0;
      candidateState.structuralCellCount = 0;
    } else {
      const auto lastBaseDuration = candidateState.lastBaseDuration;
      candidateState = {};
      candidateState.completedCycles = completedCycles;
      candidateState.lastBaseDuration = lastBaseDuration;
    }
    partStartBeat_ = beat;
    if (part.cycles >= 0) {
      ++arrangementCycle_;
      if (arrangementCycle_ >= part.cycles) {
        arrangementCycle_ = 0;
        arrangementPart_ = (arrangementPart_ + 1) % semantic.arrangement.size();
      }
    }
  }

  auto &state = *statePointer;
  const auto &timingTransforms =
      sequence->transforms[static_cast<std::size_t>(CursorLane::Sequence)];

  ArticulationAtom implicit;
  implicit.kind = ArticulationKind::Attack;
  implicit.span = sequence->nameSpan;

  SourceSpan durationSpan;
  const bool tieOnly = step && !step->atoms.empty() &&
                       std::all_of(step->atoms.begin(), step->atoms.end(),
                                   [](const ArticulationAtom &atom) {
                                     return atom.kind == ArticulationKind::Tie;
                                   });
  const double defaultDuration = sequence->separateRhythm
                                     ? sequence->rhythmSubdivisionBeats
                                     : sequence->subdivisionBeats;
  double baseDuration =
      state.lastBaseDuration > 0.0 ? state.lastBaseDuration : defaultDuration;
  if (!tieOnly) {
    baseDuration = static_cast<double>(Scalar(
        sequence->duration, state.duration,
        sequence->transforms[static_cast<std::size_t>(CursorLane::Duration)],
        cycle, seed, defaultDuration, durationSpan, 0x2000,
        sequence->alignment[static_cast<std::size_t>(CursorLane::Duration)],
        state.structuralCell, state.structuralCellCount,
        sequence
            ->alignmentSplits[static_cast<std::size_t>(CursorLane::Duration)],
        sequence
            ->alignmentSpans[static_cast<std::size_t>(CursorLane::Duration)],
        nullptr, state.freeLaneBeat + beat - partStartBeat_));
    state.lastBaseDuration = baseDuration;
  }
  const double unscaledDuration =
      baseDuration * (step ? step->durationMultiplier : 1.0);
  const double timeScale = TimeScale(timingTransforms, cycle, seed);
  const double duration = unscaledDuration * timeScale;
  const double usedUnscaledDuration =
      sequence->separateRhythm
          ? std::min(unscaledDuration,
                     std::max(0.0, sequence->pitchTimelineBeats -
                                       state.pitchTimelineBeat))
          : unscaledDuration;
  output.durationBeats = usedUnscaledDuration * timeScale;

  const bool completesNotesPass =
      sequence->separateRhythm
          ? state.pitchTimelineBeat + usedUnscaledDuration >=
                sequence->pitchTimelineBeats - 1e-9
          : !sequence->articulation.empty() &&
                state.articulation % sequence->articulation.size() == 0;
  const auto &arrangementPart = semantic.arrangement[arrangementPart_];
  bool nextPassSameSequence = false;
  if (arrangementPart.cycles < 0 ||
      arrangementCycle_ + 1 < arrangementPart.cycles) {
    nextPassSameSequence = true;
  } else {
    const auto nextPart = (arrangementPart_ + 1) % semantic.arrangement.size();
    nextPassSameSequence =
        semantic.arrangement[nextPart].sequence == sequenceIndex;
  }
  const bool continuesSameSequence =
      !completesNotesPass || nextPassSameSequence;

  const std::size_t atomCount = step ? step->atoms.size() : 1;
  bool legatoToNext = false;
  std::size_t legatoToNextVoiceCount = 0;
  bool nextPitchEntrySlides = false;
  std::size_t nextPitchTimelineIndex = std::numeric_limits<std::size_t>::max();
  if (step && !sequence->articulation.empty()) {
    const auto stepCount = sequence->articulation.size();
    std::uint64_t nextCursor = state.articulation;
    std::uint64_t nextCycle = cycle;
    for (std::size_t checked = 0; checked < stepCount; ++checked) {
      if (nextCursor % stepCount == 0) {
        if (!sequence->separateRhythm && !nextPassSameSequence)
          break;
        if (!sequence->separateRhythm)
          ++nextCycle;
      }
      const auto &nextStep =
          Cycled(sequence->articulation, nextCursor,
                 sequence->transforms[static_cast<std::size_t>(
                     sequence->separateRhythm ? CursorLane::Rhythm
                                              : CursorLane::Notes)],
                 nextCycle, seed, 0x3000);
      ++nextCursor;
      if (!StepIsPresent(nextStep, nextCycle, seed))
        continue;
      if (!nextStep.atoms.empty()) {
        const auto nextArticulationCursor = nextCursor % stepCount;
        const auto randomCursor =
            nextArticulationCursor == 0 ? stepCount : nextArticulationCursor;
        const auto nextKind =
            EffectiveArticulation(nextStep.atoms.front(),
                                  randomCursor ^ MixRandom(nextCycle), 0, seed);
        legatoToNext = nextKind == ArticulationKind::Tie ||
                       nextKind == ArticulationKind::Slide;
        if (nextKind == ArticulationKind::Slide &&
            nextStep.atoms.front().hasPitch)
          legatoToNextVoiceCount = PreparedVoiceCount(
              nextStep.atoms.front().pitch, sequence->voicing);
        if (sequence->separateRhythm && !completesNotesPass &&
            (nextKind == ArticulationKind::Attack ||
             nextKind == ArticulationKind::Slide)) {
          const double nextStepDuration =
              sequence->rhythmSubdivisionBeats * nextStep.durationMultiplier;
          const auto selection = SelectPitchTimeline(
              *sequence,
              state.pitchTimelineBeat + unscaledDuration +
                  nextStepDuration * nextStep.atoms.front().offsetFraction,
              sequence->transforms[static_cast<std::size_t>(CursorLane::Notes)],
              cycle, seed);
          if (selection.step &&
              selection.step->kind == PitchTimelineStep::Kind::Pitch &&
              selection.step->slideOnEntry) {
            nextPitchEntrySlides = true;
            nextPitchTimelineIndex = selection.logicalIndex;
            legatoToNextVoiceCount =
                PreparedVoiceCount(selection.step->pitch, sequence->voicing);
          }
        }
      }
      break;
    }
  }

  for (std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex) {
    const auto &atom = step ? step->atoms[atomIndex] : implicit;
    const auto articulationRandomKey = state.articulation ^ MixRandom(cycle);
    ArticulationKind kind =
        EffectiveArticulation(atom, articulationRandomKey, atomIndex, seed);
    // Probability is resolved at runtime, so predecessor validity must be
    // resolved here too. A missed source cannot feed a later tie; a slide with
    // no sounding source becomes a normal attack rather than an orphaned
    // legato event.
    if (kind == ArticulationKind::Tie && !state.hasSoundingPitch)
      kind = ArticulationKind::Rest;

    bool legatoFromAtom = false;
    std::size_t legatoTargetVoiceCount = 0;
    if ((kind == ArticulationKind::Attack || kind == ArticulationKind::Slide) &&
        !atom.ghost) {
      if (atomIndex + 1 < atomCount && step) {
        const auto nextKind =
            EffectiveArticulation(step->atoms[atomIndex + 1],
                                  articulationRandomKey, atomIndex + 1, seed);
        legatoFromAtom = nextKind == ArticulationKind::Tie ||
                         nextKind == ArticulationKind::Slide;
        if (nextKind == ArticulationKind::Slide &&
            step->atoms[atomIndex + 1].hasPitch)
          legatoTargetVoiceCount = PreparedVoiceCount(
              step->atoms[atomIndex + 1].pitch, sequence->voicing);
      } else {
        legatoFromAtom = legatoToNext;
        legatoTargetVoiceCount = legatoToNextVoiceCount;
      }
    }

    const double atomOffsetUnscaled = unscaledDuration * atom.offsetFraction;
    if (sequence->separateRhythm &&
        atomOffsetUnscaled >= usedUnscaledDuration - 1e-9)
      continue;
    const double atomBeat = beat + duration * atom.offsetFraction;
    const double cycleBeat = state.freeLaneBeat + atomBeat - partStartBeat_;
    double atomSpan = duration * atom.spanFraction;
    if (sequence->separateRhythm)
      atomSpan = std::min(atomSpan, beat + output.durationBeats - atomBeat);
    if (atomSpan <= 0.0)
      continue;
    PitchTimelineSelection timelineSelection;
    const PitchTimelineStep *timelineStep = nullptr;
    if (sequence->separateRhythm && (atom.kind == ArticulationKind::Attack ||
                                     atom.kind == ArticulationKind::Slide)) {
      timelineSelection = SelectPitchTimeline(
          *sequence, state.pitchTimelineBeat + atomOffsetUnscaled,
          sequence->transforms[static_cast<std::size_t>(CursorLane::Notes)],
          cycle, seed);
      timelineStep = timelineSelection.step;
      if (!timelineStep ||
          timelineStep->kind == PitchTimelineStep::Kind::Rest) {
        kind = ArticulationKind::Rest;
      } else {
        if (state.activePitchTimelineStep != timelineSelection.logicalIndex) {
          state.activePitchTimelineStep = timelineSelection.logicalIndex;
          state.pitchEntryConsumed = false;
        }
        if ((kind == ArticulationKind::Attack ||
             kind == ArticulationKind::Slide) &&
            timelineStep->slideOnEntry && !state.pitchEntryConsumed)
          kind = ArticulationKind::Slide;
        if (kind == ArticulationKind::Attack || kind == ArticulationKind::Slide)
          state.pitchEntryConsumed = true;
        if (atomIndex + 1 == atomCount && nextPitchEntrySlides &&
            nextPitchTimelineIndex != timelineSelection.logicalIndex) {
          legatoFromAtom = true;
          legatoTargetVoiceCount = legatoToNextVoiceCount;
        }
      }
    }
    if (kind == ArticulationKind::Slide && !state.hasSoundingPitch)
      kind = ArticulationKind::Attack;
    const auto structuralCell = state.structuralCell + atom.cellOffset;
    SourceSpan offsetSpan;
    bool offsetMilliseconds = false;
    const double offset = Scalar(
        sequence->offset, state.offset,
        sequence->transforms[static_cast<std::size_t>(CursorLane::Offset)],
        cycle, seed, 0.0, offsetSpan, 0x8800,
        sequence->alignment[static_cast<std::size_t>(CursorLane::Offset)],
        structuralCell, state.structuralCellCount,
        sequence->alignmentSplits[static_cast<std::size_t>(CursorLane::Offset)],
        sequence->alignmentSpans[static_cast<std::size_t>(CursorLane::Offset)],
        &offsetMilliseconds, cycleBeat);
    std::array<CvSample, CvLaneCount> cv{};
    for (std::size_t cvIndex = 0; cvIndex < cv.size(); ++cvIndex) {
      const auto lane = CvCursorLane(cvIndex);
      cv[cvIndex] = SampleCv(
          sequence->cv[cvIndex],
          sequence->transforms[static_cast<std::size_t>(lane)],
          sequence->alignment[static_cast<std::size_t>(lane)], structuralCell,
          state.structuralCellCount,
          sequence->alignmentSplits[static_cast<std::size_t>(lane)],
          sequence->alignmentSpans[static_cast<std::size_t>(lane)], atomBeat,
          cycleBeat, atomSpan, sequence->cvInterpolation[cvIndex],
          sequence->cvPower[cvIndex], cycle, seed, 0x9000 + cvIndex * 0x100);
    }
    auto applyControls = [&](RuntimeEvent &event) {
      if (offsetMilliseconds)
        event.timingOffsetMilliseconds += offset;
      else
        event.timingOffsetBeats += offset;
      event.cursors[static_cast<std::size_t>(CursorLane::Offset)] = offsetSpan;
      for (std::size_t cvIndex = 0; cvIndex < cv.size(); ++cvIndex) {
        event.cvValue[cvIndex] = cv[cvIndex].value;
        event.cvTarget[cvIndex] = cv[cvIndex].target;
        event.cvTargetBeat[cvIndex] = cv[cvIndex].targetBeat;
        event.cvInterpolation[cvIndex] = sequence->cvInterpolation[cvIndex];
        event.cvPower[cvIndex] = static_cast<float>(sequence->cvPower[cvIndex]);
        event.cvEnvelope[cvIndex] = sequence->cvEnvelope[cvIndex];
        const auto lane = CvCursorLane(cvIndex);
        event.cursors[static_cast<std::size_t>(lane)] =
            cv[cvIndex].span.valid() ? cv[cvIndex].span
                                     : sequence->cvEnvelope[cvIndex].span;
      }
    };
    std::array<SourceSpan, 5> suppressedLaneSpans{};
    const bool suppressedPitchedEvent =
        kind == ArticulationKind::Rest &&
        (atom.kind == ArticulationKind::Attack ||
         atom.kind == ArticulationKind::Slide) &&
        (atom.hasPitch || sequence->separateRhythm);
    if (suppressedPitchedEvent) {
      ++state.notes;
      const std::array<std::pair<const std::vector<ScalarItem> *, CursorLane>,
                       5>
          lanes{{{&sequence->octave, CursorLane::Octave},
                 {&sequence->velocity, CursorLane::Velocity},
                 {&sequence->gate, CursorLane::Gate},
                 {&sequence->slide, CursorLane::Slide},
                 {&sequence->ratchet, CursorLane::Ratchet}}};
      std::array<std::uint64_t *, 5> cursors{{
          &state.octave,
          &state.velocity,
          &state.gate,
          &state.slide,
          &state.ratchet,
      }};
      for (std::size_t laneIndex = 0; laneIndex < lanes.size(); ++laneIndex) {
        const auto lane = lanes[laneIndex].second;
        Scalar(*lanes[laneIndex].first, *cursors[laneIndex],
               sequence->transforms[static_cast<std::size_t>(lane)], cycle,
               seed, 0.0, suppressedLaneSpans[laneIndex],
               0xa000 + laneIndex * 0x100,
               sequence->alignment[static_cast<std::size_t>(lane)],
               structuralCell, state.structuralCellCount,
               sequence->alignmentSplits[static_cast<std::size_t>(lane)],
               sequence->alignmentSpans[static_cast<std::size_t>(lane)],
               nullptr, cycleBeat);
      }
    }
    if (kind == ArticulationKind::Tie || kind == ArticulationKind::Rest) {
      if (output.count >= output.capacity) {
        output.overflowed = true;
        break;
      }
      auto &event = output.events[output.count++];
      event = {};
      event.kind =
          kind == ArticulationKind::Tie ? EventKind::Tie : EventKind::Rest;
      event.beat = atomBeat;
      event.spanBeats = atomSpan;
      event.cursors[static_cast<std::size_t>(CursorLane::Sequence)] = partSpan;
      event.cursors[static_cast<std::size_t>(
          sequence->separateRhythm ? CursorLane::Rhythm : CursorLane::Notes)] =
          atom.span;
      event.cursors[static_cast<std::size_t>(CursorLane::Duration)] =
          durationSpan;
      if (suppressedPitchedEvent) {
        for (std::size_t laneIndex = 0; laneIndex < suppressedLaneSpans.size();
             ++laneIndex) {
          constexpr std::array<CursorLane, 5> lanes{
              CursorLane::Octave, CursorLane::Velocity, CursorLane::Gate,
              CursorLane::Slide, CursorLane::Ratchet};
          event.cursors[static_cast<std::size_t>(lanes[laneIndex])] =
              suppressedLaneSpans[laneIndex];
        }
      }
      applyControls(event);
      if (kind == ArticulationKind::Rest) {
        state.hasSoundingPitch = false;
        state.soundingVoiceCount = 0;
      }
      continue;
    }

    const auto &noteTransforms =
        sequence->transforms[static_cast<std::size_t>(CursorLane::Notes)];
    const PitchItem *selectedNote = nullptr;
    if (atom.hasPitch) {
      selectedNote = &atom.pitch;
    } else if (timelineStep &&
               timelineStep->kind == PitchTimelineStep::Kind::Pitch) {
      selectedNote = &timelineStep->pitch;
    } else {
      selectedNote = &Cycled(sequence->notes, state.notes, noteTransforms,
                             cycle, seed, 0x4000);
    }
    const auto &note = *selectedNote;
    const auto noteKey =
        sequence->separateRhythm
            ? static_cast<std::uint64_t>(timelineSelection.logicalIndex)
        : atom.hasPitch ? std::uint64_t{0}
                        : state.notes;
    const auto randomChoiceKey = sequence->separateRhythm
                                     ? noteKey ^ MixRandom(cycle)
                                 : atom.hasPitch ? cycle
                                                 : state.notes;
    ++state.notes;
    std::size_t choice = 0;
    if (note.choice == PitchItem::Choice::Alternate)
      choice = static_cast<std::size_t>(cycle % note.values.size());
    else if (note.choice == PitchItem::Choice::Random)
      choice = std::min<std::size_t>(
          note.values.size() - 1,
          static_cast<std::size_t>(
              RandomUnit(seed, randomChoiceKey, note.randomIdentity, 0x4030) *
              note.values.size()));
    PitchValue randomPitch;
    const ChordValue *chord = nullptr;
    if (note.randomDomain != PitchItem::RandomDomain::None) {
      randomPitch = SampleRandomPitch(note, *sequence, seed,
                                      noteKey ^ MixRandom(cycle), 0x4050);
    } else {
      chord = &note.values[choice];
    }
    const auto degreeTranspose =
        ModulateDegrees(noteTransforms, sequence->scale, cycle, seed, 0x4100);
    const auto degreeShift = Transpose(
        noteTransforms, TransformKind::ShiftDegree, cycle, seed, 0x4150);
    const auto semitoneTranspose = SaturatingAdd(
        Transpose(noteTransforms, TransformKind::TransposeSemitone, cycle, seed,
                  0x4200),
        TransposeToKey(noteTransforms, *sequence, cycle, seed, 0x4250));
    const auto octaveTranspose = Transpose(
        noteTransforms, TransformKind::TransposeOctave, cycle, seed, 0x4300);

    SourceSpan velocitySpan;
    SourceSpan octaveSpan;
    SourceSpan gateSpan;
    SourceSpan slideSpan;
    SourceSpan ratchetSpan;
    int sequenceOctave = sequence->scale.tonicOctave;
    if (!sequence->octave.empty()) {
      sequenceOctave = static_cast<int>(std::lround(Scalar(
          sequence->octave, state.octave,
          sequence->transforms[static_cast<std::size_t>(CursorLane::Octave)],
          cycle, seed, sequence->scale.tonicOctave, octaveSpan, 0x4400,
          sequence->alignment[static_cast<std::size_t>(CursorLane::Octave)],
          structuralCell, state.structuralCellCount,
          sequence
              ->alignmentSplits[static_cast<std::size_t>(CursorLane::Octave)],
          sequence
              ->alignmentSpans[static_cast<std::size_t>(CursorLane::Octave)],
          nullptr, cycleBeat)));
    }
    const float baseVelocity = static_cast<float>(Scalar(
        sequence->velocity, state.velocity,
        sequence->transforms[static_cast<std::size_t>(CursorLane::Velocity)],
        cycle, seed, 0.72, velocitySpan, 0x4800,
        sequence->alignment[static_cast<std::size_t>(CursorLane::Velocity)],
        structuralCell, state.structuralCellCount,
        sequence
            ->alignmentSplits[static_cast<std::size_t>(CursorLane::Velocity)],
        sequence
            ->alignmentSpans[static_cast<std::size_t>(CursorLane::Velocity)],
        nullptr, cycleBeat));
    float velocity = baseVelocity;
    float gateDefault = 0.8f;
    if (atom.gateArticulation == GateArticulation::Staccato)
      gateDefault = 0.25f;
    else if (atom.gateArticulation == GateArticulation::Tenuto)
      gateDefault = 0.95f;
    bool gateMilliseconds = false;
    bool gateNoteValue = false;
    float gate = static_cast<float>(Scalar(
        sequence->gate, state.gate,
        sequence->transforms[static_cast<std::size_t>(CursorLane::Gate)], cycle,
        seed, gateDefault, gateSpan, 0x6000,
        sequence->alignment[static_cast<std::size_t>(CursorLane::Gate)],
        structuralCell, state.structuralCellCount,
        sequence->alignmentSplits[static_cast<std::size_t>(CursorLane::Gate)],
        sequence->alignmentSpans[static_cast<std::size_t>(CursorLane::Gate)],
        &gateMilliseconds, cycleBeat, &gateNoteValue));
    bool slideMilliseconds = false;
    float slide = static_cast<float>(Scalar(
        sequence->slide, state.slide,
        sequence->transforms[static_cast<std::size_t>(CursorLane::Slide)],
        cycle, seed, sequence->glideBeats, slideSpan, 0x7000,
        sequence->alignment[static_cast<std::size_t>(CursorLane::Slide)],
        structuralCell, state.structuralCellCount,
        sequence->alignmentSplits[static_cast<std::size_t>(CursorLane::Slide)],
        sequence->alignmentSpans[static_cast<std::size_t>(CursorLane::Slide)],
        &slideMilliseconds, cycleBeat));
    float effectiveAccent = 0.f;
    if (atom.hasVelocity)
      velocity = atom.velocity;
    if (atom.hasAccent) {
      effectiveAccent = atom.accent;
      velocity = std::max(velocity, effectiveAccent);
    }
    if (atom.quiet)
      velocity *= 0.5f;
    if (atom.hasGate) {
      gate = atom.gate;
      gateMilliseconds = atom.gateMilliseconds;
      gateNoteValue = atom.gateNoteValue;
    }
    if (timelineStep && timelineStep->hasSlide) {
      slide = timelineStep->slide;
      slideMilliseconds = timelineStep->slideMilliseconds;
    }
    if (atom.hasSlide) {
      slide = atom.slide;
      slideMilliseconds = atom.slideMilliseconds;
    }
    if (atom.ghost) {
      velocity *= 0.5f;
      effectiveAccent = 0.f;
    }
    const float gateFraction =
        atom.ghost ? std::min(gateMilliseconds ? 0.1f : gate, 0.1f) : gate;
    const double laneRatchet = Scalar(
        sequence->ratchet, state.ratchet,
        sequence->transforms[static_cast<std::size_t>(CursorLane::Ratchet)],
        cycle, seed, 1.0, ratchetSpan, 0x8000,
        sequence->alignment[static_cast<std::size_t>(CursorLane::Ratchet)],
        structuralCell, state.structuralCellCount,
        sequence
            ->alignmentSplits[static_cast<std::size_t>(CursorLane::Ratchet)],
        sequence->alignmentSpans[static_cast<std::size_t>(CursorLane::Ratchet)],
        nullptr, cycleBeat);
    const std::size_t ratchets =
        atom.ratchets > 1
            ? atom.ratchets
            : std::max<std::size_t>(1, static_cast<std::size_t>(laneRatchet));

    std::array<float, MaximumPolyphony> pitches{};
    std::array<SourceSpan, MaximumPolyphony> pitchSpans{};
    std::size_t voiceCount = 0;
    if (!chord) {
      pitches[voiceCount] =
          PitchVolts(randomPitch, *sequence, sequenceOctave, degreeShift,
                     degreeTranspose, semitoneTranspose, octaveTranspose);
      pitchSpans[voiceCount] = randomPitch.span;
      ++voiceCount;
    } else if ((chord->meaning == ChordValue::Meaning::RomanSymbol ||
                chord->meaning == ChordValue::Meaning::JazzSymbol) &&
               !chord->tones.empty()) {
      const auto &rootPitch = chord->meaning == ChordValue::Meaning::RomanSymbol
                                  ? chord->romanRoot
                                  : chord->harmonicRoot;
      const float root =
          PitchVolts(rootPitch, *sequence, sequenceOctave, degreeShift,
                     degreeTranspose, semitoneTranspose, octaveTranspose);
      const int rootSemitone = BoundedSemitone(root);
      const auto realized = RealizeChordVoicing(
          *chord, sequence->voicing, rootSemitone,
          state.previousVoicingSemitones, state.previousVoicingCount);
      for (std::size_t voice = 0;
           voice < realized.count && voice < MaximumPolyphony; ++voice) {
        pitches[voiceCount] =
            static_cast<float>(realized.semitones[voice]) / 12.f;
        pitchSpans[voiceCount] = rootPitch.span;
        ++voiceCount;
      }
    } else {
      for (const auto &voice : chord->voices) {
        if (voiceCount >= MaximumPolyphony)
          break;
        pitches[voiceCount] =
            PitchVolts(voice, *sequence, sequenceOctave, degreeShift,
                       degreeTranspose, semitoneTranspose, octaveTranspose);
        pitchSpans[voiceCount] = voice.span;
        ++voiceCount;
      }
    }
    const std::size_t harmonicVoiceCount = voiceCount;
    std::array<int, MaximumPolyphony> harmonicSemitones{};
    for (std::size_t voice = 0; voice < harmonicVoiceCount; ++voice)
      harmonicSemitones[voice] = BoundedSemitone(pitches[voice]);
    if (chord && chord->hasBass && voiceCount > 0) {
      float bass =
          PitchVolts(chord->bass, *sequence, sequenceOctave, degreeShift,
                     degreeTranspose, semitoneTranspose, octaveTranspose);
      std::size_t matchingVoice = voiceCount;
      for (std::size_t voice = 0; voice < voiceCount; ++voice) {
        const double octaves = static_cast<double>(pitches[voice] - bass);
        if (std::abs(octaves - std::round(octaves)) < 1e-5) {
          matchingVoice = voice;
          break;
        }
      }
      if (matchingVoice < voiceCount) {
        for (std::size_t voice = matchingVoice; voice + 1 < voiceCount;
             ++voice) {
          pitches[voice] = pitches[voice + 1];
          pitchSpans[voice] = pitchSpans[voice + 1];
        }
        --voiceCount;
      }
      const bool explicitBassRegister =
          chord->bass.hasOctave || chord->bass.octaveOffset != 0;
      if (!explicitBassRegister && voiceCount > 0) {
        const float lowest =
            *std::min_element(pitches.begin(), pitches.begin() + voiceCount);
        while (bass >= lowest - 1e-5f)
          bass -= 1.f;
      }
      if (voiceCount < MaximumPolyphony) {
        for (std::size_t voice = voiceCount; voice > 0; --voice) {
          pitches[voice] = pitches[voice - 1];
          pitchSpans[voice] = pitchSpans[voice - 1];
        }
        pitches[0] = bass;
        pitchSpans[0] = chord->bass.span;
        ++voiceCount;
      }
    }

    if (kind == ArticulationKind::Slide &&
        state.soundingVoiceCount != voiceCount)
      kind = ArticulationKind::Attack;
    if (legatoFromAtom && legatoTargetVoiceCount != 0 &&
        legatoTargetVoiceCount != voiceCount)
      legatoFromAtom = false;

    const auto firstEmitted = output.count;
    for (std::size_t ratchet = 0; ratchet < ratchets; ++ratchet) {
      for (std::size_t voice = 0; voice < voiceCount; ++voice) {
        if (output.count >= output.capacity) {
          output.overflowed = true;
          break;
        }
        auto &event = output.events[output.count++];
        event = {};
        event.kind = kind == ArticulationKind::Slide ? EventKind::Slide
                                                     : EventKind::Attack;
        event.beat = atomBeat + atomSpan * static_cast<double>(ratchet) /
                                    static_cast<double>(ratchets);
        event.spanBeats = atomSpan / static_cast<double>(ratchets);
        event.pitchVolts = pitches[voice];
        event.velocity = std::clamp(velocity, 0.f, 1.f);
        event.accent = std::clamp(effectiveAccent, 0.f, 1.f);
        event.gateFraction = std::clamp(gateFraction, 0.f, 1.f);
        event.gateBeats = gateNoteValue ? gate : -1.f;
        event.gateMilliseconds = gateMilliseconds ? gate : -1.f;
        event.gateCapMilliseconds = atom.ghost ? 20.f : -1.f;
        event.slideBeats = slide;
        event.slideMilliseconds = slideMilliseconds ? slide : -1.f;
        event.voice = static_cast<std::uint8_t>(voice);
        event.voiceCount = static_cast<std::uint8_t>(voiceCount);
        event.cursors[static_cast<std::size_t>(CursorLane::Sequence)] =
            partSpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Notes)] =
            !sequence->separateRhythm && atom.span.valid() ? atom.span
            : timelineStep && timelineStep->span.valid()   ? timelineStep->span
            : chord && chord->span.valid()                 ? chord->span
            : pitchSpans[voice].valid()                    ? pitchSpans[voice]
                                                           : note.span;
        if (sequence->separateRhythm)
          event.cursors[static_cast<std::size_t>(CursorLane::Rhythm)] =
              atom.span;
        event.cursors[static_cast<std::size_t>(CursorLane::Octave)] =
            octaveSpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Velocity)] =
            velocitySpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Duration)] =
            durationSpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Gate)] = gateSpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Slide)] = slideSpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Ratchet)] =
            ratchetSpan;
        applyControls(event);
      }
      if (output.overflowed)
        break;
    }
    if (legatoFromAtom && output.count > firstEmitted) {
      double latestOnset = -std::numeric_limits<double>::infinity();
      for (std::size_t index = firstEmitted; index < output.count; ++index)
        latestOnset = std::max(latestOnset, output.events[index].beat);
      for (std::size_t index = firstEmitted; index < output.count; ++index) {
        auto &event = output.events[index];
        if (std::abs(event.beat - latestOnset) < 1e-7)
          event.legatoToNext = true;
      }
    }
    if (!atom.ghost && chord &&
        (chord->meaning == ChordValue::Meaning::JazzSymbol ||
         chord->meaning == ChordValue::Meaning::RomanSymbol ||
         chord->meaning == ChordValue::Meaning::ExplicitVoicing)) {
      state.previousVoicingCount = harmonicVoiceCount;
      std::copy_n(harmonicSemitones.begin(), harmonicVoiceCount,
                  state.previousVoicingSemitones.begin());
    }
    state.hasSoundingPitch = !atom.ghost;
    state.soundingVoiceCount = atom.ghost ? 0 : voiceCount;
  }
  state.structuralCell += step ? step->cellCount : 1;
  const bool completesRhythmPass =
      sequence->separateRhythm && !sequence->articulation.empty() &&
      state.articulation % sequence->articulation.size() == 0;
  if (completesRhythmPass && !completesNotesPass) {
    state.structuralCell = 0;
    state.structuralCellCount = 0;
  }
  const auto swing = Swing(timingTransforms, cycle, seed);
  if (swing.ratio != 0.5) {
    const double origin = partStartBeat_;
    for (std::size_t groupBegin = 0; groupBegin < output.count;) {
      std::size_t groupEnd = groupBegin + 1;
      while (groupEnd < output.count &&
             std::abs(output.events[groupEnd].beat -
                      output.events[groupBegin].beat) < 1e-7 &&
             std::abs(output.events[groupEnd].spanBeats -
                      output.events[groupBegin].spanBeats) < 1e-7)
        ++groupEnd;
      const double originalBeat = output.events[groupBegin].beat;
      const double originalSpan = output.events[groupBegin].spanBeats;
      const double subdivision =
          swing.subdivisionBeats > 0.0 ? swing.subdivisionBeats : originalSpan;
      const double gridPosition = (originalBeat - origin) / subdivision;
      const bool representable = std::isfinite(gridPosition) &&
                                 static_cast<long double>(gridPosition) >=
                                     std::numeric_limits<std::int64_t>::min() &&
                                 static_cast<long double>(gridPosition) <=
                                     std::numeric_limits<std::int64_t>::max();
      const auto gridIndex =
          representable ? static_cast<std::int64_t>(std::llround(gridPosition))
                        : std::int64_t{0};
      const bool aligned =
          representable && std::abs(gridPosition - gridIndex) < 1e-7;
      if (aligned) {
        const bool second = (gridIndex & 1LL) != 0;
        const double pairStart =
            origin +
            static_cast<double>(gridIndex - (second ? 1 : 0)) * subdivision;
        const double boundary = pairStart + 2.0 * subdivision * swing.ratio;
        for (std::size_t index = groupBegin; index < groupEnd; ++index) {
          auto &event = output.events[index];
          if (second)
            event.beat = boundary;
          if (std::abs(originalSpan - subdivision) < 1e-7)
            event.spanBeats = second ? 2.0 * subdivision * (1.0 - swing.ratio)
                                     : 2.0 * subdivision * swing.ratio;
        }
      }
      groupBegin = groupEnd;
    }
  }
  std::size_t onsetIndex = 0;
  double previousOnset = -std::numeric_limits<double>::infinity();
  for (std::size_t index = 0; index < output.count; ++index) {
    if (index > 0 &&
        std::abs(output.events[index].beat - previousOnset) >= 1e-7)
      ++onsetIndex;
    previousOnset = output.events[index].beat;
    const double originalBeat = output.events[index].beat;
    output.events[index].beat += output.events[index].timingOffsetBeats;
    ApplyTimingOffsets(output.events[index], timingTransforms, cycle, seed,
                       onsetIndex);
    const double beatShift = output.events[index].beat - originalBeat;
    for (auto &targetBeat : output.events[index].cvTargetBeat)
      targetBeat += beatShift;
  }
  if (sequence->separateRhythm)
    state.pitchTimelineBeat += usedUnscaledDuration;
  if (completesNotesPass) {
    const auto completedCycles = state.completedCycles + 1;
    if (continuesSameSequence) {
      state.completedCycles = completedCycles;
      state.freeLaneBeat +=
          std::max(0.0, beat + output.durationBeats - partStartBeat_);
      state.structuralCell = 0;
      state.structuralCellCount = 0;
      if (sequence->separateRhythm) {
        state.articulation = 0;
        state.pitchTimelineBeat = 0.0;
        state.activePitchTimelineStep = std::numeric_limits<std::size_t>::max();
        state.pitchEntryConsumed = false;
      }
    } else {
      const auto lastBaseDuration = state.lastBaseDuration;
      state = {};
      state.completedCycles = completedCycles;
      state.lastBaseDuration = lastBaseDuration;
    }
    partStartBeat_ = beat + output.durationBeats;
    if (arrangementPart.cycles >= 0) {
      ++arrangementCycle_;
      if (arrangementCycle_ >= arrangementPart.cycles) {
        arrangementCycle_ = 0;
        arrangementPart_ = (arrangementPart_ + 1) % semantic.arrangement.size();
      }
    }
  }
  return output;
}

} // namespace tfseq
