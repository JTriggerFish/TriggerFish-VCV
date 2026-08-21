#include "tfseq.hpp"

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

template <typename T>
const T &Cycled(const std::vector<T> &items, std::uint64_t cursor,
                const std::vector<Transform> &transforms, std::uint64_t cycle,
                std::uint64_t seed, std::uint64_t salt) noexcept {
  const auto index =
      TransformIndex(static_cast<std::size_t>(cursor), items.size(), transforms,
                     cycle, seed, salt);
  return items[index];
}

double Scalar(const std::vector<ScalarItem> &items, std::uint64_t &cursor,
              const std::vector<Transform> &transforms, std::uint64_t cycle,
              std::uint64_t seed, double fallback, SourceSpan &span,
              std::uint64_t salt) noexcept {
  if (items.empty())
    return fallback;
  const auto &item = Cycled(items, cursor++, transforms, cycle, seed, salt);
  span = item.span;
  return item.isDefault ? fallback : item.value;
}

ArticulationKind EffectiveArticulation(const ArticulationAtom &atom,
                                       std::uint64_t articulationCursor,
                                       std::size_t atomIndex,
                                       std::uint64_t seed) noexcept {
  std::uint64_t probabilityHash =
      seed ^ (articulationCursor * 0x9e3779b9ULL) ^ (atomIndex * 0x85ebca6bULL);
  probabilityHash ^= probabilityHash >> 33;
  const float random = static_cast<float>(probabilityHash & 0xffffffU) /
                       static_cast<float>(0x1000000U);
  return random < atom.probability ? atom.kind : ArticulationKind::Rest;
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

void Runtime::replaceProgram(const CompiledProgram *program,
                             double beat) noexcept {
  const auto *previousStates = states_;
  const auto previousStateCount = stateCount_;
  const auto *previousProgram = program_;
  auto *nextStates = program ? program->stateWorkspace.data() : nullptr;
  const auto nextStateCount = program ? program->stateWorkspace.size() : 0;
  if (nextStates)
    std::fill_n(nextStates, nextStateCount, SequencePlaybackState{});
  if (program && previousProgram) {
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
        nextStates[nextIndex] = previousStates[previousIndex];
        ++nextPosition;
        ++previousPosition;
      }
    }
  }

  program_ = program;
  states_ = nextStates;
  stateCount_ = nextStateCount;
  arrangementPart_ = 0;
  arrangementCycle_ = 0;
  partStartBeat_ = 0.0;
  if (!program_ || program_->semantic().arrangement.empty())
    return;
  const auto &semantic = program_->semantic();

  beat = std::max(0.0, beat);
  double arrangementDuration = 0.0;
  bool finite = true;
  for (const auto &part : semantic.arrangement) {
    if (part.cycles < 0) {
      finite = false;
      break;
    }
    arrangementDuration +=
        semantic.sequences[part.sequence].cycleBeats * part.cycles;
  }
  if (finite && arrangementDuration > 0.0) {
    const double passes = std::floor(beat / arrangementDuration);
    partStartBeat_ = passes * arrangementDuration;
    arrangementCycle_ = static_cast<int>(
        std::min(passes, static_cast<double>(std::numeric_limits<int>::max())));
  }

  for (std::size_t guard = 0; guard < semantic.arrangement.size(); ++guard) {
    const auto &part = semantic.arrangement[arrangementPart_];
    if (part.cycles < 0)
      return;
    const double duration =
        semantic.sequences[part.sequence].cycleBeats * part.cycles;
    if (beat < partStartBeat_ + duration)
      return;
    partStartBeat_ += duration;
    arrangementPart_ = (arrangementPart_ + 1) % semantic.arrangement.size();
  }
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

const Sequence *Runtime::currentSequence(double beat, SourceSpan &partSpan,
                                         std::uint64_t &cycle) noexcept {
  if (!program_ || program_->semantic().arrangement.empty())
    return nullptr;
  const auto &semantic = program_->semantic();

  for (std::size_t guard = 0; guard <= semantic.arrangement.size(); ++guard) {
    const auto &part = semantic.arrangement[arrangementPart_];
    const auto &sequence = semantic.sequences[part.sequence];
    const double elapsed = beat - partStartBeat_;
    if (part.cycles < 0 || elapsed < sequence.cycleBeats * part.cycles) {
      cycle = static_cast<std::uint64_t>(
          std::max(0.0, std::floor(elapsed / sequence.cycleBeats)));
      partSpan = part.span.valid() ? part.span : sequence.nameSpan;
      return &sequence;
    }
    partStartBeat_ += sequence.cycleBeats * part.cycles;
    arrangementPart_ = (arrangementPart_ + 1) % semantic.arrangement.size();
    ++arrangementCycle_;
  }
  return nullptr;
}

StepEvents Runtime::next(double beat) noexcept {
  StepEvents output;
  if (program_) {
    output.events = program_->stepWorkspace.data();
    output.capacity = program_->stepWorkspace.size();
  }
  SourceSpan partSpan;
  std::uint64_t cycle = 0;
  const Sequence *sequence = currentSequence(beat, partSpan, cycle);
  if (!sequence)
    return output;
  const auto &semantic = program_->semantic();
  const auto sequenceIndex =
      static_cast<std::size_t>(sequence - semantic.sequences.data());
  if (sequenceIndex >= stateCount_)
    return output;
  auto &state = states_[sequenceIndex];
  const auto seed = semantic.seed;
  const auto &timingTransforms =
      sequence->transforms[static_cast<std::size_t>(CursorLane::Sequence)];

  SourceSpan durationSpan;
  const double duration =
      static_cast<double>(Scalar(
          sequence->duration, state.duration,
          sequence->transforms[static_cast<std::size_t>(CursorLane::Duration)],
          cycle, seed, 1.f, durationSpan, 0x2000)) *
      TimeScale(timingTransforms, cycle, seed);
  output.durationBeats = duration;

  const ArticulationAtom implicit{ArticulationKind::Attack, 1, 1.f,
                                  sequence->nameSpan};
  const ArticulationStep *step = nullptr;
  if (!sequence->articulation.empty()) {
    step = &Cycled(
        sequence->articulation, state.articulation++,
        sequence
            ->transforms[static_cast<std::size_t>(CursorLane::Articulation)],
        cycle, seed, 0x3000);
  }

  const std::size_t atomCount = step ? step->atoms.size() : 1;
  bool legatoToNext = false;
  if (step && !sequence->articulation.empty()) {
    const auto &part = semantic.arrangement[arrangementPart_];
    const double nextElapsed = beat + duration - partStartBeat_;
    const bool samePart =
        part.cycles < 0 ||
        nextElapsed < sequence->cycleBeats * part.cycles - 1e-9;
    if (samePart) {
      const auto nextCycle = static_cast<std::uint64_t>(
          std::max(0.0, std::floor(nextElapsed / sequence->cycleBeats)));
      const auto &nextStep = Cycled(
          sequence->articulation, state.articulation,
          sequence
              ->transforms[static_cast<std::size_t>(CursorLane::Articulation)],
          nextCycle, seed, 0x3000);
      if (!nextStep.atoms.empty()) {
        const auto nextKind = EffectiveArticulation(
            nextStep.atoms.front(), state.articulation + 1, 0, seed);
        legatoToNext = nextKind == ArticulationKind::Tie ||
                       nextKind == ArticulationKind::Slide;
      }
    }
  }

  for (std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex) {
    const auto &atom = step ? step->atoms[atomIndex] : implicit;
    const ArticulationKind kind =
        EffectiveArticulation(atom, state.articulation, atomIndex, seed);

    const double atomBeat = beat + duration * static_cast<double>(atomIndex) /
                                       static_cast<double>(atomCount);
    const double atomSpan = duration / static_cast<double>(atomCount);
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
      event.cursors[static_cast<std::size_t>(CursorLane::Articulation)] =
          atom.span;
      event.cursors[static_cast<std::size_t>(CursorLane::Duration)] =
          durationSpan;
      continue;
    }

    const auto &noteTransforms =
        sequence->transforms[static_cast<std::size_t>(CursorLane::Notes)];
    const auto &note = Cycled(sequence->notes, state.notes++, noteTransforms,
                              cycle, seed, 0x4000);
    std::size_t choice = 0;
    if (note.choice == PitchItem::Choice::Alternate)
      choice = static_cast<std::size_t>(cycle % note.values.size());
    else if (note.choice == PitchItem::Choice::Random) {
      std::uint64_t hash = seed ^ state.notes ^ (cycle * 0x9e3779b97f4a7c15ULL);
      hash ^= hash >> 29;
      choice = static_cast<std::size_t>(hash % note.values.size());
    }
    const auto &chord = note.values[choice];
    const auto degreeTranspose =
        ModulateDegrees(noteTransforms, sequence->scale, cycle, seed, 0x4100);
    const auto degreeShift = Transpose(
        noteTransforms, TransformKind::ShiftDegree, cycle, seed, 0x4150);
    const auto semitoneTranspose = Transpose(
        noteTransforms, TransformKind::TransposeSemitone, cycle, seed, 0x4200);
    const auto octaveTranspose = Transpose(
        noteTransforms, TransformKind::TransposeOctave, cycle, seed, 0x4300);

    SourceSpan velocitySpan;
    SourceSpan octaveSpan;
    SourceSpan accentSpan;
    SourceSpan gateSpan;
    SourceSpan slideSpan;
    SourceSpan ratchetSpan;
    int sequenceOctave = sequence->scale.tonicOctave;
    if (!sequence->octave.empty()) {
      sequenceOctave = static_cast<int>(std::lround(Scalar(
          sequence->octave, state.octave,
          sequence->transforms[static_cast<std::size_t>(CursorLane::Octave)],
          cycle, seed, 0.f, octaveSpan, 0x4400)));
    }
    const float baseVelocity = static_cast<float>(Scalar(
        sequence->velocity, state.velocity,
        sequence->transforms[static_cast<std::size_t>(CursorLane::Velocity)],
        cycle, seed, 0.72, velocitySpan, 0x4800));
    const float accent = static_cast<float>(Scalar(
        sequence->accent, state.accent,
        sequence->transforms[static_cast<std::size_t>(CursorLane::Accent)],
        cycle, seed, 0.0, accentSpan, 0x5000));
    const float velocity = std::max(baseVelocity, accent);
    const float gate = static_cast<float>(
        Scalar(sequence->gate, state.gate,
               sequence->transforms[static_cast<std::size_t>(CursorLane::Gate)],
               cycle, seed, 0.8, gateSpan, 0x6000));
    const float slide = static_cast<float>(Scalar(
        sequence->slide, state.slide,
        sequence->transforms[static_cast<std::size_t>(CursorLane::Slide)],
        cycle, seed, 0.0, slideSpan, 0x7000));
    const double laneRatchet = Scalar(
        sequence->ratchet, state.ratchet,
        sequence->transforms[static_cast<std::size_t>(CursorLane::Ratchet)],
        cycle, seed, 1.0, ratchetSpan, 0x8000);
    const std::size_t ratchets =
        atom.ratchets > 1
            ? atom.ratchets
            : std::max<std::size_t>(1, static_cast<std::size_t>(laneRatchet));

    std::array<float, MaximumPolyphony> pitches{};
    std::array<SourceSpan, MaximumPolyphony> pitchSpans{};
    std::size_t voiceCount = 0;
    for (const auto &voice : chord.voices) {
      if (voiceCount >= MaximumPolyphony)
        break;
      pitches[voiceCount] =
          PitchVolts(voice, *sequence, sequenceOctave, degreeShift,
                     degreeTranspose, semitoneTranspose, octaveTranspose);
      pitchSpans[voiceCount] = voice.span;
      ++voiceCount;
    }
    if (chord.hasBass && voiceCount > 0) {
      float bass =
          PitchVolts(chord.bass, *sequence, sequenceOctave, degreeShift,
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
          chord.bass.hasOctave || chord.bass.octaveOffset != 0;
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
        pitchSpans[0] = chord.bass.span;
        ++voiceCount;
      }
    }

    for (std::size_t ratchet = 0; ratchet < ratchets; ++ratchet) {
      for (std::size_t voice = 0; voice < voiceCount; ++voice) {
        if (output.count >= output.capacity) {
          output.overflowed = true;
          break;
        }
        auto &event = output.events[output.count++];
        event = {};
        event.kind = (kind == ArticulationKind::Slide || slide > 0.f)
                         ? EventKind::Slide
                         : EventKind::Attack;
        event.beat = atomBeat + atomSpan * static_cast<double>(ratchet) /
                                    static_cast<double>(ratchets);
        event.spanBeats = atomSpan / static_cast<double>(ratchets);
        event.pitchVolts = pitches[voice];
        event.velocity = std::clamp(velocity, 0.f, 1.f);
        event.accent = std::clamp(accent, 0.f, 1.f);
        event.gateFraction = std::clamp(gate, 0.f, 1.f);
        event.slideBeats = slide > 0.f ? slide : sequence->glideBeats;
        event.voice = static_cast<std::uint8_t>(voice);
        event.voiceCount = static_cast<std::uint8_t>(voiceCount);
        event.cursors[static_cast<std::size_t>(CursorLane::Sequence)] =
            partSpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Notes)] =
            chord.span.valid()          ? chord.span
            : pitchSpans[voice].valid() ? pitchSpans[voice]
                                        : note.span;
        event.cursors[static_cast<std::size_t>(CursorLane::Octave)] =
            octaveSpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Articulation)] =
            atom.span;
        event.cursors[static_cast<std::size_t>(CursorLane::Velocity)] =
            velocitySpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Accent)] =
            accentSpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Duration)] =
            durationSpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Gate)] = gateSpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Slide)] = slideSpan;
        event.cursors[static_cast<std::size_t>(CursorLane::Ratchet)] =
            ratchetSpan;
      }
      if (output.overflowed)
        break;
    }
  }
  if (legatoToNext && output.count > 0) {
    double latestOnset = -std::numeric_limits<double>::infinity();
    for (std::size_t index = 0; index < output.count; ++index) {
      const auto &event = output.events[index];
      if (event.kind == EventKind::Attack || event.kind == EventKind::Slide)
        latestOnset = std::max(latestOnset, event.beat);
    }
    for (std::size_t index = 0; index < output.count; ++index) {
      auto &event = output.events[index];
      if ((event.kind == EventKind::Attack || event.kind == EventKind::Slide) &&
          std::abs(event.beat - latestOnset) < 1e-7)
        event.legatoToNext = true;
    }
  }
  const auto swing = Swing(timingTransforms, cycle, seed);
  if (swing.ratio != 0.5) {
    const double origin =
        partStartBeat_ + static_cast<double>(cycle) * sequence->cycleBeats;
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
    ApplyTimingOffsets(output.events[index], timingTransforms, cycle, seed,
                       onsetIndex);
  }
  return output;
}

} // namespace tfseq
