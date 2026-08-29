#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace tfseq {

inline constexpr double NormalDeviationLimit = 4.0;

namespace syntax {
struct Document;
}

constexpr std::size_t MaximumPolyphony = 16;
constexpr std::size_t CvLaneCount = 3;
// Millisecond-delayed events are prepared conservatively for clocks up to
// 1 kHz. Faster clocking remains playable but can report workspace exhaustion
// rather than allocating from the audio thread.
constexpr double PreparedMaximumClockRateHz = 1000.0;

struct SourceSpan {
  int begin = -1;
  int end = -1;
  int line = 1;
  int column = 1;

  bool valid() const noexcept { return begin >= 0 && end > begin; }
};

struct Diagnostic {
  std::string message;
  int line = 1;
  int column = 1;
};

enum class CursorLane : std::size_t {
  Sequence,
  Notes,
  Rhythm,
  Octave,
  Velocity,
  Duration,
  Gate,
  Slide,
  Ratchet,
  Offset,
  Cv1,
  Cv2,
  Cv3,
  Count
};

constexpr CursorLane CvCursorLane(const std::size_t index) noexcept {
  return static_cast<CursorLane>(static_cast<std::size_t>(CursorLane::Cv1) +
                                 index);
}

constexpr bool IsCvCursorLane(const CursorLane lane) noexcept {
  return lane >= CursorLane::Cv1 && lane <= CursorLane::Cv3;
}

constexpr std::size_t CvCursorIndex(const CursorLane lane) noexcept {
  return static_cast<std::size_t>(lane) -
         static_cast<std::size_t>(CursorLane::Cv1);
}

enum class EventKind { Attack, Slide, Tie, Rest };

enum class ArticulationKind { Attack, Slide, Tie, Rest };

enum class GateArticulation { Normal, Staccato, Tenuto };

struct PitchValue {
  bool absolute = false;
  int degree = 1;
  int accidental = 0;
  int pitchClass = 0;
  bool hasOctave = false;
  int octave = 0;
  int octaveOffset = 0;
  SourceSpan span;
};

// A chord factor retains its harmonic spelling as well as its sounding
// interval.  Voicing recipes operate on factors (3, b7, 9...) rather than on
// anonymous pitch classes so that guide tones and optional roots/fifths can
// be treated deliberately.
struct ChordTone {
  int degree = 1;
  int accidental = 0;
  int semitones = 0;
};

enum class VoicingStyle { Basic, Rootless3Notes, Rootless4Notes };

struct ChordValue {
  enum class Meaning { SinglePitch, ExplicitVoicing, JazzSymbol, RomanSymbol };
  Meaning meaning = Meaning::SinglePitch;
  // Preserve harmonic intent beside the current realization so the voicing
  // engine and future interpreters do not have to reparse display text.
  std::string jazzSymbol;
  int rootPitchClass = 0;
  PitchValue harmonicRoot;
  PitchValue romanRoot;
  std::vector<int> intervals;
  std::vector<ChordTone> tones;
  std::vector<PitchValue> voices;
  bool altered = false;
  bool hasFactorOverride = false;
  bool hasBass = false;
  PitchValue bass;
  SourceSpan span;
};

struct PitchItem {
  enum class Choice { Single, Alternate, Random };
  enum class RandomDomain { None, ScaleDegree, ChromaticSemitone };
  enum class RandomDistribution { Uniform, Normal };
  Choice choice = Choice::Single;
  RandomDomain randomDomain = RandomDomain::None;
  RandomDistribution randomDistribution = RandomDistribution::Uniform;
  // A bare `$` derives its inclusive scale-degree range from the active
  // scale. Explicit uniform bounds and normal mean/deviation use these fields.
  bool randomDefaultRange = false;
  double randomFirst = 0.0;
  double randomSecond = 0.0;
  std::uint64_t randomIdentity = 0;
  std::vector<ChordValue> values;
  SourceSpan span;
};

// When pitch and rhythm are separated, this is the slower held-value
// timeline. A value may contain one pitch or a polyphonic voicing. Rhythm hits
// sample the active entry; they never advance this timeline themselves.
struct PitchTimelineStep {
  enum class Kind { Pitch, Rest };
  Kind kind = Kind::Pitch;
  PitchItem pitch;
  double durationMultiplier = 1.0;
  bool slideOnEntry = false;
  bool hasSlide = false;
  float slide = 0.f;
  bool slideMilliseconds = false;
  SourceSpan span;
};

struct ArticulationAtom {
  struct ProbabilityGate {
    float probability = 1.f;
    std::uint64_t group = 0;
  };
  ArticulationKind kind = ArticulationKind::Attack;
  PitchItem pitch;
  bool hasPitch = false;
  std::size_t ratchets = 1;
  float probability = 1.f;
  std::uint64_t probabilityGroup = 0;
  std::vector<ProbabilityGate> enclosingProbabilityGates;
  double offsetFraction = 0.0;
  double spanFraction = 1.0;
  std::size_t cellOffset = 0;
  bool ghost = false;
  bool quiet = false;
  GateArticulation gateArticulation = GateArticulation::Normal;
  bool hasVelocity = false;
  float velocity = 0.f;
  bool hasAccent = false;
  float accent = 0.f;
  bool hasGate = false;
  float gate = 0.f;
  bool gateMilliseconds = false;
  bool gateNoteValue = false;
  bool hasSlide = false;
  float slide = 0.f;
  bool slideMilliseconds = false;
  SourceSpan span;
};

struct ArticulationStep {
  std::vector<ArticulationAtom> atoms;
  double durationMultiplier = 1.0;
  std::size_t cellCount = 1;
  float presenceProbability = 1.f;
  std::uint64_t presenceIdentity = 0;
  SourceSpan span;
};

struct ScalarItem {
  enum class RandomDistribution { None, Uniform, Normal };
  double value = 0.0;
  double randomFirst = 0.0;
  double randomSecond = 0.0;
  double randomMinimum = -std::numeric_limits<double>::infinity();
  double randomMaximum = std::numeric_limits<double>::infinity();
  std::uint64_t randomIdentity = 0;
  RandomDistribution randomDistribution = RandomDistribution::None;
  bool isDefault = false;
  bool isMilliseconds = false;
  bool isNoteValue = false;
  bool randomInteger = false;
  SourceSpan span;
};

enum class CvInterpolation { Step, Linear, Smooth, Power };

enum class CvEnvelopeMode { Ad, Ar, Adsr };
enum class CvEnvelopeComposition { Replace, Add };
enum class CvEnvelopeTimeUnit { Seconds, Beats };

struct CvEnvelopeTime {
  double value = 0.0;
  CvEnvelopeTimeUnit unit = CvEnvelopeTimeUnit::Seconds;
};

struct CvEnvelopeSpec {
  bool enabled = false;
  CvEnvelopeMode mode = CvEnvelopeMode::Ad;
  CvEnvelopeComposition composition = CvEnvelopeComposition::Replace;
  CvEnvelopeTime attack{0.005, CvEnvelopeTimeUnit::Seconds};
  CvEnvelopeTime decay{0.300, CvEnvelopeTimeUnit::Seconds};
  CvEnvelopeTime release{0.300, CvEnvelopeTimeUnit::Seconds};
  float sustain = 0.5f;
  float depth = 5.f;
  float curve = 0.f;
  bool followVelocity = false;
  float accentMultiplier = 1.f;
  SourceSpan span;
};

enum class LaneAlignment { Free, Left, Right, Edges };

enum class TransformKind {
  Reverse,
  Rotate,
  ModulateDegree,
  ShiftDegree,
  TransposeSemitone,
  TransposeKey,
  TransposeOctave,
  Fast,
  Slow,
  Swing,
  Early,
  Late,
  Rate
};

enum class TransformCondition { Always, Every, Sometimes };
enum class TransformDomain { General, Pitch, Timing, Phase };
enum class TimeUnit { Beats, Milliseconds };

struct Transform {
  TransformKind kind = TransformKind::Reverse;
  TransformDomain domain = TransformDomain::General;
  TransformCondition condition = TransformCondition::Always;
  int integer = 0;
  double number = 0.0;
  // Zero selects the event's own subdivision. A positive value supplies an
  // explicit grid in incoming-clock beats.
  double swingSubdivisionBeats = 0.0;
  bool randomAmount = false;
  TimeUnit timeUnit = TimeUnit::Beats;
  float probability = 1.f;
  int period = 1;
  SourceSpan span;
  // modulate_degree captures the scale active at its position in a pipeline,
  // preserving left-to-right composition when a later scale transform occurs.
  std::vector<int> modulationIntervals;
};

struct Scale {
  int tonicSemitone = 0;
  int tonicOctave = 4;
  std::vector<int> intervals{0, 2, 4, 5, 7, 9, 11};
};

struct Sequence {
  std::uint64_t stableId = 0;
  std::string name;
  SourceSpan nameSpan;
  double subdivisionBeats = 1.0;
  float glideBeats = 0.25f;
  Scale scale;
  // A chord key is an explicit transposition anchor and piece of harmonic
  // metadata. It never changes the meaning of an explicitly named chord root.
  bool hasKey = false;
  int keyPitchClass = 0;
  VoicingStyle voicing = VoicingStyle::Basic;
  std::vector<PitchItem> notes;
  std::vector<PitchTimelineStep> pitchTimeline;
  double pitchTimelineBeats = 0.0;
  double rhythmSubdivisionBeats = 1.0;
  // In compact form each pitched event owns its onset. With a rhythm lane,
  // pitchTimeline supplies held values and articulation supplies independent
  // attacks. Single pitches and chords use this same representation.
  bool separateRhythm = false;
  std::vector<ScalarItem> octave;
  std::vector<ArticulationStep> articulation;
  std::vector<ScalarItem> velocity;
  std::vector<ScalarItem> duration;
  std::vector<ScalarItem> gate;
  std::vector<ScalarItem> slide;
  std::vector<ScalarItem> ratchet;
  std::vector<ScalarItem> offset;
  std::array<std::vector<ScalarItem>, CvLaneCount> cv;
  std::array<CvInterpolation, CvLaneCount> cvInterpolation{
      CvInterpolation::Step, CvInterpolation::Step, CvInterpolation::Step};
  std::array<double, CvLaneCount> cvPower{1.0, 1.0, 1.0};
  std::array<CvEnvelopeSpec, CvLaneCount> cvEnvelope{};
  std::array<std::vector<Transform>,
             static_cast<std::size_t>(CursorLane::Count)>
      transforms;
  std::array<LaneAlignment, static_cast<std::size_t>(CursorLane::Count)>
      alignment{};
  std::array<std::size_t, static_cast<std::size_t>(CursorLane::Count)>
      alignmentSplits{};
  std::array<SourceSpan, static_cast<std::size_t>(CursorLane::Count)>
      alignmentSpans{};
};

struct ArrangementPart {
  std::size_t sequence = 0;
  int cycles = 1;
  SourceSpan span;
};

struct RuntimeEvent;

struct SequencePlaybackState {
  std::uint64_t notes = 0;
  std::uint64_t octave = 0;
  std::uint64_t articulation = 0;
  std::uint64_t velocity = 0;
  std::uint64_t duration = 0;
  std::uint64_t gate = 0;
  std::uint64_t slide = 0;
  std::uint64_t ratchet = 0;
  std::uint64_t offset = 0;
  std::array<std::uint64_t, CvLaneCount> cv{};
  std::uint64_t structuralCell = 0;
  std::uint64_t structuralCellCount = 0;
  std::uint64_t completedCycles = 0;
  // Score time accumulated by independently cycling lanes before the current
  // Notes pass. partStartBeat_ supplies the within-pass portion.
  double freeLaneBeat = 0.0;
  double lastBaseDuration = 0.0;
  double pitchTimelineBeat = 0.0;
  std::size_t activePitchTimelineStep = std::numeric_limits<std::size_t>::max();
  bool pitchEntryConsumed = false;
  bool hasSoundingPitch = false;
  std::size_t soundingVoiceCount = 0;
  // Contextual chord voicing is kept independently from sounding polyphony so
  // a slash bass does not participate in upper-structure voice leading.
  std::array<int, MaximumPolyphony> previousVoicingSemitones{};
  std::size_t previousVoicingCount = 0;
};

struct SemanticProgram {
  std::vector<Sequence> sequences;
  std::vector<ArrangementPart> arrangement;
  bool loopArrangement = true;
  std::uint64_t seed = 1;
};

struct CompiledProgramFactory;

// A prepared program extends the checked musical graph with bounded mutable
// storage. Only this representation may be published to the audio thread.
class CompiledProgram {
public:
  CompiledProgram(const CompiledProgram &) = delete;
  CompiledProgram &operator=(const CompiledProgram &) = delete;

  const SemanticProgram &semantic() const noexcept { return semantic_; }

  std::size_t maximumEventsPerStep = 1;
  std::size_t scheduleCapacity = 1;
  double maximumEarlyBeats = 0.0;
  double maximumEarlyMilliseconds = 0.0;
  std::vector<std::size_t> stateTransferOrder;
  // These buffers are sized while compiling on the UI thread and are the only
  // mutable part of the otherwise immutable program. Playback only overwrites
  // existing elements; it never changes capacity or allocates.
  mutable std::vector<RuntimeEvent> stepWorkspace;
  mutable std::vector<RuntimeEvent> scheduleWorkspace;
  mutable std::vector<SequencePlaybackState> stateWorkspace;
  mutable std::vector<SequencePlaybackState> snapshotStateWorkspace;

private:
  explicit CompiledProgram(SemanticProgram &&semantic)
      : semantic_(std::move(semantic)) {}

  SemanticProgram semantic_;

  friend struct CompiledProgramFactory;
};

struct CompileResult {
  std::unique_ptr<CompiledProgram> program;
  Diagnostic diagnostic;

  explicit operator bool() const noexcept { return static_cast<bool>(program); }
};

CompileResult Compile(const std::string &source);
CompileResult Compile(const syntax::Document &document);

double SchedulingLookaheadBeats(const CompiledProgram &program,
                                bool periodKnown, double periodSamples,
                                double sampleRateHz) noexcept;

struct RuntimeEvent {
  EventKind kind = EventKind::Rest;
  double beat = 0.0;
  std::uint64_t scheduleOrder = 0;
  double spanBeats = 1.0;
  float pitchVolts = 0.f;
  float velocity = 0.72f;
  float accent = 0.f;
  float gateFraction = 0.8f;
  float gateBeats = -1.f;
  float gateMilliseconds = -1.f;
  float gateCapMilliseconds = -1.f;
  float slideBeats = 0.25f;
  float slideMilliseconds = -1.f;
  std::uint8_t voice = 0;
  std::uint8_t voiceCount = 1;
  bool legatoToNext = false;
  double timingOffsetBeats = 0.0;
  double timingOffsetMilliseconds = 0.0;
  std::array<float, CvLaneCount> cvValue{};
  std::array<float, CvLaneCount> cvTarget{};
  std::array<CvInterpolation, CvLaneCount> cvInterpolation{
      CvInterpolation::Step, CvInterpolation::Step, CvInterpolation::Step};
  std::array<float, CvLaneCount> cvPower{1.f, 1.f, 1.f};
  std::array<double, CvLaneCount> cvTargetBeat{};
  std::array<CvEnvelopeSpec, CvLaneCount> cvEnvelope{};
  std::array<SourceSpan, static_cast<std::size_t>(CursorLane::Count)> cursors{};
};

struct StepEvents {
  RuntimeEvent *events = nullptr;
  std::size_t capacity = 0;
  std::size_t count = 0;
  double durationBeats = 1.0;
  bool overflowed = false;
};

class Runtime {
public:
  enum class ReplacementResult { PreservedPhase, RestartedCurrentTerm };

  void setProgram(const CompiledProgram *program) noexcept;
  bool canPreserveCurrentPhase(const CompiledProgram *program) const noexcept;
  ReplacementResult replaceProgram(const CompiledProgram *program,
                                   double beat) noexcept;
  void reset() noexcept;
  void snapshot(Runtime &destination) const noexcept;
  StepEvents next(double beat) noexcept;
  bool hasProgram() const noexcept { return program_ != nullptr; }

private:
  const CompiledProgram *program_ = nullptr;
  SequencePlaybackState *states_ = nullptr;
  std::size_t stateCount_ = 0;
  std::size_t arrangementPart_ = 0;
  int arrangementCycle_ = 0;
  double partStartBeat_ = 0.0;

  const Sequence *currentSequence(SourceSpan &partSpan,
                                  std::uint64_t &cycle) noexcept;
};

} // namespace tfseq
