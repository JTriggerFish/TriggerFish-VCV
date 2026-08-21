#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace tfseq {

namespace syntax {
struct Document;
}

constexpr std::size_t MaximumPolyphony = 16;
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
  Octave,
  Articulation,
  Velocity,
  Accent,
  Duration,
  Gate,
  Slide,
  Ratchet,
  Count
};

enum class EventKind { Attack, Slide, Tie, Rest };

enum class ArticulationKind { Attack, Slide, Tie, Rest };

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

struct ChordValue {
  enum class Meaning { SinglePitch, ExplicitVoicing, JazzSymbol };
  Meaning meaning = Meaning::SinglePitch;
  // Preserve harmonic intent beside the current default realization so a
  // future interpreter can voice or arpeggiate the chord contextually.
  std::string jazzSymbol;
  int rootPitchClass = 0;
  std::vector<PitchValue> voices;
  bool hasBass = false;
  PitchValue bass;
  SourceSpan span;
};

struct PitchItem {
  enum class Choice { Single, Alternate, Random };
  Choice choice = Choice::Single;
  std::vector<ChordValue> values;
  SourceSpan span;
};

struct ArticulationAtom {
  ArticulationKind kind = ArticulationKind::Attack;
  std::size_t ratchets = 1;
  float probability = 1.f;
  SourceSpan span;
};

struct ArticulationStep {
  std::vector<ArticulationAtom> atoms;
  SourceSpan span;
};

struct ScalarItem {
  double value = 0.0;
  bool isDefault = false;
  SourceSpan span;
};

enum class TransformKind {
  Reverse,
  Rotate,
  ModulateDegree,
  ShiftDegree,
  TransposeSemitone,
  TransposeOctave,
  Fast,
  Slow,
  Swing,
  Early,
  Late
};

enum class TransformCondition { Always, Every, Sometimes };
enum class TransformDomain { General, Pitch, Timing };
enum class TimeUnit { Beats, Milliseconds };

struct Transform {
  TransformKind kind = TransformKind::Reverse;
  TransformDomain domain = TransformDomain::General;
  TransformCondition condition = TransformCondition::Always;
  int integer = 0;
  double number = 0.0;
  // Zero selects the event's own subdivision. A positive value supplies an
  // explicit grid in incoming-clock beats (for example 1/8 or 1/16).
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
  double cycleBeats = 8.0;
  float glideBeats = 0.25f;
  Scale scale;
  std::vector<PitchItem> notes;
  std::vector<ScalarItem> octave;
  std::vector<ArticulationStep> articulation;
  std::vector<ScalarItem> velocity;
  std::vector<ScalarItem> accent;
  std::vector<ScalarItem> duration;
  std::vector<ScalarItem> gate;
  std::vector<ScalarItem> slide;
  std::vector<ScalarItem> ratchet;
  std::array<std::vector<Transform>,
             static_cast<std::size_t>(CursorLane::Count)>
      transforms;
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
  std::uint64_t accent = 0;
  std::uint64_t duration = 0;
  std::uint64_t gate = 0;
  std::uint64_t slide = 0;
  std::uint64_t ratchet = 0;
};

struct SemanticProgram {
  std::vector<Sequence> sequences;
  std::vector<ArrangementPart> arrangement;
  bool loopArrangement = true;
  bool stopped = false;
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
  double spanBeats = 1.0;
  float pitchVolts = 0.f;
  float velocity = 0.72f;
  float accent = 0.f;
  float gateFraction = 0.8f;
  float slideBeats = 0.25f;
  std::uint8_t voice = 0;
  std::uint8_t voiceCount = 1;
  bool legatoToNext = false;
  double timingOffsetMilliseconds = 0.0;
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
  void setProgram(const CompiledProgram *program) noexcept;
  void replaceProgram(const CompiledProgram *program, double beat) noexcept;
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

  const Sequence *currentSequence(double beat, SourceSpan &partSpan,
                                  std::uint64_t &cycle) noexcept;
};

} // namespace tfseq
