#pragma once

#include "tfseq.hpp"

#include <memory>
#include <string>
#include <variant>
#include <vector>

namespace tfseq::syntax {

struct Token {
  std::string text;
  SourceSpan span;
};

struct Pipeline {
  enum class Condition { Always, Every, Sometimes };
  Token operation;
  std::vector<Token> arguments;
  Condition condition = Condition::Always;
  Token conditionArgument;
  SourceSpan span;
};

enum class PatternKind {
  Atom,
  NamedPitch,
  ScaleDegree,
  JazzChord,
  RomanChord,
  ChordFactor,
  RandomPitch,
  RandomScalar,
  Event,
  Rest,
  Tie,
  Subdivision,
  CycleChoice,
  RandomChoice,
  Voicing,
  Slash
};

// Typed pattern structure. The PEG owns event, grouping, choice, chord,
// suffix, and attribute boundaries; semantic lowering converts validated leaf
// values without running another delimiter or expression parser.
struct PatternNode {
  PatternKind kind = PatternKind::Atom;
  Token atom;
  Token suffix;
  Token repeatCount;
  Token slidePrefix;
  Token dynamicPrefix;
  Token durationSuffix;
  Token ratchetCount;
  Token probability;
  bool defaultProbability = false;
  Token presenceProbability;
  bool defaultPresenceProbability = false;
  struct Attribute {
    Token name;
    Token value;
  };
  std::vector<Attribute> attributes;
  std::vector<Token> arguments;
  std::vector<PatternNode> children;
  SourceSpan span;
};

struct Pattern {
  enum class Alignment { Free, Left, Right, Edges };
  std::vector<PatternNode> steps;
  SourceSpan span;
  Alignment alignment = Alignment::Free;
  // Number of source terms before a middle ellipsis. The compiler expands
  // repetitions when it turns this into the runtime edge split.
  std::size_t alignmentSplit = 0;
};

struct Lane {
  enum class Kind { Notes, Chords, Rhythm, Percussion, Scalar, Cv, Setting };
  Token name;
  Pattern pattern;
  std::vector<Token> envelopeArguments;
  SourceSpan envelopeSpan;
  bool envelopeOnly = false;
  std::vector<Pipeline> pipelines;
  // A rhythm lane either contains an inline pattern or names a reusable
  // top-level rhythm definition. Other lane kinds leave this empty.
  Token rhythmReference;
  Kind kind = Kind::Scalar;
};

struct RhythmDefinition {
  Token name;
  Token subdivision;
  Pattern events;
  std::vector<Pipeline> pipelines;
  SourceSpan span;
};

struct SequenceDefinition {
  Token name;
  std::vector<Lane> lanes;
  std::vector<Pipeline> pipelines;
  SourceSpan span;
};

struct Expression;

struct AssignmentTerm {
  Token name;
  std::shared_ptr<Expression> group;
  int repeats = 1;
  bool hasRepeat = false;
  SourceSpan repeatSpan;
};

struct Expression {
  std::vector<AssignmentTerm> terms;
  std::vector<Pipeline> pipelines;
  SourceSpan span;
};

struct Assignment {
  Token name;
  Expression expression;
  SourceSpan span;

  bool isDerived() const noexcept {
    return expression.terms.size() == 1 && !expression.terms.front().group &&
           !expression.terms.front().hasRepeat;
  }
};

struct PlayCommand {
  Token name;
  SourceSpan span;
};

struct SeedCommand {
  Token value;
  SourceSpan span;
};

using Statement =
    std::variant<SequenceDefinition, RhythmDefinition, Assignment, PlayCommand,
                 SeedCommand>;

struct Document {
  std::vector<Statement> statements;
};

struct ParseResult {
  Document document;
  Diagnostic diagnostic;

  explicit operator bool() const noexcept { return diagnostic.message.empty(); }
};

ParseResult Parse(const std::string &source);

// Parses the smallest complete line-bounded top-level statement set that
// overlaps the selection. This uses the same PEG as Parse() and allows live
// evaluation to ignore unrelated malformed draft text.
ParseResult ParseStatementsContaining(const std::string &source,
                                       int selectionBegin,
                                       int selectionEnd);

struct SelectionDocumentResult {
  Document document;
  Diagnostic diagnostic;

  explicit operator bool() const noexcept { return diagnostic.message.empty(); }
};

// Produces a typed evaluated document by replacing only selected top-level
// statements. Retained definitions whose draft text has changed have their
// source spans invalidated so playback never highlights inactive draft text.
SelectionDocumentResult MergeSelectionDocuments(
    const Document &evaluatedDocument, const std::string &evaluatedSource,
    const Document &draftDocument, const std::string &draftSource,
    int selectionBegin, int selectionEnd);

} // namespace tfseq::syntax
