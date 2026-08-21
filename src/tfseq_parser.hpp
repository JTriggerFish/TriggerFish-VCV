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
  Subdivision,
  CycleChoice,
  RandomChoice,
  Euclidean,
  Voicing,
  Slash,
  Repeat
};

// Domain-neutral pattern syntax. Musical meaning is assigned by the semantic
// compiler, while grouping/choice/repetition is parsed exactly once here.
struct PatternNode {
  PatternKind kind = PatternKind::Atom;
  Token atom;
  Token suffix;
  Token repeatCount;
  std::vector<Token> arguments;
  std::vector<PatternNode> children;
  SourceSpan span;
};

struct Pattern {
  std::vector<PatternNode> steps;
  SourceSpan span;
};

struct Lane {
  Token name;
  Pattern pattern;
  std::vector<Pipeline> pipelines;
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

struct StopCommand {
  SourceSpan span;
};

struct SeedCommand {
  Token value;
  SourceSpan span;
};

using Statement = std::variant<SequenceDefinition, Assignment, PlayCommand,
                               StopCommand, SeedCommand>;

struct Document {
  std::vector<Statement> statements;
};

struct ParseResult {
  Document document;
  Diagnostic diagnostic;

  explicit operator bool() const noexcept { return diagnostic.message.empty(); }
};

ParseResult Parse(const std::string &source);

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
