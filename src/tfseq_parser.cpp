#include "tfseq_parser.hpp"

#include <peglib.h>

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string_view>
#include <type_traits>
#include <unordered_set>

namespace tfseq::syntax {
namespace {

constexpr const char *Grammar = R"PEG(
  Document    <- (Sequence / Play / Stop / Seed / Assignment / BlankLine)* Trailing End

  Sequence    <- H Identifier H '=' H 'sequence' H '{' H Comment? Newline
                 (Lane / BodyContinuation / BlankLine)*
                 H '}' H Comment? LineEnd SequenceContinuation*
  Play        <- H 'play' S Identifier H Comment? LineEnd
  Stop        <- H < 'stop' > H Comment? LineEnd
  Seed        <- H 'seed' S ValueText H Comment? LineEnd
  Assignment  <- H Identifier H '=' H Expression H Comment? LineEnd

  Lane         <- H Identifier S Pattern (H Pipeline)* H Comment? LineEnd
  BodyContinuation     <- H Pipeline (H Pipeline)* H Comment? LineEnd
  SequenceContinuation <- H Pipeline (H Pipeline)* H Comment? LineEnd
  Pattern        <- PatternElement (S PatternElement)*
  PatternElement <- (CycleChoice / BracketGroup / Euclidean / TopValue) Repeat?

  CycleChoice <- '<' H CycleElement (S CycleElement)* H '>'
  CycleElement <- (CycleChoice / BracketGroup / Euclidean / CycleValue) Repeat?
  CycleValue   <- (Voicing / CycleAtom) SlashSuffix?
  CycleAtom    <- < (Paren / CycleAtomChar)+ >
  CycleAtomChar <- !('[' / '<' / '(' / '/' / ']' / '>' / ')' / '!' / '|' /
                     HChar / CommentStart / Pipe / Newline) .

  BracketGroup   <- '[' H (RandomChoice / Subdivision) H ']'
  RandomChoice   <- BracketElement H '|' H BracketElement
                    (H '|' H BracketElement)*
  Subdivision    <- BracketElement (S BracketElement)*
  BracketElement <- (CycleChoice / BracketGroup / Euclidean / SlideAtom /
                     BracketValue) Repeat?
  BracketValue   <- (Voicing / BracketAtom) SlashSuffix?
  SlideAtom      <- < '>' >
  BracketAtom    <- < (Paren / BracketAtomChar)+ >
  BracketAtomChar <- !('[' / '<' / '(' / '/' / ']' / '>' / ')' / '!' / '|' /
                       HChar / CommentStart / Pipe / Newline) .

  Euclidean         <- EuclideanOperator '(' H EuclideanArgument
                       ((H ',' H / S) EuclideanArgument)* H ')'
  EuclideanOperator <- < 'x' / '.' >
  EuclideanArgument <- < (!HChar !',' !')' !Newline .)+ >
  Repeat            <- '!' RepeatCount
  RepeatCount       <- < (!HChar !CommentStart !Pipe !Newline !']' !'>' !'|'
                          !'!' .)+ >
  TopValue          <- (Voicing / TopAtom) SlashSuffix?
  TopAtom           <- < TopAtomChar+ >
  TopAtomChar       <- !('[' / '<' / '(' / '/' / ']' / ')' / '!' / '|' /
                         HChar / CommentStart / Pipe / Newline) .
  SlashSuffix       <- H '/' H SlashBass
  SlashBass         <- < (!HChar !'(' !')' !'/' !'[' !']' !'<' !'>' !'!' !'|'
                         !CommentStart !Pipe !Newline .)+ >
  Voicing           <- '(' H VoicingTone (S VoicingTone)* H VoicingSlash? H ')'
                       RegisterSuffix?
  VoicingSlash      <- '/' H VoicingTone
  VoicingTone       <- < (!HChar !'(' !')' !'/' !'[' !']' !'<' !'>' !'!' !'|'
                         !CommentStart !Pipe !Newline .)+ >
  RegisterSuffix    <- < ('@' [+-]? [0-9]+ [',]*) / [',]+ >
  Paren             <- '(' (Paren / (!')' !Newline .))* ')'

  Expression <- Concatenation (H Pipeline)*
  Concatenation <- Term (H '+' H Term)*
  Term       <- Primary (H '*' H Integer)?
  Primary    <- Identifier / '(' H Expression H ')'
  Integer    <- < [0-9]+ >

  Pipeline      <- Pipe H (EveryTransform / SometimesTransform / TransformCall)
  EveryTransform <- 'every' S ConditionArgument S TransformCall
  SometimesTransform <- 'sometimes' S ConditionArgument S TransformCall
  TransformCall <- '(' H Operation (S TransformArgument)* H ')' /
                   Operation (S TransformArgument)*
  Operation     <- < [A-Za-z_] [A-Za-z0-9_]* >
  ConditionArgument <- < (!HChar !Pipe !CommentStart !Newline !')' .)+ >
  TransformArgument <- < (!HChar !Pipe !CommentStart !Newline !')' .)+ >
  ValueText     <- < (!CommentStart !Newline .)+ >
  Identifier    <- < [A-Za-z_] [A-Za-z0-9_]* >

  ~BlankLine   <- H Comment? Newline
  ~Trailing    <- H Comment?
  ~Comment     <- CommentStart (!Newline .)*
  ~CommentStart <- '//'
  ~Pipe        <- '|>'
  ~H           <- HChar*
  ~S           <- HChar+
  ~HChar       <- [ \t]
  ~Newline     <- '\r\n' / '\n' / '\r'
  ~LineEnd     <- Newline / End
  ~End         <- !.
)PEG";

using Ast = std::shared_ptr<peg::Ast>;

SourceSpan Span(const Ast &node) {
  return {static_cast<int>(node->position),
          static_cast<int>(node->position + node->length),
          static_cast<int>(node->line), static_cast<int>(node->column)};
}

Token AstToken(const Ast &node) {
  return {node->token_to_string(), Span(node)};
}

Token TrimmedAstToken(const Ast &node) {
  std::string_view value = node->token;
  std::size_t leading = 0;
  while (leading < value.size() &&
         (value[leading] == ' ' || value[leading] == '\t'))
    ++leading;
  std::size_t trailing = value.size();
  while (trailing > leading &&
         (value[trailing - 1] == ' ' || value[trailing - 1] == '\t'))
    --trailing;
  Token token;
  token.text = std::string(value.substr(leading, trailing - leading));
  token.span = Span(node);
  token.span.begin += static_cast<int>(leading);
  token.span.end = token.span.begin + static_cast<int>(token.text.size());
  token.span.column += static_cast<int>(leading);
  return token;
}

std::vector<Ast> ChildrenNamed(const Ast &node, const std::string &name) {
  std::vector<Ast> children;
  for (const auto &child : node->nodes) {
    if (child->name == name)
      children.push_back(child);
  }
  return children;
}

Ast FirstChildNamed(const Ast &node, const std::string &name) {
  for (const auto &child : node->nodes) {
    if (child->name == name)
      return child;
  }
  return nullptr;
}

PatternNode ReadPatternElement(const Ast &element, Diagnostic &diagnostic) {
  PatternNode node;
  if (const auto choice = FirstChildNamed(element, "CycleChoice")) {
    node.kind = PatternKind::CycleChoice;
    node.span = Span(choice);
    for (const auto &child : ChildrenNamed(choice, "CycleElement"))
      node.children.push_back(ReadPatternElement(child, diagnostic));
  } else if (const auto group = FirstChildNamed(element, "BracketGroup")) {
    node.span = Span(group);
    if (const auto random = FirstChildNamed(group, "RandomChoice")) {
      node.kind = PatternKind::RandomChoice;
      for (const auto &child : ChildrenNamed(random, "BracketElement"))
        node.children.push_back(ReadPatternElement(child, diagnostic));
    } else if (const auto subdivision = FirstChildNamed(group, "Subdivision")) {
      node.kind = PatternKind::Subdivision;
      for (const auto &child : ChildrenNamed(subdivision, "BracketElement"))
        node.children.push_back(ReadPatternElement(child, diagnostic));
    }
  } else if (const auto euclidean = FirstChildNamed(element, "Euclidean")) {
    node.kind = PatternKind::Euclidean;
    node.span = Span(euclidean);
    node.atom = AstToken(FirstChildNamed(euclidean, "EuclideanOperator"));
    for (const auto &argument : ChildrenNamed(euclidean, "EuclideanArgument"))
      node.arguments.push_back(AstToken(argument));
  } else {
    Ast value = FirstChildNamed(element, "TopValue");
    if (!value)
      value = FirstChildNamed(element, "CycleValue");
    if (!value)
      value = FirstChildNamed(element, "BracketValue");

    auto readValue = [&](const Ast &source) {
      PatternNode parsed;
      if (const auto voicing = FirstChildNamed(source, "Voicing")) {
        parsed.kind = PatternKind::Voicing;
        parsed.span = Span(voicing);
        if (const auto suffix = FirstChildNamed(voicing, "RegisterSuffix"))
          parsed.suffix = AstToken(suffix);
        for (const auto &tone : ChildrenNamed(voicing, "VoicingTone")) {
          PatternNode child;
          child.kind = PatternKind::Atom;
          child.atom = AstToken(tone);
          child.span = child.atom.span;
          parsed.children.push_back(std::move(child));
        }
        if (const auto slash = FirstChildNamed(voicing, "VoicingSlash")) {
          PatternNode bass;
          bass.kind = PatternKind::Atom;
          bass.atom = AstToken(FirstChildNamed(slash, "VoicingTone"));
          bass.span = bass.atom.span;
          PatternNode joined;
          joined.kind = PatternKind::Slash;
          joined.span = Span(voicing);
          joined.children.push_back(std::move(parsed));
          joined.children.push_back(std::move(bass));
          parsed = std::move(joined);
        }
      } else {
        Ast atom = FirstChildNamed(source, "TopAtom");
        if (!atom)
          atom = FirstChildNamed(source, "CycleAtom");
        if (!atom)
          atom = FirstChildNamed(source, "BracketAtom");
        parsed.kind = PatternKind::Atom;
        parsed.atom = AstToken(atom);
        parsed.span = parsed.atom.span;
      }
      if (const auto slash = FirstChildNamed(source, "SlashSuffix")) {
        PatternNode bass;
        bass.kind = PatternKind::Atom;
        bass.atom = AstToken(FirstChildNamed(slash, "SlashBass"));
        bass.span = bass.atom.span;
        PatternNode joined;
        joined.kind = PatternKind::Slash;
        joined.span = Span(source);
        joined.children.push_back(std::move(parsed));
        joined.children.push_back(std::move(bass));
        parsed = std::move(joined);
      }
      return parsed;
    };

    if (value) {
      node = readValue(value);
    } else if (const auto slide = FirstChildNamed(element, "SlideAtom")) {
      node.kind = PatternKind::Atom;
      node.atom = AstToken(slide);
      node.span = node.atom.span;
    } else {
      diagnostic = {"internal sequencer pattern AST is incomplete",
                    static_cast<int>(element->line),
                    static_cast<int>(element->column)};
      return {};
    }
  }
  if (!diagnostic.message.empty())
    return {};

  if (const auto repeat = FirstChildNamed(element, "Repeat")) {
    PatternNode repeated;
    repeated.kind = PatternKind::Repeat;
    repeated.repeatCount = AstToken(FirstChildNamed(repeat, "RepeatCount"));
    repeated.children.push_back(std::move(node));
    repeated.span = Span(element);
    return repeated;
  }
  return node;
}

Pattern ReadPattern(const Ast &node, Diagnostic &diagnostic) {
  Pattern pattern;
  pattern.span = Span(node);
  for (const auto &element : ChildrenNamed(node, "PatternElement")) {
    pattern.steps.push_back(ReadPatternElement(element, diagnostic));
    if (!diagnostic.message.empty())
      break;
  }
  return pattern;
}

Pipeline ReadPipeline(const Ast &node) {
  Pipeline pipeline;
  pipeline.span = Span(node);
  Ast transform = FirstChildNamed(node, "TransformCall");
  if (const auto every = FirstChildNamed(node, "EveryTransform")) {
    pipeline.condition = Pipeline::Condition::Every;
    pipeline.conditionArgument =
        AstToken(FirstChildNamed(every, "ConditionArgument"));
    transform = FirstChildNamed(every, "TransformCall");
  } else if (const auto sometimes =
                 FirstChildNamed(node, "SometimesTransform")) {
    pipeline.condition = Pipeline::Condition::Sometimes;
    pipeline.conditionArgument =
        AstToken(FirstChildNamed(sometimes, "ConditionArgument"));
    transform = FirstChildNamed(sometimes, "TransformCall");
  }
  pipeline.operation = AstToken(FirstChildNamed(transform, "Operation"));
  for (const auto &argument : ChildrenNamed(transform, "TransformArgument"))
    pipeline.arguments.push_back(AstToken(argument));
  return pipeline;
}

Lane ReadLane(const Ast &node, Diagnostic &diagnostic) {
  Lane lane;
  lane.name = AstToken(FirstChildNamed(node, "Identifier"));
  const auto pattern = FirstChildNamed(node, "Pattern");
  if (pattern)
    lane.pattern = ReadPattern(pattern, diagnostic);
  for (const auto &pipeline : ChildrenNamed(node, "Pipeline"))
    lane.pipelines.push_back(ReadPipeline(pipeline));
  return lane;
}

SequenceDefinition ReadSequence(const Ast &node, Diagnostic &diagnostic) {
  SequenceDefinition sequence;
  sequence.span = Span(node);
  sequence.name = AstToken(FirstChildNamed(node, "Identifier"));
  for (const auto &child : node->nodes) {
    if (child->name == "Lane")
      sequence.lanes.push_back(ReadLane(child, diagnostic));
    else if (child->name == "BodyContinuation") {
      if (sequence.lanes.empty())
        continue;
      for (const auto &pipeline : ChildrenNamed(child, "Pipeline"))
        sequence.lanes.back().pipelines.push_back(ReadPipeline(pipeline));
    } else if (child->name == "SequenceContinuation") {
      for (const auto &pipeline : ChildrenNamed(child, "Pipeline"))
        sequence.pipelines.push_back(ReadPipeline(pipeline));
    }
  }
  return sequence;
}

Expression ReadExpression(const Ast &node);

AssignmentTerm ReadTerm(const Ast &node) {
  AssignmentTerm term;
  const auto primary = FirstChildNamed(node, "Primary");
  if (const auto name = FirstChildNamed(primary, "Identifier"))
    term.name = AstToken(name);
  else if (const auto expression = FirstChildNamed(primary, "Expression"))
    term.group = std::make_shared<Expression>(ReadExpression(expression));
  const auto repeat = FirstChildNamed(node, "Integer");
  if (repeat) {
    const auto text = repeat->token_to_string();
    char *end = nullptr;
    const auto parsed = std::strtoull(text.c_str(), &end, 10);
    term.repeats = parsed > static_cast<unsigned long long>(
                                std::numeric_limits<int>::max())
                       ? std::numeric_limits<int>::max()
                       : static_cast<int>(parsed);
    term.hasRepeat = true;
    term.repeatSpan = Span(repeat);
  }
  return term;
}

Expression ReadExpression(const Ast &node) {
  Expression expression;
  expression.span = Span(node);
  const auto concatenation = FirstChildNamed(node, "Concatenation");
  for (const auto &term : ChildrenNamed(concatenation, "Term"))
    expression.terms.push_back(ReadTerm(term));
  for (const auto &pipeline : ChildrenNamed(node, "Pipeline"))
    expression.pipelines.push_back(ReadPipeline(pipeline));
  return expression;
}

Assignment ReadAssignment(const Ast &node) {
  Assignment assignment;
  assignment.span = Span(node);
  assignment.name = AstToken(FirstChildNamed(node, "Identifier"));
  assignment.expression = ReadExpression(FirstChildNamed(node, "Expression"));
  return assignment;
}

} // namespace

ParseResult Parse(const std::string &source) {
  ParseResult result;
  peg::parser parser;
  parser.set_logger(
      [&](std::size_t line, std::size_t column, const std::string &message) {
        if (result.diagnostic.message.empty())
          result.diagnostic = {message, static_cast<int>(line),
                               static_cast<int>(column)};
      });
  if (!parser.load_grammar(Grammar)) {
    if (result.diagnostic.message.empty())
      result.diagnostic = {"internal sequencer grammar is invalid", 1, 1};
    else
      result.diagnostic.message =
          "internal sequencer grammar: " + result.diagnostic.message;
    return result;
  }
  parser.enable_ast();

  Ast ast;
  if (!parser.parse(source, ast)) {
    if (result.diagnostic.message.empty())
      result.diagnostic = {"invalid sequencer syntax", 1, 1};
    return result;
  }

  for (const auto &node : ast->nodes) {
    if (node->name == "Sequence") {
      result.document.statements.emplace_back(
          ReadSequence(node, result.diagnostic));
      if (!result.diagnostic.message.empty())
        return result;
    } else if (node->name == "Assignment") {
      result.document.statements.emplace_back(ReadAssignment(node));
    } else if (node->name == "Play") {
      result.document.statements.emplace_back(PlayCommand{
          AstToken(FirstChildNamed(node, "Identifier")), Span(node)});
    } else if (node->name == "Stop") {
      result.document.statements.emplace_back(StopCommand{Span(node)});
    } else if (node->name == "Seed") {
      result.document.statements.emplace_back(SeedCommand{
          TrimmedAstToken(FirstChildNamed(node, "ValueText")), Span(node)});
    }
  }
  return result;
}

namespace {

SourceSpan StatementSpan(const Statement &statement) {
  return std::visit([](const auto &value) { return value.span; }, statement);
}

std::string DefinitionKey(const Statement &statement) {
  if (const auto *sequence = std::get_if<SequenceDefinition>(&statement))
    return sequence->name.text;
  if (const auto *assignment = std::get_if<Assignment>(&statement))
    return assignment->name.text;
  return {};
}

bool IsTransport(const Statement &statement) {
  return std::holds_alternative<PlayCommand>(statement) ||
         std::holds_alternative<StopCommand>(statement);
}

bool IsSeed(const Statement &statement) {
  return std::holds_alternative<SeedCommand>(statement);
}

std::string StatementText(const std::string &source,
                          const Statement &statement) {
  const auto span = StatementSpan(statement);
  if (!span.valid() || span.begin < 0 || span.end < span.begin ||
      static_cast<std::size_t>(span.end) > source.size())
    return {};
  return source.substr(static_cast<std::size_t>(span.begin),
                       static_cast<std::size_t>(span.end - span.begin));
}

void Invalidate(Token &token) { token.span = {}; }
void Invalidate(SourceSpan &span) { span = {}; }

void Invalidate(Pipeline &pipeline) {
  Invalidate(pipeline.operation);
  for (auto &argument : pipeline.arguments)
    Invalidate(argument);
  Invalidate(pipeline.conditionArgument);
  pipeline.span = {};
}

void Invalidate(PatternNode &node) {
  Invalidate(node.atom);
  Invalidate(node.suffix);
  Invalidate(node.repeatCount);
  for (auto &argument : node.arguments)
    Invalidate(argument);
  for (auto &child : node.children)
    Invalidate(child);
  node.span = {};
}

void Invalidate(Pattern &pattern) {
  for (auto &step : pattern.steps)
    Invalidate(step);
  pattern.span = {};
}

void Invalidate(Expression &expression) {
  for (auto &term : expression.terms) {
    Invalidate(term.name);
    Invalidate(term.repeatSpan);
    if (term.group) {
      term.group = std::make_shared<Expression>(*term.group);
      Invalidate(*term.group);
    }
  }
  for (auto &pipeline : expression.pipelines)
    Invalidate(pipeline);
  expression.span = {};
}

void Invalidate(Statement &statement) {
  std::visit(
      [](auto &value) {
        using Value = std::decay_t<decltype(value)>;
        if constexpr (std::is_same_v<Value, SequenceDefinition>) {
          Invalidate(value.name);
          for (auto &lane : value.lanes) {
            Invalidate(lane.name);
            Invalidate(lane.pattern);
            for (auto &pipeline : lane.pipelines)
              Invalidate(pipeline);
          }
          for (auto &pipeline : value.pipelines)
            Invalidate(pipeline);
        } else if constexpr (std::is_same_v<Value, Assignment>) {
          Invalidate(value.name);
          Invalidate(value.expression);
        } else if constexpr (std::is_same_v<Value, PlayCommand>) {
          Invalidate(value.name);
        } else if constexpr (std::is_same_v<Value, SeedCommand>) {
          Invalidate(value.value);
        }
        value.span = {};
      },
      statement);
}

bool SameIdentity(const Statement &left, const Statement &right) {
  const auto leftDefinition = DefinitionKey(left);
  if (!leftDefinition.empty())
    return leftDefinition == DefinitionKey(right);
  if (IsTransport(left))
    return IsTransport(right);
  return IsSeed(left) && IsSeed(right);
}

} // namespace

SelectionDocumentResult MergeSelectionDocuments(
    const Document &evaluatedDocument, const std::string &evaluatedSource,
    const Document &draftDocument, const std::string &draftSource,
    int selectionBegin, int selectionEnd) {
  SelectionDocumentResult result;
  const int sourceSize = static_cast<int>(std::min<std::size_t>(
      draftSource.size(),
      static_cast<std::size_t>(std::numeric_limits<int>::max())));
  selectionBegin = std::clamp(selectionBegin, 0, sourceSize);
  selectionEnd = std::clamp(selectionEnd, 0, sourceSize);
  if (selectionEnd < selectionBegin)
    std::swap(selectionBegin, selectionEnd);

  std::vector<std::size_t> selected;
  std::unordered_set<std::string> replacedDefinitions;
  bool replaceTransport = false;
  bool replaceSeed = false;
  for (std::size_t index = 0; index < draftDocument.statements.size();
       ++index) {
    const auto &statement = draftDocument.statements[index];
    const auto span = StatementSpan(statement);
    const bool overlaps =
        selectionBegin == selectionEnd
            ? span.begin <= selectionBegin && selectionBegin < span.end
            : span.begin < selectionEnd && selectionBegin < span.end;
    if (!overlaps)
      continue;
    selected.push_back(index);
    const auto key = DefinitionKey(statement);
    if (!key.empty())
      replacedDefinitions.insert(key);
    replaceTransport = replaceTransport || IsTransport(statement);
    replaceSeed = replaceSeed || IsSeed(statement);
  }
  if (selected.empty()) {
    int line = 1;
    int column = 1;
    for (int index = 0; index < selectionBegin; ++index) {
      if (draftSource[static_cast<std::size_t>(index)] == '\n') {
        ++line;
        column = 1;
      } else {
        ++column;
      }
    }
    result.diagnostic = {"selection does not contain an executable statement",
                         line, column};
    return result;
  }

  std::unordered_set<std::size_t> usedDraft;
  auto appendSelectedMatching = [&](const Statement &evaluatedStatement) {
    for (const auto index : selected) {
      if (usedDraft.count(index) != 0)
        continue;
      const auto &candidate = draftDocument.statements[index];
      if (SameIdentity(evaluatedStatement, candidate)) {
        result.document.statements.push_back(candidate);
        usedDraft.insert(index);
      }
    }
  };

  bool insertedTransport = false;
  bool insertedSeed = false;
  for (const auto &evaluatedStatement : evaluatedDocument.statements) {
    const auto key = DefinitionKey(evaluatedStatement);
    if (!key.empty() && replacedDefinitions.count(key) != 0) {
      appendSelectedMatching(evaluatedStatement);
      continue;
    }
    if (replaceTransport && IsTransport(evaluatedStatement)) {
      if (!insertedTransport) {
        appendSelectedMatching(evaluatedStatement);
        insertedTransport = true;
      }
      continue;
    }
    if (replaceSeed && IsSeed(evaluatedStatement)) {
      if (!insertedSeed) {
        appendSelectedMatching(evaluatedStatement);
        insertedSeed = true;
      }
      continue;
    }

    bool usedCurrentSpans = false;
    const auto evaluatedText =
        StatementText(evaluatedSource, evaluatedStatement);
    for (std::size_t index = 0; index < draftDocument.statements.size();
         ++index) {
      if (usedDraft.count(index) != 0)
        continue;
      const auto &candidate = draftDocument.statements[index];
      if (SameIdentity(evaluatedStatement, candidate) &&
          StatementText(draftSource, candidate) == evaluatedText) {
        result.document.statements.push_back(candidate);
        usedDraft.insert(index);
        usedCurrentSpans = true;
        break;
      }
    }
    if (!usedCurrentSpans) {
      auto retained = evaluatedStatement;
      Invalidate(retained);
      result.document.statements.push_back(std::move(retained));
    }
  }

  for (const auto index : selected) {
    if (usedDraft.count(index) == 0)
      result.document.statements.push_back(draftDocument.statements[index]);
  }
  return result;
}

} // namespace tfseq::syntax
