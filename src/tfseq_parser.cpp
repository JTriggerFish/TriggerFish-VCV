#include "tfseq_parser.hpp"

#include <peglib.h>

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <memory>
#include <type_traits>
#include <unordered_set>

namespace tfseq::syntax {
namespace {

constexpr const char *Grammar = R"PEG(
  Document    <- (Sequence / Play / Stop / Seed / Assignment / BlankLine)* Trailing End

  Sequence    <- H Identifier H '=' H 'sequence' H '{' H Comment? Newline
                 (NotesLane / ScalarLane / CvLane / SettingLine /
                  BodyContinuation / BlankLine)*
                 H '}' H Comment? LineEnd SequenceContinuation*
  Play        <- H 'play' S Identifier H Comment? LineEnd
  Stop        <- H < 'stop' > H Comment? LineEnd
  Seed        <- H 'seed' S UnsignedInteger H Comment? LineEnd
  Assignment  <- H Identifier H '=' H Expression H Comment? LineEnd
                 AssignmentContinuation*

  NotesLane    <- H NotesLaneName S NotePattern (H Pipeline)* H Comment? LineEnd
  ScalarLane   <- H ScalarLaneName S ScalarPattern (H Pipeline)* H Comment? LineEnd
  CvLane       <- H CvLaneName S ScalarPattern (H Pipeline)* H Comment? LineEnd
  SettingLine  <- H SettingName S SettingValue H Comment? LineEnd
  NotesLaneName  <- < 'notes' >
  ScalarLaneName <- < 'octave' / 'velocity' / 'vel' / 'accent' /
                      'duration' / 'dur' / 'gate' / 'slide' / 'ratchet' /
                      'offset' >
  CvLaneName     <- < 'cv' PositiveInteger >
  SettingName    <- < 'cycle' / 'tonic' / 'scale' / 'glide' >
  BodyContinuation     <- H Pipeline (H Pipeline)* H Comment? LineEnd
  SequenceContinuation <- H Pipeline (H Pipeline)* H Comment? LineEnd
  AssignmentContinuation <- H Pipeline (H Pipeline)* H Comment? LineEnd

  NotePattern     <- NoteElement (S NoteElement)*
  NoteElement     <- GroupElement / NoteEvent / RestEvent / TieExtension
  GroupElement    <- GroupPrimary DurationSuffix? EuclideanSuffix?
                     ProbabilitySuffix? ReplicationSuffix?
  GroupPrimary    <- BracketGroup / Alternate
  BracketGroup    <- '[' H (RandomChoice / NotePattern) H ']'
  RandomChoice    <- NotePattern H '|' H NotePattern
                     (H '|' H NotePattern)*
  Alternate       <- '<' H NoteElement (S NoteElement)* H '>'

  NoteEvent       <- OnsetPrefix? PitchedValue DurationSuffix?
                     EuclideanSuffix? RatchetSuffix? ProbabilitySuffix?
                     ReplicationSuffix? Attributes?
  RestEvent       <- RestMark DurationSuffix? ReplicationSuffix? Attributes?
  TieExtension    <- TieMark ReplicationSuffix?
  RestMark        <- < '~' >
  TieMark         <- < '_' >

  OnsetPrefix     <- SlidePrefix? DynamicPrefix?
  SlidePrefix     <- < '>' >
  DynamicPrefix   <- < '^^' / '^' / 'x' >

  PitchedValue    <- RandomPitch / ChordValue (H SlashSuffix)? / PitchValue
  RandomPitch     <- '$' (PitchRandomDistribution? RandomArguments)?
  PitchRandomDistribution <- < 'cn' / 'c' / 'n' / 'u' >
  ChordValue      <- (ExplicitVoicing / RomanChord / JazzChord) RegisterSuffix?
  PitchValue      <- (NamedPitch / ScaleDegree) RegisterSuffix?
  ExplicitVoicing <- '(' H PitchValue (S PitchValue)+ H ')'
  SlashSuffix     <- '/' H SlashBass
  SlashBass       <- PitchValue

  RomanChord      <- < RootAccidental* (UpperRoman / LowerRoman) ChordSuffix? >
  JazzChord       <- < NamedRoot ChordSuffix >
  ChordSuffix     <- TriadQuality ChordExtension? ChordModification* /
                     ChordExtension ChordModification* / ChordModification+
  TriadQuality    <- 'sus2' / 'sus4' / 'maj' / 'min' / 'dim' / 'aug' / 'm'
  ChordExtension  <- '13' / '11' / '9' / '7' / '6' / '5'
  ChordModification <- ('b' / '#') ChordDegree / 'add' ChordDegree
  ChordDegree     <- '13' / '11' / '9' / '7' / '6' / '5' / '4' / '3' / '2' / '1'
  RootAccidental  <- 'b' / '#'
  UpperRoman      <- 'VII' / 'VI' / 'V' / 'IV' / 'III' / 'II' / 'I'
  LowerRoman      <- 'vii' / 'vi' / 'v' / 'iv' / 'iii' / 'ii' / 'i'
  NamedPitch      <- < NamedRoot >
  NamedRoot       <- [A-G] ('b' / '#')?
  ScaleDegree     <- < RootAccidental* SignedInteger >

  RegisterSuffix  <- < AbsoluteRegister RelativeRegister? / RelativeRegister >
  AbsoluteRegister <- '@' SignedInteger
  RelativeRegister <- Apostrophe+ / Comma+
  Apostrophe      <- "'"
  Comma           <- ','

  DurationSuffix  <- < Elongation Dots? / Dots >
  Elongation      <- '_' PositiveInteger / '_'+
  Dots            <- '..' / '.'
  EuclideanSuffix <- '(' H EuclideanPulses H ',' H EuclideanSteps
                     (H ',' H EuclideanRotation)? H ')'
  EuclideanPulses   <- < UnsignedInteger >
  EuclideanSteps    <- < PositiveInteger >
  EuclideanRotation <- < SignedInteger >
  RatchetSuffix   <- '*' RatchetCount
  RatchetCount    <- < PositiveInteger >
  ProbabilitySuffix <- '?' Probability?
  Probability     <- < Number >
  ReplicationSuffix <- '!' RepeatCount
  RepeatCount     <- < PositiveInteger >

  Attributes      <- '{' H Attribute (H ',' H Attribute)* H '}'
  Attribute       <- AttributeName H '=' H AttributeValue
  AttributeName   <- < Identifier >
  AttributeValue  <- < ScalarValue >

  ScalarPattern   <- RightAligned / LeftAligned / FreePattern
  RightAligned    <- Ellipsis S ScalarTerm (S ScalarTerm)*
  LeftAligned     <- ScalarTerm (S ScalarTerm)* S Ellipsis
  FreePattern     <- ScalarTerm (S ScalarTerm)*
  ScalarTerm      <- (RandomScalar / ScalarAtom / Default) ReplicationSuffix?
  RandomScalar    <- '$' ScalarRandomDistribution RandomArguments
  ScalarRandomDistribution <- < 'u' / 'n' >
  RandomArguments <- '{' H RandomArgument H ',' H RandomArgument H '}'
  RandomArgument  <- < ScalarValue >
  ScalarAtom      <- < ScalarValue >
  Ellipsis        <- '...'
  Default         <- < '.' >

  SettingValue    <- < ScalarValue / PitchValue / Identifier >
  ScalarValue     <- Number Unit?
  Unit            <- 'ms'
  Number          <- SignedInteger '/' PositiveInteger / Decimal
  Decimal         <- Sign? ([0-9]+ ('.' [0-9]+)? / '.' [0-9]+)
  SignedInteger   <- Sign? UnsignedInteger
  PositiveInteger <- [1-9] [0-9]*
  UnsignedInteger <- < [0-9]+ >
  Sign            <- '+' / '-'

  Expression <- Concatenation (H Pipeline)*
  Concatenation <- Term (H '+' H Term)*
  Term       <- Primary (H '*' H Integer)?
  Primary    <- Identifier / '(' H Expression H ')'
  Integer    <- < [1-9] [0-9]* >

  Pipeline      <- Pipe H (EveryTransform / SometimesTransform / TransformCall)
  EveryTransform <- 'every' S ConditionArgument S TransformCall
  SometimesTransform <- 'sometimes' S ConditionArgument S TransformCall
  TransformCall <- '(' H Operation (S TransformArgument)* H ')' /
                   Operation (S TransformArgument)*
  Operation     <- < [A-Za-z_] [A-Za-z0-9_]* >
  ConditionArgument <- < (!HChar !Pipe !CommentStart !Newline !')' .)+ >
  TransformArgument <- < (!HChar !Pipe !CommentStart !Newline !')' .)+ >
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

Ast FirstDescendantNamed(const Ast &node, const std::string &name) {
  if (!node)
    return nullptr;
  if (node->name == name)
    return node;
  for (const auto &child : node->nodes) {
    if (const auto found = FirstDescendantNamed(child, name))
      return found;
  }
  return nullptr;
}

PatternNode ReadPitchedValue(const Ast &pitched, Diagnostic &diagnostic) {
  auto pitchLeaf = [](const Ast &pitch) {
    PatternNode leaf;
    const auto named = FirstDescendantNamed(pitch, "NamedPitch");
    const auto degree = FirstDescendantNamed(pitch, "ScaleDegree");
    leaf.kind = named ? PatternKind::NamedPitch : PatternKind::ScaleDegree;
    leaf.atom = AstToken(named ? named : degree);
    if (const auto suffix = FirstDescendantNamed(pitch, "RegisterSuffix"))
      leaf.atom.text += AstToken(suffix).text;
    leaf.atom.span = Span(pitch);
    leaf.span = Span(pitch);
    return leaf;
  };
  PatternNode value;
  value.span = Span(pitched);
  if (const auto random = FirstChildNamed(pitched, "RandomPitch")) {
    value.kind = PatternKind::RandomPitch;
    value.span = Span(random);
    if (const auto distribution =
            FirstChildNamed(random, "PitchRandomDistribution"))
      value.atom = AstToken(distribution);
    if (const auto arguments = FirstChildNamed(random, "RandomArguments")) {
      for (const auto &argument : ChildrenNamed(arguments, "RandomArgument"))
        value.arguments.push_back(AstToken(argument));
    }
    return value;
  }
  if (const auto chord = FirstChildNamed(pitched, "ChordValue")) {
    if (const auto voicing = FirstChildNamed(chord, "ExplicitVoicing")) {
      value.kind = PatternKind::Voicing;
      value.span = Span(chord);
      for (const auto &tone : ChildrenNamed(voicing, "PitchValue")) {
        value.children.push_back(pitchLeaf(tone));
      }
      if (const auto suffix = FirstChildNamed(chord, "RegisterSuffix"))
        value.suffix = AstToken(suffix);
    } else {
      const auto roman = FirstChildNamed(chord, "RomanChord");
      value.kind = roman ? PatternKind::RomanChord : PatternKind::JazzChord;
      const auto symbol = roman ? roman : FirstChildNamed(chord, "JazzChord");
      value.atom = AstToken(symbol);
      if (const auto suffix = FirstChildNamed(chord, "RegisterSuffix"))
        value.atom.text += AstToken(suffix).text;
      value.atom.span = Span(chord);
      value.span = value.atom.span;
    }
    if (const auto slash = FirstChildNamed(pitched, "SlashSuffix")) {
      const auto bassPitch = FirstDescendantNamed(slash, "PitchValue");
      if (!bassPitch) {
        diagnostic = {"internal slash-chord AST is incomplete",
                      static_cast<int>(slash->line),
                      static_cast<int>(slash->column)};
        return {};
      }
      PatternNode bass;
      bass = pitchLeaf(bassPitch);
      PatternNode joined;
      joined.kind = PatternKind::Slash;
      joined.span = Span(pitched);
      joined.children.push_back(std::move(value));
      joined.children.push_back(std::move(bass));
      value = std::move(joined);
    }
    return value;
  }
  const auto pitch = pitched->name == "PitchValue"
                         ? pitched
                         : FirstChildNamed(pitched, "PitchValue");
  if (!pitch) {
    diagnostic = {"internal pitched-value AST is incomplete",
                  static_cast<int>(pitched->line),
                  static_cast<int>(pitched->column)};
    return {};
  }
  return pitchLeaf(pitch);
}

void ReadSuffixes(const Ast &source, PatternNode &node) {
  if (const auto duration = FirstChildNamed(source, "DurationSuffix"))
    node.durationSuffix = AstToken(duration);
  if (const auto euclidean = FirstChildNamed(source, "EuclideanSuffix"))
  {
    node.arguments.push_back(
        AstToken(FirstChildNamed(euclidean, "EuclideanPulses")));
    node.arguments.push_back(
        AstToken(FirstChildNamed(euclidean, "EuclideanSteps")));
    if (const auto rotation =
            FirstChildNamed(euclidean, "EuclideanRotation"))
      node.arguments.push_back(AstToken(rotation));
  }
  if (const auto ratchet = FirstChildNamed(source, "RatchetSuffix"))
    node.ratchetCount =
        AstToken(FirstDescendantNamed(ratchet, "RatchetCount"));
  if (const auto probability =
          FirstChildNamed(source, "ProbabilitySuffix")) {
    if (const auto amount = FirstChildNamed(probability, "Probability"))
      node.probability = AstToken(amount);
    else
      node.defaultProbability = true;
  }
  if (const auto replication =
          FirstChildNamed(source, "ReplicationSuffix"))
    node.repeatCount =
        AstToken(FirstDescendantNamed(replication, "RepeatCount"));
  if (const auto attributes = FirstChildNamed(source, "Attributes")) {
    for (const auto &attribute : ChildrenNamed(attributes, "Attribute")) {
      node.attributes.push_back(
          {AstToken(FirstChildNamed(attribute, "AttributeName")),
           AstToken(FirstChildNamed(attribute, "AttributeValue"))});
    }
  }
}

PatternNode ReadNoteElement(const Ast &element, Diagnostic &diagnostic);

std::vector<PatternNode> ReadNoteElements(const Ast &pattern,
                                          Diagnostic &diagnostic) {
  std::vector<PatternNode> result;
  for (const auto &element : ChildrenNamed(pattern, "NoteElement")) {
    result.push_back(ReadNoteElement(element, diagnostic));
    if (!diagnostic.message.empty())
      break;
  }
  return result;
}

PatternNode ReadNoteElement(const Ast &element, Diagnostic &diagnostic) {
  if (const auto event = FirstChildNamed(element, "NoteEvent")) {
    PatternNode node;
    node.kind = PatternKind::Event;
    node.span = Span(event);
    if (const auto onset = FirstChildNamed(event, "OnsetPrefix")) {
      if (const auto slide = FirstChildNamed(onset, "SlidePrefix"))
        node.slidePrefix = AstToken(slide);
      if (const auto dynamic = FirstChildNamed(onset, "DynamicPrefix"))
        node.dynamicPrefix = AstToken(dynamic);
    }
    const auto pitched = FirstChildNamed(event, "PitchedValue");
    node.children.push_back(ReadPitchedValue(pitched, diagnostic));
    ReadSuffixes(event, node);
    return node;
  }
  if (const auto rest = FirstChildNamed(element, "RestEvent")) {
    PatternNode node;
    node.kind = PatternKind::Rest;
    node.atom = AstToken(FirstDescendantNamed(rest, "RestMark"));
    node.span = Span(rest);
    ReadSuffixes(rest, node);
    return node;
  }
  if (const auto tie = FirstChildNamed(element, "TieExtension")) {
    PatternNode node;
    node.kind = PatternKind::Tie;
    node.atom = AstToken(FirstDescendantNamed(tie, "TieMark"));
    node.span = Span(tie);
    ReadSuffixes(tie, node);
    return node;
  }
  const auto group = FirstChildNamed(element, "GroupElement");
  if (!group) {
    diagnostic = {"internal note-pattern AST is incomplete",
                  static_cast<int>(element->line),
                  static_cast<int>(element->column)};
    return {};
  }
  PatternNode node;
  node.span = Span(group);
  const auto primary = FirstChildNamed(group, "GroupPrimary");
  if (const auto alternate = FirstChildNamed(primary, "Alternate")) {
    node.kind = PatternKind::CycleChoice;
    for (const auto &child : ChildrenNamed(alternate, "NoteElement"))
      node.children.push_back(ReadNoteElement(child, diagnostic));
  } else if (const auto bracket = FirstChildNamed(primary, "BracketGroup")) {
    if (const auto random = FirstChildNamed(bracket, "RandomChoice")) {
      node.kind = PatternKind::RandomChoice;
      for (const auto &branch : ChildrenNamed(random, "NotePattern")) {
        PatternNode branchNode;
        branchNode.kind = PatternKind::Subdivision;
        branchNode.span = Span(branch);
        branchNode.children = ReadNoteElements(branch, diagnostic);
        node.children.push_back(std::move(branchNode));
      }
    } else {
      node.kind = PatternKind::Subdivision;
      const auto pattern = FirstChildNamed(bracket, "NotePattern");
      node.children = ReadNoteElements(pattern, diagnostic);
    }
  }
  ReadSuffixes(group, node);
  return node;
}

Pattern ReadNotePattern(const Ast &node, Diagnostic &diagnostic) {
  Pattern pattern;
  pattern.span = Span(node);
  pattern.steps = ReadNoteElements(node, diagnostic);
  return pattern;
}

PatternNode ReadScalarTerm(const Ast &term) {
  PatternNode node;
  node.span = Span(term);
  if (const auto random = FirstChildNamed(term, "RandomScalar")) {
    node.kind = PatternKind::RandomScalar;
    node.atom =
        AstToken(FirstChildNamed(random, "ScalarRandomDistribution"));
    const auto arguments = FirstChildNamed(random, "RandomArguments");
    for (const auto &argument : ChildrenNamed(arguments, "RandomArgument"))
      node.arguments.push_back(AstToken(argument));
  } else if (const auto atom = FirstChildNamed(term, "ScalarAtom")) {
    node.kind = PatternKind::Atom;
    node.atom = AstToken(atom);
  } else {
    node.kind = PatternKind::Atom;
    node.atom = AstToken(FirstChildNamed(term, "Default"));
  }
  if (const auto replication = FirstChildNamed(term, "ReplicationSuffix"))
    node.repeatCount =
        AstToken(FirstDescendantNamed(replication, "RepeatCount"));
  return node;
}

Pattern ReadScalarPattern(const Ast &node) {
  Pattern pattern;
  pattern.span = Span(node);
  Ast body = FirstChildNamed(node, "FreePattern");
  if (const auto right = FirstChildNamed(node, "RightAligned")) {
    body = right;
    pattern.alignment = Pattern::Alignment::Right;
  } else if (const auto left = FirstChildNamed(node, "LeftAligned")) {
    body = left;
    pattern.alignment = Pattern::Alignment::Left;
  }
  for (const auto &term : ChildrenNamed(body, "ScalarTerm"))
    pattern.steps.push_back(ReadScalarTerm(term));
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
  if (node->name == "NotesLane") {
    lane.kind = Lane::Kind::Notes;
    lane.name = AstToken(FirstChildNamed(node, "NotesLaneName"));
    lane.pattern =
        ReadNotePattern(FirstChildNamed(node, "NotePattern"), diagnostic);
  } else if (node->name == "CvLane") {
    lane.kind = Lane::Kind::Cv;
    lane.name = AstToken(FirstChildNamed(node, "CvLaneName"));
    lane.pattern = ReadScalarPattern(FirstChildNamed(node, "ScalarPattern"));
  } else if (node->name == "SettingLine") {
    lane.kind = Lane::Kind::Setting;
    lane.name = AstToken(FirstChildNamed(node, "SettingName"));
    const auto value = FirstChildNamed(node, "SettingValue");
    lane.pattern.span = Span(value);
    PatternNode atom;
    atom.kind = PatternKind::Atom;
    atom.atom = AstToken(value);
    atom.span = atom.atom.span;
    lane.pattern.steps.push_back(std::move(atom));
  } else {
    lane.kind = Lane::Kind::Scalar;
    lane.name = AstToken(FirstChildNamed(node, "ScalarLaneName"));
    lane.pattern = ReadScalarPattern(FirstChildNamed(node, "ScalarPattern"));
  }
  for (const auto &pipeline : ChildrenNamed(node, "Pipeline"))
    lane.pipelines.push_back(ReadPipeline(pipeline));
  return lane;
}

SequenceDefinition ReadSequence(const Ast &node, Diagnostic &diagnostic) {
  SequenceDefinition sequence;
  sequence.span = Span(node);
  sequence.name = AstToken(FirstChildNamed(node, "Identifier"));
  for (const auto &child : node->nodes) {
    if (child->name == "NotesLane" || child->name == "ScalarLane" ||
        child->name == "CvLane" || child->name == "SettingLine")
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
  for (const auto &continuation :
       ChildrenNamed(node, "AssignmentContinuation")) {
    for (const auto &pipeline : ChildrenNamed(continuation, "Pipeline"))
      assignment.expression.pipelines.push_back(ReadPipeline(pipeline));
  }
  return assignment;
}

ParseResult ParseUsing(peg::parser &parser, const std::string &source) {
  ParseResult result;
  parser.set_logger(
      [&](std::size_t line, std::size_t column, const std::string &message) {
        if (result.diagnostic.message.empty())
          result.diagnostic = {message, static_cast<int>(line),
                               static_cast<int>(column)};
      });
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
          AstToken(FirstChildNamed(node, "UnsignedInteger")), Span(node)});
    }
  }
  return result;
}

bool LoadParser(peg::parser &parser, Diagnostic &diagnostic) {
  parser.set_logger(
      [&](std::size_t line, std::size_t column, const std::string &message) {
        if (diagnostic.message.empty())
          diagnostic = {"internal sequencer grammar: " + message,
                        static_cast<int>(line), static_cast<int>(column)};
      });
  if (!parser.load_grammar(Grammar)) {
    if (diagnostic.message.empty())
      diagnostic = {"internal sequencer grammar is invalid", 1, 1};
    return false;
  }
  parser.enable_ast();
  return true;
}

} // namespace

ParseResult Parse(const std::string &source) {
  ParseResult result;
  peg::parser parser;
  if (!LoadParser(parser, result.diagnostic))
    return result;
  return ParseUsing(parser, source);
}

ParseResult ParseStatementsContaining(const std::string &source,
                                       int selectionBegin,
                                       int selectionEnd) {
  ParseResult result;
  const int sourceSize = static_cast<int>(std::min<std::size_t>(
      source.size(), static_cast<std::size_t>(std::numeric_limits<int>::max())));
  selectionBegin = std::clamp(selectionBegin, 0, sourceSize);
  selectionEnd = std::clamp(selectionEnd, 0, sourceSize);
  if (selectionEnd < selectionBegin)
    std::swap(selectionBegin, selectionEnd);

  peg::parser parser;
  if (!LoadParser(parser, result.diagnostic))
    return result;

  std::vector<int> lineStarts{0};
  for (int index = 0; index < sourceSize; ++index) {
    if (source[static_cast<std::size_t>(index)] == '\n' &&
        index + 1 <= sourceSize)
      lineStarts.push_back(index + 1);
  }
  auto lineAt = [&](int position) {
    const auto found = std::upper_bound(lineStarts.begin(), lineStarts.end(),
                                        position);
    return static_cast<std::size_t>(
        std::max<std::ptrdiff_t>(0, found - lineStarts.begin() - 1));
  };
  const auto firstLine = lineAt(selectionBegin);
  const auto lastLine = lineAt(
      selectionEnd > selectionBegin ? selectionEnd - 1 : selectionBegin);

  std::size_t bestLines = std::numeric_limits<std::size_t>::max();
  for (std::size_t startLine = firstLine + 1; startLine-- > 0;) {
    for (std::size_t endLine = lastLine; endLine < lineStarts.size();
         ++endLine) {
      const int begin = lineStarts[startLine];
      const int end = endLine + 1 < lineStarts.size()
                          ? lineStarts[endLine + 1]
                          : sourceSize;
      if (end < selectionEnd)
        continue;

      // A continuation beginning with |> belongs to the sequence above it;
      // do not accept a shorter otherwise-valid prefix that silently drops it.
      if (endLine + 1 < lineStarts.size()) {
        int cursor = lineStarts[endLine + 1];
        while (cursor < sourceSize &&
               (source[static_cast<std::size_t>(cursor)] == ' ' ||
                source[static_cast<std::size_t>(cursor)] == '\t'))
          ++cursor;
        if (cursor + 1 < sourceSize &&
            source[static_cast<std::size_t>(cursor)] == '|' &&
            source[static_cast<std::size_t>(cursor + 1)] == '>')
          continue;
      }

      std::string candidate = source.substr(0, static_cast<std::size_t>(end));
      for (int index = 0; index < begin; ++index) {
        auto &character = candidate[static_cast<std::size_t>(index)];
        if (character != '\n' && character != '\r')
          character = ' ';
      }
      auto parsed = ParseUsing(parser, candidate);
      if (!parsed || parsed.document.statements.empty())
        continue;
      const bool overlaps = std::any_of(
          parsed.document.statements.begin(), parsed.document.statements.end(),
          [&](const Statement &statement) {
            const auto span =
                std::visit([](const auto &value) { return value.span; },
                           statement);
            return selectionBegin == selectionEnd
                       ? span.begin <= selectionBegin && selectionBegin < span.end
                       : span.begin < selectionEnd && selectionBegin < span.end;
          });
      if (!overlaps)
        continue;
      const auto lineCount = endLine - startLine + 1;
      if (lineCount < bestLines) {
        bestLines = lineCount;
        result = std::move(parsed);
      }
    }
  }
  if (bestLines != std::numeric_limits<std::size_t>::max())
    return result;

  int line = 1;
  int column = 1;
  for (int index = 0; index < selectionBegin; ++index) {
    if (source[static_cast<std::size_t>(index)] == '\n') {
      ++line;
      column = 1;
    } else {
      ++column;
    }
  }
  result.diagnostic =
      {"selected text is not a complete valid statement", line, column};
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
  Invalidate(node.slidePrefix);
  Invalidate(node.dynamicPrefix);
  Invalidate(node.durationSuffix);
  Invalidate(node.ratchetCount);
  Invalidate(node.probability);
  for (auto &attribute : node.attributes) {
    Invalidate(attribute.name);
    Invalidate(attribute.value);
  }
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
