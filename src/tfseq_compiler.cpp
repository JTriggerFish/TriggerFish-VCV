#include "tfseq.hpp"
#include "tfseq_parser.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <new>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace tfseq {

struct CompiledProgramFactory {
  static std::unique_ptr<CompiledProgram> create(SemanticProgram &&semantic) {
    return std::unique_ptr<CompiledProgram>(
        new CompiledProgram(std::move(semantic)));
  }
};

namespace {

using Token = syntax::Token;

Diagnostic Error(const SourceSpan &span, const std::string &message) {
  return {message, span.line, span.column};
}

std::uint64_t StableDefinitionId(const std::string &name) noexcept {
  std::uint64_t hash = 1469598103934665603ULL;
  for (const unsigned char character : name) {
    hash ^= character;
    hash *= 1099511628211ULL;
  }
  return hash == 0 ? 1 : hash;
}

bool ParseNumber(const std::string &text, double &value) {
  const auto slash = text.find('/');
  if (slash != std::string::npos) {
    double numerator = 0.0;
    double denominator = 0.0;
    if (!ParseNumber(text.substr(0, slash), numerator) ||
        !ParseNumber(text.substr(slash + 1), denominator) || denominator == 0.0)
      return false;
    value = numerator / denominator;
    return std::isfinite(value);
  }
  char *end = nullptr;
  value = std::strtod(text.c_str(), &end);
  return end == text.c_str() + text.size() && std::isfinite(value);
}

int NoteSemitone(char letter) {
  switch (static_cast<char>(std::toupper(static_cast<unsigned char>(letter)))) {
  case 'C':
    return 0;
  case 'D':
    return 2;
  case 'E':
    return 4;
  case 'F':
    return 5;
  case 'G':
    return 7;
  case 'A':
    return 9;
  case 'B':
    return 11;
  default:
    return -100;
  }
}

int SaturatingIntAdd(int left, int right) {
  const auto result = static_cast<long long>(left) + right;
  return static_cast<int>(
      std::clamp<long long>(result, std::numeric_limits<int>::min() + 1LL,
                            std::numeric_limits<int>::max()));
}

bool ParseIntText(const std::string &text, int &value) {
  if (text.empty())
    return false;
  errno = 0;
  char *end = nullptr;
  const auto parsed = std::strtoll(text.c_str(), &end, 10);
  if (errno == ERANGE || end != text.c_str() + text.size() ||
      parsed < std::numeric_limits<int>::min() + 1LL ||
      parsed > std::numeric_limits<int>::max())
    return false;
  value = static_cast<int>(parsed);
  return true;
}

bool ParseRegisterSuffix(std::string &text, bool &hasOctave, int &octave,
                         int &octaveOffset) {
  long long relativeOffset = 0;
  while (!text.empty() && (text.back() == '\'' || text.back() == ',')) {
    relativeOffset += text.back() == '\'' ? 1 : -1;
    if (relativeOffset < std::numeric_limits<int>::min() + 1LL ||
        relativeOffset > std::numeric_limits<int>::max())
      return false;
    text.pop_back();
  }
  octaveOffset = static_cast<int>(relativeOffset);

  const auto marker = text.rfind('@');
  if (marker == std::string::npos)
    return text.find('@') == std::string::npos;
  const std::string registerText = text.substr(marker + 1);
  int parsed = 0;
  if (!ParseIntText(registerText, parsed))
    return false;
  hasOctave = true;
  octave = parsed;
  text.resize(marker);
  return text.find('@') == std::string::npos;
}

bool ParseNoteName(const std::string &text, std::size_t &cursor,
                   int &pitchClass) {
  if (cursor >= text.size() ||
      !std::isupper(static_cast<unsigned char>(text[cursor])))
    return false;
  int semitone = NoteSemitone(text[cursor++]);
  if (semitone < 0)
    return false;
  if (cursor < text.size() && (text[cursor] == '#' || text[cursor] == 'b')) {
    semitone += text[cursor++] == '#' ? 1 : -1;
  }
  pitchClass = ((semitone % 12) + 12) % 12;
  return true;
}

bool ParseAbsoluteNote(const std::string &source, PitchValue &value) {
  std::string text = source;
  PitchValue parsed;
  if (!ParseRegisterSuffix(text, parsed.hasOctave, parsed.octave,
                           parsed.octaveOffset) ||
      !parsed.hasOctave)
    return false;
  std::size_t cursor = 0;
  if (!ParseNoteName(text, cursor, parsed.pitchClass) || cursor != text.size())
    return false;
  parsed.absolute = true;
  value = parsed;
  return true;
}

bool ParseBareNamedPitch(const Token &token, PitchValue &value) {
  std::string text = token.text;
  bool hasOctave = false;
  int octave = 0;
  int octaveOffset = 0;
  if (!ParseRegisterSuffix(text, hasOctave, octave, octaveOffset))
    return false;
  std::size_t cursor = 0;
  int pitchClass = 0;
  if (!ParseNoteName(text, cursor, pitchClass) || cursor != text.size())
    return false;
  value = {};
  value.absolute = true;
  value.pitchClass = pitchClass;
  value.hasOctave = hasOctave;
  value.octave = octave;
  value.octaveOffset = octaveOffset;
  value.span = token.span;
  return true;
}

bool ParseTonic(const std::string &source, Scale &scale) {
  std::string text = source;
  bool hasOctave = false;
  int octave = 0;
  int octaveOffset = 0;
  if (!ParseRegisterSuffix(text, hasOctave, octave, octaveOffset) ||
      octaveOffset != 0)
    return false;
  std::size_t cursor = 0;
  if (!ParseNoteName(text, cursor, scale.tonicSemitone) ||
      cursor != text.size())
    return false;
  if (hasOctave)
    scale.tonicOctave = octave;
  return true;
}

bool ParsePitchValue(const Token &token, PitchValue &value,
                     Diagnostic &diagnostic) {
  Token pitchToken = token;
  if (ParseAbsoluteNote(pitchToken.text, value)) {
    value.span = token.span;
    return true;
  }

  if (!ParseRegisterSuffix(pitchToken.text, value.hasOctave, value.octave,
                           value.octaveOffset)) {
    diagnostic =
        Error(token.span, "invalid octave suffix in '" + token.text + "'");
    return false;
  }

  std::size_t cursor = 0;
  int accidental = 0;
  while (cursor < pitchToken.text.size() &&
         (pitchToken.text[cursor] == 'b' || pitchToken.text[cursor] == '#')) {
    accidental += pitchToken.text[cursor] == '#' ? 1 : -1;
    ++cursor;
  }
  if (cursor >= pitchToken.text.size()) {
    diagnostic = Error(token.span, "expected a scale degree");
    return false;
  }
  char *end = nullptr;
  errno = 0;
  const long long degree =
      std::strtoll(pitchToken.text.c_str() + cursor, &end, 10);
  if (errno == ERANGE ||
      end != pitchToken.text.c_str() + pitchToken.text.size() ||
      degree < std::numeric_limits<int>::min() + 1LL ||
      degree > std::numeric_limits<int>::max()) {
    diagnostic = Error(token.span, "invalid scale degree '" + token.text + "'");
    return false;
  }
  value.degree = static_cast<int>(degree);
  value.accidental = accidental;
  value.span = token.span;
  return true;
}

void AddChordInterval(std::vector<int> &intervals, int semitones) {
  if (std::find(intervals.begin(), intervals.end(), semitones) ==
      intervals.end())
    intervals.push_back(semitones);
}

int NaturalChordInterval(int degree) {
  switch (degree) {
  case 1:
    return 0;
  case 2:
  case 9:
    return 14;
  case 3:
    return 4;
  case 4:
  case 11:
    return 17;
  case 5:
    return 7;
  case 6:
  case 13:
    return 21;
  case 7:
    return 11;
  default:
    return std::numeric_limits<int>::min();
  }
}

PitchValue NamedChordTone(int rootPitchClass, int interval, bool hasOctave,
                          int octave, int octaveOffset,
                          const SourceSpan &span) {
  PitchValue pitch;
  pitch.absolute = true;
  const int absoluteClass = rootPitchClass + interval;
  pitch.pitchClass = ((absoluteClass % 12) + 12) % 12;
  const int registerCarry =
      absoluteClass >= 0 ? absoluteClass / 12 : -((-absoluteClass + 11) / 12);
  pitch.hasOctave = hasOctave;
  pitch.octave = SaturatingIntAdd(octave, registerCarry);
  pitch.octaveOffset =
      hasOctave ? octaveOffset : SaturatingIntAdd(octaveOffset, registerCarry);
  pitch.span = span;
  return pitch;
}

bool ParseSlashPitch(const Token &token, PitchValue &pitch,
                     Diagnostic &diagnostic) {
  Diagnostic ignored;
  if (ParsePitchValue(token, pitch, ignored))
    return true;
  std::string text = token.text;
  if (!ParseRegisterSuffix(text, pitch.hasOctave, pitch.octave,
                           pitch.octaveOffset)) {
    diagnostic = Error(token.span, "invalid slash-bass register");
    return false;
  }
  std::size_t cursor = 0;
  if (!ParseNoteName(text, cursor, pitch.pitchClass) || cursor != text.size()) {
    diagnostic =
        Error(token.span, "invalid slash-bass note '" + token.text + "'");
    return false;
  }
  pitch.absolute = true;
  pitch.span = token.span;
  return true;
}

bool ParseJazzChord(const Token &token, ChordValue &chord,
                    Diagnostic &diagnostic) {
  std::string main = token.text;

  bool hasOctave = false;
  int octave = 0;
  int octaveOffset = 0;
  if (!ParseRegisterSuffix(main, hasOctave, octave, octaveOffset)) {
    diagnostic =
        Error(token.span, "invalid chord register in '" + token.text + "'");
    return false;
  }
  std::size_t cursor = 0;
  int root = 0;
  if (!ParseNoteName(main, cursor, root))
    return false;

  enum class Triad { Major, Minor, Diminished, Augmented, Sus2, Sus4 };
  Triad triad = Triad::Major;
  bool majorSeventh = false;
  if (main.compare(cursor, 3, "maj") == 0) {
    majorSeventh = true;
    cursor += 3;
  } else if (main.compare(cursor, 3, "min") == 0) {
    triad = Triad::Minor;
    cursor += 3;
  } else if (main.compare(cursor, 3, "dim") == 0) {
    triad = Triad::Diminished;
    cursor += 3;
  } else if (main.compare(cursor, 3, "aug") == 0) {
    triad = Triad::Augmented;
    cursor += 3;
  } else if (main.compare(cursor, 4, "sus2") == 0) {
    triad = Triad::Sus2;
    cursor += 4;
  } else if (main.compare(cursor, 4, "sus4") == 0) {
    triad = Triad::Sus4;
    cursor += 4;
  } else if (cursor < main.size() && main[cursor] == 'm') {
    triad = Triad::Minor;
    ++cursor;
  }

  int extension = 0;
  const auto extensionBegin = cursor;
  while (cursor < main.size() &&
         std::isdigit(static_cast<unsigned char>(main[cursor])))
    ++cursor;
  if (cursor > extensionBegin &&
      !ParseIntText(main.substr(extensionBegin, cursor - extensionBegin),
                    extension)) {
    diagnostic = Error(token.span, "invalid chord extension");
    return false;
  }
  if (extension != 0 && extension != 5 && extension != 6 && extension != 7 &&
      extension != 9 && extension != 11 && extension != 13) {
    diagnostic = Error(token.span,
                       "unsupported chord extension in '" + token.text + "'");
    return false;
  }

  struct Alteration {
    int degree = 0;
    int offset = 0;
    bool add = false;
  };
  std::vector<Alteration> alterations;
  while (cursor < main.size()) {
    Alteration alteration;
    if (main.compare(cursor, 3, "add") == 0) {
      alteration.add = true;
      cursor += 3;
    } else if (main[cursor] == 'b' || main[cursor] == '#') {
      alteration.offset = main[cursor++] == 'b' ? -1 : 1;
    } else {
      diagnostic =
          Error(token.span, "invalid chord quality in '" + token.text + "'");
      return false;
    }
    const auto degreeBegin = cursor;
    while (cursor < main.size() &&
           std::isdigit(static_cast<unsigned char>(main[cursor])))
      ++cursor;
    if (degreeBegin == cursor ||
        !ParseIntText(main.substr(degreeBegin, cursor - degreeBegin),
                      alteration.degree) ||
        NaturalChordInterval(alteration.degree) ==
            std::numeric_limits<int>::min()) {
      diagnostic =
          Error(token.span, "invalid chord alteration in '" + token.text + "'");
      return false;
    }
    alterations.push_back(alteration);
  }

  std::vector<int> intervals;
  switch (triad) {
  case Triad::Major:
    intervals = {0, 4, 7};
    break;
  case Triad::Minor:
    intervals = {0, 3, 7};
    break;
  case Triad::Diminished:
    intervals = {0, 3, 6};
    break;
  case Triad::Augmented:
    intervals = {0, 4, 8};
    break;
  case Triad::Sus2:
    intervals = {0, 2, 7};
    break;
  case Triad::Sus4:
    intervals = {0, 5, 7};
    break;
  }
  if (extension == 5)
    intervals = {0, 7};
  if (extension == 6)
    AddChordInterval(intervals, 9);
  if (extension >= 7) {
    const int seventh = triad == Triad::Diminished ? 9 : majorSeventh ? 11 : 10;
    AddChordInterval(intervals, seventh);
  }
  if (extension >= 9)
    AddChordInterval(intervals, 14);
  if (extension >= 11)
    AddChordInterval(intervals, 17);
  if (extension >= 13)
    AddChordInterval(intervals, 21);

  for (const auto &alteration : alterations) {
    const int natural = NaturalChordInterval(alteration.degree);
    auto found =
        std::find_if(intervals.begin(), intervals.end(), [&](int tone) {
          return ((tone - natural) % 12 + 12) % 12 == 0;
        });
    if (found != intervals.end() && !alteration.add)
      *found += alteration.offset;
    else
      AddChordInterval(intervals, natural + alteration.offset);
  }
  std::sort(intervals.begin(), intervals.end());
  intervals.erase(std::unique(intervals.begin(), intervals.end()),
                  intervals.end());
  if (intervals.size() > MaximumPolyphony) {
    diagnostic = Error(token.span, "chord exceeds Rack's 16-channel polyphony");
    return false;
  }

  chord.meaning = ChordValue::Meaning::JazzSymbol;
  chord.jazzSymbol = token.text;
  chord.rootPitchClass = root;
  chord.span = token.span;
  for (const int interval : intervals)
    chord.voices.push_back(NamedChordTone(root, interval, hasOctave, octave,
                                          octaveOffset, token.span));
  return true;
}

void ApplyChordRegister(ChordValue &chord, bool hasOctave, int octave,
                        int octaveOffset) {
  for (auto &voice : chord.voices) {
    if (hasOctave && !voice.hasOctave) {
      voice.hasOctave = true;
      voice.octave = octave;
      voice.octaveOffset = SaturatingIntAdd(voice.octaveOffset, octaveOffset);
    } else if (!hasOctave) {
      voice.octaveOffset = SaturatingIntAdd(voice.octaveOffset, octaveOffset);
    }
  }
}

bool ParseExplicitChord(const syntax::PatternNode &node, ChordValue &chord,
                        Diagnostic &diagnostic) {
  if (node.kind != syntax::PatternKind::Voicing)
    return false;
  bool hasOctave = false;
  int octave = 0;
  int octaveOffset = 0;
  std::string suffix = node.suffix.text;
  if (!ParseRegisterSuffix(suffix, hasOctave, octave, octaveOffset) ||
      !suffix.empty()) {
    diagnostic = Error(node.suffix.span.valid() ? node.suffix.span : node.span,
                       "invalid chord octave suffix");
    return false;
  }
  if (node.children.empty()) {
    diagnostic = Error(node.span, "chord cannot be empty");
    return false;
  }
  if (node.children.size() == 1) {
    const auto &tone = node.children.front().atom;
    PitchValue namedPitch;
    if (ParseBareNamedPitch(tone, namedPitch)) {
      chord.voices.push_back(namedPitch);
      chord.meaning = ChordValue::Meaning::ExplicitVoicing;
      ApplyChordRegister(chord, hasOctave, octave, octaveOffset);
      chord.span = node.span;
      return true;
    }
    Diagnostic jazzDiagnostic;
    if (ParseJazzChord(tone, chord, jazzDiagnostic)) {
      ApplyChordRegister(chord, hasOctave, octave, octaveOffset);
      chord.meaning = ChordValue::Meaning::ExplicitVoicing;
      chord.span = node.span;
      return true;
    }
  }

  for (const auto &child : node.children) {
    const Token &tone = child.atom;
    PitchValue pitch;
    Diagnostic pitchDiagnostic;
    if (!ParsePitchValue(tone, pitch, pitchDiagnostic) &&
        !ParseBareNamedPitch(tone, pitch)) {
      diagnostic = pitchDiagnostic;
      return false;
    }
    chord.voices.push_back(pitch);
  }
  if (chord.voices.empty()) {
    diagnostic = Error(node.span, "chord requires at least one pitched voice");
    return false;
  }
  if (chord.voices.size() > MaximumPolyphony) {
    diagnostic = Error(node.span, "chord exceeds Rack's 16-channel polyphony");
    return false;
  }
  chord.meaning = ChordValue::Meaning::ExplicitVoicing;
  ApplyChordRegister(chord, hasOctave, octave, octaveOffset);
  chord.span = node.span;
  return true;
}

bool ParsePitchedChoice(const Token &token, ChordValue &chord,
                        Diagnostic &diagnostic) {
  PitchValue pitch;
  Diagnostic pitchDiagnostic;
  if (ParsePitchValue(token, pitch, pitchDiagnostic)) {
    chord.voices.push_back(pitch);
    chord.span = token.span;
    return true;
  }
  Diagnostic chordDiagnostic;
  if (ParseJazzChord(token, chord, chordDiagnostic))
    return true;
  diagnostic =
      chordDiagnostic.message.empty() ? pitchDiagnostic : chordDiagnostic;
  return false;
}

std::string PatternNodeText(const syntax::PatternNode &node);

bool ParsePitchedChoice(const syntax::PatternNode &node, ChordValue &chord,
                        Diagnostic &diagnostic) {
  if (node.kind == syntax::PatternKind::Atom)
    return ParsePitchedChoice(node.atom, chord, diagnostic);
  if (node.kind == syntax::PatternKind::Voicing)
    return ParseExplicitChord(node, chord, diagnostic);
  if (node.kind != syntax::PatternKind::Slash || node.children.size() != 2 ||
      node.children[1].kind != syntax::PatternKind::Atom) {
    diagnostic = Error(node.span, "invalid pitched value");
    return false;
  }
  if (!ParsePitchedChoice(node.children[0], chord, diagnostic))
    return false;
  if (chord.hasBass) {
    diagnostic = Error(node.span, "a chord can contain only one slash bass");
    return false;
  }
  if (!ParseSlashPitch(node.children[1].atom, chord.bass, diagnostic))
    return false;
  chord.hasBass = true;
  chord.span = node.span;
  if (chord.meaning == ChordValue::Meaning::JazzSymbol)
    chord.jazzSymbol = PatternNodeText(node);
  if (chord.voices.size() + 1 > MaximumPolyphony) {
    diagnostic = Error(node.span, "chord exceeds Rack's 16-channel polyphony");
    return false;
  }
  return true;
}

bool PatternRepeat(const syntax::PatternNode &source, const std::string &label,
                   const syntax::PatternNode *&node, std::size_t &repetitions,
                   Diagnostic &diagnostic) {
  node = &source;
  repetitions = 1;
  if (source.kind != syntax::PatternKind::Repeat)
    return true;
  if (source.children.size() != 1) {
    diagnostic = Error(source.span, "invalid repeat pattern");
    return false;
  }
  errno = 0;
  char *end = nullptr;
  const auto parsed = std::strtoull(source.repeatCount.text.c_str(), &end, 10);
  if (source.repeatCount.text.empty() ||
      source.repeatCount.text.front() == '-' || errno == ERANGE ||
      parsed == 0 ||
      end != source.repeatCount.text.c_str() + source.repeatCount.text.size() ||
      parsed > std::numeric_limits<std::size_t>::max()) {
    diagnostic = Error(source.repeatCount.span,
                       label + " repetition must be a positive integer that "
                               "fits addressable memory");
    return false;
  }
  repetitions = static_cast<std::size_t>(parsed);
  node = &source.children.front();
  return true;
}

bool ParsePitchPatternNode(const syntax::PatternNode &node, PitchItem &item,
                           Diagnostic &diagnostic) {
  item.span = node.span;
  const std::vector<syntax::PatternNode> *choices = nullptr;
  if (node.kind == syntax::PatternKind::CycleChoice) {
    item.choice = PitchItem::Choice::Alternate;
    choices = &node.children;
  } else if (node.kind == syntax::PatternKind::RandomChoice) {
    item.choice = PitchItem::Choice::Random;
    choices = &node.children;
  } else if (node.kind == syntax::PatternKind::Subdivision) {
    diagnostic = Error(
        node.span,
        "note subdivisions are reserved but not executable in this version");
    return false;
  } else if (node.kind != syntax::PatternKind::Atom &&
             node.kind != syntax::PatternKind::Voicing &&
             node.kind != syntax::PatternKind::Slash) {
    diagnostic = Error(node.span, "invalid nested note pattern");
    return false;
  }

  if (!choices) {
    ChordValue value;
    if (!ParsePitchedChoice(node, value, diagnostic))
      return false;
    item.values.push_back(std::move(value));
    return true;
  }
  for (const auto &choice : *choices) {
    if (choice.kind != syntax::PatternKind::Atom &&
        choice.kind != syntax::PatternKind::Voicing &&
        choice.kind != syntax::PatternKind::Slash) {
      diagnostic =
          Error(choice.span, "nested note choices are reserved for a future "
                             "pattern evaluator");
      return false;
    }
    ChordValue value;
    if (!ParsePitchedChoice(choice, value, diagnostic))
      return false;
    item.values.push_back(std::move(value));
  }
  return true;
}

bool ParseNotes(const syntax::Pattern &pattern, Sequence &sequence,
                Diagnostic &diagnostic) {
  for (const auto &source : pattern.steps) {
    const syntax::PatternNode *musicalNode = nullptr;
    std::size_t repetitions = 1;
    if (!PatternRepeat(source, "note", musicalNode, repetitions, diagnostic))
      return false;
    PitchItem item;
    if (!ParsePitchPatternNode(*musicalNode, item, diagnostic))
      return false;
    item.span = source.span;
    for (std::size_t repetition = 0; repetition < repetitions; ++repetition)
      sequence.notes.push_back(item);
  }
  return true;
}

ArticulationAtom ParseArticulationAtom(const Token &token,
                                       Diagnostic &diagnostic) {
  ArticulationAtom atom;
  atom.span = token.span;
  if (token.text.empty()) {
    diagnostic = Error(token.span, "empty articulation");
    return atom;
  }

  const char base = token.text.front();
  bool hasProbability = false;
  bool hasRatchet = false;
  std::size_t cursor = 1;
  while (cursor < token.text.size()) {
    const char modifier = token.text[cursor++];
    if (modifier != '?' && modifier != '*') {
      diagnostic =
          Error(token.span, "unknown articulation '" + token.text + "'");
      return atom;
    }
    const auto argumentBegin = cursor;
    while (cursor < token.text.size() && token.text[cursor] != '?' &&
           token.text[cursor] != '*')
      ++cursor;
    const auto argument =
        token.text.substr(argumentBegin, cursor - argumentBegin);
    if (modifier == '?') {
      if (hasProbability) {
        diagnostic =
            Error(token.span, "articulation probability is specified twice");
        return atom;
      }
      hasProbability = true;
      atom.probability = 0.5f;
      if (!argument.empty()) {
        double parsed = 0.0;
        if (!ParseNumber(argument, parsed) || parsed < 0.0 || parsed > 1.0) {
          diagnostic = Error(token.span, "invalid articulation probability");
          return atom;
        }
        atom.probability = static_cast<float>(parsed);
      }
    } else {
      if (hasRatchet) {
        diagnostic =
            Error(token.span, "articulation ratchet is specified twice");
        return atom;
      }
      hasRatchet = true;
      double parsed = 0.0;
      if (!ParseNumber(argument, parsed) || parsed < 1.0 ||
          static_cast<long double>(parsed) >
              std::numeric_limits<std::size_t>::max() ||
          std::floor(parsed) != parsed) {
        diagnostic = Error(
            token.span, "ratchet count must be a positive addressable integer");
        return atom;
      }
      atom.ratchets = static_cast<std::size_t>(parsed);
    }
  }

  if (base == 'x' || base == '.')
    atom.kind = ArticulationKind::Attack;
  else if (base == '>')
    atom.kind = ArticulationKind::Slide;
  else if (base == '_')
    atom.kind = ArticulationKind::Tie;
  else if (base == '~')
    atom.kind = ArticulationKind::Rest;
  else
    diagnostic = Error(token.span, "unknown articulation '" + token.text + "'");
  return atom;
}

void AppendEuclidean(const Token &token, int pulses, int steps, int rotation,
                     Sequence &sequence, Diagnostic &diagnostic) {
  if (pulses < 0 || pulses > steps || steps < 1) {
    diagnostic =
        Error(token.span, "Euclidean rhythm requires 0 <= pulses <= steps");
    return;
  }
  for (int step = 0; step < steps; ++step) {
    const auto rotated =
        ((static_cast<std::int64_t>(step) + rotation) % steps + steps) % steps;
    const bool hit = (rotated * pulses) % steps < pulses;
    ArticulationStep item;
    item.span = token.span;
    item.atoms.push_back(
        {hit ? ArticulationKind::Attack : ArticulationKind::Rest, 1, 1.f,
         token.span});
    sequence.articulation.push_back(std::move(item));
  }
}

bool AppendArticulationAtoms(const syntax::PatternNode &source,
                             std::vector<ArticulationAtom> &atoms,
                             Diagnostic &diagnostic) {
  const syntax::PatternNode *node = nullptr;
  std::size_t repetitions = 1;
  if (!PatternRepeat(source, "articulation", node, repetitions, diagnostic))
    return false;
  if (node->kind != syntax::PatternKind::Atom) {
    diagnostic = Error(node->span, "nested articulation groups are not valid");
    return false;
  }
  const auto atom = ParseArticulationAtom(node->atom, diagnostic);
  if (!diagnostic.message.empty())
    return false;
  for (std::size_t repetition = 0; repetition < repetitions; ++repetition)
    atoms.push_back(atom);
  return true;
}

bool ParseEuclidean(const syntax::PatternNode &node, int &pulses, int &steps,
                    int &rotation, Diagnostic &diagnostic) {
  if (node.kind != syntax::PatternKind::Euclidean)
    return false;
  if (node.arguments.size() < 2 || node.arguments.size() > 3) {
    diagnostic = Error(node.span, "expected x(pulses,steps[,rotation])");
    return true;
  }
  auto parse = [](const Token &argument, int &value) {
    double number = 0.0;
    if (!ParseNumber(argument.text, number) || std::floor(number) != number ||
        number < std::numeric_limits<int>::min() ||
        number > std::numeric_limits<int>::max())
      return false;
    value = static_cast<int>(number);
    return true;
  };
  if (!parse(node.arguments[0], pulses) || !parse(node.arguments[1], steps) ||
      (node.arguments.size() == 3 && !parse(node.arguments[2], rotation)))
    diagnostic = Error(node.span, "expected x(pulses,steps[,rotation])");
  return true;
}

bool ParseArticulation(const syntax::Pattern &pattern, Sequence &sequence,
                       Diagnostic &diagnostic) {
  for (const auto &source : pattern.steps) {
    const syntax::PatternNode *node = nullptr;
    std::size_t repetitions = 1;
    if (!PatternRepeat(source, "articulation", node, repetitions, diagnostic))
      return false;
    if (node->kind == syntax::PatternKind::CycleChoice ||
        node->kind == syntax::PatternKind::RandomChoice) {
      diagnostic = Error(node->span,
                         "articulation choices are reserved but not executable "
                         "in this version");
      return false;
    }
    if (node->kind == syntax::PatternKind::Euclidean) {
      int pulses = 0;
      int steps = 0;
      int rotation = 0;
      const bool euclidean =
          ParseEuclidean(*node, pulses, steps, rotation, diagnostic);
      if (!diagnostic.message.empty())
        return false;
      if (euclidean) {
        const Token site{node->atom.text, node->span};
        for (std::size_t repetition = 0; repetition < repetitions; ++repetition)
          AppendEuclidean(site, pulses, steps, rotation, sequence, diagnostic);
        if (!diagnostic.message.empty())
          return false;
        continue;
      }
    }

    ArticulationStep step;
    step.span = source.span;
    if (node->kind == syntax::PatternKind::Subdivision) {
      for (const auto &child : node->children) {
        if (!AppendArticulationAtoms(child, step.atoms, diagnostic))
          return false;
      }
    } else if (!AppendArticulationAtoms(*node, step.atoms, diagnostic)) {
      return false;
    }
    for (std::size_t repetition = 0; repetition < repetitions; ++repetition)
      sequence.articulation.push_back(step);
  }
  return true;
}

bool ParseScalars(const syntax::Pattern &pattern,
                  std::vector<ScalarItem> &items, const std::string &lane,
                  Diagnostic &diagnostic) {
  for (const auto &source : pattern.steps) {
    const syntax::PatternNode *node = nullptr;
    std::size_t repetitions = 1;
    if (!PatternRepeat(source, lane, node, repetitions, diagnostic))
      return false;
    if (node->kind != syntax::PatternKind::Atom &&
        node->kind != syntax::PatternKind::Slash) {
      diagnostic =
          Error(node->span, lane + " grouped patterns are reserved but not "
                                   "executable in this version");
      return false;
    }
    Token reconstructed;
    const Token *scalarToken = &node->atom;
    if (node->kind == syntax::PatternKind::Slash) {
      reconstructed = {PatternNodeText(*node), node->span};
      scalarToken = &reconstructed;
    }
    ScalarItem item;
    item.span = source.span;
    if (scalarToken->text == ".") {
      item.isDefault = true;
    } else if (lane == "accent" && scalarToken->text == "+") {
      item.value = 0.88f;
    } else if (lane == "accent" && scalarToken->text == "++") {
      item.value = 1.f;
    } else if (lane == "slide" && scalarToken->text == ">") {
      item.value = 0.25f;
    } else {
      std::string number = scalarToken->text;
      if (!number.empty() && number.front() == '.')
        number.insert(number.begin(), '0');
      double parsed = 0.0;
      if (!ParseNumber(number, parsed)) {
        diagnostic = Error(source.span, "invalid " + lane + " value '" +
                                            scalarToken->text + "'");
        return false;
      }
      item.value = parsed;
    }
    if (!item.isDefault && !std::isfinite(item.value)) {
      diagnostic = Error(
          source.span, lane + " value is outside the supported numeric range");
      return false;
    }
    if (!item.isDefault &&
        (lane == "velocity" || lane == "vel" || lane == "accent" ||
         lane == "gate") &&
        (item.value < 0.f || item.value > 1.f)) {
      diagnostic = Error(source.span, lane + " must be from 0 to 1");
      return false;
    }
    if (!item.isDefault && (lane == "duration" || lane == "dur") &&
        item.value <= 0.f) {
      diagnostic = Error(source.span, "duration must be positive");
      return false;
    }
    if (!item.isDefault && lane == "slide" && item.value < 0.f) {
      diagnostic = Error(source.span, "slide must be non-negative");
      return false;
    }
    if (!item.isDefault && lane == "ratchet" &&
        (item.value < 1.0 ||
         static_cast<long double>(item.value) >
             std::numeric_limits<std::size_t>::max() ||
         std::floor(item.value) != item.value)) {
      diagnostic =
          Error(source.span, "ratchet must be a positive addressable integer");
      return false;
    }
    if (!item.isDefault && lane == "octave" &&
        (std::floor(item.value) != item.value ||
         item.value < std::numeric_limits<int>::min() ||
         item.value > std::numeric_limits<int>::max())) {
      diagnostic =
          Error(source.span, "octave must fit the supported integer range");
      return false;
    }
    for (std::size_t repetition = 0; repetition < repetitions; ++repetition)
      items.push_back(item);
  }
  return true;
}

enum class LaneValueKind {
  Notes,
  Articulation,
  Scalar,
  Cycle,
  Tonic,
  Scale,
  Glide
};

struct LaneSpec {
  const char *spelling;
  const char *canonical;
  LaneValueKind valueKind;
  CursorLane cursorLane;
  bool acceptsPipelines;
};

constexpr std::array<LaneSpec, 16> LaneSpecs{{
    {"notes", "notes", LaneValueKind::Notes, CursorLane::Notes, true},
    {"octave", "octave", LaneValueKind::Scalar, CursorLane::Octave, true},
    {"articulation", "articulation", LaneValueKind::Articulation,
     CursorLane::Articulation, true},
    {"art", "articulation", LaneValueKind::Articulation,
     CursorLane::Articulation, true},
    {"velocity", "velocity", LaneValueKind::Scalar, CursorLane::Velocity, true},
    {"vel", "velocity", LaneValueKind::Scalar, CursorLane::Velocity, true},
    {"accent", "accent", LaneValueKind::Scalar, CursorLane::Accent, true},
    {"duration", "duration", LaneValueKind::Scalar, CursorLane::Duration, true},
    {"dur", "duration", LaneValueKind::Scalar, CursorLane::Duration, true},
    {"gate", "gate", LaneValueKind::Scalar, CursorLane::Gate, true},
    {"slide", "slide", LaneValueKind::Scalar, CursorLane::Slide, true},
    {"ratchet", "ratchet", LaneValueKind::Scalar, CursorLane::Ratchet, true},
    {"cycle", "cycle", LaneValueKind::Cycle, CursorLane::Sequence, false},
    {"tonic", "tonic", LaneValueKind::Tonic, CursorLane::Sequence, false},
    {"scale", "scale", LaneValueKind::Scale, CursorLane::Sequence, false},
    {"glide", "glide", LaneValueKind::Glide, CursorLane::Sequence, false},
}};

const LaneSpec *FindLaneSpec(const std::string &name) {
  const auto found =
      std::find_if(LaneSpecs.begin(), LaneSpecs.end(),
                   [&](const LaneSpec &spec) { return name == spec.spelling; });
  return found == LaneSpecs.end() ? nullptr : &*found;
}

bool ParseInteger(const Token &token, int &value) {
  errno = 0;
  char *end = nullptr;
  const auto parsed = std::strtoll(token.text.c_str(), &end, 10);
  if (errno == ERANGE || end != token.text.c_str() + token.text.size() ||
      parsed < std::numeric_limits<int>::min() + 1LL ||
      parsed > std::numeric_limits<int>::max())
    return false;
  value = static_cast<int>(parsed);
  return true;
}

bool ParseSeed(const Token &token, std::uint64_t &value) {
  if (token.text.empty() || token.text.front() == '-')
    return false;
  errno = 0;
  char *end = nullptr;
  const auto parsed = std::strtoull(token.text.c_str(), &end, 10);
  if (errno == ERANGE || end != token.text.c_str() + token.text.size())
    return false;
  value = static_cast<std::uint64_t>(parsed);
  return true;
}

bool ParseTimeAmount(const Token &token, double &amount, TimeUnit &unit) {
  std::string number = token.text;
  unit = TimeUnit::Beats;
  if (number.size() > 2 && number.compare(number.size() - 2, 2, "ms") == 0) {
    unit = TimeUnit::Milliseconds;
    number.resize(number.size() - 2);
  }
  return ParseNumber(number, amount) && amount >= 0.0;
}

enum class TransformArgumentKind {
  None,
  Integer,
  PositiveNumber,
  Swing,
  Timing
};

struct TransformSpec {
  const char *spelling;
  TransformKind kind;
  TransformArgumentKind argumentKind;
  TransformDomain domain;
};

constexpr std::array<TransformSpec, 13> TransformSpecs{{
    {"rev", TransformKind::Reverse, TransformArgumentKind::None,
     TransformDomain::General},
    {"reverse", TransformKind::Reverse, TransformArgumentKind::None,
     TransformDomain::General},
    {"rotate", TransformKind::Rotate, TransformArgumentKind::Integer,
     TransformDomain::General},
    {"modulate_degree", TransformKind::ModulateDegree,
     TransformArgumentKind::Integer, TransformDomain::Pitch},
    {"shift_degree", TransformKind::ShiftDegree, TransformArgumentKind::Integer,
     TransformDomain::Pitch},
    {"transpose_semitone", TransformKind::TransposeSemitone,
     TransformArgumentKind::Integer, TransformDomain::Pitch},
    {"octave", TransformKind::TransposeOctave, TransformArgumentKind::Integer,
     TransformDomain::Pitch},
    {"transpose_octave", TransformKind::TransposeOctave,
     TransformArgumentKind::Integer, TransformDomain::Pitch},
    {"fast", TransformKind::Fast, TransformArgumentKind::PositiveNumber,
     TransformDomain::Timing},
    {"slow", TransformKind::Slow, TransformArgumentKind::PositiveNumber,
     TransformDomain::Timing},
    {"swing", TransformKind::Swing, TransformArgumentKind::Swing,
     TransformDomain::Timing},
    {"early", TransformKind::Early, TransformArgumentKind::Timing,
     TransformDomain::Timing},
    {"late", TransformKind::Late, TransformArgumentKind::Timing,
     TransformDomain::Timing},
}};

const TransformSpec *FindTransformSpec(const std::string &operation) {
  const auto found = std::find_if(
      TransformSpecs.begin(), TransformSpecs.end(),
      [&](const TransformSpec &spec) { return operation == spec.spelling; });
  return found == TransformSpecs.end() ? nullptr : &*found;
}

bool ParseTransform(const syntax::Pipeline &pipeline, Transform &transform,
                    Diagnostic &diagnostic) {
  transform.span = pipeline.span;
  if (pipeline.condition == syntax::Pipeline::Condition::Every) {
    if (!ParseInteger(pipeline.conditionArgument, transform.period) ||
        transform.period < 1) {
      diagnostic = Error(pipeline.conditionArgument.span,
                         "every requires a positive integer cycle count");
      return false;
    }
    transform.condition = TransformCondition::Every;
  } else if (pipeline.condition == syntax::Pipeline::Condition::Sometimes) {
    double probability = 0.0;
    if (!ParseNumber(pipeline.conditionArgument.text, probability) ||
        probability < 0.0 || probability > 1.0) {
      diagnostic = Error(pipeline.conditionArgument.span,
                         "sometimes requires a probability from 0 to 1");
      return false;
    }
    transform.probability = static_cast<float>(probability);
    transform.condition = TransformCondition::Sometimes;
  }

  const auto &operation = pipeline.operation.text;
  const auto &arguments = pipeline.arguments;
  const TransformSpec *spec = FindTransformSpec(operation);
  if (!spec) {
    diagnostic = Error(pipeline.operation.span,
                       "unsupported transform '" + operation + "'");
    return false;
  }
  transform.kind = spec->kind;
  transform.domain = spec->domain;
  auto requireArguments = [&](std::size_t count) {
    if (arguments.size() == count)
      return true;
    diagnostic = Error(pipeline.operation.span,
                       operation + " requires " + std::to_string(count) +
                           (count == 1 ? " argument" : " arguments"));
    return false;
  };
  auto parseIntegerArgument = [&]() {
    if (!requireArguments(1))
      return false;
    if (ParseInteger(arguments.front(), transform.integer))
      return true;
    diagnostic = Error(arguments.front().span,
                       operation + " requires an integer argument");
    return false;
  };

  if (spec->argumentKind == TransformArgumentKind::None)
    return requireArguments(0);
  if (spec->argumentKind == TransformArgumentKind::Integer)
    return parseIntegerArgument();
  if (spec->argumentKind == TransformArgumentKind::PositiveNumber) {
    if (!requireArguments(1) ||
        !ParseNumber(arguments.front().text, transform.number) ||
        transform.number <= 0.0) {
      if (diagnostic.message.empty())
        diagnostic = Error(arguments.empty() ? pipeline.operation.span
                                             : arguments.front().span,
                           operation + " requires a positive numeric factor");
      return false;
    }
    return true;
  }
  if (spec->argumentKind == TransformArgumentKind::Swing) {
    if ((arguments.size() != 1 && arguments.size() != 2) ||
        !ParseNumber(arguments.front().text, transform.number) ||
        transform.number < 0.5 || transform.number >= 1.0) {
      if (diagnostic.message.empty())
        diagnostic =
            Error(arguments.empty() ? pipeline.operation.span
                                    : arguments.front().span,
                  "swing requires RATIO or RATIO SUBDIVISION, with a ratio "
                  "from .5 up to but not including 1");
      return false;
    }
    if (arguments.size() == 2 &&
        (!ParseNumber(arguments[1].text, transform.swingSubdivisionBeats) ||
         transform.swingSubdivisionBeats <= 0.0)) {
      diagnostic =
          Error(arguments[1].span,
                "swing subdivision must be a positive incoming-clock fraction");
      return false;
    }
    return true;
  }
  if (spec->argumentKind == TransformArgumentKind::Timing) {
    std::size_t amountIndex = 0;
    if (arguments.size() == 2 && arguments.front().text == "random") {
      transform.randomAmount = true;
      amountIndex = 1;
    } else if (arguments.size() != 1) {
      diagnostic = Error(pipeline.operation.span,
                         operation + " requires AMOUNT or random AMOUNT");
      return false;
    }
    if (!ParseTimeAmount(arguments[amountIndex], transform.number,
                         transform.timeUnit)) {
      diagnostic = Error(
          arguments[amountIndex].span,
          operation + " requires a non-negative beat fraction or ms amount");
      return false;
    }
    return true;
  }
  diagnostic =
      Error(pipeline.operation.span, "invalid transform specification");
  return false;
}

std::vector<int> ScaleIntervals(const std::string &name) {
  if (name == "major" || name == "ionian")
    return {0, 2, 4, 5, 7, 9, 11};
  if (name == "minor" || name == "aeolian")
    return {0, 2, 3, 5, 7, 8, 10};
  if (name == "harmonic_minor")
    return {0, 2, 3, 5, 7, 8, 11};
  if (name == "dorian")
    return {0, 2, 3, 5, 7, 9, 10};
  if (name == "phrygian")
    return {0, 1, 3, 5, 7, 8, 10};
  if (name == "lydian")
    return {0, 2, 4, 6, 7, 9, 11};
  if (name == "mixolydian")
    return {0, 2, 4, 5, 7, 9, 10};
  if (name == "locrian")
    return {0, 1, 3, 5, 6, 8, 10};
  if (name == "minor_pentatonic")
    return {0, 3, 5, 7, 10};
  if (name == "pentatonic" || name == "major_pentatonic")
    return {0, 2, 4, 7, 9};
  if (name == "octatonic_whole_half" || name == "whole_half")
    return {0, 2, 3, 5, 6, 8, 9, 11};
  if (name == "octatonic_half_whole" || name == "half_whole")
    return {0, 1, 3, 4, 6, 7, 9, 10};
  return {};
}

std::string PatternNodeText(const syntax::PatternNode &node) {
  if (node.kind == syntax::PatternKind::Atom)
    return node.atom.text;
  if (node.kind == syntax::PatternKind::Voicing) {
    std::string text = "(";
    for (std::size_t index = 0; index < node.children.size(); ++index) {
      if (index != 0)
        text += ' ';
      text += PatternNodeText(node.children[index]);
    }
    return text + ")" + node.suffix.text;
  }
  if (node.kind == syntax::PatternKind::Slash) {
    if (node.children.size() != 2)
      return {};
    return PatternNodeText(node.children[0]) + "/" +
           PatternNodeText(node.children[1]);
  }
  if (node.kind == syntax::PatternKind::Repeat)
    return node.children.empty() ? std::string{}
                                 : PatternNodeText(node.children.front()) +
                                       "!" + node.repeatCount.text;
  if (node.kind == syntax::PatternKind::Euclidean) {
    std::string text = node.atom.text + "(";
    for (std::size_t index = 0; index < node.arguments.size(); ++index) {
      if (index != 0)
        text += ',';
      text += node.arguments[index].text;
    }
    return text + ")";
  }
  const char open = node.kind == syntax::PatternKind::CycleChoice ? '<' : '[';
  const char close = node.kind == syntax::PatternKind::CycleChoice ? '>' : ']';
  const char *separator =
      node.kind == syntax::PatternKind::RandomChoice ? "|" : " ";
  std::string text(1, open);
  for (std::size_t index = 0; index < node.children.size(); ++index) {
    if (index != 0)
      text += separator;
    text += PatternNodeText(node.children[index]);
  }
  text += close;
  return text;
}

std::string LaneText(const syntax::Lane &lane) {
  std::string text;
  for (const auto &value : lane.pattern.steps) {
    if (!text.empty())
      text += ' ';
    text += PatternNodeText(value);
  }
  return text;
}

SourceSpan LaneValueSpan(const syntax::Lane &lane) {
  if (lane.pattern.steps.empty())
    return lane.name.span;
  auto span = lane.pattern.steps.front().span;
  span.end = lane.pattern.steps.back().span.end;
  return span;
}

bool ParsePipelines(const std::vector<syntax::Pipeline> &pipelines,
                    std::vector<Transform> &transforms, CursorLane target,
                    Diagnostic &diagnostic) {
  for (const auto &pipeline : pipelines) {
    Transform transform;
    if (!ParseTransform(pipeline, transform, diagnostic))
      return false;
    const bool pitchTransform = transform.domain == TransformDomain::Pitch;
    const bool timingTransform = transform.domain == TransformDomain::Timing;
    if (target != CursorLane::Sequence &&
        ((pitchTransform && target != CursorLane::Notes) || timingTransform)) {
      diagnostic =
          Error(pipeline.operation.span,
                timingTransform ? pipeline.operation.text +
                                      " is a whole-sequence timing transform"
                                : pipeline.operation.text +
                                      " only applies to the notes lane");
      return false;
    }
    transforms.push_back(transform);
  }
  return true;
}

bool ApplySequencePipelines(const std::vector<syntax::Pipeline> &pipelines,
                            Sequence &sequence, Diagnostic &diagnostic) {
  for (const auto &pipeline : pipelines) {
    if (pipeline.operation.text == "scale") {
      if (pipeline.condition != syntax::Pipeline::Condition::Always ||
          pipeline.arguments.size() != 1) {
        diagnostic =
            Error(pipeline.operation.span,
                  "scale requires exactly one unconditional scale name");
        return false;
      }
      auto intervals = ScaleIntervals(pipeline.arguments.front().text);
      if (intervals.empty()) {
        diagnostic =
            Error(pipeline.arguments.front().span,
                  "unknown scale '" + pipeline.arguments.front().text + "'");
        return false;
      }
      sequence.scale.intervals = std::move(intervals);
      continue;
    }
    Transform transform;
    if (!ParseTransform(pipeline, transform, diagnostic))
      return false;
    if (transform.kind == TransformKind::ModulateDegree)
      transform.modulationIntervals = sequence.scale.intervals;
    const bool timingTransform = transform.domain == TransformDomain::Timing;
    const bool pitchTransform = transform.domain == TransformDomain::Pitch;
    if (timingTransform) {
      sequence.transforms[static_cast<std::size_t>(CursorLane::Sequence)]
          .push_back(transform);
    } else if (pitchTransform) {
      sequence.transforms[static_cast<std::size_t>(CursorLane::Notes)]
          .push_back(transform);
    } else {
      for (auto lane :
           {CursorLane::Notes, CursorLane::Articulation, CursorLane::Octave,
            CursorLane::Velocity, CursorLane::Accent, CursorLane::Duration,
            CursorLane::Gate, CursorLane::Slide, CursorLane::Ratchet})
        sequence.transforms[static_cast<std::size_t>(lane)].push_back(
            transform);
    }
  }
  return true;
}

bool CheckedMultiply(std::size_t left, std::size_t right, std::size_t &result) {
  if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left)
    return false;
  result = left * right;
  return true;
}

bool PrepareWorkspaces(CompiledProgram &program, Diagnostic &diagnostic) {
  const auto &sequences = program.semantic().sequences;
  std::size_t maximumEvents = 1;
  double minimumDuration = 1.0;
  double maximumEarlyBeats = 0.0;
  double maximumEarlyMilliseconds = 0.0;
  double maximumLateBeats = 0.0;
  double maximumLateMilliseconds = 0.0;
  for (const auto &sequence : sequences) {
    double fastestScale = 1.0;
    double sequenceEarlyBeats = 0.0;
    double sequenceEarlyMilliseconds = 0.0;
    double sequenceLateBeats = 0.0;
    double sequenceLateMilliseconds = 0.0;
    for (const auto &transform :
         sequence.transforms[static_cast<std::size_t>(CursorLane::Sequence)]) {
      if (transform.kind == TransformKind::Fast)
        fastestScale /= transform.number;
      if (transform.kind == TransformKind::Early) {
        if (transform.timeUnit == TimeUnit::Beats)
          sequenceEarlyBeats += transform.number;
        else
          sequenceEarlyMilliseconds += transform.number;
      } else if (transform.kind == TransformKind::Late) {
        if (transform.timeUnit == TimeUnit::Beats)
          sequenceLateBeats += transform.number;
        else
          sequenceLateMilliseconds += transform.number;
      }
    }
    if (!std::isfinite(fastestScale) || fastestScale <= 0.0) {
      diagnostic = Error(sequence.nameSpan,
                         "combined fast factors create an invalid event rate");
      return false;
    }
    maximumEarlyBeats = std::max(maximumEarlyBeats, sequenceEarlyBeats);
    maximumEarlyMilliseconds =
        std::max(maximumEarlyMilliseconds, sequenceEarlyMilliseconds);
    maximumLateBeats = std::max(maximumLateBeats, sequenceLateBeats);
    maximumLateMilliseconds =
        std::max(maximumLateMilliseconds, sequenceLateMilliseconds);
    std::size_t maximumLaneRatchet = 1;
    for (const auto &item : sequence.ratchet) {
      if (!item.isDefault)
        maximumLaneRatchet =
            std::max(maximumLaneRatchet, static_cast<std::size_t>(item.value));
    }

    std::size_t sequenceMaximum =
        sequence.articulation.empty()
            ? static_cast<std::size_t>(maximumLaneRatchet)
            : 1;
    for (const auto &step : sequence.articulation) {
      std::size_t stepEvents = 0;
      for (const auto &atom : step.atoms) {
        const auto count =
            atom.kind == ArticulationKind::Attack ||
                    atom.kind == ArticulationKind::Slide
                ? static_cast<std::size_t>(
                      atom.ratchets > 1 ? atom.ratchets : maximumLaneRatchet)
                : 1u;
        if (stepEvents > std::numeric_limits<std::size_t>::max() - count) {
          diagnostic = Error(step.span, "event density is too large");
          return false;
        }
        stepEvents += count;
      }
      sequenceMaximum = std::max(sequenceMaximum, stepEvents);
    }
    std::size_t maximumVoices = 1;
    for (const auto &item : sequence.notes) {
      for (const auto &choice : item.values) {
        maximumVoices =
            std::max(maximumVoices,
                     std::min<std::size_t>(MaximumPolyphony,
                                           choice.voices.size() +
                                               (choice.hasBass ? 1u : 0u)));
      }
    }
    if (!CheckedMultiply(sequenceMaximum, maximumVoices, sequenceMaximum)) {
      diagnostic = Error(sequence.nameSpan, "chord event density is too large");
      return false;
    }
    maximumEvents = std::max(maximumEvents, sequenceMaximum);

    for (const auto &item : sequence.duration) {
      const double value = item.isDefault ? 1.0 : item.value;
      if (value > 0.0)
        minimumDuration = std::min(minimumDuration, value * fastestScale);
    }
    if (sequence.duration.empty())
      minimumDuration = std::min(minimumDuration, fastestScale);
  }

  const double preparedEarlyBeats =
      maximumEarlyMilliseconds * 0.001 * PreparedMaximumClockRateHz;
  const double preparedLateBeats =
      maximumLateMilliseconds * 0.001 * PreparedMaximumClockRateHz;
  const double stepsAsDouble =
      std::ceil((1.0 + maximumEarlyBeats + preparedEarlyBeats +
                 maximumLateBeats + preparedLateBeats) /
                minimumDuration) +
      2.0;
  if (!std::isfinite(stepsAsDouble) ||
      stepsAsDouble >
          static_cast<double>(std::numeric_limits<std::size_t>::max())) {
    diagnostic = {"event density is too large", 1, 1};
    return false;
  }
  std::size_t scheduleCapacity = 0;
  if (!CheckedMultiply(maximumEvents, static_cast<std::size_t>(stepsAsDouble),
                       scheduleCapacity)) {
    diagnostic = {"event density is too large", 1, 1};
    return false;
  }
  program.maximumEventsPerStep = maximumEvents;
  program.scheduleCapacity = std::max<std::size_t>(1, scheduleCapacity);
  program.maximumEarlyBeats = maximumEarlyBeats;
  program.maximumEarlyMilliseconds = maximumEarlyMilliseconds;
  try {
    program.stepWorkspace.resize(program.maximumEventsPerStep);
    program.scheduleWorkspace.resize(program.scheduleCapacity);
    program.stateWorkspace.resize(sequences.size());
    program.snapshotStateWorkspace.resize(sequences.size());
    program.stateTransferOrder.resize(sequences.size());
    std::iota(program.stateTransferOrder.begin(),
              program.stateTransferOrder.end(), std::size_t{0});
    std::sort(program.stateTransferOrder.begin(),
              program.stateTransferOrder.end(),
              [&](std::size_t left, std::size_t right) {
                if (sequences[left].stableId != sequences[right].stableId)
                  return sequences[left].stableId < sequences[right].stableId;
                return sequences[left].name < sequences[right].name;
              });
  } catch (const std::bad_alloc &) {
    diagnostic = {"not enough memory for this pattern's event density", 1, 1};
    return false;
  } catch (const std::length_error &) {
    diagnostic = {"pattern event density exceeds addressable memory", 1, 1};
    return false;
  }
  return true;
}

CompileResult FinishProgram(std::unique_ptr<SemanticProgram> semantic) {
  CompileResult result;
  auto program = CompiledProgramFactory::create(std::move(*semantic));
  if (!PrepareWorkspaces(*program, result.diagnostic))
    return result;
  result.program = std::move(program);
  return result;
}

CompileResult CompileDocument(const syntax::Document &document) {
  CompileResult result;
  auto program = std::make_unique<SemanticProgram>();
  std::unordered_map<std::string, std::size_t> names;
  std::unordered_map<std::string, syntax::Assignment> assignments;
  std::string playName;
  SourceSpan playSpan;
  bool stopped = false;

  for (const auto &statement : document.statements) {
    if (const auto *play = std::get_if<syntax::PlayCommand>(&statement)) {
      playName = play->name.text;
      playSpan = play->name.span;
      stopped = false;
      continue;
    }
    if (const auto *stop = std::get_if<syntax::StopCommand>(&statement)) {
      playName.clear();
      playSpan = stop->span;
      stopped = true;
      continue;
    }
    if (const auto *seed = std::get_if<syntax::SeedCommand>(&statement)) {
      std::uint64_t parsed = 0;
      if (!ParseSeed(seed->value, parsed)) {
        result.diagnostic =
            Error(seed->value.span, "seed requires a non-negative integer");
        return result;
      }
      program->seed = parsed;
      continue;
    }
    if (const auto *assignment = std::get_if<syntax::Assignment>(&statement)) {
      if (assignments.find(assignment->name.text) != assignments.end() ||
          names.find(assignment->name.text) != names.end()) {
        result.diagnostic =
            Error(assignment->name.span,
                  "duplicate definition of '" + assignment->name.text + "'");
        return result;
      }
      assignments[assignment->name.text] = *assignment;
      continue;
    }

    const auto *definition =
        std::get_if<syntax::SequenceDefinition>(&statement);
    if (!definition)
      continue;
    Sequence sequence;
    sequence.name = definition->name.text;
    sequence.stableId = StableDefinitionId(sequence.name);
    sequence.nameSpan = definition->name.span;
    for (const auto &lane : definition->lanes) {
      const std::string laneName = lane.name.text;
      const std::string laneValue = LaneText(lane);
      const SourceSpan valueSpan = LaneValueSpan(lane);
      const LaneSpec *laneSpec = FindLaneSpec(laneName);
      if (!laneSpec) {
        result.diagnostic =
            Error(lane.name.span, "unknown sequence lane '" + laneName + "'");
        return result;
      }
      if (!laneSpec->acceptsPipelines && !lane.pipelines.empty()) {
        result.diagnostic =
            Error(lane.pipelines.front().span,
                  laneName +
                      " does not accept an inline pipeline; put whole-sequence "
                      "transforms after the closing brace");
        return result;
      }
      const CursorLane parsedLane = laneSpec->cursorLane;

      if (laneSpec->valueKind == LaneValueKind::Cycle) {
        double cycle = 0.0;
        if (!ParseNumber(laneValue, cycle) || cycle <= 0.0) {
          result.diagnostic =
              Error(valueSpan, "cycle must be a positive beat count");
          return result;
        }
        sequence.cycleBeats = cycle;
        continue;
      }
      if (laneSpec->valueKind == LaneValueKind::Tonic) {
        if (!ParseTonic(laneValue, sequence.scale)) {
          result.diagnostic = Error(
              valueSpan, "tonic must be a note name such as C, C@4, or F#@3");
          return result;
        }
        continue;
      }
      if (laneSpec->valueKind == LaneValueKind::Scale) {
        auto intervals = ScaleIntervals(laneValue);
        if (intervals.empty()) {
          result.diagnostic =
              Error(valueSpan, "unknown scale '" + laneValue + "'");
          return result;
        }
        sequence.scale.intervals = std::move(intervals);
        continue;
      }
      if (laneSpec->valueKind == LaneValueKind::Glide) {
        double glide = 0.0;
        if (!ParseNumber(laneValue, glide) || glide < 0.0 ||
            glide > std::numeric_limits<float>::max()) {
          result.diagnostic = Error(
              valueSpan,
              "glide must be a non-negative value in the supported range");
          return result;
        }
        sequence.glideBeats = static_cast<float>(glide);
        continue;
      }

      bool ok = true;
      if (laneSpec->valueKind == LaneValueKind::Notes) {
        ok = ParseNotes(lane.pattern, sequence, result.diagnostic);
      } else if (laneSpec->valueKind == LaneValueKind::Articulation) {
        ok = ParseArticulation(lane.pattern, sequence, result.diagnostic);
      } else {
        auto *items = &sequence.accent;
        if (parsedLane == CursorLane::Octave)
          items = &sequence.octave;
        else if (parsedLane == CursorLane::Velocity)
          items = &sequence.velocity;
        else if (parsedLane == CursorLane::Duration)
          items = &sequence.duration;
        else if (parsedLane == CursorLane::Gate)
          items = &sequence.gate;
        else if (parsedLane == CursorLane::Slide)
          items = &sequence.slide;
        else if (parsedLane == CursorLane::Ratchet)
          items = &sequence.ratchet;
        ok = ParseScalars(lane.pattern, *items, laneSpec->canonical,
                          result.diagnostic);
      }
      if (!ok || !ParsePipelines(
                     lane.pipelines,
                     sequence.transforms[static_cast<std::size_t>(parsedLane)],
                     parsedLane, result.diagnostic))
        return result;
    }

    if (sequence.notes.empty()) {
      result.diagnostic =
          Error(sequence.nameSpan, "sequence requires a notes lane");
      return result;
    }
    for (auto &transform :
         sequence.transforms[static_cast<std::size_t>(CursorLane::Notes)]) {
      if (transform.kind == TransformKind::ModulateDegree &&
          transform.modulationIntervals.empty())
        transform.modulationIntervals = sequence.scale.intervals;
    }
    if (!ApplySequencePipelines(definition->pipelines, sequence,
                                result.diagnostic))
      return result;

    const auto index = program->sequences.size();
    if (names.find(sequence.name) != names.end() ||
        assignments.find(sequence.name) != assignments.end()) {
      result.diagnostic = Error(sequence.nameSpan, "duplicate definition of '" +
                                                       sequence.name + "'");
      return result;
    }
    names[sequence.name] = index;
    program->sequences.push_back(std::move(sequence));
  }

  std::unordered_map<std::string, std::vector<ArrangementPart>> resolved;
  std::unordered_map<std::string, int> resolutionState;
  std::size_t generatedName = 0;
  std::function<bool(const std::string &, const SourceSpan &,
                     std::vector<ArrangementPart> &)>
      resolveName;
  std::function<bool(const syntax::Expression &, const std::string &, bool,
                     std::vector<ArrangementPart> &)>
      resolveExpression;

  auto cloneSequence = [&](std::size_t base,
                           const std::vector<syntax::Pipeline> &pipelines,
                           const std::string &name, const SourceSpan &nameSpan,
                           std::size_t &index) {
    Sequence sequence = program->sequences[base];
    sequence.name = name;
    sequence.stableId = StableDefinitionId(name);
    sequence.nameSpan = nameSpan;
    if (!ApplySequencePipelines(pipelines, sequence, result.diagnostic))
      return false;
    index = program->sequences.size();
    program->sequences.push_back(std::move(sequence));
    return true;
  };

  resolveExpression = [&](const syntax::Expression &expression,
                          const std::string &label, bool createNamedAlias,
                          std::vector<ArrangementPart> &output) {
    std::vector<ArrangementPart> parts;
    for (const auto &term : expression.terms) {
      std::vector<ArrangementPart> termParts;
      const SourceSpan termSpan =
          term.group ? term.group->span : term.name.span;
      if (term.group) {
        if (!resolveExpression(
                *term.group, label + "$group" + std::to_string(generatedName++),
                false, termParts))
          return false;
      } else if (!resolveName(term.name.text, term.name.span, termParts)) {
        return false;
      }
      if (term.repeats < 1) {
        result.diagnostic =
            Error(term.repeatSpan.valid() ? term.repeatSpan : termSpan,
                  "arrangement repeat must be a positive integer");
        return false;
      }
      if (term.hasRepeat && termParts.size() == 1) {
        const auto multiplied =
            static_cast<long long>(termParts.front().cycles) * term.repeats;
        if (multiplied > std::numeric_limits<int>::max()) {
          result.diagnostic =
              Error(term.repeatSpan, "arrangement repeat is too large");
          return false;
        }
        termParts.front().cycles = static_cast<int>(multiplied);
      } else {
        const int copies = term.hasRepeat ? term.repeats : 1;
        if (termParts.size() != 0 &&
            static_cast<std::size_t>(copies) >
                std::numeric_limits<std::size_t>::max() / termParts.size()) {
          result.diagnostic =
              Error(term.repeatSpan, "arrangement repeat is too large");
          return false;
        }
        try {
          parts.reserve(parts.size() + termParts.size() * copies);
        } catch (const std::exception &) {
          result.diagnostic =
              Error(term.repeatSpan.valid() ? term.repeatSpan : termSpan,
                    "arrangement is too large");
          return false;
        }
        for (int copy = 0; copy < copies; ++copy)
          parts.insert(parts.end(), termParts.begin(), termParts.end());
        continue;
      }
      parts.insert(parts.end(), termParts.begin(), termParts.end());
    }

    if (!expression.pipelines.empty()) {
      std::unordered_map<std::size_t, std::size_t> transformedSequences;
      for (std::size_t partIndex = 0; partIndex < parts.size(); ++partIndex) {
        const auto baseSequence = parts[partIndex].sequence;
        if (const auto transformed = transformedSequences.find(baseSequence);
            transformed != transformedSequences.end()) {
          parts[partIndex].sequence = transformed->second;
          continue;
        }
        const std::string derivedName =
            createNamedAlias && parts.size() == 1
                ? label
                : label + "$" + std::to_string(transformedSequences.size());
        std::size_t sequenceIndex = 0;
        if (!cloneSequence(baseSequence, expression.pipelines, derivedName,
                           expression.span, sequenceIndex))
          return false;
        transformedSequences.emplace(baseSequence, sequenceIndex);
        parts[partIndex].sequence = sequenceIndex;
      }
    } else if (createNamedAlias && parts.size() == 1 &&
               parts.front().cycles == 1) {
      std::size_t sequenceIndex = 0;
      if (!cloneSequence(parts.front().sequence, {}, label, expression.span,
                         sequenceIndex))
        return false;
      parts.front().sequence = sequenceIndex;
    }
    output = std::move(parts);
    return true;
  };

  resolveName = [&](const std::string &name, const SourceSpan &span,
                    std::vector<ArrangementPart> &output) {
    if (const auto sequence = names.find(name); sequence != names.end()) {
      output = {{sequence->second, 1, span}};
      return true;
    }
    if (const auto cached = resolved.find(name); cached != resolved.end()) {
      output = cached->second;
      return true;
    }
    const auto assignmentEntry = assignments.find(name);
    if (assignmentEntry == assignments.end()) {
      result.diagnostic =
          Error(span, "unknown sequence or arrangement '" + name + "'");
      return false;
    }
    if (resolutionState[name] == 1) {
      result.diagnostic =
          Error(assignmentEntry->second.name.span,
                "cyclic arrangement reference involving '" + name + "'");
      return false;
    }
    resolutionState[name] = 1;
    std::vector<ArrangementPart> parts;
    if (!resolveExpression(assignmentEntry->second.expression, name,
                           assignmentEntry->second.isDerived(), parts))
      return false;
    resolutionState[name] = 2;
    if (assignmentEntry->second.isDerived() && parts.size() == 1 &&
        parts.front().cycles == 1)
      names[name] = parts.front().sequence;
    resolved[name] = parts;
    output = std::move(parts);
    return true;
  };

  if (stopped) {
    program->stopped = true;
    return FinishProgram(std::move(program));
  }
  if (program->sequences.empty()) {
    result.diagnostic = {"program does not define a sequence", 1, 1};
    return result;
  }
  if (playName.empty())
    playName = program->sequences.front().name;

  std::vector<ArrangementPart> playedParts;
  if (!resolveName(playName, playSpan, playedParts))
    return result;
  const auto direct = names.find(playName);
  if (direct != names.end()) {
    program->arrangement.push_back(
        {direct->second, -1,
         playSpan.valid() ? playSpan
                          : program->sequences[direct->second].nameSpan});
    return FinishProgram(std::move(program));
  }

  program->arrangement = std::move(playedParts);

  return FinishProgram(std::move(program));
}

} // namespace

CompileResult Compile(const syntax::Document &document) {
  return CompileDocument(document);
}

CompileResult Compile(const std::string &source) {
  auto parsed = syntax::Parse(source);
  if (!parsed)
    return {nullptr, parsed.diagnostic};
  return Compile(parsed.document);
}

} // namespace tfseq
