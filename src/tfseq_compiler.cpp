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
#include <unordered_set>

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

std::uint64_t StableRandomIdentity(std::uint64_t definition,
                                   CursorLane lane,
                                   std::size_t position) noexcept {
  std::uint64_t hash = definition;
  auto append = [&](std::uint64_t value) {
    for (int byte = 0; byte < 8; ++byte) {
      hash ^= static_cast<unsigned char>(value & 0xffU);
      hash *= 1099511628211ULL;
      value >>= 8;
    }
  };
  append(static_cast<std::uint64_t>(lane));
  append(static_cast<std::uint64_t>(position));
  return hash == 0 ? 1 : hash;
}

void AssignRandomIdentities(Sequence &sequence) {
  for (std::size_t position = 0; position < sequence.notes.size(); ++position)
    sequence.notes[position].randomIdentity = StableRandomIdentity(
        sequence.stableId, CursorLane::Notes, position);

  std::size_t pitchedPosition = 0;
  for (auto &step : sequence.articulation) {
    for (auto &atom : step.atoms) {
      if (!atom.hasPitch)
        continue;
      atom.pitch.randomIdentity = StableRandomIdentity(
          sequence.stableId, CursorLane::Notes, pitchedPosition++);
    }
  }

  auto assignLane = [&](CursorLane lane, std::vector<ScalarItem> &items) {
    for (std::size_t position = 0; position < items.size(); ++position)
      items[position].randomIdentity =
          StableRandomIdentity(sequence.stableId, lane, position);
  };
  assignLane(CursorLane::Octave, sequence.octave);
  assignLane(CursorLane::Velocity, sequence.velocity);
  assignLane(CursorLane::Duration, sequence.duration);
  assignLane(CursorLane::Gate, sequence.gate);
  assignLane(CursorLane::Slide, sequence.slide);
  assignLane(CursorLane::Ratchet, sequence.ratchet);
  assignLane(CursorLane::Offset, sequence.offset);
  for (std::size_t cv = 0; cv < CvLaneCount; ++cv)
    assignLane(CvCursorLane(cv), sequence.cv[cv]);
}

bool IsReservedCvName(const std::string &name) noexcept {
  if (name.size() < 3 || name[0] != 'c' || name[1] != 'v' || name[2] == '0')
    return false;
  return std::all_of(name.begin() + 2, name.end(),
                     [](unsigned char value) { return std::isdigit(value); });
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

bool ParseRomanChord(const Token &token, ChordValue &chord,
                     Diagnostic &diagnostic) {
  std::string text = token.text;
  bool hasOctave = false;
  int octave = 0;
  int octaveOffset = 0;
  if (!ParseRegisterSuffix(text, hasOctave, octave, octaveOffset))
    return false;

  std::size_t cursor = 0;
  int accidental = 0;
  while (cursor < text.size() &&
         (text[cursor] == 'b' || text[cursor] == '#')) {
    accidental += text[cursor++] == '#' ? 1 : -1;
  }
  const auto numeralBegin = cursor;
  while (cursor < text.size() &&
         (text[cursor] == 'I' || text[cursor] == 'V' ||
          text[cursor] == 'i' || text[cursor] == 'v'))
    ++cursor;
  if (cursor == numeralBegin)
    return false;
  const auto numeral = text.substr(numeralBegin, cursor - numeralBegin);
  const bool lower = std::islower(static_cast<unsigned char>(numeral.front()));
  static const std::unordered_map<std::string, int> Degrees{
      {"I", 1},   {"II", 2},  {"III", 3}, {"IV", 4},
      {"V", 5},   {"VI", 6},  {"VII", 7}, {"i", 1},
      {"ii", 2},  {"iii", 3}, {"iv", 4},  {"v", 5},
      {"vi", 6},  {"vii", 7},
  };
  const auto degree = Degrees.find(numeral);
  if (degree == Degrees.end()) {
    diagnostic = Error(token.span,
                       "Roman chord degree must be I through VII");
    return false;
  }

  const std::string quality = text.substr(cursor);
  auto begins = [&](const char *prefix) { return quality.rfind(prefix, 0) == 0; };
  const bool explicitMinor = begins("min") ||
                             (begins("m") && !begins("maj"));
  const bool specialTriad =
      begins("dim") || begins("aug") || begins("sus");
  if ((!lower && explicitMinor) ||
      (lower && begins("maj"))) {
    diagnostic = Error(token.span,
                       "Roman-numeral case and chord quality contradict");
    return false;
  }
  std::string synthetic = "C";
  if (lower && quality.empty())
    synthetic += "m";
  else if (lower && !explicitMinor && !specialTriad)
    synthetic += "m" + quality;
  else
    synthetic += quality;

  ChordValue realized;
  Token syntheticToken{synthetic, token.span};
  if (!ParseJazzChord(syntheticToken, realized, diagnostic))
    return false;
  chord.meaning = ChordValue::Meaning::RomanSymbol;
  chord.jazzSymbol = token.text;
  chord.romanRoot.degree = degree->second;
  chord.romanRoot.accidental = accidental;
  chord.romanRoot.hasOctave = hasOctave;
  chord.romanRoot.octave = octave;
  chord.romanRoot.octaveOffset = octaveOffset;
  chord.romanRoot.span = token.span;
  chord.span = token.span;
  for (const auto &voice : realized.voices) {
    const int interval = voice.octaveOffset * 12 + voice.pitchClass;
    chord.intervals.push_back(interval);
  }
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
  for (const auto &child : node.children) {
    const Token &tone = child.atom;
    PitchValue pitch;
    const bool parsed = child.kind == syntax::PatternKind::NamedPitch
                            ? ParseBareNamedPitch(tone, pitch)
                            : child.kind == syntax::PatternKind::ScaleDegree
                                  ? ParsePitchValue(tone, pitch, diagnostic)
                                  : false;
    if (!parsed) {
      if (diagnostic.message.empty())
        diagnostic = Error(child.span, "invalid explicit chord tone");
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

std::string PatternNodeText(const syntax::PatternNode &node);

bool ParsePitchedChoice(const syntax::PatternNode &node, ChordValue &chord,
                        Diagnostic &diagnostic) {
  if (node.kind == syntax::PatternKind::NamedPitch) {
    PitchValue pitch;
    if (!ParseBareNamedPitch(node.atom, pitch)) {
      diagnostic = Error(node.span, "invalid named pitch");
      return false;
    }
    chord.voices.push_back(pitch);
    chord.span = node.span;
    return true;
  }
  if (node.kind == syntax::PatternKind::ScaleDegree) {
    PitchValue pitch;
    if (!ParsePitchValue(node.atom, pitch, diagnostic))
      return false;
    chord.voices.push_back(pitch);
    chord.span = node.span;
    return true;
  }
  if (node.kind == syntax::PatternKind::JazzChord)
    return ParseJazzChord(node.atom, chord, diagnostic);
  if (node.kind == syntax::PatternKind::RomanChord)
    return ParseRomanChord(node.atom, chord, diagnostic);
  if (node.kind == syntax::PatternKind::Voicing)
    return ParseExplicitChord(node, chord, diagnostic);
  if (node.kind != syntax::PatternKind::Slash || node.children.size() != 2 ||
      (node.children[1].kind != syntax::PatternKind::NamedPitch &&
       node.children[1].kind != syntax::PatternKind::ScaleDegree)) {
    diagnostic = Error(node.span, "invalid pitched value");
    return false;
  }
  if (!ParsePitchedChoice(node.children[0], chord, diagnostic))
    return false;
  if (chord.hasBass) {
    diagnostic = Error(node.span, "a chord can contain only one slash bass");
    return false;
  }
  if (node.children[1].kind == syntax::PatternKind::NamedPitch) {
    if (!ParseBareNamedPitch(node.children[1].atom, chord.bass)) {
      diagnostic = Error(node.children[1].span, "invalid slash-bass note");
      return false;
    }
  } else if (!ParsePitchValue(node.children[1].atom, chord.bass, diagnostic)) {
    return false;
  }
  chord.hasBass = true;
  chord.span = node.span;
  if (chord.meaning == ChordValue::Meaning::JazzSymbol)
    chord.jazzSymbol = PatternNodeText(node);
  const auto harmonicVoices =
      chord.meaning == ChordValue::Meaning::RomanSymbol
          ? chord.intervals.size()
          : chord.voices.size();
  if (harmonicVoices + 1 > MaximumPolyphony) {
    diagnostic = Error(node.span, "chord exceeds Rack's 16-channel polyphony");
    return false;
  }
  return true;
}

bool ParsePitchPatternNode(const syntax::PatternNode &node, PitchItem &item,
                           Diagnostic &diagnostic) {
  item.span = node.span;
  if (node.kind == syntax::PatternKind::RandomPitch) {
    const std::string distribution = node.atom.text;
    const bool chromatic = distribution == "c" || distribution == "cn";
    const bool normal = distribution == "n" || distribution == "cn";
    if (!distribution.empty() && distribution != "u" &&
        distribution != "n" && distribution != "c" &&
        distribution != "cn") {
      diagnostic = Error(node.span, "unknown random-pitch distribution");
      return false;
    }
    item.randomDomain = chromatic
                            ? PitchItem::RandomDomain::ChromaticSemitone
                            : PitchItem::RandomDomain::ScaleDegree;
    item.randomDistribution = normal
                                  ? PitchItem::RandomDistribution::Normal
                                  : PitchItem::RandomDistribution::Uniform;
    if (node.arguments.empty()) {
      if (!distribution.empty()) {
        diagnostic = Error(node.span,
                           "random-pitch distribution requires two values");
        return false;
      }
      item.randomDefaultRange = true;
      return true;
    }
    if (node.arguments.size() != 2) {
      diagnostic = Error(node.span,
                         "random pitch requires exactly two values");
      return false;
    }
    auto parseArgument = [&](const Token &token, double &value) {
      if (token.text.size() >= 2 &&
          token.text.substr(token.text.size() - 2) == "ms")
        return false;
      return ParseNumber(token.text, value);
    };
    if (!parseArgument(node.arguments[0], item.randomFirst) ||
        !parseArgument(node.arguments[1], item.randomSecond)) {
      diagnostic = Error(node.span,
                         "random-pitch values must be finite unitless numbers");
      return false;
    }
    if (normal) {
      if (item.randomSecond <= 0.0) {
        diagnostic = Error(node.arguments[1].span,
                           "normal standard deviation must be positive");
        return false;
      }
      const double spread = NormalDeviationLimit * item.randomSecond;
      if (!std::isfinite(spread) ||
          !std::isfinite(item.randomFirst - spread) ||
          !std::isfinite(item.randomFirst + spread)) {
        diagnostic = Error(node.span,
                           "normal pitch range exceeds the supported range");
        return false;
      }
    } else {
      if (std::floor(item.randomFirst) != item.randomFirst ||
          std::floor(item.randomSecond) != item.randomSecond) {
        diagnostic = Error(node.span,
                           "uniform pitch bounds must be integers");
        return false;
      }
      if (item.randomFirst > item.randomSecond) {
        diagnostic = Error(node.span,
                           "uniform pitch bounds must be low, high");
        return false;
      }
    }
    const auto minimum =
        static_cast<double>(std::numeric_limits<int>::min() + 1);
    const auto maximum = static_cast<double>(std::numeric_limits<int>::max());
    if (item.randomFirst < minimum || item.randomFirst > maximum ||
        (!normal &&
         (item.randomSecond < minimum || item.randomSecond > maximum))) {
      diagnostic = Error(node.span,
                         "random-pitch values exceed the supported range");
      return false;
    }
    return true;
  }
  ChordValue value;
  if (!ParsePitchedChoice(node, value, diagnostic))
    return false;
  item.values.push_back(std::move(value));
  return true;
}

bool ParsePositiveCount(const Token &token, const std::string &label,
                        std::size_t &count, Diagnostic &diagnostic) {
  count = 1;
  if (token.text.empty())
    return true;
  errno = 0;
  char *end = nullptr;
  const auto parsed = std::strtoull(token.text.c_str(), &end, 10);
  if (errno == ERANGE || parsed == 0 ||
      end != token.text.c_str() + token.text.size() ||
      parsed > std::numeric_limits<std::size_t>::max()) {
    diagnostic = Error(token.span, label +
                                       " must be a positive integer that fits "
                                       "addressable memory");
    return false;
  }
  count = static_cast<std::size_t>(parsed);
  return true;
}

bool ParseDurationWeight(const Token &token, double &weight,
                         Diagnostic &diagnostic) {
  weight = 1.0;
  if (token.text.empty())
    return true;
  const auto dot = token.text.find('.');
  const auto elongation = token.text.substr(0, dot);
  if (!elongation.empty()) {
    if (std::all_of(elongation.begin(), elongation.end(),
                    [](char value) { return value == '_'; })) {
      weight = static_cast<double>(elongation.size() + 1);
    } else if (elongation.front() == '_') {
      double value = 0.0;
      if (!ParseNumber(elongation.substr(1), value) || value <= 0.0) {
        diagnostic = Error(token.span, "invalid duration suffix");
        return false;
      }
      weight = value;
    }
  }
  if (dot != std::string::npos) {
    const auto dots = token.text.size() - dot;
    weight *= dots == 1 ? 1.5 : 1.75;
  }
  if (!std::isfinite(weight) || weight <= 0.0) {
    diagnostic = Error(token.span, "duration must be positive and finite");
    return false;
  }
  return true;
}

bool ParseEuclideanSuffix(const syntax::PatternNode &node, int &pulses,
                          int &steps, int &rotation,
                          Diagnostic &diagnostic) {
  if (node.arguments.empty())
    return false;
  if (node.arguments.size() < 2 || node.arguments.size() > 3) {
    diagnostic = Error(node.span, "expected (pulses,steps[,rotation])");
    return true;
  }
  auto parse = [&](const Token &argument, int &value) {
    errno = 0;
    char *end = nullptr;
    const auto parsed = std::strtoll(argument.text.c_str(), &end, 10);
    if (errno == ERANGE ||
        end != argument.text.c_str() + argument.text.size() ||
        parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max())
      return false;
    value = static_cast<int>(parsed);
    return true;
  };
  if (!parse(node.arguments[0], pulses) ||
      !parse(node.arguments[1], steps) ||
      (node.arguments.size() == 3 &&
       !parse(node.arguments[2], rotation))) {
    diagnostic = Error(node.span, "expected integer Euclidean arguments");
    return true;
  }
  if (steps < 1 || pulses < 0 || pulses > steps) {
    diagnostic = Error(node.span,
                       "Euclidean rhythm requires 0 <= pulses <= steps");
  }
  return true;
}

bool EuclideanHit(int cell, int pulses, int steps, int rotation) {
  const auto source = ((cell - rotation) % steps + steps) % steps;
  return (static_cast<std::int64_t>(source) * pulses) % steps < pulses;
}

bool ParseScalarWithUnit(const Token &token, double &value, bool &milliseconds) {
  std::string text = token.text;
  milliseconds = text.size() >= 2 && text.substr(text.size() - 2) == "ms";
  if (milliseconds)
    text.resize(text.size() - 2);
  return ParseNumber(text, value) && std::isfinite(value);
}

bool ApplyEventAttributes(const syntax::PatternNode &node,
                          ArticulationAtom &atom, double &durationWeight,
                          bool rest, Diagnostic &diagnostic) {
  std::unordered_set<std::string> seen;
  for (const auto &attribute : node.attributes) {
    if (!seen.insert(attribute.name.text).second) {
      diagnostic = Error(attribute.name.span, "duplicate event attribute '" +
                                                   attribute.name.text + "'");
      return false;
    }
    const auto &name = attribute.name.text;
    const bool hasValue = !attribute.value.text.empty();
    if (rest && name != "len") {
      diagnostic = Error(attribute.name.span,
                         "a rest accepts only the len attribute");
      return false;
    }
    if (name == "quiet" || name == "stacc" || name == "ten") {
      if (hasValue) {
        diagnostic = Error(attribute.value.span,
                           name + " is a flag and does not take a value");
        return false;
      }
      if (name == "quiet")
        atom.quiet = true;
      else if (name == "stacc")
        atom.gateArticulation = GateArticulation::Staccato;
      else
        atom.gateArticulation = GateArticulation::Tenuto;
      continue;
    }
    if (name != "len" && name != "vel" && name != "gate" &&
        name != "slide") {
      diagnostic = Error(attribute.name.span,
                         "unknown event attribute '" + name + "'");
      return false;
    }
    if (!hasValue) {
      diagnostic = Error(attribute.name.span,
                         "event attribute '" + name + "' requires a value");
      return false;
    }
    double value = 0.0;
    bool milliseconds = false;
    if (!ParseScalarWithUnit(attribute.value, value, milliseconds)) {
      diagnostic = Error(attribute.value.span, "invalid event attribute value");
      return false;
    }
    if (name == "len") {
      if (!node.durationSuffix.text.empty()) {
        diagnostic = Error(attribute.name.span,
                           "len duplicates the event duration suffix");
        return false;
      }
      if (milliseconds || value <= 0.0) {
        diagnostic = Error(attribute.value.span,
                           "len must be a positive beat value");
        return false;
      }
      durationWeight = value;
    } else if (name == "vel" && !rest) {
      if (milliseconds || value < 0.0 || value > 1.0) {
        diagnostic = Error(attribute.value.span, "vel must be from 0 to 1");
        return false;
      }
      atom.hasVelocity = true;
      atom.velocity = static_cast<float>(value);
    } else if (name == "gate" && !rest) {
      if (value < 0.0 || (!milliseconds && value > 1.0)) {
        diagnostic = Error(attribute.value.span,
                           "gate must be 0..1 or non-negative ms");
        return false;
      }
      atom.hasGate = true;
      atom.gate = static_cast<float>(value);
      atom.gateMilliseconds = milliseconds;
    } else if (name == "slide" && !rest) {
      if (node.slidePrefix.text.empty()) {
        diagnostic = Error(attribute.name.span,
                           "slide is only meaningful on a > event");
        return false;
      }
      if (value < 0.0) {
        diagnostic = Error(attribute.value.span,
                           "slide must be non-negative");
        return false;
      }
      atom.hasSlide = true;
      atom.slide = static_cast<float>(value);
      atom.slideMilliseconds = milliseconds;
    }
  }
  if (seen.count("stacc") != 0 && seen.count("ten") != 0) {
    diagnostic = Error(node.span,
                       "stacc and ten cannot be combined on one event");
    return false;
  }
  if (atom.ghost && atom.quiet) {
    diagnostic = Error(node.span, "quiet is redundant on a ghost event");
    return false;
  }
  if (atom.ghost &&
      atom.gateArticulation == GateArticulation::Staccato) {
    diagnostic = Error(node.span, "stacc is redundant on a ghost event");
    return false;
  }
  if (atom.ghost && atom.gateArticulation == GateArticulation::Tenuto) {
    diagnostic = Error(node.span, "ten contradicts a ghost event");
    return false;
  }
  return true;
}

bool ParseNodeDurationWeight(const syntax::PatternNode &node, double &weight,
                             Diagnostic &diagnostic) {
  if (!ParseDurationWeight(node.durationSuffix, weight, diagnostic))
    return false;
  for (const auto &attribute : node.attributes) {
    if (attribute.name.text != "len")
      continue;
    if (!node.durationSuffix.text.empty()) {
      diagnostic = Error(attribute.name.span,
                         "len duplicates the event duration suffix");
      return false;
    }
    bool milliseconds = false;
    if (!ParseScalarWithUnit(attribute.value, weight, milliseconds) ||
        milliseconds || weight <= 0.0) {
      diagnostic = Error(attribute.value.span,
                         "len must be a positive beat value");
      return false;
    }
  }
  return true;
}

bool ParseAtomicEvent(const syntax::PatternNode &node, Sequence &sequence,
                      ArticulationAtom &atom, double &weight,
                      Diagnostic &diagnostic) {
  atom = {};
  atom.span = node.span;
  if (!ParseDurationWeight(node.durationSuffix, weight, diagnostic))
    return false;
  if (node.kind == syntax::PatternKind::Rest) {
    atom.kind = ArticulationKind::Rest;
    return ApplyEventAttributes(node, atom, weight, true, diagnostic);
  }
  if (node.kind == syntax::PatternKind::Tie) {
    atom.kind = ArticulationKind::Tie;
    return true;
  }
  if (node.kind != syntax::PatternKind::Event || node.children.size() != 1) {
    diagnostic = Error(node.span, "invalid nested note event");
    return false;
  }
  PitchItem pitch;
  if (!ParsePitchPatternNode(node.children.front(), pitch, diagnostic))
    return false;
  pitch.span = node.span;
  atom.pitch = pitch;
  atom.hasPitch = true;
  sequence.notes.push_back(std::move(pitch));
  atom.kind = node.slidePrefix.text.empty() ? ArticulationKind::Attack
                                            : ArticulationKind::Slide;
  if (node.dynamicPrefix.text == "x") {
    if (atom.kind == ArticulationKind::Slide) {
      diagnostic = Error(node.span, "a ghost cannot also be a slide");
      return false;
    }
    atom.ghost = true;
  } else if (node.dynamicPrefix.text == "^") {
    atom.hasAccent = true;
    atom.accent = 0.88f;
  } else if (node.dynamicPrefix.text == "^^") {
    atom.hasAccent = true;
    atom.accent = 1.f;
  }
  if (!node.ratchetCount.text.empty() &&
      !ParsePositiveCount(node.ratchetCount, "ratchet count", atom.ratchets,
                          diagnostic))
    return false;
  if (atom.kind == ArticulationKind::Slide && atom.ratchets > 1) {
    diagnostic = Error(node.ratchetCount.span,
                       "a slide cannot also be ratcheted");
    return false;
  }
  atom.probability = node.defaultProbability ? 0.5f : 1.f;
  if (!node.probability.text.empty()) {
    double probability = 0.0;
    if (!ParseNumber(node.probability.text, probability) || probability < 0.0 ||
        probability > 1.0) {
      diagnostic = Error(node.probability.span,
                         "probability must be from 0 to 1");
      return false;
    }
    atom.probability = static_cast<float>(probability);
  }
  return ApplyEventAttributes(node, atom, weight, false, diagnostic);
}

bool AppendFlatGroup(const syntax::PatternNode &node, Sequence &sequence,
                     std::vector<ArticulationAtom> &atoms, double offset,
                     double span, std::uint64_t &nextProbabilityGroup,
                     Diagnostic &diagnostic) {
  std::size_t copies = 1;
  if (!ParsePositiveCount(node.repeatCount, "event replication", copies,
                          diagnostic))
    return false;
  if (copies > 1) {
    auto body = node;
    body.repeatCount = {};
    const double copySpan = span / static_cast<double>(copies);
    for (std::size_t copy = 0; copy < copies; ++copy) {
      if (!AppendFlatGroup(body, sequence, atoms,
                           offset + copySpan * static_cast<double>(copy),
                           copySpan, nextProbabilityGroup, diagnostic))
        return false;
    }
    return true;
  }

  int pulses = 0;
  int steps = 0;
  int rotation = 0;
  if (ParseEuclideanSuffix(node, pulses, steps, rotation, diagnostic)) {
    if (!diagnostic.message.empty())
      return false;
    auto body = node;
    body.arguments.clear();
    body.repeatCount = {};
    for (int cell = 0; cell < steps; ++cell) {
      const double cellOffset =
          offset + span * static_cast<double>(cell) / static_cast<double>(steps);
      const double cellSpan = span / static_cast<double>(steps);
      if (EuclideanHit(cell, pulses, steps, rotation)) {
        if (!AppendFlatGroup(body, sequence, atoms, cellOffset, cellSpan,
                             nextProbabilityGroup, diagnostic))
          return false;
      } else {
        ArticulationAtom rest;
        rest.kind = ArticulationKind::Rest;
        rest.span = node.span;
        rest.offsetFraction = cellOffset;
        rest.spanFraction = cellSpan;
        rest.cellOffset = atoms.size();
        atoms.push_back(std::move(rest));
      }
    }
    return true;
  }
  if (node.kind == syntax::PatternKind::Event ||
      node.kind == syntax::PatternKind::Rest ||
      node.kind == syntax::PatternKind::Tie) {
    double weight = 1.0;
    ArticulationAtom prototype;
    if (!ParseAtomicEvent(node, sequence, prototype, weight, diagnostic))
      return false;
    prototype.offsetFraction = offset;
    prototype.spanFraction = span;
    prototype.cellOffset = atoms.size();
    atoms.push_back(std::move(prototype));
    return true;
  }
  if (node.kind != syntax::PatternKind::Subdivision) {
    diagnostic = Error(node.span,
                       "nested random/alternate groups require equal typed "
                       "branches and are not yet executable");
    return false;
  }
  const auto firstAtom = atoms.size();
  std::vector<double> weights;
  weights.reserve(node.children.size());
  double total = 0.0;
  for (const auto &child : node.children) {
    double weight = 1.0;
    if (!ParseNodeDurationWeight(child, weight, diagnostic))
      return false;
    std::size_t copies = 1;
    if (!ParsePositiveCount(child.repeatCount, "event replication", copies,
                            diagnostic))
      return false;
    weight *= static_cast<double>(copies);
    if (!child.arguments.empty()) {
      int childPulses = 0;
      int childSteps = 0;
      int childRotation = 0;
      ParseEuclideanSuffix(child, childPulses, childSteps, childRotation,
                            diagnostic);
      if (!diagnostic.message.empty())
        return false;
      weight *= static_cast<double>(childSteps);
    }
    weights.push_back(weight);
    total += weight;
  }
  double cursor = offset;
  for (std::size_t index = 0; index < node.children.size(); ++index) {
    const double childSpan = span * weights[index] / total;
    if (!AppendFlatGroup(node.children[index], sequence, atoms, cursor,
                         childSpan, nextProbabilityGroup, diagnostic))
      return false;
    cursor += childSpan;
  }
  float groupProbability = node.defaultProbability ? 0.5f : 1.f;
  if (!node.probability.text.empty()) {
    double parsed = 0.0;
    if (!ParseNumber(node.probability.text, parsed) || parsed < 0.0 ||
        parsed > 1.0) {
      diagnostic = Error(node.probability.span,
                         "probability must be from 0 to 1");
      return false;
    }
    groupProbability = static_cast<float>(parsed);
  }
  if (groupProbability < 1.f) {
    const auto group = ++nextProbabilityGroup;
    for (std::size_t index = firstAtom; index < atoms.size(); ++index) {
      atoms[index].enclosingProbabilityGates.push_back(
          {groupProbability, group});
    }
  }
  return true;
}

std::size_t FixedVoiceCount(const PitchItem &item) {
  if (item.randomDomain != PitchItem::RandomDomain::None)
    return 1;
  std::size_t result = 0;
  for (const auto &choice : item.values) {
    const auto voices =
        (choice.meaning == ChordValue::Meaning::RomanSymbol
             ? choice.intervals.size()
             : choice.voices.size()) +
        (choice.hasBass ? 1u : 0u);
    if (voices == 0 || (result != 0 && result != voices))
      return 0;
    result = voices;
  }
  return result;
}

bool ParseNotes(const syntax::Pattern &pattern, Sequence &sequence,
                Diagnostic &diagnostic) {
  bool hasPitchedPredecessor = false;
  std::size_t predecessorVoices = 0;
  std::uint64_t nextProbabilityGroup = 0;
  auto finishStep = [&](ArticulationStep &&step) {
    step.cellCount = std::max<std::size_t>(1, step.atoms.size());
    for (const auto &atom : step.atoms) {
      if ((atom.kind == ArticulationKind::Tie ||
           atom.kind == ArticulationKind::Slide) &&
          !hasPitchedPredecessor) {
        diagnostic = Error(atom.span,
                           atom.kind == ArticulationKind::Tie
                               ? "a tie requires a preceding pitched event"
                               : "a slide requires a preceding pitched event");
        return false;
      }
      if (atom.kind == ArticulationKind::Rest)
      {
        hasPitchedPredecessor = false;
        predecessorVoices = 0;
      } else if (atom.kind == ArticulationKind::Attack ||
                 atom.kind == ArticulationKind::Slide) {
        const auto voices = atom.hasPitch ? FixedVoiceCount(atom.pitch) : 0;
        if (atom.kind == ArticulationKind::Slide &&
            (voices == 0 || voices != predecessorVoices)) {
          diagnostic = Error(
              atom.span,
              "a chord slide requires equal, statically known voice counts");
          return false;
        }
        hasPitchedPredecessor = !atom.ghost;
        predecessorVoices = atom.ghost ? 0 : voices;
      }
    }
    sequence.articulation.push_back(std::move(step));
    return true;
  };
  auto compileStep = [&](const syntax::PatternNode &body,
                         ArticulationStep &step) {
    step.span = body.span;
    if (!ParseDurationWeight(body.durationSuffix, step.durationMultiplier,
                             diagnostic))
      return false;
    if (body.kind == syntax::PatternKind::Subdivision) {
      return AppendFlatGroup(body, sequence, step.atoms, 0.0, 1.0,
                             nextProbabilityGroup, diagnostic);
    }
    if (body.kind == syntax::PatternKind::CycleChoice ||
        body.kind == syntax::PatternKind::RandomChoice) {
      PitchItem item;
      item.choice = body.kind == syntax::PatternKind::CycleChoice
                        ? PitchItem::Choice::Alternate
                        : PitchItem::Choice::Random;
      item.span = body.span;
      for (const auto &choice : body.children) {
        const syntax::PatternNode *event = &choice;
        if (choice.kind == syntax::PatternKind::Subdivision &&
            choice.children.size() == 1)
          event = &choice.children.front();
        if (event->kind != syntax::PatternKind::Event ||
            event->children.size() != 1) {
          diagnostic = Error(
              choice.span,
              "multi-event random/alternate branches are not yet executable");
          return false;
        }
        if (!event->slidePrefix.text.empty() ||
            !event->dynamicPrefix.text.empty() ||
            !event->durationSuffix.text.empty() || !event->arguments.empty() ||
            !event->ratchetCount.text.empty() ||
            !event->probability.text.empty() || event->defaultProbability ||
            !event->repeatCount.text.empty() || !event->attributes.empty()) {
          diagnostic = Error(
              event->span,
              "atomic choices currently require plain pitched alternatives");
          return false;
        }
        ChordValue value;
        if (!ParsePitchedChoice(event->children.front(), value, diagnostic))
          return false;
        item.values.push_back(std::move(value));
      }
      sequence.notes.push_back(item);
      ArticulationAtom atom;
      atom.span = body.span;
      atom.pitch = std::move(item);
      atom.hasPitch = true;
      atom.probability = body.defaultProbability ? 0.5f : 1.f;
      if (!body.probability.text.empty()) {
        double probability = 0.0;
        if (!ParseNumber(body.probability.text, probability) ||
            probability < 0.0 || probability > 1.0) {
          diagnostic = Error(body.probability.span,
                             "probability must be from 0 to 1");
          return false;
        }
        atom.probability = static_cast<float>(probability);
      }
      step.atoms.push_back(std::move(atom));
      return true;
    }
    ArticulationAtom atom;
    double ignoredWeight = 1.0;
    if (!ParseAtomicEvent(body, sequence, atom, ignoredWeight, diagnostic))
      return false;
    step.durationMultiplier = ignoredWeight;
    step.atoms.push_back(std::move(atom));
    return true;
  };

  for (const auto &source : pattern.steps) {
    std::size_t copies = 1;
    if (!ParsePositiveCount(source.repeatCount, "event replication", copies,
                            diagnostic))
      return false;
    int pulses = 0;
    int euclideanSteps = 0;
    int rotation = 0;
    const bool euclidean = ParseEuclideanSuffix(
        source, pulses, euclideanSteps, rotation, diagnostic);
    if (!diagnostic.message.empty())
      return false;
    auto body = source;
    body.repeatCount = {};
    body.arguments.clear();
    for (std::size_t copy = 0; copy < copies; ++copy) {
      const int cells = euclidean ? euclideanSteps : 1;
      for (int cell = 0; cell < cells; ++cell) {
        ArticulationStep step;
        if (!euclidean || EuclideanHit(cell, pulses, euclideanSteps, rotation)) {
          if (!compileStep(body, step))
            return false;
        } else {
          step.span = source.span;
          ArticulationAtom rest;
          rest.kind = ArticulationKind::Rest;
          rest.span = source.span;
          step.atoms.push_back(std::move(rest));
        }
        if (!finishStep(std::move(step)))
          return false;
      }
    }
  }
  return true;
}

bool ParseScalars(const syntax::Pattern &pattern,
                  std::vector<ScalarItem> &items, const std::string &lane,
                  Diagnostic &diagnostic) {
  for (const auto &source : pattern.steps) {
    std::size_t repetitions = 1;
    if (!ParsePositiveCount(source.repeatCount, lane + " repetition",
                            repetitions, diagnostic))
      return false;
    const auto *node = &source;
    if (node->kind != syntax::PatternKind::Atom &&
        node->kind != syntax::PatternKind::RandomScalar) {
      diagnostic =
          Error(node->span, lane + " grouped patterns are reserved but not "
                                   "executable in this version");
      return false;
    }
    ScalarItem item;
    item.span = source.span;
    if (node->kind == syntax::PatternKind::RandomScalar) {
      if (node->arguments.size() != 2) {
        diagnostic = Error(source.span,
                           "random scalar requires exactly two values");
        return false;
      }
      bool firstMilliseconds = false;
      bool secondMilliseconds = false;
      if (!ParseScalarWithUnit(node->arguments[0], item.randomFirst,
                               firstMilliseconds) ||
          !ParseScalarWithUnit(node->arguments[1], item.randomSecond,
                               secondMilliseconds) ||
          firstMilliseconds != secondMilliseconds) {
        diagnostic = Error(
            source.span,
            "random scalar values must be finite and use matching units");
        return false;
      }
      item.isMilliseconds = firstMilliseconds;
      if (node->atom.text == "u") {
        item.randomDistribution = ScalarItem::RandomDistribution::Uniform;
        if (item.randomFirst > item.randomSecond) {
          diagnostic = Error(source.span,
                             "uniform scalar bounds must be low, high");
          return false;
        }
      } else if (node->atom.text == "n") {
        item.randomDistribution = ScalarItem::RandomDistribution::Normal;
        if (item.randomSecond <= 0.0) {
          diagnostic = Error(node->arguments[1].span,
                             "normal standard deviation must be positive");
          return false;
        }
      } else {
        diagnostic = Error(source.span, "unknown scalar distribution");
        return false;
      }
    } else if (node->atom.text == ".") {
      item.isDefault = true;
    } else {
      std::string number = node->atom.text;
      if (number.size() >= 2 && number.substr(number.size() - 2) == "ms") {
        item.isMilliseconds = true;
        number.resize(number.size() - 2);
      }
      if (!number.empty() && number.front() == '.')
        number.insert(number.begin(), '0');
      double parsed = 0.0;
      if (!ParseNumber(number, parsed)) {
        diagnostic = Error(source.span, "invalid " + lane + " value '" +
                                            node->atom.text + "'");
        return false;
      }
      item.value = parsed;
    }
    if (!item.isDefault &&
        item.randomDistribution == ScalarItem::RandomDistribution::None &&
        !std::isfinite(item.value)) {
      diagnostic = Error(
          source.span, lane + " value is outside the supported numeric range");
      return false;
    }
    if (!item.isDefault && item.isMilliseconds && lane != "gate" &&
        lane != "slide" && lane != "offset") {
      diagnostic = Error(source.span, lane + " does not accept milliseconds");
      return false;
    }
    const bool random =
        item.randomDistribution != ScalarItem::RandomDistribution::None;
    const bool uniform =
        item.randomDistribution == ScalarItem::RandomDistribution::Uniform;
    const double validationLow = random ? item.randomFirst : item.value;
    const double validationHigh =
        random && uniform ? item.randomSecond : validationLow;
    auto setRandomDomain = [&](double minimum, double maximum) {
      item.randomMinimum = minimum;
      item.randomMaximum = maximum;
    };
    if (!item.isDefault &&
        (lane == "velocity" || lane == "vel" ||
         (lane == "gate" && !item.isMilliseconds)) &&
        (validationLow < 0.0 || validationHigh > 1.0)) {
      diagnostic = Error(source.span, lane + " must be from 0 to 1");
      return false;
    }
    if (lane == "velocity" || lane == "vel" ||
        (lane == "gate" && !item.isMilliseconds))
      setRandomDomain(0.0, 1.0);
    if (!item.isDefault && (lane == "duration" || lane == "dur") &&
        (validationLow <= 0.0 || validationHigh <= 0.0)) {
      diagnostic = Error(source.span, "duration must be positive");
      return false;
    }
    if (lane == "duration" || lane == "dur")
      setRandomDomain(1.0e-9, std::numeric_limits<double>::infinity());
    if (!item.isDefault &&
        ((lane == "gate" && item.isMilliseconds) || lane == "slide") &&
        (validationLow < 0.0 || validationHigh < 0.0)) {
      diagnostic = Error(source.span, lane + " must be non-negative");
      return false;
    }
    if ((lane == "gate" && item.isMilliseconds) || lane == "slide")
      setRandomDomain(0.0, std::numeric_limits<double>::infinity());
    if (!item.isDefault && lane == "ratchet" &&
        (validationLow < 1.0 || validationHigh < 1.0 ||
         static_cast<long double>(validationHigh) >
             std::numeric_limits<std::size_t>::max() ||
         (uniform && (std::floor(validationLow) != validationLow ||
                      std::floor(validationHigh) != validationHigh)) ||
         (!random && std::floor(item.value) != item.value))) {
      diagnostic =
          Error(source.span, "ratchet must be a positive addressable integer");
      return false;
    }
    if (lane == "ratchet") {
      item.randomInteger = true;
      setRandomDomain(
          1.0, static_cast<double>(std::numeric_limits<std::size_t>::max()));
    }
    if (!item.isDefault && lane == "octave" &&
        ((uniform && (std::floor(validationLow) != validationLow ||
                      std::floor(validationHigh) != validationHigh)) ||
         (!random && std::floor(item.value) != item.value) ||
         validationLow < std::numeric_limits<int>::min() ||
         validationHigh > std::numeric_limits<int>::max())) {
      diagnostic =
          Error(source.span, "octave must fit the supported integer range");
      return false;
    }
    if (lane == "octave") {
      item.randomInteger = true;
      setRandomDomain(static_cast<double>(std::numeric_limits<int>::min()),
                      static_cast<double>(std::numeric_limits<int>::max()));
    }
    if (item.randomDistribution == ScalarItem::RandomDistribution::Normal) {
      const double spread = NormalDeviationLimit * item.randomSecond;
      if (!std::isfinite(spread) ||
          !std::isfinite(item.randomFirst - spread) ||
          !std::isfinite(item.randomFirst + spread)) {
        diagnostic = Error(
            source.span, "normal range exceeds the supported numeric range");
        return false;
      }
    }
    for (std::size_t repetition = 0; repetition < repetitions; ++repetition)
      items.push_back(item);
  }
  return true;
}

enum class LaneValueKind {
  Notes,
  Scalar,
  Subdivision,
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

constexpr std::array<LaneSpec, 14> LaneSpecs{{
    {"notes", "notes", LaneValueKind::Notes, CursorLane::Notes, true},
    {"octave", "octave", LaneValueKind::Scalar, CursorLane::Octave, true},
    {"velocity", "velocity", LaneValueKind::Scalar, CursorLane::Velocity, true},
    {"vel", "velocity", LaneValueKind::Scalar, CursorLane::Velocity, true},
    {"duration", "duration", LaneValueKind::Scalar, CursorLane::Duration, true},
    {"dur", "duration", LaneValueKind::Scalar, CursorLane::Duration, true},
    {"gate", "gate", LaneValueKind::Scalar, CursorLane::Gate, true},
    {"slide", "slide", LaneValueKind::Scalar, CursorLane::Slide, true},
    {"ratchet", "ratchet", LaneValueKind::Scalar, CursorLane::Ratchet, true},
    {"offset", "offset", LaneValueKind::Scalar, CursorLane::Offset, true},
    {"subdiv", "subdiv", LaneValueKind::Subdivision, CursorLane::Sequence,
     false},
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

constexpr std::array<TransformSpec, 14> TransformSpecs{{
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
    {"rate", TransformKind::Rate, TransformArgumentKind::PositiveNumber,
     TransformDomain::Phase},
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
  if (node.kind == syntax::PatternKind::Atom ||
      node.kind == syntax::PatternKind::NamedPitch ||
      node.kind == syntax::PatternKind::ScaleDegree ||
      node.kind == syntax::PatternKind::JazzChord ||
      node.kind == syntax::PatternKind::RomanChord)
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
    const bool phaseTransform = transform.domain == TransformDomain::Phase;
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
    if (target == CursorLane::Sequence && phaseTransform) {
      diagnostic = Error(pipeline.operation.span,
                         "rate is lane-local; use fast or slow for a whole "
                         "sequence");
      return false;
    }
    if (target == CursorLane::Notes && phaseTransform) {
      diagnostic = Error(pipeline.operation.span,
                         "use brackets for local note density or fast/slow "
                         "for the complete sequence");
      return false;
    }
    transforms.push_back(transform);
  }
  return true;
}

bool ParseCvPipelines(const std::vector<syntax::Pipeline> &pipelines,
                      Sequence &sequence, std::size_t cvIndex,
                      CursorLane target, Diagnostic &diagnostic) {
  for (const auto &pipeline : pipelines) {
    if (pipeline.operation.text != "interp") {
      if (!ParsePipelines({pipeline},
                          sequence.transforms[static_cast<std::size_t>(target)],
                          target, diagnostic))
        return false;
      continue;
    }
    if (pipeline.condition != syntax::Pipeline::Condition::Always ||
        pipeline.arguments.empty() || pipeline.arguments.size() > 2) {
      diagnostic = Error(
          pipeline.operation.span,
          "interp requires MODE or power POSITIVE_NUMBER and is unconditional");
      return false;
    }
    const auto &mode = pipeline.arguments.front();
    if (mode.text == "step" && pipeline.arguments.size() == 1) {
      sequence.cvInterpolation[cvIndex] = CvInterpolation::Step;
    } else if (mode.text == "linear" && pipeline.arguments.size() == 1) {
      sequence.cvInterpolation[cvIndex] = CvInterpolation::Linear;
    } else if (mode.text == "smooth" && pipeline.arguments.size() == 1) {
      sequence.cvInterpolation[cvIndex] = CvInterpolation::Smooth;
    } else if (mode.text == "power" && pipeline.arguments.size() == 2) {
      double power = 0.0;
      if (!ParseNumber(pipeline.arguments[1].text, power) || power <= 0.0 ||
          !std::isfinite(power)) {
        diagnostic = Error(pipeline.arguments[1].span,
                           "interp power requires a positive finite exponent");
        return false;
      }
      sequence.cvInterpolation[cvIndex] = CvInterpolation::Power;
      sequence.cvPower[cvIndex] = power;
    } else {
      diagnostic = Error(mode.span,
                         "interp mode must be step, linear, smooth, or power P");
      return false;
    }
  }
  return true;
}

bool ValidateRomanCardinality(const Sequence &sequence,
                              Diagnostic &diagnostic) {
  for (const auto &item : sequence.notes) {
    for (const auto &choice : item.values) {
      if (choice.meaning != ChordValue::Meaning::RomanSymbol)
        continue;
      if (choice.romanRoot.degree < 1 ||
          static_cast<std::size_t>(choice.romanRoot.degree) >
              sequence.scale.intervals.size()) {
        diagnostic = Error(
            choice.romanRoot.span,
            "Roman chord degree exceeds the active scale cardinality");
        return false;
      }
    }
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
    if (transform.domain == TransformDomain::Phase) {
      diagnostic = Error(pipeline.operation.span,
                         "rate is lane-local and cannot transform a sequence");
      return false;
    }
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
            {CursorLane::Notes, CursorLane::Octave,
             CursorLane::Velocity, CursorLane::Duration,
            CursorLane::Gate, CursorLane::Slide, CursorLane::Ratchet,
             CursorLane::Offset, CursorLane::Cv1, CursorLane::Cv2,
             CursorLane::Cv3})
        sequence.transforms[static_cast<std::size_t>(lane)].push_back(
            transform);
    }
  }
  return ValidateRomanCardinality(sequence, diagnostic);
}

bool CheckedMultiply(std::size_t left, std::size_t right, std::size_t &result) {
  if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left)
    return false;
  result = left * right;
  return true;
}

struct ScalarRange {
  double minimum = 0.0;
  double maximum = 0.0;
};

ScalarRange PossibleScalarRange(const ScalarItem &item) noexcept {
  ScalarRange range{item.value, item.value};
  if (item.randomDistribution == ScalarItem::RandomDistribution::Uniform) {
    range = {item.randomFirst, item.randomSecond};
  } else if (item.randomDistribution ==
             ScalarItem::RandomDistribution::Normal) {
    const double spread = NormalDeviationLimit * item.randomSecond;
    range = {item.randomFirst - spread, item.randomFirst + spread};
    if (item.randomInteger) {
      range.minimum = std::round(range.minimum);
      range.maximum = std::round(range.maximum);
    }
  }
  range.minimum =
      std::clamp(range.minimum, item.randomMinimum, item.randomMaximum);
  range.maximum =
      std::clamp(range.maximum, item.randomMinimum, item.randomMaximum);
  return range;
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
      if (transform.kind == TransformKind::Fast ||
          transform.kind == TransformKind::Slow) {
        const double factor = transform.kind == TransformKind::Fast
                                  ? 1.0 / transform.number
                                  : transform.number;
        const bool alwaysEnabled =
            transform.condition == TransformCondition::Always ||
            (transform.condition == TransformCondition::Every &&
             transform.period == 1);
        if (alwaysEnabled || factor < 1.0)
          fastestScale *= factor;
      }
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
    double offsetEarlyBeats = 0.0;
    double offsetEarlyMilliseconds = 0.0;
    double offsetLateBeats = 0.0;
    double offsetLateMilliseconds = 0.0;
    for (const auto &item : sequence.offset) {
      if (item.isDefault)
        continue;
      const auto range = PossibleScalarRange(item);
      if (item.isMilliseconds) {
        if (range.minimum < 0.0)
          offsetEarlyMilliseconds =
              std::max(offsetEarlyMilliseconds, -range.minimum);
        if (range.maximum > 0.0)
          offsetLateMilliseconds =
              std::max(offsetLateMilliseconds, range.maximum);
      } else {
        if (range.minimum < 0.0)
          offsetEarlyBeats = std::max(offsetEarlyBeats, -range.minimum);
        if (range.maximum > 0.0)
          offsetLateBeats = std::max(offsetLateBeats, range.maximum);
      }
    }
    maximumEarlyBeats =
        std::max(maximumEarlyBeats, sequenceEarlyBeats + offsetEarlyBeats);
    maximumEarlyMilliseconds = std::max(
        maximumEarlyMilliseconds,
        sequenceEarlyMilliseconds + offsetEarlyMilliseconds);
    maximumLateBeats =
        std::max(maximumLateBeats, sequenceLateBeats + offsetLateBeats);
    maximumLateMilliseconds = std::max(
        maximumLateMilliseconds,
        sequenceLateMilliseconds + offsetLateMilliseconds);
    std::size_t maximumLaneRatchet = 1;
    for (const auto &item : sequence.ratchet) {
      if (item.isDefault)
        continue;
      const double maximum = PossibleScalarRange(item).maximum;
      if (static_cast<long double>(maximum) >
          std::numeric_limits<std::size_t>::max()) {
        diagnostic = Error(item.span, "random ratchet density is too large");
        return false;
      }
      maximumLaneRatchet =
          std::max(maximumLaneRatchet, static_cast<std::size_t>(maximum));
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
                                           (choice.meaning ==
                                                    ChordValue::Meaning::RomanSymbol
                                                ? choice.intervals.size()
                                                : choice.voices.size()) +
                                               (choice.hasBass ? 1u : 0u)));
      }
    }
    if (!CheckedMultiply(sequenceMaximum, maximumVoices, sequenceMaximum)) {
      diagnostic = Error(sequence.nameSpan, "chord event density is too large");
      return false;
    }
    maximumEvents = std::max(maximumEvents, sequenceMaximum);

    const double defaultDuration = 4.0 / sequence.subdivision;
    double shortestDurationLaneValue = defaultDuration;
    for (const auto &item : sequence.duration) {
      const double value =
          item.isDefault ? defaultDuration : PossibleScalarRange(item).minimum;
      if (value > 0.0)
        shortestDurationLaneValue =
            std::min(shortestDurationLaneValue, value);
    }
    double shortestStepMultiplier = 1.0;
    for (const auto &step : sequence.articulation) {
      shortestStepMultiplier =
          std::min(shortestStepMultiplier, step.durationMultiplier);
    }
    const double shortestPreparedStep =
        shortestDurationLaneValue * shortestStepMultiplier * fastestScale;
    if (!std::isfinite(shortestPreparedStep) || shortestPreparedStep <= 0.0) {
      diagnostic =
          Error(sequence.nameSpan,
                "combined duration and time factors create an invalid event rate");
      return false;
    }
    minimumDuration = std::min(minimumDuration, shortestPreparedStep);
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
      if (IsReservedCvName(assignment->name.text)) {
        result.diagnostic = Error(
            assignment->name.span,
            assignment->name.text +
                " is reserved for a future typed signal assignment");
        return result;
      }
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
    if (IsReservedCvName(sequence.name)) {
      result.diagnostic = Error(
          definition->name.span,
          sequence.name + " is reserved for CV lanes and signal values");
      return result;
    }
    sequence.stableId = StableDefinitionId(sequence.name);
    sequence.nameSpan = definition->name.span;
    std::array<SourceSpan, static_cast<std::size_t>(CursorLane::Count)>
        alignmentSpans{};
    std::array<syntax::Pattern::Alignment,
               static_cast<std::size_t>(CursorLane::Count)>
        alignmentModes{};
    std::unordered_set<std::string> seenLanes;
    for (const auto &lane : definition->lanes) {
      const std::string laneName = lane.name.text;
      const std::string laneValue = LaneText(lane);
      const SourceSpan valueSpan = LaneValueSpan(lane);
      LaneSpec cvSpec{};
      std::size_t cvIndex = 0;
      const bool isCv = lane.kind == syntax::Lane::Kind::Cv;
      const LaneSpec *laneSpec = FindLaneSpec(laneName);
      if (isCv) {
        errno = 0;
        char *end = nullptr;
        const auto parsed = std::strtoull(laneName.c_str() + 2, &end, 10);
        if (errno == ERANGE || end != laneName.c_str() + laneName.size() ||
             parsed < 1 || parsed > CvLaneCount) {
          result.diagnostic =
              Error(lane.name.span,
                    "this module currently provides only cv1, cv2, and cv3");
          return result;
        }
        cvIndex = static_cast<std::size_t>(parsed - 1);
        cvSpec = {laneName.c_str(), laneName.c_str(), LaneValueKind::Scalar,
                  CvCursorLane(cvIndex), true};
        laneSpec = &cvSpec;
      }
      if (!laneSpec) {
        result.diagnostic =
            Error(lane.name.span, "unknown sequence lane '" + laneName + "'");
        return result;
      }
      const std::string canonicalLane =
          isCv ? laneName : std::string(laneSpec->canonical);
      if (!seenLanes.insert(canonicalLane).second) {
        result.diagnostic =
            Error(lane.name.span, "duplicate sequence lane '" +
                                      canonicalLane + "'");
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

      if (laneSpec->valueKind == LaneValueKind::Subdivision) {
        double subdivision = 0.0;
        if (!ParseNumber(laneValue, subdivision) || subdivision <= 0.0 ||
            std::floor(subdivision) != subdivision ||
            subdivision > std::numeric_limits<int>::max()) {
          result.diagnostic =
              Error(valueSpan, "subdiv must be a positive integer");
          return result;
        }
        sequence.subdivision = static_cast<int>(subdivision);
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
      } else {
        auto *items = &sequence.velocity;
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
        else if (parsedLane == CursorLane::Offset)
          items = &sequence.offset;
        else if (IsCvCursorLane(parsedLane))
          items = &sequence.cv[CvCursorIndex(parsedLane)];
        ok = ParseScalars(lane.pattern, *items, laneSpec->canonical,
                          result.diagnostic);
        if (lane.pattern.alignment != syntax::Pattern::Alignment::Free) {
          sequence.aligned[static_cast<std::size_t>(parsedLane)] = true;
          alignmentSpans[static_cast<std::size_t>(parsedLane)] =
              lane.pattern.span;
          alignmentModes[static_cast<std::size_t>(parsedLane)] =
              lane.pattern.alignment;
        }
      }
      const bool pipelinesOk =
          isCv ? ParseCvPipelines(lane.pipelines, sequence, cvIndex, parsedLane,
                                  result.diagnostic)
               : ParsePipelines(
                     lane.pipelines,
                     sequence.transforms[static_cast<std::size_t>(parsedLane)],
                     parsedLane, result.diagnostic);
      if (!ok || !pipelinesOk)
        return result;
      if ((parsedLane == CursorLane::Offset || IsCvCursorLane(parsedLane)) &&
          lane.pattern.alignment == syntax::Pattern::Alignment::Free) {
        auto &transforms =
            sequence.transforms[static_cast<std::size_t>(parsedLane)];
        if (std::none_of(transforms.begin(), transforms.end(),
                         [](const Transform &transform) {
                           return transform.kind == TransformKind::Rate;
                         })) {
          Transform unitRate;
          unitRate.kind = TransformKind::Rate;
          unitRate.domain = TransformDomain::Phase;
          unitRate.number = 1.0;
          unitRate.span = lane.pattern.span;
          transforms.push_back(unitRate);
        }
      }
    }

    if (sequence.articulation.empty()) {
      result.diagnostic =
          Error(sequence.nameSpan, "sequence requires a notes lane");
      return result;
    }
    if (!ValidateRomanCardinality(sequence, result.diagnostic))
      return result;
    if (!sequence.ratchet.empty()) {
      for (const auto &step : sequence.articulation) {
        const auto duplicate = std::find_if(
            step.atoms.begin(), step.atoms.end(),
            [](const ArticulationAtom &atom) { return atom.ratchets > 1; });
        if (duplicate != step.atoms.end()) {
          result.diagnostic = Error(
              duplicate->span,
              "inline *N ratchet cannot be combined with a ratchet lane");
          return result;
        }
      }
    }
    std::size_t structuralCells = 0;
    for (const auto &step : sequence.articulation) {
      if (structuralCells > std::numeric_limits<std::size_t>::max() -
                                step.cellCount) {
        result.diagnostic =
            Error(step.span, "note pattern has too many structural cells");
        return result;
      }
      structuralCells += step.cellCount;
    }
    auto alignLane = [&](CursorLane lane, std::vector<ScalarItem> &items) {
      const auto index = static_cast<std::size_t>(lane);
      if (!sequence.aligned[index])
        return true;
      if (items.size() > structuralCells) {
        result.diagnostic = Error(alignmentSpans[index],
                                  "aligned lane has more values than notes");
        return false;
      }
      ScalarItem inherited;
      inherited.isDefault = true;
      inherited.span = alignmentSpans[index];
      const auto padding = structuralCells - items.size();
      const bool right = alignmentModes[index] ==
                         syntax::Pattern::Alignment::Right;
      if (right)
        items.insert(items.begin(), padding, inherited);
      else
        items.insert(items.end(), padding, inherited);
      return true;
    };
    if (!alignLane(CursorLane::Octave, sequence.octave) ||
        !alignLane(CursorLane::Velocity, sequence.velocity) ||
        !alignLane(CursorLane::Duration, sequence.duration) ||
        !alignLane(CursorLane::Gate, sequence.gate) ||
        !alignLane(CursorLane::Slide, sequence.slide) ||
        !alignLane(CursorLane::Ratchet, sequence.ratchet) ||
        !alignLane(CursorLane::Offset, sequence.offset) ||
        !alignLane(CursorLane::Cv1, sequence.cv[0]) ||
        !alignLane(CursorLane::Cv2, sequence.cv[1]) ||
        !alignLane(CursorLane::Cv3, sequence.cv[2]))
      return result;
    for (const auto lane : {CursorLane::Octave, CursorLane::Velocity,
                            CursorLane::Duration,
                            CursorLane::Gate, CursorLane::Slide,
                            CursorLane::Ratchet, CursorLane::Offset,
                            CursorLane::Cv1, CursorLane::Cv2,
                            CursorLane::Cv3}) {
      if (!sequence.aligned[static_cast<std::size_t>(lane)])
        continue;
      if (IsCvCursorLane(lane) &&
          sequence.cvInterpolation[CvCursorIndex(lane)] !=
              CvInterpolation::Step) {
        result.diagnostic = Error(
            alignmentSpans[static_cast<std::size_t>(lane)],
            "continuous interpolation is not yet supported on an "
            "edge-aligned CV lane");
        return result;
      }
      const auto &transforms =
          sequence.transforms[static_cast<std::size_t>(lane)];
      if (std::any_of(transforms.begin(), transforms.end(),
                      [](const Transform &transform) {
                        return transform.kind == TransformKind::Rate;
                      })) {
        result.diagnostic = Error(
            alignmentSpans[static_cast<std::size_t>(lane)],
            "rate cannot be combined with an edge-aligned lane");
        return result;
      }
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
    AssignRandomIdentities(sequence);

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
      } else {
        if (!resolveName(term.name.text, term.name.span, termParts))
          return false;
        // A named arrangement is an abstraction boundary for playback
        // feedback. Its expanded parts must point at the caller's term rather
        // than leaking spans from the referenced definition. Explicit groups
        // deliberately retain their inner term spans.
        for (auto &part : termParts)
          part.span = term.name.span;
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
  try {
    return CompileDocument(document);
  } catch (const std::bad_alloc &) {
    return {nullptr, {"not enough memory to compile this program", 1, 1}};
  } catch (const std::length_error &) {
    return {nullptr, {"program exceeds addressable memory", 1, 1}};
  }
}

CompileResult Compile(const std::string &source) {
  try {
    auto parsed = syntax::Parse(source);
    if (!parsed)
      return {nullptr, parsed.diagnostic};
    return Compile(parsed.document);
  } catch (const std::bad_alloc &) {
    return {nullptr, {"not enough memory to parse this program", 1, 1}};
  } catch (const std::length_error &) {
    return {nullptr, {"program exceeds addressable memory", 1, 1}};
  }
}

} // namespace tfseq
