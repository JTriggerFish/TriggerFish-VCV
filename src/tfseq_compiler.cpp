#include "tfseq.hpp"
#include "tfseq_parser.hpp"
#include "tfseq_voicing.hpp"

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

std::uint64_t StableRandomIdentity(std::uint64_t definition, CursorLane lane,
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
  for (std::size_t position = 0; position < sequence.articulation.size();
       ++position) {
    sequence.articulation[position].presenceIdentity = StableRandomIdentity(
        sequence.stableId,
        sequence.separateRhythm ? CursorLane::Rhythm : CursorLane::Notes,
        position);
  }
  for (std::size_t position = 0; position < sequence.notes.size(); ++position)
    sequence.notes[position].randomIdentity =
        StableRandomIdentity(sequence.stableId, CursorLane::Notes, position);
  for (std::size_t position = 0; position < sequence.pitchTimeline.size();
       ++position) {
    auto &step = sequence.pitchTimeline[position];
    if (step.kind == PitchTimelineStep::Kind::Pitch)
      step.pitch.randomIdentity =
          StableRandomIdentity(sequence.stableId, CursorLane::Notes, position);
  }

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

// A rhythm reference and an inline rhythm deliberately share the concise
// `rhythm VALUE` surface syntax. Identifiers which are themselves a complete
// inline rhythm would otherwise be parsed as the inline form and silently
// shadow a same-named reusable definition.
bool IsReservedRhythmSyntaxName(const std::string &name) noexcept {
  if (name == "_" || name == "x")
    return true;
  if (name.size() < 2 || name[0] != 'x' || name[1] != '_')
    return false;
  const auto suffix = name.substr(1);
  if (std::all_of(suffix.begin(), suffix.end(),
                  [](char value) { return value == '_'; }))
    return true;
  if (suffix.size() < 2 || suffix[1] < '1' || suffix[1] > '9')
    return false;
  return std::all_of(suffix.begin() + 2, suffix.end(),
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

void AddChordTone(std::vector<ChordTone> &tones, int degree, int accidental,
                  int semitones) {
  const auto found =
      std::find_if(tones.begin(), tones.end(), [&](const ChordTone &tone) {
        return tone.degree == degree;
      });
  if (found == tones.end())
    tones.push_back({degree, accidental, semitones});
  else
    *found = {degree, accidental, semitones};
}

int NaturalChordInterval(int degree) {
  switch (degree) {
  case 1:
    return 0;
  case 2:
    return 2;
  case 9:
    return 14;
  case 3:
    return 4;
  case 4:
    return 5;
  case 11:
    return 17;
  case 5:
    return 7;
  case 6:
    return 9;
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

bool ParseJazzChord(const Token &token,
                    const std::vector<syntax::PatternNode> *factorOverride,
                    ChordValue &chord, Diagnostic &diagnostic) {
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
  if (main.compare(cursor, 4, "sus2") == 0) {
    triad = Triad::Sus2;
    cursor += 4;
  } else if (main.compare(cursor, 4, "sus4") == 0) {
    triad = Triad::Sus4;
    cursor += 4;
  } else if (main.compare(cursor, 3, "sus") == 0) {
    triad = Triad::Sus4;
    cursor += 3;
  }

  struct Alteration {
    int degree = 0;
    int offset = 0;
    bool add = false;
  };
  std::vector<Alteration> alterations;
  bool altered = false;
  while (cursor < main.size()) {
    if (main.compare(cursor, 3, "alt") == 0 && cursor + 3 == main.size()) {
      altered = true;
      cursor += 3;
      break;
    }
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
  if (altered && (triad != Triad::Major || majorSeventh ||
                  (extension != 0 && extension < 7))) {
    diagnostic =
        Error(token.span, "alt requires a dominant chord such as C7alt");
    return false;
  }

  std::vector<ChordTone> tones;
  switch (triad) {
  case Triad::Major:
    tones = {{1, 0, 0}, {3, 0, 4}, {5, 0, 7}};
    break;
  case Triad::Minor:
    tones = {{1, 0, 0}, {3, -1, 3}, {5, 0, 7}};
    break;
  case Triad::Diminished:
    tones = {{1, 0, 0}, {3, -1, 3}, {5, -1, 6}};
    break;
  case Triad::Augmented:
    tones = {{1, 0, 0}, {3, 0, 4}, {5, 1, 8}};
    break;
  case Triad::Sus2:
    tones = {{1, 0, 0}, {2, 0, 2}, {5, 0, 7}};
    break;
  case Triad::Sus4:
    tones = {{1, 0, 0}, {4, 0, 5}, {5, 0, 7}};
    break;
  }
  if (extension == 5)
    tones = {{1, 0, 0}, {5, 0, 7}};
  if (extension == 6)
    AddChordTone(tones, 6, 0, 9);
  if (altered && extension == 0)
    extension = 7;
  if (extension >= 7) {
    const int seventh = triad == Triad::Diminished ? 9 : majorSeventh ? 11 : 10;
    AddChordTone(tones, 7, seventh - 11, seventh);
  }
  if (extension >= 9)
    AddChordTone(tones, 9, 0, 14);
  if (extension == 11)
    AddChordTone(tones, 11, 0, 17);
  if (extension == 13)
    AddChordTone(tones, 13, 0, 21);

  for (const auto &alteration : alterations) {
    const int degree = alteration.degree;
    const int natural = NaturalChordInterval(degree);
    auto found =
        std::find_if(tones.begin(), tones.end(), [&](const ChordTone &tone) {
          return tone.degree == degree;
        });
    if (found != tones.end() && !alteration.add)
      *found = {degree, alteration.offset, natural + alteration.offset};
    else
      AddChordTone(tones, degree, alteration.offset,
                   natural + alteration.offset);
  }
  if (altered) {
    tones.erase(std::remove_if(tones.begin(), tones.end(),
                               [](const ChordTone &tone) {
                                 return tone.degree == 5 || tone.degree == 9 ||
                                        tone.degree == 11 || tone.degree == 13;
                               }),
                tones.end());
  }

  if (factorOverride && !factorOverride->empty()) {
    std::vector<ChordTone> available = tones;
    if (altered) {
      available.push_back({9, -1, 13});
      available.push_back({9, 1, 15});
      available.push_back({5, -1, 6});
      available.push_back({5, 1, 8});
      available.push_back({11, 1, 18});
      available.push_back({13, -1, 20});
    }
    std::vector<ChordTone> selected;
    for (const auto &factorNode : *factorOverride) {
      const auto &factor = factorNode.atom;
      std::size_t factorCursor = 0;
      int accidental = 0;
      bool explicitAccidental = false;
      while (factorCursor < factor.text.size() &&
             (factor.text[factorCursor] == 'b' ||
              factor.text[factorCursor] == '#')) {
        explicitAccidental = true;
        accidental += factor.text[factorCursor++] == '#' ? 1 : -1;
      }
      int requestedDegree = 0;
      if (!ParseIntText(factor.text.substr(factorCursor), requestedDegree)) {
        diagnostic = Error(factor.span, "invalid chord factor");
        return false;
      }
      auto found = available.end();
      if (explicitAccidental) {
        const int requestedInterval =
            NaturalChordInterval(requestedDegree) + accidental;
        found = std::find_if(
            available.begin(), available.end(), [&](const ChordTone &tone) {
              return tone.degree == requestedDegree &&
                     ((tone.semitones - requestedInterval) % 12 + 12) % 12 == 0;
            });
      } else {
        found = std::find_if(available.begin(), available.end(),
                             [&](const ChordTone &tone) {
                               return tone.degree == requestedDegree;
                             });
      }
      if (found == available.end()) {
        diagnostic = Error(factor.span, "factor '" + factor.text +
                                            "' is not present in chord '" +
                                            token.text + "'");
        return false;
      }
      const int pitchClass = ((found->semitones % 12) + 12) % 12;
      const bool duplicate = std::any_of(
          selected.begin(), selected.end(), [&](const ChordTone &tone) {
            return ((tone.semitones % 12) + 12) % 12 == pitchClass;
          });
      if (duplicate) {
        diagnostic =
            Error(factor.span, "duplicate chord factor '" + factor.text + "'");
        return false;
      }
      selected.push_back(*found);
    }
    tones = std::move(selected);
    chord.hasFactorOverride = true;
  }

  std::stable_sort(tones.begin(), tones.end(),
                   [](const ChordTone &left, const ChordTone &right) {
                     return left.semitones < right.semitones;
                   });
  if (tones.size() > MaximumPolyphony) {
    diagnostic = Error(token.span, "chord exceeds Rack's 16-channel polyphony");
    return false;
  }

  chord.meaning = ChordValue::Meaning::JazzSymbol;
  chord.jazzSymbol = token.text;
  chord.rootPitchClass = root;
  chord.harmonicRoot =
      NamedChordTone(root, 0, hasOctave, octave, octaveOffset, token.span);
  chord.tones = tones;
  chord.altered = altered;
  chord.span = token.span;
  auto realizedTones = tones;
  if (altered && !chord.hasFactorOverride) {
    AddChordTone(realizedTones, 9, -1, 13);
    AddChordTone(realizedTones, 11, 1, 18);
  }
  std::stable_sort(realizedTones.begin(), realizedTones.end(),
                   [](const ChordTone &left, const ChordTone &right) {
                     return left.semitones < right.semitones;
                   });
  for (const auto &tone : realizedTones)
    chord.voices.push_back(NamedChordTone(root, tone.semitones, hasOctave,
                                          octave, octaveOffset, token.span));
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
  while (cursor < text.size() && (text[cursor] == 'b' || text[cursor] == '#')) {
    accidental += text[cursor++] == '#' ? 1 : -1;
  }
  const auto numeralBegin = cursor;
  while (cursor < text.size() && (text[cursor] == 'I' || text[cursor] == 'V' ||
                                  text[cursor] == 'i' || text[cursor] == 'v'))
    ++cursor;
  if (cursor == numeralBegin)
    return false;
  const auto numeral = text.substr(numeralBegin, cursor - numeralBegin);
  const bool lower = std::islower(static_cast<unsigned char>(numeral.front()));
  static const std::unordered_map<std::string, int> Degrees{
      {"I", 1},  {"II", 2},  {"III", 3}, {"IV", 4},  {"V", 5},
      {"VI", 6}, {"VII", 7}, {"i", 1},   {"ii", 2},  {"iii", 3},
      {"iv", 4}, {"v", 5},   {"vi", 6},  {"vii", 7},
  };
  const auto degree = Degrees.find(numeral);
  if (degree == Degrees.end()) {
    diagnostic = Error(token.span, "Roman chord degree must be I through VII");
    return false;
  }

  const std::string quality = text.substr(cursor);
  auto begins = [&](const char *prefix) {
    return quality.rfind(prefix, 0) == 0;
  };
  const bool explicitMinor = begins("min") || (begins("m") && !begins("maj"));
  const bool specialTriad = begins("dim") || begins("aug") || begins("sus");
  if ((!lower && explicitMinor) || (lower && begins("maj"))) {
    diagnostic =
        Error(token.span, "Roman-numeral case and chord quality contradict");
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
  if (!ParseJazzChord(syntheticToken, nullptr, realized, diagnostic))
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
  chord.tones = realized.tones;
  chord.altered = realized.altered;
  chord.hasFactorOverride = realized.hasFactorOverride;
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
    return ParseJazzChord(node.atom, &node.children, chord, diagnostic);
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
  const auto harmonicVoices = chord.meaning == ChordValue::Meaning::RomanSymbol
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
    if (!distribution.empty() && distribution != "u" && distribution != "n" &&
        distribution != "c" && distribution != "cn") {
      diagnostic = Error(node.span, "unknown random-pitch distribution");
      return false;
    }
    item.randomDomain = chromatic ? PitchItem::RandomDomain::ChromaticSemitone
                                  : PitchItem::RandomDomain::ScaleDegree;
    item.randomDistribution = normal ? PitchItem::RandomDistribution::Normal
                                     : PitchItem::RandomDistribution::Uniform;
    if (node.arguments.empty()) {
      if (!distribution.empty()) {
        diagnostic =
            Error(node.span, "random-pitch distribution requires two values");
        return false;
      }
      item.randomDefaultRange = true;
      return true;
    }
    if (node.arguments.size() != 2) {
      diagnostic = Error(node.span, "random pitch requires exactly two values");
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
      if (!std::isfinite(spread) || !std::isfinite(item.randomFirst - spread) ||
          !std::isfinite(item.randomFirst + spread)) {
        diagnostic =
            Error(node.span, "normal pitch range exceeds the supported range");
        return false;
      }
    } else {
      if (std::floor(item.randomFirst) != item.randomFirst ||
          std::floor(item.randomSecond) != item.randomSecond) {
        diagnostic = Error(node.span, "uniform pitch bounds must be integers");
        return false;
      }
      if (item.randomFirst > item.randomSecond) {
        diagnostic = Error(node.span, "uniform pitch bounds must be low, high");
        return false;
      }
    }
    const auto minimum =
        static_cast<double>(std::numeric_limits<int>::min() + 1);
    const auto maximum = static_cast<double>(std::numeric_limits<int>::max());
    if (item.randomFirst < minimum || item.randomFirst > maximum ||
        (!normal &&
         (item.randomSecond < minimum || item.randomSecond > maximum))) {
      diagnostic =
          Error(node.span, "random-pitch values exceed the supported range");
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

void PromoteBareNamedPitchToMajorChord(ChordValue &chord) {
  if (chord.meaning != ChordValue::Meaning::SinglePitch ||
      chord.voices.size() != 1 || !chord.voices.front().absolute)
    return;
  const auto root = chord.voices.front();
  chord.meaning = ChordValue::Meaning::JazzSymbol;
  chord.rootPitchClass = root.pitchClass;
  chord.harmonicRoot = root;
  chord.tones = {{1, 0, 0}, {3, 0, 4}, {5, 0, 7}};
  chord.voices.clear();
  for (const auto &tone : chord.tones)
    chord.voices.push_back(NamedChordTone(root.pitchClass, tone.semitones,
                                          root.hasOctave, root.octave,
                                          root.octaveOffset, chord.span));
}

void PromoteChordLaneItem(PitchItem &item) {
  if (item.randomDomain != PitchItem::RandomDomain::None)
    return;
  for (auto &choice : item.values)
    PromoteBareNamedPitchToMajorChord(choice);
}

void PromoteBareNamedChordLaneValues(Sequence &sequence) {
  for (auto &item : sequence.notes)
    PromoteChordLaneItem(item);
  for (auto &step : sequence.pitchTimeline)
    if (step.kind == PitchTimelineStep::Kind::Pitch)
      PromoteChordLaneItem(step.pitch);
  for (auto &step : sequence.articulation)
    for (auto &atom : step.atoms)
      if (atom.hasPitch)
        PromoteChordLaneItem(atom.pitch);
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
    diagnostic =
        Error(token.span, label + " must be a positive integer that fits "
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
                          int &steps, int &rotation, Diagnostic &diagnostic) {
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
  if (!parse(node.arguments[0], pulses) || !parse(node.arguments[1], steps) ||
      (node.arguments.size() == 3 && !parse(node.arguments[2], rotation))) {
    diagnostic = Error(node.span, "expected integer Euclidean arguments");
    return true;
  }
  if (steps < 1 || pulses < 0 || pulses > steps) {
    diagnostic =
        Error(node.span, "Euclidean rhythm requires 0 <= pulses <= steps");
  }
  return true;
}

bool EuclideanHit(int cell, int pulses, int steps, int rotation) {
  const auto source = ((cell - rotation) % steps + steps) % steps;
  return (static_cast<std::int64_t>(source) * pulses) % steps < pulses;
}

enum class ValueUnit { Plain, Milliseconds, NoteValue };

bool ParseNoteValue(const std::string &text, double &beats) {
  std::size_t suffixLength = 0;
  double multiplier = 1.0;
  if (text.size() > 2 && text.compare(text.size() - 2, 2, "nd") == 0) {
    suffixLength = 2;
    multiplier = 1.5;
  } else if (text.size() > 2 && text.compare(text.size() - 2, 2, "nt") == 0) {
    suffixLength = 2;
    multiplier = 2.0 / 3.0;
  } else if (text.size() > 1 && text.back() == 'n') {
    suffixLength = 1;
  } else {
    return false;
  }
  const std::string denominatorText =
      text.substr(0, text.size() - suffixLength);
  if (denominatorText.empty())
    return false;
  errno = 0;
  char *end = nullptr;
  const auto denominator = std::strtoll(denominatorText.c_str(), &end, 10);
  if (errno == ERANGE ||
      end != denominatorText.c_str() + denominatorText.size() ||
      denominator == 0 || denominator == std::numeric_limits<long long>::min())
    return false;
  const auto magnitude =
      static_cast<std::uint64_t>(denominator < 0 ? -denominator : denominator);
  if (magnitude > 128 || (magnitude & (magnitude - 1)) != 0)
    return false;
  beats = std::copysign(4.0 * multiplier / static_cast<double>(magnitude),
                        static_cast<double>(denominator));
  return std::isfinite(beats);
}

bool ParseScalarWithUnit(const Token &token, double &value, ValueUnit &unit) {
  std::string text = token.text;
  unit = ValueUnit::Plain;
  if (text.size() >= 2 && text.substr(text.size() - 2) == "ms") {
    unit = ValueUnit::Milliseconds;
    text.resize(text.size() - 2);
  } else if (ParseNoteValue(text, value)) {
    unit = ValueUnit::NoteValue;
    return true;
  }
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
      diagnostic =
          Error(attribute.name.span, "a rest accepts only the len attribute");
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
    if (name != "len" && name != "vel" && name != "gate" && name != "slide") {
      diagnostic =
          Error(attribute.name.span, "unknown event attribute '" + name + "'");
      return false;
    }
    if (!hasValue) {
      diagnostic = Error(attribute.name.span,
                         "event attribute '" + name + "' requires a value");
      return false;
    }
    double value = 0.0;
    ValueUnit unit = ValueUnit::Plain;
    if (!ParseScalarWithUnit(attribute.value, value, unit)) {
      diagnostic = Error(attribute.value.span, "invalid event attribute value");
      return false;
    }
    if (name == "len") {
      if (!node.durationSuffix.text.empty()) {
        diagnostic = Error(attribute.name.span,
                           "len duplicates the event duration suffix");
        return false;
      }
      if (unit != ValueUnit::Plain || value <= 0.0) {
        diagnostic = Error(attribute.value.span,
                           "len must be a positive span multiplier");
        return false;
      }
      durationWeight = value;
    } else if (name == "vel" && !rest) {
      if (unit != ValueUnit::Plain || value < 0.0 || value > 1.0) {
        diagnostic = Error(attribute.value.span, "vel must be from 0 to 1");
        return false;
      }
      atom.hasVelocity = true;
      atom.velocity = static_cast<float>(value);
    } else if (name == "gate" && !rest) {
      if (value < 0.0 || (unit == ValueUnit::Plain && value > 1.0)) {
        diagnostic = Error(attribute.value.span,
                           "gate must be 0..1 or a non-negative time");
        return false;
      }
      atom.hasGate = true;
      atom.gate = static_cast<float>(value);
      atom.gateMilliseconds = unit == ValueUnit::Milliseconds;
      atom.gateNoteValue = unit == ValueUnit::NoteValue;
    } else if (name == "slide" && !rest) {
      if (node.slidePrefix.text.empty()) {
        diagnostic =
            Error(attribute.name.span, "slide is only meaningful on a > event");
        return false;
      }
      if (value < 0.0) {
        diagnostic = Error(attribute.value.span, "slide must be non-negative");
        return false;
      }
      atom.hasSlide = true;
      atom.slide = static_cast<float>(value);
      atom.slideMilliseconds = unit == ValueUnit::Milliseconds;
    }
  }
  if (seen.count("stacc") != 0 && seen.count("ten") != 0) {
    diagnostic =
        Error(node.span, "stacc and ten cannot be combined on one event");
    return false;
  }
  if (atom.ghost && atom.quiet) {
    diagnostic = Error(node.span, "quiet is redundant on a ghost event");
    return false;
  }
  if (atom.ghost && atom.gateArticulation == GateArticulation::Staccato) {
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
    ValueUnit unit = ValueUnit::Plain;
    if (!ParseScalarWithUnit(attribute.value, weight, unit) ||
        unit != ValueUnit::Plain || weight <= 0.0) {
      diagnostic =
          Error(attribute.value.span, "len must be a positive span multiplier");
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
    diagnostic =
        Error(node.ratchetCount.span, "a slide cannot also be ratcheted");
    return false;
  }
  atom.probability = node.defaultProbability ? 0.5f : 1.f;
  if (!node.probability.text.empty()) {
    double probability = 0.0;
    if (!ParseNumber(node.probability.text, probability) || probability < 0.0 ||
        probability > 1.0) {
      diagnostic =
          Error(node.probability.span, "probability must be from 0 to 1");
      return false;
    }
    atom.probability = static_cast<float>(probability);
  }
  return ApplyEventAttributes(node, atom, weight, false, diagnostic);
}

bool ParsePresenceProbability(const syntax::PatternNode &node, float &value,
                              Diagnostic &diagnostic) {
  value = node.defaultPresenceProbability ? 0.5f : 1.f;
  if (node.presenceProbability.text.empty())
    return true;
  double parsed = 0.0;
  if (!ParseNumber(node.presenceProbability.text, parsed) || parsed < 0.0 ||
      parsed > 1.0) {
    diagnostic = Error(node.presenceProbability.span,
                       "presence probability must be from 0 to 1");
    return false;
  }
  value = static_cast<float>(parsed);
  return true;
}

bool HasPresenceProbability(const syntax::PatternNode &node) {
  if (!node.presenceProbability.text.empty() || node.defaultPresenceProbability)
    return true;
  return std::any_of(node.children.begin(), node.children.end(),
                     HasPresenceProbability);
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
      const double cellOffset = offset + span * static_cast<double>(cell) /
                                             static_cast<double>(steps);
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
    diagnostic =
        Error(node.span, "nested random/alternate groups require equal typed "
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
      diagnostic =
          Error(node.probability.span, "probability must be from 0 to 1");
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

std::size_t FixedVoiceCount(const PitchItem &item, const VoicingStyle style) {
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
      if (atom.kind == ArticulationKind::Rest &&
          step.presenceProbability >= 1.f) {
        hasPitchedPredecessor = false;
        predecessorVoices = 0;
      } else if (atom.kind == ArticulationKind::Attack ||
                 atom.kind == ArticulationKind::Slide) {
        const auto voices =
            atom.hasPitch ? FixedVoiceCount(atom.pitch, sequence.voicing) : 0;
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
    if (!ParsePresenceProbability(body, step.presenceProbability, diagnostic))
      return false;
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
            !event->presenceProbability.text.empty() ||
            event->defaultPresenceProbability ||
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
          diagnostic =
              Error(body.probability.span, "probability must be from 0 to 1");
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
    const auto nestedPresence =
        std::any_of(source.children.begin(), source.children.end(),
                    [](const syntax::PatternNode &child) {
                      return HasPresenceProbability(child);
                    });
    if (nestedPresence) {
      diagnostic = Error(
          source.span,
          "presence probability is only supported on top-level events and "
          "complete top-level groups");
      return false;
    }
    std::size_t copies = 1;
    if (!ParsePositiveCount(source.repeatCount, "event replication", copies,
                            diagnostic))
      return false;
    int pulses = 0;
    int euclideanSteps = 0;
    int rotation = 0;
    const bool euclidean = ParseEuclideanSuffix(source, pulses, euclideanSteps,
                                                rotation, diagnostic);
    if (!diagnostic.message.empty())
      return false;
    if (euclidean && (!source.presenceProbability.text.empty() ||
                      source.defaultPresenceProbability)) {
      diagnostic =
          Error(source.span, "presence probability cannot be combined with a "
                             "Euclidean suffix");
      return false;
    }
    auto body = source;
    body.repeatCount = {};
    body.arguments.clear();
    for (std::size_t copy = 0; copy < copies; ++copy) {
      const int cells = euclidean ? euclideanSteps : 1;
      for (int cell = 0; cell < cells; ++cell) {
        ArticulationStep step;
        if (!euclidean ||
            EuclideanHit(cell, pulses, euclideanSteps, rotation)) {
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
  if (std::none_of(sequence.articulation.begin(), sequence.articulation.end(),
                   [](const ArticulationStep &step) {
                     return step.presenceProbability >= 1.f;
                   })) {
    diagnostic = Error(pattern.span,
                       "a sequence needs at least one event without optional "
                       "presence");
    return false;
  }
  return true;
}

bool ParsePitchTimeline(const syntax::Pattern &pattern, Sequence &sequence,
                        Diagnostic &diagnostic) {
  Sequence parsed;
  parsed.scale = sequence.scale;
  if (!ParseNotes(pattern, parsed, diagnostic))
    return false;
  bool hasPitch = false;
  for (const auto &step : parsed.articulation) {
    if (step.atoms.size() != 1 || step.presenceProbability != 1.f) {
      diagnostic =
          Error(step.span,
                "with a rhythm gesture, notes/chords is a held-value timeline; "
                "use one pitch, chord, rest, or tie per top-level span");
      return false;
    }
    const auto &atom = step.atoms.front();
    if (atom.kind == ArticulationKind::Tie) {
      if (sequence.pitchTimeline.empty()) {
        diagnostic =
            Error(step.span, "a held-value timeline cannot begin with a tie");
        return false;
      }
      sequence.pitchTimeline.back().durationMultiplier +=
          step.durationMultiplier;
      continue;
    }
    if (atom.kind == ArticulationKind::Rest) {
      PitchTimelineStep value;
      value.kind = PitchTimelineStep::Kind::Rest;
      value.durationMultiplier = step.durationMultiplier;
      value.span = step.span;
      sequence.pitchTimeline.push_back(std::move(value));
      continue;
    }
    if (!atom.hasPitch || atom.ratchets != 1 || atom.probability != 1.f ||
        atom.ghost || atom.quiet ||
        atom.gateArticulation != GateArticulation::Normal || atom.hasVelocity ||
        atom.hasAccent || atom.hasGate) {
      diagnostic = Error(
          step.span,
          "with a rhythm gesture, attack probability, dynamics, gate, and "
          "ratchets belong to the rhythm or their dedicated lanes");
      return false;
    }
    PitchTimelineStep value;
    value.kind = PitchTimelineStep::Kind::Pitch;
    value.pitch = atom.pitch;
    value.durationMultiplier = step.durationMultiplier;
    value.slideOnEntry = atom.kind == ArticulationKind::Slide;
    value.hasSlide = atom.hasSlide;
    value.slide = atom.slide;
    value.slideMilliseconds = atom.slideMilliseconds;
    value.span = step.span;
    sequence.pitchTimeline.push_back(std::move(value));
    hasPitch = true;
  }
  return hasPitch && !sequence.pitchTimeline.empty();
}

bool ParseRhythm(const syntax::Pattern &pattern, Sequence &sequence,
                 Diagnostic &diagnostic) {
  Sequence parsed;
  parsed.scale = sequence.scale;
  if (!ParseNotes(pattern, parsed, diagnostic))
    return false;
  for (auto &step : parsed.articulation) {
    for (auto &atom : step.atoms) {
      if (atom.kind == ArticulationKind::Attack ||
          atom.kind == ArticulationKind::Slide) {
        atom.hasPitch = false;
        atom.pitch = {};
      }
    }
  }
  sequence.articulation = std::move(parsed.articulation);
  sequence.separateRhythm = true;
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
        diagnostic =
            Error(source.span, "random scalar requires exactly two values");
        return false;
      }
      ValueUnit firstUnit = ValueUnit::Plain;
      ValueUnit secondUnit = ValueUnit::Plain;
      if (!ParseScalarWithUnit(node->arguments[0], item.randomFirst,
                               firstUnit) ||
          !ParseScalarWithUnit(node->arguments[1], item.randomSecond,
                               secondUnit) ||
          firstUnit != secondUnit) {
        diagnostic =
            Error(source.span,
                  "random scalar values must be finite and use matching units");
        return false;
      }
      item.isMilliseconds = firstUnit == ValueUnit::Milliseconds;
      item.isNoteValue = firstUnit == ValueUnit::NoteValue;
      if (node->atom.text == "u") {
        item.randomDistribution = ScalarItem::RandomDistribution::Uniform;
        if (item.randomFirst > item.randomSecond) {
          diagnostic =
              Error(source.span, "uniform scalar bounds must be low, high");
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
      ValueUnit unit = ValueUnit::Plain;
      double parsed = 0.0;
      if (!ParseScalarWithUnit(node->atom, parsed, unit)) {
        diagnostic = Error(source.span, "invalid " + lane + " value '" +
                                            node->atom.text + "'");
        return false;
      }
      item.value = parsed;
      item.isMilliseconds = unit == ValueUnit::Milliseconds;
      item.isNoteValue = unit == ValueUnit::NoteValue;
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
    if (!item.isDefault && item.isNoteValue && lane != "duration" &&
        lane != "dur" && lane != "gate" && lane != "slide" &&
        lane != "offset") {
      diagnostic = Error(source.span, lane + " does not accept note values");
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
         (lane == "gate" && !item.isMilliseconds && !item.isNoteValue)) &&
        (validationLow < 0.0 || validationHigh > 1.0)) {
      diagnostic = Error(source.span, lane + " must be from 0 to 1");
      return false;
    }
    if (lane == "velocity" || lane == "vel" ||
        (lane == "gate" && !item.isMilliseconds && !item.isNoteValue))
      setRandomDomain(0.0, 1.0);
    if (!item.isDefault && (lane == "duration" || lane == "dur") &&
        (validationLow <= 0.0 || validationHigh <= 0.0)) {
      diagnostic = Error(source.span, "duration must be positive");
      return false;
    }
    if (lane == "duration" || lane == "dur")
      setRandomDomain(1.0e-9, std::numeric_limits<double>::infinity());
    if (!item.isDefault &&
        ((lane == "gate" && (item.isMilliseconds || item.isNoteValue)) ||
         lane == "slide") &&
        (validationLow < 0.0 || validationHigh < 0.0)) {
      diagnostic = Error(source.span, lane + " must be non-negative");
      return false;
    }
    if ((lane == "gate" && (item.isMilliseconds || item.isNoteValue)) ||
        lane == "slide")
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
      if (!std::isfinite(spread) || !std::isfinite(item.randomFirst - spread) ||
          !std::isfinite(item.randomFirst + spread)) {
        diagnostic = Error(source.span,
                           "normal range exceeds the supported numeric range");
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
  Chords,
  Rhythm,
  Scalar,
  Subdivision,
  Tonic,
  Scale,
  Key,
  Voicing,
  Glide
};

struct LaneSpec {
  const char *spelling;
  const char *canonical;
  LaneValueKind valueKind;
  CursorLane cursorLane;
  bool acceptsPipelines;
};

constexpr std::array<LaneSpec, 18> LaneSpecs{{
    {"notes", "notes", LaneValueKind::Notes, CursorLane::Notes, true},
    {"chords", "notes", LaneValueKind::Chords, CursorLane::Notes, true},
    {"rhythm", "rhythm", LaneValueKind::Rhythm, CursorLane::Rhythm, true},
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
    {"key", "key", LaneValueKind::Key, CursorLane::Sequence, false},
    {"voicing", "voicing", LaneValueKind::Voicing, CursorLane::Sequence, false},
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
  ValueUnit parsedUnit = ValueUnit::Plain;
  if (!ParseScalarWithUnit(token, amount, parsedUnit) || amount < 0.0)
    return false;
  unit = parsedUnit == ValueUnit::Milliseconds ? TimeUnit::Milliseconds
                                               : TimeUnit::Beats;
  return true;
}

enum class TransformArgumentKind {
  None,
  Integer,
  Key,
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

constexpr std::array<TransformSpec, 15> TransformSpecs{{
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
    {"transpose_key", TransformKind::TransposeKey, TransformArgumentKind::Key,
     TransformDomain::Pitch},
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
  if (spec->argumentKind == TransformArgumentKind::Key) {
    if (!requireArguments(1))
      return false;
    std::size_t cursor = 0;
    if (ParseNoteName(arguments.front().text, cursor, transform.integer) &&
        cursor == arguments.front().text.size())
      return true;
    diagnostic = Error(arguments.front().span,
                       operation + " requires a key such as C, F#, or Bb");
    return false;
  }
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
    ValueUnit subdivisionUnit = ValueUnit::Plain;
    if (arguments.size() == 2 &&
        (!ParseScalarWithUnit(arguments[1], transform.swingSubdivisionBeats,
                              subdivisionUnit) ||
         subdivisionUnit == ValueUnit::Milliseconds ||
         transform.swingSubdivisionBeats <= 0.0)) {
      diagnostic =
          Error(arguments[1].span,
                "swing subdivision must be a positive beat or note value");
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
      diagnostic =
          Error(arguments[amountIndex].span,
                operation +
                    " requires a non-negative beat, note value, or ms amount");
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
  if (name == "major_pentatonic")
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
  if (lane.envelopeOnly)
    return lane.envelopeSpan;
  if (lane.pattern.steps.empty())
    return lane.name.span;
  auto span = lane.pattern.steps.front().span;
  span.end = lane.pattern.steps.back().span.end;
  return span;
}

bool ParseEnvelopeTime(const Token &token, CvEnvelopeTime &time,
                       Diagnostic &diagnostic) {
  std::string number = token.text;
  if (number.size() > 1 && number.back() == 's' &&
      !(number.size() > 2 && number.compare(number.size() - 2, 2, "ms") == 0)) {
    number.resize(number.size() - 1);
    double value = 0.0;
    if (!ParseNumber(number, value) || value < 0.0 || !std::isfinite(value)) {
      diagnostic =
          Error(token.span, "envelope time must be a non-negative duration");
      return false;
    }
    time = {value, CvEnvelopeTimeUnit::Seconds};
    return true;
  }
  double value = 0.0;
  ValueUnit parsedUnit = ValueUnit::Plain;
  if (!ParseScalarWithUnit(token, value, parsedUnit) || value < 0.0) {
    diagnostic =
        Error(token.span, "envelope time must be a non-negative duration such "
                          "as 5ms, .3s, 16n, or 1/4 beats");
    return false;
  }
  const auto unit = parsedUnit == ValueUnit::Milliseconds
                        ? CvEnvelopeTimeUnit::Seconds
                        : CvEnvelopeTimeUnit::Beats;
  const double scale = parsedUnit == ValueUnit::Milliseconds ? 0.001 : 1.0;
  time = {value * scale, unit};
  return true;
}

bool ParseCvEnvelope(const std::vector<Token> &arguments,
                     CvEnvelopeComposition composition, const SourceSpan &span,
                     CvEnvelopeSpec &spec, Diagnostic &diagnostic) {
  if (arguments.empty()) {
    diagnostic = Error(span, "env requires ad, ar, or adsr");
    return false;
  }

  CvEnvelopeSpec parsed;
  parsed.enabled = true;
  parsed.composition = composition;
  parsed.span = span;
  std::size_t positionalCount = 0;
  if (arguments[0].text == "ad") {
    parsed.mode = CvEnvelopeMode::Ad;
  } else if (arguments[0].text == "ar") {
    parsed.mode = CvEnvelopeMode::Ar;
  } else if (arguments[0].text == "adsr") {
    parsed.mode = CvEnvelopeMode::Adsr;
    parsed.decay = {0.250, CvEnvelopeTimeUnit::Seconds};
  } else {
    diagnostic =
        Error(arguments[0].span, "envelope mode must be ad, ar, or adsr");
    return false;
  }

  bool depthSeen = false;
  bool curveSeen = false;
  bool followSeen = false;
  bool accentSeen = false;
  bool optionsStarted = false;
  auto requireValue = [&](std::size_t index, const char *option) {
    if (index + 1 < arguments.size())
      return true;
    diagnostic =
        Error(arguments[index].span, std::string(option) + " requires a value");
    return false;
  };
  for (std::size_t index = 1; index < arguments.size();) {
    const auto &argument = arguments[index];
    if (argument.text == "depth") {
      optionsStarted = true;
      if (depthSeen || !requireValue(index, "depth")) {
        if (depthSeen)
          diagnostic = Error(argument.span, "depth may be specified once");
        return false;
      }
      double value = 0.0;
      if (!ParseNumber(arguments[index + 1].text, value) ||
          !std::isfinite(value) ||
          std::abs(value) > std::numeric_limits<float>::max()) {
        diagnostic =
            Error(arguments[index + 1].span, "depth requires a finite voltage");
        return false;
      }
      parsed.depth = static_cast<float>(value);
      depthSeen = true;
      index += 2;
      continue;
    }
    if (argument.text == "curve") {
      optionsStarted = true;
      if (curveSeen || !requireValue(index, "curve")) {
        if (curveSeen)
          diagnostic = Error(argument.span, "curve may be specified once");
        return false;
      }
      double value = 0.0;
      if (!ParseNumber(arguments[index + 1].text, value) || value < -1.0 ||
          value > 1.0) {
        diagnostic =
            Error(arguments[index + 1].span, "curve must be between -1 and 1");
        return false;
      }
      parsed.curve = static_cast<float>(value);
      curveSeen = true;
      index += 2;
      continue;
    }
    if (argument.text == "follow") {
      optionsStarted = true;
      if (followSeen || !requireValue(index, "follow")) {
        if (followSeen)
          diagnostic = Error(argument.span, "follow may be specified once");
        return false;
      }
      const auto &value = arguments[index + 1];
      if (value.text != "velocity" && value.text != "vel") {
        diagnostic =
            Error(value.span, "follow currently accepts velocity or vel");
        return false;
      }
      parsed.followVelocity = true;
      followSeen = true;
      index += 2;
      continue;
    }
    if (argument.text == "accent") {
      optionsStarted = true;
      if (accentSeen || !requireValue(index, "accent")) {
        if (accentSeen)
          diagnostic = Error(argument.span, "accent may be specified once");
        return false;
      }
      double value = 0.0;
      if (!ParseNumber(arguments[index + 1].text, value) || value < 0.0 ||
          value > std::numeric_limits<float>::max()) {
        diagnostic = Error(arguments[index + 1].span,
                           "accent requires a non-negative finite multiplier");
        return false;
      }
      parsed.accentMultiplier = static_cast<float>(value);
      accentSeen = true;
      index += 2;
      continue;
    }

    if (optionsStarted) {
      diagnostic =
          Error(argument.span,
                "envelope times and sustain precede its named options");
      return false;
    }
    const std::size_t maximumPositionals =
        parsed.mode == CvEnvelopeMode::Adsr ? 4 : 2;
    if (positionalCount >= maximumPositionals) {
      diagnostic = Error(argument.span, "too many envelope parameters");
      return false;
    }
    if (argument.text != ".") {
      if (parsed.mode == CvEnvelopeMode::Adsr && positionalCount == 2) {
        double sustain = 0.0;
        if (!ParseNumber(argument.text, sustain) || sustain < 0.0 ||
            sustain > 1.0) {
          diagnostic =
              Error(argument.span, "ADSR sustain must be between 0 and 1");
          return false;
        }
        parsed.sustain = static_cast<float>(sustain);
      } else {
        CvEnvelopeTime *time = &parsed.attack;
        if (parsed.mode == CvEnvelopeMode::Ad)
          time = positionalCount == 0 ? &parsed.attack : &parsed.decay;
        else if (parsed.mode == CvEnvelopeMode::Ar)
          time = positionalCount == 0 ? &parsed.attack : &parsed.release;
        else if (positionalCount == 1)
          time = &parsed.decay;
        else if (positionalCount == 3)
          time = &parsed.release;
        if (!ParseEnvelopeTime(argument, *time, diagnostic))
          return false;
      }
    }
    ++positionalCount;
    ++index;
  }
  spec = parsed;
  return true;
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
    if ((target == CursorLane::Notes || target == CursorLane::Rhythm) &&
        phaseTransform) {
      diagnostic = Error(pipeline.operation.span,
                         "rate is not supported on notes or rhythm; use "
                         "brackets for local density or fast/slow for the "
                         "complete sequence");
      return false;
    }
    transforms.push_back(transform);
  }
  return true;
}

bool ParseCvPipelines(const std::vector<syntax::Pipeline> &pipelines,
                      Sequence &sequence, std::size_t cvIndex,
                      CursorLane target, Diagnostic &diagnostic) {
  for (std::size_t pipelineIndex = 0; pipelineIndex < pipelines.size();
       ++pipelineIndex) {
    const auto &pipeline = pipelines[pipelineIndex];
    if (pipeline.operation.text == "add") {
      if (sequence.cvEnvelope[cvIndex].enabled) {
        diagnostic = Error(pipeline.operation.span,
                           "a CV lane can contain one envelope");
        return false;
      }
      if (pipeline.condition != syntax::Pipeline::Condition::Always ||
          pipeline.arguments.empty() ||
          pipeline.arguments.front().text != "env") {
        diagnostic = Error(pipeline.operation.span,
                           "add requires an unconditional env expression");
        return false;
      }
      if (pipelineIndex + 1 != pipelines.size()) {
        diagnostic = Error(pipelines[pipelineIndex + 1].operation.span,
                           "add env must be the final CV pipeline operation");
        return false;
      }
      std::vector<Token> arguments(pipeline.arguments.begin() + 1,
                                   pipeline.arguments.end());
      if (!ParseCvEnvelope(arguments, CvEnvelopeComposition::Add, pipeline.span,
                           sequence.cvEnvelope[cvIndex], diagnostic))
        return false;
      continue;
    }
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
      diagnostic = Error(
          mode.span, "interp mode must be step, linear, smooth, or power P");
      return false;
    }
  }
  return true;
}

bool ValidateRomanCardinality(const Sequence &sequence,
                              Diagnostic &diagnostic) {
  auto validate = [&](const PitchItem &item) {
    for (const auto &choice : item.values) {
      if (choice.meaning != ChordValue::Meaning::RomanSymbol)
        continue;
      if (choice.romanRoot.degree < 1 ||
          static_cast<std::size_t>(choice.romanRoot.degree) >
              sequence.scale.intervals.size()) {
        diagnostic =
            Error(choice.romanRoot.span,
                  "Roman chord degree exceeds the active scale cardinality");
        return false;
      }
    }
    return true;
  };
  for (const auto &item : sequence.notes) {
    if (!validate(item))
      return false;
  }
  for (const auto &step : sequence.pitchTimeline)
    if (step.kind == PitchTimelineStep::Kind::Pitch && !validate(step.pitch))
      return false;
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
           {CursorLane::Notes, CursorLane::Rhythm, CursorLane::Octave,
            CursorLane::Velocity, CursorLane::Duration, CursorLane::Gate,
            CursorLane::Slide, CursorLane::Ratchet, CursorLane::Offset,
            CursorLane::Cv1, CursorLane::Cv2, CursorLane::Cv3})
        sequence.transforms[static_cast<std::size_t>(lane)].push_back(
            transform);
    }
  }
  const auto &pitchTransforms =
      sequence.transforms[static_cast<std::size_t>(CursorLane::Notes)];
  const auto keyTranspose =
      std::find_if(pitchTransforms.begin(), pitchTransforms.end(),
                   [](const Transform &transform) {
                     return transform.kind == TransformKind::TransposeKey;
                   });
  if (keyTranspose != pitchTransforms.end() && !sequence.hasKey) {
    diagnostic = Error(keyTranspose->span,
                       "transpose_key requires a written key in the sequence");
    return false;
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
    maximumEarlyMilliseconds =
        std::max(maximumEarlyMilliseconds,
                 sequenceEarlyMilliseconds + offsetEarlyMilliseconds);
    maximumLateBeats =
        std::max(maximumLateBeats, sequenceLateBeats + offsetLateBeats);
    maximumLateMilliseconds =
        std::max(maximumLateMilliseconds,
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
    auto includeVoices = [&](const PitchItem &item) {
      for (const auto &choice : item.values) {
        maximumVoices = std::max(
            maximumVoices,
            std::min<std::size_t>(
                MaximumPolyphony,
                (choice.meaning == ChordValue::Meaning::JazzSymbol ||
                         choice.meaning == ChordValue::Meaning::RomanSymbol
                     ? MaximumVoicingCount(choice, sequence.voicing)
                     : choice.voices.size()) +
                    (choice.hasBass ? 1u : 0u)));
      }
    };
    for (const auto &item : sequence.notes)
      includeVoices(item);
    for (const auto &step : sequence.pitchTimeline)
      if (step.kind == PitchTimelineStep::Kind::Pitch)
        includeVoices(step.pitch);
    if (!CheckedMultiply(sequenceMaximum, maximumVoices, sequenceMaximum)) {
      diagnostic = Error(sequence.nameSpan, "chord event density is too large");
      return false;
    }
    maximumEvents = std::max(maximumEvents, sequenceMaximum);

    const double defaultDuration = sequence.separateRhythm
                                       ? sequence.rhythmSubdivisionBeats
                                       : sequence.subdivisionBeats;
    double shortestDurationLaneValue = defaultDuration;
    for (const auto &item : sequence.duration) {
      const double value =
          item.isDefault ? defaultDuration : PossibleScalarRange(item).minimum;
      if (value > 0.0)
        shortestDurationLaneValue = std::min(shortestDurationLaneValue, value);
    }
    double shortestStepMultiplier = 1.0;
    for (const auto &step : sequence.articulation) {
      shortestStepMultiplier =
          std::min(shortestStepMultiplier, step.durationMultiplier);
    }
    const double shortestPreparedStep =
        shortestDurationLaneValue * shortestStepMultiplier * fastestScale;
    if (!std::isfinite(shortestPreparedStep) || shortestPreparedStep <= 0.0) {
      diagnostic = Error(
          sequence.nameSpan,
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
  struct ReusableRhythm {
    double subdivisionBeats = 1.0;
    std::vector<ArticulationStep> articulation;
    std::vector<Transform> transforms;
    SourceSpan span;
  };
  std::unordered_map<std::string, ReusableRhythm> rhythms;
  std::string playName;
  SourceSpan playSpan;

  std::unordered_map<std::string, SourceSpan> definitionNames;
  for (const auto &statement : document.statements) {
    std::string name;
    SourceSpan span;
    if (const auto *sequence =
            std::get_if<syntax::SequenceDefinition>(&statement)) {
      name = sequence->name.text;
      span = sequence->name.span;
    } else if (const auto *rhythm =
                   std::get_if<syntax::RhythmDefinition>(&statement)) {
      name = rhythm->name.text;
      span = rhythm->name.span;
    } else if (const auto *assignment =
                   std::get_if<syntax::Assignment>(&statement)) {
      name = assignment->name.text;
      span = assignment->name.span;
    }
    if (name.empty())
      continue;
    if (!definitionNames.emplace(name, span).second) {
      result.diagnostic = Error(span, "duplicate definition of '" + name + "'");
      return result;
    }
  }

  for (const auto &statement : document.statements) {
    const auto *definition = std::get_if<syntax::RhythmDefinition>(&statement);
    if (!definition)
      continue;
    if (IsReservedRhythmSyntaxName(definition->name.text)) {
      result.diagnostic = Error(
          definition->name.span,
          "'" + definition->name.text +
              "' is reserved because it is complete inline rhythm syntax");
      return result;
    }
    if (IsReservedCvName(definition->name.text)) {
      result.diagnostic = Error(definition->name.span,
                                definition->name.text +
                                    " is reserved for CV lanes and signals");
      return result;
    }
    ReusableRhythm rhythm;
    rhythm.span = definition->span;
    if (!definition->subdivision.text.empty()) {
      ValueUnit unit = ValueUnit::Plain;
      if (!ParseScalarWithUnit(definition->subdivision, rhythm.subdivisionBeats,
                               unit) ||
          unit != ValueUnit::NoteValue || rhythm.subdivisionBeats <= 0.0) {
        result.diagnostic = Error(
            definition->subdivision.span,
            "rhythm subdiv requires a note value such as 16n, 8nt, or 8nd");
        return result;
      }
    }
    Sequence parsed;
    if (definition->events.steps.empty()) {
      result.diagnostic =
          Error(definition->span, "a rhythm definition requires events");
      return result;
    }
    if (!ParseRhythm(definition->events, parsed, result.diagnostic)) {
      result.diagnostic.message = "in rhythm '" + definition->name.text +
                                  "': " + result.diagnostic.message;
      return result;
    }
    rhythm.articulation = std::move(parsed.articulation);
    if (!ParsePipelines(definition->pipelines, rhythm.transforms,
                        CursorLane::Rhythm, result.diagnostic))
      return result;
    rhythms.emplace(definition->name.text, std::move(rhythm));
  }

  for (const auto &statement : document.statements) {
    if (const auto *play = std::get_if<syntax::PlayCommand>(&statement)) {
      playName = play->name.text;
      playSpan = play->name.span;
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
        result.diagnostic =
            Error(assignment->name.span,
                  assignment->name.text +
                      " is reserved for a future typed signal assignment");
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
      result.diagnostic =
          Error(definition->name.span,
                sequence.name + " is reserved for CV lanes and signal values");
      return result;
    }
    sequence.stableId = StableDefinitionId(sequence.name);
    sequence.nameSpan = definition->name.span;
    std::array<SourceSpan, static_cast<std::size_t>(CursorLane::Count)>
        alignmentSpans{};
    std::unordered_set<std::string> seenLanes;
    const bool hasSeparateRhythm = std::any_of(
        definition->lanes.begin(), definition->lanes.end(),
        [](const syntax::Lane &lane) { return lane.name.text == "rhythm"; });
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
        result.diagnostic = Error(lane.name.span, "duplicate sequence lane '" +
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

      if (isCv && lane.envelopeOnly) {
        if (!lane.pipelines.empty()) {
          result.diagnostic = Error(
              lane.pipelines.front().span,
              "an envelope-only CV lane does not accept pipeline operations");
          return result;
        }
        if (!ParseCvEnvelope(lane.envelopeArguments,
                             CvEnvelopeComposition::Replace, lane.envelopeSpan,
                             sequence.cvEnvelope[cvIndex], result.diagnostic))
          return result;
        continue;
      }

      if (laneSpec->valueKind == LaneValueKind::Subdivision) {
        double subdivisionBeats = 0.0;
        ValueUnit unit = ValueUnit::Plain;
        const Token token{laneValue, valueSpan};
        if (!ParseScalarWithUnit(token, subdivisionBeats, unit) ||
            unit != ValueUnit::NoteValue || subdivisionBeats <= 0.0) {
          result.diagnostic =
              Error(valueSpan,
                    "subdiv requires a note value such as 16n, 8nt, or 8nd");
          return result;
        }
        sequence.subdivisionBeats = subdivisionBeats;
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
      if (laneSpec->valueKind == LaneValueKind::Key) {
        std::size_t cursor = 0;
        int pitchClass = 0;
        if (!ParseNoteName(laneValue, cursor, pitchClass) ||
            cursor != laneValue.size()) {
          result.diagnostic =
              Error(valueSpan, "key must be a note name such as C, F#, or Bb");
          return result;
        }
        sequence.hasKey = true;
        sequence.keyPitchClass = pitchClass;
        continue;
      }
      if (laneSpec->valueKind == LaneValueKind::Voicing) {
        if (laneValue == "basic")
          sequence.voicing = VoicingStyle::Basic;
        else if (laneValue == "rootless_3notes")
          sequence.voicing = VoicingStyle::Rootless3Notes;
        else if (laneValue == "rootless_4notes")
          sequence.voicing = VoicingStyle::Rootless4Notes;
        else {
          result.diagnostic = Error(
              valueSpan,
              "voicing must be basic, rootless_3notes, or rootless_4notes");
          return result;
        }
        continue;
      }
      if (laneSpec->valueKind == LaneValueKind::Glide) {
        double glide = 0.0;
        ValueUnit unit = ValueUnit::Plain;
        const Token token{laneValue, valueSpan};
        if (!ParseScalarWithUnit(token, glide, unit) ||
            unit == ValueUnit::Milliseconds || glide < 0.0 ||
            glide > std::numeric_limits<float>::max()) {
          result.diagnostic = Error(
              valueSpan, "glide must be a non-negative beat or note value");
          return result;
        }
        sequence.glideBeats = static_cast<float>(glide);
        continue;
      }

      bool ok = true;
      if (laneSpec->valueKind == LaneValueKind::Notes ||
          laneSpec->valueKind == LaneValueKind::Chords) {
        ok = hasSeparateRhythm
                 ? ParsePitchTimeline(lane.pattern, sequence, result.diagnostic)
                 : ParseNotes(lane.pattern, sequence, result.diagnostic);
        if (!ok && hasSeparateRhythm)
          result.diagnostic.message =
              "in held " + laneName + " lane: " + result.diagnostic.message;
        if (ok && laneSpec->valueKind == LaneValueKind::Chords)
          PromoteBareNamedChordLaneValues(sequence);
      } else if (laneSpec->valueKind == LaneValueKind::Rhythm) {
        if (!lane.rhythmReference.text.empty()) {
          const auto found = rhythms.find(lane.rhythmReference.text);
          if (found == rhythms.end()) {
            result.diagnostic =
                Error(lane.rhythmReference.span,
                      "unknown rhythm '" + lane.rhythmReference.text + "'");
            return result;
          }
          sequence.articulation = found->second.articulation;
          sequence.rhythmSubdivisionBeats = found->second.subdivisionBeats;
          sequence.transforms[static_cast<std::size_t>(CursorLane::Rhythm)] =
              found->second.transforms;
          sequence.separateRhythm = true;
        } else {
          ok = ParseRhythm(lane.pattern, sequence, result.diagnostic);
          // Inline rhythm retains the sequence's subdivision, including when
          // the setting appears later in the source block.
          sequence.rhythmSubdivisionBeats = 0.0;
        }
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
          const auto laneIndex = static_cast<std::size_t>(parsedLane);
          if (lane.pattern.alignment == syntax::Pattern::Alignment::Left) {
            sequence.alignment[laneIndex] = LaneAlignment::Left;
          } else if (lane.pattern.alignment ==
                     syntax::Pattern::Alignment::Right) {
            sequence.alignment[laneIndex] = LaneAlignment::Right;
          } else {
            sequence.alignment[laneIndex] = LaneAlignment::Edges;
            std::size_t split = 0;
            for (std::size_t term = 0; term < lane.pattern.alignmentSplit;
                 ++term) {
              std::size_t repetitions = 1;
              if (!ParsePositiveCount(lane.pattern.steps[term].repeatCount,
                                      laneName + " repetition", repetitions,
                                      result.diagnostic)) {
                ok = false;
                break;
              }
              if (split >
                  std::numeric_limits<std::size_t>::max() - repetitions) {
                result.diagnostic =
                    Error(lane.pattern.steps[term].span,
                          laneName + " edge alignment has too many values");
                ok = false;
                break;
              }
              split += repetitions;
            }
            sequence.alignmentSplits[laneIndex] = split;
          }
          sequence.alignmentSpans[static_cast<std::size_t>(parsedLane)] =
              lane.pattern.span;
          alignmentSpans[static_cast<std::size_t>(parsedLane)] =
              lane.pattern.span;
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
          Error(sequence.nameSpan,
                hasSeparateRhythm ? "sequence requires a non-empty rhythm lane"
                                  : "sequence requires a notes lane");
      return result;
    }
    if (hasSeparateRhythm && sequence.pitchTimeline.empty()) {
      result.diagnostic =
          Error(sequence.nameSpan,
                "a sequence with rhythm requires a notes or chords lane");
      return result;
    }
    if (hasSeparateRhythm) {
      sequence.pitchTimelineBeats = 0.0;
      for (const auto &value : sequence.pitchTimeline)
        sequence.pitchTimelineBeats +=
            sequence.subdivisionBeats * value.durationMultiplier;
      if (sequence.rhythmSubdivisionBeats <= 0.0)
        sequence.rhythmSubdivisionBeats = sequence.subdivisionBeats;
      const bool hasHit = std::any_of(
          sequence.articulation.begin(), sequence.articulation.end(),
          [](const ArticulationStep &step) {
            return std::any_of(step.atoms.begin(), step.atoms.end(),
                               [](const ArticulationAtom &atom) {
                                 return atom.kind == ArticulationKind::Attack ||
                                        atom.kind == ArticulationKind::Slide;
                               });
          });
      if (!hasHit) {
        result.diagnostic = Error(sequence.nameSpan,
                                  "a rhythm gesture requires at least one hit");
        return result;
      }
      const auto targetSlideTime = std::find_if(
          sequence.pitchTimeline.begin(), sequence.pitchTimeline.end(),
          [](const PitchTimelineStep &value) { return value.hasSlide; });
      const bool gestureSlideTime = std::any_of(
          sequence.articulation.begin(), sequence.articulation.end(),
          [](const ArticulationStep &step) {
            return std::any_of(
                step.atoms.begin(), step.atoms.end(),
                [](const ArticulationAtom &atom) { return atom.hasSlide; });
          });
      if (targetSlideTime != sequence.pitchTimeline.end() && gestureSlideTime) {
        result.diagnostic =
            Error(targetSlideTime->span,
                  "slide time is specified in both the held pitch timeline and "
                  "the rhythm gesture; keep one explicit value");
        return result;
      }
    }
    if (!ValidateRomanCardinality(sequence, result.diagnostic))
      return result;
    if (!sequence.ratchet.empty()) {
      for (const auto &step : sequence.articulation) {
        const auto duplicate = std::find_if(
            step.atoms.begin(), step.atoms.end(),
            [](const ArticulationAtom &atom) { return atom.ratchets > 1; });
        if (duplicate != step.atoms.end()) {
          result.diagnostic =
              Error(duplicate->span,
                    "inline *N ratchet cannot be combined with a ratchet lane");
          return result;
        }
      }
    }
    std::size_t structuralCells = 0;
    for (const auto &step : sequence.articulation) {
      if (structuralCells >
          std::numeric_limits<std::size_t>::max() - step.cellCount) {
        result.diagnostic =
            Error(step.span, "note pattern has too many structural cells");
        return result;
      }
      structuralCells += step.cellCount;
    }
    auto alignLane = [&](CursorLane lane, std::vector<ScalarItem> &items) {
      const auto index = static_cast<std::size_t>(lane);
      if (sequence.alignment[index] == LaneAlignment::Free)
        return true;
      if (items.size() > structuralCells) {
        result.diagnostic = Error(alignmentSpans[index],
                                  "aligned lane has more values than notes");
        return false;
      }
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
    for (const auto lane :
         {CursorLane::Octave, CursorLane::Velocity, CursorLane::Duration,
          CursorLane::Gate, CursorLane::Slide, CursorLane::Ratchet,
          CursorLane::Offset, CursorLane::Cv1, CursorLane::Cv2,
          CursorLane::Cv3}) {
      if (sequence.alignment[static_cast<std::size_t>(lane)] ==
          LaneAlignment::Free)
        continue;
      if (IsCvCursorLane(lane) &&
          sequence.cvInterpolation[CvCursorIndex(lane)] !=
              CvInterpolation::Step) {
        result.diagnostic =
            Error(alignmentSpans[static_cast<std::size_t>(lane)],
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
        result.diagnostic =
            Error(alignmentSpans[static_cast<std::size_t>(lane)],
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
    AssignRandomIdentities(sequence);
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
