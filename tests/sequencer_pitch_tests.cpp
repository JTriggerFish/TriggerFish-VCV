#include "tfseq.hpp"
#include "tfseq_cv.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <new>
#include <string>

namespace {
thread_local bool trackAllocations = false;
thread_local std::size_t allocations = 0;
int failures = 0;
} // namespace

void *operator new(std::size_t size) {
  if (trackAllocations)
    ++allocations;
  if (void *memory = std::malloc(size ? size : 1))
    return memory;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t size) { return ::operator new(size); }
void *operator new(std::size_t size, const std::nothrow_t &) noexcept {
  try {
    return ::operator new(size);
  } catch (...) {
    return nullptr;
  }
}
void *operator new[](std::size_t size, const std::nothrow_t &tag) noexcept {
  return ::operator new(size, tag);
}
void operator delete(void *memory) noexcept { std::free(memory); }
void operator delete[](void *memory) noexcept { std::free(memory); }
void operator delete(void *memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void *memory, std::size_t) noexcept {
  std::free(memory);
}
void operator delete(void *memory, const std::nothrow_t &) noexcept {
  std::free(memory);
}
void operator delete[](void *memory, const std::nothrow_t &) noexcept {
  std::free(memory);
}

namespace {
void check(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

std::string program(const std::string &body) {
  return "changes = sequence {\n" + body + "\n}\nplay changes\n";
}

void expectFirstChord(const std::string &body,
                      const std::array<int, 3> &expected) {
  const auto compiled = tfseq::Compile(program(body));
  check(bool(compiled), "compile: " + compiled.diagnostic.message);
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  const auto events = runtime.next(0.0);
  bool matches = events.count == expected.size();
  for (std::size_t voice = 0; voice < events.count && voice < expected.size();
       ++voice)
    matches &= std::abs(events.events[voice].pitchVolts * 12.f -
                        expected[voice]) < 1.e-4f;
  check(matches, "written key/register: " + body);
}

void keyRootsAndSettingOrder() {
  // GitHub #13, including named-note and numeric-degree equivalents.
  for (const std::string lane :
       {"chords i", "notes (1 3 5)", "chords Dm", "notes (D F A)"})
    expectFirstChord("subdiv 2n\nkey D\nscale dorian\nvoicing basic\n" + lane,
                     {2, 5, 9});

  const std::array<const char *, 12> keys{"C",  "Db", "D",  "Eb", "E",  "F",
                                          "F#", "G",  "Ab", "A",  "Bb", "B"};
  struct Mode {
    const char *name;
    int third;
    int fifth;
  };
  for (const auto mode :
       {Mode{"major", 4, 7}, Mode{"dorian", 3, 7}, Mode{"phrygian", 3, 7},
        Mode{"lydian", 4, 7}, Mode{"mixolydian", 4, 7}, Mode{"minor", 3, 7},
        Mode{"locrian", 3, 6}}) {
    for (int root = 0; root < 12; ++root) {
      const std::string key = std::string("key ") + keys[root] + "\n";
      const std::string scale = std::string("scale ") + mode.name + "\n";
      expectFirstChord(key + scale + "notes (1 3 5)",
                       {root, root + mode.third, root + mode.fifth});
      expectFirstChord(
          "chords III\n" + scale + key,
          {root + mode.third, root + mode.third + 4, root + mode.third + 7});
    }
  }
  // An explicitly written tonic remains authoritative in either order.
  for (const std::string settings :
       {"key D\ntonic E@3\n", "tonic E@3\nkey D\n"})
    expectFirstChord(settings + "scale minor\nchords i", {-8, -5, -1});
  expectFirstChord("key D\nscale dorian\nchords i@3", {-10, -7, -3});
  expectFirstChord("key D\nscale dorian\nchords i'", {14, 17, 21});
  expectFirstChord("key D\nscale dorian\nchords bVII", {11, 15, 18});
  expectFirstChord("key F#\nscale minor\nchords C", {0, 4, 7});
  expectFirstChord("scale dorian\nchords i", {0, 3, 7});

  auto compiled =
      tfseq::Compile("base = sequence {\nkey D\nscale dorian\nchords i\n}\n"
                     "shifted = base |> transpose_key E\nplay shifted\n");
  check(bool(compiled), "keyed derived sequence compiles");
  if (compiled) {
    tfseq::Runtime runtime;
    runtime.setProgram(compiled.program.get());
    const auto events = runtime.next(0.0);
    check(events.count == 3 &&
              std::abs(events.events[0].pitchVolts * 12.f - 4.f) < 1.e-4f,
          "transpose_key applies the key delta exactly once to Roman roots");
  }
}

void progressionsStayInRegister() {
  struct Progression {
    const char *chords;
    std::array<int, 3> roots;
  };
  for (const auto progression :
       {Progression{"I@4 IV@4 V@4", {0, 5, 7}},
        Progression{"Dm@4 C@5 G@4", {2, 12, 7}},
        Progression{"I IV V", {0, 5, 7}},
        Progression{"Dm9 G13 Cmaj9", {2, 7, 0}},
        Progression{"C7alt F13 Bbmaj7", {0, 5, 10}},
        Progression{"Cmaj:(3) Fmaj:(3) Gmaj:(3)", {0, 5, 7}}}) {
    for (const std::string style :
         {"basic", "rootless_3notes", "rootless_4notes"}) {
      for (const bool separateRhythm : {false, true}) {
        const auto compiled = tfseq::Compile(program(
            "subdiv 2n\nkey C\nvoicing " + style + "\nchords " +
            progression.chords + (separateRhythm ? "\nrhythm x x x" : "")));
        check(bool(compiled),
              "progression compiles: " + compiled.diagnostic.message);
        if (!compiled)
          continue;
        tfseq::Runtime runtime;
        runtime.setProgram(compiled.program.get());
        bool bounded = true;
        bool sounding = true;
        bool ordered = true;
        allocations = 0;
        trackAllocations = true;
        for (int step = 0; step < 300; ++step) {
          const auto events = runtime.next(step * 2.0);
          sounding &= events.count > 0;
          if (events.count == 0)
            continue;
          const float low = events.events[0].pitchVolts * 12.f;
          const float high = events.events[events.count - 1].pitchVolts * 12.f;
          const int target =
              progression.roots[step % 3] + (style == "basic" ? 6 : -2);
          bounded &= std::abs((low + high) * .5f - target) <= 6.001f;
          for (std::size_t voice = 1; voice < events.count; ++voice)
            ordered &= events.events[voice - 1].pitchVolts <
                       events.events[voice].pitchVolts;
        }
        trackAllocations = false;
        const std::string label = style + ": " + progression.chords +
                                  (separateRhythm ? " with rhythm" : "");
        check(sounding && ordered, "ordered sounding progression: " + label);
        check(bounded, "100 loops remain in the written register: " + label);
        check(allocations == 0,
              "playback performs no heap allocations: " + label);
      }
    }
  }
}

void registerChangesAndPlaybackState() {
  for (const std::string style :
       {"basic", "rootless_3notes", "rootless_4notes"}) {
    struct Case {
      const char *body;
      std::array<int, 3> roots;
    };
    for (const auto fixture : {Case{"chords C7@4 C7@6 C7@2", {0, 24, -24}},
                               Case{"chords C7 C7' C7,", {0, 12, -12}},
                               Case{"octave 4 6 2\nchords C7", {0, 24, -24}},
                               Case{"chords C7 |> octave 2", {24, 24, 24}}}) {
      auto compiled =
          tfseq::Compile(program("voicing " + style + "\n" + fixture.body));
      check(bool(compiled), "register-change fixture compiles");
      if (!compiled)
        continue;
      tfseq::Runtime runtime;
      runtime.setProgram(compiled.program.get());
      bool followsRegister = true;
      for (int step = 0; step < 30; ++step) {
        const auto events = runtime.next(step);
        followsRegister &= events.count > 0;
        if (events.count == 0)
          continue;
        const float midpoint =
            6.f * (events.events[0].pitchVolts +
                   events.events[events.count - 1].pitchVolts);
        const int target =
            fixture.roots[step % 3] + (style == "basic" ? 6 : -2);
        followsRegister &= std::abs(midpoint - target) <= 6.001f;
      }
      check(followsRegister,
            "written register changes override prior context: " + style);
    }
  }

  auto initial = tfseq::Compile(program("key D\nscale dorian\nchords i"));
  auto edited = tfseq::Compile(program("key E\nscale dorian\nchords i"));
  check(bool(initial) && bool(edited), "live-key fixtures compile");
  if (initial && edited) {
    tfseq::Runtime runtime;
    runtime.setProgram(initial.program.get());
    runtime.next(0.0);
    runtime.replaceProgram(edited.program.get(), 1.0);
    const auto changed = runtime.next(1.0);
    check(changed.count == 3 &&
              std::abs(changed.events[0].pitchVolts * 12.f - 4.f) < 1.e-4f,
          "a live key edit starts the new harmony in its written register");
    runtime.reset();
    const auto restarted = runtime.next(0.0);
    check(restarted.count == 3 &&
              std::abs(restarted.events[0].pitchVolts * 12.f - 4.f) < 1.e-4f,
          "reset retains the edited key");
  }

  expectFirstChord("key D\nscale dorian\nchords (D@3 F@5 A@4)", {-10, 17, 9});
  auto slash = tfseq::Compile(program("key D\nscale dorian\nchords i/A@2"));
  check(bool(slash), "keyed slash chord compiles");
  if (slash) {
    tfseq::Runtime runtime;
    runtime.setProgram(slash.program.get());
    const auto events = runtime.next(0.0);
    check(events.count == 3 &&
              std::abs(events.events[0].pitchVolts * 12.f + 15.f) < 1.e-4f,
          "a slash bass keeps its explicit register and removes the duplicate "
          "tone");
  }
}

void enharmonicSpellingsKeepTheirWrittenOctave() {
  expectFirstChord("notes (B#@4 Cb@4 E#@4)", {12, -1, 5});
  expectFirstChord("tonic C@3\nnotes (B# Cb Fb)", {0, -13, -8});
  expectFirstChord("chords B#maj@4", {12, 16, 19});
  expectFirstChord("chords Cbmaj@4", {-1, 3, 6});
  expectFirstChord("chords B#@4", {12, 16, 19});
  expectFirstChord("tonic Cb@4\nchords I", {-1, 3, 6});
  expectFirstChord("tonic B#@4\nnotes (1 3 5)", {12, 16, 19});
  expectFirstChord("notes (B# Cb E#)@4", {12, -1, 5});
}

void probabilisticSlideLookaheadMatchesLaterLoops() {
  for (const std::string pattern : {"1 >2?0.5", "1 3??0 >2?0.5"}) {
    auto compiled = tfseq::Compile(program("notes " + pattern));
    check(bool(compiled), "probabilistic slide fixture compiles");
    if (!compiled)
      continue;
    tfseq::Runtime runtime;
    runtime.setProgram(compiled.program.get());
    bool matches = true;
    bool heardSlide = false;
    bool heardRest = false;
    for (int loop = 0; loop < 100; ++loop) {
      const auto first = runtime.next(loop * 2.0);
      const bool anticipatesSlide =
          first.count == 1 && first.events[0].legatoToNext;
      const auto next = runtime.next(loop * 2.0 + 1.0);
      const bool slides =
          next.count == 1 && next.events[0].kind == tfseq::EventKind::Slide;
      const bool rests =
          next.count == 1 && next.events[0].kind == tfseq::EventKind::Rest;
      matches &= anticipatesSlide == slides;
      heardSlide |= slides;
      heardRest |= rests;
    }
    check(matches && heardSlide && heardRest,
          "slide lookahead and playback make the same probability decision on "
          "every loop");
  }
}

void separateRhythmSlideLookaheadMatchesPhraseBoundaries() {
  for (const std::string notes : {"1_3", "1_4"}) {
    auto compiled =
        tfseq::Compile("pulse = rhythm {\nsubdiv 4n\nevents x >x?0.5\n}\n" +
                       program("subdiv 4n\nnotes " + notes + "\nrhythm pulse"));
    check(bool(compiled), "independent probabilistic rhythm compiles: " +
                              compiled.diagnostic.message);
    if (!compiled)
      continue;
    tfseq::Runtime runtime;
    runtime.setProgram(compiled.program.get());
    bool matches = true;
    bool previousAttack = false;
    bool predicted = false;
    bool heardSlide = false;
    bool heardRest = false;
    for (int beat = 0; beat < 400; ++beat) {
      const auto events = runtime.next(beat);
      if (events.count != 1) {
        matches = false;
        continue;
      }
      const auto &event = events.events[0];
      if (previousAttack) {
        const bool slides = event.kind == tfseq::EventKind::Slide;
        matches &= predicted == slides;
        heardSlide |= slides;
        heardRest |= event.kind == tfseq::EventKind::Rest;
      }
      previousAttack = event.kind == tfseq::EventKind::Attack;
      predicted = event.legatoToNext;
    }
    check(matches && heardSlide && heardRest,
          "rhythm slide predictions survive rhythm and note-phrase boundaries");
  }
}

void cvPlaybackDoesNotAllocate() {
  auto compiled = tfseq::Compile(program(
      "subdiv 1n\nnotes 1\ncv1 $u{-2,2} . 4 |> interp smooth |> rate 4\n"
      "cv2 1 -1 |> interp power 2\ncv3 . 3 ."));
  check(bool(compiled), "realtime CV fixture compiles");
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  std::array<tfseq::CvLanePlayer, tfseq::CvLaneCount> players;
  bool finite = true;
  allocations = 0;
  trackAllocations = true;
  for (int frame = 0; frame < 4096; ++frame) {
    const double beat = frame / 256.0;
    if (frame % 1024 == 0) {
      const auto event = runtime.next(beat).events[0];
      for (std::size_t lane = 0; lane < players.size(); ++lane)
        players[lane].setEvent(event, lane);
    }
    for (auto &player : players)
      finite &= std::isfinite(player.process(beat));
  }
  trackAllocations = false;
  check(
      finite && allocations == 0,
      "continuous random and interpolated CV playback performs no allocations");
}
} // namespace

int main() {
  keyRootsAndSettingOrder();
  progressionsStayInRegister();
  registerChangesAndPlaybackState();
  enharmonicSpellingsKeepTheirWrittenOctave();
  probabilisticSlideLookaheadMatchesLaterLoops();
  separateRhythmSlideLookaheadMatchesPhraseBoundaries();
  cvPlaybackDoesNotAllocate();
  if (failures)
    std::cerr << failures << " sequencer pitch checks failed\n";
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
