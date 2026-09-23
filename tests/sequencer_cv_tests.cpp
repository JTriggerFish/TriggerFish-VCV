#include "tfseq_cv.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {
int failures = 0;

void check(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
  }
}

std::string program(const std::string &body) {
  return "a = sequence {\n" + body + "\n}\nplay a\n";
}

void cvPhaseIsIndependentOfNoteDensity() {
  for (const std::string notes : {"subdiv 1n\nnotes 1", "subdiv 16n\nnotes 1",
                                  "subdiv 4n\nnotes 1 ~ 2 _"}) {
    for (const std::string mode : {"step", "linear", "smooth", "power 2"}) {
      for (const int rate : {1, 4}) {
        const auto compiled =
            tfseq::Compile(program(notes + "\ncv1 0 5 |> interp " + mode +
                                   " |> rate " + std::to_string(rate)));
        check(bool(compiled),
              "CV fixture compiles: " + compiled.diagnostic.message);
        if (!compiled)
          continue;
        tfseq::Runtime runtime;
        runtime.setProgram(compiled.program.get());
        tfseq::CvLanePlayer player;
        double nextBeat = 0.0;
        bool matches = true;
        for (int frame = 0; frame < 512; ++frame) {
          const double beat = frame / 64.0;
          if (beat >= nextBeat) {
            const auto events = runtime.next(nextBeat);
            if (events.count == 1)
              player.setEvent(events.events[0], 0);
            else
              matches = false;
            nextBeat += events.durationBeats;
          }
          const double phase = std::fmod(beat * rate, 2.0);
          const bool falling = phase >= 1.0;
          double amount = phase - std::floor(phase);
          if (mode == "step")
            amount = 0.0;
          else if (mode == "smooth")
            amount = amount * amount * (3.0 - 2.0 * amount);
          else if (mode == "power 2")
            amount *= amount;
          const float expected =
              static_cast<float>(falling ? 5.0 * (1.0 - amount) : 5.0 * amount);
          matches &= std::abs(player.process(beat) - expected) < 1.e-4f;
        }
        check(matches,
              "CV follows score time through long/short notes and rests: " +
                  notes + " / " + mode + " / rate " + std::to_string(rate));
      }
    }
  }
}

void defaultsAlignmentTimingAndLifetime() {
  auto compiled = tfseq::Compile(
      program("subdiv 1n\nnotes 1\ncv1 0 . 4 . |> interp linear\n"
              "cv2 . . . |> interp smooth\ncv3 . . 3"));
  check(bool(compiled), "default CV knots compile");
  if (!compiled)
    return;
  tfseq::Runtime runtime;
  runtime.setProgram(compiled.program.get());
  auto event = runtime.next(0.0).events[0];
  tfseq::CvLanePlayer linear, empty, held;
  linear.setEvent(event, 0);
  empty.setEvent(event, 1);
  held.setEvent(event, 2);
  check(std::abs(linear.process(1.5) - 3.f) < 1.e-5f &&
            std::abs(linear.process(3.5) - 1.f) < 1.e-5f,
        "default knots are skipped in both halves of an independent CV loop");
  check(empty.process(100.0) == 0.f && held.process(100.0) == 3.f,
        "all-default CV and leading defaults remain stable");

  // Mirror the scheduler's absolute-beat and millisecond-offset conversion.
  event.beat += 10.25;
  event.cvOriginBeat += 10.25;
  linear.setEvent(event, 0);
  check(std::abs(linear.process(11.75) - 3.f) < 1.e-5f,
        "a scheduled timing shift moves the entire CV curve");
  const float last = linear.process(12.0);
  linear.detach();
  empty.detach();
  held.detach();
  compiled.program.reset();
  check(linear.process(100.0) == last,
        "detaching for live edit preserves output without retaining retired "
        "source");

  auto aligned = tfseq::Compile(program("notes 1 2 3\ncv1 2 ... 4"));
  check(bool(aligned), "aligned CV fixture compiles");
  if (aligned) {
    runtime.setProgram(aligned.program.get());
    for (int step = 0; step < 3; ++step) {
      const auto events = runtime.next(step);
      linear.setEvent(events.events[0], 0);
      const float expected = step == 0 ? 2.f : step == 2 ? 4.f : 0.f;
      check(linear.process(step + .75) == expected,
            "aligned lanes retain their structural sample-and-hold behavior");
    }
  }
}

void oversizedCvValuesAreRejected() {
  const std::string huge = "1" + std::string(40, '0');
  for (const std::string value :
       {huge, "-" + huge, "$u{0," + huge + "}", "$n{0," + huge + "}"}) {
    check(!tfseq::Compile(program("notes 1\ncv1 " + value)),
          "CV values outside finite output precision are rejected");
  }
  check(
      !tfseq::Compile(program("notes 1\ngate " + huge + "ms")) &&
          !tfseq::Compile(program("notes 1\nslide " + huge)) &&
          !tfseq::Compile(program("notes 1\ncv1 0 5 |> interp power " + huge)),
      "time and interpolation values cannot overflow output precision");

  const std::string extremeRate = "1" + std::string(308, '0');
  auto extreme = tfseq::Compile(
      program("notes 1\nvelocity .5 |> rate " + extremeRate +
              "\ncv1 $u{-1,1} 2 |> interp linear |> rate " + extremeRate));
  check(bool(extreme), "finite extreme phase-rate fixture compiles");
  if (extreme) {
    tfseq::Runtime runtime;
    runtime.setProgram(extreme.program.get());
    tfseq::CvLanePlayer player;
    for (int beat = 0; beat < 4; ++beat) {
      const auto event = runtime.next(beat).events[0];
      player.setEvent(event, 0);
      check(std::isfinite(event.cvValue[0]) &&
                std::isfinite(player.process(beat + .5)),
            "extreme phase arithmetic stays finite and safely representable");
    }
  }
}

void verySlowCvRatesKeepFiniteInterpolation() {
  const std::string tinyRate = "." + std::string(307, '0') + "1";
  for (const std::string mode : {"linear", "smooth", "power 2"}) {
    for (const std::string knots : {". . 4 0", ". . 0 4 ."}) {
      auto compiled = tfseq::Compile(program("subdiv 1n\nnotes 1\ncv1 " +
                                             knots + " |> interp " + mode +
                                             " |> rate " + tinyRate));
      check(bool(compiled), "very slow CV fixture compiles");
      if (!compiled)
        continue;
      tfseq::Runtime runtime;
      runtime.setProgram(compiled.program.get());
      const auto event = runtime.next(0.0).events[0];
      tfseq::CvLanePlayer player;
      player.setEvent(event, 0);
      check(std::isfinite(event.cvValue[0]) &&
                std::abs(player.process(0.0) - event.cvValue[0]) < 1.e-5f &&
                std::abs(player.process(1.0) - event.cvValue[0]) < 1.e-5f,
            "slow curves retain their value when knot times overflow");
      // No new note event: rendering must still advance in the lane's phase.
      const double largeBeat = 5.e307;
      tfseq::Runtime reference;
      reference.setProgram(compiled.program.get());
      reference.next(0.0);
      const float expected = reference.next(largeBeat).events[0].cvValue[0];
      check(std::isfinite(expected) &&
                std::abs(player.process(largeBeat) - expected) < 1.e-5f,
            "slow cached interpolation advances without infinite-time ratios");
    }
  }
}
} // namespace

int main() {
  cvPhaseIsIndependentOfNoteDensity();
  defaultsAlignmentTimingAndLifetime();
  oversizedCvValuesAreRejected();
  verySlowCvRatesKeepFiniteInterpolation();
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
