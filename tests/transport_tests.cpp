#include "tfseq_program_mailbox.hpp"
#include "tftransport.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>

namespace {
int failures = 0;

void check(bool condition, const char *message) {
  if (condition)
    return;
  std::cerr << "FAIL: " << message << '\n';
  ++failures;
}

struct alignas(4) TrackedProgram {
  int identity = 0;
  int *destructions = nullptr;

  ~TrackedProgram() {
    if (destructions)
      ++*destructions;
  }
};
} // namespace

int main() {
  {
    int destructions = 0;
    tfseq::ProgramMailbox<TrackedProgram, 3> programs;
    programs.publish(new TrackedProgram{1, &destructions}, 1);
    const auto first = programs.protect();
    check(tfseq::ProgramMailbox<TrackedProgram, 3>::pointer(first)->identity ==
              1,
          "the audio side protects the published program");

    programs.publish(new TrackedProgram{2, &destructions}, 2);
    check(destructions == 0,
          "publication defers deletion of a protected replacement");
    check(!programs.claim(first),
          "an exact claim rejects a concurrently replaced program");
    programs.collect();
    check(destructions == 1,
          "the UI reclaims a superseded program after hazard release");

    const auto second = programs.protect();
    auto *claimed = tfseq::ProgramMailbox<TrackedProgram, 3>::pointer(second);
    check((second & 3) == 2 && claimed->identity == 2 && programs.claim(second),
          "the audio side claims the exact tagged replacement");
    delete claimed;
    check(programs.pending() == 0 && destructions == 2,
          "a claimed program leaves mailbox ownership exactly once");
  }

  tftransport::CommandMailbox mailbox;
  tftransport::Command requested = tftransport::Command::Stop;
  check(!mailbox.consume(requested), "the command mailbox starts empty");
  mailbox.post(tftransport::Command::Pause);
  check(mailbox.consume(requested) && requested == tftransport::Command::Pause,
        "the audio side consumes a UI command exactly once");
  check(!mailbox.consume(requested),
        "a consumed command is removed from the mailbox");
  mailbox.post(tftransport::Command::Pause);
  mailbox.post(tftransport::Command::PlayFromBeginning);
  check(mailbox.consume(requested) &&
            requested == tftransport::Command::PlayFromBeginning,
        "the latest unconsumed UI command wins");

  tftransport::Engine transport;
  check(transport.state() == tftransport::State::Stopped,
        "transport initializes stopped");
  check(tftransport::PulsesPerQuarterNote == 24 &&
            tftransport::IsQuarterNotePulse(0) &&
            !tftransport::IsQuarterNotePulse(23) &&
            tftransport::IsQuarterNotePulse(24) &&
            std::abs(tftransport::BeatAtPulse(6) - 0.25) < 1e-12 &&
            tftransport::BeatPeriodFromPulseSamples(2000) == 48000.0,
        "sender and receiver share one fixed 24 PPQN timebase");
  check(tftransport::RescaleSampleCount(1000, 48000.0, 96000.0) == 2000 &&
            tftransport::RescaleSampleCount(1000, 96000.0, 48000.0) == 500 &&
            tftransport::RescaleSampleCount(0, 48000.0, 96000.0) == 0,
        "sample-rate changes preserve elapsed transport time");

  bool receiverPeriodKnown = false;
  double receiverPeriodSamples = 0.0;
  check(!tftransport::UpdateBeatPeriodEstimate(1250, false, receiverPeriodKnown,
                                               receiverPeriodSamples) &&
            !receiverPeriodKnown && receiverPeriodSamples == 0.0,
        "a skipped startup interval cannot publish a zero known period");
  check(tftransport::UpdateBeatPeriodEstimate(1250, true, receiverPeriodKnown,
                                              receiverPeriodSamples) &&
            receiverPeriodKnown && receiverPeriodSamples == 30000.0,
        "the next continuous 24-PPQN interval acquires the beat period");
  check(!tftransport::UpdateBeatPeriodEstimate(
            50000, false, receiverPeriodKnown, receiverPeriodSamples) &&
            receiverPeriodKnown && receiverPeriodSamples == 30000.0,
        "a discontinuity preserves an already established period");

  tftransport::PulseGrid grid;
  check(grid.advance() && grid.seen() && grid.pulse() == 0,
        "the first transport pulse establishes a quarter-note boundary");
  for (int pulse = 1; pulse < tftransport::PulsesPerQuarterNote; ++pulse)
    check(!grid.advance(), "intermediate transport pulses are not quarters");
  check(grid.advance() && grid.pulse() == tftransport::PulsesPerQuarterNote,
        "the independent pulse grid finds the next quarter-note boundary");
  grid.reset();
  check(grid.advance() && grid.pulse() == 0,
        "reset realigns the independent pulse grid to beat zero");

  tftransport::RunSynchronizedClockEdge synchronizedEdge;
  check(!synchronizedEdge.process(true, true, false) &&
            synchronizedEdge.pending(),
        "a clock edge waits while its RUN cable is one frame behind");
  check(synchronizedEdge.process(false, true, true) &&
            !synchronizedEdge.pending(),
        "the still-high clock edge is delivered when RUN catches up");
  check(!synchronizedEdge.process(false, true, true),
        "a synchronized clock edge is delivered only once");
  check(!synchronizedEdge.process(true, true, false) &&
            !synchronizedEdge.process(false, false, false) &&
            !synchronizedEdge.process(false, false, true),
        "a pulse that returned low while stopped is never replayed");
  synchronizedEdge.process(true, true, false);
  synchronizedEdge.reset();
  check(!synchronizedEdge.pending() &&
            !synchronizedEdge.process(false, true, true),
        "reset discards a pending RUN-synchronization edge");

  constexpr double PulsePeriodAt60Bpm = 1.0 / 24.0;
  constexpr double PulsePeriodStep = PulsePeriodAt60Bpm / 4.0;
  transport.command(tftransport::Command::Play);
  auto frame = transport.process(PulsePeriodStep, 60.0);
  check(frame.run && frame.clock && !frame.reset,
        "resuming a stopped transport starts pulse zero immediately");
  frame = transport.process(PulsePeriodStep, 60.0);
  check(frame.run && !frame.clock,
        "the clock waits for the remainder of its period");

  const double pausedPhase = transport.phase();
  transport.command(tftransport::Command::Pause);
  frame = transport.process(0.5, 60.0);
  check(!frame.run && !frame.clock && !frame.reset &&
            std::abs(transport.phase() - pausedPhase) < 1e-12,
        "pause lowers RUN and freezes clock phase without resetting");
  transport.command(tftransport::Command::Play);
  frame = transport.process(PulsePeriodStep, 60.0);
  check(frame.run && !frame.clock,
        "resume continues from the preserved fractional phase");
  transport.command(tftransport::Command::TogglePlayPause);
  check(transport.state() == tftransport::State::Paused,
        "play/pause toggles a playing transport to paused");
  transport.command(tftransport::Command::TogglePlayPause);
  check(transport.state() == tftransport::State::Playing,
        "play/pause toggles a paused transport back to playing");
  frame = transport.process(PulsePeriodStep, 60.0);
  check(frame.clock, "the resumed clock reaches its original next boundary");

  int subsequentPulses = 0;
  for (int step = 0; step < 24 * 4; ++step) {
    frame = transport.process(PulsePeriodStep, 60.0);
    subsequentPulses += frame.clock ? 1 : 0;
  }
  check(subsequentPulses == tftransport::PulsesPerQuarterNote,
        "the fixed master clock emits exactly 24 pulses per quarter note");

  transport.command(tftransport::Command::Stop);
  frame = transport.process(1.0 / 48000.0, 120.0);
  check(!frame.run && !frame.clock && frame.reset &&
            transport.state() == tftransport::State::Stopped &&
            transport.phase() == 0.0,
        "stop lowers RUN, pulses RESET, and returns to beat zero");

  transport.command(tftransport::Command::PlayFromBeginning);
  frame = transport.process(1.0 / 48000.0, 120.0);
  check(frame.reset && !frame.run && !frame.clock,
        "play from beginning sends RESET before RUN and CLOCK");
  bool started = false;
  for (int sample = 0; sample < 100; ++sample) {
    frame = transport.process(1.0 / 48000.0, 120.0);
    if (!frame.run)
      check(!frame.clock, "restart lead time contains no clock edge");
    if (frame.run) {
      check(frame.clock, "restart begins with a beat-zero clock");
      started = true;
      break;
    }
  }
  check(started, "restart completes after its reset lead time");

  if (failures != 0) {
    std::cerr << failures << " transport test(s) failed\n";
    return 1;
  }
  return 0;
}
