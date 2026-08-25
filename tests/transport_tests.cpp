#include "tftransport.hpp"

#include <cmath>
#include <iostream>

namespace {
int failures = 0;

void check(bool condition, const char *message) {
  if (condition)
    return;
  std::cerr << "FAIL: " << message << '\n';
  ++failures;
}
} // namespace

int main() {
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
