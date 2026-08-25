#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>

namespace rack::engine {
struct Module;
}

namespace tftransport {

// A fixed transport rate avoids patch-local clock conventions. Twenty-four
// pulses per quarter note is the established MIDI/DIN transport resolution;
// downstream modules retain quarter-note beat semantics while receiving enough
// timing information to lock quickly at slow tempos.
inline constexpr int PulsesPerQuarterNote = 24;
inline constexpr double BeatsPerPulse =
    1.0 / static_cast<double>(PulsesPerQuarterNote);

inline constexpr bool IsQuarterNotePulse(std::int64_t pulse) noexcept {
  return pulse % PulsesPerQuarterNote == 0;
}

inline constexpr double BeatAtPulse(std::int64_t pulse) noexcept {
  return static_cast<double>(pulse) * BeatsPerPulse;
}

inline constexpr double
BeatPeriodFromPulseSamples(std::int64_t samples) noexcept {
  return static_cast<double>(samples) * PulsesPerQuarterNote;
}

enum class State { Stopped = 0, Paused = 1, Playing = 2 };
enum class Command { PlayFromBeginning, Pause, Play, Stop, TogglePlayPause };

// UI actions publish one command without touching audio-thread transport
// state. The Transport module consumes the latest request in process().
class CommandMailbox {
public:
  void post(Command command) noexcept {
    pending_.store(static_cast<unsigned>(command) + 1,
                   std::memory_order_release);
  }

  bool consume(Command &command) noexcept {
    const unsigned pending = pending_.exchange(0, std::memory_order_acq_rel);
    if (pending == 0)
      return false;
    command = static_cast<Command>(pending - 1);
    return true;
  }

private:
  std::atomic<unsigned> pending_{0};
};

// Implemented by the Rack-facing Transport module. The source must be that
// module's RUN output, which makes the keyboard association explicit.
bool RequestModuleCommand(rack::engine::Module *target, int sourceOutputId,
                          Command command) noexcept;

struct Frame {
  bool clock = false;
  bool reset = false;
  bool run = false;
};

class Engine {
public:
  static constexpr double ResetLeadSeconds = 0.0015;

  State state() const noexcept { return state_; }
  double phase() const noexcept { return phase_; }

  void command(Command command) noexcept {
    switch (command) {
    case Command::PlayFromBeginning:
      state_ = State::Stopped;
      phase_ = 0.0;
      restartDelay_ = ResetLeadSeconds;
      restartPending_ = true;
      clockOnStart_ = false;
      resetPending_ = true;
      break;
    case Command::Pause:
      if (state_ == State::Playing)
        state_ = State::Paused;
      restartPending_ = false;
      clockOnStart_ = false;
      break;
    case Command::Play:
      restartPending_ = false;
      if (state_ == State::Stopped) {
        phase_ = 0.0;
        clockOnStart_ = true;
      }
      if (state_ != State::Playing)
        state_ = State::Playing;
      break;
    case Command::Stop:
      state_ = State::Stopped;
      phase_ = 0.0;
      restartDelay_ = 0.0;
      restartPending_ = false;
      clockOnStart_ = false;
      resetPending_ = true;
      break;
    case Command::TogglePlayPause:
      this->command(state_ == State::Playing ? Command::Pause : Command::Play);
      break;
    }
  }

  void loadState(State state) noexcept {
    state_ = state;
    phase_ = 0.0;
    restartDelay_ = 0.0;
    restartPending_ = false;
    resetPending_ = false;
    clockOnStart_ = state == State::Playing;
  }

  // CLOCK has one fixed transport resolution. Musical durations remain the
  // responsibility of downstream sequencers.
  Frame process(double sampleTime, double bpm) noexcept {
    Frame frame;
    frame.reset = resetPending_;
    resetPending_ = false;

    if (restartPending_) {
      restartDelay_ -= std::max(0.0, sampleTime);
      if (restartDelay_ <= 0.0) {
        restartPending_ = false;
        state_ = State::Playing;
        clockOnStart_ = true;
      }
    }

    frame.run = state_ == State::Playing;
    if (!frame.run)
      return frame;

    if (clockOnStart_) {
      frame.clock = true;
      clockOnStart_ = false;
    }

    const double pulsesPerSecond =
        std::clamp(bpm, 1.0, 1000.0) * PulsesPerQuarterNote / 60.0;
    phase_ += std::max(0.0, sampleTime) * pulsesPerSecond;
    if (phase_ >= 1.0) {
      phase_ -= std::floor(phase_);
      frame.clock = true;
    }
    return frame;
  }

private:
  State state_ = State::Stopped;
  double phase_ = 0.0;
  double restartDelay_ = 0.0;
  bool restartPending_ = false;
  bool resetPending_ = false;
  bool clockOnStart_ = false;
};

} // namespace tftransport
