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

// Predict whether the clock edge currently being processed is a quarter-note
// boundary in the receiver's own musical timeline. Before the first observed
// edge, that edge establishes beat zero; afterwards currentPulse is the most
// recently processed pulse.
inline constexpr bool
IsQuarterBoundaryOnClockEdge(bool clockSeen,
                             std::int64_t currentPulse) noexcept {
  return !clockSeen || IsQuarterNotePulse(currentPulse + 1);
}

inline constexpr double
BeatPeriodFromPulseSamples(std::int64_t samples) noexcept {
  return static_cast<double>(samples) * PulsesPerQuarterNote;
}

// Preserve elapsed wall-clock time when Rack changes the engine sample rate.
// Period estimates and partially measured pulse intervals must be converted
// together; converting only the period makes the first interval after an
// audio-device change appear spuriously early or late.
inline std::int64_t RescaleSampleCount(std::int64_t samples,
                                      double previousSampleRate,
                                      double nextSampleRate) noexcept {
  if (samples <= 0 || !std::isfinite(previousSampleRate) ||
      !std::isfinite(nextSampleRate) || previousSampleRate <= 0.0 ||
      nextSampleRate <= 0.0)
    return std::max<std::int64_t>(0, samples);
  return std::max<std::int64_t>(
      1, std::llround(static_cast<double>(samples) * nextSampleRate /
                      previousSampleRate));
}

// Incorporate a complete adjacent-pulse interval only when it is continuous.
// In particular, skipping the first measurement after a local/RUN pause must
// not claim that a previously unknown period is now valid: a zero-valued
// "known" period would immediately trip downstream clock-timeout guards.
inline bool UpdateBeatPeriodEstimate(std::int64_t pulseIntervalSamples,
                                     bool continuous,
                                     bool &periodKnown,
                                     double &periodSamples) noexcept {
  if (!continuous || pulseIntervalSamples <= 0)
    return false;
  const double measured =
      BeatPeriodFromPulseSamples(pulseIntervalSamples);
  periodSamples =
      periodKnown ? 0.75 * periodSamples + 0.25 * measured : measured;
  periodKnown = true;
  return true;
}

// Local mute is immediate, but local resume waits for a known master quarter
// boundary. The sequencer's score keeps following the shared clock while
// muted, so its arrangement and lane cursors remain synchronized and visible.
class QuantizedLocalPlayback {
public:
  enum class ToggleResult { Muted, ResumeQueued, ResumeCanceled };

  void reset() noexcept {
    audible_ = true;
    resumeQueued_ = false;
  }

  ToggleResult toggle() noexcept {
    if (audible_) {
      audible_ = false;
      resumeQueued_ = false;
      return ToggleResult::Muted;
    }
    resumeQueued_ = !resumeQueued_;
    return resumeQueued_ ? ToggleResult::ResumeQueued
                         : ToggleResult::ResumeCanceled;
  }

  bool applyQuarterBoundary() noexcept {
    if (!resumeQueued_)
      return false;
    audible_ = true;
    resumeQueued_ = false;
    return true;
  }

  bool applyClockEdge(bool clockSeen, std::int64_t currentPulse) noexcept {
    return IsQuarterBoundaryOnClockEdge(clockSeen, currentPulse) &&
           applyQuarterBoundary();
  }

  bool audible() const noexcept { return audible_; }
  bool resumeQueued() const noexcept { return resumeQueued_; }

private:
  bool audible_ = true;
  bool resumeQueued_ = false;
};

struct EventOutputVoltages {
  float gate = 0.f;
  float trigger = 0.f;
  float accent = 0.f;
};

// Trigger processing is deliberately unconditional: muting hides a trigger
// pulse but must not freeze it and leak the stale pulse on a later unmute.
template <typename TriggerProcessor>
inline EventOutputVoltages
ProcessEventOutputs(bool enabled, bool gateHigh, float accent,
                    TriggerProcessor processTrigger) noexcept {
  const bool triggerHigh = processTrigger();
  return {enabled && gateHigh ? 10.f : 0.f,
          enabled && triggerHigh ? 10.f : 0.f,
          enabled && gateHigh ? accent : 0.f};
}

// A transport module's RUN and CLOCK cables can become visible to different
// downstream workers one engine frame apart. Rack's Schmitt trigger consumes
// the CLOCK edge immediately, so retain that edge only while the same pulse is
// still high and deliver it when RUN catches up. A genuinely stopped or stale
// pulse is discarded as soon as CLOCK returns low.
class RunSynchronizedClockEdge {
public:
  void reset() noexcept { pending_ = false; }

  bool process(bool detectedEdge, bool clockHigh,
               bool transportRunning) noexcept {
    if (!clockHigh)
      pending_ = false;
    if (detectedEdge && !transportRunning) {
      pending_ = clockHigh;
      return false;
    }
    if (!transportRunning)
      return false;
    if (detectedEdge) {
      pending_ = false;
      return true;
    }
    if (pending_ && clockHigh) {
      pending_ = false;
      return true;
    }
    return false;
  }

  bool pending() const noexcept { return pending_; }

private:
  bool pending_ = false;
};

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
