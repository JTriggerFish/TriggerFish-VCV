#pragma once

#include <atomic>
#include <cstdint>

namespace tfseq {

// Single-publisher/single-consumer ownership handoff for immutable programs.
//
// The UI may replace an unpublished candidate while the audio thread inspects
// it. Publishing the consumer's hazard before rechecking the tagged pointer
// prevents deletion under that inspection. claim() uses compare/exchange so a
// policy decision made for one candidate is never applied to its replacement.
// Reclamation remains on the UI thread.
template <typename T, std::uintptr_t PolicyMask> class ProgramMailbox {
public:
  static_assert((alignof(T) & PolicyMask) == 0,
                "mailbox policy bits must fit the pointee alignment");

  ProgramMailbox() = default;
  ProgramMailbox(const ProgramMailbox &) = delete;
  ProgramMailbox &operator=(const ProgramMailbox &) = delete;

  ~ProgramMailbox() {
    auto *queued = pointer(pending_.exchange(0, std::memory_order_acq_rel));
    if (queued == deferred_)
      deferred_ = nullptr;
    delete queued;
    delete deferred_;
  }

  static T *pointer(const std::uintptr_t tagged) noexcept {
    return reinterpret_cast<T *>(tagged & ~PolicyMask);
  }

  // UI thread only. Ownership of candidate transfers to the mailbox.
  void publish(T *candidate, const std::uintptr_t policy = 0) noexcept {
    collect();
    const auto tagged =
        reinterpret_cast<std::uintptr_t>(candidate) | (policy & PolicyMask);
    retire(pending_.exchange(tagged, std::memory_order_acq_rel));
    collect();
  }

  // UI thread only.
  void clear() noexcept {
    collect();
    retire(pending_.exchange(0, std::memory_order_acq_rel));
    collect();
  }

  // UI thread only. Safe to call once per widget frame.
  void collect() noexcept {
    if (deferred_ && hazard_.load(std::memory_order_seq_cst) != deferred_) {
      delete deferred_;
      deferred_ = nullptr;
    }
  }

  // Audio thread only. The tagged value remains safe to inspect through
  // pointer() until claim() or release().
  std::uintptr_t protect() noexcept {
    for (;;) {
      const auto tagged = pending_.load(std::memory_order_acquire);
      auto *candidate = pointer(tagged);
      hazard_.store(candidate, std::memory_order_seq_cst);
      if (!candidate || pending_.load(std::memory_order_seq_cst) == tagged)
        return tagged;
      hazard_.store(nullptr, std::memory_order_seq_cst);
    }
  }

  // Audio thread only. A concurrent publication makes this fail instead of
  // claiming a candidate that was not the one inspected.
  bool claim(const std::uintptr_t tagged) noexcept {
    auto expected = tagged;
    const bool claimed = pending_.compare_exchange_strong(
        expected, 0, std::memory_order_acq_rel, std::memory_order_acquire);
    hazard_.store(nullptr, std::memory_order_seq_cst);
    return claimed;
  }

  // Audio thread only.
  void release() noexcept { hazard_.store(nullptr, std::memory_order_seq_cst); }

  std::uintptr_t pending() const noexcept {
    return pending_.load(std::memory_order_acquire);
  }

private:
  std::atomic<std::uintptr_t> pending_{0};
  std::atomic<T *> hazard_{nullptr};
  // A single consumer protects at most one object. UI thread only.
  T *deferred_ = nullptr;

  void retire(const std::uintptr_t tagged) noexcept {
    auto *candidate = pointer(tagged);
    if (!candidate)
      return;
    auto *hazard = hazard_.load(std::memory_order_seq_cst);
    if (deferred_ && deferred_ != hazard) {
      delete deferred_;
      deferred_ = nullptr;
    }
    if (hazard == candidate) {
      deferred_ = candidate;
      return;
    }
    delete candidate;
  }
};

} // namespace tfseq
