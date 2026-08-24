#pragma once

#include "early_reflections.hpp"

#include <chrono>
#include <condition_variable>
#include <array>
#include <atomic>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <optional>
#include <thread>

namespace tfdsp {

struct EarlyReflectionBuildRequest {
  EarlyReflectionConfig config{};
  EarlyReflectionRoom room{};
  std::array<EarlyReflectionSource, EarlyReflectionMaximumSources> sources{};
  std::size_t sourceCount{};
  EarlyReflectionMaterials materials{};
  double transitionSeconds{0.100};

  void SetSources(const std::vector<EarlyReflectionSource> &newSources) {
    if (newSources.empty() ||
        newSources.size() > EarlyReflectionMaximumSources)
      throw std::invalid_argument(
          "the ER build request must contain between 1 and 8 sources");
    sourceCount = newSources.size();
    std::copy(newSources.begin(), newSources.end(), sources.begin());
  }
};

static_assert(std::is_trivially_copyable_v<EarlyReflectionBuildRequest>,
              "audio-thread ER requests must remain fixed-size POD snapshots");

struct EarlyReflectionBuildResult {
  std::size_t sequence{};
  EarlyReflectionImpulseResponse response{};
  double buildSeconds{};
  bool publishedToConvolver{};
  bool geometryReused{};
  std::exception_ptr error{};

  bool Succeeded() const noexcept { return error == nullptr; }
};

class EarlyReflectionWorker {
public:
  using Publisher = std::function<bool(const EarlyReflectionImpulseResponse &,
                                       double transitionSeconds,
                                       std::size_t sceneSequence)>;

  explicit EarlyReflectionWorker(double maximumUpdatesPerSecond = 20.0,
                                 Publisher publisher = {});
  ~EarlyReflectionWorker();

  EarlyReflectionWorker(const EarlyReflectionWorker &) = delete;
  EarlyReflectionWorker &operator=(const EarlyReflectionWorker &) = delete;

  // Single-producer, nonblocking publication intended for Rack's audio thread.
  // Returns zero only after Stop() or if the mailbox is unexpectedly saturated.
  std::size_t Submit(const EarlyReflectionBuildRequest &request) noexcept;
  std::optional<EarlyReflectionBuildResult> TryTakeLatestResult();
  std::optional<EarlyReflectionBuildResult>
  WaitForLatestResult(std::chrono::milliseconds timeout);
  void Stop() noexcept;

  double MaximumUpdatesPerSecond() const noexcept {
    return maximumUpdatesPerSecond_;
  }

private:
  static constexpr std::size_t RequestSlotCount = 4;
  enum class RequestState : std::uint8_t { Free, Writing, Ready, Reading };
  struct RequestSlot {
    EarlyReflectionBuildRequest request{};
    std::atomic<std::size_t> sequence{};
    std::atomic<RequestState> state{RequestState::Free};
  };

  void Run() noexcept;
  bool TryTakeLatestRequest(
      std::pair<std::size_t, EarlyReflectionBuildRequest> &request) noexcept;
  bool HasReadyRequest() const noexcept;
  static bool SameGeometry(const EarlyReflectionBuildRequest &left,
                           const EarlyReflectionBuildRequest &right) noexcept;

  double maximumUpdatesPerSecond_{};
  std::chrono::steady_clock::duration minimumBuildInterval_{};
  Publisher publisher_{};
  std::array<RequestSlot, RequestSlotCount> requestSlots_{};
  std::atomic<std::size_t> nextSequence_{1};
  std::atomic<bool> stopping_{};
  std::mutex wakeMutex_{};
  std::condition_variable wakeCondition_{};
  std::mutex resultMutex_{};
  std::condition_variable resultCondition_{};
  std::optional<EarlyReflectionBuildResult> completed_{};
  std::optional<EarlyReflectionBuildRequest> cachedGeometryRequest_{};
  std::vector<EarlyReflectionPath> cachedGeometry_{};
  std::thread thread_{};
};

} // namespace tfdsp
