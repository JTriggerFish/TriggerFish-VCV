#include "early_reflections_worker.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace tfdsp {

EarlyReflectionWorker::EarlyReflectionWorker(
    const double maximumUpdatesPerSecond, Publisher publisher)
    : maximumUpdatesPerSecond_(maximumUpdatesPerSecond),
      publisher_(std::move(publisher)) {
  if (!std::isfinite(maximumUpdatesPerSecond_) ||
      maximumUpdatesPerSecond_ <= 0.0 || maximumUpdatesPerSecond_ > 20.0)
    throw std::invalid_argument(
        "ER worker update rate must lie in (0, 20] updates per second");
  minimumBuildInterval_ =
      std::chrono::duration_cast<std::chrono::steady_clock::duration>(
          std::chrono::duration<double>(1.0 / maximumUpdatesPerSecond_));
  thread_ = std::thread(&EarlyReflectionWorker::Run, this);
}

EarlyReflectionWorker::~EarlyReflectionWorker() { Stop(); }

std::size_t EarlyReflectionWorker::Submit(
    const EarlyReflectionBuildRequest &request) noexcept {
  if (stopping_.load(std::memory_order_acquire))
    return 0;

  std::size_t slotIndex = RequestSlotCount;
  for (std::size_t candidate = 0; candidate < RequestSlotCount; ++candidate) {
    RequestState expected = RequestState::Free;
    if (requestSlots_[candidate].state.compare_exchange_strong(
            expected, RequestState::Writing, std::memory_order_acq_rel,
            std::memory_order_acquire)) {
      slotIndex = candidate;
      break;
    }
  }
  if (slotIndex == RequestSlotCount)
    return 0;

  const std::size_t sequence =
      nextSequence_.fetch_add(1, std::memory_order_relaxed);
  auto &slot = requestSlots_[slotIndex];
  slot.request = request;
  slot.sequence = sequence;
  slot.state.store(RequestState::Ready, std::memory_order_release);
  const int replaced = latestRequestSlot_.exchange(
      static_cast<int>(slotIndex), std::memory_order_acq_rel);
  if (replaced >= 0) {
    RequestState expected = RequestState::Ready;
    requestSlots_[static_cast<std::size_t>(replaced)]
        .state.compare_exchange_strong(expected, RequestState::Free,
                                       std::memory_order_acq_rel,
                                       std::memory_order_acquire);
  }
  wakeCondition_.notify_one();
  return sequence;
}

std::optional<EarlyReflectionBuildResult>
EarlyReflectionWorker::TryTakeLatestResult() {
  std::lock_guard<std::mutex> lock(resultMutex_);
  if (!completed_)
    return std::nullopt;
  auto result = std::move(completed_);
  completed_.reset();
  return result;
}

std::optional<EarlyReflectionBuildResult>
EarlyReflectionWorker::WaitForLatestResult(
  const std::chrono::milliseconds timeout) {
  std::unique_lock<std::mutex> lock(resultMutex_);
  resultCondition_.wait_for(
      lock, timeout, [&] {
        return completed_.has_value() ||
               stopping_.load(std::memory_order_acquire);
      });
  if (!completed_)
    return std::nullopt;
  auto result = std::move(completed_);
  completed_.reset();
  return result;
}

void EarlyReflectionWorker::Stop() noexcept {
  if (stopping_.exchange(true, std::memory_order_acq_rel))
    return;
  latestRequestSlot_.store(InvalidRequestSlot, std::memory_order_release);
  wakeCondition_.notify_all();
  resultCondition_.notify_all();
  if (thread_.joinable())
    thread_.join();
}

bool EarlyReflectionWorker::SameGeometry(
    const EarlyReflectionBuildRequest &left,
    const EarlyReflectionBuildRequest &right) noexcept {
  if (left.diffusion != right.diffusion ||
      left.config.speedOfSound != right.config.speedOfSound ||
      left.config.responseTimeSeconds != right.config.responseTimeSeconds ||
      left.config.minimumHandoffSeconds != right.config.minimumHandoffSeconds ||
      left.config.handoffOverlapSeconds != right.config.handoffOverlapSeconds ||
      left.config.analysisSafetySeconds != right.config.analysisSafetySeconds ||
      left.room.dimensionsMetres != right.room.dimensionsMetres ||
      left.room.listenerPositionMetres != right.room.listenerPositionMetres ||
      left.sourceCount != right.sourceCount)
    return false;
  for (std::size_t source = 0; source < left.sourceCount; ++source)
    if (left.sources[source].positionMetres !=
        right.sources[source].positionMetres)
      return false;
  return true;
}

bool EarlyReflectionWorker::TryTakeLatestRequest(
    std::pair<std::size_t, EarlyReflectionBuildRequest> &request) noexcept {
  const int candidate =
      latestRequestSlot_.exchange(InvalidRequestSlot, std::memory_order_acq_rel);
  if (candidate < 0)
    return false;
  auto &slot = requestSlots_[static_cast<std::size_t>(candidate)];
  RequestState expected = RequestState::Ready;
  if (!slot.state.compare_exchange_strong(expected, RequestState::Reading,
                                          std::memory_order_acq_rel,
                                          std::memory_order_acquire))
    return false;
  request = {slot.sequence, slot.request};
  slot.state.store(RequestState::Free, std::memory_order_release);
  return true;
}

void EarlyReflectionWorker::Run() noexcept {
  using Clock = std::chrono::steady_clock;
  auto nextBuild = Clock::time_point::min();
  for (;;) {
    std::pair<std::size_t, EarlyReflectionBuildRequest> work;
    while (!TryTakeLatestRequest(work)) {
      std::unique_lock<std::mutex> lock(wakeMutex_);
      wakeCondition_.wait(lock, [&] {
        return stopping_.load(std::memory_order_acquire) ||
               latestRequestSlot_.load(std::memory_order_acquire) >= 0;
      });
      if (stopping_.load(std::memory_order_acquire))
        return;
    }
    if (Clock::now() < nextBuild) {
      std::unique_lock<std::mutex> lock(wakeMutex_);
      wakeCondition_.wait_until(lock, nextBuild, [&] {
        return stopping_.load(std::memory_order_acquire);
      });
      if (stopping_.load(std::memory_order_acquire))
        return;
      // A request submitted during the rate-limit interval supersedes the
      // snapshot already taken.
      std::pair<std::size_t, EarlyReflectionBuildRequest> newer;
      if (TryTakeLatestRequest(newer))
        work = newer;
    }

    const auto started = Clock::now();
    nextBuild = started + minimumBuildInterval_;
    EarlyReflectionBuildResult result;
    result.sequence = work.first;
    try {
      if (work.second.sourceCount == 0 ||
          work.second.sourceCount > EarlyReflectionMaximumSources)
        throw std::invalid_argument(
            "the ER build request must contain between 1 and 8 sources");
      const std::vector<EarlyReflectionSource> sources(
          work.second.sources.begin(),
          work.second.sources.begin() +
              static_cast<std::ptrdiff_t>(work.second.sourceCount));
      result.geometryReused =
          cachedGeometryRequest_ &&
          SameGeometry(*cachedGeometryRequest_, work.second);
      if (!result.geometryReused) {
        cachedGeometry_ = EnumerateEarlyReflectionPaths(
            work.second.config, work.second.room, sources, work.second.materials,
            work.second.diffusion);
        cachedGeometryRequest_ = work.second;
      }
      result.response = BuildEarlyReflectionImpulseResponse(
          work.second.config, work.second.room, sources,
          work.second.materials, work.second.convolutionLatencySamples,
          work.second.diffusion, &cachedGeometry_);
      result.buildSeconds =
          std::chrono::duration<double>(Clock::now() - started).count();

      if (latestRequestSlot_.load(std::memory_order_acquire) >= 0)
        continue;

      if (publisher_) {
        if (!publisher_(result.response, work.second.transitionSeconds))
          throw std::runtime_error(
              "no free ER convolution bank was available to the worker");
        result.publishedToConvolver = true;
      }
    } catch (...) {
      result.buildSeconds =
          std::chrono::duration<double>(Clock::now() - started).count();
      result.error = std::current_exception();
    }

    {
      std::lock_guard<std::mutex> lock(resultMutex_);
      completed_ = std::move(result);
    }
    resultCondition_.notify_all();
  }
}

} // namespace tfdsp
