// stl2step internal thread helper (D-142-1).
// SPDX-License-Identifier: MIT

#include "parallel.hpp"

#include <algorithm>
#include <atomic>
#include <thread>

namespace stl2step {
namespace detail {
namespace {

std::atomic<std::uint64_t> gSpawned{0};
thread_local int           tRequested = 0;

}  // namespace

unsigned hardwareThreadCount() {
    unsigned hw = std::thread::hardware_concurrency();
    if (hw == 0) hw = kHwFallback;
    return hw;
}

unsigned resolveThreadCount(int requested) {
    if (requested == 1) return 1;
    if (requested > 1) {
        unsigned hw = hardwareThreadCount();
        return std::min(hw, static_cast<unsigned>(requested));
    }
    // 0 (and negatives) = all cores.
    return hardwareThreadCount();
}

std::uint64_t threadsSpawnedForTest() { return gSpawned.load(); }

void resetThreadsSpawnedForTest() { gSpawned.store(0); }

void noteThreadSpawn() { gSpawned.fetch_add(1); }

int currentRequestedThreads() { return tRequested; }

void setCurrentRequestedThreads(int requested) { tRequested = requested; }

}  // namespace detail
}  // namespace stl2step
