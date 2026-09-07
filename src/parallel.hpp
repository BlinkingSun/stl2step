// stl2step internal thread helper (D-142-1). Owns every std::thread the engine
// constructs. requestedThreads <= 1 => inline, zero std::thread.
// Not installed; not part of the public API.
//
// SPDX-License-Identifier: MIT

#ifndef STL2STEP_PARALLEL_HPP
#define STL2STEP_PARALLEL_HPP

#include <cstddef>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

namespace stl2step {
namespace detail {

// Today's fallback when hardware_concurrency() returns 0 (stl2step.cpp:316,
// mesh.cpp:199). Kept in ONE place.
constexpr unsigned kHwFallback = 4;

unsigned hardwareThreadCount();
unsigned resolveThreadCount(int requested);

std::uint64_t threadsSpawnedForTest();
void          resetThreadsSpawnedForTest();
void          noteThreadSpawn();

int  currentRequestedThreads();
void setCurrentRequestedThreads(int requested);

struct ThreadBudgetScope {
    int prev;
    explicit ThreadBudgetScope(int requested)
        : prev(currentRequestedThreads()) {
        setCurrentRequestedThreads(requested);
    }
    ~ThreadBudgetScope() { setCurrentRequestedThreads(prev); }
    ThreadBudgetScope(const ThreadBudgetScope&)            = delete;
    ThreadBudgetScope& operator=(const ThreadBudgetScope&) = delete;
};

// Run workerBody(workerIndex) on nWorkers threads. nWorkers <= 1 => inline
// on the calling thread (zero std::thread). Each worker inherits the caller's
// requested-thread budget so nested pools stay consistent under reentrancy.
template <typename Fn>
void runPool(unsigned nWorkers, Fn&& workerBody) {
    if (nWorkers <= 1) {
        workerBody(0u);
        return;
    }
    const int req = currentRequestedThreads();
    std::vector<std::thread> pool;
    pool.reserve(nWorkers);
    for (unsigned w = 0; w < nWorkers; ++w) {
        pool.emplace_back([workerBody, w, req]() {
            ThreadBudgetScope inner(req);
            workerBody(w);
        });
        noteThreadSpawn();
    }
    for (auto& t : pool) t.join();
}

// Overlap bg on one worker with fg on the caller. At resolved threads <= 1
// both run inline (bg then fg) — zero std::thread.
template <typename Bg, typename Fg>
void overlap(int requested, Bg&& bg, Fg&& fg) {
    if (resolveThreadCount(requested) <= 1) {
        bg();
        fg();
        return;
    }
    const int req = currentRequestedThreads();
    noteThreadSpawn();
    std::thread t([bg = std::forward<Bg>(bg), req]() {
        ThreadBudgetScope inner(req);
        bg();
    });
    fg();
    t.join();
}

}  // namespace detail
}  // namespace stl2step

#endif  // STL2STEP_PARALLEL_HPP
