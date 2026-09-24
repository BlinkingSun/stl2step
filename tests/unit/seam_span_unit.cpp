// D-train-seams D-S8. The span partition: a residual list in, the maximal runs
// out. A run is a maximal contiguous span of residuals within sewTol that has
// two endpoints (D-S1). No engine and no threshold of its own — sewTol is the
// caller's.
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cstdio>
#include <utility>
#include <vector>

struct Runs {
    int n = 0;
    int longest = 0;
    int excluded = 0;
    std::vector<std::pair<int, int>> iv;
};

static Runs partitionResiduals(const std::vector<double>& res, double sewTol, bool closed) {
    Runs o;
    const int n = (int)res.size();
    if (n <= 0) return o;
    const int kMinRun = 2;
    std::vector<char> on((size_t)n, 0);
    for (int i = 0; i < n; i++)
        if (res[(size_t)i] <= sewTol) on[(size_t)i] = 1;
    int i = 0;
    int covered = 0;
    while (i < n) {
        if (!on[(size_t)i]) {
            i++;
            continue;
        }
        int j = i + 1;
        while (j < n && on[(size_t)j]) j++;
        const int len = j - i;
        if (len >= kMinRun) {
            o.n++;
            covered += len;
            if (len > o.longest) o.longest = len;
            o.iv.emplace_back(i, j - 1);
        }
        i = j;
    }
    if (closed && n >= kMinRun && on[0] && on[(size_t)(n - 1)]) {
        int prefix = 0;
        while (prefix < n && on[(size_t)prefix]) prefix++;
        if (prefix < n) {
            int suffix = 0;
            while (suffix < n && on[(size_t)(n - 1 - suffix)]) suffix++;
            o.iv.erase(std::remove_if(o.iv.begin(), o.iv.end(),
                                      [&](const std::pair<int, int>& r) {
                                          const bool pref = prefix >= kMinRun && r.first == 0 &&
                                                            r.second == prefix - 1;
                                          const bool suff = suffix >= kMinRun &&
                                                            r.first == n - suffix &&
                                                            r.second == n - 1;
                                          return pref || suff;
                                      }),
                       o.iv.end());
            o.n = (int)o.iv.size();
            covered = 0;
            o.longest = 0;
            for (const auto& r : o.iv) {
                const int len = r.second - r.first + 1;
                covered += len;
                if (len > o.longest) o.longest = len;
            }
            const int len = prefix + suffix;
            if (len >= kMinRun) {
                o.iv.emplace_back(n - suffix, prefix - 1);
                o.n++;
                covered += len;
                if (len > o.longest) o.longest = len;
            }
        }
    }
    o.excluded = n - covered;
    return o;
}

static int gPass = 0;
static int gFail = 0;

static void check(bool ok, const char* name) {
    if (ok) {
        ++gPass;
        std::fprintf(stderr, "PASS %s\n", name);
    } else {
        ++gFail;
        std::fprintf(stderr, "FAIL %s\n", name);
    }
}

int main() {
    const double sew = 1.0;
    {
        const Runs r = partitionResiduals({3.0, 4.0, 5.0, 9.0}, sew, false);
        check(r.n == 0 && r.longest == 0 && r.excluded == 4, "no-run");
    }
    {
        const Runs r = partitionResiduals({0.0, 0.2, 0.1, 0.0}, sew, false);
        check(r.n == 1 && r.longest == 4 && r.excluded == 0 && r.iv.size() == 1 &&
                  r.iv[0].first == 0 && r.iv[0].second == 3,
              "whole-chain");
    }
    {
        const Runs r = partitionResiduals({0.0, 0.1, 0.2, 4.0}, sew, false);
        check(r.n == 1 && r.longest == 3 && r.excluded == 1 && r.iv.size() == 1 &&
                  r.iv[0].first == 0 && r.iv[0].second == 2,
              "excluded-terminal");
    }
    {
        const Runs r = partitionResiduals({0.0, 0.1, 5.0, 0.2, 0.0}, sew, false);
        check(r.n == 2 && r.longest == 2 && r.excluded == 1 && r.iv.size() == 2 &&
                  r.iv[0].first == 0 && r.iv[0].second == 1 && r.iv[1].first == 3 &&
                  r.iv[1].second == 4,
              "excluded-interior");
    }
    {
        // Closed: the gap is interior, the two ends are one wrapped run.
        const Runs r = partitionResiduals({0.0, 0.1, 5.0, 0.2, 0.0}, sew, true);
        check(r.n == 1 && r.longest == 4 && r.excluded == 1 && r.iv.size() == 1 &&
                  r.iv[0].first == 3 && r.iv[0].second == 1,
              "closed-wrap");
    }
    std::fprintf(stderr, "seam_span_unit %d pass %d fail\n", gPass, gFail);
    return gFail == 0 ? 0 : 1;
}
