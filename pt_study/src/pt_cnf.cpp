// Parallel tempering over CNF assignments, energy = number of violated clauses.
//
// Each replica holds a full assignment of the CNF variables.  A sweep is
// n_vars single-variable Metropolis proposals per replica; after each sweep
// adjacent replicas attempt configuration swaps (even pairs on even sweeps,
// odd pairs on odd sweeps).
//
// Moves:
//   uniform  (default) variable chosen uniformly at random, Metropolis accept.
//            Satisfies detailed balance w.r.t. exp(-E/T) at every temperature.
//   focused  variable chosen from a random violated clause, Metropolis accept.
//            Faster in practice, but NOT detailed balance; replicas then do not
//            sample Boltzmann distributions and swaps are heuristic.
//
// Modes:
//   solve    stop as soon as any replica reaches E = 0; print one result line.
//   sample   (detailed-balance test, n_vars <= 24) run a fixed number of sweeps
//            and print, per replica, a histogram over all 2^n assignments.
//
// Build: g++ -O3 -march=native -std=c++17 -o pt_cnf pt_cnf.cpp
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

struct CNF {
    int n = 0;
    std::vector<std::vector<int>> clauses;       // literals, +v / -v, 1-based
    std::vector<std::vector<int>> occ_pos, occ_neg; // clause indices per variable
};

static CNF read_dimacs(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::fprintf(stderr, "cannot open %s\n", path.c_str()); std::exit(1); }
    CNF c; std::string line; std::vector<int> cur;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == 'c') continue;
        if (line[0] == 'p') { std::istringstream is(line); std::string p, cnf; int m; is >> p >> cnf >> c.n >> m; continue; }
        std::istringstream is(line); int x;
        while (is >> x) {
            if (x == 0) { c.clauses.push_back(cur); cur.clear(); }
            else { cur.push_back(x); c.n = std::max(c.n, std::abs(x)); }
        }
    }
    c.occ_pos.assign(c.n + 1, {}); c.occ_neg.assign(c.n + 1, {});
    for (int i = 0; i < (int)c.clauses.size(); ++i)
        for (int l : c.clauses[i]) (l > 0 ? c.occ_pos[l] : c.occ_neg[-l]).push_back(i);
    return c;
}

struct Replica {
    std::vector<uint8_t> x;          // x[v] in {0,1}
    std::vector<int> ntrue;          // true literals per clause
    std::vector<int> unsat;          // list of violated clauses
    std::vector<int> pos_in_unsat;   // -1 if satisfied
    long E = 0;
};

static inline bool lit_true(const Replica& r, int l) { return l > 0 ? r.x[l] : !r.x[-l]; }

static void init_replica(const CNF& c, Replica& r, std::mt19937_64& rng) {
    r.x.assign(c.n + 1, 0);
    for (int v = 1; v <= c.n; ++v) r.x[v] = rng() & 1;
    r.ntrue.assign(c.clauses.size(), 0);
    r.pos_in_unsat.assign(c.clauses.size(), -1);
    r.unsat.clear();
    for (int i = 0; i < (int)c.clauses.size(); ++i) {
        for (int l : c.clauses[i]) r.ntrue[i] += lit_true(r, l);
        if (r.ntrue[i] == 0) { r.pos_in_unsat[i] = r.unsat.size(); r.unsat.push_back(i); }
    }
    r.E = r.unsat.size();
}

// Clauses that lose / gain a true literal when v flips.
static inline int delta_E(const CNF& c, const Replica& r, int v) {
    const auto& losing = r.x[v] ? c.occ_pos[v] : c.occ_neg[v];
    const auto& gaining = r.x[v] ? c.occ_neg[v] : c.occ_pos[v];
    int brk = 0, mk = 0;
    for (int i : losing) brk += (r.ntrue[i] == 1);
    for (int i : gaining) mk += (r.ntrue[i] == 0);
    return brk - mk;
}

static inline void add_unsat(Replica& r, int i) { r.pos_in_unsat[i] = r.unsat.size(); r.unsat.push_back(i); }
static inline void del_unsat(Replica& r, int i) {
    int p = r.pos_in_unsat[i], last = r.unsat.back();
    r.unsat[p] = last; r.pos_in_unsat[last] = p; r.unsat.pop_back(); r.pos_in_unsat[i] = -1;
}

static void flip(const CNF& c, Replica& r, int v) {
    const auto& losing = r.x[v] ? c.occ_pos[v] : c.occ_neg[v];
    const auto& gaining = r.x[v] ? c.occ_neg[v] : c.occ_pos[v];
    for (int i : losing) if (--r.ntrue[i] == 0) add_unsat(r, i);
    for (int i : gaining) if (r.ntrue[i]++ == 0) del_unsat(r, i);
    r.x[v] ^= 1;
    r.E = r.unsat.size();
}

static long recount(const CNF& c, const Replica& r) {
    long e = 0;
    for (const auto& cl : c.clauses) {
        bool s = false; for (int l : cl) if (lit_true(r, l)) { s = true; break; }
        e += !s;
    }
    return e;
}

int main(int argc, char** argv) {
    std::string path, mode = "solve", move = "uniform";
    int M = 16; double Tmin = 0.25, Tmax = 3.0; uint64_t seed = 1;
    double tlimit = 60.0; long sample_sweeps = 200000, burn = 2000;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--replicas") M = std::stoi(nx());
        else if (a == "--tmin") Tmin = std::stod(nx());
        else if (a == "--tmax") Tmax = std::stod(nx());
        else if (a == "--seed") seed = std::stoull(nx());
        else if (a == "--time") tlimit = std::stod(nx());
        else if (a == "--mode") mode = nx();
        else if (a == "--move") move = nx();
        else if (a == "--sweeps") sample_sweeps = std::stol(nx());
        else if (a == "--burn") burn = std::stol(nx());
        else path = a;
    }
    CNF c = read_dimacs(path);
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::vector<double> T(M), beta(M);
    for (int k = 0; k < M; ++k) {
        T[k] = (M == 1) ? Tmin : Tmin * std::pow(Tmax / Tmin, double(k) / (M - 1));
        beta[k] = 1.0 / T[k];
    }
    // precomputed acceptance for integer dE in [1, maxdE]
    int maxdE = 1;
    for (int v = 1; v <= c.n; ++v) maxdE = std::max<int>(maxdE, c.occ_pos[v].size() + c.occ_neg[v].size());
    std::vector<std::vector<double>> acc(M, std::vector<double>(maxdE + 1));
    for (int k = 0; k < M; ++k) for (int d = 0; d <= maxdE; ++d) acc[k][d] = std::exp(-beta[k] * d);

    std::vector<Replica> R(M);
    for (auto& r : R) init_replica(c, r, rng);
    std::vector<int> slot(M); for (int k = 0; k < M; ++k) slot[k] = k; // slot[temperature] = replica
    std::vector<long> swap_try(M, 0), swap_acc(M, 0);
    const bool focused = (move == "focused");
    const bool sample = (mode == "sample");
    if (sample && c.n > 24) { std::fprintf(stderr, "sample mode needs n <= 24\n"); return 1; }
    std::vector<std::vector<long>> hist;
    if (sample) hist.assign(M, std::vector<long>(1u << c.n, 0));

    auto t0 = std::chrono::steady_clock::now();
    long sweeps = 0, flips = 0; long bestE = 1L << 60;
    std::uniform_int_distribution<int> pickv(1, std::max(1, c.n));
    int winner = -1;
    const long nsweeps_max = sample ? sample_sweeps + burn : (1L << 62);
    for (sweeps = 0; sweeps < nsweeps_max; ++sweeps) {
        for (int k = 0; k < M; ++k) {
            Replica& r = R[slot[k]];
            for (int s = 0; s < c.n; ++s) {
                int v;
                if (focused) {
                    if (r.unsat.empty()) break;
                    const auto& cl = c.clauses[r.unsat[rng() % r.unsat.size()]];
                    v = std::abs(cl[rng() % cl.size()]);
                } else v = pickv(rng);
                int dE = delta_E(c, r, v);
                if (dE <= 0 || U(rng) < acc[k][std::min(dE, maxdE)]) { flip(c, r, v); ++flips; }
                if (!sample && r.E == 0) { winner = slot[k]; break; }
            }
            bestE = std::min(bestE, r.E);
            if (winner >= 0) break;
        }
        if (winner >= 0) break;
        // replica exchange
        for (int k = sweeps & 1; k + 1 < M; k += 2) {
            Replica& a = R[slot[k]]; Replica& b = R[slot[k + 1]];
            double lr = (beta[k] - beta[k + 1]) * double(a.E - b.E);
            ++swap_try[k];
            if (lr >= 0 || U(rng) < std::exp(lr)) { std::swap(slot[k], slot[k + 1]); ++swap_acc[k]; }
        }
        if (sample && sweeps >= burn) {
            for (int k = 0; k < M; ++k) {
                const Replica& r = R[slot[k]]; uint32_t idx = 0;
                for (int v = 1; v <= c.n; ++v) idx |= uint32_t(r.x[v]) << (v - 1);
                ++hist[k][idx];
            }
        }
        if (!sample && (sweeps & 63) == 0) {
            double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            if (el > tlimit) break;
        }
    }
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    if (sample) {
        std::printf("# n %d replicas %d sweeps %ld\n", c.n, M, sample_sweeps);
        for (int k = 0; k < M; ++k) {
            std::printf("T %.6f", T[k]);
            for (long h : hist[k]) std::printf(" %ld", h);
            std::printf("\n");
        }
        std::printf("# swap_acc");
        for (int k = 0; k + 1 < M; ++k) std::printf(" %.3f", swap_try[k] ? double(swap_acc[k]) / swap_try[k] : 0.0);
        std::printf("\n");
        return 0;
    }
    bool solved = winner >= 0;
    if (solved && recount(c, R[winner]) != 0) { std::fprintf(stderr, "internal error: E mismatch\n"); return 2; }
    std::printf("result %s time %.4f sweeps %ld flips %ld bestE %ld swap_acc", solved ? "SAT" : "TIMEOUT", el, sweeps, flips, bestE);
    for (int k = 0; k + 1 < M; ++k) std::printf(" %.2f", swap_try[k] ? double(swap_acc[k]) / swap_try[k] : 0.0);
    std::printf("\n");
    if (solved) {
        std::printf("v");
        for (int v = 1; v <= c.n; ++v) std::printf(" %d", R[winner].x[v] ? v : -v);
        std::printf(" 0\n");
    }
    return 0;
}
