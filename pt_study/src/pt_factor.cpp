// Parallel tempering over the factor bits only.
//
// State: p (np bits) and q (nq bits) with the most and least significant bits
// fixed to 1, so 2d-4 free bits for d-bit primes.  Energy: Hamming distance
// between the bits of p*q and the bits of N, E = popcount(p*q XOR N).
// E = 0 iff p*q = N.  Single-bit Metropolis flips (detailed balance holds),
// adjacent-temperature swaps after every sweep.
//
// Modes: solve (default) and sample (histogram over all states, free bits <= 24).
// Build: g++ -O3 -march=native -std=c++17 -o pt_factor pt_factor.cpp
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

typedef unsigned __int128 u128;
static inline int popc(u128 x) { return __builtin_popcountll((uint64_t)x) + __builtin_popcountll((uint64_t)(x >> 64)); }

struct State { uint64_t p, q; int E; };

int main(int argc, char** argv) {
    uint64_t N = 0; int np = 0, nq = 0;
    std::string mode = "solve";
    int M = 12; double Tmin = 0.3, Tmax = 4.0; uint64_t seed = 1; double tlimit = 60;
    long sample_sweeps = 200000, burn = 2000;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--N") N = std::stoull(nx());
        else if (a == "--np") np = std::stoi(nx());
        else if (a == "--nq") nq = std::stoi(nx());
        else if (a == "--replicas") M = std::stoi(nx());
        else if (a == "--tmin") Tmin = std::stod(nx());
        else if (a == "--tmax") Tmax = std::stod(nx());
        else if (a == "--seed") seed = std::stoull(nx());
        else if (a == "--time") tlimit = std::stod(nx());
        else if (a == "--mode") mode = nx();
        else if (a == "--sweeps") sample_sweeps = std::stol(nx());
        else if (a == "--burn") burn = std::stol(nx());
    }
    if (!N || np < 2 || nq < 2 || np + nq > 62) { std::fprintf(stderr, "need --N --np --nq (np,nq >= 2)\n"); return 1; }
    // free bit positions: p bits 1..np-2, q bits 1..nq-2
    std::vector<std::pair<int, int>> fb; // (which, bit)
    for (int b = 1; b <= np - 2; ++b) fb.push_back({0, b});
    for (int b = 1; b <= nq - 2; ++b) fb.push_back({1, b});
    const int nf = fb.size();
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::vector<double> beta(M), T(M);
    for (int k = 0; k < M; ++k) {
        T[k] = (M == 1) ? Tmin : Tmin * std::pow(Tmax / Tmin, double(k) / (M - 1));
        beta[k] = 1.0 / T[k];
    }
    const int maxdE = np + nq + 2;
    std::vector<std::vector<double>> acc(M, std::vector<double>(maxdE + 1));
    for (int k = 0; k < M; ++k) for (int d = 0; d <= maxdE; ++d) acc[k][d] = std::exp(-beta[k] * d);
    auto energy = [&](uint64_t p, uint64_t q) { return popc((u128)p * q ^ (u128)N); };
    const uint64_t pfix = 1ull | (1ull << (np - 1)), qfix = 1ull | (1ull << (nq - 1));
    std::vector<State> R(M);
    for (auto& s : R) {
        s.p = pfix; s.q = qfix;
        for (auto [w, b] : fb) if (rng() & 1) (w ? s.q : s.p) |= 1ull << b;
        s.E = energy(s.p, s.q);
    }
    std::vector<int> slot(M); for (int k = 0; k < M; ++k) slot[k] = k;
    std::vector<long> st(M, 0), sa(M, 0);
    const bool sample = mode == "sample";
    if (sample && nf > 24) { std::fprintf(stderr, "sample mode needs <= 24 free bits\n"); return 1; }
    std::vector<std::vector<long>> hist;
    if (sample) hist.assign(M, std::vector<long>(1u << nf, 0));
    auto t0 = std::chrono::steady_clock::now();
    long sweeps, evals = 0; int winner = -1;
    const long smax = sample ? sample_sweeps + burn : (1L << 62);
    std::uniform_int_distribution<int> pick(0, std::max(0, nf - 1));
    for (sweeps = 0; sweeps < smax; ++sweeps) {
        for (int k = 0; k < M && winner < 0; ++k) {
            State& s = R[slot[k]];
            for (int t = 0; t < std::max(nf, 1); ++t) {
                if (nf == 0) break;
                auto [w, b] = fb[pick(rng)];
                uint64_t p2 = s.p ^ (w ? 0 : 1ull << b), q2 = s.q ^ (w ? 1ull << b : 0);
                int E2 = energy(p2, q2); ++evals;
                int dE = E2 - s.E;
                if (dE <= 0 || U(rng) < acc[k][std::min(dE, maxdE)]) { s.p = p2; s.q = q2; s.E = E2; }
                if (!sample && s.E == 0) { winner = slot[k]; break; }
            }
        }
        if (winner >= 0) break;
        for (int k = sweeps & 1; k + 1 < M; k += 2) {
            State& a = R[slot[k]]; State& b = R[slot[k + 1]];
            double lr = (beta[k] - beta[k + 1]) * double(a.E - b.E); ++st[k];
            if (lr >= 0 || U(rng) < std::exp(lr)) { std::swap(slot[k], slot[k + 1]); ++sa[k]; }
        }
        if (sample && sweeps >= burn) {
            for (int k = 0; k < M; ++k) {
                const State& s = R[slot[k]]; uint32_t idx = 0;
                for (int i = 0; i < nf; ++i) {
                    auto [w, b] = fb[i];
                    idx |= uint32_t(((w ? s.q : s.p) >> b) & 1) << i;
                }
                ++hist[k][idx];
            }
        }
        if (!sample && (sweeps & 255) == 0 &&
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() > tlimit) break;
    }
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (sample) {
        std::printf("# nfree %d replicas %d sweeps %ld\n", nf, M, sample_sweeps);
        for (int k = 0; k < M; ++k) {
            std::printf("T %.6f", T[k]);
            for (long h : hist[k]) std::printf(" %ld", h);
            std::printf("\n");
        }
        return 0;
    }
    bool ok = winner >= 0;
    if (ok && (u128)R[winner].p * R[winner].q != (u128)N) { std::fprintf(stderr, "internal error\n"); return 2; }
    std::printf("result %s time %.4f sweeps %ld evals %ld", ok ? "SAT" : "TIMEOUT", el, sweeps, evals);
    if (ok) std::printf(" p %llu q %llu", (unsigned long long)R[winner].p, (unsigned long long)R[winner].q);
    std::printf(" swap_acc");
    for (int k = 0; k + 1 < M; ++k) std::printf(" %.2f", st[k] ? double(sa[k]) / st[k] : 0.0);
    std::printf("\n");
    return 0;
}
