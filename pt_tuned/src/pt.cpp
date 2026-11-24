// Parallel tempering for the pq factoring instances, two models:
//
//   --model cnf <file.cnf>     state = assignment of the residual CNF,
//                              E = number of violated clauses
//   --model factor --N --np --nq
//                              state = free bits of p and q (MSB, LSB fixed 1),
//                              E = popcount(p*q XOR N)
//
// Local moves: Metropolis single flips.  cnf: --move uniform (detailed balance,
// default) or focused (variable from a violated clause; no detailed balance).
// A sweep = (number of free variables) proposals per replica.  After every
// sweep, adjacent temperatures attempt swaps (even/odd pairs alternating).
//
// Temperatures: --temps t1,t2,...  (ascending), or --tmin --tmax --replicas
// (geometric).
//
// Modes
//   solve   stop when any replica reaches E = 0 or after --time seconds.
//   tune    run --burn + --sweeps sweeps without stopping and report, per
//           temperature: swap acceptance with the next one, the flow fraction
//           f = n_up / (n_up + n_down) (up = last visited the coldest
//           temperature, down = last visited the hottest), <E>, var(E), the
//           fraction of samples with E = 0, 1, 2, and the number of round
//           trips (cold -> hot -> cold) completed by all replicas.
//   single  one temperature (--tmin), --sweeps sweeps after --burn: <E>,
//           var(E), flip acceptance, integrated autocorrelation time of E in
//           sweeps (Sokal automatic window, c = 6).
//
// Build: g++ -O3 -march=native -std=c++17 -o pt pt.cpp
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

typedef std::mt19937_64 RNG;
typedef unsigned __int128 u128;

// ------------------------------------------------------------------ CNF model
struct CnfModel {
    int n = 0;
    std::vector<std::vector<int>> cl, occ_pos, occ_neg;
    bool focused = false;
    int maxdE = 1;
    struct State {
        std::vector<uint8_t> x; std::vector<int> nt, unsat, pos; long E = 0;
    };
    void load(const std::string& path) {
        std::ifstream f(path);
        if (!f) { std::fprintf(stderr, "cannot open %s\n", path.c_str()); std::exit(1); }
        std::string line; std::vector<int> cur;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == 'c' || line[0] == 'p') continue;
            std::istringstream is(line); int v;
            while (is >> v) { if (!v) { cl.push_back(cur); cur.clear(); } else { cur.push_back(v); n = std::max(n, std::abs(v)); } }
        }
        occ_pos.assign(n + 1, {}); occ_neg.assign(n + 1, {});
        for (int i = 0; i < (int)cl.size(); ++i) for (int l : cl[i]) (l > 0 ? occ_pos[l] : occ_neg[-l]).push_back(i);
        for (int v = 1; v <= n; ++v) maxdE = std::max<int>(maxdE, occ_pos[v].size() + occ_neg[v].size());
    }
    int sweep_size() const { return n; }
    bool lt(const State& s, int l) const { return l > 0 ? s.x[l] : !s.x[-l]; }
    void init(State& s, RNG& rng) const {
        s.x.assign(n + 1, 0); for (int v = 1; v <= n; ++v) s.x[v] = rng() & 1;
        s.nt.assign(cl.size(), 0); s.pos.assign(cl.size(), -1); s.unsat.clear();
        for (int i = 0; i < (int)cl.size(); ++i) {
            for (int l : cl[i]) s.nt[i] += lt(s, l);
            if (!s.nt[i]) { s.pos[i] = s.unsat.size(); s.unsat.push_back(i); }
        }
        s.E = s.unsat.size();
    }
    int dE(const State& s, int v) const {
        const auto& lo = s.x[v] ? occ_pos[v] : occ_neg[v];
        const auto& ga = s.x[v] ? occ_neg[v] : occ_pos[v];
        int b = 0, m = 0;
        for (int i : lo) b += (s.nt[i] == 1);
        for (int i : ga) m += (s.nt[i] == 0);
        return b - m;
    }
    void flip(State& s, int v) const {
        const auto& lo = s.x[v] ? occ_pos[v] : occ_neg[v];
        const auto& ga = s.x[v] ? occ_neg[v] : occ_pos[v];
        for (int i : lo) if (--s.nt[i] == 0) { s.pos[i] = s.unsat.size(); s.unsat.push_back(i); }
        for (int i : ga) if (s.nt[i]++ == 0) {
            int p = s.pos[i], last = s.unsat.back(); s.unsat[p] = last; s.pos[last] = p; s.unsat.pop_back(); s.pos[i] = -1;
        }
        s.x[v] ^= 1; s.E = s.unsat.size();
    }
    // one sweep at inverse temperature with precomputed acceptance table
    // returns accepted flips; stops early (returns -1) if stop_at_zero and E hits 0
    long sweep(State& s, const std::vector<double>& acc, RNG& rng, bool stop_at_zero) const {
        std::uniform_real_distribution<double> U(0, 1);
        long a = 0;
        for (int t = 0; t < n; ++t) {
            int v;
            if (focused) {
                if (s.unsat.empty()) return stop_at_zero ? -1 : a;
                const auto& c = cl[s.unsat[rng() % s.unsat.size()]]; v = std::abs(c[rng() % c.size()]);
            } else v = 1 + rng() % n;
            int d = dE(s, v);
            if (d <= 0 || U(rng) < acc[std::min(d, maxdE)]) { flip(s, v); ++a; }
            if (stop_at_zero && s.E == 0) return -1;
        }
        return a;
    }
    bool verify(const State& s) const {
        for (const auto& c : cl) { bool ok = false; for (int l : c) if (lt(s, l)) { ok = true; break; } if (!ok) return false; }
        return true;
    }
};

// --------------------------------------------------------------- factor model
struct FactorModel {
    uint64_t N = 0; int np = 0, nq = 0; std::vector<std::pair<int, int>> fb; int maxdE = 1;
    struct State { uint64_t p, q; long E; };
    static int popc(u128 x) { return __builtin_popcountll((uint64_t)x) + __builtin_popcountll((uint64_t)(x >> 64)); }
    long energy(uint64_t p, uint64_t q) const { return popc((u128)p * q ^ (u128)N); }
    void setup() {
        for (int b = 1; b <= np - 2; ++b) fb.push_back({0, b});
        for (int b = 1; b <= nq - 2; ++b) fb.push_back({1, b});
        maxdE = np + nq + 2;
    }
    int sweep_size() const { return std::max<int>(1, fb.size()); }
    void init(State& s, RNG& rng) const {
        s.p = 1ull | (1ull << (np - 1)); s.q = 1ull | (1ull << (nq - 1));
        for (auto [w, b] : fb) if (rng() & 1) (w ? s.q : s.p) |= 1ull << b;
        s.E = energy(s.p, s.q);
    }
    long sweep(State& s, const std::vector<double>& acc, RNG& rng, bool stop_at_zero) const {
        std::uniform_real_distribution<double> U(0, 1);
        long a = 0;
        if (fb.empty()) return (stop_at_zero && s.E == 0) ? -1 : 0;
        for (size_t t = 0; t < fb.size(); ++t) {
            auto [w, b] = fb[rng() % fb.size()];
            uint64_t p2 = s.p ^ (w ? 0 : 1ull << b), q2 = s.q ^ (w ? 1ull << b : 0);
            long E2 = energy(p2, q2), d = E2 - s.E;
            if (d <= 0 || U(rng) < acc[std::min<long>(d, maxdE)]) { s.p = p2; s.q = q2; s.E = E2; ++a; }
            if (stop_at_zero && s.E == 0) return -1;
        }
        return a;
    }
    bool verify(const State& s) const { return (u128)s.p * s.q == (u128)N; }
};

// --------------------------------------------------------------- PT driver
struct Opts {
    std::string mode = "solve"; std::vector<double> T; double tlimit = 60; uint64_t seed = 1;
    long sweeps = 100000, burn = 1000;
};

static double tau_int(const std::vector<double>& x) {
    size_t n = x.size(); if (n < 100) return NAN;
    double m = 0; for (double v : x) m += v; m /= n;
    double c0 = 0; for (double v : x) c0 += (v - m) * (v - m); c0 /= n;
    if (c0 == 0) return 0.5;
    double tau = 0.5;
    for (size_t t = 1; t < n / 4; ++t) {
        double c = 0; for (size_t i = 0; i + t < n; ++i) c += (x[i] - m) * (x[i + t] - m);
        c /= (n - t); tau += c / c0;
        if ((double)t >= 6 * tau) break;
    }
    return tau;
}

template <class Model>
int run(const Model& mdl, const Opts& o) {
    RNG rng(o.seed);
    std::uniform_real_distribution<double> U(0, 1);
    const int M = o.T.size();
    std::vector<double> beta(M);
    for (int k = 0; k < M; ++k) beta[k] = 1.0 / o.T[k];
    std::vector<std::vector<double>> acc(M, std::vector<double>(mdl.maxdE + 1));
    for (int k = 0; k < M; ++k) for (int d = 0; d <= mdl.maxdE; ++d) acc[k][d] = std::exp(-beta[k] * d);
    std::vector<typename Model::State> R(M);
    for (auto& s : R) mdl.init(s, rng);
    std::vector<int> slot(M), label(M, 0); // slot[k] = replica at temperature k; label[r]: +1 up, -1 down
    for (int k = 0; k < M; ++k) slot[k] = k;
    std::vector<long> st(M, 0), sa(M, 0), nup(M, 0), ndn(M, 0), n0(M, 0), n1(M, 0), n2(M, 0), ns(M, 0);
    std::vector<double> sE(M, 0), sE2(M, 0);
    long round_trips = 0, flips = 0, acc_flips = 0;
    const bool solve = o.mode == "solve", single = o.mode == "single";
    std::vector<double> series;
    auto t0 = std::chrono::steady_clock::now();
    long sw; int winner = -1; long first_zero = -1;
    const long smax = solve ? (1L << 62) : o.burn + o.sweeps;
    for (sw = 0; sw < smax; ++sw) {
        for (int k = 0; k < M; ++k) {
            long a = mdl.sweep(R[slot[k]], acc[k], rng, solve);
            if (a < 0) { winner = slot[k]; break; }
            acc_flips += a; flips += mdl.sweep_size();
        }
        if (winner >= 0) break;
        if (M > 1) {
            for (int k = sw & 1; k + 1 < M; k += 2) {
                double lr = (beta[k] - beta[k + 1]) * double(R[slot[k]].E - R[slot[k + 1]].E);
                if (sw >= o.burn) ++st[k];
                if (lr >= 0 || U(rng) < std::exp(lr)) { std::swap(slot[k], slot[k + 1]); if (sw >= o.burn) ++sa[k]; }
            }
            // labels: coldest -> up, hottest -> down; a down replica reaching the coldest completes a round trip
            int rc = slot[0], rh = slot[M - 1];
            if (label[rc] == -1 && sw >= o.burn) ++round_trips;
            label[rc] = +1; label[rh] = -1;
        }
        if (sw >= o.burn) {
            for (int k = 0; k < M; ++k) {
                const auto& s = R[slot[k]]; long E = s.E;
                if (label[slot[k]] == 1) ++nup[k]; else if (label[slot[k]] == -1) ++ndn[k];
                sE[k] += E; sE2[k] += double(E) * E; ++ns[k];
                n0[k] += (E == 0); n1[k] += (E == 1); n2[k] += (E == 2);
                if (E == 0 && first_zero < 0) first_zero = sw;
            }
            if (single) series.push_back(R[slot[0]].E);
        }
        if (solve && (sw & 63) == 0 &&
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() > o.tlimit) break;
    }
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (solve) {
        bool ok = winner >= 0;
        if (ok && !mdl.verify(R[winner])) { std::fprintf(stderr, "verification failed\n"); return 2; }
        long minacc_k = -1; double minacc = 2;
        for (int k = 0; k + 1 < M; ++k) if (st[k]) { double a = double(sa[k]) / st[k]; if (a < minacc) { minacc = a; minacc_k = k; } }
        std::printf("result %s time %.4f sweeps %ld work %ld round_trips %ld min_swap_acc %.3f\n",
                    ok ? "SAT" : "TIMEOUT", el, sw, sw * (long)M, round_trips, minacc_k < 0 ? 0.0 : minacc);
        return 0;
    }
    if (single) {
        double m = sE[0] / ns[0], v = sE2[0] / ns[0] - m * m;
        std::printf("single T %.6f meanE %.6f varE %.6f flip_acc %.4f tau_int %.3f p0 %.6f\n",
                    o.T[0], m, v, double(acc_flips) / flips, tau_int(series), double(n0[0]) / ns[0]);
        return 0;
    }
    // tune
    std::printf("# M %d sweeps %ld burn %ld round_trips %ld time %.2f first_zero %ld\n", M, o.sweeps, o.burn, round_trips, el, first_zero);
    std::printf("# k T swap_acc_next f_up meanE varE p0 p1 p2\n");
    for (int k = 0; k < M; ++k) {
        double m = sE[k] / ns[k], v = sE2[k] / ns[k] - m * m;
        double f = (nup[k] + ndn[k]) ? double(nup[k]) / (nup[k] + ndn[k]) : NAN;
        std::printf("%d %.6f %.4f %.4f %.5f %.5f %.6f %.6f %.6f\n", k, o.T[k], (k + 1 < M && st[k]) ? double(sa[k]) / st[k] : NAN,
                    f, m, v, double(n0[k]) / ns[k], double(n1[k]) / ns[k], double(n2[k]) / ns[k]);
    }
    return 0;
}

int main(int argc, char** argv) {
    Opts o; std::string model, cnf, move = "uniform"; uint64_t N = 0; int np = 0, nq = 0;
    double tmin = 0.2, tmax = 2.0; int M = 16; std::string temps;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i]; auto nx = [&]() { return std::string(argv[++i]); };
        if (a == "--model") model = nx();
        else if (a == "--N") N = std::stoull(nx());
        else if (a == "--np") np = std::stoi(nx());
        else if (a == "--nq") nq = std::stoi(nx());
        else if (a == "--move") move = nx();
        else if (a == "--mode") o.mode = nx();
        else if (a == "--temps") temps = nx();
        else if (a == "--tmin") tmin = std::stod(nx());
        else if (a == "--tmax") tmax = std::stod(nx());
        else if (a == "--replicas") M = std::stoi(nx());
        else if (a == "--time") o.tlimit = std::stod(nx());
        else if (a == "--seed") o.seed = std::stoull(nx());
        else if (a == "--sweeps") o.sweeps = std::stol(nx());
        else if (a == "--burn") o.burn = std::stol(nx());
        else cnf = a;
    }
    if (!temps.empty()) {
        std::stringstream ss(temps); std::string t; while (std::getline(ss, t, ',')) o.T.push_back(std::stod(t));
        std::sort(o.T.begin(), o.T.end());
    } else if (o.mode == "single") o.T = {tmin};
    else for (int k = 0; k < M; ++k) o.T.push_back(M == 1 ? tmin : tmin * std::pow(tmax / tmin, double(k) / (M - 1)));
    if (o.mode == "solve") o.burn = 0;
    if (model == "cnf") { CnfModel m; m.load(cnf); m.focused = (move == "focused"); return run(m, o); }
    if (model == "factor") {
        if (!N || np < 2 || nq < 2 || np + nq > 62) { std::fprintf(stderr, "factor model needs --N --np --nq\n"); return 1; }
        FactorModel m; m.N = N; m.np = np; m.nq = nq; m.setup(); return run(m, o);
    }
    std::fprintf(stderr, "--model cnf <file> | factor --N --np --nq\n"); return 1;
}
