// Exact density of states g(E), E = popcount(p*q XOR N), over all free bits
// (MSB and LSB of p, q fixed to 1).  Prints "E g(E)" lines.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
typedef unsigned __int128 u128;
int main(int argc, char** argv) {
    uint64_t N = strtoull(argv[1], 0, 10); int np = atoi(argv[2]), nq = atoi(argv[3]);
    int fp = np - 2, fq = nq - 2;
    std::vector<unsigned long long> g(130, 0);
    for (uint64_t a = 0; a < (1ull << fp); ++a) {
        uint64_t p = 1ull | (a << 1) | (1ull << (np - 1));
        u128 P = p;
        for (uint64_t b = 0; b < (1ull << fq); ++b) {
            uint64_t q = 1ull | (b << 1) | (1ull << (nq - 1));
            u128 x = P * q ^ (u128)N;
            int e = __builtin_popcountll((uint64_t)x) + __builtin_popcountll((uint64_t)(x >> 64));
            ++g[e];
        }
    }
    for (int e = 0; e < 130; ++e) if (g[e]) printf("%d %llu\n", e, g[e]);
}
