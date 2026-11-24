# Parallel tempering vs Kissat on the pq factoring instances (Nov 2026)

## Code (src/)
* `pt_cnf.cpp`: parallel tempering over assignments of the residual CNF,
  energy = number of violated clauses.  `--move uniform` (Metropolis flip of a
  uniformly chosen variable, detailed balance) or `--move focused` (variable
  from a random violated clause; not detailed balance).  `--mode sample` writes
  per-temperature histograms over all 2^n states (n <= 24).
* `pt_factor.cpp`: parallel tempering over the factor bits only (MSB and LSB
  of p and q fixed to 1), energy = popcount(p*q XOR N).
* Build: `g++ -O3 -march=native -std=c++17 -o pt_cnf pt_cnf.cpp` (same for pt_factor).

## Detailed balance (detailed_balance/)
`python3 check_db.py f5x5.cnf f5x5.uniform.hist` compares sampled and exact
Boltzmann distributions.  5x5 residual CNF (15 vars, 52 clauses) and a random
3-CNF (12 vars, 50 clauses), 8 replicas, T = 0.3 ... 3, 2e6 sweeps.
* uniform moves: energy-distribution TV <= 0.0023, <E> to 2-3 decimals; state
  TV equal to the sampling-noise floor of 2e6 independent exact draws
  (0.0223 vs 0.0222 at T = 0.8, 0.0472 vs 0.0471 at T = 3).
* focused moves (negative control): TV 0.2-1.0; a solution is absorbing.
* pt_factor (N = 53*59): state TV <= 0.006, <E> to 3 decimals.

## Benchmark (benchmark/)
Instances: the Fig. 3 local-run set for d = 10-18 (12 symmetric pairs per d)
plus d = 5-9 built with the same script and seed (10 pairs at d = 5, all that
exist; 12 at d = 6-9).  One run per instance and method, seed 1, 300 s limit,
2 jobs on a 2-core machine.  A method was dropped after the first d with a
censored median.  Kissat 4.0.4 rerun here on the same machine; its medians are
1.1-1.5x those of the Fig. 3 run, slope 0.296 vs 0.305 decades/bit.

PT protocol, fixed for all instances: uniform moves with a geometric ladder
T = 0.15 ... 1.2 and ceil(48 sqrt(m/4698)) replicas (m = clauses); focused
moves with 16 replicas, T = 0.1 ... 0.6; factor PT with 12 replicas,
T = 0.3 ... 4.  Random walk = factor PT with 1 replica at T = 1e12.

Timing: Kissat times are process wall time (includes parsing and start-up,
about 3 ms); PT times are measured inside the program after parsing.

See summary.tsv (medians, quartiles) and fits.tsv (log10 median vs d).
