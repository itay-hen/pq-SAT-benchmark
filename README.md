**pq-SAT-benchmark**
Planted-solution SAT and Ising benchmark instances derived from integer factorization.

This repository generates structured satisfiability and optimization benchmark instances from products of two primes. Given primes p and q, the generator encodes the arithmetic constraints of
N = p q
as a Boolean constraint system and exports benchmark instances whose satisfying assignments correspond to valid factorizations of N. The known factor pair (p, q), up to the trivial swap symmetry when the two factor registers have the same width, provides a planted ground truth for solver validation.
The same residual constraint system can be exported as DIMACS CNF for SAT solvers or compiled into a quadratic Ising Hamiltonian for classical and quantum optimization benchmarks.

**Why this benchmark family?**

Many benchmark families are either random and scalable but lack a known solution, or structured and realistic but difficult to scale systematically. This project is intended to sit between those regimes:
planted ground truth: the factor pair is known by construction;
single-parameter scaling: difficulty is primarily controlled by the bit length of the prime factors;
structured, non-random constraints: the instances inherit the deterministic structure of binary multiplication;
carry structure: carries propagate column by column through the multiplication circuit and correlate the bits of p and q across the output columns;
dual SAT/Ising output: the same underlying instance can be used with SAT solvers, classical Ising optimizers, or quantum optimization hardware after embedding.
These instances are not meant as a route to factoring large integers. They are a benchmark generator: a controlled way to produce structured, planted, verifiable instances for comparing solvers and optimization platforms.

**Construction overview**
For input primes p and q, the pipeline is:
p, q
  |
  v
binary multiplication constraints
  |-- partial-product AND constraints
  |-- half-adder XOR/AND carry constraints
  |-- output-bit pinning constraints from N = p q
  |
  v
Boolean preprocessing
  |-- pin propagation
  |-- AND/XOR simplification
  |-- equivalence and complement merging
  |-- (rules applied with a worklist until a fixed point is reached)
  |
  v
residual constraint system
  |-- DIMACS CNF export
  |-- quadratic Ising export
The raw constraint system uses three elementary relation types:
c = a AND b for partial products and carries;
c = a XOR b for column sums;
pinning constraints fixing the final column outputs to the known bits of N.
The preprocessing stage repeatedly propagates fixed values, simplifies Boolean relations, merges equivalent or complementary variables, and removes constraints already implied by the known product bits.

**Scaling**
For two d-bit prime factors, the number of half-adder contractions before preprocessing is
C(d) = d^2 (d - 1)^2 / 2.
As a result, the leading-order number of Boolean variables and constraints scales as
O(d^4).
For asymmetric factors with bit lengths n_p and n_q, the corresponding contraction count is
C(n_p, n_q) = n_p n_q (n_p - 1)(n_q - 1) / 2.
The quartic count comes from using half-adders only: a column with m entries sends m - 1 carries to the next column, so the carries pile up into a long tail of d^2 - 2d columns above the most significant bit of N.

Every variable in a column k >= n_p + n_q is forced to zero (a 1 there would make the product exceed 2^(n_p+n_q)), and preprocessing removes all of those columns. What remains lives in the 2d output columns, which contain exactly

C_out(d) = d^3 - d^2 - d

contractions (for general widths, C_out = [n_p n_q (n_p + n_q) - n_p(n_p+1) - n_q(n_q+1)] / 2). The residual instances that the generator writes therefore grow as d^3; in the benchmark set of benchmark/local_run_2026-11-17 (d = 10..21, 12 pairs per d) they have 1.34-1.84 d^3 free variables and 4.70-6.43 d^3 CNF clauses, both ratios rising slowly with d (medians 1.50 d^3 and 5.26 d^3 at d = 10, 1.80 d^3 and 6.28 d^3 at d = 21).

**Changes in the November 2026 release**

* Preprocessing fix: the AND rule "output = 0 with equal inputs => input = 0" was missing, so the forced-zero high columns were left in the output. Residual instances are now about 3x smaller at d = 8 and 7x smaller at d = 16 (e.g. 2617 x 2267: 7,320 + 7,200 residual AND/XOR relations before, 1,434 + 1,315 now). The satisfying assignments are unchanged.
* New AND rule: an output equal to the complement of an input forces input = 1 and output = other input = 0.
* Preprocessing runs as a worklist to a fixed point (about 1 s at d = 19); a contradiction now raises an error instead of being ignored.
* `generate_instances.py`: d is now the bit length of each prime, as in the paper; use `--digits` for the old decimal-digit behaviour.
* Files renamed to `factorization_sat.py` and `generate_instances.py` so that the imports between scripts work.
* The sample instances in `instances/` were regenerated with the fixed code.
* Ising TSV files store fields and couplings only; with the offset E0 from the file name, H = E0 + sum h s + sum J s s is 0 at the planted state, so the stored Hamiltonian has ground-state energy -E0.

**Outputs**
The intended outputs are:
DIMACS CNF instances for SAT solvers;
metadata recording the factors, product, bit lengths, seed, and planted assignment;
Ising Hamiltonians of the form
H(s) = E0 + sum_i h_i s_i + sum_{i<j} J_ij s_i s_j,
with known planted ground state (H = 0 at the planted state, i.e. the fields and couplings alone give -E0);
benchmark data, such as seeds, solver outputs, and timing information.

## Ising compilation

By default each residual AND relation gets a quadratic AND gadget and each
XOR relation an XOR gadget with one auxiliary spin; this is the released
benchmark.  With `--gadget combined` every XOR relation whose half-adder AND
partner is also present is compiled with it as
`(1/2)(s_a + s_b - s_s - 2 s_t - 1)^2`, which needs no auxiliary spin.  The
problem, the planted solution and the CNF are unchanged; the Ising model has
fewer spins (28 instead of 40 for 11 x 13).

`analysis/` holds two scripts (they need `pip install python-sat`) that
certify the ground-state degeneracy and the exact gap of an instance:
`certify_ground_states.py 7x13 ...` and `exact_gap.py 7x13 ...`.

## Repository layout (Nov 2026)

| Path | Contents |
|---|---|
| `factorization_sat.py`, `generate_instances.py` | Instance generator (CNF and Ising), including `--gadget combined` |
| `thorough_test.py`, `verify_factorization_cnf.py`, `analysis/` | Tests, CNF verification, ground-state and gap certification |
| `benchmark/` | SAT benchmark driver (`sat_benchmark_single.cpp`), benchmark-set builder, Slurm job, analysis, and raw logs of the 2026-11-17 local run (Fig. 3) |
| `dwave/` | D-Wave instance generator, the annealing-study instances as produced by the released generator, and `fit_hardware_counts.py` (binomial fit and TTS of per-setting success counts). The D-Wave hardware data (success and read counts, run identifiers) are not in this repository. |
| `paper_support/` | Scripts and data for Fig. 1 and the old/new generator timing check |
| `pt_study/`, `pt_tuned/` | Parallel-tempering study (sources, ladders, results) |
| `solvers/kissat`, `solvers/cadical` | Git submodules pinned to kissat rel-4.0.4 and CaDiCaL rel-3.0.0 |

Clone with `git clone --recursive`, or run `git submodule update --init` after cloning.
Scripts in `benchmark/` and `paper_support/` import `factorization_sat`: `benchmark/make_benchmark_set.py` and `dwave/make_dwave_instances.py` find it in the repository root by default (`--generator`), and `paper_support/fig1_data.py` reads the root from `PQ_REPO`.
