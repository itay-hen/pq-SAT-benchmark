**pq-SAT-benchmark**
Planted-solution SAT and Ising benchmark instances derived from integer factorization.
This repository generates structured satisfiability and optimization benchmark instances from products of two primes. Given primes p and q, the generator encodes the arithmetic constraints of
N = p q
as a Boolean constraint system and exports benchmark instances whose satisfying assignments correspond to valid factorizations of N. The known factor pair (p, q), up to the trivial swap symmetry when the two factor registers have the same width, provides a planted ground truth for solver validation.
The same residual constraint system can be exported as DIMACS CNF for SAT solvers or compiled into a quadratic Ising Hamiltonian for classical and quantum optimization benchmarks.
Why this benchmark family?
Many benchmark families are either random and scalable but lack a known solution, or structured and realistic but difficult to scale systematically. This project is intended to sit between those regimes:
planted ground truth: the factor pair is known by construction;
single-parameter scaling: difficulty is primarily controlled by the bit length of the prime factors;
structured, non-random constraints: the instances inherit the deterministic structure of binary multiplication;
carry-induced long-range correlations: carries propagate through the multiplication circuit and correlate variables across long column chains;
dual SAT/Ising output: the same underlying instance can be used with SAT solvers, classical Ising optimizers, or quantum optimization hardware after embedding.
These instances are not meant as a route to factoring large integers. They are a benchmark generator: a controlled way to produce structured, planted, verifiable instances for comparing solvers and optimization platforms.
Construction overview
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
  |-- cross-clause inference
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
Scaling
For two d-bit prime factors, the number of half-adder contractions before preprocessing is
C(d) = d^2 (d - 1)^2 / 2.
As a result, the leading-order number of Boolean variables and constraints scales as
O(d^4).
For asymmetric factors with bit lengths n_p and n_q, the corresponding contraction count is
C(n_p, n_q) = n_p n_q (n_p - 1)(n_q - 1) / 2.
The quartic scaling is caused by carry cascading: contractions create carries, carries enter later columns, and those later columns require further contractions.
Outputs
The intended outputs are:
DIMACS CNF instances for SAT solvers;
metadata recording the factors, product, bit lengths, seed, and planted assignment;
Ising Hamiltonians of the form
H(s) = E0 + sum_i h_i s_i + sum_{i<j} J_ij s_i s_j,
with known planted ground state;
benchmark data, such as seeds, solver outputs, and timing information.
