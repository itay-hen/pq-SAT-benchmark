# Re-run of the solver-scaling study (Fig. 3)

Prepared Nov 2026.  This replaces the March 2026 runs in `results/` and
`results2/`, which had the following problems:

* `results/` (used for Fig. 3 of paper_v16): 90 tasks, of which 67 finished.
  Instances were drawn as 3- to 7-digit *decimal* primes, so the counts per
  bit-length were 1 to 13 (single instances at d = 8, 9, 11, 12, 16) and
  d = 15, 18, 21, 22, 25 never occurred.  The 23 unfinished tasks (d = 27 to 34)
  were ended by the 48 h Slurm limit; three of them were at d = 27, so the
  plotted d = 27 median used only the three fastest of six runs.
* `results2/`: CaDiCaL returned ERROR after 0.007 s on every instance (a
  broken binary or path), and 47 tasks did not finish.
* Both runs used CNF files from the generator before the preprocessing fix,
  which contain the forced-zero high columns.
* Solver versions were not recorded.  The source trees in this folder are
  Kissat 4.0.4 and CaDiCaL 3.0.0 (the paper said Kissat 3.0 / CaDiCaL 1.5).

## Design of the re-run

* p and q are both d-bit primes (p != q), the same number of pairs at every d.
* Instances come from the fixed generator (`../../../pq python/factorization_sat.py`).
* One Slurm task per (instance, solver): the solvers no longer share a budget.
* An explicit per-solver timeout (default 40 h) below the Slurm limit (44 h),
  so slow runs end as TIMEOUT and are recorded, not silently killed.
* Every log records the solver versions and CPU model.
* Runs that time out or are killed are treated as censored: a median is
  reported only when more than half of the runs at that d finished, and only
  those points enter the fit.

## Steps

1. On the cluster, in the folder that holds `sat_benchmark`, `kissat` and
   `cadical` (e.g. `/project2/hen_48/hen/pq/cppBenchmark`), copy this folder's
   scripts and the fixed `factorization_sat.py`.
2. Check the solver binaries first; CaDiCaL failed in `results2/`:

       ./kissat --version; ./cadical --version
       ./sat_benchmark --timeout 60 cnf/d10/<any>.cnf     # after step 3

3. Build the instance set and the task list (d = 10..28, 12 pairs per d,
   two solvers, 456 tasks; `pqSAT_rerun.job` is already set to --array=0-455%60):

       python3 make_benchmark_set.py --dmin 10 --dmax 28 --per-d 12 --seed 2026 \
           --timeout 144000 --generator .

   With seed 2026 and 12 pairs per d, the instances for d = 10..21 are
   byte-identical to those of the local run in local_run_2026-11-17/, so the
   two machines can be compared directly.

   It prints the `--array` range to put in `pqSAT_rerun.job`.
4. `mkdir -p logs && sbatch pqSAT_rerun.job`
5. When the array has finished (unfinished tasks are fine; they count as censored):

       python3 analyze_rerun.py --min-n 10

   This writes `rerun_summary.tsv`, `rerun_fits.tsv` and
   `median_runtime_vs_bits_rerun.{pdf,png}`, and prints the solver versions seen
   in the logs.  Use those numbers and the figure in the Benchmark section of the
   paper (paper_v17 marks the places with TODO(author)).

## Decision (Nov 17, 2026)

Run the extension to d = 28 on the cluster before submitting the paper,
and use the cluster data for Fig. 3 (all of d = 10..28 from one machine).
Reason: the paper's statement that d of about 18 to 28 is hard but practical
to benchmark currently rests on extrapolating the local data, which stop at
medians of about two minutes (d = 21); only the biased March data reach
d = 27.  The cost is modest (estimate below).  Keep the local run as a
cross-check of the d <= 21 points.

## Cost estimate

From the local fits (Kissat 0.315, CaDiCaL 0.341 log10 s per bit, Xeon @
2.1 GHz), the medians at d = 28 are about 1.6e4 s (Kissat) and 4.4e4 s
(CaDiCaL), and the whole d = 10..28 set is about 560 core-hours.  With 60
tasks running at once the wall time is set by the slowest single runs,
about one to two days.  The 40 h per-solver timeout leaves room for the
upper quartile at d = 28; anything slower is recorded as censored.

## Earlier estimate

With the March data (about a factor of 2 per bit, roughly 10^3 to 10^4 s at
d = 26 to 27), d = 28 needs of order 10^4 s per run or more and most of the core-hours go
to the top two or three values of d.  Lower `--dmax` or `--per-d` to cut the
cost; keep `--per-d` at 10 or more so the medians mean something.

## Files

* `make_benchmark_set.py` builds `cnf/`, `manifest.tsv`, `commands.txt`.
* `pqSAT_rerun.job` is the Slurm array script.
* `analyze_rerun.py` produces the summary, fits and figure.
* `sat_benchmark_single.cpp` is unchanged (compile with
  `g++ -std=c++17 -O2 -o sat_benchmark sat_benchmark_single.cpp`).

The kit was tested end to end with this harness on d = 8..15, including
simulated timeouts, killed tasks and a solver error.

## Local run used for Fig. 3 of paper_v17 (Nov 17, 2026)

`local_run_2026-11-17/` holds a complete run of this kit on a 2-core cloud
VM (Intel Xeon @ 2.1 GHz, two runs at a time): d = 10..21, 12 pairs per d,
seed 2026, timeout 1800 s, Kissat 4.0.4 and CaDiCaL 3.0.0 built from the
release sources.  All 288 runs were solved and passed the factorization
check (SAT F).  Fits: Kissat 0.315, CaDiCaL 0.341 log10(s) per bit.
It contains manifest.tsv, commands.txt, results/, logs/, rerun_summary.tsv,
rerun_fits.tsv and the figure.  The CNF files are not included; rerunning
make_benchmark_set.py with the same arguments regenerates them exactly:

    python3 make_benchmark_set.py --dmin 10 --dmax 21 --per-d 12 --seed 2026 \
        --timeout 1800 --generator .

The cluster run described above would extend the same protocol to d of
about 28.
