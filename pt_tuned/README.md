# PT with per-d optimized temperature ladders (Nov 19, 2026)

Follow-up to pt_study_2026-11-19.  Same test instances (d = 5-18, 12 per d,
10 at d = 5), same 300 s limit, same machine (2 cores, 2 jobs at a time).

## Code (src/)
* `pt.cpp` - one program, two models (`--model cnf <file>` / `--model factor
  --N --np --nq`), modes `solve`, `tune` (flow fraction f(T), swap acceptance,
  round trips, <E>, var E, ground-state occupation p0 per temperature) and
  `single` (tau_int of E).  Checked against exact enumeration on the 5x5
  residual CNF: <E>, var E and P0 agree to 3-4 digits at all 8 temperatures.
* `dos_factor.cpp` - exact density of states g(E) of the factor-bit model.
* `tune_ladder.py`, `tune_branch.py` - the tuning procedure below.
* `production.py`, `compare.py` - test-set runs and comparison.

## Tuning procedure (per branch, per d, on 4 training instances per d that
## are not in the test set; d = 5 reuses test pairs - only 10 exist)
1. T_max: lowest T above which every hotter grid temperature has
   single-replica tau_int(E) <= 1 sweep.
2. Shape: feedback optimization of Katzgraber, Trebst, Huse & Troyer (2006),
   equal sqrt(df) per interval, damped, 4-5 iterations, f averaged over the
   training instances.
3. For each T_min candidate, M grown (x1.4) until every adjacent swap
   acceptance >= 0.2; then p0(T_min) and training time to solution.
4. T_min adequacy:
   factor - exact equilibrium P0(T_min) >= 0.5 on every training instance
            (from g(E)); among adequate candidates, lowest training work.
   cnf    - measured p0(T_min) >= 0.5 with >= 50 round trips; if none
            qualifies, chosen by training work and flagged.
5. M = lowest cost per round trip (M x sweeps / round trips) among
   {0.7, 1, 1.4, 2} x M0 with min acceptance >= 0.15.
Final ladders: tuning/ladders/*.json; summary: tuning/ladder_summary.tsv.

## Main results (results/, figures/)
* Factor bits: tuned ladders (M = 3-8, T_max = 1.76-2.13) give a geometric-
  mean speedup of 1.4x over the untuned ladder, not significant (paired
  Wilcoxon p = 0.09, d = 13-18).  Slope 0.59 decades/bit vs 0.55 untuned and
  0.60 for an unguided random walk.  The cold end equilibrates for d <= 14
  (measured p0 matches the exact value) and not for d >= 15.
* Residual CNF, uniform moves: 1.7x geometric-mean speedup (p = 0.02); d = 8
  now 10/12 solved, median 200 s (untuned: 4/12).  Slope 1.42 decades/bit,
  unchanged.  M grows from 12 (d = 5) to 56 (d = 8).  The cold end is not
  equilibrated from d = 7 on, so T_min adequacy is not verifiable there.
* Kissat on the same machine: 0.30 decades/bit; at d = 8 it is ~10^4 times
  faster than tuned CNF PT.  Untuned focused moves (no detailed balance)
  remain ~50x faster than tuned uniform-move PT at d = 8.
d = 9 was not tuned for the CNF branch (training TTS fully censored at 60 s
and p0 unmeasurable, so neither selection criterion can be evaluated).
