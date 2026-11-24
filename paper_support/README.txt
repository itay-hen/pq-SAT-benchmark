Support data for the paper (Fig. 1 and the old/new generator check)

fig1_data.py / fig1_plot.py / fig1.json: Fig. 1 (run fig1_data.py with PQ_REPO=<path to pq-SAT-benchmark>, then fig1_plot.py). The prime draws use a seeded random.Random (seed 2026) with sympy.isprime, so fig1.json is reproducible; it was regenerated in Nov 2026 when the unseeded sympy.randprime draws were replaced.
oldnew_timing.jsonl: spot check of solver runtimes on CNFs from the old and fixed generator (same prime pairs, 9 per d, d=10-19 complete, d=20 partial). Fields: d, p, q, gen (old/new), solver (cadical153 = CaDiCaL 1.5.3, kissat404 = Kissat 4.0.4, via PySAT), vars, cls, t (s), sat.
