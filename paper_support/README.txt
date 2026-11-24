Support data for paper_v17

fig1_data.py / fig1_plot.py / fig1.json: Fig. 1 (run fig1_data.py with PQ_REPO=<path to pq-SAT-benchmark>, then fig1_plot.py).
oldnew_timing.jsonl: spot check of solver runtimes on CNFs from the old and fixed generator (same prime pairs, 9 per d, d=10-19 complete, d=20 partial). Fields: d, p, q, gen (old/new), solver (cadical153 = CaDiCaL 1.5.3, kissat404 = Kissat 4.0.4, via PySAT), vars, cls, t (s), sat.
