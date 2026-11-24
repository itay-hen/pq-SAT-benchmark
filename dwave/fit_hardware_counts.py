"""Binomial (logistic) fit of D-Wave success counts, including k = 0.

Input CSV columns:
    instance,spins,anneal_us,successes,reads[,readout_us,delay_us,programming_us]
(one row per instance, annealing time and run setting; see the Quantum
annealing section of the paper).  The original per-setting hardware counts
are not in this repository, so this script is the analysis to apply once
counts are recovered or new reads are taken; it is not a fit to them.

Model: logit p_ia = b_a + b_s (spins_i - mean spins), one intercept per
annealing time a and a common slope in the spin count, fitted by maximum
likelihood.  For each row it also reports the number of reads for 99%
success, r99 = ceil(ln 0.01 / ln(1 - p)), and
TTS99 = t_prog + r99 (t_a + t_read + t_delay).  When the fitted p is 0
(numerically), r99 and TTS99 are infinite and are reported as null.

Usage:
    python3 fit_hardware_counts.py counts.csv --out fit.json
"""
from pathlib import Path
import argparse
import csv
import json
import math

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit


def reads_for_99(p):
    """Reads needed for 99% success at per-read probability p; None if p == 0."""
    if p >= 1:
        return 1
    if p <= 0:
        return None
    denom = math.log1p(-p)
    if denom == 0.0:          # p below double resolution of 1 - p
        return None
    return int(math.ceil(math.log(0.01) / denom))


def fit(rows):
    times = sorted({float(r['anneal_us']) for r in rows})
    ix = {t: i for i, t in enumerate(times)}
    n = np.array([int(r['reads']) for r in rows])
    k = np.array([int(r['successes']) for r in rows])
    x = np.array([float(r['spins']) for r in rows])
    x0 = x.mean()
    X = np.zeros((len(rows), len(times) + 1))
    X[np.arange(len(rows)), [ix[float(r['anneal_us'])] for r in rows]] = 1
    X[:, -1] = x - x0

    def objective(th):
        z = X @ th
        p = expit(z)
        return np.sum(n * np.logaddexp(0, z) - k * z), X.T @ (n * p - k)

    result = minimize(objective, np.r_[np.repeat(-2., len(times)), 0.],
                      jac=True, method='BFGS')
    p = expit(X @ result.x)
    info = X.T @ ((n * p * (1 - p))[:, None] * X)

    tts = []
    for r, pr in zip(rows, p):
        t = float(r['anneal_us'])
        ov = float(r.get('readout_us') or 0) + float(r.get('delay_us') or 0)
        reps = reads_for_99(float(pr))
        total = None if reps is None else float(r.get('programming_us') or 0) + reps * (t + ov)
        tts.append(dict(instance=r['instance'], anneal_us=t, p_fitted=float(pr),
                        repetitions99=reps, TTS99_us=total))

    optimal = {}
    for inst in dict.fromkeys(r['instance'] for r in rows):
        finite = [a for a in tts if a['instance'] == inst and a['TTS99_us'] is not None]
        optimal[inst] = min(finite, key=lambda a: a['TTS99_us']) if finite else None

    return dict(
        converged=bool(result.success), message=str(result.message),
        center_spins=float(x0), times_us=times, parameters=result.x.tolist(),
        SE=np.sqrt(np.diag(np.linalg.pinv(info))).tolist(),
        optimal=optimal, all_TTS99=tts,
        note='Model: logit p_ia = intercept_a + slope*(spins_i-center). Assumes '
             'independent reads at each setting. This common-slope model must be '
             'checked; repeated embeddings/gauges need cluster-level uncertainty. '
             'Choosing the best time on the same small data is optimistic; '
             'validate it on held-out reads. TTS99 is null where the fitted p is 0.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    with open(a.csv) as fh:
        rows = list(csv.DictReader(fh))
    r = fit(rows)
    Path(a.out).write_text(json.dumps(r, indent=2))
    print(json.dumps(r, indent=2))
