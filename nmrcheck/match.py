"""Peak list comparison."""
import math

import numpy as np
from scipy.optimize import linear_sum_assignment

# ppm scale on which a deviation counts as large, and the window inside which
# two peaks of different species are reported as overlapping
SCALE = {"13C": 3.0, "1H": 0.3}
OVERLAP = {"13C": 0.5, "1H": 0.05}


def centres(nucleus, peaks):
    """Flatten stored peaks to a list of ppm values."""
    if nucleus == "1H":
        return [p[0] for p in peaks]
    return list(peaks)


def align(predicted, observed):
    """Pair predicted shifts with observed peaks.

    Each observed peak is used once where possible; predicted shifts left
    over (accidental overlap is common in real spectra) share the nearest
    observed peak. Returns (pairs, unexplained) where pairs is a list of
    (predicted_index, observed_index) and unexplained lists observed indices
    that no predicted shift landed on.
    """
    if not predicted or not observed:
        return [], list(range(len(observed)))
    pred = np.asarray(predicted, dtype=float)
    obs = np.asarray(observed, dtype=float)
    cost = np.abs(pred[:, None] - obs[None, :])
    rows, cols = linear_sum_assignment(cost)
    pairs = dict(zip(rows.tolist(), cols.tolist()))
    for i in range(len(pred)):
        if i not in pairs:
            pairs[i] = int(np.argmin(cost[i]))
    used = set(pairs.values())
    unexplained = [j for j in range(len(obs)) if j not in used]
    return sorted(pairs.items()), unexplained


def agreement(nucleus, predicted, observed):
    """Score in 0..1 for how well an observed list fits predicted shifts."""
    pairs, unexplained = align(predicted, observed)
    if not pairs:
        return {"score": 0.0, "mae": None, "max_dev": None,
                "unexplained": len(unexplained)}
    devs = [abs(predicted[i] - observed[j]) for i, j in pairs]
    mae = sum(devs) / len(devs)
    covered = 1.0 - len(unexplained) / len(observed)
    return {
        "score": round(math.exp(-mae / SCALE[nucleus]) * covered, 3),
        "mae": round(mae, 3),
        "max_dev": round(max(devs), 3),
        "unexplained": len(unexplained),
    }


def overlaps(nucleus, a, b):
    """Pairs (x, y) from two shift lists that fall inside the overlap window."""
    tol = OVERLAP[nucleus]
    return [(x, y) for x in a for y in b if abs(x - y) <= tol]
