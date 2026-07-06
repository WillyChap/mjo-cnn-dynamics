#!/usr/bin/env python
"""Pairwise comparison of experiment E:W bootstrap distributions (CNN vs MEAN vs CNTRL)."""
import numpy as np
z = np.load("/glade/work/wchapman/DA_ML/CESML_AI/Paper_Figures/Wang/Paper_WCD/figures/bootstrap_ew_ratio.npz", allow_pickle=True)
point = z["point"].item()
B = {k: z[k] for k in ["ERA5", "CNTRL", "MEAN", "CNN"]}

print("=== Point estimates and 90% block-bootstrap CIs (5th, 95th pct) ===")
for k in ["ERA5", "CNTRL", "MEAN", "CNN"]:
    b = B[k]
    lo, hi = np.percentile(b, [5, 95])
    print(f"  {k:6s} point={point[k]:.2f}  90% CI=[{lo:.2f}, {hi:.2f}]  mean={b.mean():.2f} sd={b.std(ddof=1):.2f}")

def compare(a, bb):
    """P(ratio_a > ratio_b) using independent bootstrap draws; two-sided p ~ 2*min(p,1-p)."""
    n = min(len(B[a]), len(B[bb]))
    ra, rb = B[a][:n], B[bb][:n]
    pg = np.mean(ra > rb)
    diff = ra - rb
    return pg, diff.mean(), np.percentile(diff, [5, 95])

print("\n=== Pairwise (independent-draw) comparisons ===")
for a, bb in [("CNN", "MEAN"), ("CNN", "CNTRL"), ("MEAN", "CNTRL"), ("CNN", "ERA5")]:
    pg, dm, dci = compare(a, bb)
    p2 = 2 * min(pg, 1 - pg)
    print(f"  {a} vs {bb}: P({a}>{bb})={pg:.3f}  diff mean={dm:+.3f} 90%CI=[{dci[0]:+.3f},{dci[1]:+.3f}]  two-sided p~{p2:.3f}")
