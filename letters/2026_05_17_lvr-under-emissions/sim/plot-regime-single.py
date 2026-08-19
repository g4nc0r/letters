#!/usr/bin/env python3
"""
Single-panel regime figure for "Loss-Versus-Rebalancing Under Emissions".

Re-plots the left panel of regime_mc.png, the break-even ratio against
volatility, as the standalone figure the note carries. Reads
regime_mc_results.json (produced by regime-mc.py); it does not re-run the
Monte Carlo, so the figure is numerically identical to the data behind the
note's fitted-exponents table.

Output: ../figures/emission-lvr-regime.png, the figure the note carries as its
Figure 1. Run regime-mc.py first if regime_mc_results.json is absent.
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SIM = Path(__file__).resolve().parent
RESULTS = SIM / "regime_mc_results.json"
OUT = SIM.parent / "figures" / "emission-lvr-regime.png"


def main():
    if not RESULTS.exists():
        sys.exit(f"missing {RESULTS}; run regime-mc.py first")

    with open(RESULTS) as fh:
        res = json.load(fh)

    managed = res["managed"]
    static = res["static"]
    fits = res["fits"]
    f_over_eps = res["params"]["f_over_eps"]

    sig = [r["sigma"] for r in managed]
    xs = np.array(sig)
    fo = [r["f_over_lvr"] for r in managed]
    eo = [r["em_over_lvr"] for r in managed]
    eo_static = [r["em_over_lvr"] for r in static]

    fF = fits["managed"]["f_over_lvr_vs_sigma"]
    fE = fits["managed"]["em_over_lvr_vs_sigma"]
    fS = fits["static"]["em_over_lvr_vs_sigma"]

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    ax.plot(sig, fo, "o-", color="C0",
            label=f"F/LVR, fee-side  (slope {fF['slope']:+.2f}, predict 0)")
    ax.plot(sig, eo, "s-", color="C3",
            label=f"em/LVR, emission-side  (slope {fE['slope']:+.2f}, predict -2)")
    ax.plot(sig, eo_static, "s--", color="C3", alpha=0.45,
            label=f"em/LVR, static position  (slope {fS['slope']:+.2f})")
    anchor = np.median([e * s ** 2 for e, s in zip(eo, sig)])
    ax.plot(xs, anchor * xs ** -2.0, color="0.45", lw=1.0,
            label=r"reference slope $\sigma^{-2}$")
    ax.axhline(f_over_eps, color="C0", ls=":", lw=1.0,
               label=r"reference $f/(-\ln(1-f))$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"annualised volatility $\sigma$")
    ax.set_ylabel("break-even ratio (income / LVR)")
    ax.set_title("(simulation)  break-even ratio against volatility")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, which="major", alpha=0.25)

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    sys.exit(main())
