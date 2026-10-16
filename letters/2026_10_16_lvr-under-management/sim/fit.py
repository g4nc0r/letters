#!/usr/bin/env python3
"""
Rate-law fit for "Loss-Versus-Rebalancing Under Management"
(../lvr-under-management.pdf).

Reads pool_days.csv, one row per pool-day:

  pool_id         pool index, 1 to 30
  width_class     tight, medium, wide or very_wide
  width_spacings  range width in tick-spacings
  period          before or after for the two weeks either side of the
                  30 July trigger change, otherwise empty
  rebalances      rebalances on the day
  f_stock         realised loss / trailing-24h stock estimate
  f_flow          realised loss / flow estimate
  mean_abs_u      mean firing depth |u|

Fits the one-variable rate law and the two-term fit with percentile
bootstrap intervals, prints the tables the note carries and writes
../figures/f_vs_rate.pdf. Deterministic under SEED.
"""
import csv
import math

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED = 20261016
B = 2000
WIDTHS = ["tight", "medium", "wide", "very_wide"]
BUCKETS = [(0, 5, "<5"), (5, 20, "5-20"), (20, 50, "20-50"),
           (50, 150, "50-150"), (150, math.inf, ">=150")]


def num(s):
    return float(s) if s else np.nan


def fit(X, y):
    return np.linalg.lstsq(X, y, rcond=None)[0]


def bootstrap(X, y, rng):
    n = len(y)
    return np.array([fit(X[i], y[i])
                     for i in (rng.integers(0, n, n) for _ in range(B))])


def ci(draws):
    return np.percentile(draws, 2.5), np.percentile(draws, 97.5)


rng = np.random.default_rng(SEED)
rows = list(csv.DictReader(open("pool_days.csv")))
pool = np.array([r["pool_id"] for r in rows])
period = np.array([r["period"] for r in rows])
wc = np.array([r["width_class"] for r in rows])
w = np.array([num(r["width_spacings"]) for r in rows])
rate = np.array([num(r["rebalances"]) for r in rows])
f = np.array([num(r["f_stock"]) for r in rows])
ff = np.array([num(r["f_flow"]) for r in rows])
u = np.array([num(r["mean_abs_u"]) for r in rows])
y = np.log(f)

print(f"pool-days {len(rows)}, pools {len(set(pool))}")
print(f"f_stock: median {np.median(f):.2f}, "
      f"IQR [{np.percentile(f, 25):.2f}, {np.percentile(f, 75):.2f}]; "
      f"f_flow median {np.nanmedian(ff):.2f}")

X1 = np.column_stack([np.ones(len(y)), np.log(rate)])
b1 = fit(X1, y)
lo, hi = ci(bootstrap(X1, y, rng)[:, 1])
law = lambda r: math.exp(b1[0] + b1[1] * math.log(r))
print(f"\nRATE LAW: log f = {b1[0]:.3f} + {b1[1]:.3f} log rate ; "
      f"slope 95% CI [{lo:.3f}, {hi:.3f}] ; "
      f"f(1)={law(1):.2f} f(10)={law(10):.2f} f(100)={law(100):.2f}")

print("\nby rate bucket: n, median f, IQR, median f_flow")
for a, b, label in BUCKETS:
    m = (rate >= a) & (rate < b)
    print(f"  {label:>7}: {m.sum():4d}  {np.median(f[m]):6.2f}  "
          f"[{np.percentile(f[m], 25):5.2f}, {np.percentile(f[m], 75):5.2f}]"
          f"  {np.nanmedian(ff[m]):5.2f}")

print("\nby width class: n, median f, median rate")
for k in WIDTHS:
    m = wc == k
    print(f"  {k:>9}: {m.sum():4d}  {np.median(f[m]):6.2f}  "
          f"{np.median(rate[m]):6.1f}")

m = ~np.isnan(u)
X3 = np.column_stack([np.ones(m.sum()), np.log(rate[m]), np.log(w[m]), u[m]])
b3 = fit(X3, y[m])
bs3 = bootstrap(X3, y[m], rng)
print(f"\nTRIGGER DECOMPOSITION (n={m.sum()}): log f = {b3[0]:.2f} "
      f"+ {b3[1]:.3f} log rate + {b3[2]:.3f} log(width/spacing) "
      f"+ {b3[3]:.3f} mean|u|")
for j, name in [(1, "rate"), (2, "width"), (3, "|u|")]:
    lo, hi = ci(bs3[:, j])
    print(f"  {name}: 95% CI [{lo:.3f}, {hi:.3f}]")
print(f"  f shift for +0.2 in |u|: x{math.exp(0.2 * b3[3]):.3f}; "
      f"for width x2: x{2 ** b3[2]:.3f}")

# The 30 July trigger change, on the pools held in both windows.
held = set(pool[period == "before"]) & set(pool[period == "after"])
for label in ["before", "after"]:
    s = (period == label) & np.isin(pool, list(held))
    resid = np.median(y[s] - (b1[0] + b1[1] * np.log(rate[s])))
    print(f"\n{label} 07-30: pool-days {s.sum()}, "
          f"median rate {np.median(rate[s]):.0f}, "
          f"median f {np.median(f[s]):.2f}, "
          f"median residual from the rate law {resid:+.2f} "
          f"(x{math.exp(resid):.2f}), mean|u| {np.nanmean(u[s]):.2f}")

colours = {"tight": "#2C5F8D", "medium": "#2CA02C",
           "wide": "#B07D2B", "very_wide": "#C0392B"}
fig, ax = plt.subplots(figsize=(5.6, 3.9))
for k in WIDTHS:
    mm = wc == k
    ax.scatter(rate[mm], f[mm], s=9, alpha=0.45, c=colours[k],
               label=k.replace("_", " "), linewidths=0)
xs = np.logspace(-0.3, 3.1, 50)
fit_label = f"$f = {math.exp(b1[0]):.2f}\\,\\mathrm{{rate}}^{{{b1[1]:.2f}}}$"
ax.plot(xs, np.exp(b1[0] + b1[1] * np.log(xs)), "k-", lw=1.4, label=fit_label)
ax.axhline(1, color="#9AA6B2", lw=0.8, ls="--")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("rebalances per pool-day")
ax.set_ylabel("$f$ = realised LVR / trailing-24h stock LVR")
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig("../figures/f_vs_rate.pdf")
print("\nwrote ../figures/f_vs_rate.pdf")
