#!/usr/bin/env python3
"""
Multi-regime Monte Carlo: the controlled demonstration of the sigma^-2 regime.

The verification surface for "Loss-Versus-Rebalancing Under Emissions"
(../lvr-under-emissions.pdf, section 8). It fixes the position geometry and the
emission rate and sweeps only sigma, so the regime claim can be demonstrated
without the confounds a field measurement carries: cross-pool heterogeneity in
kappa, and a noisy single-snapshot volatility estimate.

For each sigma on a grid spanning 20%/yr to 300%/yr it Monte-Carlos GBM
log-price paths through one representative range-clipped V3 position and
accrues three quantities:

  LVR  realised loss-versus-rebalancing by the range-clipped quadratic-
       variation method,
         dLVR = d_var_logprice / (4 * bracket) * position_value,
       bracket = 2 - sqrt(Pa/P) - sqrt(P/Pb). Accrued from realised variation,
       so no sigma estimator enters; sigma is only the path parameter.
  F    fee income a fee-collecting LP would earn. The fee world has a no-
       arbitrage band of half-width eps = -log(1-f); arbitrage is the Skorokhod
       reflection that keeps the AMM price inside it. K, the realised
       reflection, is simulated from the same path (not assumed); the fee on
       the arbitrage trade is f * (L*sqrt(P)/2) * dK. This is why the fee
       "volume tied to sigma" is the band-gated reflection, throughput
       sigma^2/(2 eps), and not the raw price total variation, which scales as
       sigma: only the band-gated quantity gives F proportional to sigma^2 and
       hence F/LVR flat.
  E    emission income, a fixed rate per unit value per unit time, independent
       of sigma. Accrued while in range.

Two position-management modes:
  managed  the position is re-centred to a fresh unit position when the price
           leaves its range, representative of an actively managed CL position.
           The +/-10% range is a 2000-tick Aerodrome range. Time-in-range
           stays ~1.
  static   the position is never re-centred, so time-in-range falls with
           sigma. This settles the note's claim that em/LVR proportional to
           sigma^-2 survives a sigma-varying time-in-range, because the
           in-range fraction cancels in the ratio.

Predicted: LVR proportional to sigma^2, F proportional to sigma^2 so F/LVR
flat at f/eps, E independent of sigma so em/LVR proportional to sigma^-2. The
convexity term V_REB - V_HODL is measured against LVR alongside, since the
note distinguishes the two benchmarks in section 9.

Self-contained; no data inputs, no network, no chain access. Deterministic
under SEED.

Outputs (this directory):
  - regime_mc_results.json : per-sigma tables (both modes) and fitted slopes.
  - regime_mc.png          : the two-panel diagnostic figure.

plot-regime-single.py re-plots the left panel of that figure as
../figures/emission-lvr-regime.png, the figure the note carries.
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
SEED = 20260510                              # programme-wide deterministic seed

WINDOW_DAYS = 30.0
YEAR_DAYS = 365.0
T_YEARS = WINDOW_DAYS / YEAR_DAYS
N_STEPS = 4000                               # fine enough to resolve the band
M_PATHS = 3000                               # Monte Carlo paths per sigma cell
SIGMA_GRID = np.geomspace(0.20, 3.00, 13)    # annualised volatility, 20%-300%
HALF_WIDTH = 0.10                            # +/-10% log-price (2000-tick range)
FEE = 0.05                                   # fee tier; see note below
EPS_BAND = -np.log(1.0 - FEE)                # Skorokhod no-arbitrage half-band
EMISSION_TARGET = 10.0                       # set emission rate so em/LVR ~ this
EMISSION_TARGET_SIGMA = 1.00                 # ... at this annualised sigma

# Note on FEE. The Skorokhod reflection law (throughput sigma^2/(2 eps)) holds
# only with scale separation sigma*sqrt(dt) << eps. A realistic 5-30 bps tier
# gives a band far narrower than any feasible time step, so the band cannot be
# resolved at a tractable step count. FEE is therefore set to a value whose
# band the discretisation resolves across the whole sigma grid. The regime
# result F/LVR = f/eps is flat in sigma for every fee tier, so the slope this
# script demonstrates is fee-tier-independent; only the flat level moves
# (f/eps -> 1 as f -> 0).


def simulate(sigma, recenter, emission_rate, rng):
    """One sigma cell. recenter=True is the managed mode, False the static mode.

    Returns per-path-averaged accruals.
    """
    dt = T_YEARS / N_STEPS
    M = M_PATHS
    h = HALF_WIDTH
    bracket_centre = 2.0 - 2.0 * np.exp(-h / 2.0)

    incr = rng.normal(0.0, sigma * np.sqrt(dt), size=(N_STEPS, M))

    logp = np.zeros(M)
    centre = np.zeros(M)
    L = np.ones(M) / (np.exp(centre / 2.0) * bracket_centre)   # unit value at P0

    Y = np.zeros(M)                          # Skorokhod mispricing in [-eps,eps]

    lvr = np.zeros(M)
    fee = np.zeros(M)
    emis = np.zeros(M)
    convexity = np.zeros(M)                  # sum of (x(P) - x_anchor) * dP
    in_range_time = np.zeros(M)
    n_recentre = np.zeros(M)

    def token0(lp, lo, hi, Lv):
        """In-range V3 token0 holding x = L (1/sqrt(P) - 1/sqrt(Pb)), clipped."""
        sqrt_p = np.exp(np.clip(lp, lo, hi) / 2.0)
        return Lv * (1.0 / sqrt_p - 1.0 / np.exp(hi / 2.0))

    x_anchor = token0(logp, centre - h, centre + h, L)   # HODL composition

    for i in range(N_STEPS):
        d = incr[i]
        logp_next = logp + d
        lo, hi = centre - h, centre + h
        in_range = (logp >= lo) & (logp < hi)
        in_range_time += in_range * dt

        # --- realised range-clipped LVR --------------------------------------
        clip_prev = np.clip(logp, lo, hi)
        clip_cur = np.clip(logp_next, lo, hi)
        d_var = (clip_cur - clip_prev) ** 2
        sqrt_pa_p = np.exp((lo - clip_prev) * 0.5)
        sqrt_p_pb = np.exp((clip_prev - hi) * 0.5)
        bracket = np.clip(2.0 - sqrt_pa_p - sqrt_p_pb, 1e-9, None)
        pos_val = L * np.exp(clip_prev / 2.0) * bracket      # V = L sqrt(P) bracket
        lvr += d_var / (4.0 * bracket) * pos_val

        # --- fee side: Skorokhod band reflection -----------------------------
        Y = Y - d
        over = np.maximum(Y - EPS_BAND, 0.0)
        Y = np.minimum(Y, EPS_BAND)
        under = np.maximum(-EPS_BAND - Y, 0.0)
        Y = np.maximum(Y, -EPS_BAND)
        dK = over + under                                    # realised reflection
        sqrt_p_clip = np.exp(clip_prev / 2.0)
        fee += np.where(in_range, FEE * (L * sqrt_p_clip / 2.0) * dK, 0.0)

        # --- emission: fixed rate, sigma-independent -------------------------
        emis += np.where(in_range, emission_rate * pos_val * dt, 0.0)

        # --- convexity term V_REB - V_HODL -----------------------------------
        dP = np.exp(logp_next) - np.exp(logp)
        x_now = token0(logp, lo, hi, L)
        convexity += np.where(in_range, (x_now - x_anchor) * dP, 0.0)

        logp = logp_next

        # --- re-centre on exit (managed mode); fresh unit position -----------
        if recenter:
            exited = (logp < lo) | (logp >= hi)
            if exited.any():
                n_recentre += exited
                centre = np.where(exited, logp, centre)
                L = np.where(exited,
                             1.0 / (np.exp(logp / 2.0) * bracket_centre), L)
                new_anchor = token0(logp, centre - h, centre + h, L)
                x_anchor = np.where(exited, new_anchor, x_anchor)

    return {
        "sigma": float(sigma),
        "lvr": float(lvr.mean()),
        "fee": float(fee.mean()),
        "emission": float(emis.mean()),
        "f_over_lvr": float((fee / np.maximum(lvr, 1e-30)).mean()),
        "em_over_lvr": float((emis / np.maximum(lvr, 1e-30)).mean()),
        "convexity_over_lvr": float(
            (np.abs(convexity) / np.maximum(lvr, 1e-30)).mean()),
        "time_in_range_frac": float((in_range_time / T_YEARS).mean()),
        "mean_recentres": float(n_recentre.mean()),
    }


def loglog_fit(xs, ys):
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    m = (xs > 0) & (ys > 0)
    b, a = np.polyfit(np.log(xs[m]), np.log(ys[m]), 1)
    yhat = a + b * np.log(xs[m])
    ss = np.sum((np.log(ys[m]) - np.log(ys[m]).mean()) ** 2)
    r2 = 1 - np.sum((np.log(ys[m]) - yhat) ** 2) / ss if ss > 0 else float("nan")
    return {"slope": float(b), "intercept": float(a), "r2": float(r2)}


def sweep(recenter, emission_rate, seed):
    rng = np.random.default_rng(seed)
    return [simulate(s, recenter, emission_rate, rng) for s in SIGMA_GRID]


def main():
    # calibrate the emission rate: em/LVR is exactly linear in it, one probe fixes it.
    probe = simulate(EMISSION_TARGET_SIGMA, True, 1.0,
                     np.random.default_rng(SEED))
    emission_rate = EMISSION_TARGET / probe["em_over_lvr"]

    managed = sweep(True, emission_rate, SEED + 1)
    static = sweep(False, emission_rate, SEED + 2)
    sig = [r["sigma"] for r in managed]

    fits = {
        "managed": {
            "lvr_vs_sigma": loglog_fit(sig, [r["lvr"] for r in managed]),
            "fee_vs_sigma": loglog_fit(sig, [r["fee"] for r in managed]),
            "emission_vs_sigma": loglog_fit(sig, [r["emission"] for r in managed]),
            "f_over_lvr_vs_sigma": loglog_fit(sig, [r["f_over_lvr"] for r in managed]),
            "em_over_lvr_vs_sigma": loglog_fit(sig, [r["em_over_lvr"] for r in managed]),
        },
        "static": {
            "lvr_vs_sigma": loglog_fit(sig, [r["lvr"] for r in static]),
            "em_over_lvr_vs_sigma": loglog_fit(sig, [r["em_over_lvr"] for r in static]),
        },
    }
    result = {
        "params": {
            "window_days": WINDOW_DAYS, "n_steps": N_STEPS, "m_paths": M_PATHS,
            "half_width_log": HALF_WIDTH, "fee": FEE, "eps_band": EPS_BAND,
            "f_over_eps": FEE / EPS_BAND,
            "emission_rate_per_year": emission_rate, "seed": SEED,
            "sigma_grid": [float(s) for s in SIGMA_GRID],
            "scale_separation_max": float(
                SIGMA_GRID[-1] * np.sqrt(T_YEARS / N_STEPS) / EPS_BAND),
        },
        "predictions": {"lvr": 2.0, "fee": 2.0, "emission": 0.0,
                        "f_over_lvr": 0.0, "em_over_lvr": -2.0},
        "fits": fits,
        "managed": managed,
        "static": static,
    }
    with open(OUT / "regime_mc_results.json", "w") as fh:
        json.dump(result, fh, indent=1)

    # --- figure ---
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.2))
    xs = np.array(sig)

    fo = [r["f_over_lvr"] for r in managed]
    eo = [r["em_over_lvr"] for r in managed]
    fF = fits["managed"]["f_over_lvr_vs_sigma"]
    fE = fits["managed"]["em_over_lvr_vs_sigma"]
    axL.plot(sig, fo, "o-", color="C0",
             label=f"F/LVR, fee-side  (slope {fF['slope']:+.2f}, predict 0)")
    axL.plot(sig, eo, "s-", color="C3",
             label=f"em/LVR, emission-side  (slope {fE['slope']:+.2f}, "
                   f"predict -2)")
    eo_static = [r["em_over_lvr"] for r in static]
    axL.plot(sig, eo_static, "s--", color="C3", alpha=0.45,
             label=f"em/LVR, static position  "
                   f"(slope {fits['static']['em_over_lvr_vs_sigma']['slope']:+.2f})")
    anchor = np.median([e * s ** 2 for e, s in zip(eo, sig)])
    axL.plot(xs, anchor * xs ** -2.0, color="0.45", lw=1.0,
             label=r"reference slope $\sigma^{-2}$")
    axL.axhline(FEE / EPS_BAND, color="C0", ls=":", lw=1.0,
                label=r"reference $f/(-\ln(1-f))$")
    axL.set_xscale("log")
    axL.set_yscale("log")
    axL.set_xlabel(r"annualised volatility $\sigma$")
    axL.set_ylabel("break-even ratio (income / LVR)")
    axL.set_title("(simulation)  break-even against volatility")
    axL.legend(fontsize=7.3, loc="best")
    axL.grid(True, which="major", alpha=0.25)

    fL = fits["managed"]["lvr_vs_sigma"]
    fFe = fits["managed"]["fee_vs_sigma"]
    fEm = fits["managed"]["emission_vs_sigma"]
    axR.plot(sig, [r["lvr"] for r in managed], "^-", color="C2",
             label=f"LVR  (slope {fL['slope']:+.2f}, predict +2)")
    axR.plot(sig, [r["fee"] for r in managed], "o-", color="C0",
             label=f"fee income  (slope {fFe['slope']:+.2f}, predict +2)")
    axR.plot(sig, [r["emission"] for r in managed], "s-", color="C3",
             label=f"emission income  (slope {fEm['slope']:+.2f}, predict 0)")
    axR.set_xscale("log")
    axR.set_yscale("log")
    axR.set_xlabel(r"annualised volatility $\sigma$")
    axR.set_ylabel("income or loss over the window (USD, unit position)")
    axR.set_title("(simulation)  LVR, fee, emission against volatility")
    axR.legend(fontsize=7.8, loc="best")
    axR.grid(True, which="major", alpha=0.25)

    fig.tight_layout()
    fig.savefig(OUT / "regime_mc.png", dpi=150)

    # --- console summary ---
    print("=== multi-regime Monte Carlo ===")
    print(f"sigma grid {SIGMA_GRID[0]:.2f}-{SIGMA_GRID[-1]:.2f} "
          f"({len(SIGMA_GRID)} pts)   paths {M_PATHS}   steps {N_STEPS}")
    print(f"position half-width +/-{HALF_WIDTH:.0%}   fee f={FEE} "
          f"(eps={EPS_BAND:.5f}, f/eps={FEE/EPS_BAND:.4f})")
    print(f"emission rate {emission_rate:.3f}/yr   "
          f"scale-separation sigma*sqrt(dt)/eps max "
          f"{result['params']['scale_separation_max']:.3f}")
    print()
    print(f"{'sigma':>7s} | {'LVR':>9s} {'F/LVR':>8s} {'em/LVR':>9s} "
          f"{'cvx/LVR':>8s} {'recentr':>8s} | {'em/LVR':>9s} {'TIR':>6s} (static)")
    for rm, rs in zip(managed, static):
        print(f"{rm['sigma']:7.3f} | {rm['lvr']:9.5f} {rm['f_over_lvr']:8.3f} "
              f"{rm['em_over_lvr']:9.3f} {rm['convexity_over_lvr']:8.3f} "
              f"{rm['mean_recentres']:8.1f} | {rs['em_over_lvr']:9.3f} "
              f"{rs['time_in_range_frac']:6.3f}")
    print()
    print("log-log fits  quantity = c * sigma^b :")
    m = fits["managed"]
    print(f"  managed  LVR       b = {m['lvr_vs_sigma']['slope']:+.3f}  "
          f"R2={m['lvr_vs_sigma']['r2']:.4f}   [predict +2]")
    print(f"  managed  fee       b = {m['fee_vs_sigma']['slope']:+.3f}  "
          f"R2={m['fee_vs_sigma']['r2']:.4f}   [predict +2]")
    print(f"  managed  emission  b = {m['emission_vs_sigma']['slope']:+.3f}  "
          f"R2={m['emission_vs_sigma']['r2']:.4f}   [predict  0]")
    print(f"  managed  F/LVR     b = {m['f_over_lvr_vs_sigma']['slope']:+.3f}  "
          f"R2={m['f_over_lvr_vs_sigma']['r2']:.4f}   [predict  0]")
    print(f"  managed  em/LVR    b = {m['em_over_lvr_vs_sigma']['slope']:+.3f}  "
          f"R2={m['em_over_lvr_vs_sigma']['r2']:.4f}   [predict -2]")
    print(f"  static   em/LVR    b = "
          f"{fits['static']['em_over_lvr_vs_sigma']['slope']:+.3f}  "
          f"R2={fits['static']['em_over_lvr_vs_sigma']['r2']:.4f}   "
          f"[predict -2; time-in-range varies]")
    print()
    print("outputs written: regime_mc_results.json, regime_mc.png")


if __name__ == "__main__":
    sys.exit(main())
