# Simulation output

Captured run of `regime-mc.py`, the verification surface for the note
*Loss-Versus-Rebalancing Under Emissions* (`../lvr-under-emissions.pdf`).
Deterministic under `SEED = 20260510`.

```
=== multi-regime Monte Carlo ===
sigma grid 0.20-3.00 (13 pts)   paths 3000   steps 4000
position half-width +/-10%   fee f=0.05 (eps=0.05129, f/eps=0.9748)
emission rate 25.659/yr   scale-separation sigma*sqrt(dt)/eps max 0.265

  sigma |       LVR    F/LVR    em/LVR  cvx/LVR  recentr |    em/LVR    TIR (static)
  0.200 |   0.00843    0.715   249.532    0.952      0.2 |   249.669  0.976
  0.251 |   0.01323    0.800   158.949    0.905      0.4 |   158.880  0.942
  0.314 |   0.02077    0.867   101.197    0.726      0.6 |   101.158  0.885
  0.394 |   0.03263    0.894    64.436    0.539      1.1 |    64.442  0.810
  0.493 |   0.05127    0.897    41.005    0.458      1.8 |    41.079  0.720
  0.618 |   0.08043    0.918    26.132    0.362      2.9 |    26.209  0.629
  0.775 |   0.12619    0.924    16.653    0.302      4.6 |    16.702  0.551
  0.971 |   0.19798    0.927    10.615    0.239      7.3 |    10.683  0.457
  1.216 |   0.31056    0.911     6.766    0.195     11.1 |     6.825  0.380
  1.524 |   0.48670    0.910     4.318    0.160     17.5 |     4.363  0.313
  1.910 |   0.76151    0.894     2.759    0.126     26.9 |     2.798  0.263
  2.394 |   1.18983    0.880     1.766    0.103     41.4 |     1.790  0.213
  3.000 |   1.85374    0.863     1.133    0.084     63.1 |     1.152  0.170

log-log fits  quantity = c * sigma^b :
  managed  LVR       b = +1.993  R2=1.0000   [predict +2]
  managed  fee       b = +2.037  R2=0.9987   [predict +2]
  managed  emission  b = -0.000  R2=0.7815   [predict  0]
  managed  F/LVR     b = +0.044  R2=0.2905   [predict  0]
  managed  em/LVR    b = -1.994  R2=1.0000   [predict -2]
  static   em/LVR    b = -1.987  R2=1.0000   [predict -2; time-in-range varies]

outputs written: regime_mc_results.json, regime_mc.png
```
