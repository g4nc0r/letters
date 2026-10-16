# Simulation output

Captured run of `fit.py`, the verification surface for the note
*Loss-Versus-Rebalancing Under Management* (`../lvr-under-management.pdf`).
Deterministic under `SEED = 20261016`.

```
pool-days 1453, pools 30
f_stock: median 3.36, IQR [1.71, 7.83]; f_flow median 2.27

RATE LAW: log f = 0.263 + 0.356 log rate ; slope 95% CI [0.309, 0.403] ; f(1)=1.30 f(10)=2.95 f(100)=6.70

by rate bucket: n, median f, IQR, median f_flow
       <5:  287    2.02  [ 0.86,  4.32]   3.15
     5-20:  394    2.60  [ 1.40,  6.04]   2.05
    20-50:  455    3.99  [ 2.11,  7.73]   2.21
   50-150:  213    4.78  [ 2.51,  8.89]   2.45
    >=150:  104   14.12  [ 8.18, 21.96]   2.19

by width class: n, median f, median rate
      tight:  243   16.63    48.0
     medium:  714    2.94    29.0
       wide:  418    2.61     7.0
  very_wide:   78    2.86     3.0

TRIGGER DECOMPOSITION (n=1453): log f = -0.40 + 0.356 log rate + 0.118 log(width/spacing) + 0.638 mean|u|
  rate: 95% CI [0.312, 0.403]
  width: 95% CI [0.057, 0.183]
  |u|: 95% CI [0.468, 0.926]
  f shift for +0.2 in |u|: x1.136; for width x2: x1.086

before 07-30: pool-days 70, median rate 40, median f 4.00, median residual from the rate law -0.12 (x0.89), mean|u| 0.65

after 07-30: pool-days 63, median rate 13, median f 2.22, median residual from the rate law -0.43 (x0.65), mean|u| 0.67

wrote ../figures/f_vs_rate.pdf
```
