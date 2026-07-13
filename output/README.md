# Simulation Scenarios: S1–S16 Overview

This file maps every scenario label used (S1–S16) in the aggregation script 
to the data-generating process (DGP) and design settings behind it.

The 16 scenarios are a **4 x 4 crossing**: four base DGPs (S1–S4) evaluated
under four (N, K) design cells. The four base DGPs are chosen to switch the 
two sufficient conditions of the median–mean equivalence proposition on and 
off one at a time. The DGP surface follows the non-linear treatment-effect 
setting of Wager and Athey (2018).

---

## 1. Quick-reference map

| Scenario | Base DGP | N     | K   | Design role                    |
|----------|----------|-------|-----|--------------------------------|
| **S1**   | S1       | 1000  | 10  | Baseline cell                  |
| **S2**   | S2       | 1000  | 10  | Baseline cell                  |
| **S3**   | S3       | 1000  | 10  | Baseline cell                  |
| **S4**   | S4       | 1000  | 10  | Baseline cell                  |
| **S5**   | S1       | 2000  | 10  | Sample-size check (N=2000)     |
| **S6**   | S2       | 2000  | 10  | Sample-size check (N=2000)     |
| **S7**   | S3       | 2000  | 10  | Sample-size check (N=2000)     |
| **S8**   | S4       | 2000  | 10  | Sample-size check (N=2000)     |
| **S9**   | S1       | 1000  | 20  | High-dimension check (K=20)    |
| **S10**  | S2       | 1000  | 20  | High-dimension check (K=20)    |
| **S11**  | S3       | 1000  | 20  | High-dimension check (K=20)    |
| **S12**  | S4       | 1000  | 20  | High-dimension check (K=20)    |
| **S13**  | S1       | 1000  | 5   | Low-dimension check (K=5)      |
| **S14**  | S2       | 1000  | 5   | Low-dimension check (K=5)      |
| **S15**  | S3       | 1000  | 5   | Low-dimension check (K=5)      |
| **S16**  | S4       | 1000  | 5   | Low-dimension check (K=5)      |


---

## 2. The four base DGPs (S1–S4)

All four share the same target CATE surface `tau(x)` and the same covariate
generation. They differ only in the outcome noise, the treatment-effect
structure, and the individual-level component `U_i`. Because `E[U_i] = 0`
throughout, the population CATE `tau(x)` is identical across S1–S4. The
median treatment effect `tau_med(x) = tau(x) + med(U_i)` coincides with
`tau(x)` in S1, S2, and S3, and separates from it only in S4.

Conditions (i) and (ii) below refer to the two sufficient conditions of the
median–mean equivalence proposition: (i) leafwise location shift, and
(ii) within-leaf symmetry of the individual treatment effect.

### S1 — Gaussian baseline
- **Noise:** `epsilon_i ~ N(0, 1)`
- **Treatment effect:** smooth, non-sparse `tau(x)`
- **ITE component:** `U_i = 0`
- **Conditions:** (i) holds, (ii) holds

### S2 — Heavy-tailed noise
- **Noise:** `epsilon_i ~ t_3 / sqrt(3)` (Student-t with 3 df)
- **Treatment effect:** smooth, non-sparse `tau(x)` (same as S1)
- **ITE component:** `U_i = 0`
- **Conditions:** (i) holds, (ii) holds

### S3 — Sparse extreme-responder region
- **Noise:** `epsilon_i ~ N(0, 1)`
- **Treatment effect:** sparse `tau(x)` with an extreme-responder boundary
- **ITE component:** `U_i = 0`
- **Conditions:** (i) violated, (ii) violated (jointly)

### S4 — Skewed individual treatment effect
- **Noise:** `epsilon_i ~ N(0, 1)`
- **Treatment effect:** smooth, non-sparse `tau(x)` (same as S1)
- **ITE component:** skewed, mean-zero `U_i`, independent of `X_i`
- **Conditions:** (i) holds, (ii) violated

---


