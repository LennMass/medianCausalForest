# medianCausalForest

Replication materials for *Median-based Splitting Rules for the Causal Tree*, which introduces a **Median
Squared Deviation (MSD)** splitting criterion for honest causal trees within the causal forest framework
([Athey and Imbens, 2016](https://www.pnas.org/doi/10.1073/pnas.1510489113); [Wager and Athey, 2018](https://www.pnas.org/doi/10.1073/pnas.1510489113])).
Split selection is anchored on the **Hodges-Lehmann**
location-shift estimator, giving robust estimation of conditional average
treatment effects (CATEs) under heavy-tailed and skewed outcome distributions.

## Overview

Honest causal trees [(Athey and Imbens, 2016)](https://www.pnas.org/doi/10.1073/pnas.1510489113) select splits by maximising
estimated treatment-effect heterogeneity under a mean-based objective. When
outcomes are heavy-tailed or skewed, that objective becomes unstable and CATE
estimates might degrade. This repository replaces the mean-based criterion with the
Median Squared Deviation (MSD), evaluated around Hodges-Lehmann leaf effects, and
provides the simulation study and empirical applications reported in the paper.
MSD is the primary, theoretically grounded criterion. Median Absolute Deviation (MAD) and Least Median Squares (LMS) are included
as secondary heuristic variants.

## Repository structure

```
medianCausalForest/
├── README.md
├── LICENSE
├── .gitignore
├── dependencies.R          install required R packages
├── core/                   core functions 
│   ├──
│   └── 
├── simulation/             simulation files
│   ├──
│   └── 
├── application/
│   ├── progresa/           Progresa application (de la O, 2013)
│   └── actg175/            ACTG 175 application (Hammer et al., 1996)
└── output/
    ├── figures/             generated figures
    └── tables/              generated tables
```

## Requirements

- R (>= 4.1).
- Packages: `causalTree`, `htetree`, `grf`, `dplyr`, `ggplot2`, `patchwork`,
  `kableExtra`. The ACTG 175 data are commonly accessed through (`speff2trial`)[https://cran.r-project.org/web/packages/speff2trial/index.html].
  The Progresa data is accessed from the [repository](https://github.com/ghoshadi/RRE)
  that corresponds to [Ghosh et al. (2026)](https://arxiv.org/pdf/2111.15524). 
- `causalTree` is not on CRAN. It is installed from GitHub by `dependencies.R`.



## Outputs



## License


