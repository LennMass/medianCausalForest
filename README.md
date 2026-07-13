# medianCausalForest

Replication materials for *Median-based Splitting Rules for the Causal Tree*, which introduces a **Median
Squared Deviation (MSD)** splitting criterion for honest causal trees within the causal forest framework
([Athey and Imbens, 2016](https://www.pnas.org/doi/10.1073/pnas.1510489113); [Wager and Athey, 2018](https://www.pnas.org/doi/10.1073/pnas.1510489113])).
Split selection is anchored on the **Hodges-Lehmann**
location-shift estimator, giving robust estimation of conditional average
treatment effects (CATEs) under heavy-tailed and skewed outcome distributions.

## Overview

Honest causal trees [(Athey and Imbens, 2016)](https://www.pnas.org/doi/10.1073/pnas.1510489113) select splits by maximizing
estimated treatment effect heterogeneity under a mean-based objective. When
outcomes are heavy-tailed or skewed, that objective becomes unstable and CATE
estimates might degrade. This repository replaces the mean-based criterion with the
Median Squared Deviation (MSD), evaluated around Hodges-Lehmann leaf effects, and
provides the simulation study and empirical applications reported in the paper.
MSD is the primary criterion that we implement, besides including the secondary criterion variants of 
Median Absolute Deviation (MAD) and Least Median Squares (LMS). 

## Repository structure

```
medianCausalForest/
├── README.md
├── LICENSE
├── .gitignore
├── 00_setup.R                 install required R packages once
├── 01_config.R                load required R packages
├── core/                      core functions 
│   └── robustCausalTree/      local causal forest package that integrates median-based splitting rules
├── simulation/                simulation files
│   ├── raw_results/
│   └── aggrgeated_results/
├── application/               application files
│   ├── progresa/              
│   └── actg175/               
└── output/                    generated results
│   ├── application_plots/                     
│   ├── application_tables/                     
│   ├── sim_plots/                     
│   └── sim_tables/                     
```

## Requirements

- R (>= 4.3.0). Tested under R version 4.5.3.
- Necessary packages are installed with `00_setup.R` and loaded with `01_config.R`.
- The ACTG 175 data are commonly accessed through [`speff2trial`](https://cran.r-project.org/web/packages/speff2trial/index.html).
  The Progresa data are accessed from the public [repository](https://github.com/ghoshadi/RRE)
  that corresponds to [Ghosh et al. (2026)](https://arxiv.org/pdf/2111.15524).





