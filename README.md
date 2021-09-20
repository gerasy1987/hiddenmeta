
# hiddenmeta

![R](https://github.com/gerasy1987/hiddenmeta/workflows/R/badge.svg)
[![codecov](https://codecov.io/gh/gerasy1987/hiddenmeta/branch/main/graph/badge.svg?token=ZG9A64Q0A1)](https://codecov.io/gh/gerasy1987/hiddenmeta)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/gerasy1987/hiddenmeta/blob/main/LICENSE)

### Installation

``` r
install.packages("DeclareDesign")

devtools::install_github("gerasy1987/hiddenmeta", build_vignettes = TRUE)
```

### Overview

We rely on [`DeclareDesign`](https://declaredesign.org/) package for the
simulations. The key idea of the simulation workflow is that we simulate
the design as a set of steps. When concatenated the steps produce an
instance of a design. When replicated these instances let us assess the
overall performance of a design.

Fully specified designs include information on:

1.  **M**: The background model—for instance, background beliefs about
    network structures.
2.  **I**: The inquiries, i.e. what exactly is the quantity we want to
    estimate. The **estimands** should be readable from the background
    model.
3.  **D**: What is the **data** strategy? For instance how will sampling
    be undertaken?
4.  **A**: What is the **answer** strategy? How will estimation be done?

Each of these MIDA steps can be defined using “handlers” from the
`hiddenmeta` package and can be adjusted to reflect details of
individual study designs.

### Individual designs

The individual study design simulation involves:

1.  Selecting parameters
2.  Using these to declare a design
3.  Diagnosing the design

### Meta analysis

Simulated study level results feed into the meta-analysis:

1.  Conduct multi-study design for as many sampling-estimator pairs in
    each study as possible, then diagnose the multi-study design.
    Calculate average (across simulations) estimand and bias of
    sampling-estimator for each of the studies and estimator sampling
    strategies. These will serve as population we will be drawing
    population-sampling-estimator triads

2.  Estimands include:

    -   Average estimand by inquiry label (within study)
    -   Average bias of specific sampling-estimator pair (across
        studies) compared to truth
    -   Average relative bias of sampling-estimator pair (across
        studies) compared to “gold standard”
    -   Ratio of average bias to costs of sampling-estimator pair

3.  Sampling consists of drawing population-sampling-estimator triads
    presuming that each study uses at least two sampling strategies at a
    time

4.  Once we draw sample we use Stan model to estimate study-specific
    estimands and sampling-estimator specific errors (biases, study
    level prevalence or size, cost-effectiveness)

------------------------------------------------------------------------

### Getting started

-   To familiarize yourself with the `hiddenmeta` workflow please read
    [this
    vignette](https://gsyunyaev.com/hiddenmeta/articles/hiddenmeta-base.html)

------------------------------------------------------------------------

This project was generously supported by a grant from the [African
Programming & Research Initiative to End Slavery
(APRIES)](https://apries.uga.edu/prevalenceforum/).
