README
================

# Bayesian Matrix Completion for Hypothesis Testing

This form documents the artifacts associated with the article (i.e., the
data and code supporting the computational findings) and describes how
to reproduce the findings.

# Part 1: Data

## Abstract

The US Environmental Protection Agency (EPA) has developed programs
called Toxicity Forecaster (ToxCast) and Toxicology Testing in the 21st
Century (Tox21) that involve thousands of chemicals tested across
numerous biological endpoints. All data from ToxCast/Tox21 are available
online at
<https://epa.figshare.com/articles/ToxCast/Tox21_Database_invitroDB_for_Mac_Users/6062620>.
We engage a subset of chemicals and assay endpoints relevant to two
biomedical outcomes: neurodevelopmental disorders and obesity. Selection
criteria for the chemicals and assay endpoints are discussed in the
Supplemental Material of the manuscript. We provide processed datasets
at <https://github.com/jinbora0720/BMC>.

## Description

### File format(s)

-   Preprocessed data that are used to produce results presented in the
    manuscript: .RDS
-   Original data from the EPA’s website: .sql

### Data dictionary

-   Available at the following URL:
    <https://github.com/jinbora0720/BMC/blob/master/data_dic.md>

# Part 2: Code

## Abstract

Code for Bayesian Matrix Completion (BMC) for Hypothesis Testing is
available at <https://github.com/jinbora0720/BMC>. This repository has
multiple R script files to reproduce results presented in the
manuscript. In order to fit BMC, one can consult `bmc.R`. The main
function `bmc.R` is able to fit data with/without heteroscedastic
noises, using the adaptive/fixed number of factors in modeling
*γ*<sub>*i**j*</sub> and *t*<sub>*i**j*</sub>. The function also
provides an option for simpler structures in *γ*<sub>*i**j*</sub>,
namely, BMC<sub>0</sub>, BMC<sub>*i*</sub>, and BMC<sub>*j*</sub> as
presented in the manuscript. These simpler structures can be used by
specifying “bmc0”, “bmci”, “bmcj” for `gamma_simpler` argument in `bmc`
as in `bmc(..., gamma_simpler="bmc0")`. Under “bmc0”, “bmci”, and
“bmcj”, *t*<sub>*i**j*</sub> follows a simple Bernoulli distribution
with a global probability parameter.

## Description

### Code format(s)

-   Script files
    -   R
    -   .cpp files

### Supporting software requirements

#### Version of primary software used

R version 4.1.0, 4.0.3

#### Libraries and dependencies used by the code

<!-- sessionInfo() -->

MASS 7.3-54, Rcpp 1.0.6, tidyverse 1.3.1, splines 4.1.0

# Part 3: Reproducibility workflow

## Scope

Code to reproduce figures about data structure or observations is
available at
<https://github.com/jinbora0720/BMC/blob/master/R/eda_plots.R>.

The provided workflow reproduces selected tables and figures in the
paper regarding the analysis, that is, Figure 2-6, Figure S4-S11, and
Table 1-2 and Table S1-S3.

## Workflow

### Format(s)

-   Master code files

### Instructions

Reproducible code:

-   `neuro.R` runs BMC and produces Figure 5 and S8, regarding
    neurodevelopmental disorders.
    -   `obese_active_aenm.RDS` contains obesity-relevant assay
        endpoints likely to be activated by top 5 chemicals (Triclosan,
        BPA, 2,4,5-Trichlorophenol, DDT, and p,p’-DDE). This is saved in
        advance by running `obese.R`.
    -   Run in R version 4.1.0 under Ubuntu 20.04.2 LTS
-   `obese.R` runs BMC and produces obesity-related results: Figure 4,
    6, S7, S9, and S11.
    -   `neuro_active_aenm.RDS` contains assay endpoints related to
        neurodevelopmental disorders, which are likely to be activated
        by the top 5 chemicals. This is saved in advance by running
        `neuro.R`.
    -   Run in R version 4.1.0 under Ubuntu 20.04.2 LTS
-   `sim1.R` conducts 30 simulation studies when BMC is the true data
    generating process. It runs BMC (`sim1_BMC_cluster.R`),
    BMC<sub>*i*</sub> (`sim1_BMCi_cluter.R`), tcpl, and ZIPLL
    (`sim1_tcplZIPLL.R`). Due to computational time, each model was run
    separately with parallelisation using ten cores and one node.
    Results are scraped by running `sim1_collectresult_cluster.R`. The R
    file `sim1.R` produces Figure 2-3, S4-S6, and Table 1.
    -   ZIPLL and tcpl run in R version 4.1.0 under Ubuntu 20.04.2 LTS
    -   BMC’s run in R version 4.0.3 under CentOS Linux 8
-   `sim2.R` conducts 50 simulation studies when BMC is not the true
    data generating process. It produces Table S1.
    -   Run in R version 4.1.0 under Ubuntu 20.04.2 LTS
-   `sim3.R` makes Table 2 with results from `sim3_m5J5.R`,
    `sim3_m5J6.R`, `sim3_m5J20.R`, `sim3_m5J50.R`, and `sim3_m5J100.R`.
    Simulation 3 shows how BMC adjusts for multiplicity.
    -   Run in R version 4.1.0 under Ubuntu 20.04.2 LTS
-   `sim4_highcor.R` and `sim4_weakcor.R` conduct 30 simulation studies
    with varying missingness under highly correlated chemical structure
    vs. weakly correlated structure, respectively. They produce Figure
    S10, and Table S3-S4.

Source code:

-   `bmc.R` is the main function to run bmc, which requires

    -   `bmc_sampler.R`
    -   `bmc_sampler.cpp`
    -   `samplerBits.cpp`.

-   `simulate_data.R` consists of functions to simulate data by BMC
    assumptions and create missingness in data.

-   `metadata.R` produces informative variables to process data
    (e.g. K\_ij: a matrix of number of observations for each (i,j) cell,
    idx\_j: a list of chemical indices with at least one observation for
    each assay endpoint, etc.).

-   `dosres_plot.R` consists of functions to draw dose-response plots
    and residual versus fitted values plots.

-   `msf.cpp` is for varimax rotation of latent factors.

-   `list_mean.cpp` is to compute mean across list elements.

-   `simdata.R` simulates data with ZIPLL assumptions.

-   `ZIPLL.R` is needed to properly run ZIPLL code.

### Expected run-time

Approximate time needed to reproduce the analyses on a standard desktop
machine is &gt; 8 hours.

### Additional information

Additional dependencies: Rcpp 1.0.6, tcpl 2.0.2, ROCR 1.0-11, ZIPLL 1.0
(available at <https://github.com/AnderWilson/ZIPLL>)
