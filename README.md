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

  - Preprocessed data that are used to produce results presented in the
    manuscript: .rda
  - Original data from the EPA’s website: .sql

### Data dictionary

  - Available at the following URL:
    <https://github.com/jinbora0720/BMC/blob/master/data_dic.md>

# Part 2: Code

## Abstract

Code for Bayesian Matrix Completion (BMC) for Hypothesis Testing is
available at <https://github.com/jinbora0720/BMC>. This repository has
multiple R script files to reproduce results presented in the
manuscript. In order to fit BMC, one can consult `bmc.R` and
`ex_runbmc.R`. The main function `bmc.R` is able to fit data
with/without (default) heteroscedastic noises, using the adaptive
(default)/fixed number of factors in modeling \(\gamma_{ij}\). The
function also provides an option for simpler structures in
\(\gamma_{ij}\), namely, BMC\(_0\), BMC\(_i\), and BMC\(_j\) as
presented in the manuscript. These simpler structures can be used by
specifying “bmc0”, “bmci”, “bmcj” for `simpler` argument in `bmc` as in
`bmc(..., simpler="bmc0")`.

## Description

### Code format(s)

  - Script files
      - R
      - .cpp files

### Supporting software requirements

#### Version of primary software used

R version 3.6.0, 3.6.1, 4.0.0, and 4.0.2

#### Libraries and dependencies used by the code

MASS 7.3-51.6, Rcpp 1.0.5, tidyverse 1.3.0, splines 4.0.2

# Part 3: Reproducibility workflow

## Scope

Code to reproduce figures about data structure or observations is
available at
<https://github.com/jinbora0720/BMC/blob/master/R/eda_plots.R>.

The provided workflow reproduces selected tables and figures in the
paper regarding the analysis, that is, Figure 2-5, Figure S5-S12, and
Table 1 and Table S1.

## Workflow

### Format(s)

  - Master code files

### Instructions

Reproducible code:

  - `neuro.R` runs BMC and produces Figure S10 and S11, regarding
    neurodevelopmental disorders.
      - `obese_active_aenm.rds` contains obesity-relevant assay
        endpoints likely to be activated by top 5 chemicals (Triclosan,
        BPA, 2,4,5-Trichlorophenol, DDT, and p,p’-DDE). This is saved in
        advance by running `obese.R`.
      - Run in R version 3.6.1 under Fedora 29 (MATE-Compiz)
  - `obese.R` runs BMC and produces obesity-related results: Figure 4-5,
    S8, S9, and S12.
      - `neuro_active_aenm.rds` contains assay endpoints related to
        neurodevelopmental disorders, which are likely to be activated
        by the top 5 chemicals. This is saved in advance by running
        `neuro.R`.
      - Run in R version 3.6.1 under Fedora 29 (MATE-Compiz)
  - `sim1.R` conducts 50 simulation studies when BMC is the true data
    generating process. It runs BMC, BMC\(_0\), BMC\(_i\), BMC\(_j\),
    tcpl, and ZIPLL. Due to computational time, each model was run
    separately with parallelisation using four cores and three nodes. It
    produces Figure 2-3, S5-S7, and Table 1.
      - ZIPLL and tcpl run in R version 4.0.2 under macOS Catalina
        10.15.6
      - BMC’s run in R version 4.0.0 under Red Hat Enterprise Linux
        Server 7.8 (Maipo)
  - `sim2.R` conducts 50 simulation studies when BMC is not the true
    data generating process. It produces Table S1.
      - tcpl run in R version 4.0.2 under macOS Catalina 10.15.6
      - BMC’s and ZIPLL run in R version 3.6.1 under Fedora 29
        (MATE-Compiz)

Source code:

  - `simdata.R` simulates data with ZIPLL assumptions.

  - `simulate_data.R` consists of functions to simulate data by BMC
    assumptions and create missingness in data.

  - `dosres_plot.R` consists of functions to draw dose-response plots
    and residual versus fitted values plots.

  - `metadata.R` produces informative variables to process data
    (e.g. K\_ij: a matrix of number of observations for each (i,j)
    cell, idx\_j: a list of chemical indices with at least one
    observation for each assay endpoint, etc.).

  - `bmc.R` is the main function to run bmc, which requires
    
      - `bmc_sampler.R`
      - `bmc_sampler2.cpp`
      - `samplerBits.cpp`.

### Expected run-time

Approximate time needed to reproduce the analyses on a standard desktop
machine is \> 8 hours.

### Additional information

Additional dependencies: Rcpp 1.0.5, tcpl 2.0.2, ROCR 1.0-11, ZIPLL 1.0
(available at <https://github.com/AnderWilson/ZIPLL>) (HDInterval 0.2.2,
gridExtra 2.3, ggpubr 0.4.0, reshape2 1.4.4)
