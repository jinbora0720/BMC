Data Dictionary
================

## Data Structure

`datalist.rds` is a list containing `neuro_data`, `obese_data`, and
`hit_vec`.

  - `neuro_data` is the processed ToxCast/Tox21 data with 30 chemicals
    and 131 assay endpoints relevant to neurodevelopmental disorders.
  - `obese_data` is the processed ToxCast/Tox21 data with 30 chemicals
    and 271 obesity-related assay endpoints.
  - `hit_vec` includes hit-call values based on invitroDBv2 of EPA.

## Variables

For `neuro_data` and `obese_data`,

  - `aeid`: assay endpoint ID
  - `aenm`: assay endpoint name
  - `casn`: CAS registry number
  - `code`: C followed by `casn` without hyphens
  - `chnm`: chemical name
  - `logc`: log base 10 concentration
  - `resp`: response value.

`hit_vec` has

  - `code`: C followed by `casn` without hyphens
  - `aenm`: assay endpoint name
  - `hitc`: hit-call value where 1 if active, 0 if inactive, -1 if
    cannot be determined, and NA if not tested.

<!-- end list -->

``` r
names(datalist)
```

    ## [1] "neuro_data" "obese_data" "hit_vec"

``` r
neuro_data <- datalist$neuro_data
obese_data <- datalist$neuro_data
hit_vec <- datalist$hit_vec
```
