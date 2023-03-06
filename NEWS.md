# hiddenmeta 0.1.8

## Bug Fixes

-   Lots of them
-   Adjusted vignettes
-   Allowed RDS sample to sample new seeds

# hiddenmeta 0.1.9

## Enhancements

-   Added function that allows to read study parameters from the Google spreadsheet
-   Added drafts of several vignettes: On reading study parameters, and on simple meta-analysis simulation

## Bug Fixes

-   Fixed meta-estimation handler to allow specification of custom Stan code
-   Added option to provide more controls over Stan fit call

# hiddenmeta 0.2.0

## Enhancements

-   Significantly reworked RDS sampling procedure
-   Full support of import and diagnosis of study designs from Google Spreadsheets

## Bug Fixes

-   Numerous bug fixes


# hiddenmeta 0.3.0

## Enhancements

-   Study level population, sampling and inquiries functions are rewritten using `data.table` instead of `tidyverse` for significant speed gains
-   `sample_rds` allows to specify RDS+ sampling with many heterogeneous seeds
-   introducing `get_study_est_linktrace` implementing Link-Tracing population size estimation on RDS+ sample

## Bug Fixes

-   Numerous bug fixes

# hiddenmeta 0.3.2

## Enhancements

-   Add Multiple Systems population size estimator

## Bug Fixes

-   Fix imports to avoid conflicts between packages
-   Unify naming of sample prefix arguments to be `prefix`
-   Change RDS+ to LTS (Link-Tracing Sample)
