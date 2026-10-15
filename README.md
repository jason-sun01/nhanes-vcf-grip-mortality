# Vertebral fracture, low grip strength, and all-cause mortality in US adults

This repository contains R code and tabular outputs supporting the NHANES 2013-2014 linked mortality study.

## Data sources
Download the public-use NHANES 2013-2014 files DEMO_H, DXXVFA_H, MGX_H, BMX_H, SMQ_H, ALQ_H, PAQ_H, DIQ_H, BPQ_H, and MCQ_H, together with the 2019 Public-use Linked Mortality File. Place all files in data_raw. Original NHANES files are not redistributed here.

## Software
R 4.6.1. Required packages: survival, survey, splines, haven, foreign, dplyr, readr, ggplot2, nnet, ragg, and patchwork. The optional svglite or Cairo package can be used for editable SVG export.

## Reproduction
Run scripts from the repository root in numeric order. Generated analytic data are written to data_processed, tabular outputs to results, and publication figures to figures. Scripts 05-07 reproduce the 20-imputation analysis, proportional-hazards diagnostics and exploratory interaction/influence analyses, and income-to-poverty-ratio non-proportionality sensitivity analyses reported in the article. Script 08 fits a four-knot restricted cubic spline for continuous grip strength and limited exploratory sex and age interaction models, then exports the dose-response and subgroup figures.

## Definitions
Grade >=2 VCF means at least one T4-L4 vertebra with fracture grade 2 or 3. Low grip strength means maximum single-trial grip strength below 26 kg in men or below 16 kg in women.

## Notes
The models retain NHANES MEC weights, strata, and PSUs. Both exposures are entered together. Sequential models use the fixed Model 4 complete-case cohort. Primary estimates are independently validated with survey::svycoxph.

## Data availability
The underlying data are publicly available from NCHS. No restricted-use or identifiable data are included.

## Licence
The analysis code is released under the MIT License.

## Citation
The current Zenodo record is the previously released version 1.0.1: Sun S, Zhu L, Xiong H, Cui H, Yang X (2026). Code for vertebral fracture, low grip strength, and all-cause mortality in US adults. Zenodo. https://doi.org/10.5281/zenodo.22930336. This local enhancement bundle adds script 08 and its outputs; create a new GitHub/Zenodo release before resubmission and then update the version-specific citation if Zenodo assigns a new DOI. The associated article DOI will be added when available.

