# Vertebral fracture, low grip strength, and all-cause mortality in US adults

This repository contains R code and tabular outputs supporting the NHANES 2013-2014 linked mortality study.

## Data sources
Download the public-use NHANES 2013-2014 files DEMO_H, DXXVFA_H, MGX_H, BMX_H, SMQ_H, ALQ_H, PAQ_H, DIQ_H, BPQ_H, and MCQ_H, together with the 2019 Public-use Linked Mortality File. Place all files in data_raw. Original NHANES files are not redistributed here.

## Software
R 4.6.1. Required packages: survival, survey, splines, haven, foreign, dplyr, readr, ggplot2.

## Reproduction
Run scripts from the repository root in numeric order. Generated analytic data are written to data_processed and outputs to results.

## Definitions
Grade >=2 VCF means at least one T4-L4 vertebra with fracture grade 2 or 3. Low grip strength means maximum single-trial grip strength below 26 kg in men or below 16 kg in women.

## Notes
The models retain NHANES MEC weights, strata, and PSUs. Both exposures are entered together. Sequential models use the fixed Model 4 complete-case cohort. Primary estimates are independently validated with survey::svycoxph.

## Data availability
The underlying data are publicly available from NCHS. No restricted-use or identifiable data are included.

## Licence
The analysis code is released under the MIT License.

## Citation
Repository: https://github.com/jason-sun01/nhanes-vcf-grip-mortality. The article DOI and archived release DOI will be added when available.
