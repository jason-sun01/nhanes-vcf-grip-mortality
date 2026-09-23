# NHANES 2013-2014 VCF + low grip + 2019 mortality feasibility analysis
# Packages: haven, dplyr, readr
library(haven)
library(dplyr)
library(readr)

project_root <- normalizePath(".", winslash="/", mustWork=TRUE)
data_dir <- file.path(project_root, "data_raw")
dir.create(file.path(project_root, "data_processed"), showWarnings=FALSE)
dir.create(file.path(project_root, "results"), showWarnings=FALSE)
vfa  <- read_xpt(file.path(data_dir, "DXXVFA_H.xpt"))
grip <- read_xpt(file.path(data_dir, "MGX_H.xpt"))
demo <- read_xpt(file.path(data_dir, "DEMO_H.xpt"))
mort <- read_fwf(
  file.path(data_dir, "NHANES_2013_2014_MORT_2019_PUBLIC.dat"),
  fwf_cols(SEQN=c(1,6), ELIGSTAT=c(15,15), MORTSTAT=c(16,16),
           UCOD_LEADING=c(17,19), DIABETES=c(20,20), HYPERTEN=c(21,21),
           PERMTH_INT=c(43,45), PERMTH_EXM=c(46,48)),
  col_types="iiiiiiii", na=c("", "."))

fx <- c(paste0("DXXT",4:12,"FX"), paste0("DXXL",1:4,"FX"))
trials <- c("MGXH1T1","MGXH2T1","MGXH1T2","MGXH2T2","MGXH1T3","MGXH2T3")

vfa2 <- vfa %>%
  mutate(across(all_of(fx), ~ifelse(abs(.x)<1e-12, 0, .x)),
         vcf_grade2=as.integer(if_any(all_of(fx), ~.x>=2)),
         n_vertebrae_read=rowSums(!is.na(pick(all_of(fx))))) %>%
  filter(DXXVFAST %in% c(1,2)) %>%
  select(SEQN,DXXVFAST,vcf_grade2,n_vertebrae_read)

grip2 <- grip %>%
  mutate(grip_max=do.call(pmax, c(pick(all_of(trials)), na.rm=TRUE)),
         grip_max=ifelse(is.infinite(grip_max), NA, grip_max)) %>%
  select(SEQN,grip_max)

d <- demo %>% filter(RIDAGEYR>=40) %>%
  select(SEQN,RIDAGEYR,RIAGENDR,WTMEC2YR,SDMVPSU,SDMVSTRA) %>%
  inner_join(vfa2,by="SEQN") %>% inner_join(grip2,by="SEQN") %>%
  inner_join(mort,by="SEQN") %>% filter(ELIGSTAT==1,!is.na(grip_max)) %>%
  mutate(low_grip=as.integer((RIAGENDR==1 & grip_max<26) |
                             (RIAGENDR==2 & grip_max<16)),
         phenotype=case_when(
           vcf_grade2==0 & low_grip==0 ~ "No VCF + normal grip",
           vcf_grade2==0 & low_grip==1 ~ "No VCF + low grip",
           vcf_grade2==1 & low_grip==0 ~ "Grade>=2 VCF + normal grip",
           TRUE ~ "Grade>=2 VCF + low grip"),
         followup_years=PERMTH_EXM/12)

result <- d %>% group_by(phenotype) %>% summarise(
  N=n(), deaths=sum(MORTSTAT), death_percent=100*mean(MORTSTAT),
  mean_followup_years=mean(followup_years),
  median_followup_years=median(followup_years),
  total_person_years=sum(followup_years),
  deaths_per_1000_person_years=1000*deaths/total_person_years,
  .groups="drop")
print(result)

# Demographic-adjusted screening Cox models (unweighted; not the final NHANES model)
library(survival)
d <- d %>% mutate(race=factor(RIDRETH1), age10=RIDAGEYR/10)
fit_vcf <- coxph(Surv(followup_years,MORTSTAT) ~ vcf_grade2+age10+factor(RIAGENDR)+race, data=d)
fit_grip <- coxph(Surv(followup_years,MORTSTAT) ~ low_grip+age10+factor(RIAGENDR)+race, data=d)
fit_mutual <- coxph(Surv(followup_years,MORTSTAT) ~ vcf_grade2+low_grip+age10+factor(RIAGENDR)+race, data=d)
summary(fit_vcf); summary(fit_grip); summary(fit_mutual)

write.csv(d, file.path(project_root, "data_processed", "nhanes_vcf_grip_analysis_dataset.csv"), row.names=FALSE)
write.csv(result, file.path(project_root, "results", "feasibility_by_phenotype.csv"), row.names=FALSE)
