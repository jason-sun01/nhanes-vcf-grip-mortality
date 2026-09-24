library(survival)

d <- read.csv("data_processed/nhanes_full_covariate_analysis_dataset.csv", na.strings=c("NA", ""))
d$race <- factor(d$RIDRETH1)
d$education <- factor(d$education, levels=1:5)
d$smoking <- factor(d$smoking, levels=c("Never", "Former", "Current"))
d$alcohol <- factor(d$alcohol,
                    levels=c("Lifetime <12 drinks", "Former/light lifetime", "Past-year drinker"))

# Keep the exact Model 4 complete-case cohort.
base_vars <- c("time", "event", "vcf_grade2", "low_grip", "age10", "female", "race",
               "education", "INDFMPIR", "bmi", "smoking", "alcohol", "phys_active",
               "diabetes", "hypertension", "cvd", "cancer", "WTMEC2YR", "SDMVSTRA", "SDMVPSU")
cc <- d[complete.cases(d[, base_vars]), ]

design_fit <- function(formula, data, label) {
  z <- data
  z$w_norm <- z$WTMEC2YR / mean(z$WTMEC2YR)
  fit <- coxph(formula, data=z, weights=w_norm, ties="breslow", x=TRUE)
  db <- residuals(fit, type="dfbeta", weighted=TRUE)
  if (is.null(dim(db))) db <- matrix(db, ncol=1)
  colnames(db) <- names(coef(fit))
  agg <- rowsum(db, interaction(z$SDMVSTRA, z$SDMVPSU, drop=TRUE), reorder=FALSE)
  stratum <- as.numeric(sub("\\..*$", "", rownames(agg)))
  V <- matrix(0, ncol(db), ncol(db))
  for (h in unique(stratum)) {
    u <- agg[stratum == h, , drop=FALSE]
    m <- nrow(u)
    if (m > 1) {
      uc <- sweep(u, 2, colMeans(u))
      V <- V + m/(m-1) * crossprod(uc)
    }
  }
  b <- coef(fit)
  se <- sqrt(diag(V))
  data.frame(analysis=label, N=length(unique(z$SEQN)), deaths=sum(z$event), term=names(b),
             HR=exp(b), CI_low=exp(b-1.96*se), CI_high=exp(b+1.96*se),
             p=2*pnorm(-abs(b/se)), row.names=NULL)
}

common_cov <- "+age10+female+race+education+bmi+smoking+alcohol+phys_active+diabetes+hypertension+cvd+cancer"

# Primary specification retained for a direct same-sample comparison.
f_primary <- as.formula(paste("Surv(time,event)~vcf_grade2+low_grip+INDFMPIR", common_cov))
primary <- design_fit(f_primary, cc, "Primary: continuous PIR")

# Sensitivity 1: allow a different baseline hazard across weighted-sample PIR quartiles.
q <- quantile(cc$INDFMPIR, probs=c(0, .25, .5, .75, 1), na.rm=TRUE)
q <- unique(q)
cc$pir_quartile <- cut(cc$INDFMPIR, breaks=q, include.lowest=TRUE, ordered_result=TRUE)
f_strata <- as.formula(paste("Surv(time,event)~vcf_grade2+low_grip+strata(pir_quartile)", common_cov))
pir_strata <- design_fit(f_strata, cc, "PIR-quartile-stratified baseline hazard")

# Sensitivity 2: permit PIR effects to differ across prespecified follow-up
# intervals (0-2, 2-4, and >4 years), avoiding the risk-set expansion of tt().
split <- survSplit(Surv(time, event) ~ ., data=cc, cut=c(2, 4),
                   start="tstart", end="time", event="event", episode="period")
split$period <- factor(split$period)
f_piece <- as.formula(paste(
  "Surv(tstart,time,event)~vcf_grade2+low_grip+INDFMPIR+INDFMPIR:period+strata(period)",
  common_cov))
pir_piece <- design_fit(f_piece, split, "Piecewise PIR effects (0-2, 2-4, >4 y)")

all_results <- rbind(primary, pir_strata, pir_piece)
main_terms <- all_results[all_results$term %in% c("vcf_grade2", "low_grip", "INDFMPIR",
                                                   "INDFMPIR:period2", "INDFMPIR:period3"), ]
write.csv(all_results, "results/nhanes_PIR_nonproportional_sensitivity_all_terms.csv", row.names=FALSE)
write.csv(main_terms, "results/nhanes_PIR_nonproportional_sensitivity_key_results.csv", row.names=FALSE)

cat("Model 4 complete-case cohort:", nrow(cc), "participants and", sum(cc$event), "deaths\n")
print(main_terms)

