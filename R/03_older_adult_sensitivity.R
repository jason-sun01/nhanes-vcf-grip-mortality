project_root <- normalizePath(".", winslash="/", mustWork=TRUE)
out_dir <- file.path(project_root, "results")
base <- file.path(project_root, "data_raw")
dir.create(out_dir, showWarnings=FALSE)
library(survival)

find_xpt <- function(name) {
  z <- list.files(base, pattern = paste0("^", name, "$"), recursive = TRUE,
    full.names = TRUE, ignore.case = TRUE)
  if (length(z) != 1) stop("Expected exactly one file: ", name)
  z
}
read_keep <- function(name, vars) {
  x <- foreign::read.xport(find_xpt(name))
  x[, intersect(vars, names(x)), drop = FALSE]
}

d <- read.csv(file.path(project_root, "data_processed", "nhanes_vcf_grip_analysis_dataset.csv"))
demo <- foreign::read.xport(find_xpt("DEMO_H.xpt"))[, c("SEQN", "DMDEDUC2", "INDFMPIR")]
bmx <- read_keep("BMX_H.xpt", c("SEQN", "BMXBMI"))
smq <- read_keep("SMQ_H.xpt", c("SEQN", "SMQ020", "SMQ040"))
alq <- read_keep("ALQ_H.xpt", c("SEQN", "ALQ101", "ALQ110"))
paq <- read_keep("PAQ_H.xpt", c("SEQN", "PAQ650", "PAQ665"))
diq <- read_keep("DIQ_H.xpt", c("SEQN", "DIQ010"))
bpq <- read_keep("BPQ_H.xpt", c("SEQN", "BPQ020"))
mcq <- read_keep("MCQ_H.xpt", c("SEQN", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F", "MCQ220"))
for (x in list(demo, bmx, smq, alq, paq, diq, bpq, mcq)) d <- merge(d, x, by = "SEQN", all.x = TRUE)

d$time <- d$PERMTH_EXM / 12
d$event <- as.integer(d$MORTSTAT)
d$age10 <- d$RIDAGEYR / 10
d$female <- as.integer(d$RIAGENDR == 2)
d$race <- factor(d$RIDRETH1)
d$education <- factor(ifelse(d$DMDEDUC2 %in% 1:5, d$DMDEDUC2, NA))
d$bmi <- ifelse(d$BMXBMI > 10 & d$BMXBMI < 100, d$BMXBMI, NA)
d$smoking <- NA_character_
d$smoking[d$SMQ020 == 2] <- "Never"
d$smoking[d$SMQ020 == 1 & d$SMQ040 == 3] <- "Former"
d$smoking[d$SMQ020 == 1 & d$SMQ040 %in% c(1, 2)] <- "Current"
d$smoking <- relevel(factor(d$smoking), "Never")
d$alcohol <- NA_character_
d$alcohol[d$ALQ101 == 1] <- "Past-year drinker"
d$alcohol[d$ALQ101 == 2 & d$ALQ110 == 1] <- "Former/light lifetime"
d$alcohol[d$ALQ101 == 2 & d$ALQ110 == 2] <- "Lifetime <12 drinks"
d$alcohol <- relevel(factor(d$alcohol), "Lifetime <12 drinks")
d$phys_active <- ifelse(d$PAQ650 == 1 | d$PAQ665 == 1, 1,
  ifelse(d$PAQ650 == 2 & d$PAQ665 == 2, 0, NA))
d$diabetes <- ifelse(d$DIQ010 == 1, 1, ifelse(d$DIQ010 %in% c(2, 3), 0, NA))
d$hypertension <- ifelse(d$BPQ020 == 1, 1, ifelse(d$BPQ020 == 2, 0, NA))
heartvars <- c("MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F")
d$cvd <- apply(d[, heartvars], 1, function(x) if (any(x == 1, na.rm = TRUE)) 1 else if (all(x == 2)) 0 else NA)
d$cancer <- ifelse(d$MCQ220 == 1, 1, ifelse(d$MCQ220 == 2, 0, NA))

design_cox <- function(formula, data, label) {
  vars <- all.vars(formula)
  keep <- complete.cases(data[, unique(c(vars, "WTMEC2YR", "SDMVSTRA", "SDMVPSU"))])
  z <- data[keep, ]
  z$w_norm <- z$WTMEC2YR / mean(z$WTMEC2YR)
  fit <- coxph(formula, data = z, weights = w_norm, ties = "breslow", x = TRUE)
  db <- residuals(fit, type = "dfbeta", weighted = TRUE)
  agg <- rowsum(db, interaction(z$SDMVSTRA, z$SDMVPSU, drop = TRUE), reorder = FALSE)
  stratum <- as.numeric(sub("\\..*$", "", rownames(agg)))
  V <- matrix(0, ncol(db), ncol(db))
  for (h in unique(stratum)) {
    u <- agg[stratum == h, , drop = FALSE]
    m <- nrow(u)
    if (m > 1) {
      uc <- sweep(u, 2, colMeans(u))
      V <- V + m / (m - 1) * crossprod(uc)
    }
  }
  b <- coef(fit)
  se <- sqrt(diag(V))
  list(
    result = data.frame(subgroup = label, N = nrow(z), deaths = sum(z$event),
      term = names(b), HR = exp(b), CI_low = exp(b - 1.96 * se),
      CI_high = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)), row.names = NULL),
    data = z
  )
}

f4 <- Surv(time, event) ~ vcf_grade2 + low_grip + age10 + female + race +
  education + INDFMPIR + bmi + smoking + alcohol + phys_active +
  diabetes + hypertension + cvd + cancer

thresholds <- c(40, 50, 60, 65)
fits <- lapply(thresholds, function(a) design_cox(f4, d[d$RIDAGEYR >= a, ], paste0("Age >=", a, " years")))
all_terms <- do.call(rbind, lapply(fits, `[[`, "result"))
main <- all_terms[all_terms$term %in% c("vcf_grade2", "low_grip"), ]
counts <- do.call(rbind, Map(function(x, a) {
  z <- x$data
  data.frame(subgroup = paste0("Age >=", a, " years"), N = nrow(z), deaths = sum(z$event),
    vcf_N = sum(z$vcf_grade2 == 1), vcf_deaths = sum(z$event[z$vcf_grade2 == 1]),
    low_grip_N = sum(z$low_grip == 1), low_grip_deaths = sum(z$event[z$low_grip == 1]),
    joint_N = sum(z$vcf_grade2 == 1 & z$low_grip == 1),
    joint_deaths = sum(z$event[z$vcf_grade2 == 1 & z$low_grip == 1]))
}, fits, thresholds))

write.csv(all_terms, file.path(out_dir, "older_adult_sensitivity_all_terms.csv"), row.names = FALSE)
write.csv(main, file.path(out_dir, "older_adult_sensitivity_main_exposures.csv"), row.names = FALSE)
write.csv(counts, file.path(out_dir, "older_adult_sensitivity_counts.csv"), row.names = FALSE)
print(main)
print(counts)
