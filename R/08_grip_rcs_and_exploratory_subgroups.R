# Continuous grip-strength dose-response and limited exploratory subgroup analyses
# NHANES 2013-2014 VFA + grip strength + 2019 linked mortality
#
# This script uses the same Model 4 complete-case cohort and the same
# Taylor-linearized, design-based variance estimator that was independently
# validated against survey::svycoxph in R/02_primary_and_validation_models.R.

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(dplyr)
  library(ragg)
})

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "results")
fig_dir <- file.path(project_root, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

d <- read.csv(
  file.path(project_root, "data_processed", "nhanes_full_covariate_analysis_dataset.csv"),
  na.strings = c("NA", "")
)

# Restore analysis factor levels used in the primary model.
d$race <- factor(d$RIDRETH1)
d$education <- factor(d$education, levels = 1:5)
d$smoking <- factor(d$smoking, levels = c("Never", "Former", "Current"))
d$alcohol <- factor(
  d$alcohol,
  levels = c("Lifetime <12 drinks", "Former/light lifetime", "Past-year drinker")
)
d$older65 <- as.integer(d$RIDAGEYR >= 65)

primary_vars <- c(
  "time", "event", "vcf_grade2", "low_grip", "grip_max", "age10", "female",
  "race", "education", "INDFMPIR", "bmi", "smoking", "alcohol", "phys_active",
  "diabetes", "hypertension", "cvd", "cancer", "WTMEC2YR", "SDMVSTRA", "SDMVPSU"
)
cc <- d[complete.cases(d[, primary_vars]), ]
stopifnot(nrow(cc) == 2623L, sum(cc$event) == 220L)

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

# Harrell-style restricted cubic spline basis. With four knots, the basis has
# one linear and two nonlinear terms. Scaling affects coefficients but not fit.
rcs_basis <- function(x, knots) {
  stopifnot(length(knots) >= 4L, all(diff(knots) > 0))
  k <- length(knots)
  tp3 <- function(z, a) pmax(z - a, 0)^3
  den <- (knots[k] - knots[1])^2
  nonlinear <- do.call(cbind, lapply(seq_len(k - 2L), function(j) {
    tp3(x, knots[j]) -
      tp3(x, knots[k - 1L]) * (knots[k] - knots[j]) / (knots[k] - knots[k - 1L]) +
      tp3(x, knots[k]) * (knots[k - 1L] - knots[j]) / (knots[k] - knots[k - 1L])
  })) / den
  out <- cbind(x, nonlinear)
  colnames(out) <- c("grip_rcs_linear", paste0("grip_rcs_nl", seq_len(k - 2L)))
  out
}

# Fit a weighted Cox model and estimate the covariance by Taylor linearization
# after aggregating weighted dfbeta contributions at the NHANES PSU level.
design_cox <- function(formula, data, label) {
  needed <- unique(c(all.vars(formula), "WTMEC2YR", "SDMVSTRA", "SDMVPSU"))
  z <- data[complete.cases(data[, needed]), ]
  z$w_norm <- z$WTMEC2YR / mean(z$WTMEC2YR)
  fit <- coxph(formula, data = z, weights = w_norm, ties = "breslow", x = TRUE)

  db <- residuals(fit, type = "dfbeta", weighted = TRUE)
  if (is.null(dim(db))) db <- matrix(db, ncol = 1L)
  colnames(db) <- names(coef(fit))

  psu_group <- interaction(z$SDMVSTRA, z$SDMVPSU, drop = TRUE)
  agg <- rowsum(db, psu_group, reorder = FALSE)
  psu_stratum <- tapply(z$SDMVSTRA, psu_group, function(x) x[1])
  stratum <- as.numeric(psu_stratum[rownames(agg)])

  V <- matrix(0, ncol(db), ncol(db), dimnames = list(colnames(db), colnames(db)))
  for (h in unique(stratum)) {
    u <- agg[stratum == h, , drop = FALSE]
    m <- nrow(u)
    if (m > 1L) {
      uc <- sweep(u, 2L, colMeans(u))
      V <- V + m / (m - 1) * crossprod(uc)
    }
  }

  design_df <- nrow(unique(z[, c("SDMVSTRA", "SDMVPSU")])) -
    length(unique(z$SDMVSTRA))
  list(fit = fit, V = V, data = z, label = label, design_df = design_df)
}

extract_term <- function(obj, term) {
  b <- coef(obj$fit)[term]
  se <- sqrt(obj$V[term, term])
  data.frame(
    analysis = obj$label,
    N = nrow(obj$data),
    deaths = sum(obj$data$event),
    term = term,
    beta = unname(b),
    SE = unname(se),
    HR = exp(b),
    CI_low = exp(b - 1.96 * se),
    CI_high = exp(b + 1.96 * se),
    p = 2 * pnorm(-abs(b / se)),
    row.names = NULL
  )
}

linear_contrast <- function(obj, weights, label) {
  b <- coef(obj$fit)
  coef_names <- intersect(names(b), rownames(obj$V))
  b <- b[coef_names]
  V <- obj$V[coef_names, coef_names, drop = FALSE]
  missing_terms <- setdiff(names(weights), coef_names)
  if (length(missing_terms)) {
    stop("Contrast term(s) absent from fitted model: ", paste(missing_terms, collapse = ", "))
  }
  a <- setNames(rep(0, length(b)), names(b))
  a[names(weights)] <- weights
  est <- sum(a * b)
  se <- sqrt(drop(t(a) %*% V %*% a))
  data.frame(
    subgroup = label,
    beta = est,
    SE = se,
    HR = exp(est),
    CI_low = exp(est - 1.96 * se),
    CI_high = exp(est + 1.96 * se),
    p = 2 * pnorm(-abs(est / se)),
    row.names = NULL
  )
}

wald_test <- function(obj, terms) {
  b <- coef(obj$fit)[terms]
  V <- obj$V[terms, terms, drop = FALSE]
  q <- length(terms)
  wald <- drop(t(b) %*% solve(V, b))
  F_stat <- wald / q
  data.frame(
    test_df = q,
    denominator_df = obj$design_df,
    Wald_chisq = wald,
    F = F_stat,
    p = pf(F_stat, q, obj$design_df, lower.tail = FALSE)
  )
}

# -------------------------------------------------------------------------
# 1. Restricted cubic spline for maximum single-trial grip strength
# -------------------------------------------------------------------------

knots <- weighted_quantile(cc$grip_max, cc$WTMEC2YR, c(0.05, 0.35, 0.65, 0.95))
reference_grip <- weighted_quantile(cc$grip_max, cc$WTMEC2YR, 0.50)
plot_limits <- weighted_quantile(cc$grip_max, cc$WTMEC2YR, c(0.025, 0.975))

cc_basis <- rcs_basis(cc$grip_max, knots)
cc <- cbind(cc, cc_basis)

f_rcs <- Surv(time, event) ~ vcf_grade2 + grip_rcs_linear + grip_rcs_nl1 +
  grip_rcs_nl2 + age10 + female + race + education + INDFMPIR + bmi + smoking +
  alcohol + phys_active + diabetes + hypertension + cvd + cancer

rcs_fit <- design_cox(f_rcs, cc, "Four-knot restricted cubic spline")
rcs_terms <- c("grip_rcs_linear", "grip_rcs_nl1", "grip_rcs_nl2")
rcs_overall <- wald_test(rcs_fit, rcs_terms)
rcs_overall$test <- "Overall association"
rcs_nonlinear <- wald_test(rcs_fit, c("grip_rcs_nl1", "grip_rcs_nl2"))
rcs_nonlinear$test <- "Nonlinearity"
rcs_tests <- bind_rows(rcs_overall, rcs_nonlinear) %>%
  select(test, everything())

grid <- data.frame(grip_max = seq(plot_limits[1], plot_limits[2], length.out = 240L))
grid_basis <- rcs_basis(grid$grip_max, knots)
ref_basis <- drop(rcs_basis(reference_grip, knots))
coef_rcs <- coef(rcs_fit$fit)[rcs_terms]
V_rcs <- rcs_fit$V[rcs_terms, rcs_terms, drop = FALSE]

contrast_matrix <- sweep(grid_basis, 2L, ref_basis, "-")
grid$logHR <- drop(contrast_matrix %*% coef_rcs)
grid$SE <- sqrt(rowSums((contrast_matrix %*% V_rcs) * contrast_matrix))
grid$HR <- exp(grid$logHR)
grid$CI_low <- exp(grid$logHR - 1.96 * grid$SE)
grid$CI_high <- exp(grid$logHR + 1.96 * grid$SE)
grid$reference_grip <- reference_grip

rcs_metadata <- data.frame(
  item = c(
    "N", "Deaths", "Design degrees of freedom", "Reference grip (kg)",
    "Knot 1 (kg; weighted 5th percentile)",
    "Knot 2 (kg; weighted 35th percentile)",
    "Knot 3 (kg; weighted 65th percentile)",
    "Knot 4 (kg; weighted 95th percentile)",
    "Plot lower limit (kg; weighted 2.5th percentile)",
    "Plot upper limit (kg; weighted 97.5th percentile)"
  ),
  value = c(
    nrow(cc), sum(cc$event), rcs_fit$design_df, reference_grip,
    knots, plot_limits
  )
)

write.csv(grid, file.path(out_dir, "grip_rcs_curve_estimates.csv"), row.names = FALSE)
write.csv(rcs_tests, file.path(out_dir, "grip_rcs_wald_tests.csv"), row.names = FALSE)
write.csv(rcs_metadata, file.path(out_dir, "grip_rcs_metadata.csv"), row.names = FALSE)

palette <- c(
  signal_blue = "#2C7FB8",
  signal_fill = "#A6CEE3",
  neutral_dark = "#222222",
  neutral_mid = "#777777"
)

p_rcs <- ggplot(grid, aes(x = grip_max, y = HR)) +
  geom_ribbon(
    aes(ymin = CI_low, ymax = CI_high),
    fill = palette[["signal_fill"]], alpha = 0.55, linewidth = 0
  ) +
  geom_line(colour = palette[["signal_blue"]], linewidth = 0.85) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35,
             colour = palette[["neutral_mid"]]) +
  geom_vline(xintercept = reference_grip, linetype = "dotted", linewidth = 0.35,
             colour = palette[["neutral_dark"]]) +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 4, 8),
    labels = c("0.5", "1", "2", "4", "8"),
    expand = expansion(mult = c(0.03, 0.09))
  ) +
  scale_x_continuous(
    breaks = seq(15, 65, by = 10),
    limits = range(grid$grip_max),
    expand = expansion(mult = c(0, 0))
  ) +
  annotate(
    "label", x = min(grid$grip_max) + 1.5, y = max(grid$CI_high) / 1.25,
    hjust = 0, vjust = 1, size = 2.35, linewidth = 0,
    label = paste0(
      "Overall P = ", formatC(
        rcs_tests$p[rcs_tests$test == "Overall association"],
        digits = 3, format = "f"
      ),
      "\nP for nonlinearity = ", formatC(
        rcs_tests$p[rcs_tests$test == "Nonlinearity"],
        digits = 3, format = "f"
      )
    )
  ) +
  labs(
    x = "Maximum single-trial grip strength (kg)",
    y = "Hazard ratio for all-cause mortality"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8.2, base_family = "Helvetica") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = 8.4, colour = "black"),
    axis.text = element_text(size = 7.5, colour = "black"),
    plot.margin = margin(6, 9, 6, 7),
    panel.grid = element_blank()
  )

# -------------------------------------------------------------------------
# 2. Limited exploratory sex and age subgroup analyses
# -------------------------------------------------------------------------

base_cov_without_sex <- paste(
  "+ age10 + race + education + INDFMPIR + bmi + smoking + alcohol + phys_active",
  "+ diabetes + hypertension + cvd + cancer"
)
f_sex <- as.formula(paste(
  "Surv(time,event) ~ vcf_grade2*female + low_grip*female",
  base_cov_without_sex
))
sex_fit <- design_cox(f_sex, cc, "Sex interaction model")

base_cov_age <- paste(
  "+ age10 + female + race + education + INDFMPIR + bmi + smoking + alcohol",
  "+ phys_active + diabetes + hypertension + cvd + cancer"
)
f_age <- as.formula(paste(
  "Surv(time,event) ~ vcf_grade2*older65 + low_grip*older65",
  base_cov_age
))
age_fit <- design_cox(f_age, cc, "Age interaction model")

interaction_name <- function(obj, term_a, term_b) {
  candidates <- c(paste0(term_a, ":", term_b), paste0(term_b, ":", term_a))
  found <- intersect(candidates, names(coef(obj$fit)))
  if (length(found) != 1L) {
    stop("Could not uniquely identify interaction between ", term_a, " and ", term_b)
  }
  found
}

sex_vcf_int <- interaction_name(sex_fit, "vcf_grade2", "female")
sex_grip_int <- interaction_name(sex_fit, "low_grip", "female")
age_vcf_int <- interaction_name(age_fit, "vcf_grade2", "older65")
age_grip_int <- interaction_name(age_fit, "low_grip", "older65")

f_primary <- Surv(time, event) ~ vcf_grade2 + low_grip + age10 + female + race +
  education + INDFMPIR + bmi + smoking + alcohol + phys_active + diabetes +
  hypertension + cvd + cancer
primary_fit <- design_cox(f_primary, cc, "Overall")

overall_vcf <- linear_contrast(primary_fit, c(vcf_grade2 = 1), "Overall")
overall_grip <- linear_contrast(primary_fit, c(low_grip = 1), "Overall")

sex_vcf <- bind_rows(
  linear_contrast(sex_fit, c(vcf_grade2 = 1), "Men"),
  linear_contrast(
    sex_fit,
    setNames(c(1, 1), c("vcf_grade2", sex_vcf_int)),
    "Women"
  )
)
sex_grip <- bind_rows(
  linear_contrast(sex_fit, c(low_grip = 1), "Men"),
  linear_contrast(
    sex_fit,
    setNames(c(1, 1), c("low_grip", sex_grip_int)),
    "Women"
  )
)

age_vcf <- bind_rows(
  linear_contrast(age_fit, c(vcf_grade2 = 1), "Age 40-64 years"),
  linear_contrast(
    age_fit,
    setNames(c(1, 1), c("vcf_grade2", age_vcf_int)),
    "Age >=65 years"
  )
)
age_grip <- bind_rows(
  linear_contrast(age_fit, c(low_grip = 1), "Age 40-64 years"),
  linear_contrast(
    age_fit,
    setNames(c(1, 1), c("low_grip", age_grip_int)),
    "Age >=65 years"
  )
)

subgroup <- bind_rows(
  mutate(overall_vcf, exposure = "Grade >=2 VCF", modifier = "Overall"),
  mutate(overall_grip, exposure = "FNIH low grip strength", modifier = "Overall"),
  mutate(sex_vcf, exposure = "Grade >=2 VCF", modifier = "Sex"),
  mutate(sex_grip, exposure = "FNIH low grip strength", modifier = "Sex"),
  mutate(age_vcf, exposure = "Grade >=2 VCF", modifier = "Age"),
  mutate(age_grip, exposure = "FNIH low grip strength", modifier = "Age")
)

count_for <- function(label, data) {
  data.frame(
    subgroup = label,
    N = nrow(data),
    deaths = sum(data$event),
    vcf_N = sum(data$vcf_grade2 == 1),
    vcf_deaths = sum(data$event[data$vcf_grade2 == 1]),
    low_grip_N = sum(data$low_grip == 1),
    low_grip_deaths = sum(data$event[data$low_grip == 1])
  )
}
counts <- bind_rows(
  count_for("Overall", cc),
  count_for("Men", cc[cc$female == 0, ]),
  count_for("Women", cc[cc$female == 1, ]),
  count_for("Age 40-64 years", cc[cc$older65 == 0, ]),
  count_for("Age >=65 years", cc[cc$older65 == 1, ])
)
subgroup <- left_join(subgroup, counts, by = "subgroup")

interaction_p <- data.frame(
  exposure = rep(c("Grade >=2 VCF", "FNIH low grip strength"), each = 2),
  modifier = rep(c("Sex", "Age"), 2),
  interaction_term = c(
    sex_vcf_int, age_vcf_int,
    sex_grip_int, age_grip_int
  ),
  p_interaction = c(
    extract_term(sex_fit, sex_vcf_int)$p,
    extract_term(age_fit, age_vcf_int)$p,
    extract_term(sex_fit, sex_grip_int)$p,
    extract_term(age_fit, age_grip_int)$p
  )
)
subgroup <- left_join(subgroup, interaction_p, by = c("exposure", "modifier"))

write.csv(subgroup, file.path(out_dir, "exploratory_subgroup_estimates.csv"), row.names = FALSE)
write.csv(interaction_p, file.path(out_dir, "exploratory_interaction_pvalues.csv"), row.names = FALSE)
write.csv(counts, file.path(out_dir, "exploratory_subgroup_event_counts.csv"), row.names = FALSE)

subgroup_order <- c("Overall", "Men", "Women", "Age 40-64 years", "Age >=65 years")
subgroup$subgroup <- factor(subgroup$subgroup, levels = rev(subgroup_order))
subgroup$exposure <- factor(
  subgroup$exposure,
  levels = c("Grade >=2 VCF", "FNIH low grip strength")
)
levels(subgroup$exposure) <- c("Grade 2 or higher VCF", "FNIH low grip strength")
levels(subgroup$subgroup) <- c(
  "Age 65 years or older", "Age 40-64 years", "Women", "Men", "Overall"
)

pd <- position_dodge(width = 0.52)
p_subgroup <- ggplot(
  subgroup,
  aes(x = HR, y = subgroup, colour = exposure, shape = exposure)
) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35,
             colour = palette[["neutral_mid"]]) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y", width = 0.15, position = pd, linewidth = 0.55
  ) +
  geom_point(position = pd, size = 2.25, stroke = 0.65) +
  scale_x_log10(
    limits = c(0.1, 10),
    breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
    labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8")
  ) +
  scale_colour_manual(values = c(
    "Grade 2 or higher VCF" = "#2C7FB8",
    "FNIH low grip strength" = "#D95F0E"
  )) +
  scale_shape_manual(values = c(
    "Grade 2 or higher VCF" = 16,
    "FNIH low grip strength" = 17
  )) +
  labs(
    x = "Hazard ratio for all-cause mortality (95% CI)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 8.2, base_family = "Helvetica") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.title = element_text(size = 8.4, colour = "black"),
    axis.text = element_text(size = 7.5, colour = "black"),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 7.2),
    plot.margin = margin(6, 9, 6, 7),
    panel.grid = element_blank()
  )

save_plot <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  # The standard PDF device preserves searchable/editable Helvetica text and
  # avoids transformed 1-point glyph runs produced by some Windows Cairo builds.
  grDevices::pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in,
    family = "Helvetica", onefile = TRUE, useDingbats = FALSE
  )
  print(plot)
  grDevices::dev.off()

  # A Cairo PDF is also retained for journals that prefer Cairo-embedded fonts.
  grDevices::cairo_pdf(
    paste0(stem, "_cairo.pdf"), width = width_in, height = height_in,
    family = "Arial", onefile = TRUE
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_in, height = height_in,
    units = "in", res = dpi, compression = "lzw"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300
  )
  print(plot)
  grDevices::dev.off()

  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(
      paste0(stem, ".svg"), width = width_in, height = height_in,
      bg = "white"
    )
    print(plot)
    grDevices::dev.off()
  } else if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::CairoSVG(
      paste0(stem, ".svg"), width = width_in, height = height_in,
      bg = "white", onefile = TRUE
    )
    print(plot)
    grDevices::dev.off()
  }
}

save_plot(
  p_rcs,
  file.path(fig_dir, "Figure3_grip_strength_rcs"),
  width_mm = 183, height_mm = 120
)
save_plot(
  p_subgroup,
  file.path(fig_dir, "FigureS5_exploratory_subgroups"),
  width_mm = 183, height_mm = 120
)

session <- capture.output(sessionInfo())
writeLines(session, file.path(out_dir, "analysis_08_session_info.txt"))

cat("Complete-case cohort:", nrow(cc), "participants and", sum(cc$event), "deaths\n")
cat("RCS knots (kg):", paste(sprintf("%.1f", knots), collapse = ", "), "\n")
cat("RCS reference (kg):", sprintf("%.1f", reference_grip), "\n")
print(rcs_tests)
print(interaction_p)
print(subgroup[, c(
  "exposure", "modifier", "subgroup", "N", "deaths", "HR", "CI_low", "CI_high",
  "p", "p_interaction"
)])

