project_root <- normalizePath(".", winslash="/", mustWork=TRUE)
dir.create(file.path(project_root, "results"), showWarnings=FALSE)
evalue_rr <- function(rr) rr + sqrt(rr * (rr - 1))
primary <- data.frame(
  exposure=c("Grade >=2 vertebral compression fracture", "FNIH-defined low grip strength"),
  HR=c(1.61,1.74), CI_low=c(1.05,1.20), CI_high=c(2.47,2.52)
)
primary$E_value_point <- evalue_rr(primary$HR)
primary$E_value_CI_closest_to_null <- evalue_rr(primary$CI_low)
write.csv(primary,file.path(project_root,"results","E_values_primary_analysis.csv"),row.names=FALSE)
print(primary)
